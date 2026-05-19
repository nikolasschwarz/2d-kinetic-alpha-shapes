#include "SegmentBuilder.hpp"

#include "KineticDelaunayFlipEvent.hpp"
#include "KineticDelaunayRadiusEvent.hpp"
#include "Logger.hpp"
#include "SegmentBuilderCrossingCallback.hpp"
#include "SegmentBuilderFlipCallback.hpp"
#include "SegmentBuilderRadiusCallback.hpp"
#include "SegmentBuilderSectionCallback.hpp"
#include "SegmentBuilderSubdivisionCallback.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iomanip>
#include <limits>
#include <glm/gtx/exterior_product.hpp>
#include <map>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <optional>
#include <sstream>
#include <stdexcept>

using namespace kinDS;

namespace
{
void interpolateFlexibleVerticesAlongEdge(
  VoronoiMesh& mesh, std::vector<int>& flex, size_t anchor_old_vertex, size_t anchor_new_vertex)
{
  if (flex.empty())
  {
    return;
  }
  auto& verts = mesh.getVertices();
  if (anchor_old_vertex >= verts.size() || anchor_new_vertex >= verts.size())
  {
    return;
  }
  const glm::dvec3 p0 = verts[anchor_old_vertex];
  const glm::dvec3 p1 = verts[anchor_new_vertex];
  const double z0 = p0.z;
  const double z1 = p1.z;
  const double denom = z1 - z0;
  const size_t k = flex.size();
  for (size_t j = 0; j < k; ++j)
  {
    const int fj = flex[j];
    if (fj < 0)
    {
      continue;
    }
    const size_t fju = static_cast<size_t>(fj);
    if (fju >= verts.size())
    {
      continue;
    }
    const double fz = verts[fju].z;
    double s = 0.0;
    if (std::abs(denom) > 1e-18)
    {
      s = (fz - z0) / denom;
      if (s < 0.0)
      {
        s = 0.0;
      }
      else if (s > 1.0)
      {
        s = 1.0;
      }
    }
    else
    {
      // Anchors share z: cannot parameterize by z; fall back to uniform spacing along the segment in xy.
      s = static_cast<double>(j + 1) / static_cast<double>(k + 1);
    }
    const double x = p0.x + s * (p1.x - p0.x);
    const double y = p0.y + s * (p1.y - p0.y);
    mesh.replaceVertex(fju, glm::dvec3(x, y, fz));
  }
}

const char* boundaryEventTypeToString(SegmentBuilder::BoundaryEventType event_type)
{
  switch (event_type)
  {
  case SegmentBuilder::BoundaryEventType::Init:
    return "init";
  case SegmentBuilder::BoundaryEventType::Section:
    return "section_event";
  case SegmentBuilder::BoundaryEventType::Radius:
    return "radius_event";
  case SegmentBuilder::BoundaryEventType::Crossing:
    return "crossing_event";
  case SegmentBuilder::BoundaryEventType::Subdivision:
    return "subdivision_event";
  default:
    return "unknown_event";
  }
}

const char* boundarySegmentActionToString(SegmentBuilder::BoundarySegmentAction segment_action)
{
  switch (segment_action)
  {
  case SegmentBuilder::BoundarySegmentAction::NewSegment:
    return "new_segment";
  case SegmentBuilder::BoundarySegmentAction::SegmentCompleted:
    return "segment_completed";
  case SegmentBuilder::BoundarySegmentAction::SegmentRemoved:
    return "segment_removed";
  case SegmentBuilder::BoundarySegmentAction::SegmentRemapped:
    return "segment_remapped";
  default:
    return "unknown_action";
  }
}

std::string makeBoundaryMeshMetadata(
  SegmentBuilder::BoundaryEventType event_type, SegmentBuilder::BoundarySegmentAction segment_action)
{
  std::ostringstream o;
  o << "{\"event_type\":\"" << boundaryEventTypeToString(event_type) << "\",\"segment_action\":\""
    << boundarySegmentActionToString(segment_action) << "\"}";
  return o.str();
}

std::string appendBoundaryVertexPosMetadata(const std::string& base_meta, const char* pos)
{
  if (base_meta.empty() || base_meta.back() != '}')
  {
    return base_meta;
  }
  std::string out = base_meta;
  out.pop_back();
  out += ",\"pos\":\"";
  out += pos;
  out += "\"}";
  return out;
}

std::optional<std::pair<size_t, size_t>> closingMeshIntersectionListIndices(const KineticDelaunay& kd,
  KineticDelaunay::CrossingData::EdgeIntersectionRef ref)
{
  const auto& crossing_data = kd.getCrossingData();
  const auto& ir = *ref;
  const size_t d_edge = ir.delaunay_edge_id;
  const size_t v_edge = ir.voronoi_edge_id;
  if (d_edge >= crossing_data.delaunay_edge_intersections.size()
    || v_edge >= crossing_data.voronoi_edge_intersections.size())
  {
    return std::nullopt;
  }
  size_t d_list_idx = 0;
  for (auto d_it = crossing_data.delaunay_edge_intersections[d_edge].begin();
    d_it != crossing_data.delaunay_edge_intersections[d_edge].end(); ++d_it, ++d_list_idx)
  {
    const KineticDelaunay::CrossingData::EdgeIntersectionRef r = *d_it;
    if (&(*r) != &ir)
    {
      continue;
    }
    size_t v_list_idx = 0;
    for (auto v_it = crossing_data.voronoi_edge_intersections[v_edge].begin();
      v_it != crossing_data.voronoi_edge_intersections[v_edge].end(); ++v_it, ++v_list_idx)
    {
      if (&(**v_it) == &ir)
      {
        return std::make_pair(d_list_idx, v_list_idx);
      }
    }
    break;
  }
  return std::nullopt;
}

void logClosingMeshIntersectionVertex(const KineticDelaunay& kd, KineticDelaunay::CrossingData::EdgeIntersectionRef ref)
{
  const auto& ir = *ref;
  const auto idx = closingMeshIntersectionListIndices(kd, ref);
  if (idx.has_value())
  {
    KINDS_DEBUG("Adding vertex at intersection of Delaunay edge " << ir.delaunay_edge_id << " and Voronoi edge "
                                                                 << ir.voronoi_edge_id << " with intersection indices ("
                                                                 << idx->first << "," << idx->second << ")");
  }
  else
  {
    KINDS_DEBUG("Adding vertex at intersection of Delaunay edge " << ir.delaunay_edge_id << " and Voronoi edge "
                                                                 << ir.voronoi_edge_id << " with intersection indices (?,?)");
  }
}

void logClosingMeshStripEndpointVertex(const KineticDelaunay& kd, size_t segment_index,
  const std::optional<KineticDelaunay::CrossingData::EdgeIntersectionRef>& endpoint_crossing,
  size_t closing_voronoi_edge_id)
{
  if (endpoint_crossing.has_value())
  {
    logClosingMeshIntersectionVertex(kd, endpoint_crossing.value());
  }
  else if (closing_voronoi_edge_id == static_cast<size_t>(-1))
  {
    KINDS_DEBUG("Polygon walk: using mesh vertex on segment " << segment_index << " of Voronoi edge ? (strip circumcenter).");
  }
  else
  {
    KINDS_DEBUG("Polygon walk: using mesh vertex on segment " << segment_index << " of Voronoi edge "
                                                             << closing_voronoi_edge_id << " (strip circumcenter).");
  }
}

void logExtractionClosingMeshStripEndpoint(const KineticDelaunay& kd, double t, size_t extraction_index,
  const std::optional<KineticDelaunay::CrossingData::EdgeIntersectionRef>& endpoint_crossing,
  size_t closing_voronoi_edge_id)
{
  if (endpoint_crossing.has_value())
  {
    const auto ref = endpoint_crossing.value();
    const auto& ir = *ref;
    const auto idx = closingMeshIntersectionListIndices(kd, ref);
    glm::dvec2 xy(0.0, 0.0);
    const bool have_xy = tryComputeCrossingIntersectionPosition2D(kd, endpoint_crossing, t, xy);
    if (idx.has_value())
    {
      if (have_xy)
      {
        KINDS_DEBUG("Closing mesh extraction: strip endpoint at Delaunay edge " << ir.delaunay_edge_id << " / Voronoi edge "
                                                                               << ir.voronoi_edge_id << ", list indices ("
                                                                               << idx->first << "," << idx->second
                                                                               << "), 2D=(" << xy.x << "," << xy.y << ")");
      }
      else
      {
        KINDS_DEBUG("Closing mesh extraction: strip endpoint at Delaunay edge " << ir.delaunay_edge_id << " / Voronoi edge "
                                                                               << ir.voronoi_edge_id << ", list indices ("
                                                                               << idx->first << "," << idx->second << ")");
      }
    }
    else if (have_xy)
    {
      KINDS_DEBUG("Closing mesh extraction: strip endpoint at Delaunay edge " << ir.delaunay_edge_id << " / Voronoi edge "
                                                                             << ir.voronoi_edge_id
                                                                             << ", list indices (?,?), 2D=(" << xy.x
                                                                             << "," << xy.y << ")");
    }
    else
    {
      KINDS_DEBUG("Closing mesh extraction: strip endpoint at Delaunay edge " << ir.delaunay_edge_id << " / Voronoi edge "
                                                                           << ir.voronoi_edge_id << ", list indices (?,?)");
    }
  }
  else if (closing_voronoi_edge_id == static_cast<size_t>(-1))
  {
    KINDS_DEBUG("Closing mesh extraction: strip endpoint on extraction segment " << extraction_index
                                                                                << " of Voronoi edge ? (circumcenter, "
                                                                                   "no crossing ref).");
  }
  else
  {
    KINDS_DEBUG("Closing mesh extraction: strip endpoint on extraction segment " << extraction_index << " of Voronoi edge "
                                                                                << closing_voronoi_edge_id
                                                                                << " (circumcenter, no crossing ref).");
  }
}

std::string makeRegularStripVertexMetadataJson(double kinetic_time, size_t voronoi_edge_id, size_t even_half_edge_id,
  int strand_even_origin, int strand_odd_origin, SegmentBuilder::BoundaryEventType event_type,
  SegmentBuilder::BoundarySegmentAction segment_action,
  const std::optional<KineticDelaunay::CrossingData::EdgeIntersectionRef>& crossing, const char* pos, const char* op)
{
  std::ostringstream o;
  o << std::setprecision(17);
  o << "{\"mesh_type\":\"regular\",\"event_type\":\"" << boundaryEventTypeToString(event_type) << "\",\"segment_action\":\""
    << boundarySegmentActionToString(segment_action) << "\",\"time\":" << kinetic_time << ",\"voronoi_edge_id\":" << voronoi_edge_id
    << ",\"even_half_edge_id\":" << even_half_edge_id << ",\"strand_even_origin\":" << strand_even_origin
    << ",\"strand_odd_origin\":" << strand_odd_origin;
  if (crossing.has_value())
  {
    const auto& ir = *crossing.value();
    o << ",\"crossing_delaunay_edge_id\":" << ir.delaunay_edge_id << ",\"crossing_voronoi_edge_id\":" << ir.voronoi_edge_id
      << ",\"crossing_delaunay_param\":" << ir.delaunay_edge_param;
  }
  o << ",\"pos\":\"" << pos << "\"";
  if (op != nullptr && op[0] != '\0')
  {
    o << ",\"op\":\"" << op << "\"";
  }
  o << "}";
  return o.str();
}

std::string makeRegularStripFaceMetadataJson(double kinetic_time, size_t voronoi_edge_id, size_t even_half_edge_id,
  int strand_even_origin, int strand_odd_origin, SegmentBuilder::BoundaryEventType event_type,
  SegmentBuilder::BoundarySegmentAction segment_action, const char* op)
{
  std::ostringstream o;
  o << std::setprecision(17);
  o << "{\"mesh_type\":\"regular\",\"event_type\":\"" << boundaryEventTypeToString(event_type) << "\",\"segment_action\":\""
    << boundarySegmentActionToString(segment_action) << "\",\"time\":" << kinetic_time << ",\"voronoi_edge_id\":" << voronoi_edge_id
    << ",\"even_half_edge_id\":" << even_half_edge_id << ",\"strand_even_origin\":" << strand_even_origin
    << ",\"strand_odd_origin\":" << strand_odd_origin << ",\"op\":\"" << op << "\"}";
  return o.str();
}
} // namespace

static bool raySegmentIntersection(
  const glm::dvec2& C, const glm::dvec2& D, const glm::dvec2& A, const glm::dvec2& B, double& t_out)
{
  glm::dvec2 E = B - A;
  glm::dvec2 AC = A - C;

  double det = glm::cross(D, E);
  if (std::abs(det) < 1e-12)
    return false;

  double t = glm::cross(AC, E) / det;
  double u = glm::cross(AC, D) / det;

  if (t >= 0.0 && u >= 0.0 && u <= 1.0)
  {
    t_out = t;
    return true;
  }
  return false;
}

static std::vector<double> rayCast(
  const std::vector<BoundaryPoint>& polygon, const glm::dvec2& origin, const glm::dvec2& dir)
{
  double lenCP = glm::length(dir);

  if (lenCP < 1e-12)
    return {};

  std::vector<double> hits;
  const size_t n = polygon.size();
  for (size_t i = 0; i < n; ++i)
  {
    const glm::dvec2& A = polygon[i].p;
    const glm::dvec2& B = polygon[(i + 1) % n].p;

    double t;
    if (raySegmentIntersection(origin, dir, A, B, t))
    {
      hits.emplace_back(t);
    }
  }

  return hits;
}

static std::vector<double> rayCast(
  const std::vector<glm::dvec2>& polygon, const glm::dvec2& origin, const glm::dvec2& dir)
{
  double lenCP = glm::length(dir);

  if (lenCP < 1e-12)
    return {};

  std::vector<double> hits;
  const size_t n = polygon.size();
  for (size_t i = 0; i < n; ++i)
  {
    const glm::dvec2& A = polygon[i];
    const glm::dvec2& B = polygon[(i + 1) % n];

    double t;
    if (raySegmentIntersection(origin, dir, A, B, t))
    {
      hits.emplace_back(t);
    }
  }

  return hits;
}

static double relativeDistanceFromCenter(
  const std::vector<BoundaryPoint>& polygon, const glm::dvec2& center, const glm::dvec2& point)
{
  glm::dvec2 dir = point - center;
  auto hits = rayCast(polygon, center, dir);

  if (hits.empty())
    return std::numeric_limits<double>::quiet_NaN();

  double t_max = -std::numeric_limits<double>::infinity();

  for (double t : hits)
  {
    t_max = std::max(t_max, t);
  }

  if (t_max < 0)
  {
    return std::numeric_limits<double>::quiet_NaN();
  }

  // |B - C| = t_max * |D|
  return 1.0 / t_max;
}

static bool isInside(const std::vector<BoundaryPoint>& polygon, const glm::dvec2& center, const glm::dvec2& point)
{
  double rel_dist = relativeDistanceFromCenter(polygon, center, point);
  if (std::isnan(rel_dist))
  {
    throw std::runtime_error("Center lies outside of polygon");
  }
  return rel_dist <= 1.0;
}

[[nodiscard]] glm::dvec3 kinDS::SegmentBuilder::computeVoronoiVertex(size_t half_edge_id, double t) const
{
  return kin_del.computeVoronoiVertexClampedInfinity(half_edge_id, t);
}

bool clampVoronoiVertices(glm::dvec3& left_vertex, glm::dvec3& right_vertex,
  const std::vector<BoundaryPoint>& boundary_points, const glm::dvec2& centroid)
{
  // don't do this for now
  return true;

  bool left_inside = isInside(boundary_points, centroid, glm::dvec2 { left_vertex[0], left_vertex[1] });
  bool right_inside = isInside(boundary_points, centroid, glm::dvec2 { right_vertex[0], right_vertex[1] });

  if (left_inside && !right_inside)
  {
    KINDS_DEBUG(
      "Clamping right vertex: (" << right_vertex[0] << ", " << right_vertex[1] << ", " << right_vertex[2] << ")");
    // clamp right to the boundary
    glm::dvec2 origin { left_vertex[0], left_vertex[1] };
    glm::dvec2 dir { right_vertex[0] - left_vertex[0], right_vertex[1] - left_vertex[1] };
    auto hits = rayCast(boundary_points, origin, dir);

    double s_min = std::numeric_limits<double>::infinity();

    for (double s : hits)
    {
      if (s >= 0)
      {
        s_min = std::min(s_min, s);
      }
    }

    auto right_vertex_2d = origin + s_min * dir;
    right_vertex = glm::dvec3 { right_vertex_2d[0], right_vertex_2d[1], right_vertex[2] };
    KINDS_DEBUG(
      "Clamped right vertex to: (" << right_vertex[0] << ", " << right_vertex[1] << ", " << right_vertex[2] << ")");
  }
  else if (!left_inside && right_inside)
  {
    KINDS_DEBUG("Clamping left vertex: (" << left_vertex[0] << ", " << left_vertex[1] << ", " << left_vertex[2] << ")");
    // clamp left to the boundary
    glm::dvec2 origin { right_vertex[0], right_vertex[1] };
    glm::dvec2 dir { left_vertex[0] - right_vertex[0], left_vertex[1] - right_vertex[1] };
    auto hits = rayCast(boundary_points, origin, dir);

    double s_min = std::numeric_limits<double>::infinity();

    for (double s : hits)
    {
      if (s >= 0)
      {
        s_min = std::min(s_min, s);
      }
    }

    auto left_vertex_2d = origin + s_min * dir;
    left_vertex = glm::dvec3 { left_vertex_2d[0], left_vertex_2d[1], left_vertex[2] };
    KINDS_DEBUG(
      "Clamped left vertex to: (" << left_vertex[0] << ", " << left_vertex[1] << ", " << left_vertex[2] << ")");
  }
  else if (!left_inside && !right_inside)
  {
    // I'm not sure yet if this will work, especially for connections, but for now we will just discard it
    return false;
  }
  return true;
}

void kinDS::SegmentBuilder::meshletDiagnosticLogLine(const char* tag, size_t half_edge_id, double t, const char* extra_note) const
{
  const size_t even = half_edge_id & ~1;
  const size_t dual_edge = even / 2;
  size_t pi = static_cast<size_t>(-1);
  if (even < half_edge_index_to_segment_mesh_pair_index.size())
  {
    pi = half_edge_index_to_segment_mesh_pair_index[even];
  }
  size_t nv = 0;
  size_t nt = 0;
  size_t nstrip = 0;
  size_t nclosed = 0;
  const bool pi_ok = (pi < meshes.size()) && (pi < segment_mesh_pair_last_left_and_right_vertex.size());
  if (pi_ok)
  {
    nv = meshes[pi].getVertices().size();
    nt = meshes[pi].getTriangleCount();
    nstrip = segment_mesh_pair_last_left_and_right_vertex[pi].size();
    for (const auto& s : segment_mesh_pair_last_left_and_right_vertex[pi])
    {
      if (s.mesh_start_vertex_id >= 0 && s.mesh_end_vertex_id >= 0)
      {
        ++nclosed;
      }
    }
  }
  std::ostringstream oss;
  oss << "meshlet_diag " << tag << " he=" << half_edge_id << " even=" << even << " dual_edge=" << dual_edge << " t=" << t;
  if (even < half_edge_index_to_segment_mesh_pair_index.size() && !pi_ok)
  {
    oss << " pair=INVALID(raw_pi=" << pi << ")";
  }
  else
  {
    oss << " pair=" << (pi_ok ? static_cast<long long>(pi) : -1);
  }
  oss << " verts=" << nv << " tris=" << nt << " strips=" << nstrip << " closed_strips=" << nclosed;
  if (extra_note != nullptr && extra_note[0] != '\0')
  {
    oss << ' ' << extra_note;
  }
  KINDS_DEBUG(oss.str());
}

void kinDS::SegmentBuilder::meshletDiagnosticWarnIfUnexpectedEmptyAfterStartNewMesh(size_t half_edge_even, double t,
  bool initial_left_inside, const VoronoiMesh& mesh, const std::list<MeshingData>& strips) const
{
  const size_t dual_edge = half_edge_even / 2;
  size_t with_any_mesh_id = 0;
  for (const auto& s : strips)
  {
    if (s.mesh_start_vertex_id >= 0 || s.mesh_end_vertex_id >= 0)
    {
      ++with_any_mesh_id;
    }
  }
  const size_t nv = mesh.getVertices().size();
  if (nv == 0 && initial_left_inside)
  {
    KINDS_WARNING("meshlet_diag unexpected_empty after startNewMesh: dual_edge=" << dual_edge << " t=" << t
      << " — left circumcenter face was inside alpha-shape but mesh has no vertices (strips=" << strips.size() << ").");
  }
  if (nv == 0 && with_any_mesh_id > 0)
  {
    KINDS_WARNING("meshlet_diag unexpected_empty after startNewMesh: dual_edge=" << dual_edge << " t=" << t
      << " — strip entries reference mesh vertex ids but mesh has no vertices (strips=" << strips.size() << ").");
  }
  if (nv > 0 && strips.empty())
  {
    KINDS_WARNING("meshlet_diag inconsistent after startNewMesh: dual_edge=" << dual_edge << " t=" << t
      << " — mesh has " << nv << " vertices but strip list is empty.");
  }
}

SegmentBuilder::RegularMeshStripIntervalEndpoints SegmentBuilder::regularMeshStripIntervalFromMeshingData(
  const MeshingData& segment, size_t even_half_edge_id, size_t odd_half_edge_id)
{
  RegularMeshStripIntervalEndpoints interval;
  if (segment.start_half_edge_id < 0)
  {
    interval.start_open_voronoi_half_edge_id = even_half_edge_id;
  }
  else
  {
    interval.start_crossing = segment.start_crossing;
    interval.start_crossed_inside_half_edge_id = segment.start_half_edge_id;
  }
  if (segment.end_half_edge_id < 0)
  {
    interval.end_open_voronoi_half_edge_id = odd_half_edge_id;
  }
  else
  {
    interval.end_crossing = segment.end_crossing;
    interval.end_crossed_inside_half_edge_id = segment.end_half_edge_id;
  }
  return interval;
}

glm::dvec3 SegmentBuilder::regularMeshStripIntervalEndpointPositionAt(const RegularMeshStripIntervalEndpoints& interval,
  bool at_start, size_t even_half_edge_id, size_t odd_half_edge_id, size_t voronoi_edge_id, double t,
  const RadiusBoundaryTransitionShiftContext* boundary_transition_shift) const
{
  const char* endpoint_label = at_start ? "start" : "end";
  if (at_start)
  {
    if (interval.start_open_voronoi_half_edge_id.has_value())
    {
      const glm::dvec3 p = computeVoronoiVertex(interval.start_open_voronoi_half_edge_id.value(), t);
      if (!std::isfinite(p.x) || !std::isfinite(p.y) || !std::isfinite(p.z) || (p.x == 0.0 && p.y == 0.0))
      {
        KINDS_WARNING("regularMeshStripIntervalEndpointPositionAt(" << endpoint_label
          << "): degenerate circumcenter endpoint for voronoi_edge=" << voronoi_edge_id << " t=" << t << " -> (" << p.x
          << ", " << p.y << ", " << p.z << ").");
      }
      return p;
    }
    if (interval.start_crossing.has_value())
    {
      const KineticDelaunay::CrossingData::EdgeIntersectionRef orig_ref = interval.start_crossing.value();
      const glm::dvec3 old_chord_pos
        = closingMeshVoronoiDelaunayCrossingPosition(t, orig_ref->voronoi_edge_id, orig_ref->delaunay_edge_id);
      KineticDelaunay::CrossingData::EdgeIntersectionRef position_ref = orig_ref;
      if (auto neighbor_opt = neighborIntersectionOnTargetAlongVoronoiEdge(orig_ref, voronoi_edge_id, boundary_transition_shift);
        neighbor_opt.has_value())
      {
        position_ref = neighbor_opt.value();
      }
      glm::dvec2 xy;
      if (tryComputeCrossingIntersectionPosition2D(kin_del, position_ref, t, xy))
      {
        const glm::dvec3 p(xy, t);
        if (!std::isfinite(p.x) || !std::isfinite(p.y) || !std::isfinite(p.z) || (p.x == 0.0 && p.y == 0.0))
        {
          KINDS_WARNING("regularMeshStripIntervalEndpointPositionAt(" << endpoint_label
            << "): degenerate CrossingData endpoint for voronoi_edge=" << voronoi_edge_id << " t=" << t << " -> (" << p.x
            << ", " << p.y << ", " << p.z << ").");
        }
        if (position_ref != orig_ref)
        {
          logRadiusBoundaryTransitionVertexShift(
            "finishRegularMeshStripInterval_endpoint", t, orig_ref, position_ref, old_chord_pos, p);
        }
        return p;
      }
      if (boundary_transition_shift != nullptr)
      {
        const glm::dvec3 p = closingMeshVoronoiDelaunayCrossingPosition(
          t, position_ref->voronoi_edge_id, position_ref->delaunay_edge_id);
        if (position_ref != orig_ref)
        {
          logRadiusBoundaryTransitionVertexShift(
            "finishRegularMeshStripInterval_endpoint", t, orig_ref, position_ref, old_chord_pos, p);
        }
        return p;
      }
      const size_t delaunay_edge_id = static_cast<size_t>(interval.start_crossed_inside_half_edge_id) / 2;
      return closingMeshVoronoiDelaunayCrossingPosition(t, voronoi_edge_id, delaunay_edge_id);
    }
  }
  else
  {
    if (interval.end_open_voronoi_half_edge_id.has_value())
    {
      const glm::dvec3 p = computeVoronoiVertex(interval.end_open_voronoi_half_edge_id.value(), t);
      if (!std::isfinite(p.x) || !std::isfinite(p.y) || !std::isfinite(p.z) || (p.x == 0.0 && p.y == 0.0))
      {
        KINDS_WARNING("regularMeshStripIntervalEndpointPositionAt(" << endpoint_label
          << "): degenerate circumcenter endpoint for voronoi_edge=" << voronoi_edge_id << " t=" << t << " -> (" << p.x
          << ", " << p.y << ", " << p.z << ").");
      }
      return p;
    }
    if (interval.end_crossing.has_value())
    {
      const KineticDelaunay::CrossingData::EdgeIntersectionRef orig_ref = interval.end_crossing.value();
      const glm::dvec3 old_chord_pos
        = closingMeshVoronoiDelaunayCrossingPosition(t, orig_ref->voronoi_edge_id, orig_ref->delaunay_edge_id);
      KineticDelaunay::CrossingData::EdgeIntersectionRef position_ref = orig_ref;
      if (auto neighbor_opt = neighborIntersectionOnTargetAlongVoronoiEdge(orig_ref, voronoi_edge_id, boundary_transition_shift);
        neighbor_opt.has_value())
      {
        position_ref = neighbor_opt.value();
      }
      glm::dvec2 xy;
      if (tryComputeCrossingIntersectionPosition2D(kin_del, position_ref, t, xy))
      {
        const glm::dvec3 p(xy, t);
        if (position_ref != orig_ref)
        {
          logRadiusBoundaryTransitionVertexShift(
            "finishRegularMeshStripInterval_endpoint", t, orig_ref, position_ref, old_chord_pos, p);
        }
        return p;
      }
      if (boundary_transition_shift != nullptr)
      {
        const glm::dvec3 p = closingMeshVoronoiDelaunayCrossingPosition(
          t, position_ref->voronoi_edge_id, position_ref->delaunay_edge_id);
        if (position_ref != orig_ref)
        {
          logRadiusBoundaryTransitionVertexShift(
            "finishRegularMeshStripInterval_endpoint", t, orig_ref, position_ref, old_chord_pos, p);
        }
        return p;
      }
      const size_t delaunay_edge_id = static_cast<size_t>(interval.end_crossed_inside_half_edge_id) / 2;
      return closingMeshVoronoiDelaunayCrossingPosition(t, voronoi_edge_id, delaunay_edge_id);
    }
  }
  return glm::dvec3(0.0);
}

void SegmentBuilder::finishRegularMeshStripInterval(VoronoiMesh& mesh, const std::vector<BoundaryPoint>& boundary_polygon,
  const glm::dvec2& centroid, size_t even_half_edge_id, size_t voronoi_edge_id, double t, size_t strand_vertex_id,
  int strand_even_origin_i, int strand_odd_origin_i, BoundaryEventType event_type, BoundarySegmentAction segment_action,
  const RegularMeshStripIntervalEndpoints& interval, size_t last_start_vertex_index, size_t last_end_vertex_index,
  const std::string& finish_face_metadata, const RadiusBoundaryTransitionShiftContext* boundary_transition_shift,
  size_t& out_start_vertex_index, size_t& out_end_vertex_index)
{
  const size_t odd_half_edge_id = even_half_edge_id + 1;
  const glm::dvec3 new_start_pos
    = regularMeshStripIntervalEndpointPositionAt(interval, true, even_half_edge_id, odd_half_edge_id, voronoi_edge_id, t,
      boundary_transition_shift);
  const glm::dvec3 new_end_pos = regularMeshStripIntervalEndpointPositionAt(interval, false, even_half_edge_id,
    odd_half_edge_id, voronoi_edge_id, t, boundary_transition_shift);

  const auto& graph = kin_del.getGraph();
  const std::optional<size_t> start_vv = interval.start_open_voronoi_half_edge_id.has_value()
    ? std::optional<size_t>(graph.getHalfEdges()[even_half_edge_id].face)
    : std::nullopt;
  const std::optional<size_t> end_vv = interval.end_open_voronoi_half_edge_id.has_value()
    ? std::optional<size_t>(graph.getHalfEdges()[odd_half_edge_id].face)
    : std::nullopt;

  const std::string meta_start = makeRegularStripVertexMetadataJson(t, voronoi_edge_id, even_half_edge_id,
    strand_even_origin_i, strand_odd_origin_i, event_type, segment_action, interval.start_crossing, "left", nullptr);
  const std::string meta_end = makeRegularStripVertexMetadataJson(t, voronoi_edge_id, even_half_edge_id,
    strand_even_origin_i, strand_odd_origin_i, event_type, segment_action, interval.end_crossing, "right", nullptr);

  out_start_vertex_index
    = addMeshletVertex(mesh, boundary_polygon, centroid, new_start_pos, strand_vertex_id, t, start_vv, meta_start);
  out_end_vertex_index
    = addMeshletVertex(mesh, boundary_polygon, centroid, new_end_pos, strand_vertex_id, t, end_vv, meta_end);

  if (last_start_vertex_index == last_end_vertex_index)
  {
    addMeshletTriangle(mesh, out_start_vertex_index, last_end_vertex_index, out_end_vertex_index, finish_face_metadata);
  }
  else if (mesh.getVertices()[last_start_vertex_index][2] < mesh.getVertices()[last_end_vertex_index][2])
  {
    addMeshletTriangle(mesh, last_start_vertex_index, last_end_vertex_index, out_start_vertex_index, finish_face_metadata);
    addMeshletTriangle(mesh, out_start_vertex_index, last_end_vertex_index, out_end_vertex_index, finish_face_metadata);
  }
  else
  {
    addMeshletTriangle(mesh, last_start_vertex_index, last_end_vertex_index, out_end_vertex_index, finish_face_metadata);
    addMeshletTriangle(mesh, last_start_vertex_index, out_end_vertex_index, out_start_vertex_index, finish_face_metadata);
  }
}

/**
 * @brief Extends every active strip on one Voronoi edge to time @p t (emit quads as two triangles per strip).
 *
 * @details Expects @ref startNewMesh to have run earlier on this edge so @c last_segments holds seeded corners.
 * Per strip, the “before” state is (@c last_left, @c last_right) = (@c mesh_start_vertex_id, @c mesh_end_vertex_id).
 * @ref finishRegularMeshStripInterval places new corners at the interval endpoints at @p t, connects them to the old
 * corners, then overwrites @c mesh_start_vertex_id / @c mesh_end_vertex_id so a later finish continues from the new front.
 * Strips with either corner unset (@c -1) are skipped (e.g. open strip never fully seeded).
 */
void kinDS::SegmentBuilder::finishMesh(size_t he_id, double t, const std::vector<BoundaryPoint>& boundary_points,
  BoundaryEventType event_type, BoundarySegmentAction segment_action,
  const RadiusBoundaryTransitionShiftContext* boundary_transition_shift)
{
  if (kin_del.getGraph().isInfinite(he_id) && kin_del.computeBoundaryOnTheFly())
  {
    meshletDiagnosticLogLine("finish_mesh_skip", he_id, t, "reason=infinite_boundary_on_the_fly");
    return;
  }

  size_t segment_mesh_pair_index = half_edge_index_to_segment_mesh_pair_index[he_id];
  if (segment_mesh_pair_index >= meshes.size() || segment_mesh_pair_index >= segment_mesh_pair_last_left_and_right_vertex.size())
  {
    KINDS_WARNING("meshlet_diag finish_mesh: invalid pair index he=" << he_id << " t=" << t << " pair=" << segment_mesh_pair_index);
    return;
  }

  meshletDiagnosticLogLine("finish_mesh_enter", he_id, t, "");

  // Live strip state for this dual edge (same list @ref startNewMesh populated).
  auto& last_segments = segment_mesh_pair_last_left_and_right_vertex[segment_mesh_pair_index];
  if (last_segments.empty())
  {
    meshletDiagnosticLogLine("finish_mesh_noop", he_id, t, "last_segments empty (no extension)");
    return;
  }

  VoronoiMesh& mesh = meshes[segment_mesh_pair_index];
  const size_t even_id = he_id & ~1;
  const size_t odd_id = even_id + 1;
  const size_t voronoi_edge_id = even_id / 2;
  auto& he = kin_del.getGraph().getHalfEdges()[even_id];
  const auto& twin_he_finish = kin_del.getGraph().getHalfEdges()[odd_id];
  const int strand_even_origin_i = static_cast<int>(he.origin);
  const int strand_odd_origin_i = static_cast<int>(twin_he_finish.origin);
  glm::dvec2 centroid = polygonCentroid(boundary_points);

  const std::string finish_face_meta = makeRegularStripFaceMetadataJson(t, voronoi_edge_id, even_id, strand_even_origin_i,
    strand_odd_origin_i, event_type, segment_action, "finish_strip");

  if (mesh.getVertices().empty())
  {
    size_t closed_strip_refs = 0;
    for (const auto& s : last_segments)
    {
      if (s.mesh_start_vertex_id >= 0 && s.mesh_end_vertex_id >= 0)
      {
        ++closed_strip_refs;
      }
    }
    if (closed_strip_refs > 0)
    {
      KINDS_WARNING("meshlet_diag finish_mesh_enter: mesh has no vertices but " << closed_strip_refs
        << " strip(s) have mesh_start/end set (expected non-empty mesh) dual_edge=" << voronoi_edge_id << " he=" << he_id
        << " t=" << t << ".");
    }
  }

  // Strand site used as @c origin for addMeshletVertex on this edge (finite endpoint of the even directed edge).
  size_t v = he.origin;
  if (v == -1)
  {
    v = kin_del.getGraph().destination(even_id);
  }

  // Snapshot interval topology from current MeshingData (crossing refs refreshed after any CrossingData rebuild).
  std::vector<RegularMeshStripIntervalEndpoints> finish_intervals;
  finish_intervals.reserve(last_segments.size());
  for (auto& segment : last_segments)
  {
    if (segment.mesh_start_vertex_id < 0 || segment.mesh_end_vertex_id < 0)
    {
      continue;
    }
    refreshMeshingDataCrossingRefs(segment, voronoi_edge_id);
    finish_intervals.push_back(regularMeshStripIntervalFromMeshingData(segment, even_id, odd_id));
  }

  const size_t processable_strips = finish_intervals.size();
  size_t loops_ran = 0;
  size_t tris_before = mesh.getTriangleCount();

  // Parallel walk: one finish_interval per processable strip, in list order.
  auto interval_it = finish_intervals.begin();
  for (auto& segment : last_segments)
  {
    if (segment.mesh_start_vertex_id < 0 || segment.mesh_end_vertex_id < 0)
    {
      continue;
    }

    ++loops_ran;
    const RegularMeshStripIntervalEndpoints& interval = *interval_it++;
    const size_t last_left = static_cast<size_t>(segment.mesh_start_vertex_id);
    const size_t last_right = static_cast<size_t>(segment.mesh_end_vertex_id);
    size_t new_start_vertex_index = 0;
    size_t new_end_vertex_index = 0;
    finishRegularMeshStripInterval(mesh, boundary_points, centroid, even_id, voronoi_edge_id, t, v, strand_even_origin_i,
      strand_odd_origin_i, event_type, segment_action, interval, last_left, last_right, finish_face_meta,
      boundary_transition_shift, new_start_vertex_index, new_end_vertex_index);
    // After finish: strip front is at the new corners; next event’s finish will use these as last_left / last_right.
    segment.mesh_start_vertex_id = static_cast<int>(new_start_vertex_index);
    segment.mesh_end_vertex_id = static_cast<int>(new_end_vertex_index);
  }

  const size_t tris_after = mesh.getTriangleCount();
  const long long d_tris = static_cast<long long>(tris_after) - static_cast<long long>(tris_before);
  {
    std::ostringstream oss;
    oss << "tris_before=" << tris_before << " tris_after=" << tris_after << " d_tris=" << d_tris
        << " processable_strips=" << processable_strips << " loops_ran=" << loops_ran;
    meshletDiagnosticLogLine("finish_mesh_exit", he_id, t, oss.str().c_str());
  }
  if (processable_strips > 0 && tris_after == tris_before)
  {
    KINDS_WARNING("meshlet_diag finish_mesh: expected triangle extension (processable_strips=" << processable_strips
      << ") but triangle count unchanged for dual_edge=" << voronoi_edge_id << " he=" << he_id << " t=" << t << ".");
  }
  // Strips may exist in the list but be “open” on one side (only one seeded corner) — nothing to extrude this step.
  if (!last_segments.empty() && processable_strips == 0)
  {
    KINDS_DEBUG("meshlet_diag finish_mesh: strips present but none had both mesh_start/end set — no extension for he="
      << he_id << " dual_edge=" << voronoi_edge_id << " t=" << t << ".");
  }
}

SegmentBuilder::SegmentBuilder(
  KineticDelaunay& kin_del, std::vector<std::pair<size_t, double>> subdivisions, bool create_transformed_mesh)
  : kin_del(kin_del)
  , section_callback_(std::make_unique<SegmentBuilderSectionCallback>(*this))
  , flip_callback_(std::make_unique<SegmentBuilderFlipCallback>(*this))
  , radius_callback_(std::make_unique<SegmentBuilderRadiusCallback>(*this))
  , crossing_callback_(std::make_unique<SegmentBuilderCrossingCallback>(*this))
  , subdivision_callback_(std::make_unique<SegmentBuilderSubdivisionCallback>(*this))
  , create_transformed_mesh(create_transformed_mesh)
  , boundary_mesh({ "bark", "interior" })
{
  assert(std::is_sorted(subdivisions.begin(), subdivisions.end(),
    [](const auto& a, const auto& b) { return a.second < b.second; }));
  kin_del.setSubdivisionSchedule(std::move(subdivisions));
}

SegmentBuilder::SegmentBuilder(KineticDelaunay& kin_del, bool create_transformed_mesh)
  : kin_del(kin_del)
  , section_callback_(std::make_unique<SegmentBuilderSectionCallback>(*this))
  , flip_callback_(std::make_unique<SegmentBuilderFlipCallback>(*this))
  , radius_callback_(std::make_unique<SegmentBuilderRadiusCallback>(*this))
  , crossing_callback_(std::make_unique<SegmentBuilderCrossingCallback>(*this))
  , subdivision_callback_(std::make_unique<SegmentBuilderSubdivisionCallback>(*this))
  , create_transformed_mesh(create_transformed_mesh)
  , boundary_mesh({ "bark", "interior" })
{
  kin_del.setSubdivisionSchedule({});
}

SegmentBuilder::~SegmentBuilder() = default;

static glm::dvec2 intersectSegments(
  const glm::dvec2& p, const glm::dvec2& p2, const glm::dvec2& q, const glm::dvec2& q2)
{
  glm::dvec2 r = p2 - p;
  glm::dvec2 s = q2 - q;

  double rxs = r.x * s.y - r.y * s.x;
  double qpxr = (q - p).x * r.y - (q - p).y * r.x;

  if (rxs == 0.0 && qpxr == 0.0)
    return glm::dvec2(std::numeric_limits<double>::quiet_NaN()); // collinear

  if (rxs == 0.0 && qpxr != 0.0)
    return glm::dvec2(std::numeric_limits<double>::infinity()); // parallel

  double t = ((q - p).x * s.y - (q - p).y * s.x) / rxs;
  double u = ((q - p).x * r.y - (q - p).y * r.x) / rxs;

  if (t < 0 || t > 1 || u < 0 || u > 1)
  {
    KINDS_WARNING("Segments do not intersect within the segment bounds");
  }

  return p + t * r;
}

std::vector<SegmentBuilder::RegularMeshStripIntervalEndpoints> SegmentBuilder::collectRegularMeshStripIntervalsOnVoronoiEdge(
  size_t even_half_edge_id, size_t voronoi_edge_id, size_t left_containing_tri_id) const
{
  std::vector<RegularMeshStripIntervalEndpoints> intervals;
  const size_t odd_half_edge_id = even_half_edge_id + 1;
  const auto& graph = kin_del.getGraph();
  const size_t left_voronoi_vertex_id = graph.getHalfEdges()[even_half_edge_id].face;
  const size_t right_voronoi_vertex_id = graph.getHalfEdges()[odd_half_edge_id].face;

  bool inside = kin_del.getFaceInside(left_containing_tri_id);
  std::optional<RegularMeshStripIntervalEndpoints> open_interval;

  if (inside || !kin_del.computeBoundaryOnTheFly())
  {
    RegularMeshStripIntervalEndpoints interval;
    interval.start_open_voronoi_half_edge_id = even_half_edge_id;
    open_interval = interval;
  }

  if (kin_del.computeBoundaryOnTheFly())
  {
    const KineticDelaunay::CrossingData& crossing_data = kin_del.getCrossingData();
    if (voronoi_edge_id >= crossing_data.voronoi_edge_intersections.size())
    {
      std::ostringstream oss;
      oss << "collectRegularMeshStripIntervalsOnVoronoiEdge: voronoi_edge_intersections has no slot for voronoi_edge_id="
          << voronoi_edge_id << " (even_half_edge_id=" << even_half_edge_id << ", voronoi_vertices=[" << left_voronoi_vertex_id
          << "," << right_voronoi_vertex_id << "], left_containing_tri_id=" << left_containing_tri_id
          << ", voronoi_edge_intersections.size=" << crossing_data.voronoi_edge_intersections.size() << ")";
      throw std::runtime_error(oss.str());
    }

    const auto& v_intersections = crossing_data.voronoi_edge_intersections[voronoi_edge_id];
    std::vector<std::pair<size_t, KineticDelaunay::CrossingData::EdgeIntersectionRef>> directed_crossings_with_ref;
    directed_crossings_with_ref.reserve(v_intersections.size());
    size_t current_face_id = left_containing_tri_id;
    size_t step_index = 0;
    for (const KineticDelaunay::CrossingData::EdgeIntersectionRef& iref : v_intersections)
    {
      const size_t d = iref->delaunay_edge_id;
      const size_t he0 = d * 2;
      const size_t he1 = he0 + 1;
      size_t crossed_he_id = static_cast<size_t>(-1);
      if (graph.getHalfEdges()[he0].face == current_face_id)
      {
        crossed_he_id = he0;
      }
      else if (graph.getHalfEdges()[he1].face == current_face_id)
      {
        crossed_he_id = he1;
      }
      else
      {
        std::ostringstream oss;
        oss << "collectRegularMeshStripIntervalsOnVoronoiEdge: cannot orient crossing along voronoi_edge_intersections at step "
            << step_index << " (voronoi_edge_id=" << voronoi_edge_id << ", even_half_edge_id=" << even_half_edge_id
            << ", delaunay_edge_id=" << d << ", current_face_id=" << current_face_id << ", face(he0)="
            << graph.getHalfEdges()[he0].face << ", face(he1)=" << graph.getHalfEdges()[he1].face
            << ", left_containing_tri_id=" << left_containing_tri_id << ", voronoi_vertices=[" << left_voronoi_vertex_id
            << "," << right_voronoi_vertex_id << "])";
        throw std::runtime_error(oss.str());
      }
      directed_crossings_with_ref.emplace_back(crossed_he_id, iref);
      current_face_id = graph.getHalfEdges()[crossed_he_id ^ 1].face;
      ++step_index;
    }

    for (const auto& entry : directed_crossings_with_ref)
    {
      const size_t crossed_he_id = entry.first;
      const auto& crossing_ref = entry.second;
      const size_t next_face_id = graph.getHalfEdges()[crossed_he_id ^ 1].face;
      const bool next_inside = kin_del.getFaceInside(next_face_id);

      if (inside != next_inside)
      {
        if (inside)
        {
          if (open_interval.has_value())
          {
            open_interval->end_crossing = crossing_ref;
            open_interval->end_crossed_inside_half_edge_id = static_cast<int>(crossed_he_id);
            intervals.push_back(open_interval.value());
            open_interval.reset();
          }
        }
        else
        {
          RegularMeshStripIntervalEndpoints interval;
          interval.start_crossing = crossing_ref;
          interval.start_crossed_inside_half_edge_id = static_cast<int>(crossed_he_id ^ 1);
          open_interval = interval;
        }
      }

      inside = next_inside;
    }
  }

  if (inside || !kin_del.computeBoundaryOnTheFly())
  {
    if (open_interval.has_value())
    {
      open_interval->end_open_voronoi_half_edge_id = odd_half_edge_id;
      intervals.push_back(open_interval.value());
    }
  }

  return intervals;
}

SegmentBuilder::MeshingData SegmentBuilder::meshRegularStripInterval(VoronoiMesh& mesh,
  const std::vector<BoundaryPoint>& boundary_polygon, const glm::dvec2& centroid, size_t even_half_edge_id,
  size_t voronoi_edge_id, double t, int strand_even_origin_i, int strand_odd_origin_i, BoundaryEventType event_type,
  BoundarySegmentAction segment_action, const RegularMeshStripIntervalEndpoints& interval,
  const RadiusBoundaryTransitionShiftContext* boundary_transition_shift)
{
  const auto& graph = kin_del.getGraph();

  auto crossing_position = [&](KineticDelaunay::CrossingData::EdgeIntersectionRef crossing_ref) -> glm::dvec3
  {
    KineticDelaunay::CrossingData::EdgeIntersectionRef position_ref = crossing_ref;
    const glm::dvec3 old_chord_pos
      = closingMeshVoronoiDelaunayCrossingPosition(t, crossing_ref->voronoi_edge_id, crossing_ref->delaunay_edge_id);
    if (auto neighbor_opt
      = neighborIntersectionOnTargetAlongVoronoiEdge(crossing_ref, voronoi_edge_id, boundary_transition_shift);
      neighbor_opt.has_value())
    {
      position_ref = neighbor_opt.value();
    }
    const glm::dvec3 intersection_point = closingMeshVoronoiDelaunayCrossingPosition(
      t, position_ref->voronoi_edge_id, position_ref->delaunay_edge_id);
    if (position_ref != crossing_ref)
    {
      logRadiusBoundaryTransitionVertexShift(
        "meshRegularStripInterval_crossing", t, crossing_ref, position_ref, old_chord_pos, intersection_point);
    }
    return intersection_point;
  };

  int start_half_edge_id = -1;
  size_t start_vertex_index = 0;
  if (interval.start_crossing.has_value())
  {
    start_half_edge_id = interval.start_crossed_inside_half_edge_id;
    const std::string meta_start = makeRegularStripVertexMetadataJson(t, voronoi_edge_id, even_half_edge_id,
      strand_even_origin_i, strand_odd_origin_i, event_type, segment_action,
      std::optional<KineticDelaunay::CrossingData::EdgeIntersectionRef>(interval.start_crossing), "left", nullptr);
    start_vertex_index = addMeshletVertex(mesh, boundary_polygon, centroid, crossing_position(interval.start_crossing.value()),
      static_cast<size_t>(-1), t, std::nullopt, meta_start);
  }
  else if (interval.start_open_voronoi_half_edge_id.has_value())
  {
    const size_t open_he = interval.start_open_voronoi_half_edge_id.value();
    const glm::dvec3 pos = computeVoronoiVertex(open_he, t);
    const size_t voronoi_vertex_id = graph.getHalfEdges()[open_he].face;
    const int origin = graph.getHalfEdges()[open_he].origin;
    const std::string meta_start = makeRegularStripVertexMetadataJson(t, voronoi_edge_id, even_half_edge_id,
      strand_even_origin_i, strand_odd_origin_i, event_type, segment_action, std::nullopt, "left", nullptr);
    start_vertex_index
      = addMeshletVertex(mesh, boundary_polygon, centroid, pos, origin, t, std::optional<size_t>(voronoi_vertex_id), meta_start);
  }

  int end_half_edge_id = -1;
  size_t end_vertex_index = 0;
  if (interval.end_crossing.has_value())
  {
    end_half_edge_id = interval.end_crossed_inside_half_edge_id;
    const std::string meta_end = makeRegularStripVertexMetadataJson(t, voronoi_edge_id, even_half_edge_id,
      strand_even_origin_i, strand_odd_origin_i, event_type, segment_action,
      std::optional<KineticDelaunay::CrossingData::EdgeIntersectionRef>(interval.end_crossing), "right", nullptr);
    end_vertex_index = addMeshletVertex(mesh, boundary_polygon, centroid, crossing_position(interval.end_crossing.value()),
      static_cast<size_t>(-1), t, std::nullopt, meta_end);
  }
  else if (interval.end_open_voronoi_half_edge_id.has_value())
  {
    const size_t open_he = interval.end_open_voronoi_half_edge_id.value();
    const glm::dvec3 pos = computeVoronoiVertex(open_he, t);
    const size_t voronoi_vertex_id = graph.getHalfEdges()[open_he].face;
    const int origin = graph.getHalfEdges()[open_he].origin;
    const std::string meta_end = makeRegularStripVertexMetadataJson(t, voronoi_edge_id, even_half_edge_id,
      strand_even_origin_i, strand_odd_origin_i, event_type, segment_action, std::nullopt, "right", nullptr);
    end_vertex_index
      = addMeshletVertex(mesh, boundary_polygon, centroid, pos, origin, t, std::optional<size_t>(voronoi_vertex_id), meta_end);
  }

  MeshingData segment { static_cast<int>(start_vertex_index), static_cast<int>(end_vertex_index), start_half_edge_id,
    end_half_edge_id };
  segment.start_crossing = interval.start_crossing;
  segment.end_crossing = interval.end_crossing;
  return segment;
}

/**
 * @brief (Re)builds the regular strip meshlet(s) on one undirected Voronoi / Delaunay dual edge at time @p t.
 *
 * @details One call corresponds to “open” or “re-open” meshing on this edge. Each inside-alpha interval along the edge
 * becomes one @ref MeshingData strip stored in @c segment_mesh_pair_last_left_and_right_vertex[pair]. A strip’s state
 * after this function (before any @ref finishMesh):
 *   - @c mesh_start_vertex_id / @c mesh_end_vertex_id: mesh indices of the **left** and **right** strip corners seeded
 *     at interval endpoints (circumcenter or boundary crossing).
 *   - @c start_half_edge_id / @c end_half_edge_id: inside-oriented Delaunay half-edge at each boundary endpoint, or @c -1
 *     for an open circumcenter end.
 *   - @c start_crossing / @c end_crossing: iterators into @c CrossingData when the endpoint is a stored crossing (refreshed
 *     after collection).
 * No quads are emitted here; @ref finishMesh advances the corners and adds triangles between old and new positions.
 *
 * @param half_edge_id Directed Delaunay half-edge on the strand side used for lookup (even or odd of the dual edge).
 * @param reuse_existing_pair_and_mesh If true, append vertices into the mesh already tied to this edge and replace the
 *   strip list; if false, allocate a new @c SegmentMeshPair and mesh when none exists yet.
 */
void SegmentBuilder::startNewMesh(size_t half_edge_id, double t, bool reuse_existing_pair_and_mesh,
  BoundaryEventType event_type, BoundarySegmentAction segment_action,
  const RadiusBoundaryTransitionShiftContext* boundary_transition_shift)
{
  if (kin_del.getGraph().isInfinite(half_edge_id) && kin_del.computeBoundaryOnTheFly())
  {
    meshletDiagnosticLogLine("start_new_mesh_skip", half_edge_id, t, "reason=infinite_boundary_on_the_fly");
    return;
  }

  // Canonical even/odd pair for the undirected dual edge: even = left Voronoi vertex, odd = right Voronoi vertex.
  size_t even_id = half_edge_id & ~1;
  size_t odd_id = even_id + 1;
  {
    std::ostringstream oss;
    oss << "reuse_existing_pair=" << (reuse_existing_pair_and_mesh ? "true" : "false");
    meshletDiagnosticLogLine("start_new_mesh_begin", half_edge_id, t, oss.str().c_str());
  }

  const auto& graph = kin_del.getGraph();
  const auto& he = graph.getHalfEdges()[even_id];
  const auto& twin_he = graph.getHalfEdges()[odd_id];

  // --- SegmentMeshPair: book-keeping linking this Voronoi-edge mesh to the two incident strand segments (sites). ---
  size_t segment_mesh_pair_index = static_cast<size_t>(-1);
  bool created_new_pair = false;
  if (reuse_existing_pair_and_mesh && even_id < half_edge_index_to_segment_mesh_pair_index.size())
  {
    segment_mesh_pair_index = half_edge_index_to_segment_mesh_pair_index[even_id];
  }

  MeshStructure::SegmentMeshPair segment_mesh_pair;
  segment_mesh_pair.segment_index0 = he.origin == -1 ? -1 : strand_to_segment_indices[he.origin].back();
  segment_mesh_pair.segment_index1 = twin_he.origin == -1 ? -1 : strand_to_segment_indices[twin_he.origin].back();
  if (segment_mesh_pair_index == static_cast<size_t>(-1))
  {
    created_new_pair = true;
    segment_mesh_pair_index = segment_mesh_pairs.size();
    half_edge_index_to_segment_mesh_pair_index[even_id] = segment_mesh_pair_index;
    half_edge_index_to_segment_mesh_pair_index[odd_id] = segment_mesh_pair_index;
    segment_mesh_pairs.push_back(segment_mesh_pair);
  }
  else
  {
    half_edge_index_to_segment_mesh_pair_index[even_id] = segment_mesh_pair_index;
    half_edge_index_to_segment_mesh_pair_index[odd_id] = segment_mesh_pair_index;
    segment_mesh_pairs[segment_mesh_pair_index] = segment_mesh_pair;
  }

  KINDS_DEBUG("Using mesh pair index " << segment_mesh_pair_index << " for half-edge " << even_id << " (segment indices "
    << segment_mesh_pair.segment_index0 << ", " << segment_mesh_pair.segment_index1 << ")");

  // Component boundary + centroid for alpha-shape checks and texture coordinates on new vertices.
  size_t vertex = std::max(he.origin, twin_he.origin);
  size_t component_id = kin_del.component_data.component_map[vertex];

  std::vector<bool> he_visited(graph.getHalfEdges().size(), false);
  updateBoundary(t, he_visited, component_id);

  auto& boundary_polygon = kin_del.component_data.component_boundaries[component_id][0];
  auto& centroid = kin_del.component_data.component_centroids[component_id];

  // --- Target mesh: new pair → local mesh registered at end; reuse → append to existing VoronoiMesh for this pair. ---
  VoronoiMesh mesh_local;
  const bool reuse_in_place = !created_new_pair && segment_mesh_pair_index < meshes.size();
  VoronoiMesh& mesh = reuse_in_place ? meshes[segment_mesh_pair_index] : mesh_local;
  const size_t vertex_count_before = mesh.getVertexCount();

  const size_t voronoi_edge_id = even_id / 2;
  const int strand_even_origin_i = static_cast<int>(he.origin);
  const int strand_odd_origin_i = static_cast<int>(twin_he.origin);

  // Walk the dual edge from the left circumcenter: inside intervals are strips bounded by crossings or open ends.
  const size_t left_voronoi_vertex_id = kin_del.getGraph().getHalfEdges()[even_id].face;
  const size_t left_containing_tri_id = kin_del.getCrossingDataContainingTriId(left_voronoi_vertex_id);
  const size_t right_voronoi_vertex_id = kin_del.getGraph().getHalfEdges()[odd_id].face;
  const bool initial_left_inside = kin_del.getFaceInside(left_containing_tri_id);

  const std::vector<RegularMeshStripIntervalEndpoints> strip_intervals
    = collectRegularMeshStripIntervalsOnVoronoiEdge(even_id, voronoi_edge_id, left_containing_tri_id);
  const size_t voronoi_edge_crossing_count = kin_del.computeBoundaryOnTheFly()
    && voronoi_edge_id < kin_del.getCrossingData().voronoi_edge_intersections.size()
    ? kin_del.getCrossingData().voronoi_edge_intersections[voronoi_edge_id].size()
    : 0;
  const bool ending_inside = !strip_intervals.empty() && strip_intervals.back().end_open_voronoi_half_edge_id.has_value();

  // Replace strip list for this pair: previous MeshingData (if any) is discarded; new strips hold seeded corners only.
  if (segment_mesh_pair_index >= segment_mesh_pair_last_left_and_right_vertex.size())
  {
    segment_mesh_pair_last_left_and_right_vertex.resize(segment_mesh_pair_index + 1);
  }
  segment_mesh_pair_last_left_and_right_vertex[segment_mesh_pair_index].clear();
  auto& segments_for_pair = segment_mesh_pair_last_left_and_right_vertex[segment_mesh_pair_index];

  for (const RegularMeshStripIntervalEndpoints& interval : strip_intervals)
  {
    segments_for_pair.emplace_back(meshRegularStripInterval(mesh, boundary_polygon, centroid, even_id, voronoi_edge_id, t,
      strand_even_origin_i, strand_odd_origin_i, event_type, segment_action, interval, boundary_transition_shift));
  }

  for (auto& seg : segments_for_pair)
  {
    refreshMeshingDataCrossingRefs(seg, voronoi_edge_id);
  }

  meshletDiagnosticWarnIfUnexpectedEmptyAfterStartNewMesh(even_id, t, initial_left_inside, mesh, segments_for_pair);

  const size_t vertex_count_after = mesh.getVertexCount();
  const size_t vertices_added = vertex_count_after - vertex_count_before;

  if (reuse_in_place)
  {
    if (!std::isfinite(mesh.getCreationKineticTime()))
    {
      mesh.setCreationKineticTime(t);
    }
  }
  else
  {
    registerMeshletWithSuffix(std::move(mesh_local), std::string("_voronoi") + std::to_string(voronoi_edge_id), t);
  }

  {
    std::ostringstream oss;
    oss << "created_new_pair=" << (created_new_pair ? "true" : "false") << " pair_idx=" << segment_mesh_pair_index;
    meshletDiagnosticLogLine("start_new_mesh_end", even_id, t, oss.str().c_str());
  }

  KINDS_DEBUG("startNewMesh summary: vertices_added=" << vertices_added << " mesh_vertex_count_after=" << vertex_count_after
    << " mesh_vertex_count_before=" << vertex_count_before << " half_edge_id=" << half_edge_id << " even_id=" << even_id
    << " voronoi_edge_id=" << voronoi_edge_id << " t=" << t << " pair_idx=" << segment_mesh_pair_index
    << " created_new_pair=" << (created_new_pair ? "true" : "false") << " reuse_in_place=" << (reuse_in_place ? "true" : "false")
    << " boundary_on_the_fly=" << (kin_del.computeBoundaryOnTheFly() ? "true" : "false") << " strip_segments="
    << segments_for_pair.size() << " voronoi_edge_crossing_list_size=" << voronoi_edge_crossing_count
    << " initial_left_inside=" << (initial_left_inside ? "true" : "false") << " ending_inside=" << (ending_inside ? "true" : "false")
    << " component_id=" << component_id << " voronoi_vertices=[" << left_voronoi_vertex_id << "," << right_voronoi_vertex_id << "]"
    << " segment_indices=[" << segment_mesh_pair.segment_index0 << "," << segment_mesh_pair.segment_index1 << "]");

  if (segment_mesh_pair_index >= segment_mesh_pair_last_left_and_right_vertex.size())
  {
    KINDS_WARNING("startNewMesh: strip storage out of bounds for pair index " << segment_mesh_pair_index);
  }

  assert(segment_mesh_pairs.size() == meshes.size());
}

std::string SegmentBuilder::composeBoundaryMetadata(BoundaryEventType event_type, BoundarySegmentAction segment_action)
{
  return makeBoundaryMeshMetadata(event_type, segment_action);
}

std::string SegmentBuilder::composeRegularStripVertexMetadata(double kinetic_time, size_t voronoi_edge_id,
  size_t even_half_edge_id, int strand_even_origin, int strand_odd_origin, BoundaryEventType event_type,
  BoundarySegmentAction segment_action,
  const std::optional<KineticDelaunay::CrossingData::EdgeIntersectionRef>& crossing, const char* pos, const char* op)
{
  return makeRegularStripVertexMetadataJson(kinetic_time, voronoi_edge_id, even_half_edge_id, strand_even_origin,
    strand_odd_origin, event_type, segment_action, crossing, pos, op);
}

std::string SegmentBuilder::composeRegularStripFaceMetadata(double kinetic_time, size_t voronoi_edge_id,
  size_t even_half_edge_id, int strand_even_origin, int strand_odd_origin, BoundaryEventType event_type,
  BoundarySegmentAction segment_action, const char* op)
{
  return makeRegularStripFaceMetadataJson(kinetic_time, voronoi_edge_id, even_half_edge_id, strand_even_origin,
    strand_odd_origin, event_type, segment_action, op);
}

void SegmentBuilder::logRadiusBoundaryTransitionVertexShift(const char* context, double t,
  KineticDelaunay::CrossingData::EdgeIntersectionRef from_ref, KineticDelaunay::CrossingData::EdgeIntersectionRef to_ref,
  const glm::dvec3& old_pos, const glm::dvec3& new_pos) const
{
  if (from_ref == to_ref)
  {
    return;
  }
  KINDS_DEBUG("Radius boundary transition [" << context << "]: vertex shifted from Delaunay edge "
                                              << from_ref->delaunay_edge_id << " (Voronoi edge "
                                              << from_ref->voronoi_edge_id << ") to Delaunay edge "
                                              << to_ref->delaunay_edge_id << " (Voronoi edge " << to_ref->voronoi_edge_id
                                              << ") at t=" << t << "; old=(" << old_pos.x << "," << old_pos.y << "," << old_pos.z
                                              << ") new=(" << new_pos.x << "," << new_pos.y << "," << new_pos.z << ")");
}

std::optional<KineticDelaunay::CrossingData::EdgeIntersectionRef> SegmentBuilder::neighborIntersectionOnTargetAlongVoronoiEdge(
  KineticDelaunay::CrossingData::EdgeIntersectionRef crossing_ref, size_t voronoi_edge_id,
  const RadiusBoundaryTransitionShiftContext* boundary_transition_shift) const
{
  if (!radius_boundary_transition_shift_enabled || boundary_transition_shift == nullptr
    || !boundary_transition_shift->roles_valid)
  {
    return std::nullopt;
  }

  const size_t d_edge = crossing_ref->delaunay_edge_id;
  const bool on_source = (d_edge == boundary_transition_shift->source_delaunay_edges[0]
    || d_edge == boundary_transition_shift->source_delaunay_edges[1]);
  if (!on_source)
  {
    return std::nullopt;
  }

  const size_t target_d_edge = boundary_transition_shift->target_delaunay_edge;
  const auto& crossing_data = kin_del.getCrossingData();
  if (voronoi_edge_id >= crossing_data.voronoi_edge_intersections.size())
  {
    return std::nullopt;
  }

  const auto& v_list = crossing_data.voronoi_edge_intersections[voronoi_edge_id];
  for (auto it = v_list.begin(); it != v_list.end(); ++it)
  {
    if (*it != crossing_ref)
    {
      continue;
    }
    if (it != v_list.begin())
    {
      auto pit = std::prev(it);
      if ((*pit)->delaunay_edge_id == target_d_edge)
      {
        return *pit;
      }
    }
    auto nit = std::next(it);
    if (nit != v_list.end() && (*nit)->delaunay_edge_id == target_d_edge)
    {
      return *nit;
    }
    KINDS_WARNING("neighborIntersectionOnTargetAlongVoronoiEdge: source crossing has no adjacent crossing on target "
                  "delaunay_edge="
                    << target_d_edge << " along voronoi_edge=" << voronoi_edge_id << " (source delaunay_edge=" << d_edge
                    << ").");
    return std::nullopt;
  }

  KINDS_WARNING("neighborIntersectionOnTargetAlongVoronoiEdge: crossing iterator not in voronoi_edge=" << voronoi_edge_id
                                                                                                       << " list.");
  return std::nullopt;
}

glm::dvec3 SegmentBuilder::crossingPositionWithRadiusBoundaryTransitionShift(double t,
  KineticDelaunay::CrossingData::EdgeIntersectionRef orig_ref,
  const RadiusBoundaryTransitionShiftContext* boundary_transition_shift) const
{
  KineticDelaunay::CrossingData::EdgeIntersectionRef use_ref = orig_ref;
  if (auto neighbor_opt
    = neighborIntersectionOnTargetAlongVoronoiEdge(orig_ref, orig_ref->voronoi_edge_id, boundary_transition_shift);
    neighbor_opt.has_value())
  {
    use_ref = neighbor_opt.value();
  }
  glm::dvec2 xy;
  if (tryComputeCrossingIntersectionPosition2D(kin_del, use_ref, t, xy))
  {
    return glm::dvec3(xy, t);
  }
  return closingMeshVoronoiDelaunayCrossingPosition(t, use_ref->voronoi_edge_id, use_ref->delaunay_edge_id);
}

bool SegmentBuilder::delaunayUndirectedEdgeHasVertex(
  const HalfEdgeDelaunayGraph& graph, size_t delaunay_edge_id, size_t vertex_id)
{
  const size_t he_even = 2 * delaunay_edge_id;
  if (he_even + 1 >= graph.getHalfEdges().size())
  {
    return false;
  }
  const int a = graph.getHalfEdges()[he_even].origin;
  const int b = graph.getHalfEdges()[he_even + 1].origin;
  return (a >= 0 && static_cast<size_t>(a) == vertex_id) || (b >= 0 && static_cast<size_t>(b) == vertex_id);
}

std::optional<size_t> SegmentBuilder::oppositeFiniteDelaunayVertexOnUndirectedEdge(
  const HalfEdgeDelaunayGraph& graph, size_t delaunay_edge_id, size_t site_vertex_id)
{
  const size_t he_even = 2 * delaunay_edge_id;
  if (he_even + 1 >= graph.getHalfEdges().size())
  {
    return std::nullopt;
  }
  const int a = graph.getHalfEdges()[he_even].origin;
  const int b = graph.getHalfEdges()[he_even + 1].origin;
  if (a >= 0 && static_cast<size_t>(a) == site_vertex_id)
  {
    if (b >= 0)
    {
      return static_cast<size_t>(b);
    }
    return std::nullopt;
  }
  if (b >= 0 && static_cast<size_t>(b) == site_vertex_id)
  {
    if (a >= 0)
    {
      return static_cast<size_t>(a);
    }
    return std::nullopt;
  }
  return std::nullopt;
}

std::optional<glm::dvec3> SegmentBuilder::radiusTransitionInterpolatedSitePosition(double t, size_t site_vertex_id,
  size_t strip_delaunay_edge_id, const RadiusBoundaryTransitionShiftContext* boundary_transition_shift) const
{
  if (!radius_boundary_transition_shift_enabled || boundary_transition_shift == nullptr
    || !boundary_transition_shift->roles_valid)
  {
    return std::nullopt;
  }
  const auto& graph = kin_del.getGraph();
  const size_t s0 = boundary_transition_shift->source_delaunay_edges[0];
  const size_t s1 = boundary_transition_shift->source_delaunay_edges[1];
  if (!delaunayUndirectedEdgeHasVertex(graph, s0, site_vertex_id)
    || !delaunayUndirectedEdgeHasVertex(graph, s1, site_vertex_id))
  {
    return std::nullopt;
  }
  const size_t d_strip = strip_delaunay_edge_id;
  if (d_strip != s0 && d_strip != s1)
  {
    return std::nullopt;
  }
  const size_t d_other = (d_strip == s0) ? s1 : s0;

  const auto anchor_old_and_new = [&](size_t delaunay_edge_id) -> std::optional<std::pair<glm::dvec3, glm::dvec3>>
  {
    const size_t he_even = 2 * delaunay_edge_id;
    if (he_even + 1 >= graph.getHalfEdges().size())
    {
      return std::nullopt;
    }
    const std::vector<KineticDelaunay::CrossingData::EdgeIntersectionRef> refs
      = getBoundaryIntersectionsInBoundaryOrder(delaunay_edge_id);
    if (!refs.empty())
    {
      const int o0 = graph.getHalfEdges()[he_even].origin;
      const int o1 = graph.getHalfEdges()[he_even + 1].origin;
      std::optional<KineticDelaunay::CrossingData::EdgeIntersectionRef> ref;
      if (o0 >= 0 && static_cast<size_t>(o0) == site_vertex_id)
      {
        ref = refs.front();
      }
      else if (o1 >= 0 && static_cast<size_t>(o1) == site_vertex_id)
      {
        ref = refs.back();
      }
      if (!ref.has_value())
      {
        return std::nullopt;
      }
      const glm::dvec3 p_old
        = closingMeshVoronoiDelaunayCrossingPosition(t, ref.value()->voronoi_edge_id, ref.value()->delaunay_edge_id);
      const glm::dvec3 p_new = crossingPositionWithRadiusBoundaryTransitionShift(t, ref.value(), boundary_transition_shift);
      return std::make_pair(p_old, p_new);
    }
    if (!delaunayUndirectedEdgeHasVertex(graph, delaunay_edge_id, site_vertex_id))
    {
      return std::nullopt;
    }
    if (auto v_opp = oppositeFiniteDelaunayVertexOnUndirectedEdge(graph, delaunay_edge_id, site_vertex_id);
      v_opp.has_value())
    {
      const glm::dvec2 xy = kin_del.getPointAt(t, v_opp.value());
      const glm::dvec3 p(xy, t);
      return std::make_pair(p, p);
    }
    return std::nullopt;
  };

  const auto a0 = anchor_old_and_new(d_strip);
  const auto a1 = anchor_old_and_new(d_other);
  if (!a0.has_value() || !a1.has_value())
  {
    return std::nullopt;
  }

  const glm::dvec3 p0_old = a0.value().first;
  const glm::dvec3 p1_old = a1.value().first;
  const glm::dvec2 site_xy = kin_del.getPointAt(t, site_vertex_id);
  const double d0 = glm::length(glm::dvec2(p0_old.x, p0_old.y) - site_xy);
  const double d1 = glm::length(glm::dvec2(p1_old.x, p1_old.y) - site_xy);
  const double sum = d0 + d1;
  if (!(sum > 1e-30))
  {
    return std::nullopt;
  }

  const glm::dvec3 p0_new = a0.value().second;
  const glm::dvec3 p1_new = a1.value().second;

  const double w0 = d1 / sum;
  const double w1 = d0 / sum;
  const glm::dvec3 out(w0 * p0_new.x + w1 * p1_new.x, w0 * p0_new.y + w1 * p1_new.y, t);
  KINDS_DEBUG("Radius boundary transition corner site: cell=" << site_vertex_id << " strip_d=" << d_strip << " other_d="
                                                                << d_other << " w0=" << w0 << " w1=" << w1 << " d0=" << d0
                                                                << " d1=" << d1 << " t=" << t << " out=(" << out.x << ","
                                                                << out.y << ")");
  return out;
}

size_t SegmentBuilder::resolveIntersectionMeshPairIndex(size_t voronoi_cell_id,
  std::optional<KineticDelaunay::CrossingData::EdgeIntersectionRef> start_intersection,
  std::optional<KineticDelaunay::CrossingData::EdgeIntersectionRef> end_intersection, double event_time) const
{
  if (!start_intersection.has_value() && !end_intersection.has_value())
  {
    const std::string msg
      = "resolveIntersectionMeshPairIndex: both start_intersection and end_intersection are null.";
    KINDS_ERROR(msg);
    throw std::runtime_error(msg);
  }

  if(start_intersection.has_value() && end_intersection.has_value())
  {
    if (start_intersection.value()->delaunay_edge_id != end_intersection.value()->delaunay_edge_id)
    {
      std::ostringstream oss;
      oss << "resolveIntersectionMeshPairIndex: start/end intersection Delaunay edge mismatch (start="
          << start_intersection.value()->delaunay_edge_id << ", end=" << end_intersection.value()->delaunay_edge_id
          << ", voronoi_cell_id=" << voronoi_cell_id << ", event_time=" << event_time << ").";
      KINDS_ERROR(oss.str());
      throw std::runtime_error(oss.str());
    }

    // We assume correct order here, check that
    auto& start_value = start_intersection.value();
    auto& end_value = end_intersection.value();

    if(start_value->next_segment_mesh_pair_index == end_value->prev_segment_mesh_pair_index)
    {
      return start_value->next_segment_mesh_pair_index;
    }

    // Try to recover, perhaps they are swapped, but issue a warning since this should not happen
    if(start_value->prev_segment_mesh_pair_index == end_value->next_segment_mesh_pair_index)
    {
      KINDS_WARNING("resolveIntersectionMeshPairIndex: start/end intersection mesh pair index mismatch but recoverable by swapping start/end (start_next="
        << start_value->next_segment_mesh_pair_index << ", start_prev=" << start_value->prev_segment_mesh_pair_index
        << ", end_next=" << end_value->next_segment_mesh_pair_index << ", end_prev=" << end_value->prev_segment_mesh_pair_index
        << ", voronoi_cell_id=" << voronoi_cell_id << ", event_time=" << event_time << ").");
      return start_value->prev_segment_mesh_pair_index;
    }

    std::ostringstream oss;
    oss << "resolveIntersectionMeshPairIndex: start/end intersection mesh pair index mismatch (start_next="
        << start_value->next_segment_mesh_pair_index << ", end_prev=" << end_value->prev_segment_mesh_pair_index
        << ", voronoi_cell_id=" << voronoi_cell_id << ", event_time=" << event_time << ").";
    KINDS_ERROR(oss.str());
    throw std::runtime_error(oss.str());
  }

  // Now handle cases where one is null
  if(start_intersection.has_value())
  {
    // Verify that this is indeed the last intersection in the delaunay edge
    size_t delaunay_edge_id = start_intersection.value()->delaunay_edge_id;
    auto d_ref = start_intersection.value()->delaunay_ref;
    auto next_ref = std::next(d_ref);
    if(next_ref != kin_del.getCrossingData().delaunay_edge_intersections[delaunay_edge_id].end())
    {
      std::ostringstream oss;
      oss << "resolveIntersectionMeshPairIndex: start_intersection is not the last one on its Delaunay edge (total_on_edge=" << kin_del.getCrossingData().delaunay_edge_intersections[delaunay_edge_id].size()
          << ", voronoi_cell_id=" << voronoi_cell_id << ", event_time=" << event_time << ").";
      KINDS_ERROR(oss.str());
      throw std::runtime_error(oss.str());
    }

    // check if value is valid
    if(start_intersection.value()->next_segment_mesh_pair_index == static_cast<size_t>(-1))
    {
      // print the other value for debugging
      size_t prev_index = start_intersection.value()->prev_segment_mesh_pair_index;
      std::ostringstream oss;
      oss << "resolveIntersectionMeshPairIndex: start_intersection has invalid next_segment_mesh_pair_index (voronoi_cell_id=" << voronoi_cell_id
          << ", event_time=" << event_time << ", prev_segment_mesh_pair_index=" << prev_index << ").";
      KINDS_ERROR(oss.str());
      throw std::runtime_error(oss.str());
    }

    return start_intersection.value()->next_segment_mesh_pair_index;
  }
  else // if(end_intersection.has_value())
  {
    // Verify that this is indeed the first intersection in the delaunay edge
    size_t delaunay_edge_id = end_intersection.value()->delaunay_edge_id;
    auto d_ref = end_intersection.value()->delaunay_ref;
    if(d_ref != kin_del.getCrossingData().delaunay_edge_intersections[delaunay_edge_id].begin())
    {
      std::ostringstream oss;
      oss << "resolveIntersectionMeshPairIndex: end_intersection is not the first one on its Delaunay edge (voronoi_cell_id=" << voronoi_cell_id << ", event_time=" << event_time << ").";
      KINDS_ERROR(oss.str());
      throw std::runtime_error(oss.str());
    }

    // check if value is valid
    if(end_intersection.value()->prev_segment_mesh_pair_index == static_cast<size_t>(-1))
    {
      // print the other value for debugging
      size_t next_index = end_intersection.value()->next_segment_mesh_pair_index;
      std::ostringstream oss;
      oss << "resolveIntersectionMeshPairIndex: end_intersection has invalid prev_segment_mesh_pair_index (voronoi_cell_id=" << voronoi_cell_id
          << ", event_time=" << event_time << ", next_segment_mesh_pair_index=" << next_index << ").";
      KINDS_ERROR(oss.str());
      throw std::runtime_error(oss.str());
    }

    return end_intersection.value()->prev_segment_mesh_pair_index;
  }
}

size_t SegmentBuilder::startNewMeshFromIntersections(size_t voronoi_cell_id, double t,
  std::optional<KineticDelaunay::CrossingData::EdgeIntersectionRef> start_intersection,
  std::optional<KineticDelaunay::CrossingData::EdgeIntersectionRef> end_intersection, bool reuse_existing_pair_and_mesh,
  BoundaryEventType event_type, BoundarySegmentAction segment_action, bool force_single_seed_vertex,
  const RadiusBoundaryTransitionShiftContext* boundary_transition_shift)
{
  if (!start_intersection.has_value() && !end_intersection.has_value())
  {
    throw std::runtime_error("startNewMeshFromIntersections requires at least one intersection reference.");
  }

  const auto& graph = kin_del.getGraph();

  // Interval semantics (start/end of the strip) and mesh topology are passed explicitly; only `SegmentBuilder`
  // state is captured via `[this]`.
  const auto finish_for_delaunay_edge = [this, boundary_transition_shift](const HalfEdgeDelaunayGraph& graph,
    size_t delaunay_edge_id, size_t voronoi_cell_id, double t,
    const std::optional<KineticDelaunay::CrossingData::EdgeIntersectionRef>& interval_start_crossing,
    const std::optional<KineticDelaunay::CrossingData::EdgeIntersectionRef>& interval_end_crossing,
    bool reuse_existing_pair_and_mesh, BoundaryEventType event_type, BoundarySegmentAction segment_action,
    bool force_single_seed_vertex) -> size_t {
    // Directed half-edges for this undirected Delaunay edge: even index is one orientation, odd is the twin.
    const size_t he_even = 2 * delaunay_edge_id;
    const size_t he_odd = he_even + 1;
    if (he_odd >= graph.getHalfEdges().size())
    {
      return static_cast<size_t>(-1);
    }

    // `MeshingData` stores the *inside* boundary half-edge for each endpoint that is an actual crossing (not an open site).
    // On a boundary Delaunay edge, one directed half-edge points outward from the meshed component; we store the opposite.
    int inside_boundary_he_id = -1;
    if (kin_del.isOnComponentBoundary(he_even))
    {
      const bool even_is_outside = kin_del.isOnComponentBoundaryOutside(he_even);
      inside_boundary_he_id = static_cast<int>(even_is_outside ? he_odd : he_even);
    }

    size_t intersection_pair_index = static_cast<size_t>(-1);
    bool created_new_pair = false;
    if (reuse_existing_pair_and_mesh)
    {
      // Pair lookup keys off the same interval tuple the caller passed (cell, start/end refs, time).
      intersection_pair_index
        = resolveIntersectionMeshPairIndex(voronoi_cell_id, interval_start_crossing, interval_end_crossing, t);
      if (intersection_pair_index != static_cast<size_t>(-1) && intersection_pair_index >= intersection_segment_mesh_pairs.size())
      {
        intersection_pair_index = static_cast<size_t>(-1);
      }
      if (intersection_pair_index == static_cast<size_t>(-1))
      {
        return static_cast<size_t>(-1);
      }
    }

    const auto& he = graph.getHalfEdges()[he_even];
    const auto& twin_he = graph.getHalfEdges()[he_odd];
    const size_t owner_segment_id = (voronoi_cell_id < strand_to_segment_indices.size() && !strand_to_segment_indices[voronoi_cell_id].empty())
      ? strand_to_segment_indices[voronoi_cell_id].back()
      : static_cast<size_t>(-1);
    MeshStructure::SegmentMeshPair segment_mesh_pair;
    segment_mesh_pair.segment_index0 = he.origin == -1 ? -1 : strand_to_segment_indices[he.origin].back();
    segment_mesh_pair.segment_index1 = twin_he.origin == -1 ? -1 : strand_to_segment_indices[twin_he.origin].back();
    if (intersection_pair_index == static_cast<size_t>(-1))
    {
      created_new_pair = true;
      intersection_pair_index = intersection_segment_mesh_pairs.size();
      intersection_segment_mesh_pairs.push_back(segment_mesh_pair);
    }
    else
    {
      intersection_segment_mesh_pairs[intersection_pair_index] = segment_mesh_pair;
    }

    if (intersection_pair_index < intersection_mesh_pair_metadata.size())
    {
      intersection_mesh_pair_metadata[intersection_pair_index]
        = MeshStructure::IntersectionMeshPairMetadata { voronoi_cell_id, owner_segment_id,
            interval_start_crossing.has_value() ? interval_start_crossing.value()->delaunay_edge_id : static_cast<size_t>(-1),
            interval_end_crossing.has_value() ? interval_end_crossing.value()->delaunay_edge_id : static_cast<size_t>(-1) };
    }
    else
    {
      intersection_mesh_pair_metadata.resize(intersection_pair_index + 1);
      intersection_mesh_pair_metadata[intersection_pair_index]
        = MeshStructure::IntersectionMeshPairMetadata { voronoi_cell_id, owner_segment_id,
            interval_start_crossing.has_value() ? interval_start_crossing.value()->delaunay_edge_id : static_cast<size_t>(-1),
            interval_end_crossing.has_value() ? interval_end_crossing.value()->delaunay_edge_id : static_cast<size_t>(-1) };
    }

    const size_t component_vertex = std::max(he.origin, twin_he.origin);
    const size_t component_id = kin_del.component_data.component_map[component_vertex];
    std::vector<bool> he_visited(graph.getHalfEdges().size(), false);
    updateBoundary(t, he_visited, component_id);
    auto& boundary_polygon = kin_del.component_data.component_boundaries[component_id][0];
    auto& centroid = kin_del.component_data.component_centroids[component_id];

    const bool reuse_in_place = !created_new_pair && intersection_pair_index < intersection_meshes.size();
    VoronoiMesh mesh_local;
    VoronoiMesh& mesh = reuse_in_place ? intersection_meshes[intersection_pair_index] : mesh_local;

    auto intersection_or_cell_position = [&](std::optional<KineticDelaunay::CrossingData::EdgeIntersectionRef> ref)
    {
      if (ref.has_value())
      {
        const KineticDelaunay::CrossingData::EdgeIntersectionRef orig_ref = ref.value();
        const glm::dvec3 old_pos
          = closingMeshVoronoiDelaunayCrossingPosition(t, orig_ref->voronoi_edge_id, orig_ref->delaunay_edge_id);
        const glm::dvec3 new_pos = crossingPositionWithRadiusBoundaryTransitionShift(t, orig_ref, boundary_transition_shift);
        if (boundary_transition_shift != nullptr && new_pos != old_pos)
        {
          KineticDelaunay::CrossingData::EdgeIntersectionRef use_ref = orig_ref;
          if (auto neighbor_opt = neighborIntersectionOnTargetAlongVoronoiEdge(
                orig_ref, orig_ref->voronoi_edge_id, boundary_transition_shift);
            neighbor_opt.has_value())
          {
            use_ref = neighbor_opt.value();
            logRadiusBoundaryTransitionVertexShift(
              "startNewMeshFromIntersections_interval", t, orig_ref, use_ref, old_pos, new_pos);
          }
        }
        return new_pos;
      }
      if (auto site_shifted
        = radiusTransitionInterpolatedSitePosition(t, voronoi_cell_id, delaunay_edge_id, boundary_transition_shift);
        site_shifted.has_value())
      {
        const glm::dvec2 p_site = kin_del.getPointAt(t, voronoi_cell_id);
        KINDS_DEBUG("Radius boundary transition [startNewMeshFromIntersections_site]: cell=" << voronoi_cell_id
                                                                                              << " strip_d=" << delaunay_edge_id
                                                                                              << " t=" << t << " old=("
                                                                                              << p_site.x << "," << p_site.y
                                                                                              << ") new=(" << site_shifted->x
                                                                                              << "," << site_shifted->y << ")");
        return site_shifted.value();
      }
      const glm::dvec2 p = kin_del.getPointAt(t, voronoi_cell_id);
      return glm::dvec3(p.x, p.y, t);
    };

    const glm::dvec3 start_pos = intersection_or_cell_position(interval_start_crossing);
    const glm::dvec3 end_pos = intersection_or_cell_position(interval_end_crossing);
    // Open interval ends are tagged with the owning Voronoi cell id so vertex export can distinguish site vertices from crossings.
    const std::optional<size_t> start_voronoi_vertex
      = interval_start_crossing.has_value() ? std::nullopt : std::optional<size_t>(voronoi_cell_id);
    const std::optional<size_t> end_voronoi_vertex
      = interval_end_crossing.has_value() ? std::nullopt : std::optional<size_t>(voronoi_cell_id);

    std::ostringstream boundary_start_meta_stream;
    boundary_start_meta_stream << "{\"event_type\":\"" << boundaryEventTypeToString(event_type) << "\",\"segment_action\":\""
                               << boundarySegmentActionToString(segment_action) << "\",\"voronoi_cell_id\":" << voronoi_cell_id
                               << ",\"time\":" << t << ",\"start_ref\":";
    if (interval_start_crossing.has_value())
    {
      boundary_start_meta_stream << "{\"delaunay_edge_id\":" << interval_start_crossing.value()->delaunay_edge_id
                                 << ",\"voronoi_edge_id\":" << interval_start_crossing.value()->voronoi_edge_id
                                 << ",\"delaunay_param\":" << interval_start_crossing.value()->delaunay_edge_param << "}";
    }
    else
    {
      boundary_start_meta_stream << "null";
    }
    boundary_start_meta_stream << ",\"end_ref\":";
    if (interval_end_crossing.has_value())
    {
      boundary_start_meta_stream << "{\"delaunay_edge_id\":" << interval_end_crossing.value()->delaunay_edge_id
                                 << ",\"voronoi_edge_id\":" << interval_end_crossing.value()->voronoi_edge_id
                                 << ",\"delaunay_param\":" << interval_end_crossing.value()->delaunay_edge_param << "}";
    }
    else
    {
      boundary_start_meta_stream << "null";
    }
    boundary_start_meta_stream << "}";
    const std::string boundary_start_meta = boundary_start_meta_stream.str();
    // "left"/"right" label the two mesh strip endpoints in interval order (start of interval vs end), not world L/R.
    const std::string boundary_start_meta_left = appendBoundaryVertexPosMetadata(boundary_start_meta, "left");
    const std::string boundary_start_meta_right = appendBoundaryVertexPosMetadata(boundary_start_meta, "right");
    const std::string boundary_start_meta_uniform = appendBoundaryVertexPosMetadata(boundary_start_meta, "uniform");
    const double dx = start_pos.x - end_pos.x;
    const double dy = start_pos.y - end_pos.y;
    const double dz = start_pos.z - end_pos.z;
    const bool same_endpoint = force_single_seed_vertex || (dx * dx + dy * dy + dz * dz) <= 1e-20;

    // First vertex = interval start (left tag); second = interval end (right tag). Uniform when the strip collapses to one vertex.
    const std::string& start_vertex_meta = same_endpoint ? boundary_start_meta_uniform : boundary_start_meta_left;
    const size_t start_vertex_index
      = addMeshletVertex(mesh, boundary_polygon, centroid, start_pos, voronoi_cell_id, t, start_voronoi_vertex, start_vertex_meta,
        same_endpoint ? std::optional<glm::dvec3>(glm::dvec3(1.0, 0.0, 1.0))
                      : std::optional<glm::dvec3>(glm::dvec3(1.0, 0.0, 0.0)));
    const size_t end_vertex_index = same_endpoint
      ? start_vertex_index
      : addMeshletVertex(mesh, boundary_polygon, centroid, end_pos, voronoi_cell_id, t, end_voronoi_vertex,
        boundary_start_meta_right, glm::dvec3(0.0, 0.0, 1.0));

    std::list<MeshingData> local_segments;
    MeshingData seg { static_cast<int>(start_vertex_index), static_cast<int>(end_vertex_index),
      interval_start_crossing.has_value() ? inside_boundary_he_id : -1,
      interval_end_crossing.has_value() ? inside_boundary_he_id : -1 };
    seg.start_crossing = interval_start_crossing;
    seg.end_crossing = interval_end_crossing;
    if (reuse_in_place)
    {
      if (intersection_pair_index >= intersection_mesh_pair_last_left_and_right_vertex.size())
      {
        intersection_mesh_pair_last_left_and_right_vertex.resize(intersection_pair_index + 1);
      }
      auto& segments_for_pair = intersection_mesh_pair_last_left_and_right_vertex[intersection_pair_index];
      segments_for_pair.clear();
      segments_for_pair.emplace_back(std::move(seg));
    }
    else
    {
      local_segments.emplace_back(std::move(seg));
    }

    // Prev/next on `CrossingData` follow list order on the Delaunay edge for `[ref,ref]`; one-null cases use `voronoi_cell_id`
    // vs even-half-edge origin inside `writeIntersectionPairLinks` (see there). Arguments mirror interval start/end crossings.
    writeIntersectionPairLinks(intersection_pair_index, voronoi_cell_id, interval_start_crossing, interval_end_crossing);

    if (reuse_in_place)
    {
      mesh.setCreationKineticTime(t);
    }
    else
    {
      intersection_meshes.push_back(std::move(mesh_local));
      intersection_meshlet_export_suffixes.push_back(std::string("_intersection_d") + std::to_string(delaunay_edge_id));
      intersection_mesh_pair_last_left_and_right_vertex.emplace_back(std::move(local_segments));
      intersection_meshes.back().setCreationKineticTime(t);
    }

    return intersection_pair_index;
  };

  // Both crossings `[ref, ref]` — must agree on the same Delaunay edge.
  if (start_intersection.has_value() && end_intersection.has_value())
  {
    if (start_intersection.value()->delaunay_edge_id != end_intersection.value()->delaunay_edge_id)
    {
      return static_cast<size_t>(-1);
    }
    return finish_for_delaunay_edge(graph, start_intersection.value()->delaunay_edge_id, voronoi_cell_id, t, start_intersection,
      end_intersection, reuse_existing_pair_and_mesh, event_type, segment_action, force_single_seed_vertex);
  }

  // Crossing at start, open site at end `[ref, null]`.
  if (start_intersection.has_value())
  {
    return finish_for_delaunay_edge(graph, start_intersection.value()->delaunay_edge_id, voronoi_cell_id, t, start_intersection,
      end_intersection, reuse_existing_pair_and_mesh, event_type, segment_action, force_single_seed_vertex);
  }

  // Open site at start, crossing at end `[null, ref]`.
  return finish_for_delaunay_edge(graph, end_intersection.value()->delaunay_edge_id, voronoi_cell_id, t, start_intersection,
    end_intersection, reuse_existing_pair_and_mesh, event_type, segment_action, force_single_seed_vertex);
}

void SegmentBuilder::finishMeshFromIntersections(size_t voronoi_cell_id, double t,
  std::optional<KineticDelaunay::CrossingData::EdgeIntersectionRef> start_intersection,
  std::optional<KineticDelaunay::CrossingData::EdgeIntersectionRef> end_intersection, BoundaryEventType event_type,
  BoundarySegmentAction segment_action, const RadiusBoundaryTransitionShiftContext* boundary_transition_shift)
{
  if (!start_intersection.has_value() && !end_intersection.has_value())
  {
    throw std::runtime_error("finishMeshFromIntersections requires at least one intersection reference.");
  }

  const size_t intersection_pair_index
    = resolveIntersectionMeshPairIndex(voronoi_cell_id, start_intersection, end_intersection, t);

  if (intersection_pair_index == static_cast<size_t>(-1) || intersection_pair_index >= intersection_meshes.size()
    || intersection_pair_index >= intersection_mesh_pair_last_left_and_right_vertex.size())
  {
    std::ostringstream oss;
    oss << "finishMeshFromIntersections: failed to resolve mesh pair for extension (cell=" << voronoi_cell_id << ", t=" << t
        << ", start_ref=";
    if (start_intersection.has_value())
    {
      oss << "{d_edge=" << start_intersection.value()->delaunay_edge_id << ", v_edge=" << start_intersection.value()->voronoi_edge_id
          << ", param=" << start_intersection.value()->delaunay_edge_param << "}";
    }
    else
    {
      oss << "null";
    }
    oss << ", end_ref=";
    if (end_intersection.has_value())
    {
      oss << "{d_edge=" << end_intersection.value()->delaunay_edge_id << ", v_edge=" << end_intersection.value()->voronoi_edge_id
          << ", param=" << end_intersection.value()->delaunay_edge_param << "}";
    }
    else
    {
      oss << "null";
    }
    oss << ").";
    KINDS_ERROR(oss.str());
    //throw std::runtime_error(oss.str());
    return;
  }

  auto& segs = intersection_mesh_pair_last_left_and_right_vertex[intersection_pair_index];
  if (segs.empty())
  {
    std::ostringstream oss;
    oss << "finishMeshFromIntersections: missing segment state for pair " << intersection_pair_index << " (cell=" << voronoi_cell_id
        << ", t=" << t << ").";
    throw std::runtime_error(oss.str());
  }
  auto& seg = segs.front();
  if (seg.mesh_start_vertex_id < 0 || seg.mesh_end_vertex_id < 0)
  {
    std::ostringstream oss;
    oss << "finishMeshFromIntersections: missing last-left/last-right indices for pair " << intersection_pair_index
        << " (mesh_start_vertex_id=" << seg.mesh_start_vertex_id << ", mesh_end_vertex_id=" << seg.mesh_end_vertex_id
        << ", cell=" << voronoi_cell_id << ", t=" << t << ").";
    throw std::runtime_error(oss.str());
  }

  auto& mesh = intersection_meshes[intersection_pair_index];

  size_t component_id = kin_del.component_data.component_map[voronoi_cell_id];
  std::vector<bool> he_visited(kin_del.getGraph().getHalfEdges().size(), false);
  updateBoundary(t, he_visited, component_id);
  auto& boundary_polygon = kin_del.component_data.component_boundaries[component_id][0];
  auto centroid = kin_del.component_data.component_centroids[component_id];

  const size_t strip_delaunay_edge_id_for_site
    = start_intersection.has_value() ? start_intersection.value()->delaunay_edge_id : end_intersection.value()->delaunay_edge_id;

  auto endpoint_position_at_t = [&](bool at_start) -> glm::dvec3
  {
    const auto& input_ref = at_start ? start_intersection : end_intersection;
    if (!input_ref.has_value())
    {
      if (auto site_shifted = radiusTransitionInterpolatedSitePosition(
            t, voronoi_cell_id, strip_delaunay_edge_id_for_site, boundary_transition_shift);
        site_shifted.has_value())
      {
        const glm::dvec2 p_site = kin_del.getPointAt(t, voronoi_cell_id);
        KINDS_DEBUG("Radius boundary transition [finishMeshFromIntersections_site]: cell=" << voronoi_cell_id
                                                                                           << " strip_d="
                                                                                           << strip_delaunay_edge_id_for_site
                                                                                           << " t=" << t << " old=(" << p_site.x
                                                                                           << "," << p_site.y << ") new=("
                                                                                           << site_shifted->x << ","
                                                                                           << site_shifted->y << ")");
        return site_shifted.value();
      }
      const glm::dvec2 p = kin_del.getPointAt(t, voronoi_cell_id);
      return glm::dvec3(p.x, p.y, t);
    }

    const KineticDelaunay::CrossingData::EdgeIntersectionRef orig_ref = input_ref.value();
    const glm::dvec3 old_pos
      = closingMeshVoronoiDelaunayCrossingPosition(t, orig_ref->voronoi_edge_id, orig_ref->delaunay_edge_id);
    const glm::dvec3 new_pos = crossingPositionWithRadiusBoundaryTransitionShift(t, orig_ref, boundary_transition_shift);
    if (boundary_transition_shift != nullptr && new_pos != old_pos)
    {
      KineticDelaunay::CrossingData::EdgeIntersectionRef use_ref = orig_ref;
      if (auto neighbor_opt
        = neighborIntersectionOnTargetAlongVoronoiEdge(orig_ref, orig_ref->voronoi_edge_id, boundary_transition_shift);
        neighbor_opt.has_value())
      {
        use_ref = neighbor_opt.value();
        logRadiusBoundaryTransitionVertexShift("finishMeshFromIntersections_interval", t, orig_ref, use_ref, old_pos, new_pos);
      }
    }
    return new_pos;
  };

  const glm::dvec3 new_start_pos = endpoint_position_at_t(true);
  const glm::dvec3 new_end_pos = endpoint_position_at_t(false);
  const std::optional<size_t> start_vv = !start_intersection.has_value() ? std::optional<size_t>(voronoi_cell_id) : std::nullopt;
  const std::optional<size_t> end_vv = !end_intersection.has_value() ? std::optional<size_t>(voronoi_cell_id) : std::nullopt;
  std::ostringstream boundary_finish_meta_stream;
  boundary_finish_meta_stream << "{\"event_type\":\"" << boundaryEventTypeToString(event_type) << "\",\"segment_action\":\""
                              << boundarySegmentActionToString(segment_action) << "\",\"voronoi_cell_id\":" << voronoi_cell_id
                              << ",\"time\":" << t << ",\"start_ref\":";
  if (start_intersection.has_value())
  {
    boundary_finish_meta_stream << "{\"delaunay_edge_id\":" << start_intersection.value()->delaunay_edge_id
                                << ",\"voronoi_edge_id\":" << start_intersection.value()->voronoi_edge_id
                                << ",\"delaunay_param\":" << start_intersection.value()->delaunay_edge_param << "}";
  }
  else
  {
    boundary_finish_meta_stream << "null";
  }
  boundary_finish_meta_stream << ",\"end_ref\":";
  if (end_intersection.has_value())
  {
    boundary_finish_meta_stream << "{\"delaunay_edge_id\":" << end_intersection.value()->delaunay_edge_id
                                << ",\"voronoi_edge_id\":" << end_intersection.value()->voronoi_edge_id
                                << ",\"delaunay_param\":" << end_intersection.value()->delaunay_edge_param << "}";
  }
  else
  {
    boundary_finish_meta_stream << "null";
  }
  boundary_finish_meta_stream << "}";
  const std::string boundary_finish_meta = boundary_finish_meta_stream.str();
  const std::string boundary_finish_meta_left = appendBoundaryVertexPosMetadata(boundary_finish_meta, "left");
  const std::string boundary_finish_meta_right = appendBoundaryVertexPosMetadata(boundary_finish_meta, "right");
  const size_t new_start_vertex_index
    = addMeshletVertex(mesh, boundary_polygon, centroid, new_start_pos, voronoi_cell_id, t, start_vv, boundary_finish_meta_left,
      glm::dvec3(1.0, 0.0, 0.0));
  const size_t new_end_vertex_index
    = addMeshletVertex(mesh, boundary_polygon, centroid, new_end_pos, voronoi_cell_id, t, end_vv, boundary_finish_meta_right,
      glm::dvec3(0.0, 0.0, 1.0));

  size_t ordered_new_start_vertex_index = new_start_vertex_index;
  size_t ordered_new_end_vertex_index = new_end_vertex_index;
  const int old_fixed_start_id = seg.mesh_start_vertex_id;
  const int old_fixed_end_id = seg.mesh_end_vertex_id;
  const size_t last_left = static_cast<size_t>(old_fixed_start_id);
  const size_t last_right = static_cast<size_t>(old_fixed_end_id);
  const size_t eff_left = intersectionStripEffectiveVertexIndex(seg, true);
  const size_t eff_right = intersectionStripEffectiveVertexIndex(seg, false);
  if (last_left < mesh.getVertices().size() && last_right < mesh.getVertices().size())
  {
    const glm::dvec3& prev_start = mesh.getVertices()[last_left];
    const glm::dvec3& prev_end = mesh.getVertices()[last_right];
    const glm::dvec3& curr_start = mesh.getVertices()[ordered_new_start_vertex_index];
    const glm::dvec3& curr_end = mesh.getVertices()[ordered_new_end_vertex_index];
    const auto sqd = [](const glm::dvec3& a, const glm::dvec3& b) {
      const glm::dvec3 d = a - b;
      return d.x * d.x + d.y * d.y + d.z * d.z;
    };
    const double direct_cost = sqd(prev_start, curr_start) + sqd(prev_end, curr_end);
    const double swapped_cost = sqd(prev_start, curr_end) + sqd(prev_end, curr_start);
    if (swapped_cost + 1e-18 < direct_cost)
    {
      KINDS_WARNING("finishMeshFromIntersections: detected swapped endpoint continuity for pair "
                    << intersection_pair_index << " (cell=" << voronoi_cell_id << ", t=" << t
                    << "), auto-correcting start/end assignment.");
      std::swap(ordered_new_start_vertex_index, ordered_new_end_vertex_index);
    }
  }

  // Single new anchor (collapsed interval / one mesh vertex for both ends): both flex chains meet that point.
  const bool uniform_finish_targets = (ordered_new_start_vertex_index == ordered_new_end_vertex_index);
  const size_t flex_interp_target = ordered_new_start_vertex_index;
  interpolateFlexibleVerticesAlongEdge(mesh, seg.flexible_left_vertex_ids, static_cast<size_t>(old_fixed_start_id),
    uniform_finish_targets ? flex_interp_target : ordered_new_start_vertex_index);
  interpolateFlexibleVerticesAlongEdge(mesh, seg.flexible_right_vertex_ids, static_cast<size_t>(old_fixed_end_id),
    uniform_finish_targets ? flex_interp_target : ordered_new_end_vertex_index);
  seg.flexible_left_vertex_ids.clear();
  seg.flexible_right_vertex_ids.clear();

  const int inside_boundary_he_id = (seg.start_half_edge_id >= 0) ? seg.start_half_edge_id : seg.end_half_edge_id;
  const auto& verts_after = mesh.getVertices();
  bool start_apex_first = false;
  if (eff_left != eff_right)
  {
    const double zl = verts_after[eff_left].z;
    const double zr = verts_after[eff_right].z;
    constexpr double z_eps = 1e-12;
    if (std::abs(zl - zr) > z_eps)
    {
      start_apex_first = zl < zr;
    }
    else
    {
      // Flexible placeholders often share z with strip corners; pick the diagonal from projected order along eff_l–eff_r.
      const glm::dvec2 pl(verts_after[eff_left].x, verts_after[eff_left].y);
      const glm::dvec2 pr(verts_after[eff_right].x, verts_after[eff_right].y);
      const glm::dvec2 ps(verts_after[ordered_new_start_vertex_index].x, verts_after[ordered_new_start_vertex_index].y);
      const glm::dvec2 pe(verts_after[ordered_new_end_vertex_index].x, verts_after[ordered_new_end_vertex_index].y);
      const glm::dvec2 d = pr - pl;
      const double strip_len2 = glm::dot(d, d);
      if (strip_len2 > 1e-24)
      {
        start_apex_first = glm::dot(ps - pl, d) < glm::dot(pe - pl, d);
      }
      else
      {
        start_apex_first = zl < zr;
      }
    }
  }
  if (eff_left == eff_right)
  {
    addBoundaryIntervalTriangleOriented(
      mesh, ordered_new_start_vertex_index, eff_right, ordered_new_end_vertex_index, inside_boundary_he_id, t,
      boundary_finish_meta);
  }
  else if (start_apex_first)
  {
    addBoundaryIntervalTriangleOriented(
      mesh, eff_left, eff_right, ordered_new_start_vertex_index, inside_boundary_he_id, t, boundary_finish_meta);
    addBoundaryIntervalTriangleOriented(
      mesh, ordered_new_start_vertex_index, eff_right, ordered_new_end_vertex_index, inside_boundary_he_id, t,
      boundary_finish_meta);
  }
  else
  {
    addBoundaryIntervalTriangleOriented(
      mesh, eff_left, eff_right, ordered_new_end_vertex_index, inside_boundary_he_id, t, boundary_finish_meta);
    addBoundaryIntervalTriangleOriented(
      mesh, eff_left, ordered_new_end_vertex_index, ordered_new_start_vertex_index, inside_boundary_he_id, t,
      boundary_finish_meta);
  }

  seg.mesh_start_vertex_id = static_cast<int>(ordered_new_start_vertex_index);
  seg.mesh_end_vertex_id = static_cast<int>(ordered_new_end_vertex_index);
}

std::vector<KineticDelaunay::CrossingData::EdgeIntersectionRef> SegmentBuilder::getBoundaryIntersectionsInBoundaryOrder(
  size_t delaunay_edge_id) const
{
  const auto& crossing_data = kin_del.getCrossingData();
  if (delaunay_edge_id >= crossing_data.delaunay_edge_intersections.size())
  {
    std::ostringstream oss;
    oss << "getBoundaryIntersectionsInBoundaryOrder: Delaunay edge out of bounds (edge=" << delaunay_edge_id
        << ", crossing_data_size=" << crossing_data.delaunay_edge_intersections.size() << ").";
    throw std::runtime_error(oss.str());
  }

  std::vector<KineticDelaunay::CrossingData::EdgeIntersectionRef> refs;
  const auto& d_intersections = crossing_data.delaunay_edge_intersections[delaunay_edge_id];
  for (const auto& ref : d_intersections)
  {
    refs.push_back(ref);
  }
  return refs;
}

bool SegmentBuilder::writeOneNullIntersectionPairLinkByNullVertex(
  size_t intersection_pair_index, size_t null_vertex_id, KineticDelaunay::CrossingData::EdgeIntersectionRef ref,
  bool interval_is_ref_to_null)
{
  const size_t d_edge_id = ref->delaunay_edge_id;
  const size_t he_even = 2 * d_edge_id;
  const size_t he_odd = he_even + 1;
  const auto& graph = kin_del.getGraph();
  if (he_odd >= graph.getHalfEdges().size())
  {
    throw std::runtime_error("writeOneNullIntersectionPairLinkByNullVertex: Delaunay half-edge out of bounds.");
  }

  const int even_origin = graph.getHalfEdges()[he_even].origin;
  const bool write_prev = even_origin >= 0 && null_vertex_id == static_cast<size_t>(even_origin);
  if (write_prev)
  {
    ref->prev_segment_mesh_pair_index = intersection_pair_index;
  }
  else
  {
    ref->next_segment_mesh_pair_index = intersection_pair_index;
  }
  KINDS_DEBUG("one-null link write " << (interval_is_ref_to_null ? "[ref,null]" : "[null,ref]") << ": de=" << d_edge_id
                                     << " pair=" << intersection_pair_index << " null_vertex=" << null_vertex_id
                                     << " even_origin=" << even_origin << " wrote=" << (write_prev ? "prev" : "next"));
  return write_prev;
}

void SegmentBuilder::writeIntersectionPairLinks(size_t intersection_pair_index, size_t voronoi_cell_id,
  std::optional<KineticDelaunay::CrossingData::EdgeIntersectionRef> start_intersection,
  std::optional<KineticDelaunay::CrossingData::EdgeIntersectionRef> end_intersection)
{
  const auto edge_intersection_ref_desc
    = [](const std::optional<KineticDelaunay::CrossingData::EdgeIntersectionRef>& ref) -> std::string {
      if (!ref.has_value())
      {
        return "null";
      }
      return "{d_edge=" + std::to_string(ref.value()->delaunay_edge_id) + ", v_edge="
        + std::to_string(ref.value()->voronoi_edge_id) + ", param=" + std::to_string(ref.value()->delaunay_edge_param) + "}";
    };

  std::string wrote_pair_index_to;

  // Keep prev/next links consistent with the actual order in delaunay_edge_intersections:
  // if y follows x in list order -> x.next = pair, y.prev = pair.
  if (start_intersection.has_value() && end_intersection.has_value())
  {
    const size_t d_edge_id = start_intersection.value()->delaunay_edge_id;
    if (d_edge_id < kin_del.getCrossingData().delaunay_edge_intersections.size())
    {
      const auto& d_list = kin_del.getCrossingData().delaunay_edge_intersections[d_edge_id];
      bool seen_start = false;
      bool seen_end = false;
      for (auto it = d_list.begin(); it != d_list.end(); ++it)
      {
        if (*it == start_intersection.value())
        {
          seen_start = true;
          break;
        }
        if (*it == end_intersection.value())
        {
          seen_end = true;
          break;
        }
      }
      if (seen_start)
      {
        start_intersection.value()->next_segment_mesh_pair_index = intersection_pair_index;
        end_intersection.value()->prev_segment_mesh_pair_index = intersection_pair_index;
        wrote_pair_index_to = "start_ref->next_segment_mesh_pair_index, end_ref->prev_segment_mesh_pair_index";
      }
      else if (seen_end)
      {
        start_intersection.value()->prev_segment_mesh_pair_index = intersection_pair_index;
        end_intersection.value()->next_segment_mesh_pair_index = intersection_pair_index;
        wrote_pair_index_to = "start_ref->prev_segment_mesh_pair_index, end_ref->next_segment_mesh_pair_index";
      }
      else
      {
        start_intersection.value()->next_segment_mesh_pair_index = intersection_pair_index;
        end_intersection.value()->prev_segment_mesh_pair_index = intersection_pair_index;
        wrote_pair_index_to
          = "start_ref->next_segment_mesh_pair_index, end_ref->prev_segment_mesh_pair_index (list lookup fallback)";
      }
    }
    else
    {
      start_intersection.value()->next_segment_mesh_pair_index = intersection_pair_index;
      end_intersection.value()->prev_segment_mesh_pair_index = intersection_pair_index;
      wrote_pair_index_to
        = "start_ref->next_segment_mesh_pair_index, end_ref->prev_segment_mesh_pair_index (missing crossing list fallback)";
    }
  }
  else
  {
    std::ostringstream w;
    if (start_intersection.has_value())
    {
      const bool prev = writeOneNullIntersectionPairLinkByNullVertex(
        intersection_pair_index, voronoi_cell_id, start_intersection.value(), true);
      w << "start_ref->" << (prev ? "prev_segment_mesh_pair_index" : "next_segment_mesh_pair_index");
    }
    if (end_intersection.has_value())
    {
      const bool prev = writeOneNullIntersectionPairLinkByNullVertex(
        intersection_pair_index, voronoi_cell_id, end_intersection.value(), false);
      if (start_intersection.has_value())
      {
        w << ", ";
      }
      w << "end_ref->" << (prev ? "prev_segment_mesh_pair_index" : "next_segment_mesh_pair_index");
    }
    wrote_pair_index_to = w.str();
  }

  std::ostringstream log;
  log << "writeIntersectionPairLinks: pair=" << intersection_pair_index << " cell=" << voronoi_cell_id << " start_ref="
      << edge_intersection_ref_desc(start_intersection) << " end_ref=" << edge_intersection_ref_desc(end_intersection)
      << " wrote_pair_index_to=" << wrote_pair_index_to;
  KINDS_DEBUG(log.str());
}

size_t SegmentBuilder::determineVoronoiCellForBoundaryIntersectionInterval(size_t delaunay_edge_id,
  std::optional<KineticDelaunay::CrossingData::EdgeIntersectionRef> start_intersection,
  std::optional<KineticDelaunay::CrossingData::EdgeIntersectionRef> end_intersection) const
{
  if (!start_intersection.has_value() && !end_intersection.has_value())
  {
    throw std::runtime_error(
      "determineVoronoiCellForBoundaryIntersectionInterval requires at least one intersection reference.");
  }

  const auto& graph = kin_del.getGraph();
  const size_t he_even = 2 * delaunay_edge_id;
  const size_t he_odd = he_even + 1;
  if (he_odd >= graph.getHalfEdges().size())
  {
    std::ostringstream oss;
    oss << "determineVoronoiCellForBoundaryIntersectionInterval: Delaunay edge out of bounds (edge=" << delaunay_edge_id
        << ", half_edge_count=" << graph.getHalfEdges().size() << ").";
    throw std::runtime_error(oss.str());
  }

  if (!start_intersection.has_value() && end_intersection.has_value())
  {
    return static_cast<size_t>(graph.getHalfEdges()[he_even].origin);
  }
  if (start_intersection.has_value() && !end_intersection.has_value())
  {
    return static_cast<size_t>(graph.getHalfEdges()[he_odd].origin);
  }

  const size_t ve0 = start_intersection.value()->voronoi_edge_id;
  const size_t ve1 = end_intersection.value()->voronoi_edge_id;
  const size_t ve0_he_even = 2 * ve0;
  const size_t ve0_he_odd = ve0_he_even + 1;
  const size_t ve1_he_even = 2 * ve1;
  const size_t ve1_he_odd = ve1_he_even + 1;
  if (ve0_he_odd >= graph.getHalfEdges().size() || ve1_he_odd >= graph.getHalfEdges().size())
  {
    std::ostringstream oss;
    oss << "determineVoronoiCellForBoundaryIntersectionInterval: Voronoi edge out of bounds (ve0=" << ve0
        << ", ve1=" << ve1 << ", half_edge_count=" << graph.getHalfEdges().size() << ").";
    throw std::runtime_error(oss.str());
  }

  const size_t a0 = static_cast<size_t>(graph.getHalfEdges()[ve0_he_even].origin);
  const size_t b0 = static_cast<size_t>(graph.getHalfEdges()[ve0_he_odd].origin);
  const size_t a1 = static_cast<size_t>(graph.getHalfEdges()[ve1_he_even].origin);
  const size_t b1 = static_cast<size_t>(graph.getHalfEdges()[ve1_he_odd].origin);
  if (a0 == a1 || a0 == b1)
  {
    return a0;
  }
  if (b0 == a1 || b0 == b1)
  {
    return b0;
  }

  std::ostringstream oss;
  oss << "determineVoronoiCellForBoundaryIntersectionInterval: no shared Voronoi cell between intersection Voronoi edges "
      << ve0 << " and " << ve1 << " on Delaunay edge " << delaunay_edge_id << ".";
  throw std::runtime_error(oss.str());
}

size_t kinDS::SegmentBuilder::registerMeshletWithSuffix(VoronoiMesh&& mesh, std::string suffix, double creation_kinetic_time)
{
  if (std::isfinite(creation_kinetic_time))
  {
    mesh.setCreationKineticTime(creation_kinetic_time);
  }
  const size_t index = meshes.size();
  meshes.push_back(std::move(mesh));
  meshlet_export_suffixes.push_back(std::move(suffix));
  return index;
}


void kinDS::SegmentBuilder::completeBoundaryMeshSection(size_t he_id, size_t new_left, size_t new_right)
{
  const auto& last_left_and_right = boundary_mesh_last_left_and_right_vertex[he_id];
  if (last_left_and_right.first != -1)
  {
    // distinguish the case that we have previously flipped an infinite edge that became a boundary edge
    if (half_edge_to_boundary_vertex_index[he_id] == -1)
    {
      // no edge flip
      addBoundaryTriangle(last_left_and_right.first, new_right, new_left);
      if (last_left_and_right.second != -1)
      {
        addBoundaryTriangle(last_left_and_right.second, new_right, last_left_and_right.first);
      }
    }
    else
    {
      assert(last_left_and_right.second != -1);
      // he_id was previously flipped, it's corresponding vertex is no longer part of the boundary
      addBoundaryTriangle(last_left_and_right.second, new_right, half_edge_to_boundary_vertex_index[he_id]);
      addBoundaryTriangle(new_left, last_left_and_right.first, half_edge_to_boundary_vertex_index[he_id]);
      addBoundaryTriangle(new_left, half_edge_to_boundary_vertex_index[he_id], new_right);

      // reset the half-edge to boundary vertex index
      half_edge_to_boundary_vertex_index[he_id] = -1;
    }
  }
  else
  {
    assert(last_left_and_right.second == -1);
  }
}

size_t kinDS::SegmentBuilder::addBoundaryTriangle(size_t u, size_t v, size_t w)
{
  // check bounds
  if (u >= boundary_mesh.getVertexCount() || v >= boundary_mesh.getVertexCount() || w >= boundary_mesh.getVertexCount())
  {
    KINDS_ERROR("Vertex index out of boundary mesh range.");
    return -1;
  }

  if (u >= boundary_mesh_raw_uvs.size() || v >= boundary_mesh_raw_uvs.size() || w >= boundary_mesh_raw_uvs.size())
  {
    KINDS_ERROR("Vertex index out of raw uv range.");
    return -1;
  }

  // get raw UVs
  glm::dvec3 uv_u = { boundary_mesh_raw_uvs[u][0], boundary_mesh_raw_uvs[u][1], 0.0 };
  glm::dvec3 uv_v = { boundary_mesh_raw_uvs[v][0], boundary_mesh_raw_uvs[v][1], 0.0 };
  glm::dvec3 uv_w = { boundary_mesh_raw_uvs[w][0], boundary_mesh_raw_uvs[w][1], 0.0 };

  // output UVs
  /*KINDS_DEBUG("Adding boundary triangle with raw UVs: u(" + std::to_string(uv_u[0]) + ", " + std::to_string(uv_u[1])
     +
                "), v(" + std::to_string(uv_v[0]) + ", " + std::to_string(uv_v[1]) + "), w(" + std::to_string(uv_w[0]) +
                ", " + std::to_string(uv_w[1]) + ")");*/

  // adjust UVs to avoid seams, first coordinate is the angle normalized to [-0.5, 0.5]
  // As a heuristic, we take the first angle and adjust the others such that they have less than 0.5 difference
  double base_angle = uv_u[0];
  double& angle_v = uv_v[0];
  double diff_v = angle_v - base_angle;
  double adjustment = std::round(diff_v);
  angle_v -= adjustment;

  double& angle_w = uv_w[0];
  double diff_w = angle_w - base_angle;
  adjustment = std::round(diff_w);
  angle_w -= adjustment;

  uv_u[0] *= uv_circum_factor;
  uv_v[0] *= uv_circum_factor;
  uv_w[0] *= uv_circum_factor;
  uv_u[1] *= uv_height_factor;
  uv_v[1] *= uv_height_factor;
  uv_w[1] *= uv_height_factor;

  // add adjusted UVs
  size_t uv_index_u = boundary_mesh.addUV(uv_u);
  size_t uv_index_v = boundary_mesh.addUV(uv_v);
  size_t uv_index_w = boundary_mesh.addUV(uv_w);

  /*KINDS_DEBUG("UVs after adjustment: u(" + std::to_string(uv_u[0]) + ", " + std::to_string(uv_u[1]) + "), v(" +
                std::to_string(uv_v[0]) + ", " + std::to_string(uv_v[1]) + "), w(" + std::to_string(uv_w[0]) + ", " +
                std::to_string(uv_w[1]) + ")");*/
  return boundary_mesh.addTriangle(u, v, w, uv_index_u, uv_index_v, uv_index_w, 0);
}

size_t kinDS::SegmentBuilder::addBoundaryVertex(glm::dvec3 vertex, glm::dvec2 centroid, size_t strand_id, double t)
{
  double angle = std::atan2(centroid[1] - vertex[1], centroid[0] - vertex[0]);

  glm::dvec2 raw_uv { angle / (2.0 * glm::pi<double>()), vertex[2] };
  if (create_transformed_mesh)
  {
    vertex = kin_del.getStrandTree().transformToObjectSpace(vertex, strand_id, t);
  }

  size_t index = boundary_mesh.addVertex(vertex);
  boundary_vertex_to_strand_id.push_back(strand_id);
  boundary_mesh_raw_uvs.resize(index + 1, glm::dvec2 {});
  boundary_mesh_raw_uvs[index] = raw_uv;

  return index;
}

size_t kinDS::SegmentBuilder::addMeshletTriangle(
  VoronoiMesh& mesh, size_t u, size_t v, size_t w, const std::string& metadata, int material_id)
{
  if (mesh.getMaterialNames().empty())
  {
    mesh.setMaterialNames({ MeshletExportMaterialNames[static_cast<size_t>(RegularMeshletMaterialId)] });
  }
  return mesh.addTriangle(u, v, w, u, v, w, material_id, metadata);
}

size_t kinDS::SegmentBuilder::addBoundaryIntervalTriangleOriented(
  VoronoiMesh& mesh, size_t u, size_t v, size_t w, int inside_boundary_he_id, double t, const std::string& metadata)
{
  (void)t;
  if (mesh.getMaterialNames().empty())
  {
    mesh.setMaterialNames({ MeshletExportMaterialNames[static_cast<size_t>(BoundaryIntervalMeshletMaterialId)] });
  }
  const int boundary_material_id = 0;
  if (inside_boundary_he_id < 0 || static_cast<size_t>(inside_boundary_he_id) >= kin_del.getGraph().getHalfEdges().size())
  {
    return addMeshletTriangle(mesh, u, v, w, metadata, boundary_material_id);
  }

  // `inside_boundary_he_id` is the inside-directed boundary half-edge; its twin is the outside one on the same Delaunay edge.
  const size_t outside_he = static_cast<size_t>(inside_boundary_he_id) ^ 1u;
  if (outside_he >= kin_del.getGraph().getHalfEdges().size())
  {
    return addMeshletTriangle(mesh, u, v, w, metadata, boundary_material_id);
  }

  if ((outside_he & 1u) != 0u)
  {
    std::swap(v, w);
  }
  return addMeshletTriangle(mesh, u, v, w, metadata, boundary_material_id);
}

void kinDS::SegmentBuilder::warnIfVoronoiVertexOutsideAlphaShape(
  const char* context, size_t voronoi_vertex_id, const glm::dvec3& position) const
{
  const size_t containing_tri_id = kin_del.getCrossingDataContainingTriId(voronoi_vertex_id);
  if (kin_del.getFaceInside(containing_tri_id))
  {
    return;
  }
  KINDS_WARNING("SegmentBuilder: " << context << " - Voronoi vertex " << voronoi_vertex_id << " (containing Delaunay triangle "
    << containing_tri_id << ") is outside the alpha-shape; position (" << position.x << ", " << position.y << ", "
    << position.z << ").");
}

size_t kinDS::SegmentBuilder::addMeshletVertex(VoronoiMesh& mesh, const std::vector<BoundaryPoint>& boundary_polygon,
  const glm::dvec2& centroid, glm::dvec3 vertex, size_t strand_id, double t,
  std::optional<size_t> meshlet_voronoi_vertex_for_alpha_check, const std::string& metadata,
  const std::optional<glm::dvec3>& debug_color)
{
  const auto warn_degenerate_or_non_finite = [&](const glm::dvec3& p, const char* stage)
  {
    if (!std::isfinite(p.x) || !std::isfinite(p.y) || !std::isfinite(p.z))
    {
      KINDS_WARNING("addMeshletVertex(" << stage << "): non-finite vertex for strand " << strand_id << " at t=" << t
                                        << " -> (" << p.x << ", " << p.y << ", " << p.z << ").");
      return;
    }
    if (p.x == 0.0 && p.y == 0.0)
    {
      KINDS_WARNING("addMeshletVertex(" << stage << "): degenerate XY vertex (0,0,z) for strand " << strand_id
                                        << " at t=" << t << " -> (" << p.x << ", " << p.y << ", " << p.z << ").");
    }
  };

  warn_degenerate_or_non_finite(vertex, "input");
  if (create_transformed_mesh)
  {
    vertex = kin_del.getStrandTree().transformToObjectSpace(vertex, strand_id, t);
  }
  warn_degenerate_or_non_finite(vertex, "stored");
  if (meshlet_voronoi_vertex_for_alpha_check.has_value())
  {
    warnIfVoronoiVertexOutsideAlphaShape("addMeshletVertex", meshlet_voronoi_vertex_for_alpha_check.value(), vertex);
  }
  size_t index = mesh.addVertex(vertex, metadata, debug_color.has_value() ? debug_color.value() : glm::dvec3(1.0));
  double rel_dist = relativeDistanceFromCenter(boundary_polygon, centroid, glm::dvec2 { vertex[0], vertex[1] });

  /*if (rel_dist > 1.0 + std::numeric_limits<double>::epsilon()) {
    KINDS_WARNING("Adding vertex that is too far outside, relative distance: " << rel_dist);
  }*/

  // TODO: this can be simplified to not use trigonometric functions
  double angle = std::atan2(centroid[1] - vertex[1], centroid[0] - vertex[0]);
  double u = 0.5 + texture_diameter * rel_dist * 0.5 * std::cos(angle);
  double v = 0.5 + texture_diameter * rel_dist * 0.5 * std::sin(angle);
  mesh.addUV(u, v, vertex[2] * uv_height_factor);
  return index;
}

size_t SegmentBuilder::intersectionStripEffectiveVertexIndex(const MeshingData& seg, bool left_side) const
{
  if (left_side)
  {
    if (!seg.flexible_left_vertex_ids.empty())
    {
      const int v = seg.flexible_left_vertex_ids.back();
      if (v >= 0)
      {
        return static_cast<size_t>(v);
      }
    }
    if (seg.mesh_start_vertex_id < 0)
    {
      throw std::runtime_error("intersectionStripEffectiveVertexIndex: invalid mesh_start_vertex_id.");
    }
    return static_cast<size_t>(seg.mesh_start_vertex_id);
  }
  if (!seg.flexible_right_vertex_ids.empty())
  {
    const int v = seg.flexible_right_vertex_ids.back();
    if (v >= 0)
    {
      return static_cast<size_t>(v);
    }
  }
  if (seg.mesh_end_vertex_id < 0)
  {
    throw std::runtime_error("intersectionStripEffectiveVertexIndex: invalid mesh_end_vertex_id.");
  }
  return static_cast<size_t>(seg.mesh_end_vertex_id);
}

void SegmentBuilder::addFlexibleVertexToIntersectionMesh(VoronoiMesh& mesh, MeshingData& seg, bool flexible_on_left_side,
  const std::vector<BoundaryPoint>& boundary_polygon, const glm::dvec2& centroid, size_t strand_id, double t,
  const std::string& metadata)
{
  if (!intersection_strip_flexible_vertices_enabled)
  {
    return;
  }
  if (seg.mesh_start_vertex_id < 0 || seg.mesh_end_vertex_id < 0)
  {
    return;
  }
  // Wedge base must match crossing/finish: latest flex on a side wins over fixed mesh_start/mesh_end. Snapshot before
  // appending the new placeholder so the base is (eff_l, eff_r), not always (mesh_start, mesh_end).
  const size_t eff_l = intersectionStripEffectiveVertexIndex(seg, true);
  const size_t eff_r = intersectionStripEffectiveVertexIndex(seg, false);
  if (eff_l == eff_r)
  {
    return;
  }

  // XY must not be degenerate on the strip line or (0,0): addBoundaryIntervalTriangleOriented uses mesh positions for winding.
  const glm::dvec3 placeholder(centroid.x, centroid.y, t);
  const std::string vertex_meta = appendBoundaryVertexPosMetadata(metadata, flexible_on_left_side ? "left" : "right");
  const size_t idx = addMeshletVertex(mesh, boundary_polygon, centroid, placeholder, strand_id, t, std::nullopt, vertex_meta,
    std::nullopt);
  if (flexible_on_left_side)
  {
    seg.flexible_left_vertex_ids.push_back(static_cast<int>(idx));
  }
  else
  {
    seg.flexible_right_vertex_ids.push_back(static_cast<int>(idx));
  }

  const int inside_he = seg.start_half_edge_id >= 0 ? seg.start_half_edge_id : seg.end_half_edge_id;
  if (inside_he < 0)
  {
    return;
  }
  addBoundaryIntervalTriangleOriented(mesh, eff_l, eff_r, idx, inside_he, t, vertex_meta);
}

void SegmentBuilder::applyIntersectionStripOneSidedFixedVertex(VoronoiMesh& mesh, MeshingData& seg, bool fixed_start_side,
  size_t new_fixed_vertex_index, int inside_half_edge_id,
  const std::optional<KineticDelaunay::CrossingData::EdgeIntersectionRef>& new_crossing_for_updated_side,
  const std::vector<BoundaryPoint>& boundary_polygon, const glm::dvec2& centroid, size_t strand_id, double t,
  bool keep_strip_alive)
{
  const int old_fixed = fixed_start_side ? seg.mesh_start_vertex_id : seg.mesh_end_vertex_id;
  if (old_fixed < 0)
  {
    return;
  }
  std::vector<int>& flex = fixed_start_side ? seg.flexible_left_vertex_ids : seg.flexible_right_vertex_ids;
  interpolateFlexibleVerticesAlongEdge(mesh, flex, static_cast<size_t>(old_fixed), new_fixed_vertex_index);
  flex.clear();
  if (!keep_strip_alive)
  {
    return;
  }
  if (fixed_start_side)
  {
    seg.mesh_start_vertex_id = static_cast<int>(new_fixed_vertex_index);
    seg.start_half_edge_id = inside_half_edge_id;
    seg.start_crossing = new_crossing_for_updated_side;
  }
  else
  {
    seg.mesh_end_vertex_id = static_cast<int>(new_fixed_vertex_index);
    seg.end_half_edge_id = inside_half_edge_id;
    seg.end_crossing = new_crossing_for_updated_side;
  }
  addFlexibleVertexToIntersectionMesh(mesh, seg, !fixed_start_side, boundary_polygon, centroid, strand_id, t,
    "{\"intersection_flexible_placeholder\":true}");
}

void SegmentBuilder::applyIntersectionStripUniformClosureVertex(VoronoiMesh& mesh, MeshingData& seg, size_t closure_vertex_index)
{
  if (seg.mesh_start_vertex_id < 0 || seg.mesh_end_vertex_id < 0)
  {
    return;
  }
  interpolateFlexibleVerticesAlongEdge(
    mesh, seg.flexible_left_vertex_ids, static_cast<size_t>(seg.mesh_start_vertex_id), closure_vertex_index);
  interpolateFlexibleVerticesAlongEdge(
    mesh, seg.flexible_right_vertex_ids, static_cast<size_t>(seg.mesh_end_vertex_id), closure_vertex_index);
  seg.flexible_left_vertex_ids.clear();
  seg.flexible_right_vertex_ids.clear();
}

void kinDS::SegmentBuilder::addVoronoiTriangulationToBoundaryMesh(double t, bool invert_orientation, double offset)
{
  auto& graph = kin_del.getGraph();

  size_t index_offset = boundary_mesh.getVertices().size();
  std::vector<double> relative_center_distances;
  // add all vertices
  for (size_t i = 0; i < graph.getVertexCount(); i++)
  {
    glm::dvec2 vertex = kin_del.getPointAt(t, i);

    auto component_index = kin_del.component_data.component_map[i];
    auto& boundary_polygon = kin_del.component_data.component_boundaries[component_index][0];
    auto& centroid = kin_del.component_data.component_centroids[component_index];

    addBoundaryVertex(glm::dvec3 { vertex[0], vertex[1], t + offset }, centroid, i, t);
    // KINDS_DEBUG("New raw uv: " << raw_uv[0] << ", " << raw_uv[1] << " for vertex: " << vertex_index);

    double rel_dist = relativeDistanceFromCenter(boundary_polygon, centroid, vertex);
    relative_center_distances.push_back(rel_dist);
  }
  // add all triangles
  // for (const auto& triangle : graph.getFaces()) {
  for (size_t face_index = 0; face_index < graph.getFaces().size(); face_index++)
  {
    const auto& triangle = graph.getFaces()[face_index];
    const auto& he_ids = triangle.half_edges;
    auto vertices = graph.adjacentTriangleVertices(triangle.half_edges[0]);

    // check if on component boundary and update last left and right vertices accordingly
    // store last left and right
    for (size_t i = 0; i < 3; i++)
    {
      if (kin_del.isOnComponentBoundaryOutside(he_ids[i]))
      {
        completeBoundaryMeshSection(he_ids[i], index_offset + vertices[i], index_offset + vertices[(i + 1) % 3]);
        // KINDS_DEBUG("Assigning boundary last vertices for he_id " << he_ids[i]);
        boundary_mesh_last_left_and_right_vertex[he_ids[i]]
          = std::make_pair(index_offset + vertices[i], index_offset + vertices[(i + 1) % 3]);
      }
    }
    // skip faces that are outside
    if (!kin_del.getFaceInside(face_index))
    {
      continue;
    }

    // check for infinite vertices
    if (vertices[0] == -1 || vertices[1] == -1 || vertices[2] == -1)
    {
      continue; // skip triangles with infinite vertices
    }

    if (invert_orientation)
    {
      std::swap(vertices[1], vertices[2]);
    }

    // next, get the angles from the raw UVs and convert back to cartesian coordinates centered at (0.5, 0.5)
    size_t uv_indices[3];

    for (size_t i = 0; i < 3; i++)
    {
      double rel_dist = relative_center_distances[vertices[i]];
      double angle = boundary_mesh_raw_uvs[index_offset + vertices[i]][0] * 2.0 * glm::pi<double>();
      double u = 0.5 + texture_diameter * rel_dist * 0.5 * std::cos(angle);
      double v = 0.5 + texture_diameter * rel_dist * 0.5 * std::sin(angle);
      uv_indices[i] = boundary_mesh.addUV(u, v, 0.0);
    }

    // as an exception, we directly add the triangle here to have access to the UV indices
    boundary_mesh.addTriangle(index_offset + vertices[0], index_offset + vertices[1], index_offset + vertices[2],
      uv_indices[0], uv_indices[1], uv_indices[2], 1);
  }
}

std::vector<BoundaryPoint> kinDS::SegmentBuilder::traceConvexHull(double t) const
{
  const auto& graph = kin_del.getGraph();
  std::vector<BoundaryPoint> convex_hull_points;
  for (HalfEdgeDelaunayGraph::ConvexHullEdgeIterator it = graph.boundaryEdgesBegin(), end = graph.boundaryEdgesEnd();
    it != end; ++it)
  {
    size_t he_id = *it;
    size_t strand_index = graph.getHalfEdges()[he_id].origin;

    glm::dvec2 convex_hull_point = kin_del.getPointAt(t, strand_index);
    convex_hull_points.push_back({ strand_index, he_id, convex_hull_point });
  }

  return convex_hull_points;
}

void kinDS::SegmentBuilder::advanceBoundaryMesh(
  double t, const std::vector<BoundaryPoint>& boundary_points, const glm::dvec2& centroid)
{
  std::vector<size_t> new_vertex_indices;

  for (size_t i = 0; i < boundary_points.size(); i++)
  {
    glm::dvec2 boundary_point = boundary_points[i].p;

    new_vertex_indices.push_back(addBoundaryVertex(
      glm::dvec3 { boundary_point[0], boundary_point[1], t }, centroid, boundary_points[i].vertex_id, t));
  }

  for (size_t i = 0; i < boundary_points.size(); i++)
  {
    size_t he_id = boundary_points[i].he_id;
    size_t left_vertex_index = new_vertex_indices[i];
    size_t right_vertex_index = new_vertex_indices[(i + 1) % boundary_points.size()];
    auto& left_and_right = boundary_mesh_last_left_and_right_vertex[he_id];
    completeBoundaryMeshSection(he_id, left_vertex_index, right_vertex_index);

    // KINDS_DEBUG("Assigning boundary last vertices for he_id " << he_id);
    left_and_right.first = left_vertex_index;
    left_and_right.second = right_vertex_index;
  }
}

void kinDS::SegmentBuilder::updateBoundary(double t, std::vector<bool>& visited, size_t component_index)
{
  if (kin_del.component_data.component_last_updated[component_index] != t)
  {
    kin_del.component_data.component_boundaries[component_index]
      = kin_del.extractComponentBoundaries(kin_del.component_data.components[component_index], t, visited);
    kin_del.component_data.component_centroids[component_index]
      = polygonCentroid(kin_del.component_data.component_boundaries[component_index][0]);
    kin_del.component_data.component_last_updated[component_index] = t;
  }
}

void kinDS::SegmentBuilder::updateBoundaries(double t)
{
  std::vector<bool> visited(kin_del.getGraph().getHalfEdges().size(), false);

  for (size_t component_index = 0; component_index < kin_del.component_data.components.size(); component_index++)
  {
    updateBoundary(t, visited, component_index);
  }
}

void kinDS::SegmentBuilder::advanceBoundaryMeshes(double t)
{
  updateBoundaries(t);

  for (size_t component_index = 0; component_index < kin_del.component_data.components.size(); component_index++)
  {
    auto& boundaries = kin_del.component_data.component_boundaries[component_index];
    auto& centroid = kin_del.component_data.component_centroids[component_index];

    for (auto& boundary_points : boundaries)
    {
      advanceBoundaryMesh(t, boundary_points, centroid);
    }
  }
}

std::vector<glm::dvec3> SegmentBuilder::computeClampedVoronoiVertices(
  size_t strand_id, double t, const std::vector<BoundaryPoint>& boundary_polygon, const glm::dvec2& centroid)
{
  auto& graph = kin_del.getGraph();
  std::vector<glm::dvec3> voronoi_vertices;

  // we just create a triangle fan because the Voronoi cell is convex
  // iterate over all segment indices of the strand
  for (HalfEdgeDelaunayGraph::IncidentEdgeIterator it = graph.incidentEdgesBegin(strand_id),
                                                   end = graph.incidentEdgesEnd(strand_id);
    it != end; ++it)
  {
    glm::dvec3 voronoi_vertex = computeVoronoiVertex(*it, t);
    voronoi_vertices.emplace_back(voronoi_vertex);
    // addMeshletVertex(mesh, boundary_polygon, centroid, voronoi_vertex);
  }

  std::vector<glm::dvec3> clamped_voronoi_vertices;
  // go through the vertices and clamp them
  for (size_t i = 0; i < voronoi_vertices.size(); i++)
  {
    auto left = voronoi_vertices[i];
    auto right = voronoi_vertices[(i + 1) % voronoi_vertices.size()];

    bool inside = clampVoronoiVertices(left, right, boundary_polygon, centroid);

    if (inside)
    {
      clamped_voronoi_vertices.push_back(left);
      clamped_voronoi_vertices.push_back(right);
    }
  }

  voronoi_vertices.clear();
  // Put back into original vector and remove duplicates
  voronoi_vertices.emplace_back(clamped_voronoi_vertices.front());
  for (size_t i = 1; i < clamped_voronoi_vertices.size() - 1; i++)
  {
    if (clamped_voronoi_vertices[i] != voronoi_vertices.back())
    {
      voronoi_vertices.push_back(clamped_voronoi_vertices[i]);
    }
  }
  if (clamped_voronoi_vertices.front() != clamped_voronoi_vertices.back())
  {
    voronoi_vertices.push_back(clamped_voronoi_vertices.back());
  }

  return voronoi_vertices;
}

size_t kinDS::SegmentBuilder::closingMeshCountStrandIncidentEdges(size_t strand_id) const
{
  auto& graph = kin_del.getGraph();
  size_t num = 0;
  for (HalfEdgeDelaunayGraph::IncidentEdgeIterator it = graph.incidentEdgesBegin(strand_id),
                                                   end = graph.incidentEdgesEnd(strand_id);
    it != end; ++it)
  {
    ++num;
  }
  return num;
}

int kinDS::SegmentBuilder::closingMeshAppendVertex(VoronoiMesh& mesh,
  const std::vector<BoundaryPoint>& boundary_polygon, const glm::dvec2& centroid, size_t strand_id, double t,
  const glm::dvec3& position, std::optional<size_t> voronoi_vertex_for_alpha_check)
{
  const size_t vertex_id
    = addMeshletVertex(mesh, boundary_polygon, centroid, position, strand_id, t, voronoi_vertex_for_alpha_check);
  return static_cast<int>(vertex_id);
}

std::optional<KineticDelaunay::CrossingData::EdgeIntersectionRef> kinDS::SegmentBuilder::closingMeshFindVoronoiEdgeIntersection(
  size_t voronoi_edge_id, size_t crossed_delaunay_he_id) const
{
  const size_t d = crossed_delaunay_he_id / 2;
  const auto& crossing_data = kin_del.getCrossingData();
  if (voronoi_edge_id >= crossing_data.voronoi_edge_intersections.size())
  {
    return std::nullopt;
  }
  for (const auto& ref : crossing_data.voronoi_edge_intersections[voronoi_edge_id])
  {
    if (ref->delaunay_edge_id == d)
    {
      return ref;
    }
  }
  return std::nullopt;
}

void kinDS::SegmentBuilder::refreshMeshingDataCrossingRefs(MeshingData& seg, size_t voronoi_edge_id) const
{
  if (seg.start_half_edge_id >= 0)
  {
    seg.start_crossing
      = closingMeshFindVoronoiEdgeIntersection(voronoi_edge_id, static_cast<size_t>(seg.start_half_edge_id));
  }
  else
  {
    seg.start_crossing.reset();
  }
  if (seg.end_half_edge_id >= 0)
  {
    seg.end_crossing
      = closingMeshFindVoronoiEdgeIntersection(voronoi_edge_id, static_cast<size_t>(seg.end_half_edge_id));
  }
  else
  {
    seg.end_crossing.reset();
  }
}

void kinDS::SegmentBuilder::refreshStripCrossingRefsIncidentToVoronoiVertex(size_t voronoi_vertex_id)
{
  auto& graph = kin_del.getGraph();
  const auto& face = graph.getFaces()[voronoi_vertex_id];
  std::unordered_set<size_t> refreshed_pairs;
  for (size_t vhe : face.half_edges)
  {
    const size_t even_he = vhe & ~1;
    if (even_he >= half_edge_index_to_segment_mesh_pair_index.size())
    {
      continue;
    }
    const size_t pair_idx = half_edge_index_to_segment_mesh_pair_index[even_he];
    if (pair_idx >= segment_mesh_pair_last_left_and_right_vertex.size())
    {
      continue;
    }
    if (!refreshed_pairs.insert(pair_idx).second)
    {
      continue;
    }
    const size_t strip_ve = even_he / 2;
    for (auto& seg : segment_mesh_pair_last_left_and_right_vertex[pair_idx])
    {
      refreshMeshingDataCrossingRefs(seg, strip_ve);
    }
  }
}

void kinDS::SegmentBuilder::refreshCrossingRefsForAllStrips()
{
  for (size_t pair_idx = 0; pair_idx < segment_mesh_pair_last_left_and_right_vertex.size(); ++pair_idx)
  {
    size_t voronoi_edge_id = 0;
    bool found = false;
    for (size_t he = 0; he < half_edge_index_to_segment_mesh_pair_index.size(); ++he)
    {
      if (half_edge_index_to_segment_mesh_pair_index[he] == pair_idx)
      {
        voronoi_edge_id = (he & ~1) / 2;
        found = true;
        break;
      }
    }
    if (!found)
    {
      continue;
    }
    for (auto& seg : segment_mesh_pair_last_left_and_right_vertex[pair_idx])
    {
      refreshMeshingDataCrossingRefs(seg, voronoi_edge_id);
    }
  }
}

glm::dvec3 kinDS::SegmentBuilder::closingMeshVoronoiDelaunayCrossingPosition(
  double t, size_t voronoi_edge_id, size_t delaunay_edge_id) const
{
  auto& graph = kin_del.getGraph();
  const size_t voronoi_he0 = 2 * voronoi_edge_id;
  const size_t voronoi_he1 = voronoi_he0 + 1;

  const glm::dvec3 left_vertex = computeVoronoiVertex(voronoi_he0, t);
  const glm::dvec3 right_vertex = computeVoronoiVertex(voronoi_he1, t);
  const glm::dvec2 left2(left_vertex.x, left_vertex.y);
  const glm::dvec2 right2(right_vertex.x, right_vertex.y);

  const size_t d_he0 = 2 * delaunay_edge_id;
  const size_t d_he1 = d_he0 + 1;
  const size_t d_a = graph.getHalfEdges()[d_he0].origin;
  const size_t d_b = graph.getHalfEdges()[d_he1].origin;
  if (d_a == size_t(-1) || d_b == size_t(-1))
  {
    KINDS_WARNING("closingMeshVoronoiDelaunayCrossingPosition: degenerate Delaunay edge for voronoi_edge=" << voronoi_edge_id << " delaunay_edge=" << delaunay_edge_id << " t=" << t);
    return glm::dvec3(std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN(), t);
  }

  const glm::dvec2 d0 = kin_del.getPointAt(t, d_a);
  const glm::dvec2 d1 = kin_del.getPointAt(t, d_b);
  const glm::dvec2 p = intersectSegments(left2, right2, d0, d1);
  return glm::dvec3(p, t);
}

auto kinDS::SegmentBuilder::extractSegmentsForVoronoiEdge(double t, int incident_edge_index, size_t voronoi_edge_id,
  const std::function<int(const glm::dvec3&, std::optional<size_t>)>& track_vertex, bool reverse) -> std::vector<MeshingData>
{
  std::vector<MeshingData> out_segments;
  auto& graph = kin_del.getGraph();

  const size_t voronoi_he_even = 2 * voronoi_edge_id;
  const size_t voronoi_he_odd = voronoi_he_even + 1;

  if (graph.isInfinite(voronoi_he_even))
  {
    return out_segments;
  }

  // `reverse`: when false, walk even circumcenter → odd; when true, odd → even (strand-side-first in closingMesh).
  const size_t voronoi_he_strand_side = reverse ? voronoi_he_odd : voronoi_he_even;
  const size_t voronoi_he_other_side = voronoi_he_strand_side ^ 1;
  const bool strand_at_voronoi_even_he = !reverse;

  const glm::dvec3 strand_cm_vertex = computeVoronoiVertex(voronoi_he_strand_side, t);
  const glm::dvec3 other_cm_vertex = computeVoronoiVertex(voronoi_he_other_side, t);
  const glm::dvec2 strand2 = glm::dvec2(strand_cm_vertex);
  const glm::dvec2 other2 = glm::dvec2(other_cm_vertex);
  const glm::dvec2 vor_dir = other2 - strand2;
  const double vor_len2 = glm::dot(vor_dir, vor_dir);
  if (vor_len2 < 1e-14)
  {
    return out_segments;
  }

  const size_t strand_face_id = graph.getHalfEdges()[voronoi_he_strand_side].face;
  const size_t strand_containing_tri_id = kin_del.getCrossingDataContainingTriId(strand_face_id);
  bool inside = kin_del.getFaceInside(strand_containing_tri_id);
  if (inside)
  {
    out_segments.push_back(MeshingData { track_vertex(strand_cm_vertex, graph.getHalfEdges()[voronoi_he_strand_side].face),
      -1, -1, -1 });
    out_segments.back().closing_incident_edge_index = incident_edge_index;
    out_segments.back().closing_voronoi_edge_id = voronoi_edge_id;
    out_segments.back().closing_strand_at_voronoi_even_he = strand_at_voronoi_even_he;
  }

  if (!kin_del.computeBoundaryOnTheFly())
  {
    if (inside && !out_segments.empty())
    {
      out_segments.back().mesh_end_vertex_id
        = track_vertex(other_cm_vertex, graph.getHalfEdges()[voronoi_he_other_side].face);
    }
    for (auto& seg : out_segments)
    {
      refreshMeshingDataCrossingRefs(seg, voronoi_edge_id);
    }
    return out_segments;
  }

  const auto crossed_edges_params = kin_del.computeCrossedHalfEdges(strand_containing_tri_id, other2, strand2, t);
  for (size_t crossed_he_id : crossed_edges_params.first)
  {
    const size_t next_face_id = graph.getHalfEdges()[crossed_he_id ^ 1].face;
    const bool next_inside = kin_del.getFaceInside(next_face_id);

    if (inside != next_inside)
    {
      const size_t da = graph.getHalfEdges()[crossed_he_id].origin;
      const size_t db = graph.getHalfEdges()[crossed_he_id ^ 1].origin;
      if (da != size_t(-1) && db != size_t(-1))
      {
        const glm::dvec2 d0 = kin_del.getPointAt(t, da);
        const glm::dvec2 d1 = kin_del.getPointAt(t, db);
        const glm::dvec2 p = intersectSegments(strand2, other2, d0, d1);

        if (std::isfinite(p.x) && std::isfinite(p.y))
        {
          const int mesh_vertex_id = track_vertex(glm::dvec3(p, t), std::nullopt);
          const auto intersection_it = closingMeshFindVoronoiEdgeIntersection(voronoi_edge_id, crossed_he_id);

          if (inside)
          {
            if (!out_segments.empty())
            {
              auto& segment = out_segments.back();
              segment.mesh_end_vertex_id = mesh_vertex_id;
              segment.end_half_edge_id = static_cast<int>(crossed_he_id);
              segment.end_crossing = intersection_it;
            }
          }
          else
          {
            MeshingData segment { mesh_vertex_id, -1, static_cast<int>(crossed_he_id ^ 1), -1 };
            segment.start_crossing = intersection_it;
            segment.closing_incident_edge_index = incident_edge_index;
            segment.closing_voronoi_edge_id = voronoi_edge_id;
            segment.closing_strand_at_voronoi_even_he = strand_at_voronoi_even_he;
            out_segments.push_back(segment);
          }
        }
      }
    }

    inside = next_inside;
  }

  if (inside && !out_segments.empty())
  {
    out_segments.back().mesh_end_vertex_id
      = track_vertex(other_cm_vertex, graph.getHalfEdges()[voronoi_he_other_side].face);
  }

  for (auto& seg : out_segments)
  {
    refreshMeshingDataCrossingRefs(seg, voronoi_edge_id);
  }

  return out_segments;
}

auto kinDS::SegmentBuilder::closingMeshExtractRawSegmentsForVoronoiEdge(size_t strand_id, double t, VoronoiMesh& mesh,
  const std::vector<BoundaryPoint>& boundary_polygon, const glm::dvec2& centroid, int incident_edge_index,
  size_t incident_he, const std::function<int(const glm::dvec3&, std::optional<size_t>)>& track_vertex)
  -> std::vector<MeshingData>
{
  auto& graph = kin_del.getGraph();

  // Undirected Delaunay edge id (matches `CrossingData` / `computeEdgeIntersections`: voronoi_edge_id == he/2).
  const size_t voronoi_edge_id = incident_he / 2;
  const size_t voronoi_he_even = 2 * voronoi_edge_id;
  const size_t voronoi_he_odd = voronoi_he_even + 1;

  if (graph.isInfinite(voronoi_he_even))
  {
    return {};
  }

  const int origin_even = graph.getHalfEdges()[voronoi_he_even].origin;
  const int origin_odd = graph.getHalfEdges()[voronoi_he_odd].origin;
  const auto finite_origin_is_strand = [&](int o) { return o >= 0 && static_cast<size_t>(o) == strand_id; };
  if (!finite_origin_is_strand(origin_even) && !finite_origin_is_strand(origin_odd))
  {
    KINDS_ERROR("Closing mesh: Voronoi edge for dual Delaunay edge_id "
      << voronoi_edge_id << " (half-edges " << voronoi_he_even << ", " << voronoi_he_odd << ", incident_he "
      << incident_he << ") has no finite half-edge with origin strand_id " << strand_id << " (origins " << origin_even
      << ", " << origin_odd << ")");
  }

  // Orient the inside trace from the Voronoi half-edge whose origin is the strand (matches how we leave the strand
  // along the dual). `CrossingData::computeEdgeIntersections` still stores the Voronoi list in even->odd order; when
  // the strand is on the odd half-edge we reverse traversal of `delaunay_edge_intersections` during the boundary
  // walk.
  const bool strand_at_voronoi_even_he = finite_origin_is_strand(origin_even);
  const bool reverse = !strand_at_voronoi_even_he;

  return extractSegmentsForVoronoiEdge(t, incident_edge_index, voronoi_edge_id, track_vertex, reverse);
}

SegmentBuilder::ClosingMeshOrderedIndex kinDS::SegmentBuilder::closingMeshBuildOrderedIndex(
  std::list<MeshingData>& closing_segments)
{
  ClosingMeshOrderedIndex out;
  out.ordered_segments.reserve(closing_segments.size());
  for (auto& seg : closing_segments)
  {
    if (seg.mesh_start_vertex_id == -1 || seg.mesh_end_vertex_id == -1)
    {
      continue;
    }
    out.ordered_segments.push_back(&seg);
  }

  out.start_crossing_to_segment.reserve(out.ordered_segments.size());
  for (size_t i = 0; i < out.ordered_segments.size(); ++i)
  {
    MeshingData* s = out.ordered_segments[i];
    if (s->start_crossing.has_value())
    {
      const auto ref = s->start_crossing.value();
      const auto ins = out.start_crossing_to_segment.insert({ref, i});
      if (!ins.second)
      {
        KINDS_ERROR("Closing mesh: duplicate start_crossing: ordered_seg[" << ins.first->second << "] and ordered_seg["
                                                                             << i << "] share the same start crossing.");
      }
    }
  }
  return out;
}

void kinDS::SegmentBuilder::closingMeshValidateOrderedSegmentGeometry(
  double t, const VoronoiMesh& mesh, const std::vector<MeshingData*>& ordered_segments)
{
  constexpr double k_closing_cap_geom_eps = 1e-5;
  auto mesh_vertex_xy = [&](int mesh_vid) -> glm::dvec2
  {
    const glm::dvec3& v = mesh.getVertices()[static_cast<size_t>(mesh_vid)];
    return glm::dvec2(v.x, v.y);
  };

  for (size_t si = 0; si < ordered_segments.size(); ++si)
  {
    const MeshingData* s = ordered_segments[si];
    if (s->closing_voronoi_edge_id == static_cast<size_t>(-1))
    {
      continue;
    }
    const size_t ve = s->closing_voronoi_edge_id;
    {
      const std::string ctxs = std::string("ordered_seg[") + std::to_string(si)
        + "] incident=" + std::to_string(s->closing_incident_edge_index);
      validateClosingCapCrossingRef(kin_del, (ctxs + " start_ref").c_str(), s->start_crossing, ve, s->start_half_edge_id);
      validateClosingCapCrossingRef(kin_del, (ctxs + " end_ref").c_str(), s->end_crossing, ve, s->end_half_edge_id);
    }

    const glm::dvec3 L3 = kin_del.computeVoronoiVertexClampedInfinity(2 * ve, t);
    const glm::dvec3 R3 = kin_del.computeVoronoiVertexClampedInfinity(2 * ve + 1, t);
    const glm::dvec2 L(L3.x, L3.y);
    const glm::dvec2 R(R3.x, R3.y);
    const glm::dvec2 axis = R - L;
    const double ax2 = glm::dot(axis, axis);
    if (ax2 < 1e-24)
    {
      KINDS_ERROR("Closing mesh ordered_seg[" << si << "]: degenerate Voronoi edge " << ve << " (zero length).");
      continue;
    }
    const glm::dvec2 ps = mesh_vertex_xy(s->mesh_start_vertex_id);
    const glm::dvec2 pe = mesh_vertex_xy(s->mesh_end_vertex_id);
    const double axis_len = glm::length(axis);
    const auto perp_dist
      = [&](const glm::dvec2& p) -> double { return std::abs((p.x - L.x) * axis.y - (p.y - L.y) * axis.x) / axis_len; };
    if (perp_dist(ps) > k_closing_cap_geom_eps || perp_dist(pe) > k_closing_cap_geom_eps)
    {
      KINDS_ERROR("Closing mesh ordered_seg["
        << si << "]: mesh_start/end not collinear with canonical Voronoi "
        << "edge " << ve << " (even->odd circumcenters); check segment direction vs inside walk.");
    }
    const double ts = glm::dot(ps - L, axis) / ax2;
    const double te = glm::dot(pe - L, axis) / ax2;
    // mesh_start is always the strand-side circumcenter; along the canonical even->odd axis that is the smaller t iff
    // the strand Voronoi half-edge is even.
    if (s->closing_strand_at_voronoi_even_he)
    {
      if (ts > te + 1e-8)
      {
        KINDS_ERROR("Closing mesh ordered_seg[" << si << "]: mesh vertices reversed along Voronoi edge " << ve
                                                << " (t_start=" << ts << " > t_end=" << te
                                                << "); expected mesh_start nearer voronoi_he_even.");
      }
    }
    else if (te > ts + 1e-8)
    {
      KINDS_ERROR("Closing mesh ordered_seg[" << si << "]: mesh vertices reversed along Voronoi edge " << ve
                                                << " (t_start=" << ts << " > t_end=" << te
                                                << " along even->odd); expected mesh_start nearer "
                                                << "voronoi_he_odd (strand-side).");
    }

    glm::dvec2 ip;
    if (s->start_crossing.has_value()
      && tryComputeCrossingIntersectionPosition2D(kin_del, s->start_crossing, t, ip))
    {
      if (glm::distance(ip, ps) > k_closing_cap_geom_eps)
      {
        KINDS_ERROR("Closing mesh ordered_seg[" << si << "]: start_ref 2D position does not match mesh_start.");
      }
    }
    if (s->end_crossing.has_value()
      && tryComputeCrossingIntersectionPosition2D(kin_del, s->end_crossing, t, ip))
    {
      if (glm::distance(ip, pe) > k_closing_cap_geom_eps)
      {
        KINDS_ERROR("Closing mesh ordered_seg[" << si << "]: end_ref 2D position does not match mesh_end.");
      }
    }
  }
}

SegmentBuilder::ClosingMeshPolygonsTraceResult kinDS::SegmentBuilder::closingMeshTraceCapPolygons(size_t strand_id,
  double t, size_t num_incident_edges, VoronoiMesh& mesh, const std::vector<BoundaryPoint>& boundary_polygon,
  const glm::dvec2& centroid, std::vector<size_t> mesh_vertex_ids, const std::vector<MeshingData*>& ordered_segments,
  const std::unordered_map<KineticDelaunay::CrossingData::EdgeIntersectionRef, size_t, ClosingMeshCrossingIteratorHash,
    ClosingMeshCrossingIteratorEq>& start_crossing_to_segment)
{
  ClosingMeshPolygonsTraceResult result;
  result.mesh_vertex_ids = std::move(mesh_vertex_ids);
  result.segment_used.assign(ordered_segments.size(), false);

  auto& graph = kin_del.getGraph();
  const auto& crossing_data = kin_del.getCrossingData();

  auto& cap_vertex_ids = result.mesh_vertex_ids;
  auto add_mesh_vertex = [&](const glm::dvec3& v) -> int
  {
    const int id = closingMeshAppendVertex(mesh, boundary_polygon, centroid, strand_id, t, v);
    cap_vertex_ids.push_back(static_cast<size_t>(id));
    return id;
  };

  auto polygon_push = [&](std::vector<size_t>& poly, size_t vertex_id)
  {
    if (poly.empty() || poly.back() != vertex_id)
    {
      poly.push_back(vertex_id);
    }
  };

  // Only connect along a Delaunay boundary edge to intersections whose Voronoi edge is dual to an edge incident to
  // this strand; otherwise the first CrossingData hit on that Delaunay edge can belong to another cell and we
  // incorrectly skip the strand's Voronoi continuation.
  auto voronoi_edge_incident_to_strand = [&](size_t voronoi_edge_id) -> bool
  {
    const size_t he0 = 2 * voronoi_edge_id;
    const size_t he1 = he0 + 1;
    if (he1 >= graph.getHalfEdges().size())
    {
      return false;
    }
    const int o0 = graph.getHalfEdges()[he0].origin;
    const int o1 = graph.getHalfEdges()[he1].origin;
    auto finite_origin_is_strand = [&](int o) { return o >= 0 && static_cast<size_t>(o) == strand_id; };
    return finite_origin_is_strand(o0) || finite_origin_is_strand(o1);
  };

  // Dual graph: Voronoi half-edge `.face` is the incident Voronoi vertex (circumcenter) id — stable across strips;
  // `mesh_*_vertex_id` are mesh buffer indices and must not be used to match circumcenter joins.
  auto closing_strip_strand_voronoi_half_edge = [](const MeshingData& s) -> size_t
  {
    if (s.closing_voronoi_edge_id == static_cast<size_t>(-1))
    {
      return static_cast<size_t>(-1);
    }
    const size_t ve = s.closing_voronoi_edge_id;
    return s.closing_strand_at_voronoi_even_he ? (2 * ve) : (2 * ve + 1);
  };
  auto closing_strip_voronoi_vertex_at_strand_endpoint = [&](const MeshingData& s) -> size_t
  {
    const size_t he_strand = closing_strip_strand_voronoi_half_edge(s);
    if (he_strand == static_cast<size_t>(-1) || he_strand >= graph.getHalfEdges().size())
    {
      return static_cast<size_t>(-1);
    }
    return graph.getHalfEdges()[he_strand].face;
  };
  auto closing_strip_voronoi_vertex_at_other_endpoint = [&](const MeshingData& s) -> size_t
  {
    const size_t he_strand = closing_strip_strand_voronoi_half_edge(s);
    if (he_strand == static_cast<size_t>(-1) || he_strand >= graph.getHalfEdges().size())
    {
      return static_cast<size_t>(-1);
    }
    return graph.getHalfEdges()[he_strand ^ 1].face;
  };

  for (;;)
  {
    size_t seed_index = static_cast<size_t>(-1);
    for (size_t i = 0; i < ordered_segments.size(); ++i)
    {
      if (!result.segment_used[i])
      {
        seed_index = i;
        break;
      }
    }
    if (seed_index == static_cast<size_t>(-1))
    {
      break;
    }

    std::vector<size_t> polygon;
    MeshingData* seed_seg = ordered_segments[seed_index];
    {
      const size_t ve_seed = seed_seg->closing_voronoi_edge_id;
      if (ve_seed == static_cast<size_t>(-1))
      {
        KINDS_DEBUG("Starting polygon walk from ordered segment " << seed_index << " of Voronoi edge ?");
      }
      else
      {
        KINDS_DEBUG("Starting polygon walk from ordered segment " << seed_index << " of Voronoi edge " << ve_seed);
      }
    }

    size_t current_segment_index = seed_index;
    size_t last_walk_segment_index = seed_index;
    size_t guard = 0;

    std::vector<int> cap_start_hits(ordered_segments.size(), 0);
    std::vector<int> cap_end_handed_off(ordered_segments.size(), 0);
    std::vector<size_t> walk_segment_order;
    walk_segment_order.reserve(ordered_segments.size());
    bool exited_fail = false;
    bool exited_closed_loop = false;
    std::optional<KineticDelaunay::CrossingData::EdgeIntersectionRef> walk_closure_crossing_ref;

    auto record_cap_handoff = [&](size_t from_idx, size_t to_idx)
    {
      cap_end_handed_off[from_idx]++;
      cap_start_hits[to_idx]++;
    };

    while (guard++ < ordered_segments.size() * 4 && !exited_fail)
    {
      bool closed_walk_at_seed = false;
      MeshingData* current_segment = ordered_segments[current_segment_index];
      last_walk_segment_index = current_segment_index;
      result.segment_used[current_segment_index] = true;
      walk_segment_order.push_back(current_segment_index);

      polygon_push(polygon, static_cast<size_t>(current_segment->mesh_start_vertex_id));
      logClosingMeshStripEndpointVertex(
        kin_del, current_segment_index, current_segment->start_crossing, current_segment->closing_voronoi_edge_id);
      polygon_push(polygon, static_cast<size_t>(current_segment->mesh_end_vertex_id));
      logClosingMeshStripEndpointVertex(
        kin_del, current_segment_index, current_segment->end_crossing, current_segment->closing_voronoi_edge_id);

      const bool ends_at_boundary_crossing
        = current_segment->end_crossing.has_value() && current_segment->end_half_edge_id >= 0;

      if (!ends_at_boundary_crossing)
      {
        const size_t joint_vv = closing_strip_voronoi_vertex_at_other_endpoint(*current_segment);
        if (joint_vv == static_cast<size_t>(-1))
        {
          KINDS_ERROR("Closing mesh: ordered_seg[" << current_segment_index
            << "] has no boundary end_crossing but closing_voronoi_edge_id is unset; cannot resolve circumcenter Voronoi "
               "vertex for join.");
          exited_fail = true;
          break;
        }
        size_t next_j = static_cast<size_t>(-1);
        for (size_t j = 0; j < ordered_segments.size(); ++j)
        {
          if (j == current_segment_index || result.segment_used[j])
          {
            continue;
          }
          if (closing_strip_voronoi_vertex_at_strand_endpoint(*ordered_segments[j]) == joint_vv)
          {
            next_j = j;
            break;
          }
        }
        if (next_j != static_cast<size_t>(-1))
        {
          record_cap_handoff(current_segment_index, next_j);
          KINDS_DEBUG("Polygon walk: Voronoi-vertex join from ordered_seg[" << current_segment_index << "] to ordered_seg["
            << next_j << "] at Voronoi vertex (circumcenter) id " << joint_vv
            << " (strip end at circumcenter, no boundary end_crossing).");
          current_segment_index = next_j;
          continue;
        }
        const size_t seed_start_vv = closing_strip_voronoi_vertex_at_strand_endpoint(*seed_seg);
        if (joint_vv == seed_start_vv && current_segment_index != seed_index
          && seed_start_vv != static_cast<size_t>(-1))
        {
          walk_closure_crossing_ref.reset();
          record_cap_handoff(current_segment_index, seed_index);
          KINDS_DEBUG("Polygon walk: Voronoi-vertex join closes loop to seed ordered_seg[" << seed_index
            << "] at Voronoi vertex (circumcenter) id " << joint_vv << ".");
          exited_closed_loop = true;
          break;
        }
        KINDS_ERROR("Closing mesh: ordered_seg[" << current_segment_index
          << "] ends at circumcenter Voronoi vertex id " << joint_vv
          << " with no boundary end_crossing, but no unused strip has that id as its strand-side Voronoi vertex (and it "
             "does not match seed's strand-side Voronoi vertex); cannot continue cap walk.");
        exited_fail = true;
        break;
      }

      const int exit_he_id = current_segment->end_half_edge_id;
      std::optional<KineticDelaunay::CrossingData::EdgeIntersectionRef> current_ref_opt = current_segment->end_crossing;

      const size_t start_boundary_he = kin_del.isOnComponentBoundaryOutside(static_cast<size_t>(exit_he_id))
        ? static_cast<size_t>(exit_he_id)
        : static_cast<size_t>(exit_he_id) ^ 1;

      size_t boundary_he = start_boundary_he;
      bool found_next_segment = false;
      size_t boundary_guard = 0;

      constexpr double k_strand_corner_dist_eps = 1e-9;
      auto append_strand_corner_if_needed = [&]()
      {
        const glm::dvec2 xy = kin_del.getPointAt(t, strand_id);
        const glm::dvec3 strand_pos(xy, t);
        if (!polygon.empty())
        {
          const glm::dvec3& last = mesh.getVertices()[polygon.back()];
          if (glm::distance(glm::dvec2(last.x, last.y), xy) <= k_strand_corner_dist_eps)
          {
            return;
          }
        }
        const int vid = add_mesh_vertex(strand_pos);
        polygon_push(polygon, static_cast<size_t>(vid));
        KINDS_DEBUG("Adding vertex at strand " << strand_id);
      };

      while (boundary_guard++ < graph.getHalfEdges().size() * 2)
      {
        const int he_origin = graph.getHalfEdges()[boundary_he].origin;
        if (he_origin >= 0 && static_cast<size_t>(he_origin) == strand_id)
        {
          append_strand_corner_if_needed();
        }

        const auto& d_intersections = crossing_data.delaunay_edge_intersections[boundary_he / 2];
        // `delaunay_edge_intersections[e]` is ordered along the even Delaunay half-edge (2*e).
        // For a boundary edge, select the inside-oriented half-edge; traverse forward iff that inside half-edge is even,
        // otherwise traverse backward.
        const size_t inside_boundary_he
          = kin_del.isOnComponentBoundaryOutside(boundary_he) ? (boundary_he ^ 1) : boundary_he;
        const bool effective_list_forward = ((inside_boundary_he % 2) != 0);

        // `traced_boundary_intervals` record endpoints in canonical interval order (even-half-edge / increasing list
        // index). The walk may list crossings backward when `effective_list_forward` is false — then swap walk order
        // (a_walk -> b_walk) to (start_intersection, end_intersection).
        const auto append_traced_boundary_interval
          = [&](const std::optional<KineticDelaunay::CrossingData::EdgeIntersectionRef>& a_walk,
              const std::optional<KineticDelaunay::CrossingData::EdgeIntersectionRef>& b_walk) {
            if (effective_list_forward)
            {
              result.traced_boundary_intervals.push_back(BoundaryIntersectionInterval { strand_id, a_walk, b_walk });
            }
            else
            {
              result.traced_boundary_intervals.push_back(BoundaryIntersectionInterval { strand_id, b_walk, a_walk });
            }
          };

        // Walk every crossing along this directed Delaunay half-edge after `current_ref_opt`: each crossing gets a
        // mesh vertex. Strand-incident crossings that match an ordered segment start_crossing hand off to that segment;
        // others only advance along the same Delaunay edge (no skipping to the corner).
        bool exit_boundary_chain_to_voronoi = false;
        if (!d_intersections.empty())
        {
          std::optional<KineticDelaunay::CrossingData::EdgeIntersectionRef> last_ref_on_edge = std::nullopt;
          if (current_ref_opt.has_value() && current_ref_opt.value()->delaunay_edge_id == (boundary_he / 2))
          {
            last_ref_on_edge = current_ref_opt;
          }

          auto find_next_intersection = [&](std::optional<KineticDelaunay::CrossingData::EdgeIntersectionRef> ref_opt)
            -> std::list<KineticDelaunay::CrossingData::EdgeIntersectionRef>::const_iterator {
            if (d_intersections.empty())
            {
              return d_intersections.end();
            }
            if (!ref_opt.has_value())
            {
              if (effective_list_forward)
              {
                return d_intersections.begin();
              }
              return std::prev(d_intersections.end());
            }
            auto it = std::find(d_intersections.begin(), d_intersections.end(), ref_opt.value());
            if (it == d_intersections.end())
            {
              return d_intersections.end();
            }
            if (effective_list_forward)
            {
              return std::next(it);
            }
            if (it == d_intersections.begin())
            {
              return d_intersections.end();
            }
            return std::prev(it);
          };

          for (size_t hop = 0; hop <= d_intersections.size(); ++hop)
          {
            const auto next_it = find_next_intersection(current_ref_opt);
            if (next_it == d_intersections.end())
            {
              break;
            }

            const KineticDelaunay::CrossingData::EdgeIntersectionRef candidate_ref = *next_it;
            const std::string candidate_log = formatCrossingIntersectionForLog(kin_del, std::make_optional(candidate_ref));

            append_traced_boundary_interval(last_ref_on_edge, candidate_ref);
            last_ref_on_edge = candidate_ref;

            const bool strand_inc = voronoi_edge_incident_to_strand(candidate_ref->voronoi_edge_id);
            const auto seg_it_start
              = strand_inc ? start_crossing_to_segment.find(candidate_ref) : start_crossing_to_segment.end();
            size_t next_segment_index = static_cast<size_t>(-1);
            if (seg_it_start != start_crossing_to_segment.end())
            {
              next_segment_index = seg_it_start->second;
            }
            if (next_segment_index == current_segment_index && next_segment_index != seed_index)
            {
              KINDS_DEBUG("Skipping intersection " << candidate_log
                                                 << " during walk: handoff would target current segment's own start.");
              KINDS_ERROR("Closing mesh: Delaunay handoff from ordered_seg[" << current_segment_index
                << "] targets its own start_crossing (only the seed loop may close onto seed).");
              exited_fail = true;
              exit_boundary_chain_to_voronoi = true;
              break;
            }
            // Close the loop only when Delaunay walk reaches the seed strip's *start* crossing (end of prior strip hands
            // off to start of next; the seed's start receives the last handoff).
            const bool closes_polygon_at_seed_crossing = (next_segment_index == seed_index);

            // Walking a new boundary half-edge resets current_ref_opt to nullopt; the first listed crossing on that edge
            // can be the seed segment's start crossing — already emitted as mesh_start. Closing the loop must not add
            // another mesh vertex there.
            if (!closes_polygon_at_seed_crossing)
            {
              const glm::dvec3 inter_pos = closingMeshVoronoiDelaunayCrossingPosition(
                t, candidate_ref->voronoi_edge_id, candidate_ref->delaunay_edge_id);
              if (std::isfinite(inter_pos.x) && std::isfinite(inter_pos.y))
              {
                const int nv = add_mesh_vertex(inter_pos);
                polygon_push(polygon, static_cast<size_t>(nv));
                logClosingMeshIntersectionVertex(kin_del, candidate_ref);
              }
              else
              {
                KINDS_DEBUG("Skipping intersection " << candidate_log
                                                   << " during walk: computed intersection position is non-finite.");
                KINDS_WARNING("Closing mesh: crossing on Delaunay edge has non-finite 2D position (boundary_he="
                  << boundary_he << ", crossing=" << candidate_log << ").");
              }
            }
            else
            {
              KINDS_DEBUG("Skipping intersection " << candidate_log
                                                 << " during walk: closes loop at seed start crossing; "
                                                    "seed start vertex already on polygon.");
            }

            if (strand_inc)
            {
              if (next_segment_index != static_cast<size_t>(-1))
              {
                if (closes_polygon_at_seed_crossing)
                {
                  walk_closure_crossing_ref = candidate_ref;
                  record_cap_handoff(current_segment_index, seed_index);
                  found_next_segment = true;
                  closed_walk_at_seed = true;
                  exit_boundary_chain_to_voronoi = true;
                  break;
                }
                record_cap_handoff(current_segment_index, next_segment_index);
                found_next_segment = true;
                if (next_segment_index == seed_index)
                {
                  closed_walk_at_seed = true;
                }
                else
                {
                  current_segment_index = next_segment_index;
                }
                exit_boundary_chain_to_voronoi = true;
                break;
              }
              KINDS_ERROR("Closing mesh: strand-incident Voronoi crossing on Delaunay boundary matches no ordered segment "
                          "start_crossing at boundary_he="
                          << boundary_he << ".");
              KINDS_DEBUG("Skipping intersection " << candidate_log
                                                 << " during walk: strand-incident crossing has no matching segment "
                                                    "start_crossing.");
              exited_fail = true;
              exit_boundary_chain_to_voronoi = true;
              break;
            }

            current_ref_opt = candidate_ref;
          }
        }

        if (exit_boundary_chain_to_voronoi)
        {
          break;
        }

        const int corner_vertex_id = graph.destination(boundary_he);
        if (!d_intersections.empty())
        {
          std::optional<KineticDelaunay::CrossingData::EdgeIntersectionRef> tail_ref = std::nullopt;
          if (current_ref_opt.has_value() && current_ref_opt.value()->delaunay_edge_id == (boundary_he / 2))
          {
            tail_ref = current_ref_opt;
          }
          if (tail_ref.has_value())
          {
            append_traced_boundary_interval(tail_ref, std::nullopt);
          }
        }
        if (corner_vertex_id >= 0)
        {
          if (static_cast<size_t>(corner_vertex_id) != strand_id)
          {
            KINDS_ERROR("Closing mesh: Delaunay boundary walk (half-edge " << boundary_he << ") reached corner vertex "
                                                                          << corner_vertex_id << " (expected strand "
                                                                          << strand_id << ").");
            exited_fail = true;
            const glm::dvec2 corner = kin_del.getPointAt(t, static_cast<size_t>(corner_vertex_id));
            const int cv = add_mesh_vertex(glm::dvec3(corner, t));
            polygon_push(polygon, static_cast<size_t>(cv));
            break;
          }
          else
          {
            append_strand_corner_if_needed();
          }
        }

        boundary_he = kin_del.nextOnComponentBoundaryId(boundary_he);
        current_ref_opt.reset();
      }

      if (exited_fail)
      {
        break;
      }
      if (closed_walk_at_seed)
      {
        exited_closed_loop = true;
        break;
      }
      if (!found_next_segment)
      {
        KINDS_ERROR("Closing mesh: Delaunay boundary walk from end of ordered_seg[" << current_segment_index
          << "] did not reach any strip's start_crossing (nor close the loop to ordered_seg[" << seed_index << "]).");
        exited_fail = true;
        break;
      }
    }

    if (exited_closed_loop)
    {
      for (size_t idx : walk_segment_order)
      {
        if (cap_start_hits[idx] != 1)
        {
          KINDS_ERROR("Closing mesh: ordered_seg[" << idx << "] start: expected exactly one incoming handoff (Delaunay or "
                                                     "Voronoi-vertex join), got "
                                                     << cap_start_hits[idx] << ".");
          exited_fail = true;
        }
        if (cap_end_handed_off[idx] != 1)
        {
          KINDS_ERROR("Closing mesh: ordered_seg[" << idx << "] end: expected exactly one outgoing handoff (Delaunay or "
                                                     "Voronoi-vertex join), got "
                                                     << cap_end_handed_off[idx] << ".");
          exited_fail = true;
        }
      }
    }
    else if (!exited_fail && !walk_segment_order.empty())
    {
      KINDS_ERROR("Closing mesh: polygon walk for seed ordered_seg[" << seed_index
        << "] ended without closing the loop (each strip end must hand off via a Delaunay boundary walk or a Voronoi "
           "circumcenter join matching dual Voronoi vertex ids to another strip's strand end or seed's strand-side "
           "vertex).");
      exited_fail = true;
    }

    {
      const MeshingData* end_seg = ordered_segments[last_walk_segment_index];
      const size_t ve_end = end_seg->closing_voronoi_edge_id;
      std::ostringstream head;
      head << "Ended polygon walk at ordered segment " << last_walk_segment_index << " of Voronoi edge ";
      if (ve_end == static_cast<size_t>(-1))
      {
        head << "?";
      }
      else
      {
        head << ve_end;
      }
      head << ".";
      if (exited_closed_loop && walk_closure_crossing_ref.has_value())
      {
        glm::dvec2 xy(0.0, 0.0);
        const bool have_xy = tryComputeCrossingIntersectionPosition2D(kin_del, walk_closure_crossing_ref, t, xy);
        const auto& ir = *walk_closure_crossing_ref.value();
        const glm::dvec3 geom
          = closingMeshVoronoiDelaunayCrossingPosition(t, ir.voronoi_edge_id, ir.delaunay_edge_id);
        head << " Loop closed at crossing " << formatCrossingIntersectionForLog(kin_del, walk_closure_crossing_ref) << ".";
        if (have_xy)
        {
          head << " CrossingData 2D=(" << xy.x << "," << xy.y << ")";
        }
        else
        {
          head << " CrossingData 2D=(unavailable)";
        }
        if (std::isfinite(geom.x) && std::isfinite(geom.y))
        {
          head << "; Voronoi-Delaunay chord intersection=(" << geom.x << "," << geom.y << ").";
        }
        else
        {
          head << "; Voronoi-Delaunay chord intersection=(non-finite).";
        }
      }
      else if (exited_closed_loop)
      {
        head << " Loop closed at seed mesh_start (Voronoi circumcenter join; no Delaunay closure crossing).";
      }
      else if (!exited_fail)
      {
        head << " (walk did not close to seed start crossing).";
      }
      KINDS_DEBUG(head.str());
    }

    if (polygon.size() >= 3 && exited_closed_loop && !exited_fail)
    {
      double area2 = 0.0;
      for (size_t i = 0; i < polygon.size(); ++i)
      {
        const glm::dvec3& p0 = mesh.getVertices()[polygon[i]];
        const glm::dvec3& p1 = mesh.getVertices()[polygon[(i + 1) % polygon.size()]];
        area2 += p0.x * p1.y - p1.x * p0.y;
      }
      if (area2 < 0.0)
      {
        std::reverse(polygon.begin() + 1, polygon.end());
      }
      result.polygons.push_back(std::move(polygon));
    }
  }

  return result;
}

void kinDS::SegmentBuilder::closingMeshLogUnmatchedOrderedSegments(
  size_t strand_id, double t, const std::vector<MeshingData*>& ordered_segments, const std::vector<bool>& segment_used)
{
  size_t unmatched_count = 0;
  for (size_t i = 0; i < segment_used.size(); ++i)
  {
    if (!segment_used[i])
    {
      ++unmatched_count;
    }
  }
  if (unmatched_count > 0)
  {
    KINDS_WARNING("Closing mesh: " << unmatched_count << " of " << segment_used.size()
                                   << " ordered Voronoi segment(s) not reached by polygon tracing (strand " << strand_id
                                   << " t=" << t << ").");
    for (size_t i = 0; i < segment_used.size(); ++i)
    {
      if (!segment_used[i])
      {
        const MeshingData* s = ordered_segments[i];
        KINDS_WARNING("Closing mesh: unmatched ordered segment [" << i << "] incident_edge=" << s->closing_incident_edge_index
                                                                  << " voronoi_edge="
                                                                  << (s->closing_voronoi_edge_id == static_cast<size_t>(-1)
                                                                        ? -1
                                                                        : static_cast<int>(s->closing_voronoi_edge_id))
                                                                  << ".");
      }
    }
  }
}

void kinDS::SegmentBuilder::closingMeshTriangulatePolygonsFan(
  VoronoiMesh& mesh, const std::vector<std::vector<size_t>>& polygons)
{
  for (const auto& polygon : polygons)
  {
    for (size_t i = 2; i < polygon.size(); ++i)
    {
      addMeshletTriangle(mesh, polygon[0], polygon[i - 1], polygon[i]);
    }
  }
}

size_t kinDS::SegmentBuilder::createClosingMesh(size_t strand_id, double t,
  const std::vector<BoundaryPoint>& boundary_polygon, const glm::dvec2& centroid,
  std::vector<BoundaryIntersectionInterval>* traced_boundary_intervals)
{
  KINDS_DEBUG("createClosingMesh strand " << strand_id << " t=" << t);
  const size_t num_incident_edges = closingMeshCountStrandIncidentEdges(strand_id);

  if (traced_boundary_intervals != nullptr)
  {
    traced_boundary_intervals->clear();
  }

  segment_mesh_pairs.push_back(MeshStructure::SegmentMeshPair {});

  VoronoiMesh mesh;
  std::vector<size_t> mesh_vertex_ids;
  mesh_vertex_ids.reserve(32);
  auto track_vertex = [&](const glm::dvec3& pos, std::optional<size_t> vv = std::nullopt) -> int
  {
    const int id = closingMeshAppendVertex(mesh, boundary_polygon, centroid, strand_id, t, pos, vv);
    mesh_vertex_ids.push_back(static_cast<size_t>(id));
    return id;
  };

  auto& graph = kin_del.getGraph();
  // Raw segments keyed by undirected Voronoi (dual) edge id; one bucket per incident finite edge.
  std::map<size_t, std::vector<MeshingData>> segments_by_voronoi_edge;
  std::list<MeshingData> closing_segments;

  // Trace every Voronoi edge incident to this strand and keep only portions that lie inside.
  // Use the same half-edge walk as `startNewMesh` (`computeCrossedHalfEdges`) so face inside/outside stays in sync;
  // then attach optional `CrossingData` edge-intersection iterators for Delaunay boundary tracing.
  int incident_edge_index = -1;
  for (HalfEdgeDelaunayGraph::IncidentEdgeIterator it = graph.incidentEdgesBegin(strand_id),
                                                   end = graph.incidentEdgesEnd(strand_id);
    it != end; ++it)
  {
    ++incident_edge_index;
    const size_t incident_he = *it;
    const size_t voronoi_edge_id = incident_he / 2;
    std::vector<MeshingData> bucket = closingMeshExtractRawSegmentsForVoronoiEdge(strand_id, t, mesh, boundary_polygon,
      centroid, incident_edge_index, incident_he, track_vertex);
    segments_by_voronoi_edge[voronoi_edge_id] = bucket;
    closing_segments.insert(closing_segments.end(), std::make_move_iterator(bucket.begin()),
      std::make_move_iterator(bucket.end()));
  }

  {
    size_t extraction_index = 0;
    for (const auto& seg : closing_segments)
    {
      const int ve_log
        = seg.closing_voronoi_edge_id == static_cast<size_t>(-1) ? -1 : static_cast<int>(seg.closing_voronoi_edge_id);
      KINDS_DEBUG("Closing mesh extraction, segment " << extraction_index << " (incident_edge="
                                                       << seg.closing_incident_edge_index << " voronoi_edge=" << ve_log
                                                       << " mesh_v " << seg.mesh_start_vertex_id << "->"
                                                       << seg.mesh_end_vertex_id << "):");
      logExtractionClosingMeshStripEndpoint(kin_del, t, extraction_index, seg.start_crossing, seg.closing_voronoi_edge_id);
      logExtractionClosingMeshStripEndpoint(kin_del, t, extraction_index, seg.end_crossing, seg.closing_voronoi_edge_id);
      ++extraction_index;
    }
  }

  ClosingMeshOrderedIndex index_data = closingMeshBuildOrderedIndex(closing_segments);

  closingMeshValidateOrderedSegmentGeometry(t, mesh, index_data.ordered_segments);

  ClosingMeshPolygonsTraceResult trace
    = closingMeshTraceCapPolygons(strand_id, t, num_incident_edges, mesh, boundary_polygon, centroid,
      std::move(mesh_vertex_ids), index_data.ordered_segments, index_data.start_crossing_to_segment);
  if (traced_boundary_intervals != nullptr)
  {
    *traced_boundary_intervals = trace.traced_boundary_intervals;
  }

  closingMeshLogUnmatchedOrderedSegments(strand_id, t, index_data.ordered_segments, trace.segment_used);
  closingMeshTriangulatePolygonsFan(mesh, trace.polygons);

  const size_t index = registerMeshletWithSuffix(std::move(mesh), std::string("_strand") + std::to_string(strand_id), t);
  segment_mesh_pair_last_left_and_right_vertex.emplace_back();
  segment_mesh_pair_last_left_and_right_vertex.back() = std::move(closing_segments);

  return index;
}

void kinDS::SegmentBuilder::accumulateSegmentProperties()
{
  // Iterate through all pairs and accumulate properties
  for (size_t pair_id = 0; pair_id < segment_mesh_pairs.size(); ++pair_id)
  {
    auto& pair = segment_mesh_pairs[pair_id];
    if (pair.segment_index0 != -1)
    {
      // make sure there is space left
      if (segment_properties[pair.segment_index0].neighbor_count >= MeshStructure::SegmentProperties::MAX_NEIGHBORS)
      {
        KINDS_ERROR("Exceeded maximum number of neighbors for segment: "
          << segment_properties[pair.segment_index0].neighbor_count
          << " >= " << MeshStructure::SegmentProperties::MAX_NEIGHBORS);
        // throw std::runtime_error("Exceeded maximum number of neighbors for segment.");
      }
      else
      {
        segment_properties[pair.segment_index0]
          .mesh_pair_indices[segment_properties[pair.segment_index0].neighbor_count] = pair_id; // add mesh pair index
        segment_properties[pair.segment_index0].neighbor_indices[segment_properties[pair.segment_index0].neighbor_count]
          = pair.segment_index1; // add neighbor
        segment_properties[pair.segment_index0].neighbor_count++;
      }
    }

    if (pair.segment_index1 != -1)
    {
      // make sure there is space left
      if (segment_properties[pair.segment_index1].neighbor_count >= MeshStructure::SegmentProperties::MAX_NEIGHBORS)
      {
        KINDS_ERROR("Exceeded maximum number of neighbors for segment: "
          << segment_properties[pair.segment_index1].neighbor_count
          << " >= " << MeshStructure::SegmentProperties::MAX_NEIGHBORS);
        // throw std::runtime_error("Exceeded maximum number of neighbors for segment.");
      }
      else
      {
        segment_properties[pair.segment_index1]
          .mesh_pair_indices[segment_properties[pair.segment_index1].neighbor_count] = pair_id; // add mesh pair index
        segment_properties[pair.segment_index1].neighbor_indices[segment_properties[pair.segment_index1].neighbor_count]
          = pair.segment_index0; // add neighbor
        segment_properties[pair.segment_index1].neighbor_count++;
      }
    }
  }
}

void kinDS::SegmentBuilder::splitComponent(
  size_t component_id, const std::vector<std::vector<size_t>>& new_components, double t)
{
  if (new_components.empty())
  {
    return;
  }

  // Update component data
  std::vector<size_t> component_ids(new_components.size(), -1);

  component_ids[0] = component_id;
  size_t new_component_id = kin_del.component_data.components.size();
  for (size_t i = 1; i < new_components.size(); i++)
  {
    component_ids[i] = new_component_id;
    new_component_id++;
  }

  size_t new_size = kin_del.component_data.components.size() + new_components.size() - 1;
  kin_del.component_data.components.resize(new_size);
  kin_del.component_data.component_boundaries.resize(new_size);
  kin_del.component_data.component_centroids.resize(new_size);

  std::vector<bool> he_visited(kin_del.getGraph().getHalfEdges().size(), false);

  for (size_t i = 0; i < new_components.size(); i++)
  {
    size_t cid = component_ids[i];

    for (size_t v : new_components[i])
    {
      kin_del.component_data.component_map[v] = cid;
    }

    kin_del.component_data.components[cid] = new_components[i];
    kin_del.component_data.component_boundaries[cid]
      = kin_del.extractComponentBoundaries(new_components[i], t, he_visited);
    kin_del.component_data.component_centroids[cid]
      = polygonCentroid(kin_del.component_data.component_boundaries[cid][0]);
    kin_del.component_data.component_last_updated[cid] = t;
  }
}

void SegmentBuilder::init()
{
  kin_del.registerEventCallbacks(
    section_callback_.get(), flip_callback_.get(), radius_callback_.get(), crossing_callback_.get());
  kin_del.registerSubdivisionEventCallback(subdivision_callback_.get());

  auto& graph = kin_del.getGraph();

  size_t strand_count = graph.getVertexCount();
  strand_to_segment_indices.resize(strand_count);
  half_edge_index_to_segment_mesh_pair_index.resize(graph.getHalfEdges().size(), -1);
  corner_to_cutoff_mesh_indices.resize(graph.getHalfEdges().size(), -1);

  // Initialize the strand geometries at t = 0.0
  double t = 0.0; // TODO: might be customized later

  // We need a ruled surface for each half-edge in the graph except those having the infinite vertex as
  // origin
  size_t half_edge_count = graph.getHalfEdges().size();

  // initialize segment mesh properties for each strand
  for (size_t strand_id = 0; strand_id < strand_count; ++strand_id)
  {
    size_t new_segment_id = segment_properties.size();
    MeshStructure::SegmentProperties properties;
    segment_properties.push_back(properties);
    strand_to_segment_indices[strand_id].push_back(new_segment_id);

    size_t component_index = kin_del.component_data.component_map[strand_id];
    const auto& component = kin_del.component_data.components[component_index];
    const auto& boundary_polygon = kin_del.component_data.component_boundaries[component_index][0];
    const auto& centroid = kin_del.component_data.component_centroids[component_index];
    // create a closing mesh
    size_t closing_mesh_index = createClosingMesh(strand_id, t, boundary_polygon, centroid);
    MeshStructure::SegmentMeshPair& segment_mesh_pair = segment_mesh_pairs[new_segment_id];
    segment_mesh_pair.segment_index0 = -1;
    segment_mesh_pair.segment_index1 = strand_to_segment_indices[strand_id].back();
  }

  // now go through all half-edges and create a segment mesh pair
  for (size_t i = 0; i < half_edge_count; i += 2)
  {
    startNewMesh(i, t);
  }

  // Create boundary interval meshes at t=0 from precomputed crossing data.
  for (size_t i = 0; i < half_edge_count; i += 2)
  {
    if (!kin_del.isOnComponentBoundary(i))
    {
      continue;
    }

    const size_t d_edge_id = i / 2;
    const auto& d_intersections = kin_del.getCrossingData().delaunay_edge_intersections[d_edge_id];
    if (d_intersections.empty())
    {
      continue;
    }

    const std::vector<KineticDelaunay::CrossingData::EdgeIntersectionRef> refs
      = getBoundaryIntersectionsInBoundaryOrder(d_edge_id);

    {
      const size_t first_cell
        = determineVoronoiCellForBoundaryIntersectionInterval(d_edge_id, std::nullopt, refs.front());
      startNewMeshFromIntersections(first_cell, t, std::nullopt, refs.front(), false, BoundaryEventType::Init,
        BoundarySegmentAction::NewSegment);
    }

    for (size_t k = 0; k + 1 < refs.size(); ++k)
    {
      const size_t mid_cell
        = determineVoronoiCellForBoundaryIntersectionInterval(d_edge_id, refs[k], refs[k + 1]);
      startNewMeshFromIntersections(
        mid_cell, t, refs[k], refs[k + 1], false, BoundaryEventType::Init, BoundarySegmentAction::NewSegment);
    }

    {
      const size_t last_cell = determineVoronoiCellForBoundaryIntersectionInterval(d_edge_id, refs.back(), std::nullopt);
      startNewMeshFromIntersections(last_cell, t, refs.back(), std::nullopt, false, BoundaryEventType::Init,
        BoundarySegmentAction::NewSegment);
    }
  }

  // initialize boundary mesh
  boundary_mesh_last_left_and_right_vertex.resize(half_edge_count, std::make_pair(-1, -1));
  half_edge_to_boundary_vertex_index.resize(half_edge_count, -1);
  addVoronoiTriangulationToBoundaryMesh(t, false, -0.01);
}

void SegmentBuilder::finalize(double t)
{
  updateBoundaries(t);

  // Finalize the segments by finishing all meshes
  auto& graph = kin_del.getGraph();
  size_t half_edge_count = graph.getHalfEdges().size();

  for (size_t i = 0; i < half_edge_count; i += 2)
  {
    auto vertex = graph.getHalfEdges()[i].origin;

    // fall back for infinite vertices
    if (vertex == -1)
    {
      vertex = graph.destination(i);
    }

    size_t component_index = kin_del.component_data.component_map[vertex];
    auto& boundary_points = kin_del.component_data.component_boundaries[component_index][0];

    finishMesh(i, t, boundary_points);
  }

  // Finalize boundary-interval meshes once more at the final time for all boundary Delaunay-edge sections.
  for (size_t i = 0; i < half_edge_count; i += 2)
  {
    if (!kin_del.isOnComponentBoundary(i))
    {
      continue;
    }
    const size_t d_edge_id = i / 2;
    const auto& d_intersections = kin_del.getCrossingData().delaunay_edge_intersections[d_edge_id];
    if (d_intersections.empty())
    {
      continue;
    }

    const std::vector<KineticDelaunay::CrossingData::EdgeIntersectionRef> refs
      = getBoundaryIntersectionsInBoundaryOrder(d_edge_id);

    {
      const size_t first_cell = determineVoronoiCellForBoundaryIntersectionInterval(d_edge_id, std::nullopt, refs.front());
      finishMeshFromIntersections(first_cell, t, std::nullopt, refs.front(), BoundaryEventType::Section,
        BoundarySegmentAction::SegmentCompleted);
    }
    for (size_t k = 0; k + 1 < refs.size(); ++k)
    {
      const size_t mid_cell = determineVoronoiCellForBoundaryIntersectionInterval(d_edge_id, refs[k], refs[k + 1]);
      finishMeshFromIntersections(
        mid_cell, t, refs[k], refs[k + 1], BoundaryEventType::Section, BoundarySegmentAction::SegmentCompleted);
    }
    {
      const size_t last_cell = determineVoronoiCellForBoundaryIntersectionInterval(d_edge_id, refs.back(), std::nullopt);
      finishMeshFromIntersections(
        last_cell, t, refs.back(), std::nullopt, BoundaryEventType::Section, BoundarySegmentAction::SegmentCompleted);
    }
  }

  // finalize closing meshes
  for (size_t strand_id = 0; strand_id < graph.getVertexCount(); ++strand_id)
  {
    // create a closing mesh
    size_t component_index = kin_del.component_data.component_map[strand_id];
    auto& boundary_points = kin_del.component_data.component_boundaries[component_index][0];
    auto& centroid = kin_del.component_data.component_centroids[component_index];
    size_t closing_mesh_index = createClosingMesh(strand_id, t, boundary_points, centroid);
    MeshStructure::SegmentMeshPair& segment_mesh_pair = segment_mesh_pairs[closing_mesh_index];
    segment_mesh_pair.segment_index0 = strand_to_segment_indices[strand_id].back();
    segment_mesh_pair.segment_index1 = -1;
  }

  accumulateSegmentProperties();

  addVoronoiTriangulationToBoundaryMesh(t, true, 0.01);

  // compute normals
  for (auto& meshlet : meshes)
  {
    meshlet.computeNormals(NormalMode::PerTriangleCorner);
  }

  auto remap1 = boundary_mesh.mergeDuplicateVertices();
  boundary_mesh.removeDegenerateTriangles();
  auto remap2 = boundary_mesh.removeIsolatedVertices();
  boundary_mesh.computeNormals(NormalMode::PerTriangleCorner);

  // Update boundary vertex to strand id mapping
  std::vector<size_t> new_boundary_vertex_to_strand_id;
  new_boundary_vertex_to_strand_id.resize(boundary_mesh.getVertices().size());
  for (size_t old_index = 0; old_index < boundary_vertex_to_strand_id.size(); ++old_index)
  {
    size_t new_index = remap2[remap1[old_index]];
    if (new_index != size_t(-1))
    {
      new_boundary_vertex_to_strand_id[new_index] = boundary_vertex_to_strand_id[old_index];
    }
  }

  boundary_vertex_to_strand_id.swap(new_boundary_vertex_to_strand_id);

  finalized = true; // Set the finalized flag to true
}

std::vector<VoronoiMesh> kinDS::SegmentBuilder::extractMeshes() const { return meshes; }

std::pair<std::vector<VoronoiMesh>, std::vector<std::vector<int>>> kinDS::SegmentBuilder::extractSegmentMeshlets(
  bool merge_by_segment) const
{
  if (!merge_by_segment)
  {
    std::vector<VoronoiMesh> raw_meshlets = meshes;
    std::vector<std::vector<int>> raw_neighbors;
    raw_neighbors.reserve(raw_meshlets.size());

    for (size_t meshlet_index = 0; meshlet_index < raw_meshlets.size(); ++meshlet_index)
    {
      std::vector<int> neighbors_for_meshlet;
      neighbors_for_meshlet.assign(raw_meshlets[meshlet_index].getTriangleCount(), -1);
      raw_neighbors.push_back(std::move(neighbors_for_meshlet));
    }

    return std::make_pair(raw_meshlets, raw_neighbors);
  }

  std::vector<VoronoiMesh> meshlets;
  std::vector<std::vector<int>> neighbor_segments; // accessed as [segment_id][triangle_index]
  for (size_t segment_id = 0; segment_id < segment_properties.size(); ++segment_id)
  {
    VoronoiMesh segment_mesh(MeshletExportMaterialNames);
    std::vector<int> neighbor_segments_for_meshlet;
    const auto& properties = segment_properties[segment_id];
    double earliest_creation = std::numeric_limits<double>::quiet_NaN();
    auto append_oriented_mesh = [&](VoronoiMesh mesh, size_t seg0, size_t seg1) {
      const double mesh_ct = mesh.getCreationKineticTime();
      if (std::isfinite(mesh_ct))
      {
        if (!std::isfinite(earliest_creation) || mesh_ct < earliest_creation)
        {
          earliest_creation = mesh_ct;
        }
      }
      if (seg0 != segment_id)
      {
        mesh.flipOrientation();
      }
      const int neighbor = (seg0 == segment_id) ? static_cast<int>(seg1) : static_cast<int>(seg0);
      neighbor_segments_for_meshlet.insert(neighbor_segments_for_meshlet.end(), mesh.getTriangleCount(), neighbor);
      segment_mesh += mesh;
    };

    // Regular meshlets referenced through segment_properties.
    for (size_t neighbor_index = 0; neighbor_index < properties.neighbor_count; ++neighbor_index)
    {
      size_t mesh_pair_index = properties.mesh_pair_indices[neighbor_index];
      VoronoiMesh mesh = meshes[mesh_pair_index];
      append_oriented_mesh(mesh, segment_mesh_pairs[mesh_pair_index].segment_index0, segment_mesh_pairs[mesh_pair_index].segment_index1);
    }

    // Boundary-interval meshlets are tracked in a separate pair/index space and must be merged explicitly.
    // Unlike regular Voronoi-edge strips, each boundary interval belongs to exactly one segment (its Voronoi cell),
    // and should not be orientation-flipped during merged export.
    for (size_t pair_idx = 0; pair_idx < intersection_segment_mesh_pairs.size(); ++pair_idx)
    {
      if (pair_idx >= intersection_meshes.size() || pair_idx >= intersection_mesh_pair_metadata.size())
      {
        continue;
      }

      const size_t owner_segment_id = intersection_mesh_pair_metadata[pair_idx].owner_segment_id;
      if (owner_segment_id == static_cast<size_t>(-1))
      {
        continue;
      }
      if (owner_segment_id != segment_id)
      {
        continue;
      }

      VoronoiMesh mesh = intersection_meshes[pair_idx];
      const double mesh_ct = mesh.getCreationKineticTime();
      if (std::isfinite(mesh_ct))
      {
        if (!std::isfinite(earliest_creation) || mesh_ct < earliest_creation)
        {
          earliest_creation = mesh_ct;
        }
      }
      neighbor_segments_for_meshlet.insert(neighbor_segments_for_meshlet.end(), mesh.getTriangleCount(), -1);
      segment_mesh += mesh;
    }
    if (std::isfinite(earliest_creation))
    {
      segment_mesh.setCreationKineticTime(earliest_creation);
    }
    neighbor_segments.push_back(neighbor_segments_for_meshlet);
    segment_mesh.mergeDuplicateVertices(1e-4);
    meshlets.push_back(segment_mesh);
  }

  return std::make_pair(meshlets, neighbor_segments);
}

std::vector<std::string> kinDS::SegmentBuilder::extractSegmentMeshletExportSuffixes(bool merge_by_segment) const
{
  if (!merge_by_segment)
  {
    return meshlet_export_suffixes;
  }

  std::vector<std::string> suffixes;
  suffixes.reserve(segment_properties.size());
  for (size_t segment_id = 0; segment_id < segment_properties.size(); ++segment_id)
  {
    suffixes.push_back(std::string("_segment") + std::to_string(segment_id));
  }
  return suffixes;
}

std::vector<VoronoiMesh> kinDS::SegmentBuilder::extractBoundaryIntervalMeshlets() const
{
  return intersection_meshes;
}

std::vector<std::string> kinDS::SegmentBuilder::extractBoundaryIntervalMeshletExportSuffixes() const
{
  std::vector<std::string> out = intersection_meshlet_export_suffixes;
  if (out.size() < intersection_mesh_pair_metadata.size())
  {
    out.resize(intersection_mesh_pair_metadata.size());
  }
  const size_t count = std::min(out.size(), intersection_mesh_pair_metadata.size());
  for (size_t i = 0; i < count; ++i)
  {
    const auto& md = intersection_mesh_pair_metadata[i];
    std::string suffix = out[i];
    suffix += std::string("_cell") + std::to_string(md.voronoi_cell_id);
    if (md.start_delaunay_edge_id != static_cast<size_t>(-1))
    {
      suffix += std::string("_de0") + std::to_string(md.start_delaunay_edge_id);
    }
    if (md.end_delaunay_edge_id != static_cast<size_t>(-1))
    {
      suffix += std::string("_de1") + std::to_string(md.end_delaunay_edge_id);
    }
    out[i] = std::move(suffix);
  }
  return out;
}

const VoronoiMesh& kinDS::SegmentBuilder::getBoundaryMesh() const { return boundary_mesh; }

const std::vector<size_t>& kinDS::SegmentBuilder::getBoundaryVertexToStrandId() const
{
  return boundary_vertex_to_strand_id;
}

const std::vector<std::vector<size_t>>& kinDS::SegmentBuilder::getStrandToSegmentIndices() const
{
  return strand_to_segment_indices;
}

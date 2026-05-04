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
#include <limits>
#include <glm/gtx/exterior_product.hpp>
#include <map>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <optional>
#include <sstream>

using namespace kinDS;

namespace
{
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

void kinDS::SegmentBuilder::finishMesh(size_t he_id, double t, const std::vector<BoundaryPoint>& boundary_points)
{
  // check if half-edge is infinite
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

  auto& last_segments = segment_mesh_pair_last_left_and_right_vertex[segment_mesh_pair_index];
  if (last_segments.empty())
  {
    meshletDiagnosticLogLine("finish_mesh_noop", he_id, t, "last_segments empty (no extension)");
    // this can happen if the first Voronoi vertex was discarded because it was outside of the boundary, in this case we
    // just don't create any triangles
    return;
  }

  VoronoiMesh& mesh = meshes[segment_mesh_pair_index];
  const size_t even_id = he_id & ~1;
  const size_t odd_id = even_id + 1;
  const size_t voronoi_edge_id = even_id / 2;
  auto& he = kin_del.getGraph().getHalfEdges()[even_id];
  glm::dvec2 centroid = polygonCentroid(boundary_points);

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

  size_t v = he.origin;
  if (v == -1)
  {
    v = kin_del.getGraph().destination(even_id);
  }

  auto endpoint_position_at_t = [&](const MeshingData& segment, bool at_start) -> glm::dvec3
  {
    const char* endpoint_label = at_start ? "start" : "end";
    const int ep_he = at_start ? segment.start_half_edge_id : segment.end_half_edge_id;
    if (ep_he < 0)
    {
      // Case 1: Open strip endpoint (no boundary half-edge). Use canonical Voronoi circumcenter on this edge.
      const glm::dvec3 p = computeVoronoiVertex(at_start ? even_id : odd_id, t);
      if (!std::isfinite(p.x) || !std::isfinite(p.y) || !std::isfinite(p.z) || (p.x == 0.0 && p.y == 0.0))
      {
        KINDS_WARNING("finishMesh endpoint_position_at_t(" << endpoint_label
          << "): degenerate circumcenter endpoint for voronoi_edge=" << voronoi_edge_id << " t=" << t << " -> (" << p.x
          << ", " << p.y << ", " << p.z << ").");
      }
      return p;
    }
    const auto& ref = at_start ? segment.start_crossing : segment.end_crossing;
    glm::dvec2 xy;
    if (ref.has_value() && tryComputeCrossingIntersectionPosition2D(kin_del, ref, t, xy))
    {
      // Case 2: Crossing iterator resolved in CrossingData; use its current geometric 2D position.
      const glm::dvec3 p(xy, t);
      if (!std::isfinite(p.x) || !std::isfinite(p.y) || !std::isfinite(p.z) || (p.x == 0.0 && p.y == 0.0))
      {
        KINDS_WARNING("finishMesh endpoint_position_at_t(" << endpoint_label
          << "): degenerate CrossingData endpoint for voronoi_edge=" << voronoi_edge_id << " t=" << t << " -> (" << p.x
          << ", " << p.y << ", " << p.z << ").");
      }
      return p;
    }
    // Case 3: Crossing ref missing/stale for this boundary half-edge. Fall back to geometric Voronoi-Delaunay chord
    // intersection using the endpoint Delaunay edge id.
    const size_t delaunay_edge_id = static_cast<size_t>(ep_he) / 2;
    const glm::dvec3 p = closingMeshVoronoiDelaunayCrossingPosition(t, voronoi_edge_id, delaunay_edge_id);
    if (!std::isfinite(p.x) || !std::isfinite(p.y) || !std::isfinite(p.z) || (p.x == 0.0 && p.y == 0.0))
    {
      KINDS_WARNING("finishMesh endpoint_position_at_t(" << endpoint_label
        << "): degenerate fallback chord endpoint for voronoi_edge=" << voronoi_edge_id
        << " d_edge=" << delaunay_edge_id << " t=" << t << " -> (" << p.x << ", " << p.y << ", " << p.z << ").");
    }
    return p;
  };

  size_t tris_before = mesh.getTriangleCount();
  size_t processable_strips = 0;
  for (const auto& segment : last_segments)
  {
    if (segment.mesh_start_vertex_id >= 0 && segment.mesh_end_vertex_id >= 0)
    {
      ++processable_strips;
    }
  }
  size_t loops_ran = 0;

  for (auto& segment : last_segments)
  {
    if (segment.mesh_start_vertex_id < 0 || segment.mesh_end_vertex_id < 0)
    {
      continue;
    }

    ++loops_ran;
    refreshMeshingDataCrossingRefs(segment, voronoi_edge_id);

    const glm::dvec3 new_start_pos = endpoint_position_at_t(segment, true);
    const glm::dvec3 new_end_pos = endpoint_position_at_t(segment, false);

    const auto& graph_finish = kin_del.getGraph();
    const std::optional<size_t> start_vv = segment.start_half_edge_id < 0
      ? std::optional<size_t>(graph_finish.getHalfEdges()[even_id].face)
      : std::nullopt;
    const std::optional<size_t> end_vv = segment.end_half_edge_id < 0
      ? std::optional<size_t>(graph_finish.getHalfEdges()[odd_id].face)
      : std::nullopt;

    const size_t new_start_vertex_index = addMeshletVertex(mesh, boundary_points, centroid, new_start_pos, v, t, start_vv);
    const size_t new_end_vertex_index = addMeshletVertex(mesh, boundary_points, centroid, new_end_pos, v, t, end_vv);

    const size_t last_left = static_cast<size_t>(segment.mesh_start_vertex_id);
    const size_t last_right = static_cast<size_t>(segment.mesh_end_vertex_id);
    if (last_left == last_right)
    {
      addMeshletTriangle(mesh, new_start_vertex_index, last_right, new_end_vertex_index);
    }
    else if (mesh.getVertices()[last_left][2] < mesh.getVertices()[last_right][2])
    {
      addMeshletTriangle(mesh, last_left, last_right, new_start_vertex_index);
      addMeshletTriangle(mesh, new_start_vertex_index, last_right, new_end_vertex_index);
    }
    else
    {
      addMeshletTriangle(mesh, last_left, last_right, new_end_vertex_index);
      addMeshletTriangle(mesh, last_left, new_end_vertex_index, new_start_vertex_index);
    }

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

void SegmentBuilder::startNewMesh(size_t half_edge_id, double t, bool reuse_existing_pair_and_mesh)
{
  // check if half-edge is infinite
  if (kin_del.getGraph().isInfinite(half_edge_id) && kin_del.computeBoundaryOnTheFly())
  {
    meshletDiagnosticLogLine("start_new_mesh_skip", half_edge_id, t, "reason=infinite_boundary_on_the_fly");
    return;
  }
  
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

  size_t vertex = std::max(he.origin, twin_he.origin);
  size_t component_id = kin_del.component_data.component_map[vertex];

  std::vector<bool> he_visited(graph.getHalfEdges().size(), false);
  updateBoundary(t, he_visited, component_id);

  auto& boundary_polygon = kin_del.component_data.component_boundaries[component_id][0];
  auto& centroid = kin_del.component_data.component_centroids[component_id];

  // New pair: build into a fresh mesh and register. Reuse: append strip vertices onto the existing mesh; strip
  // metadata below is cleared and rebuilt — mesh vertex indices from addMeshletVertex are absolute (append order).
  VoronoiMesh mesh_local;
  const bool reuse_in_place = !created_new_pair && segment_mesh_pair_index < meshes.size();
  VoronoiMesh& mesh = reuse_in_place ? meshes[segment_mesh_pair_index] : mesh_local;

  glm::dvec3 left_vertex = computeVoronoiVertex(even_id, t);
  glm::dvec3 right_vertex = computeVoronoiVertex(odd_id, t);

  const size_t voronoi_edge_id = even_id / 2;

  // Track how the Voronoi edge between these two vertices intersects with the boundary, we will need this information
  // to correctly build the boundary mesh
  size_t left_voronoi_vertex_id = kin_del.getGraph().getHalfEdges()[even_id].face;
  size_t left_containing_tri_id = kin_del.getCrossingDataContainingTriId(left_voronoi_vertex_id);

  size_t right_voronoi_vertex_id = kin_del.getGraph().getHalfEdges()[odd_id].face;

  // Now go through all faces
  const bool initial_left_inside = kin_del.getFaceInside(left_containing_tri_id);
  bool inside = initial_left_inside;

  if (segment_mesh_pair_index < segment_mesh_pair_last_left_and_right_vertex.size())
  {
    segment_mesh_pair_last_left_and_right_vertex[segment_mesh_pair_index].clear();
  }
  else
  {
    segment_mesh_pair_last_left_and_right_vertex.emplace_back();
  }
  auto& segments_for_pair = segment_mesh_pair_last_left_and_right_vertex[segment_mesh_pair_index];
  // If already inside, the left Voronoi vertex is the first one to be added to the mesh
  if (inside || !kin_del.computeBoundaryOnTheFly())
  {
    size_t vertex_index = addMeshletVertex(
      mesh, boundary_polygon, centroid, left_vertex, he.origin, t, std::optional<size_t>(left_voronoi_vertex_id));
    segments_for_pair.emplace_back(MeshingData { static_cast<int>(vertex_index), -1, -1, -1 });
  }

  if (kin_del.computeBoundaryOnTheFly())
  {
    // KINDS_DEBUG("Determining edge crossings of voronoi edge from " << left_voronoi_vertex_id << " to " <<
    // right_voronoi_vertex_id);
    // TODO: I don't think we need this anymore, we can use the edge intersection list instead
    std::vector<size_t> crossed_half_edges
      = kin_del.computeCrossedHalfEdges(left_containing_tri_id, right_vertex, left_vertex, t).first;
    // KINDS_DEBUG("Number of crossed half edges: " << crossed_half_edges.size() << " with voronoi vertex ids: " <<
    // left_voronoi_vertex_id << " and " << right_voronoi_vertex_id);
    for (size_t crossed_he_id : crossed_half_edges)
    {
      // KINDS_DEBUG("Crossed half edge: " << crossed_he_id);
      size_t prev_face_id = graph.getHalfEdges()[crossed_he_id].face;
      size_t next_face_id = graph.getHalfEdges()[crossed_he_id ^ 1].face;

      bool next_inside = kin_del.getFaceInside(next_face_id);

      if (inside != next_inside)
      {
        // Compute the intersection point of the Voronoi edge with the boundary corresponding to this half-edge
        glm::dvec2 intersection_point
          = intersectSegments(glm::dvec2(left_vertex[0], left_vertex[1]), glm::dvec2(right_vertex[0], right_vertex[1]),
            kin_del.getPointAt(t, graph.getHalfEdges()[crossed_he_id].origin),
            kin_del.getPointAt(t, graph.getHalfEdges()[crossed_he_id ^ 1].origin));

        // Add the intersection point as a vertex to the mesh
        size_t vertex_index
          = addMeshletVertex(mesh, boundary_polygon, centroid, glm::dvec3(intersection_point, t), -1, t);

        if (inside)
        {
          // leaving the inside, add the intersection point as second of the last added pair
          auto& back_seg = segments_for_pair.back();
          back_seg.mesh_end_vertex_id = static_cast<int>(vertex_index);
          back_seg.end_half_edge_id = static_cast<int>(crossed_he_id);
          back_seg.end_crossing = closingMeshFindVoronoiEdgeIntersection(voronoi_edge_id, crossed_he_id);
        }
        else
        {
          // Use the twin half-edge as it is located inside the component.
          MeshingData seg { static_cast<int>(vertex_index), -1, static_cast<int>(crossed_he_id ^ 1), -1 };
          seg.start_crossing = closingMeshFindVoronoiEdgeIntersection(voronoi_edge_id, crossed_he_id);
          segments_for_pair.emplace_back(std::move(seg));
        }
      }

      // current_face_id = next_face_id;
      inside = next_inside;
    }
  }

  if (inside || !kin_del.computeBoundaryOnTheFly())
  {
    // If we end inside, the right Voronoi vertex is the last one to be added to the mesh
    size_t vertex_index = addMeshletVertex(
      mesh, boundary_polygon, centroid, right_vertex, twin_he.origin, t, std::optional<size_t>(right_voronoi_vertex_id));
    segments_for_pair.back().mesh_end_vertex_id = static_cast<int>(vertex_index);
  }

  for (auto& seg : segments_for_pair)
  {
    refreshMeshingDataCrossingRefs(seg, voronoi_edge_id);
  }

  meshletDiagnosticWarnIfUnexpectedEmptyAfterStartNewMesh(even_id, t, initial_left_inside, mesh, segments_for_pair);

  if (reuse_in_place)
  {
    mesh.setCreationKineticTime(t);
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

  assert(segment_mesh_pairs.size() == segment_mesh_pair_last_left_and_right_vertex.size());
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

size_t kinDS::SegmentBuilder::addMeshletTriangle(VoronoiMesh& mesh, size_t u, size_t v, size_t w)
{
  return mesh.addTriangle(u, v, w, u, v, w); // For meshlets, the UVs are assigned per vertex so the indices match
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
  std::optional<size_t> meshlet_voronoi_vertex_for_alpha_check)
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
  size_t index = mesh.addVertex(vertex);
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

        // Walk every crossing along this directed Delaunay half-edge after `current_ref_opt`: each crossing gets a
        // mesh vertex. Strand-incident crossings that match an ordered segment start_crossing hand off to that segment;
        // others only advance along the same Delaunay edge (no skipping to the corner).
        bool exit_boundary_chain_to_voronoi = false;
        if (!d_intersections.empty())
        {
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

size_t kinDS::SegmentBuilder::createClosingMesh(
  size_t strand_id, double t, const std::vector<BoundaryPoint>& boundary_polygon, const glm::dvec2& centroid)
{
  KINDS_DEBUG("createClosingMesh strand " << strand_id << " t=" << t);
  const size_t num_incident_edges = closingMeshCountStrandIncidentEdges(strand_id);

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
    VoronoiMesh segment_mesh;
    std::vector<int> neighbor_segments_for_meshlet;
    const auto& properties = segment_properties[segment_id];
    double earliest_creation = std::numeric_limits<double>::quiet_NaN();
    for (size_t neighbor_index = 0; neighbor_index < properties.neighbor_count; ++neighbor_index)
    {
      size_t mesh_pair_index = properties.mesh_pair_indices[neighbor_index];
      VoronoiMesh mesh = meshes[mesh_pair_index];
      const double mesh_ct = mesh.getCreationKineticTime();
      if (std::isfinite(mesh_ct))
      {
        if (!std::isfinite(earliest_creation) || mesh_ct < earliest_creation)
        {
          earliest_creation = mesh_ct;
        }
      }
      if (segment_mesh_pairs[mesh_pair_index].segment_index0 != segment_id)
      {
        mesh.flipOrientation();
      }
      // Append the mesh to the segment mesh
      neighbor_segments_for_meshlet.insert(
        neighbor_segments_for_meshlet.end(), mesh.getTriangleCount(), properties.neighbor_indices[neighbor_index]);
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

const VoronoiMesh& kinDS::SegmentBuilder::getBoundaryMesh() const { return boundary_mesh; }

const std::vector<size_t>& kinDS::SegmentBuilder::getBoundaryVertexToStrandId() const
{
  return boundary_vertex_to_strand_id;
}

const std::vector<std::vector<size_t>>& kinDS::SegmentBuilder::getStrandToSegmentIndices() const
{
  return strand_to_segment_indices;
}

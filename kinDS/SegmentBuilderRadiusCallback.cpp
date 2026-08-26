#include "SegmentBuilderRadiusCallback.hpp"

#include "SegmentBuilder.hpp"
#include "DebugExportFormatting.hpp"
#include "KineticDelaunayCrossingEvent.hpp"
#include "SegmentBuilderVisualDebug.hpp"
#include "VisualDebugHighlight.hpp"
#include "Logger.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <optional>
#include <sstream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace kinDS
{
namespace
{
std::string radiusDebugNumberLiteral(double value)
{
  std::ostringstream o;
  o << std::setprecision(kDebugExportTimePrecision) << std::showpoint << value;
  return o.str();
}

std::string radiusDebugSvgEscape(const std::string& text)
{
  std::string out;
  out.reserve(text.size());
  for (char ch : text)
  {
    switch (ch)
    {
    case '&':
      out += "&amp;";
      break;
    case '<':
      out += "&lt;";
      break;
    case '>':
      out += "&gt;";
      break;
    case '"':
      out += "&quot;";
      break;
    default:
      out += ch;
      break;
    }
  }
  return out;
}

struct RadiusRingWalkDebugVertex
{
  glm::dvec2 xy {};
  std::string source = "unknown";
  std::optional<size_t> strand_id {};
  std::optional<size_t> voronoi_vertex_id {};
  std::optional<size_t> delaunay_edge_id {};
  std::optional<size_t> voronoi_edge_id {};
  /// Set when this walked vertex is known to be wrong for the traced cell (e.g. wrong site).
  bool incorrect = false;
  std::string note {};
};

std::filesystem::path makeRadiusRingWalkFailDebugPath(const KineticDelaunay& kin_del, EventTime occurrence_time,
  std::optional<size_t> runtime_branch_id, size_t delaunay_face_id, size_t strand_cell_id, size_t debug_counter,
  const char* extension)
{
  const std::string filename = formatDebugExportTimeToken(occurrence_time) + "_radius_ring_walk_FAIL_face"
    + std::to_string(delaunay_face_id) + "_strand" + std::to_string(strand_cell_id) + "_"
    + std::to_string(debug_counter) + extension;
  const std::string branch_folder = runtime_branch_id.has_value()
    ? ("branch" + std::to_string(runtime_branch_id.value()))
    : kVisualDebugUnresolvedBranchFolder;
  std::filesystem::path filepath = std::filesystem::path(branch_folder) / filename;
  if (const std::optional<std::filesystem::path>& output_root = kin_del.getVisualDebugOutputRoot();
    output_root.has_value())
  {
    filepath = *output_root / filepath;
  }
  if (filepath.has_parent_path())
  {
    std::filesystem::create_directories(filepath.parent_path());
  }
  return filepath;
}

void appendRadiusRingWalkVertexSourceFields(std::ostream& out, const RadiusRingWalkDebugVertex& vert)
{
  out << " source=" << vert.source;
  if (vert.strand_id.has_value())
  {
    out << " strand_id=" << vert.strand_id.value();
  }
  if (vert.voronoi_vertex_id.has_value())
  {
    out << " voronoi_vertex_id=" << vert.voronoi_vertex_id.value();
  }
  if (vert.delaunay_edge_id.has_value())
  {
    out << " delaunay_edge_id=" << vert.delaunay_edge_id.value();
  }
  if (vert.voronoi_edge_id.has_value())
  {
    out << " voronoi_edge_id=" << vert.voronoi_edge_id.value();
  }
  if (vert.incorrect)
  {
    out << " status=incorrect";
  }
  if (!vert.note.empty())
  {
    out << " note=\"" << vert.note << '"';
  }
}

void writeRadiusRingWalkFailDebugTxt(const std::filesystem::path& filepath, EventTime occurrence_time,
  std::optional<size_t> runtime_branch_id, size_t delaunay_face_id, size_t strand_cell_id,
  const std::string& fail_reason, const std::vector<RadiusRingWalkDebugVertex>& ring,
  const std::vector<std::string>& incorrect_vertices, const std::vector<std::string>& unmatched_vertices)
{
  std::ofstream out(filepath);
  if (!out)
  {
    KINDS_WARNING("radius ring walk FAIL: failed to open debug TXT " << filepath.generic_string()
      << " runtime_branch="
      << (runtime_branch_id.has_value() ? std::to_string(runtime_branch_id.value()) : "unresolved")
      << " strand_cell=" << strand_cell_id << " face=" << delaunay_face_id << " t=" << occurrence_time);
    return;
  }

  out << "# tag=radius_ring_walk_FAIL\n";
  out << "# occurrence_time=" << occurrence_time << '\n';
  out << "# runtime_branch_id="
      << (runtime_branch_id.has_value() ? std::to_string(runtime_branch_id.value()) : "unresolved") << '\n';
  out << "# delaunay_face_id=" << delaunay_face_id << '\n';
  out << "# strand_cell_id=" << strand_cell_id << '\n';
  out << "# fail_reason=" << fail_reason << '\n';
  out << "# ring_vertex_count=" << ring.size() << '\n';
  out << "# incorrect_count=" << incorrect_vertices.size() << '\n';
  out << "# unmatched_count=" << unmatched_vertices.size() << "\n\n";

  if (!incorrect_vertices.empty())
  {
    out << "## incorrectly_found_vertices\n";
    for (const std::string& entry : incorrect_vertices)
    {
      out << entry << '\n';
    }
    out << '\n';
  }

  if (!unmatched_vertices.empty())
  {
    out << "## unmatched_vertices\n";
    for (const std::string& entry : unmatched_vertices)
    {
      out << entry << '\n';
    }
    out << '\n';
  }

  out << "## ring=walked vertex_count=" << ring.size() << '\n';
  for (size_t i = 0; i < ring.size(); ++i)
  {
    for (size_t j = 0; j < i; ++j)
    {
      if (ring[i].xy.x == ring[j].xy.x && ring[i].xy.y == ring[j].xy.y)
      {
        out << "# duplicate_xy: index=" << i << " matches index=" << j << '\n';
      }
    }
  }

  out << "# columns: index profile_x profile_y <source meta...>\n";
  for (size_t i = 0; i < ring.size(); ++i)
  {
    const RadiusRingWalkDebugVertex& vert = ring[i];
    out << i << ' ' << radiusDebugNumberLiteral(vert.xy.x) << ' ' << radiusDebugNumberLiteral(vert.xy.y);
    appendRadiusRingWalkVertexSourceFields(out, vert);
    out << '\n';
  }
  out << '\n';

  KINDS_WARNING("radius ring walk FAIL: wrote debug TXT to " << filepath.generic_string()
    << " runtime_branch="
    << (runtime_branch_id.has_value() ? std::to_string(runtime_branch_id.value()) : "unresolved")
    << " strand_cell=" << strand_cell_id << " face=" << delaunay_face_id << " t=" << occurrence_time);
}

/// Debug dump for radius cell ring-walk failure (before any ear-clip). Filename encodes ring_walk_FAIL.
/// Also writes the common event-style kinetic SVG with the affected face/sites/edges/VVs highlighted.
void writeRadiusRingWalkFailDebug(KineticDelaunay& kin_del, EventTime occurrence_time,
  std::optional<size_t> runtime_branch_id, size_t delaunay_face_id, size_t strand_cell_id,
  const std::string& fail_reason, std::vector<RadiusRingWalkDebugVertex> ring,
  std::vector<std::string> incorrect_vertices, std::vector<std::string> unmatched_vertices)
{
  // Annotate walked sites / Voronoi verts that clearly do not belong to this cell/face.
  const auto& crossing_data = kin_del.getCrossingData();
  for (RadiusRingWalkDebugVertex& vert : ring)
  {
    if (vert.source == "site" && vert.strand_id.has_value() && vert.strand_id.value() != strand_cell_id)
    {
      vert.incorrect = true;
      if (vert.note.empty())
      {
        vert.note = "site strand != traced cell";
      }
      incorrect_vertices.push_back("site strand_id=" + std::to_string(vert.strand_id.value())
        + " (expected cell " + std::to_string(strand_cell_id) + ")");
    }
    else if (vert.source == "Voronoi vertex" && vert.voronoi_vertex_id.has_value())
    {
      const size_t vv = vert.voronoi_vertex_id.value();
      if (!crossing_data.isVoronoiVertexRegistered(vv)
        || crossing_data.getContainingTriId(vv) != delaunay_face_id)
      {
        vert.incorrect = true;
        if (vert.note.empty())
        {
          vert.note = "Voronoi vertex not contained in affected triangle";
        }
        incorrect_vertices.push_back(
          "voronoi_vertex_id=" + std::to_string(vv) + " (not in face " + std::to_string(delaunay_face_id) + ")");
      }
    }
  }

  // Cell site should appear on a successful closed ring; list it as unmatched when absent.
  bool ring_contains_cell_site = false;
  for (const RadiusRingWalkDebugVertex& vert : ring)
  {
    if (vert.source == "site" && vert.strand_id.has_value() && vert.strand_id.value() == strand_cell_id)
    {
      ring_contains_cell_site = true;
      break;
    }
  }
  if (!ring_contains_cell_site)
  {
    unmatched_vertices.push_back("site strand_id=" + std::to_string(strand_cell_id) + " (not reached by walk)");
  }

  // Deduplicate diagnostic lines while preserving order.
  auto unique_preserve = [](std::vector<std::string>& entries)
  {
    std::unordered_set<std::string> seen;
    std::vector<std::string> out;
    out.reserve(entries.size());
    for (std::string& entry : entries)
    {
      if (seen.insert(entry).second)
      {
        out.push_back(std::move(entry));
      }
    }
    entries = std::move(out);
  };
  unique_preserve(incorrect_vertices);
  unique_preserve(unmatched_vertices);

  static size_t debug_counter = 0;
  ++debug_counter;

  writeRadiusRingWalkFailDebugTxt(
    makeRadiusRingWalkFailDebugPath(
      kin_del, occurrence_time, runtime_branch_id, delaunay_face_id, strand_cell_id, debug_counter, ".txt"),
    occurrence_time, runtime_branch_id, delaunay_face_id, strand_cell_id, fail_reason, ring, incorrect_vertices,
    unmatched_vertices);

  // Common event-style kinetic SVG with the affected triangle / walked sites / edges / VVs highlighted.
  {
    const auto& graph = kin_del.getGraph();
    VisualDebugHighlight highlight;
    if (graph.isLiveFace(delaunay_face_id))
    {
      highlight.addDelaunayTriangle(graph, delaunay_face_id);
    }
    highlight.delaunay_vertices.insert(strand_cell_id);
    for (const RadiusRingWalkDebugVertex& vert : ring)
    {
      if (vert.strand_id.has_value())
      {
        highlight.delaunay_vertices.insert(vert.strand_id.value());
      }
      if (vert.voronoi_vertex_id.has_value())
      {
        highlight.voronoi_vertices.insert(vert.voronoi_vertex_id.value());
      }
      if (vert.delaunay_edge_id.has_value())
      {
        const size_t de = vert.delaunay_edge_id.value();
        highlight.label_crossings_on_delaunay_edges.insert(de);
        highlight.addUndirectedDelaunayEdge(graph, 2 * de);
      }
      if (vert.voronoi_edge_id.has_value())
      {
        const size_t ve = vert.voronoi_edge_id.value();
        highlight.voronoi_edges.insert(ve);
        highlight.label_crossings_on_voronoi_edges.insert(ve);
      }
      if (vert.delaunay_edge_id.has_value() && vert.voronoi_edge_id.has_value())
      {
        highlight.crossing_intersection_keys.insert(
          (static_cast<uint64_t>(vert.delaunay_edge_id.value()) << 32) | vert.voronoi_edge_id.value());
      }
    }
    const std::string event_descriptor = "radius_ring_walk_FAIL_face" + std::to_string(delaunay_face_id) + "_strand"
      + std::to_string(strand_cell_id);
    writeSegmentBuilderVisualDebugSvg(
      true, kin_del, graph, occurrence_time, "error", event_descriptor, highlight, runtime_branch_id);
  }

  const std::filesystem::path filepath = makeRadiusRingWalkFailDebugPath(
    kin_del, occurrence_time, runtime_branch_id, delaunay_face_id, strand_cell_id, debug_counter, ".svg");

  double min_x = std::numeric_limits<double>::infinity();
  double min_y = std::numeric_limits<double>::infinity();
  double max_x = -std::numeric_limits<double>::infinity();
  double max_y = -std::numeric_limits<double>::infinity();
  for (const RadiusRingWalkDebugVertex& vert : ring)
  {
    if (!std::isfinite(vert.xy.x) || !std::isfinite(vert.xy.y))
    {
      continue;
    }
    min_x = std::min(min_x, vert.xy.x);
    min_y = std::min(min_y, vert.xy.y);
    max_x = std::max(max_x, vert.xy.x);
    max_y = std::max(max_y, vert.xy.y);
  }
  if (!std::isfinite(min_x) || !std::isfinite(min_y) || !std::isfinite(max_x) || !std::isfinite(max_y))
  {
    min_x = -1.0;
    min_y = -1.0;
    max_x = 1.0;
    max_y = 1.0;
  }

  constexpr double pad = 0.05;
  const double span_x = std::max(max_x - min_x, 1e-6);
  const double span_y = std::max(max_y - min_y, 1e-6);
  min_x -= pad * span_x;
  min_y -= pad * span_y;
  max_x += pad * span_x;
  max_y += pad * span_y;
  const double width = max_x - min_x;
  const double height = max_y - min_y;
  auto svg_x = [&](double x) { return x - min_x; };
  auto svg_y = [&](double y) { return max_y - y; };

  std::ostringstream poly_points;
  for (size_t i = 0; i < ring.size(); ++i)
  {
    if (!std::isfinite(ring[i].xy.x) || !std::isfinite(ring[i].xy.y))
    {
      continue;
    }
    if (poly_points.tellp() > 0)
    {
      poly_points << ' ';
    }
    poly_points << svg_x(ring[i].xy.x) << ',' << svg_y(ring[i].xy.y);
  }

  std::ofstream out(filepath);
  if (!out)
  {
    KINDS_WARNING("radius ring walk FAIL: failed to open debug SVG " << filepath.generic_string()
      << " runtime_branch="
      << (runtime_branch_id.has_value() ? std::to_string(runtime_branch_id.value()) : "unresolved")
      << " strand_cell=" << strand_cell_id << " face=" << delaunay_face_id << " t=" << occurrence_time);
    return;
  }

  out << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n";
  out << "<svg xmlns=\"http://www.w3.org/2000/svg\" viewBox=\"0 0 " << width << ' ' << height
      << "\" width=\"" << width * 100.0 << "\" height=\"" << height * 100.0 << "\">\n";
  out << "<rect x=\"0\" y=\"0\" width=\"" << width << "\" height=\"" << height << "\" fill=\"#f8f8f8\"/>\n";
  if (poly_points.tellp() > 0)
  {
    out << "<polyline points=\"" << poly_points.str()
        << "\" fill=\"none\" stroke=\"#cc0000\" stroke-width=\"" << span_x * 0.002 << "\"/>\n";
  }

  for (size_t i = 0; i < ring.size(); ++i)
  {
    const RadiusRingWalkDebugVertex& vert = ring[i];
    if (!std::isfinite(vert.xy.x) || !std::isfinite(vert.xy.y))
    {
      continue;
    }
    const double cx = svg_x(vert.xy.x);
    const double cy = svg_y(vert.xy.y);
    const double r = span_x * 0.008;
    const char* fill = vert.incorrect ? "#ff0000" : "#ffaa00";
    out << "<circle cx=\"" << cx << "\" cy=\"" << cy << "\" r=\"" << r
        << "\" fill=\"" << fill << "\" stroke=\"#222\" stroke-width=\"" << span_x * 0.0008 << "\"/>\n";

    std::ostringstream label;
    label << "i=" << i << " " << vert.source << " (" << radiusDebugNumberLiteral(vert.xy.x) << ","
          << radiusDebugNumberLiteral(vert.xy.y) << ")";
    if (vert.voronoi_vertex_id.has_value())
    {
      label << " vv=" << vert.voronoi_vertex_id.value();
    }
    if (vert.strand_id.has_value())
    {
      label << " strand=" << vert.strand_id.value();
    }
    if (vert.delaunay_edge_id.has_value())
    {
      label << " de=" << vert.delaunay_edge_id.value();
    }
    if (vert.voronoi_edge_id.has_value())
    {
      label << " ve=" << vert.voronoi_edge_id.value();
    }
    if (vert.incorrect)
    {
      label << " INCORRECT";
    }
    out << "<text x=\"" << (cx + r * 1.5) << "\" y=\"" << (cy - r * 1.5) << "\" font-size=\"" << span_x * 0.02
        << "\" fill=\"#111\">" << radiusDebugSvgEscape(label.str()) << "</text>\n";
  }

  out << "<text x=\"" << span_x * 0.02 << "\" y=\"" << span_y * 0.05 << "\" font-size=\"" << span_x * 0.025
      << "\" fill=\"#111\">radius ring walk FAIL face=" << delaunay_face_id << " strand=" << strand_cell_id
      << " verts=" << ring.size() << "</text>\n";
  out << "<text x=\"" << span_x * 0.02 << "\" y=\"" << span_y * 0.09 << "\" font-size=\"" << span_x * 0.02
      << "\" fill=\"#444\">" << radiusDebugSvgEscape(fail_reason) << "</text>\n";
  out << "</svg>\n";

  KINDS_WARNING("radius ring walk FAIL: wrote debug SVG to " << filepath.generic_string()
    << " runtime_branch="
    << (runtime_branch_id.has_value() ? std::to_string(runtime_branch_id.value()) : "unresolved")
    << " strand_cell=" << strand_cell_id << " face=" << delaunay_face_id << " t=" << occurrence_time);
}

void trySetOwnedInteriorVoronoiVertexForRadiusShift(RadiusBoundaryTransitionShiftContext& ctx,
  const KineticDelaunay& kin_del, size_t triangle_face_id)
{
  const std::vector<size_t> contained = kin_del.getCrossingDataVoronoiVerticesInTri(triangle_face_id);
  if (contained.empty())
  {
    return;
  }

  // Prefer the triangle's own circumcenter when it is contained; otherwise any contained VV is enough to
  // reject shifting and take the traced triangle-cap path.
  if (kin_del.isCrossingDataVoronoiVertexRegistered(triangle_face_id)
    && kin_del.getCrossingDataContainingTriId(triangle_face_id) == triangle_face_id)
  {
    ctx.interior_voronoi_vertex_id = triangle_face_id;
    return;
  }
  ctx.interior_voronoi_vertex_id = contained.front();
}

bool triangleSpansFutureBranchesOfPendingSplit(const KineticDelaunay& kin_del, size_t triangle_face_id)
{
  const auto& graph = kin_del.getGraph();
  if (!graph.isLiveFace(triangle_face_id))
  {
    return false;
  }

  const std::array<int, 3> triangle_vertices = graph.getTriangleVertexIndices(triangle_face_id);
  std::array<size_t, 3> finite_vertices {};
  size_t finite_vertex_count = 0;
  for (int vertex : triangle_vertices)
  {
    if (vertex < 0 || kin_del.isDummyBoundary(static_cast<size_t>(vertex)))
    {
      return false;
    }
    finite_vertices[finite_vertex_count++] = static_cast<size_t>(vertex);
  }
  if (finite_vertex_count != 3)
  {
    return false;
  }

  bool spans_pending_split = false;
  kin_del.visitPendingBranchSplits(
    [&](size_t, const PendingBranchSplit& split)
    {
      if (spans_pending_split || split.frozen_parent_strands.empty())
      {
        return;
      }

      const auto belongs_to_split = [&](size_t strand_id)
      {
        return std::find(split.frozen_parent_strands.begin(), split.frozen_parent_strands.end(), strand_id)
          != split.frozen_parent_strands.end();
      };
      if (!belongs_to_split(finite_vertices[0]) || !belongs_to_split(finite_vertices[1])
        || !belongs_to_split(finite_vertices[2]))
      {
        return;
      }

      // notePendingBranchSplit assigns split-off strands their child runtime branch ids immediately,
      // before the graph cut. These are therefore the branch ids the vertices will retain after the split.
      const size_t branch0 = kin_del.getRuntimeBranchIdForStrand(finite_vertices[0]);
      const size_t branch1 = kin_del.getRuntimeBranchIdForStrand(finite_vertices[1]);
      const size_t branch2 = kin_del.getRuntimeBranchIdForStrand(finite_vertices[2]);
      spans_pending_split = branch0 != branch1 || branch0 != branch2;
    });
  return spans_pending_split;
}

void classifyRadiusBoundaryTransitionShiftContext(RadiusBoundaryTransitionShiftContext& ctx,
  size_t pre_boundary_edge_count, const std::array<bool, 3>& pre_is_boundary_edge,
  const std::array<size_t, 3>& pre_face_he, size_t post_boundary_edge_count,
  const std::array<bool, 3>& post_is_boundary_edge, const std::array<size_t, 3>& post_face_he,
  const KineticDelaunay& kin_del, size_t pre_face_id, size_t post_face_id, double t)
{
  if (pre_boundary_edge_count == 2 && post_boundary_edge_count == 1)
  {
    size_t out = 0;
    for (size_t i = 0; i < 3; ++i)
    {
      if (pre_is_boundary_edge[i])
      {
        ctx.source_delaunay_edges[out++] = pre_face_he[i] / 2;
      }
    }
    for (size_t i = 0; i < 3; ++i)
    {
      if (post_is_boundary_edge[i])
      {
        ctx.target_delaunay_edge = post_face_he[i] / 2;
        break;
      }
    }
    if (out == 2)
    {
      ctx.roles_valid = true;
      trySetOwnedInteriorVoronoiVertexForRadiusShift(ctx, kin_del, pre_face_id);
      KINDS_DEBUG("Radius boundary transition 2->1: source_delaunay_edges=("
                  << ctx.source_delaunay_edges[0] << "," << ctx.source_delaunay_edges[1] << ") target_delaunay_edge="
                  << ctx.target_delaunay_edge << " interior_voronoi_vertex="
                  << (ctx.interior_voronoi_vertex_id.has_value()
                       ? std::to_string(ctx.interior_voronoi_vertex_id.value())
                       : std::string("none"))
                  << " pre_face=" << pre_face_id << " post_face=" << post_face_id << " t=" << t);
    }
  }
  else if (pre_boundary_edge_count == 1 && post_boundary_edge_count == 2)
  {
    size_t out = 0;
    for (size_t i = 0; i < 3; ++i)
    {
      if (post_is_boundary_edge[i])
      {
        ctx.source_delaunay_edges[out++] = post_face_he[i] / 2;
      }
    }
    for (size_t i = 0; i < 3; ++i)
    {
      if (pre_is_boundary_edge[i])
      {
        ctx.target_delaunay_edge = pre_face_he[i] / 2;
        break;
      }
    }
    if (out == 2)
    {
      ctx.roles_valid = true;
      trySetOwnedInteriorVoronoiVertexForRadiusShift(ctx, kin_del, post_face_id);
      KINDS_DEBUG("Radius boundary transition 1->2: source_delaunay_edges=("
                  << ctx.source_delaunay_edges[0] << "," << ctx.source_delaunay_edges[1] << ") target_delaunay_edge="
                  << ctx.target_delaunay_edge << " interior_voronoi_vertex="
                  << (ctx.interior_voronoi_vertex_id.has_value()
                       ? std::to_string(ctx.interior_voronoi_vertex_id.value())
                       : std::string("none"))
                  << " pre_face=" << pre_face_id << " post_face=" << post_face_id << " t=" << t);
    }
  }
}
} // namespace

void SegmentBuilderRadiusCallback::beforeEvent(KineticDelaunay::Event& e)
{
  auto* radius = dynamic_cast<KineticDelaunay::RadiusEvent*>(&e);
  if (!radius)
  {
    return;
  }
  segment_builder_.clearRadiusShiftedSiteCache();
  radius_pre_finished_one_edge_meshlet_ids_.clear();
  SegmentBuilder::ScopedMetadataCallbackPhase callback_phase(segment_builder_, "before");
  auto& graph = segment_builder_.kin_del.getGraph();
  const auto radius_vertices = graph.adjacentTriangleVertices(radius->half_edge_id);
  const size_t runtime_branch_id
    = segment_builder_.kin_del.getRuntimeBranchIdForHalfEdge(radius->half_edge_id);

  size_t face_id = graph.halfEdge(radius->half_edge_id).face;
  bool is_inside = segment_builder_.kin_del.getFaceInside(face_id);
  // Label by the sign-change direction carried by the event, not the mutable pre-toggle state.
  const std::string radius_transition_tag = radius->target_inside ? "add" : "remove";

  writeSegmentBuilderVisualDebugSvg(segment_builder_.visual_debug, segment_builder_.kin_del, graph,
    radius->occurrence_time, "before",
    "radius_he" + std::to_string(radius->half_edge_id) + "_" + radius_transition_tag,
    VisualDebugHighlight::forRadius(graph, radius->half_edge_id), runtime_branch_id,
    /*separation_offset_segments=*/nullptr, /*seam_outlines=*/nullptr, /*explicit_runtime_branch_ids=*/nullptr,
    radius->creation_time);
  const auto& face_half_edges = graph.face(face_id).half_edges;
  const double t = radius->occurrence_time;

  // Intersection dParams are only refreshed on crossing/flip topology updates; recompute them at event time
  // in kinetic space before boundary-interval finishing uses stale ordering or shift interpolation.
  segment_builder_.kin_del.refreshTriangleDelaunayEdgeIntersectionParams(face_id, t);

  std::array<bool, 3> is_boundary_edge {};
  size_t boundary_edge_count = 0;
  for (size_t i = 0; i < 3; ++i)
  {
    const size_t he_id = face_half_edges[i];
    is_boundary_edge[i] = segment_builder_.kin_del.isOnComponentBoundary(he_id);
    if (is_boundary_edge[i])
    {
      ++boundary_edge_count;
    }
  }

  // Pre-flip 2-boundary case: finishing intersection strips on the two boundary edges before the flip; neighbor remap
  // uses the third (internal) Delaunay edge as target — same as post-flip 2->1 classification after topology updates.
  RadiusBoundaryTransitionShiftContext radius_pre_finish_shift_ctx {};
  const RadiusBoundaryTransitionShiftContext* radius_finish_shift_arg = nullptr;
  bool radius_event_will_create_mixed_branch_split = false;
  if (segment_builder_.radius_pending_split_triangle_fallback_enabled
    && segment_builder_.radius_boundary_transition_shift_enabled)
  {
    std::vector<bool> target_face_inside = segment_builder_.kin_del.getFacesInside();
    if (face_id < target_face_inside.size())
    {
      target_face_inside[face_id] = radius->target_inside;
      radius_event_will_create_mixed_branch_split
        = !segment_builder_.kin_del.checkForSplit(radius_vertices, target_face_inside, radius->occurrence_time)
             .empty();
    }
  }
  const bool pre_triangle_spans_pending_split
    = segment_builder_.radius_pending_split_triangle_fallback_enabled
    && segment_builder_.radius_boundary_transition_shift_enabled
    && (triangleSpansFutureBranchesOfPendingSplit(segment_builder_.kin_del, face_id)
      || radius_event_will_create_mixed_branch_split);
  if (segment_builder_.radius_boundary_transition_shift_enabled && boundary_edge_count == 2)
  {
    size_t out = 0;
    size_t internal_corner_i = static_cast<size_t>(-1);
    for (size_t i = 0; i < 3; ++i)
    {
      if (is_boundary_edge[i])
      {
        radius_pre_finish_shift_ctx.source_delaunay_edges[out++] = face_half_edges[i] / 2;
      }
      else
      {
        internal_corner_i = i;
      }
    }
    if (out == 2 && internal_corner_i < 3)
    {
      radius_pre_finish_shift_ctx.target_delaunay_edge = face_half_edges[internal_corner_i] / 2;
      radius_pre_finish_shift_ctx.roles_valid = true;
      trySetOwnedInteriorVoronoiVertexForRadiusShift(radius_pre_finish_shift_ctx, segment_builder_.kin_del, face_id);
      if (!pre_triangle_spans_pending_split
        && radiusBoundaryTransitionShiftApplicable(
          segment_builder_.radius_boundary_transition_shift_enabled, radius_pre_finish_shift_ctx))
      {
        radius_finish_shift_arg = &radius_pre_finish_shift_ctx;
      }
      KINDS_DEBUG("Radius beforeEvent: finishMeshFromIntersections shift context (pre-flip 2 boundary) sources=("
                  << radius_pre_finish_shift_ctx.source_delaunay_edges[0] << ","
                  << radius_pre_finish_shift_ctx.source_delaunay_edges[1] << ") target_internal_delaunay_edge="
                  << radius_pre_finish_shift_ctx.target_delaunay_edge
                  << " interior_voronoi_vertex="
                  << (radius_pre_finish_shift_ctx.interior_voronoi_vertex_id.has_value()
                       ? std::to_string(radius_pre_finish_shift_ctx.interior_voronoi_vertex_id.value())
                       : std::string("none"))
                  << " shift_applicable="
                  << (radius_finish_shift_arg != nullptr ? "true" : "false")
                  << " pending_split_mixed_future_branches="
                  << (pre_triangle_spans_pending_split ? "true" : "false")
                  << " radius_event_creates_split="
                  << (radius_event_will_create_mixed_branch_split ? "true" : "false")
                  << " face=" << face_id << " t=" << t);
    }
  }

  // Before the radius topology update, finish all active boundary-interval meshes on
  // boundary Delaunay edges of the affected triangle.
  std::unordered_set<size_t> processed_boundary_he_even;
  for (size_t he_id : face_half_edges)
  {
    const size_t he_even = he_id & ~1;
    if (!processed_boundary_he_even.insert(he_even).second)
    {
      continue;
    }
    if (!segment_builder_.kin_del.isOnComponentBoundary(he_even))
    {
      continue;
    }

    const size_t d_edge_id = he_even / 2;
    const auto& d_intersections = segment_builder_.kin_del.getCrossingData().delaunay_edge_intersections[d_edge_id];
    if (d_intersections.empty())
    {
      continue;
    }

    const std::vector<KineticDelaunay::CrossingData::EdgeIntersectionRef> refs
      = segment_builder_.getBoundaryIntersectionsInBoundaryOrder(d_edge_id);

    const bool buffer_one_edge_finishes = (boundary_edge_count == 1);
    auto finish_and_maybe_buffer
      = [&](size_t cell, std::optional<KineticDelaunay::CrossingData::EdgeIntersectionRef> start,
          std::optional<KineticDelaunay::CrossingData::EdgeIntersectionRef> end)
    {
      const size_t mesh_id
        = segment_builder_.resolveIntersectionMeshPairIndex(cell, start, end, t);
      segment_builder_.finishMeshFromIntersections(cell, t, start, end, SegmentBuilder::BoundaryEventType::Radius,
        SegmentBuilder::BoundarySegmentAction::SegmentCompleted, radius_finish_shift_arg, mesh_id);
      if (buffer_one_edge_finishes && mesh_id != static_cast<size_t>(-1))
      {
        radius_pre_finished_one_edge_meshlet_ids_.push_back(mesh_id);
      }
    };

    {
      const size_t first_cell
        = segment_builder_.determineVoronoiCellForBoundaryIntersectionInterval(d_edge_id, std::nullopt, refs.front());
      finish_and_maybe_buffer(first_cell, std::nullopt, refs.front());
    }
    for (size_t k = 0; k + 1 < refs.size(); ++k)
    {
      const size_t mid_cell
        = segment_builder_.determineVoronoiCellForBoundaryIntersectionInterval(d_edge_id, refs[k], refs[k + 1]);
      if (mid_cell == static_cast<size_t>(-1))
      {
        continue;
      }
      finish_and_maybe_buffer(mid_cell, refs[k], refs[k + 1]);
    }
    {
      const size_t last_cell
        = segment_builder_.determineVoronoiCellForBoundaryIntersectionInterval(d_edge_id, refs.back(), std::nullopt);
      finish_and_maybe_buffer(last_cell, refs.back(), std::nullopt);
    }
  }

  if (radius_finish_shift_arg != nullptr)
  {
    // Complementary target mid-interval is not finished here (internal for 2-boundary pre-flip).
    // Queue happens in afterEvent when the target mid is started, or via existing-mid lookup.
    const auto& graph_for_site = segment_builder_.kin_del.getGraph();
    const size_t s0 = radius_finish_shift_arg->source_delaunay_edges[0];
    const size_t s1 = radius_finish_shift_arg->source_delaunay_edges[1];
    std::optional<size_t> shared_site;
    const int a0 = graph_for_site.halfEdge(2 * s0).origin;
    const int a1 = graph_for_site.halfEdge(2 * s0 + 1).origin;
    const int b0 = graph_for_site.halfEdge(2 * s1).origin;
    const int b1 = graph_for_site.halfEdge(2 * s1 + 1).origin;
    for (int a : { a0, a1 })
    {
      if (a < 0)
      {
        continue;
      }
      if (a == b0 || a == b1)
      {
        shared_site = static_cast<size_t>(a);
        break;
      }
    }
    if (shared_site.has_value())
    {
      // Ensure the shifted site is computed once during finish of the two source one-nulls.
      (void)segment_builder_.getOrFillRadiusShiftedSiteCache(t, shared_site.value(), radius_finish_shift_arg);
    }
  }

  radius_pre_boundary_edge_count_ = boundary_edge_count;
  radius_pre_face_id_ = face_id;
  for (size_t i = 0; i < 3; ++i)
  {
    radius_pre_face_he_[i] = face_half_edges[i];
    radius_pre_is_boundary_edge_[i] = is_boundary_edge[i];
  }

  switch (boundary_edge_count)
  {
  case 0:
  {
    size_t vertices[3];
    for (size_t i = 0; i < 3; ++i)
    {
      vertices[i] = graph.halfEdge(face_half_edges[i]).origin;

      if (vertices[i] == size_t(-1))
      {
        KINDS_ERROR("Boundary triangle was turned (t=" << radius->occurrence_time
                                                      << "); that is impossible and will be ignored!");
        break;
      }
    }
    glm::dvec2 p0 = segment_builder_.kin_del.getPointAt(radius->occurrence_time, vertices[0], false, false);
    glm::dvec2 p1 = segment_builder_.kin_del.getPointAt(radius->occurrence_time, vertices[1], false, false);
    glm::dvec2 p2 = segment_builder_.kin_del.getPointAt(radius->occurrence_time, vertices[2], false, false);
    glm::dvec2 new_point = (p0 + p1 + p2) / 3.0;

    size_t new_vertex_index = segment_builder_.boundary_mesh.getVertices().size();
    segment_builder_.addBoundaryVertex(glm::dvec3 { new_point[0], new_point[1], radius->occurrence_time },
      glm::dvec2 { 0.0, 0.0 }, vertices[0], radius->occurrence_time, false);

    if (is_inside)
    {
      for (size_t i = 0; i < 3; ++i)
      {
        segment_builder_.boundary_mesh_last_left_and_right_vertex[face_half_edges[i]]
          = std::make_pair(new_vertex_index, new_vertex_index);
      }
    }
    else
    {
      for (size_t i = 0; i < 3; ++i)
      {
        size_t outer_he_id = face_half_edges[i] ^ 1;
        segment_builder_.boundary_mesh_last_left_and_right_vertex[outer_he_id]
          = std::make_pair(new_vertex_index, new_vertex_index);
      }
    }

    break;
  }
  case 1:
  {
    size_t boundary_he_index = 0;
    for (size_t i = 0; i < 3; ++i)
    {
      if (is_boundary_edge[i])
      {
        boundary_he_index = i;
        break;
      }
    }

    size_t inner_he_id = face_half_edges[boundary_he_index];
    size_t outer_he_id = inner_he_id ^ 1;
    size_t opposite_vertex = graph.triangleOppositeVertex(inner_he_id);

    glm::dvec2 opposite_point = segment_builder_.kin_del.getPointAt(radius->occurrence_time, opposite_vertex, false, false);
    size_t u = graph.halfEdge(inner_he_id).origin;
    glm::dvec2 p_u = segment_builder_.kin_del.getPointAt(radius->occurrence_time, u, false, false);
    size_t v = graph.halfEdge(outer_he_id).origin;
    glm::dvec2 p_v = segment_builder_.kin_del.getPointAt(radius->occurrence_time, v, false, false);

    glm::dvec2 new_boundary_vertex = (opposite_point + p_u + p_v) / 3.0;

    size_t component_id = segment_builder_.kin_del.component_data.component_map[v];
    auto& boundary_polygon = segment_builder_.kin_del.component_data.component_boundaries[component_id][0];
    auto centroid = polygonCentroid(boundary_polygon);

    size_t boundary_he_id = outer_he_id;
    if (!is_inside)
    {
      boundary_he_id = inner_he_id;
    }

    const auto& boundary_last_vertices = segment_builder_.boundary_mesh_last_left_and_right_vertex[boundary_he_id];

    size_t new_boundary_vertex_index = segment_builder_.addBoundaryVertex(
      glm::dvec3 { new_boundary_vertex[0], new_boundary_vertex[1], radius->occurrence_time }, centroid, opposite_vertex, radius->occurrence_time, false);

    size_t index
      = segment_builder_.addBoundaryTriangle(boundary_last_vertices.first, boundary_last_vertices.second, new_boundary_vertex_index);

    if (index == size_t(-1))
    {
      KINDS_DEBUG("\ninner_he_id: " << inner_he_id << "\nouter_he_id: " << outer_he_id
                                    << "\nlast_vertices: " << boundary_last_vertices.first << ", "
                                    << boundary_last_vertices.second << "\ninner last vertices: "
                                    << segment_builder_.boundary_mesh_last_left_and_right_vertex[inner_he_id].first << ", "
                                    << segment_builder_.boundary_mesh_last_left_and_right_vertex[inner_he_id].second << "\nu: " << u
                                    << ", v: " << v << ", opposite: " << opposite_vertex);
    }

    size_t he1_id = graph.halfEdge(inner_he_id).next;
    size_t he2_id = graph.halfEdge(he1_id).next;

    if (!is_inside)
    {
      he1_id = he1_id ^ 1;
      he2_id = he2_id ^ 1;
    }

    segment_builder_.boundary_mesh_last_left_and_right_vertex[he1_id]
      = std::make_pair(boundary_last_vertices.first, new_boundary_vertex_index);
    segment_builder_.boundary_mesh_last_left_and_right_vertex[he2_id]
      = std::make_pair(new_boundary_vertex_index, boundary_last_vertices.second);
    segment_builder_.boundary_mesh_last_left_and_right_vertex[boundary_he_id] = std::make_pair(-1, -1);

    break;
  }
  case 2:
  {
    size_t non_boundary_he_index = 0;
    for (size_t i = 0; i < 3; ++i)
    {
      if (!is_boundary_edge[i])
      {
        non_boundary_he_index = i;
        break;
      }
    }

    size_t inner_he_id = face_half_edges[non_boundary_he_index];
    size_t outer_he_id = inner_he_id ^ 1;

    size_t opposite_vertex = graph.triangleOppositeVertex(inner_he_id);

    glm::dvec2 opposite_point = segment_builder_.kin_del.getPointAt(radius->occurrence_time, opposite_vertex, false, false);
    size_t u = graph.halfEdge(inner_he_id).origin;
    glm::dvec2 p_u = segment_builder_.kin_del.getPointAt(radius->occurrence_time, u, false, false);
    size_t v = graph.halfEdge(outer_he_id).origin;
    glm::dvec2 p_v = segment_builder_.kin_del.getPointAt(radius->occurrence_time, v, false, false);
    glm::dvec2 old_boundary_vertex = (opposite_point + p_u + p_v) / 3.0;

    size_t component_id = segment_builder_.kin_del.component_data.component_map[v];
    auto& boundary_polygon = segment_builder_.kin_del.component_data.component_boundaries[component_id][0];
    auto centroid = polygonCentroid(boundary_polygon);

    size_t old_boundary_vertex_index = segment_builder_.boundary_mesh.getVertices().size();
    segment_builder_.addBoundaryVertex(
      glm::dvec3 { old_boundary_vertex[0], old_boundary_vertex[1], radius->occurrence_time }, centroid, opposite_vertex, radius->occurrence_time, false);

    size_t he1_id = graph.halfEdge(inner_he_id).next;
    size_t he2_id = graph.halfEdge(he1_id).next;
    if (is_inside)
    {
      he1_id = he1_id ^ 1;
      he2_id = he2_id ^ 1;
    }

    size_t index = segment_builder_.addBoundaryTriangle(segment_builder_.boundary_mesh_last_left_and_right_vertex[he1_id].first,
      segment_builder_.boundary_mesh_last_left_and_right_vertex[he1_id].second, old_boundary_vertex_index);
    if (index == size_t(-1))
    {
      KINDS_DEBUG("he1_id: " << he1_id << "\ntwin: " << (he1_id ^ 1)
                             << "\nlast_vertices: " << segment_builder_.boundary_mesh_last_left_and_right_vertex[he1_id].first << ", "
                             << segment_builder_.boundary_mesh_last_left_and_right_vertex[he1_id].second
                             << "\ntwin last vertices: " << segment_builder_.boundary_mesh_last_left_and_right_vertex[he1_id ^ 1].first
                             << ", " << segment_builder_.boundary_mesh_last_left_and_right_vertex[he1_id ^ 1].second
                             << "\nopposite: " << opposite_vertex);
    }

    index = segment_builder_.addBoundaryTriangle(segment_builder_.boundary_mesh_last_left_and_right_vertex[he2_id].first,
      segment_builder_.boundary_mesh_last_left_and_right_vertex[he2_id].second, old_boundary_vertex_index);
    if (index == size_t(-1))
    {
      KINDS_DEBUG("he2_id: " << he2_id << "\ntwin: " << (he2_id ^ 1)
                             << "\nlast_vertices: " << segment_builder_.boundary_mesh_last_left_and_right_vertex[he2_id].first << ", "
                             << segment_builder_.boundary_mesh_last_left_and_right_vertex[he2_id].second
                             << "\ntwin last vertices: " << segment_builder_.boundary_mesh_last_left_and_right_vertex[he2_id ^ 1].first
                             << ", " << segment_builder_.boundary_mesh_last_left_and_right_vertex[he2_id ^ 1].second
                             << "\nopposite: " << opposite_vertex);
    }

    segment_builder_.half_edge_to_boundary_vertex_index[outer_he_id] = old_boundary_vertex_index;

    if (!is_inside)
    {
      segment_builder_.boundary_mesh_last_left_and_right_vertex[outer_he_id]
        = std::make_pair(segment_builder_.boundary_mesh_last_left_and_right_vertex[he1_id].first,
          segment_builder_.boundary_mesh_last_left_and_right_vertex[he2_id].second);
    }
    else
    {
      segment_builder_.boundary_mesh_last_left_and_right_vertex[inner_he_id]
        = std::make_pair(segment_builder_.boundary_mesh_last_left_and_right_vertex[he2_id].second,
          segment_builder_.boundary_mesh_last_left_and_right_vertex[he1_id].first);
    }

    segment_builder_.boundary_mesh_last_left_and_right_vertex[he1_id] = std::make_pair(-1, -1);
    segment_builder_.boundary_mesh_last_left_and_right_vertex[he2_id] = std::make_pair(-1, -1);
    break;
  }
  case 3:
  {
    size_t vertices[3];
    for (size_t i = 0; i < 3; ++i)
    {
      vertices[i] = graph.halfEdge(face_half_edges[i]).origin;
    }
    glm::dvec2 p0 = segment_builder_.kin_del.getPointAt(radius->occurrence_time, vertices[0], false, false);
    glm::dvec2 p1 = segment_builder_.kin_del.getPointAt(radius->occurrence_time, vertices[1], false, false);
    glm::dvec2 p2 = segment_builder_.kin_del.getPointAt(radius->occurrence_time, vertices[2], false, false);
    glm::dvec2 new_point = (p0 + p1 + p2) / 3.0;

    size_t new_vertex_index = segment_builder_.addBoundaryVertex(glm::dvec3 { new_point[0], new_point[1], radius->occurrence_time },
      glm::dvec2 { 0.0, 0.0 }, vertices[0], radius->occurrence_time, false);

    for (size_t i = 0; i < 3; ++i)
    {
      const auto& last_vertices = segment_builder_.boundary_mesh_last_left_and_right_vertex[face_half_edges[i]];
      segment_builder_.addBoundaryTriangle(last_vertices.first, last_vertices.second, new_vertex_index);
    }

    break;
  }
  }

  // Wedge caps on boundary-interval strip meshes (intersection meshlets) at triangle corners — not `boundary_mesh`.
  // Triangle Delaunay edges are handled by the switch above; skip those. Resolve strip via one-null interval on the
  // incident boundary edge (even vs odd Delaunay origin), same pairing as `resolveIntersectionMeshPairIndex`.
  //
  // `mesh_start` / `mesh_end` follow interval order (see `startNewMeshFromIntersections`): for `[null, ref]` the open
  // site is interval *start* → update `mesh_start_vertex_id`; for `[ref, null]` the open site is interval *end* →
  // update `mesh_end_vertex_id`. That matches even-origin vs odd-origin from `determineVoronoiCellForBoundaryIntersectionInterval`.
  const auto& crossing_data = segment_builder_.kin_del.getCrossingData();
  std::unordered_set<size_t> radius_tri_edge_evens;
  for (size_t i = 0; i < 3; ++i)
  {
    radius_tri_edge_evens.insert(face_half_edges[i] & ~static_cast<size_t>(1));
  }
  std::unordered_set<size_t> processed_extra_boundary_edge_evens;
  for (size_t ti = 0; ti < 3; ++ti)
  {
    const int corner = graph.halfEdge(face_half_edges[ti]).origin;
    if (corner < 0)
    {
      continue;
    }
    const size_t corner_u = static_cast<size_t>(corner);
    const size_t corner_component = segment_builder_.kin_del.component_data.component_map[corner_u];

    for (auto it = graph.incidentEdgesBegin(corner_u); it != graph.incidentEdgesEnd(corner_u); ++it)
    {
      const size_t inc_he = *it;
      const size_t inc_even = inc_he & ~static_cast<size_t>(1);
      if (radius_tri_edge_evens.count(inc_even) != 0u)
      {
        continue;
      }
      if (!segment_builder_.kin_del.isOnComponentBoundary(inc_even))
      {
        continue;
      }
      if (!processed_extra_boundary_edge_evens.insert(inc_even).second)
      {
        continue;
      }

      const size_t d_edge_id = inc_even >> 1;
      const auto& d_inters = crossing_data.delaunay_edge_intersections[d_edge_id];
      if (d_inters.empty())
      {
        continue;
      }

      const std::vector<KineticDelaunay::CrossingData::EdgeIntersectionRef> refs
        = segment_builder_.getBoundaryIntersectionsInBoundaryOrder(d_edge_id);
      if (refs.empty())
      {
        continue;
      }

      const size_t he_even_corner = d_edge_id * 2u;
      const size_t he_odd_corner = he_even_corner | 1u;
      const int even_origin = graph.halfEdge(he_even_corner).origin;
      const int odd_origin = graph.halfEdge(he_odd_corner).origin;

      size_t pair_idx = static_cast<size_t>(-1);
      bool update_start_endpoint = false;

      if (corner == even_origin && refs.front()->prev_segment_mesh_pair_index != static_cast<size_t>(-1))
      {
        if (refs.front()->delaunay_ref.has_value() && *refs.front()->delaunay_ref == d_inters.begin())
        {
          pair_idx = refs.front()->prev_segment_mesh_pair_index;
          // `[null, first]` — open site is interval start (left / mesh_start).
          update_start_endpoint = true;
        }
      }
      else if (corner == odd_origin && refs.back()->next_segment_mesh_pair_index != static_cast<size_t>(-1))
      {
        const auto d_last_ref = refs.back()->delaunay_ref;
        if (d_last_ref.has_value() && *d_last_ref != d_inters.end() && std::next(*d_last_ref) == d_inters.end())
        {
          pair_idx = refs.back()->next_segment_mesh_pair_index;
          // `[last, null]` — open site is interval end (right / mesh_end).
          update_start_endpoint = false;
        }
      }

      if (pair_idx != static_cast<size_t>(-1))
      {
        pair_idx = segment_builder_.preferLiveBoundaryMeshPair(pair_idx);
      }

      if (pair_idx == static_cast<size_t>(-1) || pair_idx >= segment_builder_.intersection_meshes.size()
        || pair_idx >= segment_builder_.intersection_mesh_pair_last_left_and_right_vertex.size())
      {
        continue;
      }

      auto& segs = segment_builder_.intersection_mesh_pair_last_left_and_right_vertex[pair_idx];
      if (segs.empty())
      {
        continue;
      }
      auto& seg = segs.front();
      if (seg.mesh_start_vertex_id < 0 || seg.mesh_end_vertex_id < 0)
      {
        continue;
      }

      const bool boundary_even_he_is_outside = segment_builder_.kin_del.isOnComponentBoundaryOutside(he_even_corner);
      const int inside_boundary_he_id
        = boundary_even_he_is_outside ? static_cast<int>(he_odd_corner) : static_cast<int>(he_even_corner);

      std::vector<bool> he_visited(graph.halfEdgeSlotCount(), false);
      segment_builder_.updateBoundary(t, he_visited, corner_component);
      auto& boundary_polygon = segment_builder_.kin_del.component_data.component_boundaries[corner_component][0];
      const auto centroid = polygonCentroid(boundary_polygon);

      auto& mesh = segment_builder_.intersection_meshes[pair_idx];
      // Mesh site placement must never use getPointAt / radiusTransitionInterpolatedSitePosition XY:
      // those are profile/local (and the latter may return a synthetic blend or VV circumcenter).
      // Unshifted Delaunay-frame site comes only from computeMeshSiteVertexPosition.
      // If the corner is the radius-transition junction with an interior VV, snap to that VV instead.
      glm::dvec3 corner_pos
        = segment_builder_.computeMeshSiteVertexPosition(glm::dvec3(0.0, 0.0, t), corner_u, t);
      std::optional<size_t> corner_snap_voronoi_vertex_id;
      if (radius_finish_shift_arg != nullptr)
      {
        if (auto site_shifted = segment_builder_.radiusTransitionInterpolatedSitePosition(
              t, corner_u, d_edge_id, radius_finish_shift_arg);
          site_shifted.has_value() && site_shifted->snap_voronoi_vertex_id.has_value())
        {
          corner_snap_voronoi_vertex_id = site_shifted->snap_voronoi_vertex_id;
          // Intentionally ignore site_shifted->position (local VV / blend). addMeshletVertex places
          // via barycentric transfer when snap_voronoi_vertex_id is set.
        }
      }
      const std::string radius_corner_meta = segment_builder_.composeBoundaryMetadata(
        SegmentBuilder::BoundaryEventType::Radius, SegmentBuilder::BoundarySegmentAction::SegmentRemapped);
      // Placement is a site corner, optionally snapped to the interior VV of a radius boundary-transition shift.
      const char* corner_source
        = corner_snap_voronoi_vertex_id.has_value() ? "Voronoi vertex" : "site";
      const std::string vertex_meta = segment_builder_.store_mesh_metadata
        ? SegmentBuilder::MetadataBuilder::fromObject(radius_corner_meta)
            .addString("source", corner_source)
            .addString("pos", update_start_endpoint ? "left" : "right")
            .build()
        : std::string {};
      const glm::dvec3 vertex_color
        = update_start_endpoint ? glm::dvec3(1.0, 0.0, 0.0) : glm::dvec3(0.0, 0.0, 1.0);
      const size_t eff_l = segment_builder_.intersectionStripEffectiveVertexIndex(seg, true);
      const size_t eff_r = segment_builder_.intersectionStripEffectiveVertexIndex(seg, false);
      const size_t new_vid = segment_builder_.addMeshletVertex(
        mesh, boundary_polygon, centroid, corner_pos, corner_u, t, false, corner_snap_voronoi_vertex_id, vertex_meta,
        vertex_color);
      segment_builder_.addBoundaryIntervalTriangleOriented(
        mesh, eff_l, eff_r, new_vid, inside_boundary_he_id, t, radius_corner_meta, pair_idx);
      segment_builder_.applyIntersectionStripOneSidedFixedVertex(mesh, seg, update_start_endpoint, new_vid,
        inside_boundary_he_id, std::nullopt, boundary_polygon, centroid, corner_u, t, true, true);
    }
  }
}

void SegmentBuilderRadiusCallback::afterEvent(KineticDelaunay::Event& e)
{
  auto* radius = dynamic_cast<KineticDelaunay::RadiusEvent*>(&e);
  if (!radius)
  {
    return;
  }
  segment_builder_.clearRadiusShiftedSiteCache();
  SegmentBuilder::ScopedMetadataCallbackPhase callback_phase(segment_builder_, "after");
  if (segment_builder_.diagnostics)
  {
    segment_builder_.logDiagnosticsMonitoredFaceInsideState(radius->occurrence_time, "radius_event");
  }
  auto& graph = segment_builder_.kin_del.getGraph();
  auto vertices = graph.adjacentTriangleVertices(radius->half_edge_id);

  // Outside→inside can glue previously separate kinetic components across this triangle. Fold them into the
  // smallest component id before any split bookkeeping runs on the (now possibly merged) parent.
  if (radius->target_inside)
  {
    segment_builder_.kin_del.consolidateComponentsAtTriangle(vertices, radius->occurrence_time);
  }

  size_t component_id = segment_builder_.kin_del.component_data.component_map[vertices[0]];
  auto split = segment_builder_.kin_del.checkForSplit(vertices, radius->occurrence_time);
  segment_builder_.splitComponent(component_id, split, radius->occurrence_time);

  std::unordered_set<size_t> components_to_refresh;
  for (int vertex : vertices)
  {
    if (vertex < 0)
    {
      continue;
    }
    const size_t strand_id = static_cast<size_t>(vertex);
    if (strand_id < segment_builder_.kin_del.component_data.component_map.size())
    {
      components_to_refresh.insert(segment_builder_.kin_del.component_data.component_map[strand_id]);
    }
  }

  std::vector<bool> visited(segment_builder_.kin_del.getGraph().halfEdgeSlotCount(), false);
  for (size_t refresh_component_id : components_to_refresh)
  {
    if (refresh_component_id >= segment_builder_.kin_del.component_data.components.size()
      || segment_builder_.kin_del.component_data.components[refresh_component_id].empty())
    {
      continue;
    }
    segment_builder_.kin_del.component_data.component_boundaries[refresh_component_id]
      = segment_builder_.kin_del.extractComponentBoundaries(
        segment_builder_.kin_del.component_data.components[refresh_component_id], radius->occurrence_time, visited,
        false);
    if (!segment_builder_.kin_del.component_data.component_boundaries[refresh_component_id].empty()
      && !segment_builder_.kin_del.component_data.component_boundaries[refresh_component_id][0].empty())
    {
      segment_builder_.kin_del.component_data.component_centroids[refresh_component_id] = polygonCentroid(
        segment_builder_.kin_del.component_data.component_boundaries[refresh_component_id][0]);
    }
    segment_builder_.kin_del.component_data.component_last_updated[refresh_component_id] = radius->occurrence_time;
  }

  // Keep the same sign-change direction tag as beforeEvent.
  const std::string radius_transition_tag = radius->target_inside ? "add" : "remove";

  writeSegmentBuilderVisualDebugSvg(segment_builder_.visual_debug, segment_builder_.kin_del, graph,
    radius->occurrence_time, "after",
    "radius_he" + std::to_string(radius->half_edge_id) + "_" + radius_transition_tag,
    VisualDebugHighlight::forRadius(graph, radius->half_edge_id),
    segment_builder_.kin_del.getRuntimeBranchIdForHalfEdge(radius->half_edge_id),
    /*separation_offset_segments=*/nullptr, /*seam_outlines=*/nullptr, /*explicit_runtime_branch_ids=*/nullptr,
    radius->creation_time);

  auto triangle_he_ids = graph.getTriangleHalfEdgeIndices(radius->half_edge_id);
  std::unordered_set<size_t> affected_delaunay_edges;
  for (size_t tri_he_id : triangle_he_ids)
  {
    affected_delaunay_edges.insert(tri_he_id / 2);
  }

  const auto& crossing_data = segment_builder_.kin_del.getCrossingData();
  const size_t affected_face_id = graph.halfEdge(radius->half_edge_id).face;
  const EventTime occurrence = radius->occurrence_time;
  const double t = occurrence.real_time;
  const size_t runtime_branch_id = segment_builder_.kin_del.getRuntimeBranchIdForHalfEdge(radius->half_edge_id);
  const bool new_inside_state = segment_builder_.kin_del.getFaceInside(affected_face_id);
  const bool inside_to_outside = !new_inside_state; // post-flip outside <=> inside->outside transition
  const auto affected_face_he = graph.face(affected_face_id).half_edges;

  std::array<bool, 3> post_is_boundary_edge {};
  size_t post_boundary_edge_count = 0;
  for (size_t i = 0; i < 3; ++i)
  {
    post_is_boundary_edge[i] = segment_builder_.kin_del.isOnComponentBoundary(affected_face_he[i]);
    if (post_is_boundary_edge[i])
    {
      ++post_boundary_edge_count;
    }
  }

  segment_builder_.invalidateStaleAlphaBoundaryMeshLinksOnTriangleEdgesLeftBoundary(affected_face_he);

  RadiusBoundaryTransitionShiftContext radius_boundary_shift_ctx {};
  if (segment_builder_.radius_boundary_transition_shift_enabled)
  {
    classifyRadiusBoundaryTransitionShiftContext(radius_boundary_shift_ctx, radius_pre_boundary_edge_count_,
      radius_pre_is_boundary_edge_, radius_pre_face_he_, post_boundary_edge_count, post_is_boundary_edge,
      affected_face_he, segment_builder_.kin_del, radius_pre_face_id_, affected_face_id, t);
  }

  const bool triangle_spans_pending_split
    = segment_builder_.radius_pending_split_triangle_fallback_enabled
    && triangleSpansFutureBranchesOfPendingSplit(segment_builder_.kin_del, affected_face_id);
  const bool use_radius_boundary_shift
    = !triangle_spans_pending_split
    && radiusBoundaryTransitionShiftApplicable(
      segment_builder_.radius_boundary_transition_shift_enabled, radius_boundary_shift_ctx);
  const RadiusBoundaryTransitionShiftContext* radius_boundary_shift_arg
    = use_radius_boundary_shift ? &radius_boundary_shift_ctx : nullptr;
  // only_adjacent_segment filters strips by crossing topology; requires shift context for neighbor remap.
  const bool radius_only_adjacent_strips = use_radius_boundary_shift;

  if (segment_builder_.radius_boundary_transition_shift_enabled && !use_radius_boundary_shift)
  {
    KINDS_DEBUG("Radius: falling back to traced Voronoi-cell / triangle-cap meshlets (shift preconditions unmet:"
                << " roles_valid=" << (radius_boundary_shift_ctx.roles_valid ? "true" : "false")
                << " interior_voronoi_vertex="
                << (radius_boundary_shift_ctx.interior_voronoi_vertex_id.has_value()
                     ? std::to_string(radius_boundary_shift_ctx.interior_voronoi_vertex_id.value())
                     : std::string("none"))
                << " pre_boundary_edges=" << radius_pre_boundary_edge_count_
                << " post_boundary_edges=" << post_boundary_edge_count
                << " pending_split_mixed_future_branches="
                << (triangle_spans_pending_split ? "true" : "false") << " face=" << affected_face_id << " t=" << t);
  }

  auto edge_endpoints = [&](size_t d_edge_id) -> std::array<int, 2>
  {
    const size_t he_even = 2 * d_edge_id;
    return { graph.halfEdge(he_even).origin, graph.destination(he_even) };
  };

  auto voronoi_vertex_on_edge = [&](size_t d_edge_id, bool even_side) -> size_t
  {
    const size_t he = 2 * d_edge_id + (even_side ? 0 : 1);
    return static_cast<size_t>(graph.halfEdge(he).face);
  };

  auto intersection_position = [&](KineticDelaunay::CrossingData::EdgeIntersectionRef ref) -> glm::dvec3
  {
    // Trace/triangulate in Delaunay space; addMeshletVertex converts to mesh space for stored vertices.
    return getCrossingCoordsInDelaunaySpace(segment_builder_.kin_del, ref, t);
  };

  struct RadiusTraceVertex
  {
    glm::dvec3 position {};
    std::optional<KineticDelaunay::CrossingData::EdgeIntersectionRef> intersection {};
    std::optional<size_t> voronoi_vertex_id {};
    std::optional<size_t> site_strand_id {};
  };

  auto make_intersection_trace_vertex = [&](KineticDelaunay::CrossingData::EdgeIntersectionRef ref) -> RadiusTraceVertex
  {
    return RadiusTraceVertex { intersection_position(ref), ref, std::nullopt, std::nullopt };
  };
  auto make_voronoi_trace_vertex = [&](size_t vv, const glm::dvec3& pos) -> RadiusTraceVertex
  {
    return RadiusTraceVertex { pos, std::nullopt, vv, std::nullopt };
  };
  auto make_site_trace_vertex = [&](size_t site_id, const glm::dvec3& pos) -> RadiusTraceVertex
  {
    return RadiusTraceVertex { pos, std::nullopt, std::nullopt, site_id };
  };

  auto dual_triangle_edges = [&](size_t voronoi_vertex_id) -> std::array<size_t, 3>
  {
    const auto& tri = graph.face(voronoi_vertex_id);
    return { tri.half_edges[0] / 2, tri.half_edges[1] / 2, tri.half_edges[2] / 2 };
  };

  auto choose_voronoi_vertex_towards_inside = [&](size_t d_edge_id, std::optional<size_t> prev_vv) -> std::optional<size_t>
  {
    const size_t vv0 = voronoi_vertex_on_edge(d_edge_id, true);
    const size_t vv1 = voronoi_vertex_on_edge(d_edge_id, false);
    const auto vv_inside = [&](size_t vv) -> bool
    {
      return crossing_data.isVoronoiVertexRegistered(vv)
        && crossing_data.getContainingTriId(vv) == affected_face_id;
    };
    const bool vv0_inside = vv_inside(vv0);
    const bool vv1_inside = vv_inside(vv1);
    if (!vv0_inside && !vv1_inside)
    {
      return std::nullopt;
    }
    if (vv0_inside && !vv1_inside)
    {
      return vv0;
    }
    if (vv1_inside && !vv0_inside)
    {
      return vv1;
    }
    if (prev_vv.has_value())
    {
      if (vv0 != prev_vv.value())
      {
        return vv0;
      }
      if (vv1 != prev_vv.value())
      {
        return vv1;
      }
    }
    return vv0;
  };

  auto next_dual_edge_for_cell = [&](size_t vv, size_t current_d_edge, size_t cell_id) -> std::optional<size_t>
  {
    const auto tri_edges = dual_triangle_edges(vv);
    for (size_t candidate_d : tri_edges)
    {
      if (candidate_d == current_d_edge)
      {
        continue;
      }
      const auto ep = edge_endpoints(candidate_d);
      if (ep[0] >= 0 && static_cast<size_t>(ep[0]) == cell_id)
      {
        return candidate_d;
      }
      if (ep[1] >= 0 && static_cast<size_t>(ep[1]) == cell_id)
      {
        return candidate_d;
      }
    }
    return std::nullopt;
  };

  auto first_triangle_intersection_on_voronoi_edge
    = [&](size_t d_edge_id, size_t from_vv,
        std::optional<KineticDelaunay::CrossingData::EdgeIntersectionRef> exclude_ref)
    -> std::optional<KineticDelaunay::CrossingData::EdgeIntersectionRef> {
    const auto& v_list = crossing_data.voronoi_edge_intersections[d_edge_id];
    const size_t vv_even = voronoi_vertex_on_edge(d_edge_id, true);
    const bool forward = (from_vv == vv_even);
    if (forward)
    {
      for (auto it = v_list.begin(); it != v_list.end(); ++it)
      {
        auto ref = *it;
        if (exclude_ref.has_value() && ref == exclude_ref.value())
        {
          continue;
        }
        if (affected_delaunay_edges.find(ref->delaunay_edge_id) != affected_delaunay_edges.end())
        {
          return ref;
        }
      }
      return std::nullopt;
    }
    for (auto it = v_list.rbegin(); it != v_list.rend(); ++it)
    {
      auto ref = *it;
      if (exclude_ref.has_value() && ref == exclude_ref.value())
      {
        continue;
      }
      if (affected_delaunay_edges.find(ref->delaunay_edge_id) != affected_delaunay_edges.end())
      {
        return ref;
      }
    }
    return std::nullopt;
  };

  auto other_triangle_intersection_on_same_voronoi_edge
    = [&](size_t d_edge_id, KineticDelaunay::CrossingData::EdgeIntersectionRef exclude_ref)
    -> std::optional<KineticDelaunay::CrossingData::EdgeIntersectionRef> {
    const auto& v_list = crossing_data.voronoi_edge_intersections[d_edge_id];
    for (auto ref : v_list)
    {
      if (ref == exclude_ref)
      {
        continue;
      }
      if (affected_delaunay_edges.find(ref->delaunay_edge_id) != affected_delaunay_edges.end())
      {
        return ref;
      }
    }
    return std::nullopt;
  };

  auto edge_index_in_affected_triangle = [&](size_t delaunay_edge_id) -> std::optional<size_t>
  {
    for (size_t i = 0; i < 3; ++i)
    {
      if ((affected_face_he[i] / 2) == delaunay_edge_id)
      {
        return i;
      }
    }
    return std::nullopt;
  };

  auto triangle_boundary_chain_vertices = [&](size_t from_edge, size_t to_edge) -> std::vector<size_t>
  {
    std::vector<size_t> out;
    const auto from_idx_opt = edge_index_in_affected_triangle(from_edge);
    const auto to_idx_opt = edge_index_in_affected_triangle(to_edge);
    if (!from_idx_opt.has_value() || !to_idx_opt.has_value())
    {
      return out;
    }

    const size_t from_idx = from_idx_opt.value();
    const size_t to_idx = to_idx_opt.value();
    size_t idx = from_idx;
    while (idx != to_idx)
    {
      const size_t next_he = affected_face_he[idx];
      const size_t v = graph.destination(next_he);
      if (v != size_t(-1))
      {
        out.push_back(v);
      }
      idx = (idx + 1) % 3;
    }
    return out;
  };

  auto choose_triangle_boundary_chain = [&](size_t from_edge, size_t to_edge, size_t cell_id) -> std::vector<size_t>
  {
    const auto cw = triangle_boundary_chain_vertices(from_edge, to_edge);
    const auto ccw = triangle_boundary_chain_vertices(to_edge, from_edge);

    auto score = [&](const std::vector<size_t>& path) -> std::pair<size_t, size_t>
    {
      size_t contains_cell = 0;
      for (size_t v : path)
      {
        if (v == cell_id)
        {
          contains_cell = 1;
          break;
        }
      }
      // Prefer chain containing the cell site; tie-break by fewer corner vertices.
      return { contains_cell, std::numeric_limits<size_t>::max() - path.size() };
    };

    const auto s_cw = score(cw);
    const auto s_ccw = score(ccw);
    if (s_ccw > s_cw)
    {
      return ccw;
    }
    return cw;
  };

  std::vector<KineticDelaunay::CrossingData::EdgeIntersectionRef> seed_intersections;
  seed_intersections.reserve(32);
  for (size_t d_edge_id : affected_delaunay_edges)
  {
    const auto& d_list = crossing_data.delaunay_edge_intersections[d_edge_id];
    for (auto ref : d_list)
    {
      seed_intersections.push_back(ref);
    }
  }

  struct IntersectionKey
  {
    size_t voronoi_edge_id = 0;
    size_t delaunay_edge_id = 0;
    uint64_t delaunay_param_bits = 0;

    bool operator==(const IntersectionKey& other) const noexcept
    {
      return voronoi_edge_id == other.voronoi_edge_id && delaunay_edge_id == other.delaunay_edge_id
        && delaunay_param_bits == other.delaunay_param_bits;
    }
  };
  struct IntersectionKeyHash
  {
    size_t operator()(const IntersectionKey& k) const noexcept
    {
      size_t h = std::hash<size_t> {}(k.voronoi_edge_id);
      h ^= std::hash<size_t> {}(k.delaunay_edge_id) + 0x9e3779b97f4a7c15ull + (h << 6) + (h >> 2);
      h ^= std::hash<uint64_t> {}(k.delaunay_param_bits) + 0x9e3779b97f4a7c15ull + (h << 6) + (h >> 2);
      return h;
    }
  };
  auto ref_key = [](KineticDelaunay::CrossingData::EdgeIntersectionRef ref) -> IntersectionKey
  {
    IntersectionKey k;
    k.voronoi_edge_id = ref->voronoi_edge_id;
    k.delaunay_edge_id = ref->delaunay_edge_id;
    std::memcpy(&k.delaunay_param_bits, &ref->delaunay_edge_param, sizeof(uint64_t));
    return k;
  };

  auto edge_touches_cell = [&](size_t voronoi_edge_id, size_t cell_id) -> bool
  {
    const auto ep = edge_endpoints(voronoi_edge_id);
    if (ep[0] >= 0 && static_cast<size_t>(ep[0]) == cell_id)
    {
      return true;
    }
    if (ep[1] >= 0 && static_cast<size_t>(ep[1]) == cell_id)
    {
      return true;
    }
    return false;
  };

  auto phase2_intersection_on_triangle_edge
    = [&](size_t he_id, bool traverse_he_forward,
        std::optional<KineticDelaunay::CrossingData::EdgeIntersectionRef> start_ref)
    -> std::optional<KineticDelaunay::CrossingData::EdgeIntersectionRef> {
    const size_t d_edge_id = he_id / 2;
    const auto& d_list = crossing_data.delaunay_edge_intersections[d_edge_id];
    if (d_list.empty())
    {
      return std::nullopt;
    }

    // Delaunay lists are sorted by parameter along the even half-edge direction.
    const bool traverse_even_forward = (he_id % 2 == 0) ? traverse_he_forward : !traverse_he_forward;

    if (!start_ref.has_value() || start_ref.value()->delaunay_edge_id != d_edge_id)
    {
      return traverse_even_forward ? d_list.front() : d_list.back();
    }

    std::vector<KineticDelaunay::CrossingData::EdgeIntersectionRef> refs;
    refs.reserve(d_list.size());
    for (auto ref : d_list)
    {
      refs.push_back(ref);
    }
    size_t start_i = size_t(-1);
    for (size_t i = 0; i < refs.size(); ++i)
    {
      if (refs[i] == start_ref.value())
      {
        start_i = i;
        break;
      }
    }
    if (start_i == size_t(-1))
    {
      return traverse_even_forward ? refs.front() : refs.back();
    }
    if (traverse_even_forward)
    {
      if (start_i + 1 < refs.size())
      {
        return refs[start_i + 1];
      }
      return std::nullopt;
    }
    if (start_i > 0)
    {
      return refs[start_i - 1];
    }
    return std::nullopt;
  };

  auto relative_intersection_id = [&](const auto& refs, KineticDelaunay::CrossingData::EdgeIntersectionRef target) -> size_t {
    size_t i = 0;
    for (auto ref : refs)
    {
      if (ref == target)
      {
        return i;
      }
      ++i;
    }
    return size_t(-1);
  };

  auto log_intersection = [&](const char* tag, KineticDelaunay::CrossingData::EdgeIntersectionRef ref) {
    const auto& d_refs = crossing_data.delaunay_edge_intersections[ref->delaunay_edge_id];
    const auto& v_refs = crossing_data.voronoi_edge_intersections[ref->voronoi_edge_id];
    const size_t rel_d = relative_intersection_id(d_refs, ref);
    const size_t rel_v = relative_intersection_id(v_refs, ref);
    KINDS_DEBUG("Radius trace " << tag << " face=" << affected_face_id << " t=" << t << " de=" << ref->delaunay_edge_id
                                << " ve=" << ref->voronoi_edge_id << " rel_i=(" << rel_d << "," << rel_v << ")"
                                << " de_param=" << ref->delaunay_edge_param);
  };

  // One representative intersection per touched cell.
  std::unordered_map<size_t, KineticDelaunay::CrossingData::EdgeIntersectionRef> representative_by_cell;
  for (auto ref : seed_intersections)
  {
    const auto ep = edge_endpoints(ref->voronoi_edge_id);
    if (ep[0] >= 0)
    {
      representative_by_cell.emplace(static_cast<size_t>(ep[0]), ref);
    }
    if (ep[1] >= 0)
    {
      representative_by_cell.emplace(static_cast<size_t>(ep[1]), ref);
    }
  }

  std::unordered_set<size_t> encountered_voronoi_edges_all;
  std::vector<SegmentBuilder::MeshingData> radius_operated_finished_strips;
  std::vector<SegmentBuilder::MeshingData> radius_operated_started_strips;

  for (const auto& entry : representative_by_cell)
  {
    const size_t cell_id = entry.first;
    auto seed_ref = entry.second;
    std::unordered_set<size_t> encountered_voronoi_edges;
    encountered_voronoi_edges.insert(seed_ref->voronoi_edge_id);
    std::vector<RadiusTraceVertex> polygon;
    bool success = false;
    const size_t max_steps = 128;
    std::string terminal_fail_reason;

    std::string fail_reason;
    polygon.clear();
    KINDS_DEBUG("Starting new trace for cell " << cell_id << ":");
    polygon.push_back(make_intersection_trace_vertex(seed_ref));
    {
      const auto& p = polygon.back().position;
      const auto& d_refs = crossing_data.delaunay_edge_intersections[seed_ref->delaunay_edge_id];
      const auto& v_refs = crossing_data.voronoi_edge_intersections[seed_ref->voronoi_edge_id];
      const size_t rel_d = relative_intersection_id(d_refs, seed_ref);
      const size_t rel_v = relative_intersection_id(v_refs, seed_ref);
      KINDS_DEBUG("Added point: seed-intersection x=" << p.x << " y=" << p.y << " t=" << p.z << " de="
                                                      << seed_ref->delaunay_edge_id << " ve=" << seed_ref->voronoi_edge_id
                                                      << " rel_i=(" << rel_d << "," << rel_v << ")");
    }
    KINDS_DEBUG("Radius trace start face=" << affected_face_id << " cell=" << cell_id << " t=" << t);
    log_intersection("seed", seed_ref);

    auto current_ref = seed_ref;
    bool closed_to_seed = false;
    bool wrong_direction = false;
    size_t step_count = 0;

    while (!closed_to_seed && !wrong_direction && step_count < max_steps)
    {
      ++step_count;

      // Phase 1: try both directions independently.
      std::optional<KineticDelaunay::CrossingData::EdgeIntersectionRef> end_ref = std::nullopt;
      std::string phase1_wrong_direction_reason;
      bool phase1_non_direction_failure = false;
      size_t current_d_edge = current_ref->voronoi_edge_id;
      encountered_voronoi_edges.insert(current_d_edge);
      KINDS_DEBUG("Radius trace phase1 start face=" << affected_face_id << " cell=" << cell_id << " step=" << step_count
                                                    << " current_dual_edge=" << current_d_edge);
      log_intersection("phase1-current-intersection", current_ref);

      for (size_t phase1_attempt = 0; phase1_attempt < 2 && !end_ref.has_value() && !phase1_non_direction_failure; ++phase1_attempt)
      {
        const bool phase1_reverse = (phase1_attempt == 1);
        KINDS_DEBUG("Radius trace phase1 attempt face=" << affected_face_id << " cell=" << cell_id
                                                        << " attempt=" << (phase1_attempt + 1) << "/2 reverse_orientation="
                                                        << phase1_reverse);
        std::optional<size_t> prev_vv = std::nullopt;
        if (phase1_reverse)
        {
          prev_vv = voronoi_vertex_on_edge(current_d_edge, true);
        }
        bool phase1_wrong_direction_local = false;
        std::string phase1_fail_local;

        for (size_t hop = 0; hop < max_steps; ++hop)
        {
          auto vv_opt = choose_voronoi_vertex_towards_inside(current_d_edge, prev_vv);
          if (!vv_opt.has_value())
          {
            KINDS_DEBUG("Radius trace phase1 no inward Voronoi vertex face=" << affected_face_id << " cell=" << cell_id
                                                                              << " ve=" << current_d_edge);
            auto paired_ref_opt = other_triangle_intersection_on_same_voronoi_edge(current_d_edge, current_ref);
            if (paired_ref_opt.has_value())
            {
              end_ref = paired_ref_opt.value();
              encountered_voronoi_edges.insert(end_ref.value()->voronoi_edge_id);
              log_intersection("phase1-paired-end", end_ref.value());
            }
            else
            {
              phase1_fail_local = "phase1: no inward Voronoi vertex and no paired triangle intersection on dual edge "
                + std::to_string(current_d_edge);
            }
            break;
          }

          const size_t vv = vv_opt.value();
          KINDS_DEBUG("Radius trace phase1 visit voronoi-vertex face=" << affected_face_id << " cell=" << cell_id << " vv=" << vv
                                                                        << " via_ve=" << current_d_edge
                                                                        << " containing_tri=" << crossing_data.getContainingTriId(vv));
          if (crossing_data.getContainingTriId(vv) != affected_face_id)
          {
            phase1_wrong_direction_local = true;
            phase1_fail_local = "phase1: Voronoi vertex " + std::to_string(vv) + " not in affected triangle "
              + std::to_string(affected_face_id) + " (wrong direction)";
            break;
          }

          const glm::dvec3 vv_pos = segment_builder_.computeVoronoiVertex(graph.face(vv).half_edges[0], t);
          if (std::isfinite(vv_pos.x) && std::isfinite(vv_pos.y))
          {
            polygon.push_back(make_voronoi_trace_vertex(vv, glm::dvec3 { vv_pos.x, vv_pos.y, t }));
            const auto& p = polygon.back().position;
            KINDS_DEBUG("Added point: phase1-voronoi-vertex dual_vv=" << vv << " via_ve=" << current_d_edge << " x=" << p.x
                                                                       << " y=" << p.y << " t=" << p.z);
          }

          const auto tri_edges = dual_triangle_edges(vv);
          std::vector<size_t> candidates;
          for (size_t candidate_d : tri_edges)
          {
            if (candidate_d == current_d_edge)
            {
              continue;
            }
            const auto ep = edge_endpoints(candidate_d);
            if ((ep[0] >= 0 && static_cast<size_t>(ep[0]) == cell_id) || (ep[1] >= 0 && static_cast<size_t>(ep[1]) == cell_id))
            {
              candidates.push_back(candidate_d);
            }
          }
          if (candidates.empty())
          {
            phase1_wrong_direction_local = true;
            phase1_fail_local = "phase1: no incident dual edge from Voronoi vertex " + std::to_string(vv) + " for cell "
              + std::to_string(cell_id) + " (wrong direction)";
            break;
          }
          const size_t next_d = phase1_reverse ? candidates.back() : candidates.front();
          encountered_voronoi_edges.insert(next_d);
          KINDS_DEBUG("Radius trace phase1 select-next-ve face=" << affected_face_id << " cell=" << cell_id << " vv=" << vv
                                                                 << " current_ve=" << current_d_edge << " next_ve=" << next_d
                                                                 << " candidate_count=" << candidates.size());

          auto next_ref_opt = first_triangle_intersection_on_voronoi_edge(next_d, vv, current_ref);
          if (!next_ref_opt.has_value())
          {
            // No intersection on this Voronoi edge: continue VV→VV if the other endpoint is still
            // inside the affected triangle.
            const size_t vv_even = voronoi_vertex_on_edge(next_d, true);
            const size_t vv_odd = voronoi_vertex_on_edge(next_d, false);
            const size_t other_vv = (vv_even == vv) ? vv_odd : vv_even;
            if (crossing_data.isVoronoiVertexRegistered(other_vv)
              && crossing_data.getContainingTriId(other_vv) == affected_face_id)
            {
              KINDS_DEBUG("Radius trace phase1 VV→VV continue face=" << affected_face_id << " cell=" << cell_id
                                                                      << " from_vv=" << vv << " to_vv=" << other_vv
                                                                      << " via_ve=" << next_d);
              prev_vv = vv;
              current_d_edge = next_d;
              continue;
            }

            phase1_fail_local = "phase1: no triangle-edge intersection on dual edge " + std::to_string(next_d)
              + " from Voronoi vertex " + std::to_string(vv) + " and other VV " + std::to_string(other_vv)
              + " not in affected triangle";
            break;
          }

          end_ref = next_ref_opt.value();
          encountered_voronoi_edges.insert(end_ref.value()->voronoi_edge_id);
          log_intersection("phase1-end-intersection", end_ref.value());
          break;
        }

        if (!end_ref.has_value())
        {
          if (phase1_wrong_direction_local)
          {
            phase1_wrong_direction_reason = phase1_fail_local;
          }
          else
          {
            phase1_non_direction_failure = true;
            fail_reason = phase1_fail_local.empty() ? "phase1: Voronoi walk ended without a next triangle-edge intersection"
                                                    : phase1_fail_local;
          }
        }
      }

      if (!end_ref.has_value())
      {
        if (!phase1_non_direction_failure)
        {
          wrong_direction = true;
          fail_reason = phase1_wrong_direction_reason.empty() ? "phase1: wrong direction in both attempts"
                                                              : phase1_wrong_direction_reason;
        }
        break;
      }

      if (end_ref.value() == seed_ref)
      {
        KINDS_DEBUG("Radius trace close-to-seed in phase1 face=" << affected_face_id << " cell=" << cell_id);
        closed_to_seed = true;
        break;
      }
      polygon.push_back(make_intersection_trace_vertex(end_ref.value()));
      {
        const auto added_ref = end_ref.value();
        const auto& p = polygon.back().position;
        const auto& d_refs = crossing_data.delaunay_edge_intersections[added_ref->delaunay_edge_id];
        const auto& v_refs = crossing_data.voronoi_edge_intersections[added_ref->voronoi_edge_id];
        const size_t rel_d = relative_intersection_id(d_refs, added_ref);
        const size_t rel_v = relative_intersection_id(v_refs, added_ref);
        KINDS_DEBUG("Added point: phase1-end-intersection x=" << p.x << " y=" << p.y << " t=" << p.z << " de="
                                                               << added_ref->delaunay_edge_id << " ve="
                                                               << added_ref->voronoi_edge_id << " rel_i=(" << rel_d << ","
                                                               << rel_v << ")");
      }

      // Phase 2: try both directions independently.
      KINDS_DEBUG("Radius trace phase2 start face=" << affected_face_id << " cell=" << cell_id << " end_de="
                                                    << end_ref.value()->delaunay_edge_id << " seed_de="
                                                    << seed_ref->delaunay_edge_id);
      log_intersection("phase2-end-intersection", end_ref.value());
      const auto end_idx_opt = edge_index_in_affected_triangle(end_ref.value()->delaunay_edge_id);
      const auto seed_idx_opt = edge_index_in_affected_triangle(seed_ref->delaunay_edge_id);
      if (!end_idx_opt.has_value() || !seed_idx_opt.has_value())
      {
        fail_reason = "phase2: delaunay edge not on affected triangle (end_edge="
          + std::to_string(end_ref.value()->delaunay_edge_id) + " seed_edge=" + std::to_string(seed_ref->delaunay_edge_id) + ")";
        break;
      }

      bool phase2_advanced = false;
      bool phase2_closed = false;
      std::string phase2_wrong_direction_reason;
      std::string phase2_non_direction_reason;

      for (size_t phase2_attempt = 0; phase2_attempt < 2 && !phase2_advanced && !phase2_closed; ++phase2_attempt)
      {
        const bool phase2_reverse = (phase2_attempt == 1);
        KINDS_DEBUG("Radius trace phase2 attempt face=" << affected_face_id << " cell=" << cell_id
                                                        << " attempt=" << (phase2_attempt + 1) << "/2 reverse_orientation="
                                                        << phase2_reverse);

        size_t idx = end_idx_opt.value();
        bool local_wrong_direction = false;
        std::string local_fail_reason;
        std::vector<RadiusTraceVertex> local_polygon = polygon;
        auto phase2_current_ref = end_ref.value();

        while (true)
        {
          const size_t he_id = affected_face_he[idx];
          const size_t d_edge_id = he_id / 2;
          const bool traverse_he_forward = !phase2_reverse;
          KINDS_DEBUG("Radius trace phase2 inspect triangle-edge face=" << affected_face_id << " cell=" << cell_id
                                                                         << " idx=" << idx << " he=" << he_id << " de=" << d_edge_id
                                                                         << " traverse_he_forward=" << traverse_he_forward);
          const auto start_ref = (phase2_current_ref->delaunay_edge_id == d_edge_id)
            ? std::optional<KineticDelaunay::CrossingData::EdgeIntersectionRef>(phase2_current_ref)
            : std::nullopt;
          if (start_ref.has_value())
          {
            log_intersection("phase2-inspect-start-intersection", start_ref.value());
          }
          auto boundary_intersection_opt = phase2_intersection_on_triangle_edge(he_id, traverse_he_forward, start_ref);
          if (boundary_intersection_opt.has_value())
          {
            auto boundary_ref = boundary_intersection_opt.value();
            encountered_voronoi_edges.insert(boundary_ref->voronoi_edge_id);
            log_intersection("phase2-boundary-intersection", boundary_ref);
            if (!edge_touches_cell(boundary_ref->voronoi_edge_id, cell_id))
            {
              local_wrong_direction = true;
              local_fail_reason = "phase2: boundary intersection uses voronoi_edge " + std::to_string(boundary_ref->voronoi_edge_id)
                + " that does not touch cell " + std::to_string(cell_id) + " (wrong direction)";
              break;
            }

            if (boundary_ref == seed_ref)
            {
              phase2_closed = true;
              break;
            }
            if (boundary_ref == end_ref.value())
            {
              local_fail_reason = "phase2: returned to the same intersection it started from (de="
                + std::to_string(boundary_ref->delaunay_edge_id) + " ve=" + std::to_string(boundary_ref->voronoi_edge_id) + ")";
              KINDS_ERROR("Radius strand meshlet: " << local_fail_reason << " face=" << affected_face_id << " cell=" << cell_id
                                                    << " t=" << occurrence);
              break;
            }

            local_polygon.push_back(make_intersection_trace_vertex(boundary_ref));
            {
              const auto& p = local_polygon.back().position;
              const auto& d_refs = crossing_data.delaunay_edge_intersections[boundary_ref->delaunay_edge_id];
              const auto& v_refs = crossing_data.voronoi_edge_intersections[boundary_ref->voronoi_edge_id];
              const size_t rel_d = relative_intersection_id(d_refs, boundary_ref);
              const size_t rel_v = relative_intersection_id(v_refs, boundary_ref);
              KINDS_DEBUG("Added point: phase2-boundary-intersection x=" << p.x << " y=" << p.y << " t=" << p.z << " de="
                                                                          << boundary_ref->delaunay_edge_id << " ve="
                                                                          << boundary_ref->voronoi_edge_id << " rel_i=("
                                                                          << rel_d << "," << rel_v << ")");
            }
            phase2_current_ref = boundary_ref;
            current_ref = boundary_ref;
            log_intersection("phase2-next-current", current_ref);
            polygon = std::move(local_polygon);
            phase2_advanced = true;
            break;
          }

          const bool traverse_even_forward = (he_id % 2 == 0) ? traverse_he_forward : !traverse_he_forward;
          const size_t he_even = 2 * d_edge_id;
          const size_t tri_vertex_id = traverse_even_forward ? graph.destination(he_even) : graph.halfEdge(he_even).origin;
          KINDS_DEBUG("Radius trace phase2 visit site (triangle-vertex) face=" << affected_face_id << " cell=" << cell_id
                                                                                << " site=" << tri_vertex_id << " on_he=" << he_id
                                                                                << " on_de=" << d_edge_id
                                                                                << " traverse_even_forward=" << traverse_even_forward);
          if (tri_vertex_id != size_t(-1))
          {
            if (tri_vertex_id != cell_id)
            {
              local_wrong_direction = true;
              local_fail_reason = "phase2: triangle corner vertex " + std::to_string(tri_vertex_id) + " != cell "
                + std::to_string(cell_id) + " (wrong direction)";
              break;
            }
            const glm::dvec2 p2 = segment_builder_.kin_del.getPointInDelaunaySpace(tri_vertex_id, t);
            local_polygon.push_back(make_site_trace_vertex(tri_vertex_id, glm::dvec3 { p2.x, p2.y, t }));
            {
              const auto& p = local_polygon.back().position;
              KINDS_DEBUG("Added point: phase2-triangle-vertex site=" << tri_vertex_id << " on_he="
                                                                       << affected_face_he[idx] << " on_de=" << d_edge_id
                                                                       << " x=" << p.x << " y=" << p.y << " t=" << p.z);
            }
          }
          // We must inspect the seed edge too: even after reaching it, there can still be
          // a directed step along that edge before we hit the seed intersection.
          if (idx == seed_idx_opt.value())
          {
            break;
          }
          idx = phase2_reverse ? ((idx + 2) % 3) : ((idx + 1) % 3);
        }

        if (phase2_closed)
        {
          // Commit points accumulated in this successful phase2 attempt before closing.
          polygon = std::move(local_polygon);
          closed_to_seed = true;
          break;
        }
        if (phase2_advanced)
        {
          break;
        }

        if (local_wrong_direction)
        {
          phase2_wrong_direction_reason = local_fail_reason;
        }
        else
        {
          if (phase2_non_direction_reason.empty())
          {
            phase2_non_direction_reason
              = local_fail_reason.empty()
              ? "phase2: directed triangle-edge traversal (including seed edge) found no valid next boundary intersection (considering half-edge parity and forward/backward walk direction)"
              : local_fail_reason;
          }
        }
      }

      if (!closed_to_seed && !phase2_advanced)
      {
        if (!phase2_non_direction_reason.empty())
        {
          fail_reason = phase2_non_direction_reason;
        }
        else
        {
          wrong_direction = true;
          fail_reason = phase2_wrong_direction_reason.empty() ? "phase2: wrong direction in both attempts"
                                                              : phase2_wrong_direction_reason;
        }
        break;
      }
    }

    if (!closed_to_seed && !wrong_direction && step_count >= max_steps)
    {
      fail_reason = "traversal step limit reached (" + std::to_string(max_steps) + " steps)";
    }

    if (!wrong_direction && closed_to_seed && polygon.size() >= 3)
    {
      success = true;
    }
    else if (!wrong_direction)
    {
      if (fail_reason.empty())
      {
        if (closed_to_seed && polygon.size() < 3)
        {
          fail_reason = "polygon closed to seed but only " + std::to_string(polygon.size()) + " vertices (need >= 3)";
        }
        else
        {
          fail_reason = "traversal incomplete: not closed to seed without direction error";
        }
      }
      terminal_fail_reason = fail_reason;
      KINDS_ERROR("Radius strand meshlet: non-retryable failure face=" << affected_face_id << " cell=" << cell_id
        << " t=" << occurrence << " reason=" << terminal_fail_reason);
    }
    else
    {
      terminal_fail_reason = fail_reason;
      KINDS_ERROR("Radius strand meshlet: wrong direction after phase retries face=" << affected_face_id << " cell=" << cell_id
        << " t=" << occurrence << " reason=" << terminal_fail_reason);
    }

    auto emit_radius_cell_mesh = [&](const std::vector<RadiusTraceVertex>& poly)
    {
      if (poly.empty())
      {
        return;
      }
      std::string radius_triangulation_meta;
      if (segment_builder_.store_mesh_metadata)
      {
        radius_triangulation_meta = SegmentBuilder::MetadataBuilder()
                                      .addString("event_type", "radius_event")
                                      .addString("mesh_type", "bark")
                                      .addString("segment_action", "new_segment")
                                      .addDouble("time", t)
                                      .addDouble("t", t)
                                      .addSize("delaunay_face_id", affected_face_id)
                                      .addSize("strand_cell_id", cell_id)
                                      .addString("op", "radius_strand_cell_triangulation")
                                      .build();
      }
      VoronoiMesh mesh;
      segment_builder_.configureMeshletStorage(mesh);
      std::vector<size_t> ids;
      ids.reserve(poly.size());
      const auto& boundary_polygon = segment_builder_.kin_del.component_data.component_boundaries[component_id][0];
      const auto& centroid = segment_builder_.kin_del.component_data.component_centroids[component_id];
      for (const auto& vert : poly)
      {
        const char* source = vert.voronoi_vertex_id.has_value()
          ? "Voronoi vertex"
          : (vert.intersection.has_value() ? "intersection" : "site");
        std::string radius_vertex_meta;
        if (segment_builder_.store_mesh_metadata)
        {
          SegmentBuilder::MetadataBuilder builder;
          builder.addString("event_type", "radius_event")
            .addString("source", source)
            .addString("mesh_type", "bark")
            .addString("segment_action", "new_segment")
            .addDouble("time", t)
            .addDouble("t", t)
            .addSize("delaunay_face_id", affected_face_id)
            .addSize("strand_cell_id", cell_id);
          if (vert.voronoi_vertex_id.has_value())
          {
            builder.addSize("voronoi_vertex_id", vert.voronoi_vertex_id.value());
          }
          radius_vertex_meta = builder.build();
        }
        const size_t strand_for_vertex = vert.site_strand_id.value_or(cell_id);
        SegmentBuilder::MeshletVertexRuntimeInfo runtime_info {};
        // Exact Delaunay-plane XY used to construct the polygon ring (stored for profile/UV; fan triangulation
        // does not re-query geometry).
        runtime_info.triangulation_plane_xy = glm::dvec2(vert.position.x, vert.position.y);
        if (vert.intersection.has_value())
        {
          runtime_info.position_intersection = vert.intersection;
          runtime_info.conceptual_intersection = vert.intersection;
        }
        ids.push_back(segment_builder_.addMeshletVertex(mesh, boundary_polygon, centroid, vert.position, strand_for_vertex,
          t, false, vert.voronoi_vertex_id, radius_vertex_meta, std::nullopt, runtime_info));
      }
      double signed_area2 = 0.0;
      for (size_t i = 0; i < poly.size(); ++i)
      {
        const glm::dvec2 p0(poly[i].position.x, poly[i].position.y);
        const glm::dvec2 p1(poly[(i + 1) % poly.size()].position.x, poly[(i + 1) % poly.size()].position.y);
        signed_area2 += p0.x * p1.y - p1.x * p0.y;
      }
      constexpr double area_eps = 1e-12;
      const bool polygon_ccw = signed_area2 > area_eps;
      const bool orient_upwards = (std::abs(signed_area2) <= area_eps)
        ? inside_to_outside
        : (inside_to_outside ? polygon_ccw : !polygon_ccw);
      // Radius traced cell rings are convex — fan triangulation only (no ear-clip / plane geometry).
      // Non-shift radius cell fans close the bark surface of the transitioning Delaunay triangle.
      segment_builder_.fanTriangulateConvexPolygon(mesh, ids, radius_triangulation_meta,
        SegmentBuilder::BoundaryIntervalMeshletMaterialId, orient_upwards);
      if (ids.size() < 3)
      {
        KINDS_WARNING("Radius: traced cell polygon has fewer than three vertices; no triangles emitted for cell "
          << cell_id << " face=" << affected_face_id << " t=" << occurrence << " runtime_branch=" << runtime_branch_id
          << "; " << segment_builder_.formatStrandBranchLogInfo(cell_id, t));
      }
      size_t owner_segment_id = static_cast<size_t>(-1);
      if (cell_id < segment_builder_.strand_to_segment_indices.size()
        && !segment_builder_.strand_to_segment_indices[cell_id].empty())
      {
        owner_segment_id = segment_builder_.strand_to_segment_indices[cell_id].back();
      }
      const size_t stored_segment_pair_index = segment_builder_.segment_mesh_pairs.size();
      segment_builder_.segment_mesh_pairs.push_back(
        MeshStructure::SegmentMeshPair { owner_segment_id, static_cast<size_t>(-1), 0, 0, 1 });
      std::string suffix = std::string("_delaunay") + std::to_string(affected_face_id) + "_strand"
        + std::to_string(cell_id);
      KINDS_DEBUG("Radius: stored extracted strand meshlet segment segment_mesh_pairs_index=" << stored_segment_pair_index
                                                                                              << " cell_id=" << cell_id
                                                                                              << " owner_segment_id=" << owner_segment_id
                                                                                              << " delaunay_face=" << affected_face_id
                                                                                              << " t=" << t << " polygon_vertices=" << poly.size()
                                                                                              << " meshlet_suffix=" << suffix);
      segment_builder_.registerMeshletWithSuffix(std::move(mesh), std::move(suffix), t);
      segment_builder_.segment_mesh_pair_last_left_and_right_vertex.emplace_back();
    };

    if (!success)
    {
      // Incomplete / wrong-direction rings must not reach ear-clip; dump ring-walk debug SVG/TXT instead.
      std::vector<RadiusRingWalkDebugVertex> debug_ring;
      debug_ring.reserve(polygon.size());
      std::vector<std::string> incorrect_vertices;
      std::vector<std::string> unmatched_vertices;
      for (const RadiusTraceVertex& vert : polygon)
      {
        RadiusRingWalkDebugVertex debug_vert;
        debug_vert.xy = glm::dvec2(vert.position.x, vert.position.y);
        if (vert.voronoi_vertex_id.has_value())
        {
          debug_vert.source = "Voronoi vertex";
          debug_vert.voronoi_vertex_id = vert.voronoi_vertex_id;
        }
        else if (vert.intersection.has_value())
        {
          debug_vert.source = "intersection";
          debug_vert.delaunay_edge_id = vert.intersection.value()->delaunay_edge_id;
          debug_vert.voronoi_edge_id = vert.intersection.value()->voronoi_edge_id;
        }
        else
        {
          debug_vert.source = "site";
          debug_vert.strand_id = vert.site_strand_id.value_or(cell_id);
        }
        debug_ring.push_back(std::move(debug_vert));
      }

      const std::string reason = terminal_fail_reason.empty() ? fail_reason : terminal_fail_reason;
      if (reason.find("no ") != std::string::npos || reason.find("found no") != std::string::npos
        || reason.find("incomplete") != std::string::npos || reason.find("step limit") != std::string::npos)
      {
        unmatched_vertices.push_back("missing_step: " + reason);
      }
      if (segment_builder_.shouldDumpErrorFiles())
      {
        writeRadiusRingWalkFailDebug(segment_builder_.kin_del, occurrence, runtime_branch_id, affected_face_id, cell_id, reason,
          std::move(debug_ring), std::move(incorrect_vertices), std::move(unmatched_vertices));
      }
      encountered_voronoi_edges_all.insert(encountered_voronoi_edges.begin(), encountered_voronoi_edges.end());
      continue;
    }

    if (!use_radius_boundary_shift)
    {
      emit_radius_cell_mesh(polygon);
    }
    encountered_voronoi_edges_all.insert(encountered_voronoi_edges.begin(), encountered_voronoi_edges.end());
  }


  if (segment_builder_.kin_del.computeBoundaryOnTheFly())
  {
    // Close/extend existing strip meshes at the current kinetic time, then rebuild crossing iterators for the whole
    // meshing state (full `computeEdgeIntersections` invalidates all `EdgeIntersectionRef`s), then seed fresh strips.
    for (size_t voronoi_edge_id : encountered_voronoi_edges_all)
    {
      const size_t he_even = 2 * voronoi_edge_id;
      if (he_even >= graph.halfEdgeSlotCount())
      {
        continue;
      }
      if (graph.isInfinite(he_even))
      {
        continue;
      }
      const auto& he = graph.halfEdge(he_even);
      const auto& twin_he = graph.halfEdge(he_even ^ 1);
      const size_t vertex = std::max(he.origin, twin_he.origin);
      const size_t component_id = segment_builder_.kin_del.component_data.component_map[vertex];
      std::vector<bool> he_visited(graph.halfEdgeSlotCount(), false);
      segment_builder_.updateBoundary(t, he_visited, component_id);
      auto& boundary_polygon = segment_builder_.kin_del.component_data.component_boundaries[component_id][0];
      std::vector<SegmentBuilder::MeshingData> finished_strips
        = segment_builder_.finishMesh(he_even, t, boundary_polygon, SegmentBuilder::BoundaryEventType::Radius,
            SegmentBuilder::BoundarySegmentAction::SegmentCompleted, radius_boundary_shift_arg,
            radius_only_adjacent_strips);
      radius_operated_finished_strips.insert(radius_operated_finished_strips.end(), finished_strips.begin(), finished_strips.end());
    }
  }

  auto for_each_triangle_edge_crossing_with_boundary_change
    = [&](bool became_boundary, const auto& fn)
  {
    for (size_t i = 0; i < 3; ++i)
    {
      const size_t he = radius_pre_face_he_[i];
      const size_t d_edge = he / 2;
      const size_t neighbor_face = graph.halfEdge(he ^ 1).face;
      // IO↔II only: the neighboring triangle stays inside, so the intersection moves between
      // the alpha boundary and the interior (not OO↔OI).
      if (!segment_builder_.kin_del.getFaceInside(neighbor_face))
      {
        continue;
      }
      const bool pre_boundary = radius_pre_is_boundary_edge_[i];
      const bool post_boundary = segment_builder_.kin_del.isOnComponentBoundary(he);
      if (became_boundary)
      {
        if (pre_boundary || !post_boundary)
        {
          continue;
        }
      }
      else if (!pre_boundary || post_boundary)
      {
        continue;
      }
      if (d_edge >= crossing_data.delaunay_edge_intersections.size())
      {
        continue;
      }
      for (KineticDelaunay::CrossingData::EdgeIntersectionRef ref : crossing_data.delaunay_edge_intersections[d_edge])
      {
        fn(ref);
      }
    }
  };

  if (!use_radius_boundary_shift && segment_builder_.kin_del.computeBoundaryOnTheFly())
  {
    // Interior → boundary: finishMesh just created the top triangle; split it now.
    // Neighboring crossings on the same Voronoi edge (triangle between them flipped) are both inserted,
    // in order along the top edge so the second split uses the remaining sub-triangle.
    std::vector<KineticDelaunay::CrossingData::EdgeIntersectionRef> became_boundary_crossings;
    for_each_triangle_edge_crossing_with_boundary_change(true, [&](KineticDelaunay::CrossingData::EdgeIntersectionRef ref)
    {
      if (encountered_voronoi_edges_all.find(ref->voronoi_edge_id) != encountered_voronoi_edges_all.end())
      {
        became_boundary_crossings.push_back(ref);
      }
    });
    segment_builder_.applyOrQueueRadiusNonShiftInteriorCrossingSplits(became_boundary_crossings, t, true);
  }

  for (size_t voronoi_edge_id : encountered_voronoi_edges_all)
  {
    const size_t he_even = 2 * voronoi_edge_id;
    if (he_even >= segment_builder_.kin_del.getGraph().halfEdgeSlotCount())
    {
      continue;
    }
    std::vector<SegmentBuilder::MeshingData> started_strips = segment_builder_.startNewMesh(he_even, t, true,
      SegmentBuilder::BoundaryEventType::Radius, SegmentBuilder::BoundarySegmentAction::NewSegment, radius_boundary_shift_arg,
      radius_only_adjacent_strips);
    radius_operated_started_strips.insert(radius_operated_started_strips.end(), started_strips.begin(), started_strips.end());
  }

  if (!use_radius_boundary_shift && segment_builder_.kin_del.computeBoundaryOnTheFly())
  {
    // Boundary → interior: reuse startNewMesh seeded two vertices and no triangle yet; split in finalize.
    // Symmetric neighboring pair: both crossings lie on the upcoming seed edge and are queued in order.
    std::vector<KineticDelaunay::CrossingData::EdgeIntersectionRef> left_boundary_crossings;
    for_each_triangle_edge_crossing_with_boundary_change(false, [&](KineticDelaunay::CrossingData::EdgeIntersectionRef ref)
    {
      if (encountered_voronoi_edges_all.find(ref->voronoi_edge_id) != encountered_voronoi_edges_all.end())
      {
        left_boundary_crossings.push_back(ref);
      }
    });
    segment_builder_.applyOrQueueRadiusNonShiftInteriorCrossingSplits(left_boundary_crossings, t, false);
  }

  std::unordered_set<size_t> triangle_voronoi_edge_ids;
  for (size_t i = 0; i < 3; ++i)
  {
    triangle_voronoi_edge_ids.insert(affected_face_he[i] / 2);
  }

  std::vector<KineticDelaunay::CrossingData::EdgeIntersectionRef> radius_strip_endpoint_intersections_off_triangle;
  std::unordered_set<IntersectionKey, IntersectionKeyHash> off_triangle_strip_endpoint_keys;
  auto collect_strip_endpoint_intersections_off_triangle = [&](const std::vector<SegmentBuilder::MeshingData>& strips)
  {
    for (const SegmentBuilder::MeshingData& strip : strips)
    {
      const std::optional<KineticDelaunay::CrossingData::EdgeIntersectionRef>* endpoints[]
        = { &strip.start_crossing, &strip.end_crossing };
      for (const auto* endpoint : endpoints)
      {
        if (!endpoint->has_value())
        {
          continue;
        }
        const KineticDelaunay::CrossingData::EdgeIntersectionRef ref = endpoint->value();
        if (triangle_voronoi_edge_ids.find(ref->voronoi_edge_id) != triangle_voronoi_edge_ids.end())
        {
          continue;
        }
        if (off_triangle_strip_endpoint_keys.insert(ref_key(ref)).second)
        {
          radius_strip_endpoint_intersections_off_triangle.push_back(ref);
        }
      }
    }
  };
  collect_strip_endpoint_intersections_off_triangle(radius_operated_finished_strips);
  collect_strip_endpoint_intersections_off_triangle(radius_operated_started_strips);

  KINDS_DEBUG("Radius: strip endpoint intersections off updated triangle voronoi edges count="
              << radius_strip_endpoint_intersections_off_triangle.size() << " face=" << affected_face_id << " t=" << t);

  for (const KineticDelaunay::CrossingData::EdgeIntersectionRef& endpoint_ref :
    radius_strip_endpoint_intersections_off_triangle)
  {
    if (endpoint_ref->prev_segment_mesh_pair_index != static_cast<size_t>(-1))
    {
      segment_builder_.extendIntersectionMeshAtSharedCrossing(endpoint_ref->prev_segment_mesh_pair_index, endpoint_ref, false,
        t, SegmentBuilder::BoundaryEventType::Radius, SegmentBuilder::BoundarySegmentAction::SegmentRemapped,
        radius_boundary_shift_arg, true);
    }
    if (endpoint_ref->next_segment_mesh_pair_index != static_cast<size_t>(-1))
    {
      segment_builder_.extendIntersectionMeshAtSharedCrossing(endpoint_ref->next_segment_mesh_pair_index, endpoint_ref, true, t,
        SegmentBuilder::BoundaryEventType::Radius, SegmentBuilder::BoundarySegmentAction::SegmentRemapped,
        radius_boundary_shift_arg, true);
    }
  }

  // After the radius topology update, start/reseed all boundary-interval meshes on
  // boundary Delaunay edges of the affected (now updated) triangle.

  std::optional<size_t> radius_shift_shared_site;
  const RadiusBoundaryTransitionShiftContext* roles_for_shared_site
    = radius_boundary_shift_arg != nullptr
    ? radius_boundary_shift_arg
    : (radius_boundary_shift_ctx.roles_valid ? &radius_boundary_shift_ctx : nullptr);
  if (roles_for_shared_site != nullptr)
  {
    const size_t s0 = roles_for_shared_site->source_delaunay_edges[0];
    const size_t s1 = roles_for_shared_site->source_delaunay_edges[1];
    const int a0 = graph.halfEdge(2 * s0).origin;
    const int a1 = graph.halfEdge(2 * s0 + 1).origin;
    const int b0 = graph.halfEdge(2 * s1).origin;
    const int b1 = graph.halfEdge(2 * s1 + 1).origin;
    for (int a : { a0, a1 })
    {
      if (a < 0)
      {
        continue;
      }
      if (a == b0 || a == b1)
      {
        radius_shift_shared_site = static_cast<size_t>(a);
        break;
      }
    }
  }

  const bool radius_shift_two_to_one
    = radius_pre_boundary_edge_count_ == 2 && post_boundary_edge_count == 1;
  const bool radius_shift_one_to_two
    = radius_pre_boundary_edge_count_ == 1 && post_boundary_edge_count == 2;
  if (radius_boundary_shift_arg != nullptr && radius_shift_shared_site.has_value() && radius_shift_one_to_two)
  {
    // 1->2: complementary mid was finished in beforeEvent on the old boundary (target_d). Split its top
    // triangle now using the shifted corner site. Pair ids were buffered then — afterEvent clears
    // CrossingData links on the edge that left the alpha boundary, so do not re-resolve via prev/next.
    segment_builder_.maybeQueueRadiusComplementarySplitForExistingMid(
      t, radius_shift_shared_site.value(), radius_boundary_shift_arg, radius_pre_finished_one_edge_meshlet_ids_);
  }
  radius_pre_finished_one_edge_meshlet_ids_.clear();

  std::unordered_set<size_t> started_boundary_he_even;
  for (size_t he_id : affected_face_he)
  {
    const size_t he_even = he_id & ~1;
    if (!started_boundary_he_even.insert(he_even).second)
    {
      continue;
    }
    if (!segment_builder_.kin_del.isOnComponentBoundary(he_even))
    {
      continue;
    }

    const size_t d_edge_id = he_even / 2;
    const auto& d_intersections = segment_builder_.kin_del.getCrossingData().delaunay_edge_intersections[d_edge_id];
    if (d_intersections.empty())
    {
      continue;
    }

    const std::vector<KineticDelaunay::CrossingData::EdgeIntersectionRef> refs
      = segment_builder_.getBoundaryIntersectionsInBoundaryOrder(d_edge_id);

    KINDS_DEBUG("Radius: extracted boundary intersection segment for reseed d_edge_id=" << d_edge_id << " he_even=" << he_even
                                                                                        << " ordered_crossing_count=" << refs.size()
                                                                                        << " t=" << t);

    auto format_crossing = [](KineticDelaunay::CrossingData::EdgeIntersectionRef r) {
      std::ostringstream o;
      o << "de=" << r->delaunay_edge_id << " ve=" << r->voronoi_edge_id << " param=" << r->delaunay_edge_param;
      return o.str();
    };

    {
      const size_t first_cell
        = segment_builder_.determineVoronoiCellForBoundaryIntersectionInterval(d_edge_id, std::nullopt, refs.front());
      KINDS_DEBUG("Radius: passing extracted segment interval [null,first_crossing] to startNewMeshFromIntersections voronoi_cell="
                   << first_cell << " t=" << t << " first_crossing={" << format_crossing(refs.front()) << "}");
      segment_builder_.startNewMeshFromIntersections(
        first_cell, t, std::nullopt, refs.front(), false, SegmentBuilder::BoundaryEventType::Radius,
        SegmentBuilder::BoundarySegmentAction::NewSegment, false, radius_boundary_shift_arg);
    }
    for (size_t k = 0; k + 1 < refs.size(); ++k)
    {
      const size_t mid_cell
        = segment_builder_.determineVoronoiCellForBoundaryIntersectionInterval(d_edge_id, refs[k], refs[k + 1]);
      if (mid_cell == static_cast<size_t>(-1))
      {
        continue;
      }
      KINDS_DEBUG("Radius: passing extracted segment interval [crossing,crossing] to startNewMeshFromIntersections k=" << k
                   << " voronoi_cell=" << mid_cell << " t=" << t << " start={" << format_crossing(refs[k]) << "} end={"
                   << format_crossing(refs[k + 1]) << "}");
      const size_t mid_pair = segment_builder_.startNewMeshFromIntersections(
        mid_cell, t, refs[k], refs[k + 1], false, SegmentBuilder::BoundaryEventType::Radius,
        SegmentBuilder::BoundarySegmentAction::NewSegment, false, radius_boundary_shift_arg);
      if (radius_boundary_shift_arg != nullptr && radius_shift_shared_site.has_value() && radius_shift_two_to_one
        && d_edge_id == radius_boundary_shift_arg->target_delaunay_edge)
      {
        segment_builder_.maybeQueueRadiusComplementarySplitForStartedMid(
          mid_pair, refs[k], refs[k + 1], t, radius_shift_shared_site.value(), radius_boundary_shift_arg);
      }
    }
    {
      const size_t last_cell
        = segment_builder_.determineVoronoiCellForBoundaryIntersectionInterval(d_edge_id, refs.back(), std::nullopt);
      KINDS_DEBUG("Radius: passing extracted segment interval [last_crossing,null] to startNewMeshFromIntersections voronoi_cell="
                   << last_cell << " t=" << t << " last_crossing={" << format_crossing(refs.back()) << "}");
      segment_builder_.startNewMeshFromIntersections(
        last_cell, t, refs.back(), std::nullopt, false, SegmentBuilder::BoundaryEventType::Radius,
        SegmentBuilder::BoundarySegmentAction::NewSegment, false, radius_boundary_shift_arg);
    }
  }
}
} // namespace kinDS


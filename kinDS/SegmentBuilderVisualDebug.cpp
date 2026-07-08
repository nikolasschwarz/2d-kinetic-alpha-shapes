#include "SegmentBuilderVisualDebug.hpp"

#include "HalfEdgeDelaunayGraphToSVG.hpp"
#include "KineticDelaunay.hpp"
#include "KineticDelaunayCrossingEvent.hpp"

#include <cmath>
#include <filesystem>
#include <limits>
#include <unordered_map>
#include <unordered_set>

namespace kinDS
{
namespace
{
size_t inputBranchSectionIndex(double t, size_t tree_height)
{
  const size_t lower_index = static_cast<size_t>(std::floor(t));
  const double frac = t - static_cast<double>(lower_index);
  size_t section_index = lower_index;
  if (frac > std::numeric_limits<double>::epsilon())
  {
    section_index = lower_index + 1;
  }
  if (tree_height > 0 && section_index >= tree_height)
  {
    section_index = tree_height - 1;
  }
  return section_index;
}

std::unordered_map<size_t, size_t> buildSiteInputBranchLabels(
  KineticDelaunay& kin_del, const HalfEdgeDelaunayGraph& graph, double occurrence_time)
{
  std::unordered_map<size_t, size_t> labels;
  const size_t vertex_count = graph.getVertexCount();
  const size_t section_index = inputBranchSectionIndex(occurrence_time, kin_del.getStrandTree().getHeight());

  for (size_t strand_id = 0; strand_id < vertex_count; ++strand_id)
  {
    if (kin_del.isDummyBoundary(strand_id))
    {
      continue;
    }

    labels[strand_id] = kin_del.getStrandTree().getBranchIndex(strand_id, section_index);
  }

  return labels;
}

void noteRuntimeBranchForStrand(
  std::unordered_set<size_t>& runtime_branches, const KineticDelaunay& kin_del, size_t strand_id)
{
  if (kin_del.isDummyBoundary(strand_id))
  {
    return;
  }

  const auto& runtime_branch_map = kin_del.getRuntimeBranchMap();
  if (strand_id >= runtime_branch_map.size())
  {
    return;
  }

  runtime_branches.insert(kin_del.unsplitRuntimeBranchId(runtime_branch_map[strand_id]));
}

std::optional<size_t> uniqueBranchOrNull(const std::unordered_set<size_t>& branches)
{
  if (branches.size() == 1)
  {
    return *branches.begin();
  }
  return std::nullopt;
}

std::optional<size_t> inferEventRuntimeBranchFromHighlight(
  const HalfEdgeDelaunayGraph& graph, const KineticDelaunay& kin_del, const VisualDebugHighlight& highlight)
{
  std::unordered_set<size_t> runtime_branches;

  for (size_t strand_id : highlight.delaunay_vertices)
  {
    noteRuntimeBranchForStrand(runtime_branches, kin_del, strand_id);
  }

  for (size_t he_id : highlight.directed_half_edges)
  {
    const int origin = graph.halfEdge(he_id).origin;
    if (origin >= 0)
    {
      noteRuntimeBranchForStrand(runtime_branches, kin_del, static_cast<size_t>(origin));
    }
    const int dest = graph.destination(he_id);
    if (dest >= 0)
    {
      noteRuntimeBranchForStrand(runtime_branches, kin_del, static_cast<size_t>(dest));
    }
  }

  for (size_t voronoi_vertex_id : highlight.voronoi_vertices)
  {
    if (voronoi_vertex_id >= graph.faceSlotCount() || !graph.isLiveFace(voronoi_vertex_id))
    {
      continue;
    }
    const auto verts = graph.getTriangleVertexIndices(voronoi_vertex_id);
    for (int v : verts)
    {
      if (v >= 0)
      {
        noteRuntimeBranchForStrand(runtime_branches, kin_del, static_cast<size_t>(v));
      }
    }
  }

  for (size_t voronoi_edge_id : highlight.voronoi_edges)
  {
    const size_t he_even = 2 * voronoi_edge_id;
    if (he_even + 1 >= graph.halfEdgeSlotCount())
    {
      continue;
    }
    const int origin_even = graph.halfEdge(he_even).origin;
    const int origin_odd = graph.halfEdge(he_even + 1).origin;
    if (origin_even >= 0)
    {
      noteRuntimeBranchForStrand(runtime_branches, kin_del, static_cast<size_t>(origin_even));
    }
    if (origin_odd >= 0)
    {
      noteRuntimeBranchForStrand(runtime_branches, kin_del, static_cast<size_t>(origin_odd));
    }
  }

  return uniqueBranchOrNull(runtime_branches);
}

/// Chronologically sortable phase token: "before" gets a '!' prefix so it sorts ahead of "after" (which is otherwise
/// alphabetically first) within the same event.
std::string chronologicalPhaseToken(const char* phase)
{
  return std::string(phase) == "before" ? "!before" : phase;
}

std::string visualDebugSvgRelativePath(double occurrence_time, const char* phase, const std::string& event_descriptor,
  std::optional<size_t> runtime_branch_id, const std::optional<std::filesystem::path>& output_root)
{
  const std::string basename = "t" + std::to_string(occurrence_time) + "_segmentbuilder_"
    + chronologicalPhaseToken(phase) + "_" + event_descriptor + ".svg";
  const std::string relative = runtime_branch_id.has_value()
    ? ("branch" + std::to_string(runtime_branch_id.value()) + "/" + basename)
    : ("branch0/" + basename);
  if (!output_root.has_value())
  {
    return relative;
  }
  return (*output_root / relative).generic_string();
}

std::unordered_set<size_t> collectLiveStrandIds(const HalfEdgeDelaunayGraph& graph)
{
  std::unordered_set<size_t> strand_ids;
  auto add_vertex = [&](int vertex)
  {
    if (vertex >= 0)
    {
      strand_ids.insert(static_cast<size_t>(vertex));
    }
  };

  for (size_t he_id : graph.liveDelaunayEdges())
  {
    add_vertex(graph.halfEdge(he_id).origin);
    add_vertex(graph.destination(he_id));
  }

  return strand_ids;
}

std::unordered_set<size_t> collectActiveRuntimeBranches(
  const KineticDelaunay& kin_del, const std::unordered_set<size_t>& live_strand_ids)
{
  std::unordered_set<size_t> runtime_branches;
  const auto& runtime_branch_map = kin_del.getRuntimeBranchMap();
  for (size_t strand_id : live_strand_ids)
  {
    if (strand_id < runtime_branch_map.size())
    {
      runtime_branches.insert(kin_del.unsplitRuntimeBranchId(runtime_branch_map[strand_id]));
    }
  }
  return runtime_branches;
}

std::unordered_set<size_t> collectStrandIdsForRuntimeBranch(
  const KineticDelaunay& kin_del, const std::unordered_set<size_t>& live_strand_ids, size_t runtime_branch_id)
{
  // @p runtime_branch_id is an unsplit branch id: gather its own strands plus those of any pending split-off children,
  // keeping only the live ones.
  std::unordered_set<size_t> strand_ids;
  for (size_t strand_id : kin_del.collectUnsplitRuntimeBranchStrands(runtime_branch_id))
  {
    if (live_strand_ids.count(strand_id) != 0)
    {
      strand_ids.insert(strand_id);
    }
  }
  return strand_ids;
}

std::pair<std::vector<glm::dvec2>, std::unordered_set<size_t>> buildVisualDebugStrandPositions(
  KineticDelaunay& kin_del, const HalfEdgeDelaunayGraph& graph, double occurrence_time,
  const std::unordered_set<size_t>& strand_ids)
{
  std::vector<glm::dvec2> points(graph.getVertexCount(), glm::dvec2(0.0));
  std::unordered_set<size_t> positioned_strands;
  positioned_strands.reserve(strand_ids.size());

  for (size_t strand_id : strand_ids)
  {
    if (kin_del.isDummyBoundary(strand_id))
    {
      continue;
    }
    points[strand_id] = kin_del.getPointAt(strand_id, occurrence_time);
    positioned_strands.insert(strand_id);
  }

  return { std::move(points), std::move(positioned_strands) };
}

std::unordered_map<size_t, glm::dvec3> buildSiteWorldPositions(
  KineticDelaunay& kin_del, double occurrence_time, const std::unordered_set<size_t>& strand_ids)
{
  std::unordered_map<size_t, glm::dvec3> world_positions;
  world_positions.reserve(strand_ids.size());

  for (size_t strand_id : strand_ids)
  {
    if (kin_del.isDummyBoundary(strand_id))
    {
      continue;
    }

    const glm::dvec2 profile_pos = kin_del.getPointAt(strand_id, occurrence_time);
    const size_t representative_strand_id = kin_del.representativeStrandIdForRuntimeBranch(strand_id);
    glm::dvec3 profile_pos_3d(profile_pos, occurrence_time);
    world_positions[strand_id]
      = kin_del.getStrandTree().transformToObjectSpace(profile_pos_3d, representative_strand_id, occurrence_time);
  }

  return world_positions;
}

std::unordered_map<size_t, glm::dvec3> buildVoronoiVertexWorldPositions(
  KineticDelaunay& kin_del, const HalfEdgeDelaunayGraph& graph, double occurrence_time)
{
  std::unordered_map<size_t, glm::dvec3> world_positions;

  for (size_t voronoi_vertex_id : graph.liveFaces())
  {
    const auto& face = graph.face(voronoi_vertex_id);
    size_t reference_strand_id = static_cast<size_t>(-1);
    for (size_t he_id : face.half_edges)
    {
      const int origin = graph.halfEdge(he_id).origin;
      if (origin >= 0 && !kin_del.isDummyBoundary(static_cast<size_t>(origin)))
      {
        reference_strand_id = static_cast<size_t>(origin);
        break;
      }
    }
    if (reference_strand_id == static_cast<size_t>(-1))
    {
      continue;
    }

    const glm::dvec3 profile_pos
      = kin_del.computeVoronoiVertexClampedInfinity(face.half_edges[0], occurrence_time);
    const size_t representative_strand_id = kin_del.representativeStrandIdForRuntimeBranch(reference_strand_id);
    glm::dvec3 world_pos_input = profile_pos;
    world_positions[voronoi_vertex_id]
      = kin_del.getStrandTree().transformToObjectSpace(world_pos_input, representative_strand_id, occurrence_time);
  }

  return world_positions;
}

void writeVisualDebugSvgFile(const std::string& relative_path, const std::vector<glm::dvec2>& points,
  const HalfEdgeDelaunayGraph& graph, KineticDelaunay& kin_del,
  const std::vector<size_t>& containing_tri_ids,
  const std::vector<HalfEdgeDelaunayGraphToSVG::IntersectionDebugInfo>& intersection_debug_data,
  const VisualDebugHighlight& highlight, const std::unordered_set<size_t>* highlighted_strands,
  const std::unordered_set<size_t>& positioned_strands,
  const std::unordered_map<size_t, size_t>& site_input_branch_labels,
  const std::unordered_map<size_t, glm::dvec3>& site_world_positions,
  const std::unordered_map<size_t, glm::dvec3>& voronoi_vertex_world_positions,
  const std::vector<HalfEdgeDelaunayGraphToSVG::SeparationOffsetSegment>* separation_offset_segments = nullptr,
  const std::vector<std::vector<glm::dvec2>>* seam_outlines = nullptr)
{
  const std::filesystem::path filepath(relative_path);
  if (filepath.has_parent_path())
  {
    std::filesystem::create_directories(filepath.parent_path());
  }
  HalfEdgeDelaunayGraphToSVG::write(points, graph, relative_path, 0.1, &kin_del.getFacesInside(), true,
    &containing_tri_ids, &intersection_debug_data, &highlight, highlighted_strands, &site_input_branch_labels,
    &positioned_strands, &site_world_positions, &voronoi_vertex_world_positions, separation_offset_segments,
    seam_outlines);
}
} // namespace

void writeSegmentBuilderVisualDebugSvg(bool visual_debug, KineticDelaunay& kin_del, const HalfEdgeDelaunayGraph& graph,
  double occurrence_time, const char* phase, const std::string& event_descriptor,
  const VisualDebugHighlight& highlight, std::optional<size_t> event_runtime_branch_id,
  const std::vector<HalfEdgeDelaunayGraphToSVG::SeparationOffsetSegment>* separation_offset_segments,
  const std::vector<std::vector<glm::dvec2>>* seam_outlines)
{
  if (!visual_debug)
  {
    return;
  }

  // Debug SVGs are always grouped by unsplit branch, so map any pending split-off child id back to its parent.
  if (event_runtime_branch_id.has_value())
  {
    event_runtime_branch_id = kin_del.unsplitRuntimeBranchId(event_runtime_branch_id.value());
  }

  const auto& containing_tri_ids = kin_del.getCrossingData().getContainingTriIds();
  const auto intersection_debug_data = kin_del.getCrossingIntersectionDebugData();
  const std::unordered_map<size_t, size_t> site_input_branch_labels
    = buildSiteInputBranchLabels(kin_del, graph, occurrence_time);
  const std::unordered_set<size_t> live_strand_ids = collectLiveStrandIds(graph);
  const std::unordered_map<size_t, glm::dvec3> live_site_world_positions
    = buildSiteWorldPositions(kin_del, occurrence_time, live_strand_ids);
  const std::unordered_map<size_t, glm::dvec3> voronoi_vertex_world_positions
    = buildVoronoiVertexWorldPositions(kin_del, graph, occurrence_time);
  const std::optional<std::filesystem::path>& output_root = kin_del.getVisualDebugOutputRoot();
  const std::unordered_set<size_t> active_runtime_branches
    = collectActiveRuntimeBranches(kin_del, live_strand_ids);

  const bool per_branch_svgs = active_runtime_branches.size() > 1;

  auto write_for_runtime_branch = [&](size_t runtime_branch_id)
  {
    const std::unordered_set<size_t> branch_strands
      = collectStrandIdsForRuntimeBranch(kin_del, live_strand_ids, runtime_branch_id);
    const std::unordered_map<size_t, glm::dvec3> branch_site_world_positions
      = buildSiteWorldPositions(kin_del, occurrence_time, branch_strands);
    const auto [points, positioned_strands]
      = buildVisualDebugStrandPositions(kin_del, graph, occurrence_time, branch_strands);
    if (positioned_strands.empty())
    {
      return;
    }
    const std::string filename = visualDebugSvgRelativePath(
      occurrence_time, phase, event_descriptor, runtime_branch_id, output_root);
    writeVisualDebugSvgFile(filename, points, graph, kin_del, containing_tri_ids, intersection_debug_data, highlight,
      &branch_strands, positioned_strands, site_input_branch_labels, branch_site_world_positions,
      voronoi_vertex_world_positions, separation_offset_segments, seam_outlines);
  };

  if (!per_branch_svgs)
  {
    std::optional<size_t> runtime_branch_id = event_runtime_branch_id;
    if (!runtime_branch_id.has_value())
    {
      runtime_branch_id = inferEventRuntimeBranchFromHighlight(graph, kin_del, highlight);
    }
    if (!runtime_branch_id.has_value() && active_runtime_branches.size() == 1)
    {
      runtime_branch_id = *active_runtime_branches.begin();
    }

    if (runtime_branch_id.has_value())
    {
      write_for_runtime_branch(runtime_branch_id.value());
      return;
    }

    const auto [points, positioned_strands]
      = buildVisualDebugStrandPositions(kin_del, graph, occurrence_time, live_strand_ids);
    if (positioned_strands.empty())
    {
      return;
    }
    const std::string filename
      = visualDebugSvgRelativePath(occurrence_time, phase, event_descriptor, std::nullopt, output_root);
    writeVisualDebugSvgFile(filename, points, graph, kin_del, containing_tri_ids, intersection_debug_data, highlight,
      nullptr, positioned_strands, site_input_branch_labels, live_site_world_positions, voronoi_vertex_world_positions,
      separation_offset_segments, seam_outlines);
    return;
  }

  std::optional<size_t> runtime_branch_id = event_runtime_branch_id;
  if (!runtime_branch_id.has_value())
  {
    runtime_branch_id = inferEventRuntimeBranchFromHighlight(graph, kin_del, highlight);
  }

  if (runtime_branch_id.has_value())
  {
    write_for_runtime_branch(runtime_branch_id.value());
    return;
  }

  for (size_t branch_id : active_runtime_branches)
  {
    write_for_runtime_branch(branch_id);
  }
}

} // namespace kinDS

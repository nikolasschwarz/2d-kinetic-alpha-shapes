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

  runtime_branches.insert(runtime_branch_map[strand_id]);
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

std::string visualDebugSvgRelativePath(double occurrence_time, const char* phase, const std::string& event_descriptor,
  std::optional<size_t> runtime_branch_id, const std::optional<std::filesystem::path>& output_root)
{
  const std::string basename = "t" + std::to_string(occurrence_time) + "_segmentbuilder_" + phase + "_"
    + event_descriptor + ".svg";
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
      runtime_branches.insert(runtime_branch_map[strand_id]);
    }
  }
  return runtime_branches;
}

std::unordered_set<size_t> collectStrandIdsForRuntimeBranch(
  const KineticDelaunay& kin_del, const std::unordered_set<size_t>& live_strand_ids, size_t runtime_branch_id)
{
  std::unordered_set<size_t> strand_ids;
  const auto& runtime_branch_map = kin_del.getRuntimeBranchMap();
  for (size_t strand_id : live_strand_ids)
  {
    if (strand_id < runtime_branch_map.size() && runtime_branch_map[strand_id] == runtime_branch_id)
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

void writeVisualDebugSvgFile(const std::string& relative_path, const std::vector<glm::dvec2>& points,
  const HalfEdgeDelaunayGraph& graph, KineticDelaunay& kin_del,
  const std::vector<size_t>& containing_tri_ids,
  const std::vector<HalfEdgeDelaunayGraphToSVG::IntersectionDebugInfo>& intersection_debug_data,
  const VisualDebugHighlight& highlight, const std::unordered_set<size_t>* highlighted_strands,
  const std::unordered_set<size_t>& positioned_strands,
  const std::unordered_map<size_t, size_t>& site_input_branch_labels)
{
  const std::filesystem::path filepath(relative_path);
  if (filepath.has_parent_path())
  {
    std::filesystem::create_directories(filepath.parent_path());
  }
  HalfEdgeDelaunayGraphToSVG::write(points, graph, relative_path, 0.1, &kin_del.getFacesInside(), true,
    &containing_tri_ids, &intersection_debug_data, &highlight, highlighted_strands, &site_input_branch_labels,
    &positioned_strands);
}
} // namespace

void writeSegmentBuilderVisualDebugSvg(bool visual_debug, KineticDelaunay& kin_del, const HalfEdgeDelaunayGraph& graph,
  double occurrence_time, const char* phase, const std::string& event_descriptor,
  const VisualDebugHighlight& highlight, std::optional<size_t> event_runtime_branch_id)
{
  if (!visual_debug)
  {
    return;
  }

  const auto& containing_tri_ids = kin_del.getCrossingData().getContainingTriIds();
  const auto intersection_debug_data = kin_del.getCrossingIntersectionDebugData();
  const std::unordered_map<size_t, size_t> site_input_branch_labels
    = buildSiteInputBranchLabels(kin_del, graph, occurrence_time);
  const std::unordered_set<size_t> live_strand_ids = collectLiveStrandIds(graph);
  const std::optional<std::filesystem::path>& output_root = kin_del.getVisualDebugOutputRoot();
  const std::unordered_set<size_t> active_runtime_branches
    = collectActiveRuntimeBranches(kin_del, live_strand_ids);

  const bool per_branch_svgs = active_runtime_branches.size() > 1;

  auto write_for_runtime_branch = [&](size_t runtime_branch_id)
  {
    const std::unordered_set<size_t> branch_strands
      = collectStrandIdsForRuntimeBranch(kin_del, live_strand_ids, runtime_branch_id);
    const auto [points, positioned_strands]
      = buildVisualDebugStrandPositions(kin_del, graph, occurrence_time, branch_strands);
    if (positioned_strands.empty())
    {
      return;
    }
    const std::string filename = visualDebugSvgRelativePath(
      occurrence_time, phase, event_descriptor, runtime_branch_id, output_root);
    writeVisualDebugSvgFile(filename, points, graph, kin_del, containing_tri_ids, intersection_debug_data, highlight,
      &branch_strands, positioned_strands, site_input_branch_labels);
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
      nullptr, positioned_strands, site_input_branch_labels);
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

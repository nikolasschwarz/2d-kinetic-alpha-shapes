#include "SegmentBuilderVisualDebug.hpp"

#include "HalfEdgeDelaunayGraphToSVG.hpp"
#include "KineticDelaunay.hpp"
#include "KineticDelaunayCrossingEvent.hpp"

#include <cmath>
#include <filesystem>
#include <unordered_map>
#include <unordered_set>

namespace kinDS
{
namespace
{
std::unordered_map<size_t, size_t> buildSiteRuntimeBranchLabelsIfDiff(
  KineticDelaunay& kin_del, const HalfEdgeDelaunayGraph& graph, double occurrence_time)
{
  std::unordered_map<size_t, size_t> labels;
  const auto& branch_indices = kin_del.getStrandTree().getBranchIndices();
  const auto& component_map = kin_del.component_data.component_map;
  const size_t vertex_count = graph.getVertexCount();

  size_t section_index = static_cast<size_t>(std::ceil(occurrence_time));
  const size_t height = kin_del.getStrandTree().getHeight();
  if (height > 0 && section_index >= height)
  {
    section_index = height - 1;
  }

  for (size_t strand_id = 0; strand_id < vertex_count; ++strand_id)
  {
    if (kin_del.isDummyBoundary(strand_id))
    {
      continue;
    }
    if (strand_id >= component_map.size())
    {
      continue;
    }

    const size_t runtime_branch = component_map[strand_id];
    size_t data_structure_branch = runtime_branch;
    if (strand_id < branch_indices.size() && section_index < branch_indices[strand_id].size())
    {
      data_structure_branch = branch_indices[strand_id][section_index];
    }

    if (data_structure_branch != runtime_branch)
    {
      labels[strand_id] = runtime_branch;
    }
  }

  return labels;
}

void noteBranchForSite(std::unordered_set<size_t>& branches, const std::vector<size_t>& component_map, size_t strand_id)
{
  if (strand_id < component_map.size())
  {
    branches.insert(component_map[strand_id]);
  }
}

std::optional<size_t> uniqueBranchOrNull(const std::unordered_set<size_t>& branches)
{
  if (branches.size() == 1)
  {
    return *branches.begin();
  }
  return std::nullopt;
}

std::optional<size_t> inferEventBranchFromHighlight(
  const HalfEdgeDelaunayGraph& graph, const KineticDelaunay& kin_del, const VisualDebugHighlight& highlight)
{
  const auto& component_map = kin_del.component_data.component_map;
  std::unordered_set<size_t> branches;

  for (size_t strand_id : highlight.delaunay_vertices)
  {
    noteBranchForSite(branches, component_map, strand_id);
  }

  for (size_t he_id : highlight.directed_half_edges)
  {
    const int origin = graph.halfEdge(he_id).origin;
    if (origin >= 0)
    {
      noteBranchForSite(branches, component_map, static_cast<size_t>(origin));
    }
    const int dest = graph.destination(he_id);
    if (dest >= 0)
    {
      noteBranchForSite(branches, component_map, static_cast<size_t>(dest));
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
        noteBranchForSite(branches, component_map, static_cast<size_t>(v));
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
      noteBranchForSite(branches, component_map, static_cast<size_t>(origin_even));
    }
    if (origin_odd >= 0)
    {
      noteBranchForSite(branches, component_map, static_cast<size_t>(origin_odd));
    }
  }

  return uniqueBranchOrNull(branches);
}

std::string visualDebugSvgRelativePath(
  double occurrence_time, const char* phase, const std::string& event_descriptor, std::optional<size_t> branch_id)
{
  const std::string basename = "t" + std::to_string(occurrence_time) + "_segmentbuilder_" + phase + "_"
    + event_descriptor + ".svg";
  if (!branch_id.has_value())
  {
    return basename;
  }
  return "branch" + std::to_string(branch_id.value()) + "/" + basename;
}

void writeVisualDebugSvgFile(const std::string& relative_path, const std::vector<glm::dvec2>& points,
  const HalfEdgeDelaunayGraph& graph, KineticDelaunay& kin_del,
  const std::vector<size_t>& containing_tri_ids,
  const std::vector<HalfEdgeDelaunayGraphToSVG::IntersectionDebugInfo>& intersection_debug_data,
  const VisualDebugHighlight& highlight, const std::unordered_set<size_t>* component_strands,
  const std::unordered_map<size_t, size_t>& site_runtime_branch_if_diff)
{
  const std::filesystem::path filepath(relative_path);
  if (filepath.has_parent_path())
  {
    std::filesystem::create_directories(filepath.parent_path());
  }
  HalfEdgeDelaunayGraphToSVG::write(points, graph, relative_path, 0.1, &kin_del.getFacesInside(), true,
    &containing_tri_ids, &intersection_debug_data, &highlight, component_strands, &site_runtime_branch_if_diff);
}
} // namespace

void writeSegmentBuilderVisualDebugSvg(bool visual_debug, KineticDelaunay& kin_del, const HalfEdgeDelaunayGraph& graph,
  double occurrence_time, const char* phase, const std::string& event_descriptor,
  const VisualDebugHighlight& highlight, std::optional<size_t> event_branch_id)
{
  if (!visual_debug)
  {
    return;
  }

  const std::vector<glm::dvec2> points = kin_del.getPointsAt(occurrence_time);
  const auto& containing_tri_ids = kin_del.getCrossingData().getContainingTriIds();
  const auto intersection_debug_data = kin_del.getCrossingIntersectionDebugData();
  const auto& components = kin_del.component_data.components;
  const std::unordered_map<size_t, size_t> site_runtime_branch_if_diff
    = buildSiteRuntimeBranchLabelsIfDiff(kin_del, graph, occurrence_time);

  const bool per_branch_svgs = kin_del.isGraphRetriangulatedForComponents() && components.size() > 1;

  if (!per_branch_svgs)
  {
    const std::string filename = visualDebugSvgRelativePath(occurrence_time, phase, event_descriptor, std::nullopt);
    writeVisualDebugSvgFile(filename, points, graph, kin_del, containing_tri_ids, intersection_debug_data, highlight,
      nullptr, site_runtime_branch_if_diff);
    return;
  }

  std::optional<size_t> branch_id = event_branch_id;
  if (!branch_id.has_value())
  {
    branch_id = inferEventBranchFromHighlight(graph, kin_del, highlight);
  }

  if (branch_id.has_value())
  {
    if (branch_id.value() >= components.size())
    {
      return;
    }
    std::unordered_set<size_t> component_strands(components[branch_id.value()].begin(), components[branch_id.value()].end());
    const std::string filename
      = visualDebugSvgRelativePath(occurrence_time, phase, event_descriptor, branch_id);
    writeVisualDebugSvgFile(filename, points, graph, kin_del, containing_tri_ids, intersection_debug_data, highlight,
      &component_strands, site_runtime_branch_if_diff);
    return;
  }

  for (size_t component_id = 0; component_id < components.size(); ++component_id)
  {
    std::unordered_set<size_t> component_strands(components[component_id].begin(), components[component_id].end());
    const std::string filename
      = visualDebugSvgRelativePath(occurrence_time, phase, event_descriptor, component_id);
    writeVisualDebugSvgFile(filename, points, graph, kin_del, containing_tri_ids, intersection_debug_data, highlight,
      &component_strands, site_runtime_branch_if_diff);
  }
}

} // namespace kinDS

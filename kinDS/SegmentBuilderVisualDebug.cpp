#include "SegmentBuilderVisualDebug.hpp"

#include "HalfEdgeDelaunayGraphToSVG.hpp"
#include "KineticDelaunay.hpp"
#include "KineticDelaunayCrossingEvent.hpp"

#include <cmath>
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
} // namespace

void writeSegmentBuilderVisualDebugSvg(bool visual_debug, KineticDelaunay& kin_del, const HalfEdgeDelaunayGraph& graph,
  double occurrence_time, const char* phase, const std::string& event_descriptor,
  const VisualDebugHighlight& highlight)
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
    const std::string filename = "t" + std::to_string(occurrence_time) + "_segmentbuilder_" + phase + "_"
      + event_descriptor + ".svg";
    HalfEdgeDelaunayGraphToSVG::write(points, graph, filename, 0.1, &kin_del.getFacesInside(), true,
      &containing_tri_ids, &intersection_debug_data, &highlight, nullptr, &site_runtime_branch_if_diff);
    return;
  }

  for (size_t component_id = 0; component_id < components.size(); ++component_id)
  {
    std::unordered_set<size_t> component_strands(components[component_id].begin(), components[component_id].end());
    const std::string filename = "t" + std::to_string(occurrence_time) + "_branch" + std::to_string(component_id)
      + "_segmentbuilder_" + phase + "_" + event_descriptor + ".svg";
    HalfEdgeDelaunayGraphToSVG::write(points, graph, filename, 0.1, &kin_del.getFacesInside(), true,
      &containing_tri_ids, &intersection_debug_data, &highlight, &component_strands, &site_runtime_branch_if_diff);
  }
}

} // namespace kinDS

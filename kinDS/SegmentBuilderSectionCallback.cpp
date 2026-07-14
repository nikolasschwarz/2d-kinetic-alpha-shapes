#include "SegmentBuilderSectionCallback.hpp"

#include "KineticDelaunayCrossingEvent.hpp"
#include "SegmentBuilder.hpp"
#include "SegmentBuilderVisualDebug.hpp"

namespace kinDS
{
void SegmentBuilderSectionCallback::beforeEvent(KineticDelaunay::Event& e)
{
  auto* section = dynamic_cast<KineticDelaunay::SectionEvent*>(&e);
  if (!section)
  {
    return;
  }
  SegmentBuilder::ScopedMetadataCallbackPhase callback_phase(segment_builder_, "before");
  const size_t index = section->section_id;

  const double t = static_cast<double>(index);
  auto& graph = segment_builder_.kin_del.getGraph();
  writeSegmentBuilderVisualDebugSvg(segment_builder_.visual_debug, segment_builder_.kin_del, graph, t, "before",
    "section_" + std::to_string(index),
    VisualDebugHighlight::forSectionBoundary(
      graph, [&](size_t even_he) { return segment_builder_.kin_del.isOnComponentBoundary(even_he); }));

  segment_builder_.advanceBoundaryMeshes(t);

  const size_t live_edge_count = graph.liveDelaunayEdgeCount();
  segment_builder_.parallel_for(live_edge_count,
    [&](size_t live_index)
  {
    const size_t i = graph.liveDelaunayEdgeId(live_index);
    // use the origin of the half edge to obtain the correct component
    auto vertex = graph.halfEdge(i).origin;

    // fall back for infinite vertices
    if (vertex == -1)
    {
      vertex = graph.destination(i);
    }
    size_t component_index = segment_builder_.kin_del.component_data.component_map[vertex];
    auto& boundary_points = segment_builder_.kin_del.component_data.component_boundaries[component_index][0];

    segment_builder_.finishMesh(i, t, boundary_points, SegmentBuilder::BoundaryEventType::Section,
      SegmentBuilder::BoundarySegmentAction::SegmentCompleted);

    // Advance boundary-interval meshes on boundary Delaunay edges using the same interval decomposition as init():
    // [null, first], [k, k+1], [last, null].
    if (!segment_builder_.kin_del.isOnComponentBoundary(i))
    {
      return;
    }
    const size_t d_edge_id = i / 2;
    const auto& d_intersections = segment_builder_.kin_del.getCrossingData().delaunay_edge_intersections[d_edge_id];
    if (d_intersections.empty())
    {
      return;
    }

    const std::vector<KineticDelaunay::CrossingData::EdgeIntersectionRef> refs
      = segment_builder_.getBoundaryIntersectionsInBoundaryOrder(d_edge_id);

    {
      const size_t first_cell
        = segment_builder_.determineVoronoiCellForBoundaryIntersectionInterval(d_edge_id, std::nullopt, refs.front());
      segment_builder_.finishMeshFromIntersections(
        first_cell, t, std::nullopt, refs.front(), SegmentBuilder::BoundaryEventType::Section,
        SegmentBuilder::BoundarySegmentAction::SegmentCompleted);
    }
    for (size_t k = 0; k + 1 < refs.size(); ++k)
    {
      const size_t mid_cell
        = segment_builder_.determineVoronoiCellForBoundaryIntersectionInterval(d_edge_id, refs[k], refs[k + 1]);
      if (mid_cell == static_cast<size_t>(-1))
      {
        continue;
      }
      segment_builder_.finishMeshFromIntersections(
        mid_cell, t, refs[k], refs[k + 1], SegmentBuilder::BoundaryEventType::Section,
        SegmentBuilder::BoundarySegmentAction::SegmentCompleted);
    }
    {
      const size_t last_cell
        = segment_builder_.determineVoronoiCellForBoundaryIntersectionInterval(d_edge_id, refs.back(), std::nullopt);
      segment_builder_.finishMeshFromIntersections(
        last_cell, t, refs.back(), std::nullopt, SegmentBuilder::BoundaryEventType::Section,
        SegmentBuilder::BoundarySegmentAction::SegmentCompleted);
    }
  });

  segment_builder_.createClosingCapsForInputBranchesFinishingAtSection(t);

  for (size_t input_branch_id : segment_builder_.kin_del.inputBranchesFinishingAtSection(t))
  {
    segment_builder_.addDelaunayTriangulationToBoundaryMesh(t, input_branch_id, true, 0.01);
  }
}

void SegmentBuilderSectionCallback::afterEvent(KineticDelaunay::Event& e)
{
  auto* section = dynamic_cast<KineticDelaunay::SectionEvent*>(&e);
  if (!section)
  {
    return;
  }
  SegmentBuilder::ScopedMetadataCallbackPhase callback_phase(segment_builder_, "after");

  const size_t index = section->section_id;
  auto& graph = segment_builder_.kin_del.getGraph();
  const double t = static_cast<double>(index);
  if (segment_builder_.diagnostics)
  {
    segment_builder_.kin_del.validateAllFaceInsideStatesAtTime(t, "section_event");
    segment_builder_.logDiagnosticsMonitoredFaceInsideState(t, "section_event");
  }
  writeSegmentBuilderVisualDebugSvg(segment_builder_.visual_debug, segment_builder_.kin_del, graph, t, "after",
    "section_" + std::to_string(index),
    VisualDebugHighlight::forSectionBoundary(
      graph, [&](size_t even_he) { return segment_builder_.kin_del.isOnComponentBoundary(even_he); }));
}
} // namespace kinDS


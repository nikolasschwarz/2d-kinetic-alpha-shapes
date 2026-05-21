#include "SegmentBuilderSectionCallback.hpp"

#include "SegmentBuilder.hpp"
#include "KineticDelaunayCrossingEvent.hpp"
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
  const size_t index = section->section_id;

  auto& graph = segment_builder_.kin_del.getGraph();

  const double t = static_cast<double>(index);
  writeSegmentBuilderVisualDebugSvg(segment_builder_.visual_debug, segment_builder_.kin_del, graph, t, "before",
    "section_" + std::to_string(index),
    VisualDebugHighlight::forSectionBoundary(
      graph, [&](size_t even_he) { return segment_builder_.kin_del.isOnComponentBoundary(even_he); }));

  segment_builder_.advanceBoundaryMeshes(index);

  size_t half_edge_count = graph.getHalfEdges().size();
  for (size_t i = 0; i < half_edge_count; i += 2)
  {
    // use the origin of the half edge to obtain the correct component
    auto vertex = graph.getHalfEdges()[i].origin;

    // fall back for infinite vertices
    if (vertex == -1)
    {
      vertex = graph.destination(i);
    }
    size_t component_index = segment_builder_.kin_del.component_data.component_map[vertex];
    auto& boundary_points = segment_builder_.kin_del.component_data.component_boundaries[component_index][0];

    segment_builder_.finishMesh(i, index, boundary_points, SegmentBuilder::BoundaryEventType::Section,
      SegmentBuilder::BoundarySegmentAction::SegmentCompleted);

    // Advance boundary-interval meshes on boundary Delaunay edges using the same interval decomposition as init():
    // [null, first], [k, k+1], [last, null].
    if (!segment_builder_.kin_del.isOnComponentBoundary(i))
    {
      continue;
    }
    const size_t d_edge_id = i / 2;
    const auto& d_intersections = segment_builder_.kin_del.getCrossingData().delaunay_edge_intersections[d_edge_id];
    if (d_intersections.empty())
    {
      continue;
    }

    const std::vector<KineticDelaunay::CrossingData::EdgeIntersectionRef> refs
      = segment_builder_.getBoundaryIntersectionsInBoundaryOrder(d_edge_id);

    {
      const size_t first_cell
        = segment_builder_.determineVoronoiCellForBoundaryIntersectionInterval(d_edge_id, std::nullopt, refs.front());
      segment_builder_.finishMeshFromIntersections(
        first_cell, index, std::nullopt, refs.front(), SegmentBuilder::BoundaryEventType::Section,
        SegmentBuilder::BoundarySegmentAction::SegmentCompleted);
    }
    for (size_t k = 0; k + 1 < refs.size(); ++k)
    {
      const size_t mid_cell
        = segment_builder_.determineVoronoiCellForBoundaryIntersectionInterval(d_edge_id, refs[k], refs[k + 1]);
      segment_builder_.finishMeshFromIntersections(
        mid_cell, index, refs[k], refs[k + 1], SegmentBuilder::BoundaryEventType::Section,
        SegmentBuilder::BoundarySegmentAction::SegmentCompleted);
    }
    {
      const size_t last_cell
        = segment_builder_.determineVoronoiCellForBoundaryIntersectionInterval(d_edge_id, refs.back(), std::nullopt);
      segment_builder_.finishMeshFromIntersections(
        last_cell, index, refs.back(), std::nullopt, SegmentBuilder::BoundaryEventType::Section,
        SegmentBuilder::BoundarySegmentAction::SegmentCompleted);
    }
  }

}

void SegmentBuilderSectionCallback::afterEvent(KineticDelaunay::Event& e)
{
  auto* section = dynamic_cast<KineticDelaunay::SectionEvent*>(&e);
  if (!section)
  {
    return;
  }

  const size_t index = section->section_id;
  auto& graph = segment_builder_.kin_del.getGraph();
  const double t = static_cast<double>(index);
  writeSegmentBuilderVisualDebugSvg(segment_builder_.visual_debug, segment_builder_.kin_del, graph, t, "after",
    "section_" + std::to_string(index),
    VisualDebugHighlight::forSectionBoundary(
      graph, [&](size_t even_he) { return segment_builder_.kin_del.isOnComponentBoundary(even_he); }));
}
} // namespace kinDS


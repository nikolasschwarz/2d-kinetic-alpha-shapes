#include "SegmentBuilderSectionCallback.hpp"

#include "SegmentBuilder.hpp"
#include "KineticDelaunayCrossingEvent.hpp"

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

  // Check if we need to insert a subdivision before handling this event
  while (segment_builder_.subdivision_index < segment_builder_.subdivisions.size()
    && segment_builder_.subdivisions[segment_builder_.subdivision_index].second <= index)
  {
    segment_builder_.insertSubdivision(segment_builder_.subdivisions[segment_builder_.subdivision_index].first,
      segment_builder_.subdivisions[segment_builder_.subdivision_index].second);
    segment_builder_.subdivision_index++;
  }

  auto& graph = segment_builder_.kin_del.getGraph();

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

    segment_builder_.finishMesh(i, index, boundary_points);
  }

  // Visual debug: export SVG with Voronoi vertex labels at section boundaries.
  if (segment_builder_.visual_debug)
  {
    double t = static_cast<double>(index);
    std::vector<glm::dvec2> points = segment_builder_.kin_del.getPointsAt(t);

    std::string filename
      = "t" + std::to_string(t) + "_segmentbuilder_between_section_" + std::to_string(index) + ".svg";
    const auto& containing_tri_ids = segment_builder_.kin_del.getCrossingData().getContainingTriIds();
    auto intersection_debug_data = segment_builder_.kin_del.getCrossingIntersectionDebugData();
    HalfEdgeDelaunayGraphToSVG::write(points, graph, filename, 0.1, &segment_builder_.kin_del.getFacesInside(), true,
      &containing_tri_ids, &intersection_debug_data);
    KINDS_INFO("SegmentBuilder wrote SVG: " << filename);
  }
}

void SegmentBuilderSectionCallback::afterEvent(KineticDelaunay::Event& e)
{
  (void)e;
  // No SegmentBuilder hook needed after section processing currently.
}
} // namespace kinDS


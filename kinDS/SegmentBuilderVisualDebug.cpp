#include "SegmentBuilderVisualDebug.hpp"

#include "HalfEdgeDelaunayGraphToSVG.hpp"
#include "KineticDelaunay.hpp"
#include "KineticDelaunayCrossingEvent.hpp"

namespace kinDS
{
void writeSegmentBuilderVisualDebugSvg(bool visual_debug, KineticDelaunay& kin_del, const HalfEdgeDelaunayGraph& graph,
  double occurrence_time, const char* phase, const std::string& event_descriptor,
  const VisualDebugHighlight& highlight)
{
  if (!visual_debug)
  {
    return;
  }

  const std::vector<glm::dvec2> points = kin_del.getPointsAt(occurrence_time);
  const std::string filename = "t" + std::to_string(occurrence_time) + "_segmentbuilder_" + phase + "_"
    + event_descriptor + ".svg";
  const auto& containing_tri_ids = kin_del.getCrossingData().getContainingTriIds();
  const auto intersection_debug_data = kin_del.getCrossingIntersectionDebugData();
  HalfEdgeDelaunayGraphToSVG::write(points, graph, filename, 0.1, &kin_del.getFacesInside(), true,
    &containing_tri_ids, &intersection_debug_data, &highlight);
}

} // namespace kinDS

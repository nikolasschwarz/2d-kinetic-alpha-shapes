#pragma once

#include "HalfEdgeDelaunayGraph.hpp"
#include "VisualDebugHighlight.hpp"

#include <string>

namespace kinDS
{
class KineticDelaunay;

/// Writes `t{time}_segmentbuilder_{phase}_{event_descriptor}.svg` when @p visual_debug is true.
void writeSegmentBuilderVisualDebugSvg(bool visual_debug, KineticDelaunay& kin_del, const HalfEdgeDelaunayGraph& graph,
  double occurrence_time, const char* phase, const std::string& event_descriptor,
  const VisualDebugHighlight& highlight);

} // namespace kinDS

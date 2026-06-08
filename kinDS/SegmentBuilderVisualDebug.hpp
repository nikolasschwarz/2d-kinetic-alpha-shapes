#pragma once

#include "HalfEdgeDelaunayGraph.hpp"
#include "VisualDebugHighlight.hpp"

#include <string>

namespace kinDS
{
class KineticDelaunay;

/// Writes one SVG per connected component (runtime branch) when @p visual_debug is true:
/// `t{time}_branch{id}_segmentbuilder_{phase}_{event_descriptor}.svg`.
void writeSegmentBuilderVisualDebugSvg(bool visual_debug, KineticDelaunay& kin_del, const HalfEdgeDelaunayGraph& graph,
  double occurrence_time, const char* phase, const std::string& event_descriptor,
  const VisualDebugHighlight& highlight);

} // namespace kinDS

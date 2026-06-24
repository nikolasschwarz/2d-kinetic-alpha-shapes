#pragma once

#include "HalfEdgeDelaunayGraph.hpp"
#include "VisualDebugHighlight.hpp"

#include <string>

namespace kinDS
{
class KineticDelaunay;

/// Writes debug SVG snapshots when @p visual_debug is true.
/// After a runtime component split, per-branch files (`t{time}_branch{id}_...`) are only written once the
/// Delaunay graph has been retriangulated at the next section; until then a single unsplit snapshot is exported.
void writeSegmentBuilderVisualDebugSvg(bool visual_debug, KineticDelaunay& kin_del, const HalfEdgeDelaunayGraph& graph,
  double occurrence_time, const char* phase, const std::string& event_descriptor,
  const VisualDebugHighlight& highlight);

} // namespace kinDS

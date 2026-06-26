#pragma once

#include "HalfEdgeDelaunayGraph.hpp"
#include "VisualDebugHighlight.hpp"

#include <optional>
#include <string>

namespace kinDS
{
class KineticDelaunay;

/// Writes debug SVG snapshots when @p visual_debug is true.
/// After a runtime component split, exports go under `branch{id}/` and only the branch affected by the event is
/// written when @p event_branch_id is set or can be inferred from @p highlight. Section-wide events with no single
/// branch export one file per branch folder.
void writeSegmentBuilderVisualDebugSvg(bool visual_debug, KineticDelaunay& kin_del, const HalfEdgeDelaunayGraph& graph,
  double occurrence_time, const char* phase, const std::string& event_descriptor,
  const VisualDebugHighlight& highlight, std::optional<size_t> event_branch_id = std::nullopt);

} // namespace kinDS

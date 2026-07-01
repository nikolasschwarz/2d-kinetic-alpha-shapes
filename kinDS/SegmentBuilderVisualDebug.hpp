#pragma once

#include "HalfEdgeDelaunayGraph.hpp"
#include "VisualDebugHighlight.hpp"

#include <optional>
#include <string>

namespace kinDS
{
class KineticDelaunay;

/// Writes debug SVG snapshots when @p visual_debug is true.
/// Exports go under `branch{runtime_branch_id}/` using @ref KineticDelaunay::getRuntimeBranchIdForStrand (not
/// inside-face @ref ComponentData::component_map). When @p event_runtime_branch_id is set or can be inferred from
/// @p highlight, only that runtime branch is written; otherwise one file is written per active runtime branch.
void writeSegmentBuilderVisualDebugSvg(bool visual_debug, KineticDelaunay& kin_del, const HalfEdgeDelaunayGraph& graph,
  double occurrence_time, const char* phase, const std::string& event_descriptor,
  const VisualDebugHighlight& highlight, std::optional<size_t> event_runtime_branch_id = std::nullopt);

} // namespace kinDS

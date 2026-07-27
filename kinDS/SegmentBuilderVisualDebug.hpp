#pragma once

#include "HalfEdgeDelaunayGraph.hpp"
#include "HalfEdgeDelaunayGraphToSVG.hpp"
#include "VisualDebugHighlight.hpp"
#include "DebugExportFormatting.hpp"

#include <cstddef>
#include <optional>
#include <string>
#include <vector>

namespace kinDS
{
class KineticDelaunay;

/// Folder for debug SVGs when runtime-branch resolution fails (all live strands are exported).
inline constexpr const char* kVisualDebugUnresolvedBranchFolder = "branchX";

/// Writes debug SVG snapshots when @p visual_debug is true.
/// Exports go under `branch{runtime_branch_id}/` when a unique runtime branch is resolved (preferred id, highlight
/// inference, or a single active branch). When no unique branch is known, the SVG is written once under
/// @ref kVisualDebugUnresolvedBranchFolder with all live strands.
///
/// By default, pending split-off child runtime branches collapse into the parent (unsplit) folder. When
/// @ref KineticDelaunay::visualDebugSeparatePendingSplits is enabled, pending children keep their own folders and
/// strand subsets from the moment the radius event notes the split.
///
/// Filename times use @ref formatDebugExportTimeToken / @ref kDebugExportTimePrecision.
/// Occurrence time is the leading `t…` token; when @p creation_time is set, another `t…` token is
/// appended before the extension so duplicate events at the same occurrence time stay distinct.
///
/// If @p explicit_runtime_branch_ids is non-null and non-empty, SVGs are written only for those runtime branch ids
/// (duplicates allowed when endpoints span multiple branches). Failure to position any of them falls back to
/// @ref kVisualDebugUnresolvedBranchFolder.
void writeSegmentBuilderVisualDebugSvg(bool visual_debug, KineticDelaunay& kin_del, const HalfEdgeDelaunayGraph& graph,
  double occurrence_time, const char* phase, const std::string& event_descriptor,
  const VisualDebugHighlight& highlight, std::optional<size_t> event_runtime_branch_id = std::nullopt,
  const std::vector<HalfEdgeDelaunayGraphToSVG::SeparationOffsetSegment>* separation_offset_segments = nullptr,
  const std::vector<std::vector<glm::dvec2>>* seam_outlines = nullptr,
  const std::vector<size_t>* explicit_runtime_branch_ids = nullptr,
  std::optional<double> creation_time = std::nullopt);

} // namespace kinDS

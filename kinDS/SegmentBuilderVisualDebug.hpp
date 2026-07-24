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
/// inference, or a single active branch). When no unique branch is known (e.g. section events), one SVG is written
/// per active runtime branch. @ref kVisualDebugUnresolvedBranchFolder is used only when a preferred/explicit branch
/// cannot be positioned, or when no active branch has strands.
///
/// By default, pending split-off child runtime branches collapse into the parent (unsplit) folder. When
/// @ref KineticDelaunay::visualDebugSeparatePendingSplits is enabled, pending children keep their own folders and
/// strand subsets from the moment the radius event notes the split.
///
/// Filename times use @ref formatDebugExportTimeToken / @ref kDebugExportTimePrecision.
///
/// If @p explicit_runtime_branch_ids is non-null and non-empty, SVGs are written only for those runtime branch ids
/// (no fan-out to every other active branch).
void writeSegmentBuilderVisualDebugSvg(bool visual_debug, KineticDelaunay& kin_del, const HalfEdgeDelaunayGraph& graph,
  double occurrence_time, const char* phase, const std::string& event_descriptor,
  const VisualDebugHighlight& highlight, std::optional<size_t> event_runtime_branch_id = std::nullopt,
  const std::vector<HalfEdgeDelaunayGraphToSVG::SeparationOffsetSegment>* separation_offset_segments = nullptr,
  const std::vector<std::vector<glm::dvec2>>* seam_outlines = nullptr,
  const std::vector<size_t>* explicit_runtime_branch_ids = nullptr);

} // namespace kinDS

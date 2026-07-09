#include "SegmentBuilderSeparationCallback.hpp"

#include "KineticDelaunaySeparationEvent.hpp"
#include "SegmentBuilder.hpp"
#include "SegmentBuilderVisualDebug.hpp"
#include "VisualDebugHighlight.hpp"

using namespace kinDS;

namespace
{
std::string separationEventDescriptor(size_t parent_component_id, size_t iteration)
{
  return "separation_parent" + std::to_string(parent_component_id) + "_iter" + std::to_string(iteration);
}

std::optional<size_t> runtimeBranchForPendingSplitParent(const KineticDelaunay& kin_del, size_t parent_component_id)
{
  const std::optional<PendingBranchSplit> split = kin_del.getPendingBranchSplit(parent_component_id);
  if (!split.has_value() || split->split_component_ids.empty())
  {
    return std::nullopt;
  }

  const size_t retained_component_id = split->split_component_ids.front();
  if (retained_component_id >= kin_del.component_data.components.size())
  {
    return std::nullopt;
  }

  const std::vector<size_t>& retained_strands = kin_del.component_data.components[retained_component_id];
  for (size_t strand_id : retained_strands)
  {
    if (!kin_del.isDummyBoundary(strand_id))
    {
      return kin_del.getRuntimeBranchIdForStrand(strand_id);
    }
  }

  return std::nullopt;
}
} // namespace

void SegmentBuilderSeparationCallback::writeSeparationVisualDebugSvg(
  const KineticDelaunay::SeparationEvent& separation, const char* phase) const
{
  KineticDelaunay& kin_del = segment_builder_.kin_del;
  const HalfEdgeDelaunayGraph& graph = kin_del.getGraph();
  const size_t parent_component_id = separation.parent_component_id;
  const double t = separation.occurrence_time;

  const std::optional<PendingBranchSplit> split = kin_del.getPendingBranchSplit(parent_component_id);
  const size_t iteration = split ? split->separation_iteration : 0;
  const VisualDebugHighlight highlight = kin_del.buildSeparationRecomputeHighlight(parent_component_id);
  const std::vector<HalfEdgeDelaunayGraphToSVG::SeparationOffsetSegment> offset_segments
    = kin_del.collectSeparationOffsetSegments(parent_component_id, t);

  const std::vector<std::vector<BoundaryPoint>> seam_outline_points
    = kin_del.collectPendingSplitBranchOutlines(parent_component_id, t);
  std::vector<std::vector<glm::dvec2>> seam_outlines;
  seam_outlines.reserve(seam_outline_points.size());
  for (const std::vector<BoundaryPoint>& outline : seam_outline_points)
  {
    std::vector<glm::dvec2> loop;
    loop.reserve(outline.size());
    for (const BoundaryPoint& bp : outline)
    {
      loop.push_back(bp.p);
    }
    seam_outlines.push_back(std::move(loop));
  }

  writeSegmentBuilderVisualDebugSvg(segment_builder_.visual_debug, kin_del, graph, t, phase,
    separationEventDescriptor(parent_component_id, iteration), highlight,
    runtimeBranchForPendingSplitParent(kin_del, parent_component_id), &offset_segments, &seam_outlines);
}

void SegmentBuilderSeparationCallback::beforeEvent(KineticDelaunay::Event& e)
{
  auto* separation = dynamic_cast<KineticDelaunay::SeparationEvent*>(&e);
  if (!separation)
  {
    return;
  }
  SegmentBuilder::ScopedMetadataCallbackPhase callback_phase(segment_builder_, "before");

  writeSeparationVisualDebugSvg(*separation, "before");
}

void SegmentBuilderSeparationCallback::afterEvent(KineticDelaunay::Event& e)
{
  auto* separation = dynamic_cast<KineticDelaunay::SeparationEvent*>(&e);
  if (!separation)
  {
    return;
  }
  SegmentBuilder::ScopedMetadataCallbackPhase callback_phase(segment_builder_, "after");

  writeSeparationVisualDebugSvg(*separation, "after");
}

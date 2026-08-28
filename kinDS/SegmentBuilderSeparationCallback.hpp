#pragma once

#include "KineticDelaunay.hpp"
#include "KineticDelaunaySeparationEvent.hpp"

#include <cstddef>
#include <optional>
#include <vector>

namespace kinDS
{
class SegmentBuilder;

class SegmentBuilderSeparationCallback final : public KineticDelaunay::EventCallback
{
 public:
  explicit SegmentBuilderSeparationCallback(SegmentBuilder& segment_builder)
    : segment_builder_(segment_builder)
  {
  }

  void beforeEvent(KineticDelaunay::Event& e) override;
  void afterEvent(KineticDelaunay::Event& e) override;

  /// Visual-debug SVG for a pending/active separation (SeparationEvent or infinitesimal activation).
  void writeSeparationVisualDebugSvg(size_t parent_component_id, EventTime occurrence_time, const char* phase,
    std::optional<EventTime> creation_time = std::nullopt) const;

 private:
  SegmentBuilder& segment_builder_;
  /// Cached across before/after when the graph cut clears the pending split mid-event.
  mutable size_t cached_split_parent_component_id_ = static_cast<size_t>(-1);
  mutable size_t cached_split_parent_branch_id_ = static_cast<size_t>(-1);
  mutable std::vector<size_t> cached_split_child_branch_ids_;

  void writeSeparationVisualDebugSvg(const KineticDelaunay::SeparationEvent& separation, const char* phase) const;
};
} // namespace kinDS

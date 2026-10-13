#pragma once

#include "KineticDelaunay.hpp"
#include "KineticDelaunayRadiusEvent.hpp"

#include <array>
#include <cstddef>
#include <vector>

namespace kinDS
{
class SegmentBuilder;

class SegmentBuilderRadiusCallback final : public KineticDelaunay::EventCallback
{
 public:
  explicit SegmentBuilderRadiusCallback(SegmentBuilder& segment_builder)
    : segment_builder_(segment_builder)
  {
  }

  void beforeEvent(KineticDelaunay::Event& e) override;
  void afterEvent(KineticDelaunay::Event& e) override;

 private:
  SegmentBuilder& segment_builder_;

  /// Pre-flip triangle state for radius boundary-transition shift (used between @c beforeEvent and @c afterEvent).
  size_t radius_pre_boundary_edge_count_ = 0;
  size_t radius_pre_face_id_ = 0;
  std::array<size_t, 3> radius_pre_face_he_ {};
  std::array<bool, 3> radius_pre_is_boundary_edge_ {};
  /// Boundary-interval pair indices finished in @c beforeEvent on the triangle's single pre-flip boundary edge
  /// (1→2 candidate pool for @ref SegmentBuilder::maybeQueueRadiusComplementarySplitForExistingMid).
  std::vector<size_t> radius_pre_finished_one_edge_meshlet_ids_;
};
} // namespace kinDS

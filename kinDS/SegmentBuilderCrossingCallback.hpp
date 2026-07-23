#pragma once

#include "KineticDelaunay.hpp"
#include "KineticDelaunayCrossingEvent.hpp"
#include <string>
#include <vector>

namespace kinDS
{
class SegmentBuilder;

class SegmentBuilderCrossingCallback final : public KineticDelaunay::EventCallback
{
 public:
  explicit SegmentBuilderCrossingCallback(SegmentBuilder& segment_builder)
    : segment_builder_(segment_builder)
  {
  }

  void beforeEvent(KineticDelaunay::Event& e) override;
  void afterEvent(KineticDelaunay::Event& e) override;

 private:
  struct CrossingEdgeSnapshotEntry
  {
    size_t voronoi_edge_id = static_cast<size_t>(-1);
    size_t prev_pair_idx = static_cast<size_t>(-1);
    size_t next_pair_idx = static_cast<size_t>(-1);
  };

  std::vector<CrossingEdgeSnapshotEntry> crossing_edge_snapshot_;
  size_t crossing_edge_snapshot_delaunay_edge_id_ = static_cast<size_t>(-1);
  /// @c (prev, next, ...) mesh-pair links on the crossed Delaunay edge, captured in @ref beforeEvent.
  std::string crossing_edge_links_before_;
  SegmentBuilder& segment_builder_;
};
} // namespace kinDS


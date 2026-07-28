#pragma once

#include "KineticDelaunay.hpp"

#include <cmath>

namespace kinDS
{
/**
 * \brief Scheduled mesh event: apply branch separation geometry for a pending component split.
 *
 * Does not change the Delaunay triangulation; triggers @ref KineticDelaunay::EventCallback hooks (typically
 * @ref SegmentBuilderSeparationCallback).
 */
class KineticDelaunay::SeparationEvent final : public KineticDelaunay::Event
{
 public:
  static constexpr uint32_t scheduled_iteration_dispatch_order = 15u;
  /// Runs after @ref RadiusEvent (30) at the same occurrence time when enqueued from @ref notePendingBranchSplit.
  static constexpr uint32_t after_radius_dispatch_order = 35u;

  static double scheduledOccurrenceTime(double split_time)
  {
    return 0.5 * (split_time + std::floor(split_time) + 1.0);
  }

  KineticDelaunay* kd_;
  size_t parent_component_id = static_cast<size_t>(-1);
  double split_time = 0.0;

  SeparationEvent(KineticDelaunay* kd, double occurrence_time, size_t parent_component_id, double split_time,
    double creation_time, uint32_t dispatch_order = scheduled_iteration_dispatch_order)
    : KineticDelaunay::Event(occurrence_time, creation_time, dispatch_order)
    , kd_(kd)
    , parent_component_id(parent_component_id)
    , split_time(split_time)
  {
  }

  void handleEvent() override;
};

class KineticDelaunay::SeparationEventManager final : public KineticDelaunay::EventManager
{
 public:
  explicit SeparationEventManager(KineticDelaunay* kd)
    : kd_(kd)
  {
  }

  /// Separation events are enqueued when a pending branch split is recorded; nothing incremental.
  void computeEvents(double /*t*/, size_t /*event_id*/) override { }

 private:
  KineticDelaunay* kd_;
};
} // namespace kinDS

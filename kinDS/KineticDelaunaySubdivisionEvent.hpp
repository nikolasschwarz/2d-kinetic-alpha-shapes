#pragma once

#include "KineticDelaunay.hpp"

#include <stdexcept>

namespace kinDS
{
/**
 * \brief Scheduled mesh event: insert an extra longitudinal segment on a strand at a given kinetic time.
 *
 * Does not change the Delaunay triangulation; only triggers @ref KineticDelaunay::EventCallback hooks (typically
 * @ref SegmentBuilderSubdivisionCallback).
 */
class KineticDelaunay::SubdivisionEvent final : public KineticDelaunay::Event
{
 public:
  KineticDelaunay* kd_;
  size_t strand_id;

  SubdivisionEvent(KineticDelaunay* kd, double occurrence_time, size_t strand_id, double creation_time)
    : KineticDelaunay::Event(occurrence_time, creation_time, 0u)
    , kd_(kd)
    , strand_id(strand_id)
  {
  }

  void handleEvent() override;
};

class KineticDelaunay::SubdivisionEventManager final : public KineticDelaunay::EventManager
{
 public:
  explicit SubdivisionEventManager(KineticDelaunay* kd)
    : kd_(kd)
  {
  }

  /// Subdivision events are enqueued once in @ref KineticDelaunay::compute (after section scheduling); nothing incremental.
  void computeEvents(double /*t*/, size_t /*event_id*/) override { }

 private:
  KineticDelaunay* kd_;
};
} // namespace kinDS

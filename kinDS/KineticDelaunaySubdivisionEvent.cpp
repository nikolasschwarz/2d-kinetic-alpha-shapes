#include "KineticDelaunaySubdivisionEvent.hpp"

#include "Logger.hpp"

#include <stdexcept>

using namespace kinDS;

void KineticDelaunay::SubdivisionEvent::handleEvent()
{
  auto* kd = kd_;
  if (!kd)
  {
    throw std::runtime_error("SubdivisionEvent has no KineticDelaunay pointer");
  }

  KINDS_DEBUG("SubdivisionEvent at t=" << occurrence_time << " strand_id=" << strand_id << " queue_seq=" << queue_sequence_);

  const auto topology_lock = kd->lockSharedTopology();
  if (!kd->eventInActiveRuntimeScope(*this))
  {
    KINDS_DEBUG("SubdivisionEvent SKIP strand_id=" << strand_id << " not owned by active runtime branch");
    kd->exportRuntimeScopeSkipSvg(*this);
    return;
  }

  auto* handler = kd->subdivision_event_manager_->getCallback();
  if (handler)
  {
    handler->beforeEvent(*this);
  }

  kd->validateSitesInsideConvexHull("SubdivisionEvent", occurrence_time);
  kd->validateDelaunayVoronoiIntersections("SubdivisionEvent", occurrence_time, eventId(), creation_time);

  if (handler)
  {
    handler->afterEvent(*this);
  }
}

#include "KineticDelaunaySeparationEvent.hpp"

#include "Logger.hpp"

#include <stdexcept>

using namespace kinDS;

void KineticDelaunay::SeparationEvent::handleEvent()
{
  auto* kd = kd_;
  if (!kd)
  {
    throw std::runtime_error("SeparationEvent has no KineticDelaunay pointer");
  }

  const std::optional<PendingBranchSplit> split = kd->getPendingBranchSplit(parent_component_id);
  const size_t iteration = split ? split->separation_iteration : 0;
  KINDS_DEBUG("SeparationEvent handle parent_component_id=" << parent_component_id << " event_time=" << occurrence_time
                                                            << " iteration=" << iteration << " split_time=" << split_time
                                                            << " queue_seq=" << queue_sequence_);

  auto* handler = kd->separation_event_manager_->getCallback();
  if (handler)
  {
    handler->beforeEvent(*this);
  }

  kd->handleSeparationEventAtTime(parent_component_id, occurrence_time);

  if (handler)
  {
    handler->afterEvent(*this);
  }

  kd->validateSitesInsideConvexHull("SeparationEvent:afterEvent", occurrence_time);
}

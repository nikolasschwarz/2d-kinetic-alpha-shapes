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
  if (!split.has_value() || split->on_hiatus || split->split_component_ids.size() < 2)
  {
    KINDS_DEBUG("SeparationEvent SKIP parent_component_id=" << parent_component_id << " event_time=" << occurrence_time
                                                            << " split_time=" << split_time
                                                            << " (pending gone, on hiatus, or incomplete)");
    return;
  }

  KINDS_DEBUG("SeparationEvent handle parent_component_id=" << parent_component_id << " event_time=" << occurrence_time
                                                            << " epoch=" << split->infinitesimal_epoch
                                                            << " split_time=" << split_time
                                                            << " queue_seq=" << queue_sequence_);

  auto* handler = kd->separation_event_manager_->getCallback();
  if (handler)
  {
    handler->beforeEvent(*this);
  }

  // apply_cut_now=true: finalize convex cut (or start infinitesimal if seams are no longer convex).
  kd->handleSeparationEventAtTime(parent_component_id, occurrence_time);

  if (handler)
  {
    handler->afterEvent(*this);
  }

  kd->validateSitesInsideConvexHull("SeparationEvent:afterEvent", occurrence_time);
}

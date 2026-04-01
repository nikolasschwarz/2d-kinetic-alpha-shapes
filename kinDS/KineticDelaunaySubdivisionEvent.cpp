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

  auto* handler = kd->subdivision_event_manager_->getCallback();
  if (handler)
  {
    handler->beforeEvent(*this);
    handler->afterEvent(*this);
  }
}

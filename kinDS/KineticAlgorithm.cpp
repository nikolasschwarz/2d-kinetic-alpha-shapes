#include "KineticAlgorithm.hpp"

#include <optional>

using namespace kinDS;

void KineticAlgorithm::clear()
{
  while (!events_.empty())
  {
    events_.pop();
  }
  next_queue_sequence_ = 0;
}

void KineticAlgorithm::processEvents(std::optional<double> end_time, Statistics* statistics)
{
  while (!events_.empty())
  {
    auto event = events_.top();
    if (end_time.has_value() && !(event->occurrence_time.real_time < *end_time))
    {
      // Queue is ordered by occurrence_time; everything left is at/after the stop time.
      clear();
      break;
    }
    events_.pop();
    if (statistics != nullptr)
    {
      Statistics::EventListRow row;
      row.event_id = event->eventId();
      row.occurrence_t = event->occurrence_time.real_time;
      row.occurrence_infinitesimal_t = event->occurrence_time.infinitesimal_time;
      row.type = event->eventType();
      event->populateEventListFields(row);
      statistics->recordEvent(std::move(row));
    }
    event->handleEvent();
  }
}

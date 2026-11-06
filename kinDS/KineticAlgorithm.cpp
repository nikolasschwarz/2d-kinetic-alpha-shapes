#include "KineticAlgorithm.hpp"
#include "Logger.hpp"

#include <exception>
#include <optional>

using namespace kinDS;

void EventScheduler::clear()
{
  while (!events_.empty())
  {
    events_.pop();
  }
  next_queue_sequence_ = 0;
}

std::vector<std::shared_ptr<EventScheduler::Event>> EventScheduler::takeAll()
{
  std::vector<std::shared_ptr<Event>> out;
  out.reserve(events_.size());
  while (!events_.empty())
  {
    out.push_back(events_.top());
    events_.pop();
  }
  return out;
}

void EventScheduler::processEvents(
  EndTimeQuery get_end_time, Statistics* statistics, std::optional<size_t> runtime_branch_id,
  EventScheduler::MeshingErrorHandler on_meshing_error)
{
  last_processed_occurrence_.reset();
  while (!events_.empty())
  {
    auto event = events_.top();
    const std::optional<double> end_time = get_end_time ? get_end_time() : std::nullopt;
    if (end_time.has_value() && !(event->occurrence_time.real_time < *end_time))
    {
      KINDS_INFO("EventScheduler stop at end_time=" << *end_time << " next_event_t=" << event->occurrence_time
                                                    << " discarded=" << events_.size());
      // Queue is ordered by occurrence_time; everything left is at/after the stop time.
      clear();
      break;
    }
    events_.pop();
    last_processed_occurrence_ = event->occurrence_time;
    if (statistics != nullptr)
    {
      Statistics::EventListRow row;
      row.event_id = event->eventId();
      row.occurrence_t = event->occurrence_time.real_time;
      row.occurrence_infinitesimal_t = event->occurrence_time.infinitesimal_time;
      row.type = event->eventType();
      event->populateEventListFields(row);
      if (!row.runtime_branch_id.has_value())
      {
        row.runtime_branch_id = runtime_branch_id;
      }
      statistics->recordEvent(std::move(row));
    }
    try
    {
      event->handleEvent();
    }
    catch (const std::exception& ex)
    {
      if (!on_meshing_error || !on_meshing_error(*event, ex))
      {
        throw;
      }
    }
  }
}

void EventScheduler::processEvents(
  std::optional<double> end_time, Statistics* statistics, std::optional<size_t> runtime_branch_id,
  EventScheduler::MeshingErrorHandler on_meshing_error)
{
  processEvents([end_time]() { return end_time; }, statistics, runtime_branch_id, std::move(on_meshing_error));
}

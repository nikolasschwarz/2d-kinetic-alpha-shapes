#pragma once

#include <cstddef>
#include <cstdint>
#include <memory>
#include <optional>
#include <queue>
#include <tuple>
#include <vector>

namespace kinDS
{
class KineticAlgorithm
{
 public:
  class Event
  {
   public:
    virtual ~Event() = default;

    double occurrence_time;
    double creation_time;
    /// Lower value runs before higher when @p occurrence_time ties (e.g. strand subdivision before flip).
    uint32_t queue_dispatch_order_;
    /// Monotonic id assigned in @ref enqueueEvent for stable ordering when time and dispatch order tie.
    uint64_t queue_sequence_ = 0;

    virtual void handleEvent() = 0;

    double getTime() const { return occurrence_time; }

    void assignQueueSequence(uint64_t sequence) { queue_sequence_ = sequence; }

   protected:
    Event(double occurrence_time, double creation_time, uint32_t queue_dispatch_order = 10u)
      : occurrence_time(occurrence_time)
      , creation_time(creation_time)
      , queue_dispatch_order_(queue_dispatch_order)
    {
    }
  };

  class EventCallback
  {
   public:
    virtual ~EventCallback() = default;
    virtual void beforeEvent(Event& e) { }
    virtual void afterEvent(Event& e) { }
  };

  class EventManager
  {
   public:
    virtual ~EventManager() = default;

    virtual void computeEvents(double t, size_t event_id) = 0;
    void setCallback(EventCallback* callback) { callback_ = callback; }
    EventCallback* getCallback() const { return callback_; }

   private:
    EventCallback* callback_ = nullptr;
  };

 private:
  struct EventPtrCompare
  {
    bool operator()(const std::shared_ptr<Event>& a, const std::shared_ptr<Event>& b) const
    {
      return std::tie(a->occurrence_time, a->queue_dispatch_order_, a->queue_sequence_)
        > std::tie(b->occurrence_time, b->queue_dispatch_order_, b->queue_sequence_);
    }
  };

  using EventQueue = std::priority_queue<std::shared_ptr<Event>, std::vector<std::shared_ptr<Event>>, EventPtrCompare>;
  EventQueue events_;
  uint64_t next_queue_sequence_ = 0;

 public:
  void enqueueEvent(const std::shared_ptr<Event>& event)
  {
    event->assignQueueSequence(next_queue_sequence_++);
    events_.push(event);
  }
  bool empty() const { return events_.empty(); }
  void clear();
  /// Process queued events in time order. If @p end_time is set, discard the remainder once the next
  /// event has @c occurrence_time >= @p end_time (that event is not executed).
  void processEvents(std::optional<double> end_time = std::nullopt);
};
} // namespace kinDS

#pragma once

#include "Statistics.hpp"

#include <cstddef>
#include <cstdint>
#include <memory>
#include <optional>
#include <ostream>
#include <queue>
#include <tuple>
#include <vector>

namespace kinDS
{
/// Kinetic time plus optional virtual / infinitesimal coordinate.
/// Lexicographic order: (real_time, infinitesimal_time).
/// Plain kinetic times use infinitesimal_time == 0.
struct EventTime
{
  double real_time = 0.0;
  double infinitesimal_time = 0.0;

  EventTime() = default;
  EventTime(double real_time, double infinitesimal_time = 0.0)
    : real_time(real_time)
    , infinitesimal_time(infinitesimal_time)
  {
  }

  /// Conservatively allow passing EventTime where only kinetic/real time is needed.
  operator double() const { return real_time; }

  friend bool operator==(const EventTime& a, const EventTime& b)
  {
    return a.real_time == b.real_time && a.infinitesimal_time == b.infinitesimal_time;
  }
  friend bool operator!=(const EventTime& a, const EventTime& b) { return !(a == b); }
  friend bool operator<(const EventTime& a, const EventTime& b)
  {
    return std::tie(a.real_time, a.infinitesimal_time) < std::tie(b.real_time, b.infinitesimal_time);
  }
  friend bool operator>(const EventTime& a, const EventTime& b) { return b < a; }
  friend bool operator<=(const EventTime& a, const EventTime& b) { return !(b < a); }
  friend bool operator>=(const EventTime& a, const EventTime& b) { return !(a < b); }

  friend std::ostream& operator<<(std::ostream& os, const EventTime& t)
  {
    os << t.real_time;
    if (t.infinitesimal_time != 0.0)
    {
      os << "+inf*" << t.infinitesimal_time;
    }
    return os;
  }
};

/// Explicit schedule mode for @ref KineticAlgorithm::EventManager::computeEvents.
/// When set, managers build frozen+virtual site trajectories and enqueue by infinitesimal time
/// with roots in (@ref min_infinitesimal_t, +∞) at frozen real time @c t.
struct InfinitesimalComputeContext
{
  double min_infinitesimal_t = 0.0;
  size_t parent_component_id = static_cast<size_t>(-1);
};

class KineticAlgorithm
{
 public:
  class Event
  {
   public:
    virtual ~Event() = default;

    EventTime occurrence_time;
    EventTime creation_time;
    /// Epoch of the pending split that produced this event; stale after finalize / re-activate.
    uint64_t infinitesimal_epoch_ = 0;
    /// Secondary order when @ref occurrence_time ties (e.g. strand subdivision before flip).
    uint32_t queue_dispatch_order_;
    /// Monotonic id assigned in @ref enqueueEvent for stable ordering when higher keys tie.
    uint64_t queue_sequence_ = 0;

    virtual void handleEvent() = 0;
    virtual KineticEventType eventType() const = 0;

    double getTime() const { return occurrence_time.real_time; }

    void assignQueueSequence(uint64_t sequence) { queue_sequence_ = sequence; }

   protected:
    Event(EventTime occurrence_time, EventTime creation_time, uint32_t queue_dispatch_order = 10u)
      : occurrence_time(occurrence_time)
      , creation_time(creation_time)
      , queue_dispatch_order_(queue_dispatch_order)
    {
    }

    /// Convenience for ordinary kinetic events (infinitesimal_time = 0).
    Event(double occurrence_time, double creation_time, uint32_t queue_dispatch_order = 10u)
      : Event(EventTime(occurrence_time), EventTime(creation_time), queue_dispatch_order)
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

    /// @param infinitesimal When set, use virtual/infinitesimal polynomials and time; otherwise kinetic.
    virtual void computeEvents(double t, size_t event_id,
      std::optional<InfinitesimalComputeContext> infinitesimal = std::nullopt)
      = 0;
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
      return std::tie(a->occurrence_time.real_time, a->occurrence_time.infinitesimal_time, a->queue_dispatch_order_,
               a->queue_sequence_)
        > std::tie(b->occurrence_time.real_time, b->occurrence_time.infinitesimal_time, b->queue_dispatch_order_,
          b->queue_sequence_);
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
  /// When @p statistics is set, each dequeued event is recorded before @c handleEvent.
  void processEvents(std::optional<double> end_time = std::nullopt, Statistics* statistics = nullptr);
};
} // namespace kinDS

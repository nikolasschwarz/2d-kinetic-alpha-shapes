#pragma once

#include "Statistics.hpp"

#include <cstddef>
#include <cstdint>
#include <exception>
#include <functional>
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

/// Priority queue of kinetic events with enqueue-order sequence ids.
class EventScheduler
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
    /// Parent component of that pending split. Epochs are per split, so concurrent separations share epoch values.
    size_t infinitesimal_parent_component_ = static_cast<size_t>(-1);
    /// Secondary order when @ref occurrence_time ties (e.g. strand subdivision before flip).
    uint32_t queue_dispatch_order_;
    /// Monotonic id assigned in @ref enqueueEvent (enqueue order). Used for queue tie-breaks and debug filenames.
    uint64_t queue_sequence_ = 0;

    virtual void handleEvent() = 0;
    virtual KineticEventType eventType() const = 0;
    /// Fill type-specific columns for the meshing event-list CSV (optional fields).
    virtual void populateEventListFields(Statistics::EventListRow& row) const { (void)row; }

    double getTime() const { return occurrence_time.real_time; }

    /// Consecutive event id assigned at enqueue (@ref queue_sequence_).
    uint64_t eventId() const { return queue_sequence_; }

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
  std::optional<EventTime> last_processed_occurrence_;

 public:
  void enqueueEvent(const std::shared_ptr<Event>& event)
  {
    event->assignQueueSequence(next_queue_sequence_++);
    events_.push(event);
  }
  bool empty() const { return events_.empty(); }
  size_t size() const { return events_.size(); }
  std::shared_ptr<Event> peek() const { return events_.empty() ? nullptr : events_.top(); }
  void clear();
  /// Pop every pending event (does not run @c handleEvent). Used by failsafe to keep subdivisions.
  [[nodiscard]] std::vector<std::shared_ptr<Event>> takeAll();
  std::optional<EventTime> lastProcessedOccurrence() const { return last_processed_occurrence_; }
  /// After @c handleEvent throws, return true to keep processing this queue (failsafe already rewound
  /// to the next section). Return false to rethrow. Empty handler rethrows.
  using MeshingErrorHandler = std::function<bool(const Event& failed, const std::exception& ex)>;
  /// Recomputed each loop iteration so a parent remainder can keep going after exclusive end grows.
  using EndTimeQuery = std::function<std::optional<double>()>;
  /// Process queued events in time order. If @p get_end_time returns a cutoff, discard the remainder
  /// once the next event has @c occurrence_time >= that cutoff (that event is not executed).
  /// When @p statistics is set, each dequeued event is recorded before @c handleEvent.
  /// @p runtime_branch_id is stored on event-list rows when the event did not already fill it.
  /// Meshing errors unwind @c handleEvent completely; @p on_meshing_error may rebuild the queue, then
  /// this loop continues as if that event had returned.
  void processEvents(EndTimeQuery get_end_time, Statistics* statistics = nullptr,
    std::optional<size_t> runtime_branch_id = std::nullopt, MeshingErrorHandler on_meshing_error = {});
  void processEvents(std::optional<double> end_time = std::nullopt, Statistics* statistics = nullptr,
    std::optional<size_t> runtime_branch_id = std::nullopt, MeshingErrorHandler on_meshing_error = {});
};

/// Facade that routes enqueue/process to an active @ref EventScheduler (branch-local or root).
class KineticAlgorithm
{
 public:
  using Event = EventScheduler::Event;

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

  KineticAlgorithm() = default;
  KineticAlgorithm(const KineticAlgorithm&) = delete;
  KineticAlgorithm& operator=(const KineticAlgorithm&) = delete;
  ~KineticAlgorithm() { unbindIfOwnedByThisThread(); }

  EventScheduler& rootScheduler() { return root_scheduler_; }
  const EventScheduler& rootScheduler() const { return root_scheduler_; }

  /// Queue this thread enqueues into. Falls back to @ref rootScheduler until @ref setActiveScheduler.
  EventScheduler& activeScheduler() { return *activeSchedulerPtr(); }
  const EventScheduler& activeScheduler() const { return *activeSchedulerPtr(); }

  void bindRootScheduler() { setActiveScheduler(&root_scheduler_); }

  void setActiveScheduler(EventScheduler* scheduler)
  {
    ThreadBinding& binding = threadBinding();
    binding.owner = this;
    binding.scheduler = scheduler != nullptr ? scheduler : &root_scheduler_;
  }

  void unbindIfOwnedByThisThread()
  {
    ThreadBinding& binding = threadBinding();
    if (binding.owner == this)
    {
      binding.owner = nullptr;
      binding.scheduler = nullptr;
    }
  }

  void enqueueEvent(const std::shared_ptr<Event>& event) { activeScheduler().enqueueEvent(event); }
  bool empty() const { return activeScheduler().empty(); }
  void clear() { activeScheduler().clear(); }
  void processEvents(std::optional<double> end_time = std::nullopt, Statistics* statistics = nullptr)
  {
    activeScheduler().processEvents(end_time, statistics);
  }

 private:
  struct ThreadBinding
  {
    const KineticAlgorithm* owner = nullptr;
    EventScheduler* scheduler = nullptr;
  };

  static ThreadBinding& threadBinding()
  {
    thread_local ThreadBinding binding;
    return binding;
  }

  EventScheduler* activeSchedulerPtr() const
  {
    const ThreadBinding& binding = threadBinding();
    if (binding.owner == this && binding.scheduler != nullptr)
    {
      return binding.scheduler;
    }
    return &root_scheduler_;
  }

  mutable EventScheduler root_scheduler_ {};
};
} // namespace kinDS

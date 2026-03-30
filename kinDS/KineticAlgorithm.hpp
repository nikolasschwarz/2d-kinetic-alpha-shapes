#pragma once

#include <cstddef>
#include <memory>
#include <queue>
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

    virtual void handleEvent() = 0;

    double getTime() const { return occurrence_time; }

   protected:
    Event(double occurrence_time, double creation_time)
      : occurrence_time(occurrence_time)
      , creation_time(creation_time)
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
      return a->occurrence_time > b->occurrence_time;
    }
  };

  using EventQueue = std::priority_queue<std::shared_ptr<Event>, std::vector<std::shared_ptr<Event>>, EventPtrCompare>;
  EventQueue events_;

 public:
  void enqueueEvent(const std::shared_ptr<Event>& event) { events_.push(event); }
  bool empty() const { return events_.empty(); }
  void clear();
  void processEvents();
};
} // namespace kinDS

#pragma once

#include "KineticDelaunay.hpp"
#include "ProgressBar.hpp"

namespace kinDS
{
class KineticDelaunay::SectionEvent final : public KineticDelaunay::Event
{
 public:
  size_t section_id;
  glm::dvec2 position;

  SectionEvent(KineticDelaunay* kd, double t, size_t section_id, double creation_time, glm::dvec2 position)
    : KineticDelaunay::Event(kd, t, creation_time)
    , section_id(section_id)
    , position(position)
  {
  }

  void handleEvent() override;
};

class KineticDelaunay::SectionEventManager final : public KineticDelaunay::EventManager
{
 public:
  explicit SectionEventManager(KineticDelaunay* kd)
    : EventManager(kd)
  {
  }

  void computeEvents(double t, size_t event_id) override;

  void updateProgress(size_t section_index);
  void finishProgressIfNeeded(size_t section_index);

  void resetProgress()
  {
    progress_bar_.reset();
    section_count_ = 0;
  }

 private:
  size_t section_count_ = 0;
  std::unique_ptr<ProgressBar> progress_bar_;
};
}
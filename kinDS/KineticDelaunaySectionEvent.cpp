#include "KineticDelaunaySectionEvent.hpp"

#include <stdexcept>

using namespace kinDS;

void KineticDelaunay::SectionEventManager::computeEvents(double t, size_t event_id)
{
  const size_t section_count = kd_->getSectionCount();

  if (event_id != static_cast<size_t>(-1))
  {
    KINDS_ERROR("SectionEventManager::computeEvents only supports full scheduling.");
    return;
  }

  resetProgress();
  section_count_ = section_count;

  // Ensure the section counter starts from 0 for full compute runs.
  kd_->sections_advanced = 0;

  // Clear any queued events from previous runs.
  kd_->kinetic_algorithm_->clear();

  progress_bar_
    = std::make_unique<ProgressBar>(0, section_count, "Computing Kinetic Voronoi Sections", ProgressBar::Display::Absolute);

  // Create events for every section starting at t=0.0.
  for (size_t i = 0; i < section_count; ++i)
  {
    const double time = t + static_cast<double>(i);
    kd_->kinetic_algorithm_->enqueueEvent(std::make_shared<SectionEvent>(kd_, time, i, time, glm::dvec2 { 0.0, 0.0 }));
  }
}

void KineticDelaunay::SectionEventManager::updateProgress(size_t section_index)
{
  if (!progress_bar_)
    return;
  progress_bar_->Update(section_index);
}

void KineticDelaunay::SectionEventManager::finishProgressIfNeeded(size_t section_index)
{
  if (!progress_bar_)
    return;

  if (section_index + 1 >= section_count_)
  {
    progress_bar_->Finish();
  }
}

void KineticDelaunay::SectionEvent::handleEvent()
{
  auto* kd = kd_;
  if (!kd)
  {
    throw std::runtime_error("SectionEvent has no KineticDelaunay pointer");
  }

  const size_t section_index = section_id;
  const size_t section_count = kd->getSectionCount();

  assert(section_index < section_count);

  // In-order processing is required for correct triangulation state updates.
  assert(section_index == kd->sections_advanced);

  auto* event_handler = kd->section_event_manager_->getCallback();

  // This used to happen in KineticDelaunay::compute() / external loops (when section_index != 0).
  if (section_index != 0 && event_handler)
  {
    event_handler->beforeEvent(*this);
  }

  // Retire finished input branches before graph cuts so vertex positions are not queried past strand data.
  kd->retireFinishedInputBranches(static_cast<double>(section_index));

  // Replace the old `advanceOneSection()` logic.
  if (kd->component_data.components.size() > kd->prev_component_count)
  {
    kd->applyPendingComponentGraphSplit(static_cast<double>(section_index));
  }

  kd->precomputeStep(static_cast<double>(section_index));

  kd->section_event_manager_->updateProgress(section_index);
  kd->sections_advanced++;
  kd->section_event_manager_->finishProgressIfNeeded(section_index);

  if (section_index != 0 && event_handler)
  {
    event_handler->afterEvent(*this);
  }
}
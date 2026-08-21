#include "KineticDelaunaySectionEvent.hpp"

#include <algorithm>
#include <stdexcept>
#include <vector>

using namespace kinDS;

void KineticDelaunay::SectionEventManager::computeEvents(double /*t*/, size_t event_id)
{
  const size_t section_count = kd_->getSectionCount();

  if (event_id != static_cast<size_t>(-1))
  {
    KINDS_ERROR("SectionEventManager::computeEvents only supports full scheduling.");
    return;
  }

  resetProgress();

  const size_t start_section = kd_->getStartSection();
  const size_t end_section = kd_->getEndSection();
  if (section_count == 0 || start_section > section_count || end_section < start_section
    || end_section > section_count)
  {
    KINDS_ERROR("SectionEventManager::computeEvents: invalid section range start="
      << start_section << " end=" << end_section << " count=" << section_count);
    return;
  }

  // Section events open the next kinetic interval. The stop/finalize time is @ref getEndSection, so
  // schedule only [start_section, end_section) — the event at end_section must not run.
  section_count_ = end_section > start_section ? end_section - start_section : 0;

  // Align the in-order section counter with the first scheduled section.
  kd_->sections_advanced = start_section;

  // Clear any queued events from previous runs.
  kd_->kinetic_algorithm_->clear();

  if (section_count_ > 0)
  {
    progress_bar_ = std::make_unique<ProgressBar>(
      0, section_count_, "Computing Kinetic Voronoi Sections", ProgressBar::Display::Absolute);
  }

  for (size_t i = start_section; i < end_section; ++i)
  {
    const double time = static_cast<double>(i);
    kd_->kinetic_algorithm_->enqueueEvent(std::make_shared<SectionEvent>(kd_, time, i, time, glm::dvec2 { 0.0, 0.0 }));
  }
}

void KineticDelaunay::SectionEventManager::updateProgress(size_t section_index)
{
  if (!progress_bar_)
    return;
  const size_t start_section = kd_->getStartSection();
  if (section_index < start_section)
  {
    return;
  }
  progress_bar_->Update(section_index - start_section + 1);
}

void KineticDelaunay::SectionEventManager::finishProgressIfNeeded(size_t section_index)
{
  if (!progress_bar_)
    return;

  // Last scheduled section event is end_section - 1.
  if (kd_->getEndSection() > 0 && section_index + 1 >= kd_->getEndSection())
  {
    progress_bar_->Finish();
  }
}

void KineticDelaunay::SectionEventManager::finishProgress()
{
  if (progress_bar_)
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

  // Bootstrap section (normally 0, or @ref getStartSection when starting mid-tree): init already set up meshes.
  const bool is_bootstrap_section = section_index == kd->getStartSection();
  if (!is_bootstrap_section && event_handler)
  {
    event_handler->beforeEvent(*this);
  }

  // Retire finished input branches before graph cuts so vertex positions are not queried past strand data.
  kd->retireFinishedInputBranches(static_cast<double>(section_index));
  kd->post_split_frame_transitions_.expireBeforeHeight(section_index);

  if (kd->collectStatistics())
  {
    const auto [strand_count, branch_count] = kd->countLiveStrandsAndBranches();
    kd->statistics().setSectionTopology(section_index, strand_count, branch_count);
  }

  // Look ahead to height section_index+1: if strands in a live component already sit on multiple input
  // branches there, hold the common frame through that endpoint before scheduling this section's events.
  kd->registerUpcomingPostSplitFrameTransitions(section_index);

  // Force any still-pending splits at the section boundary, one parent *runtime branch* at a time.
  // Hiatus pending splits are not finalizable (seam check / cut would be ill-formed) — leave them alone.
  if (kd->component_data.components.size() > kd->prev_component_count)
  {
    std::vector<size_t> pending_runtime_parents;
    kd->visitPendingBranchSplits(
      [&](size_t, const PendingBranchSplit& split)
      {
        if (split.on_hiatus || split.parent_runtime_branch == KineticDelaunay::RuntimeBranchData::no_branch)
        {
          return;
        }
        if (kd->isPendingRuntimeBranchOnHiatus(split.parent_runtime_branch))
        {
          return;
        }
        if (std::find(pending_runtime_parents.begin(), pending_runtime_parents.end(), split.parent_runtime_branch)
          == pending_runtime_parents.end())
        {
          pending_runtime_parents.push_back(split.parent_runtime_branch);
        }
      });
    for (size_t parent_runtime_branch_id : pending_runtime_parents)
    {
      kd->applyPendingRuntimeBranchSplit(static_cast<double>(section_index), parent_runtime_branch_id);
    }
  }

  kd->precomputeStep(static_cast<double>(section_index));

  kd->section_event_manager_->updateProgress(section_index);
  kd->sections_advanced++;
  kd->section_event_manager_->finishProgressIfNeeded(section_index);

  if (!is_bootstrap_section && event_handler)
  {
    event_handler->afterEvent(*this);
  }
  kd->validateSitesInsideConvexHull("SectionEvent:afterEvent", occurrence_time);
}

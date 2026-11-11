#include "KineticDelaunaySectionEvent.hpp"

#include <algorithm>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

using namespace kinDS;

void KineticDelaunay::SectionEventManager::computeEvents(double /*t*/, size_t event_id,
  std::optional<InfinitesimalComputeContext> /*infinitesimal*/)
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

  // Bind this instance's root queue before clear(). Thread-local scheduler pointers from a previous
  // TreeMesher / KineticDelaunay in this process would otherwise pop freed shared_ptr events.
  kd_->kinetic_algorithm_->bindRootScheduler();
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

void KineticDelaunay::SectionEventManager::updateProgress()
{
  if (!progress_bar_)
  {
    return;
  }
  const size_t start_section = kd_->getStartSection();
  const size_t end_section = kd_->getEndSection();
  const std::optional<size_t> lowest = kd_->lowestNonRetiredSchedulerSection();
  if (!lowest.has_value() || *lowest >= end_section)
  {
    progress_bar_->Finish();
    progress_bar_.reset();
    return;
  }
  const size_t value = *lowest >= start_section ? *lowest - start_section : 0;
  progress_bar_->Update(value);
}

void KineticDelaunay::SectionEventManager::finishProgress()
{
  if (progress_bar_)
  {
    progress_bar_->Finish();
    progress_bar_.reset();
  }
}

void KineticDelaunay::SectionEvent::handleEvent()
{
  auto* kd = kd_;
  if (!kd)
  {
    throw std::runtime_error("SectionEvent has no KineticDelaunay pointer");
  }
  // Retires, cuts and forks on the shared graph and runtime-branch tables.
  const auto topology_lock = kd->lockSharedTopology();
  if (!kd->eventInActiveRuntimeScope(*this))
  {
    kd->exportRuntimeScopeSkipSvg(*this);
    return;
  }

  const size_t section_index = section_id;
  const size_t section_count = kd->getSectionCount();

  if (const BranchScheduler* scheduler = kd->activeBranchScheduler(); scheduler != nullptr)
  {
    if (section_index >= scheduler->handoff_section
      || (scheduler->independent && section_index >= kd->exclusiveEndSectionForScheduler(*scheduler)))
    {
      KINDS_INFO("SectionEvent skip runtime_branch="
        << (scheduler->runtime_branch_id == BranchScheduler::no_branch
              ? std::string("root")
              : std::to_string(scheduler->runtime_branch_id))
        << " section=" << section_index << " handoff=" << scheduler->handoff_section
        << " exclusive_end=" << kd->exclusiveEndSectionForScheduler(*scheduler)
        << " clearing queue size=" << kd->kinetic_algorithm_->activeScheduler().size());
      kd->kinetic_algorithm_->activeScheduler().clear();
      return;
    }
    if (scheduler->independent && !scheduler->owned_strands.empty())
    {
      bool any_profile = false;
      for (size_t strand_id : scheduler->owned_strands)
      {
        if (kd->strandHasProfileAt(strand_id, static_cast<double>(section_index)))
        {
          any_profile = true;
          break;
        }
      }
      if (!any_profile)
      {
        KINDS_INFO("SectionEvent skip runtime_branch=" << scheduler->runtime_branch_id
                                                       << " section=" << section_index
                                                       << " (no owned strand profile), clearing queue size="
                                                       << kd->kinetic_algorithm_->activeScheduler().size());
        kd->kinetic_algorithm_->activeScheduler().clear();
        return;
      }
    }
  }

  assert(section_index < section_count);

  // This thread's scheduler must run its own section events in order. Other branches may be at
  // different heights; their cursors live on their schedulers, not in a global counter.
  assert(section_index == kd->getActiveSectionsAdvanced());

  auto* event_handler = kd->section_event_manager_->getCallback();

  // Bootstrap section (normally 0, or @ref getStartSection when starting mid-tree): init already set up meshes.
  const bool is_bootstrap_section = section_index == kd->getStartSection()
    && runtime_branch_id == BranchScheduler::no_branch;
  if (!is_bootstrap_section && event_handler)
  {
    event_handler->beforeEvent(*this);
  }

  // Retire finished input branches before graph cuts so vertex positions are not queried past strand data.
  kd->retireFinishedInputBranches(static_cast<double>(section_index), runtime_branch_id);
  // Lookups check their own height window, so expiry is only cleanup. Independent schedulers sit at different
  // heights; expiring here would drop transitions another branch still needs and rebuild the table under it.
  if (runtime_branch_id == BranchScheduler::no_branch)
  {
    kd->post_split_frame_transitions_.expireBeforeHeight(section_index);
  }

  if (kd->collectStatistics())
  {
    const auto [strand_count, branch_count] = kd->countLiveStrandsAndBranches(section_index, runtime_branch_id);
    std::optional<size_t> stats_branch_id;
    if (runtime_branch_id != BranchScheduler::no_branch)
    {
      stats_branch_id = runtime_branch_id;
    }
    kd->statistics().setSectionTopology(section_index, stats_branch_id, strand_count, branch_count);
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
        if (runtime_branch_id != BranchScheduler::no_branch
          && split.parent_runtime_branch != runtime_branch_id)
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

  if (BranchScheduler* scheduler = kd->activeBranchScheduler();
    scheduler != nullptr && scheduler->independent)
  {
    kd->continueIndependentSchedulerFromSection(*scheduler, section_index + 1);
  }

  if (!is_bootstrap_section)
  {
    kd->scheduleRadiusEventsForCutoffRegimeChanges(static_cast<double>(section_index));
  }

  kd->precomputeStepForRuntimeBranch(static_cast<double>(section_index), runtime_branch_id);

  kd->mutableActiveSectionsAdvanced()++;
  kd->section_event_manager_->updateProgress();

  kd->validateSitesInsideConvexHull("SectionEvent", occurrence_time);
  kd->validateDelaunayVoronoiIntersections("SectionEvent", occurrence_time, eventId(), creation_time);

  if (!is_bootstrap_section && event_handler)
  {
    event_handler->afterEvent(*this);
  }
}

#pragma once

#include "KineticAlgorithm.hpp"

#include <cstddef>
#include <cstdint>
#include <limits>
#include <unordered_set>
#include <vector>

namespace kinDS
{

/// Per-runtime-branch event schedule after a physical graph cut.
/// Serial / pre-fork mode uses a single scheduler with @ref runtime_branch_id == no_branch.
struct BranchScheduler
{
  static constexpr size_t no_branch = (std::numeric_limits<size_t>::max)();

  size_t runtime_branch_id = no_branch;
  EventScheduler queue {};
  size_t sections_advanced = 0;
  /// When true this scheduler is independent of the pre-fork global queue (post-cut).
  bool independent = false;
  /// Section from which forked schedulers own this queue's branches; its section event drains the queue.
  size_t handoff_section = no_branch;
  /// Owned strand ids (copied from runtime branch strand list at fork).
  std::vector<size_t> owned_strands {};
  /// Live even half-edge ids owned by this branch (filled at fork / refreshed after local topology ops).
  std::unordered_set<size_t> owned_even_half_edges {};
  /// Live face ids owned by this branch.
  std::unordered_set<size_t> owned_faces {};
};

} // namespace kinDS

#pragma once

#include <array>
#include <chrono>
#include <cstddef>
#include <filesystem>
#include <optional>
#include <string>
#include <vector>

namespace kinDS
{
/// Kinetic Delaunay event kinds counted by @ref Statistics.
enum class KineticEventType : size_t
{
  Subdivision = 0,
  Section,
  Separation,
  Flip,
  Radius,
  Crossing,
  Count
};

inline constexpr size_t kineticEventTypeCount = static_cast<size_t>(KineticEventType::Count);

const char* kineticEventTypeName(KineticEventType type);

/**
 * @brief Per-section and total wall-clock / event-count statistics for one tree-meshing run.
 *
 * Section attribution uses @c floor(occurrence_time). Wall time for section @c i is the span from the first
 * dequeued event belonging to @c i until the first event of a later section (or @ref endRun / finalize).
 * Strand/branch snapshots are optional and recorded explicitly after section retirement.
 */
class Statistics
{
 public:
  struct SectionStats
  {
    size_t section_id = 0;
    double runtime_seconds = 0.0;
    std::array<size_t, kineticEventTypeCount> event_counts {};
    /// Live non-dummy strands after retirement at this section; unset until snapshotted.
    std::optional<size_t> strand_count {};
    /// Alive runtime branches after retirement at this section; unset until snapshotted.
    std::optional<size_t> branch_count {};
  };

  void reset();

  /// Start a new collection window (call before processing kinetic events).
  void beginRun();

  /// Close the open section timer (call after the last timed work, including finalize if attributed here).
  void endRun();

  bool isCollecting() const { return run_active_; }

  /// Record a dequeued kinetic event (updates section timer + counts).
  void onEvent(KineticEventType type, double occurrence_time);

  /// Snapshot live strand/branch counts for @p section_id (after phased-out strands are retired).
  void setSectionTopology(size_t section_id, size_t strand_count, size_t branch_count);

  /// Add wall time to the current open section (and totals), e.g. @c SegmentBuilder::finalize.
  void addWallTimeSeconds(double seconds);

  bool empty() const { return sections_.empty(); }
  const std::vector<SectionStats>& sections() const { return sections_; }
  const SectionStats& totals() const { return totals_; }

  /// CSV: @c section_id,runtime_s,strand_count,branch_count,<event types...>; final @c total row.
  bool writeCsv(const std::filesystem::path& path) const;

 private:
  using Clock = std::chrono::steady_clock;

  void closeOpenSection(Clock::time_point now);
  void openSection(size_t section_id, Clock::time_point now);
  SectionStats& ensureSection(size_t section_id);

  bool run_active_ = false;
  bool section_open_ = false;
  size_t current_section_id_ = 0;
  Clock::time_point section_started_ {};
  std::vector<SectionStats> sections_;
  SectionStats totals_ {};
};
} // namespace kinDS

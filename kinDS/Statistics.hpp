#pragma once

#include <array>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <iosfwd>
#include <limits>
#include <memory>
#include <mutex>
#include <optional>
#include <string>
#include <unordered_map>
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
  Failsafe,
  Count
};

inline constexpr size_t kineticEventTypeCount = static_cast<size_t>(KineticEventType::Count);

const char* kineticEventTypeName(KineticEventType type);

/**
 * @brief Per-(section, runtime-branch) and total wall-clock / event-count statistics for one tree-meshing run.
 *
 * Section attribution uses @c floor(occurrence_time). Each runtime branch has its own open-section timer:
 * wall time for (@c i, branch) is the span from that branch's first dequeued event in @c i until its first
 * event of a later section (or @ref endRun). Strand/branch snapshots are optional: per-row values record
 * live topology of that runtime branch after retirement; the totals row uses the full input strand inventory
 * and all distinct input-tree branch IDs. All mutators and CSV writers are internally synchronized.
 */
class Statistics
{
 public:
  struct SectionStats
  {
    size_t section_id = 0;
    /// Blank on the pre-fork / serial scheduler (@c BranchScheduler::no_branch) and on the totals row.
    std::optional<size_t> runtime_branch_id {};
    double runtime_seconds = 0.0;
    std::array<size_t, kineticEventTypeCount> event_counts {};
    /// Per row: live non-dummy strands after retirement (scoped to the runtime branch when set).
    /// Totals: full input strand inventory.
    std::optional<size_t> strand_count {};
    /// Per unscoped row: alive runtime branches after retirement. Per scoped row: distinct input-tree
    /// branches among those live strands. Totals: all input-tree branches (alive + retired).
    std::optional<size_t> branch_count {};
  };

  /// One dequeued kinetic event for the companion event-list CSV.
  struct EventListRow
  {
    uint64_t event_id = 0;
    double occurrence_t = 0.0;
    double occurrence_infinitesimal_t = 0.0;
    KineticEventType type = KineticEventType::Count;
    std::optional<size_t> runtime_branch_id {};
    std::optional<size_t> half_edge_id {};
    std::optional<size_t> delaunay_edge_id {};
    std::optional<size_t> voronoi_vertex_id {};
    std::optional<size_t> section_id {};
    std::optional<size_t> strand_id {};
    std::optional<size_t> parent_component_id {};
    std::optional<double> split_time {};
    std::optional<bool> target_inside {};
    /// Set for radius events after meshing decides shift vs traced-cell fallback; blank for other types.
    std::optional<bool> radius_shift {};
  };

  Statistics();
  ~Statistics();
  Statistics(const Statistics&) = delete;
  Statistics& operator=(const Statistics&) = delete;

  void reset();

  /// Start a new collection window (call before processing kinetic events).
  void beginRun();

  /// Close the open section timer (call after the last timed work, including finalize if attributed here).
  void endRun();

  bool isCollecting() const;

  /// Record a dequeued kinetic event (updates that branch's section timer + counts). Prefer @ref recordEvent
  /// when the full event object is available so the event-list CSV gets typed columns.
  void onEvent(KineticEventType type, double occurrence_time, std::optional<size_t> runtime_branch_id = {});

  /// Record counts and append an event-list row (type-specific fields filled via @p row).
  void recordEvent(EventListRow row);

  /// Fill @c radius_shift on the event-list row matching @p event_id (no-op if not collecting / not found).
  void setRadiusShift(uint64_t event_id, bool shifted);

  /// Snapshot live strand/branch counts for (@p section_id, @p runtime_branch_id) after retirement.
  void setSectionTopology(size_t section_id, std::optional<size_t> runtime_branch_id, size_t strand_count,
    size_t branch_count);

  /// Totals-row topology: full strand inventory and all input-tree branches (alive + retired).
  void setTotalsTopology(size_t strand_count, size_t branch_count);

  /// Totals-row extras written only on the @c total CSV row (blank on per-section rows).
  void setTotalsAlpha(double alpha);
  void setTotalsMeshCounts(size_t triangle_count, size_t vertex_count);
  /// Optional failure note for the totals row (e.g. meshing exception during a parameter sweep).
  void setTotalsFailure(std::string message);
  /// Append a recoverable-failure note; keeps prior messages so failsafe skips are not silent.
  void appendTotalsFailure(const std::string& message);

  /// Optional experiment tag inserted into the timestamped CSV stem
  /// (@c meshing_statistics_<tag>_<sectionCount>_YYYYMMDD_…). Spaces should already be underscores.
  void setFilenameExperimentTag(std::string tag);

  /// After @ref beginRun, open a stable incremental CSV (header only) and keep it open for the run.
  /// On Windows the handle uses @c FILE_SHARE_READ so other processes (LibreOffice) can only open
  /// read-only; our threads append through @ref mutex_. One row is written each time a section closes.
  void startIncrementalCsv(const std::filesystem::path& base_path);

  /// Add wall time to the totals row only (not branch-local), e.g. @c SegmentBuilder::finalize.
  void addWallTimeSeconds(double seconds);

  bool empty() const;
  const std::vector<SectionStats>& sections() const { return sections_; }
  const SectionStats& totals() const { return totals_; }
  const std::vector<EventListRow>& eventList() const { return event_list_; }
  const std::optional<double>& totalsAlpha() const { return totals_alpha_; }
  const std::optional<size_t>& totalsTriangleCount() const { return totals_triangle_count_; }
  const std::optional<size_t>& totalsVertexCount() const { return totals_vertex_count_; }
  const std::optional<std::string>& totalsFailure() const { return totals_failure_; }

  /// Insert a local timestamp before the extension so repeated writes never collide.
  /// @c meshing_statistics.csv becomes @c meshing_statistics_YYYYMMDD_HHMMSS_mmm.csv.
  /// When @p experiment_tag / section count are set on this instance, @ref writeCsv builds a richer stem.
  static std::filesystem::path timestampedCsvPath(const std::filesystem::path& path);

  /// Companion event-list base path beside a statistics CSV
  /// (@c meshing_statistics.csv → @c meshing_event_list.csv).
  static std::filesystem::path eventListCsvPathBeside(const std::filesystem::path& statistics_csv_path);

  /// CSV: @c section_id,runtime_branch_id,runtime_s,strand_count,branch_count,segment_count,<event types...>,
  /// @c alpha,triangle_count,vertex_count,failure;
  /// @c segment_count is @c strand_count + subdivision events when strand_count is set (else blank).
  /// Rows are sorted by @c section_id then @c runtime_branch_id (blank branch id first).
  /// Final @c total row (alpha / mesh counts / failure filled only there when set).
  /// Writes to a timestamped sibling of @p path so existing files are never overwritten.
  bool writeCsv(const std::filesystem::path& path) const;

  /// CSV: one row per kinetic event in dequeue order (blank cells when not applicable).
  /// Writes to a timestamped sibling of @p path.
  bool writeEventListCsv(const std::filesystem::path& path) const;

 private:
  using Clock = std::chrono::steady_clock;

  struct OpenSectionTimer
  {
    size_t section_id = 0;
    std::optional<size_t> runtime_branch_id {};
    Clock::time_point started {};
  };

  static size_t branchTimerKey(std::optional<size_t> runtime_branch_id);
  static bool sectionRowLess(const SectionStats& a, const SectionStats& b);

  void resetLocked();
  void onEventLocked(KineticEventType type, double occurrence_time, std::optional<size_t> runtime_branch_id);
  void closeOpenBranch(size_t timer_key, Clock::time_point now);
  void openBranchSection(size_t section_id, std::optional<size_t> runtime_branch_id, Clock::time_point now);
  SectionStats& ensureSection(size_t section_id, std::optional<size_t> runtime_branch_id);
  std::filesystem::path statisticsCsvStemPath(const std::filesystem::path& path) const;
  void writeCsvHeader(std::ostream& out) const;
  void writeCsvSectionRow(std::ostream& out, const SectionStats& row) const;
  void appendIncrementalSectionRow(const SectionStats& row);
  void closeIncrementalCsvLocked();
  bool writeIncrementalCsvLocked(const std::string& text);

  mutable std::mutex mutex_;
  bool run_active_ = false;
  std::unordered_map<size_t, OpenSectionTimer> open_by_branch_;
  std::vector<SectionStats> sections_;
  SectionStats totals_ {};
  std::vector<EventListRow> event_list_;
  std::optional<double> totals_alpha_ {};
  std::optional<size_t> totals_triangle_count_ {};
  std::optional<size_t> totals_vertex_count_ {};
  std::optional<std::string> totals_failure_ {};
  std::string filename_experiment_tag_ {};
  std::filesystem::path incremental_csv_path_ {};
  struct IncrementalCsvFile;
  std::unique_ptr<IncrementalCsvFile> incremental_csv_file_;
};
} // namespace kinDS

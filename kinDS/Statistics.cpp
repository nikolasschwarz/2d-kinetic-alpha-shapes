#include "Statistics.hpp"

#include "Logger.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstring>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <limits>
#include <mutex>
#include <ostream>
#include <sstream>
#include <string>
#include <vector>

#ifdef _WIN32
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#else
#include <errno.h>
#include <fcntl.h>
#include <sys/file.h>
#include <unistd.h>
#endif

namespace kinDS
{
struct Statistics::IncrementalCsvFile
{
#ifdef _WIN32
  HANDLE handle = INVALID_HANDLE_VALUE;
#else
  int fd = -1;
#endif

  ~IncrementalCsvFile() { close(); }

  IncrementalCsvFile() = default;
  IncrementalCsvFile(const IncrementalCsvFile&) = delete;
  IncrementalCsvFile& operator=(const IncrementalCsvFile&) = delete;
  IncrementalCsvFile(IncrementalCsvFile&&) = delete;
  IncrementalCsvFile& operator=(IncrementalCsvFile&&) = delete;

  void close()
  {
#ifdef _WIN32
    if (handle != INVALID_HANDLE_VALUE)
    {
      CloseHandle(handle);
      handle = INVALID_HANDLE_VALUE;
    }
#else
    if (fd >= 0)
    {
      flock(fd, LOCK_UN);
      ::close(fd);
      fd = -1;
    }
#endif
  }

  bool isOpen() const
  {
#ifdef _WIN32
    return handle != INVALID_HANDLE_VALUE;
#else
    return fd >= 0;
#endif
  }

  bool open(const std::filesystem::path& path)
  {
    close();
#ifdef _WIN32
    handle = CreateFileW(path.c_str(), GENERIC_WRITE, FILE_SHARE_READ, nullptr, CREATE_ALWAYS,
      FILE_ATTRIBUTE_NORMAL, nullptr);
    return handle != INVALID_HANDLE_VALUE;
#else
    fd = ::open(path.c_str(), O_WRONLY | O_CREAT | O_TRUNC | O_CLOEXEC, 0644);
    if (fd < 0)
    {
      return false;
    }
    if (flock(fd, LOCK_EX) != 0)
    {
      ::close(fd);
      fd = -1;
      return false;
    }
    return true;
#endif
  }

  bool append(const std::string& text)
  {
    if (text.empty())
    {
      return true;
    }
#ifdef _WIN32
    if (handle == INVALID_HANDLE_VALUE)
    {
      return false;
    }
    const char* cursor = text.data();
    DWORD remaining = static_cast<DWORD>(text.size());
    while (remaining > 0)
    {
      DWORD written = 0;
      if (!WriteFile(handle, cursor, remaining, &written, nullptr))
      {
        return false;
      }
      cursor += written;
      remaining -= written;
    }
    return FlushFileBuffers(handle) != 0;
#else
    if (fd < 0)
    {
      return false;
    }
    const char* cursor = text.data();
    size_t remaining = text.size();
    while (remaining > 0)
    {
      const ssize_t written = ::write(fd, cursor, remaining);
      if (written < 0)
      {
        if (errno == EINTR)
        {
          continue;
        }
        return false;
      }
      cursor += static_cast<size_t>(written);
      remaining -= static_cast<size_t>(written);
    }
    return ::fsync(fd) == 0;
#endif
  }
};

Statistics::Statistics() = default;
Statistics::~Statistics() = default;

const char* kineticEventTypeName(KineticEventType type)
{
  switch (type)
  {
  case KineticEventType::Subdivision:
    return "subdivision";
  case KineticEventType::Section:
    return "section";
  case KineticEventType::Separation:
    return "separation";
  case KineticEventType::Flip:
    return "flip";
  case KineticEventType::Radius:
    return "radius";
  case KineticEventType::Crossing:
    return "crossing";
  case KineticEventType::Failsafe:
    return "failsafe";
  case KineticEventType::Count:
    break;
  }
  return "unknown";
}

size_t Statistics::branchTimerKey(std::optional<size_t> runtime_branch_id)
{
  return runtime_branch_id.value_or((std::numeric_limits<size_t>::max)());
}

bool Statistics::sectionRowLess(const SectionStats& a, const SectionStats& b)
{
  if (a.section_id != b.section_id)
  {
    return a.section_id < b.section_id;
  }
  if (a.runtime_branch_id.has_value() != b.runtime_branch_id.has_value())
  {
    return !a.runtime_branch_id.has_value();
  }
  return a.runtime_branch_id.has_value() && *a.runtime_branch_id < *b.runtime_branch_id;
}

void Statistics::resetLocked()
{
  run_active_ = false;
  open_by_branch_.clear();
  sections_.clear();
  totals_ = {};
  event_list_.clear();
  totals_alpha_.reset();
  totals_triangle_count_.reset();
  totals_vertex_count_.reset();
  totals_failure_.reset();
  filename_experiment_tag_.clear();
  closeIncrementalCsvLocked();
}

void Statistics::reset()
{
  std::lock_guard<std::mutex> lock(mutex_);
  resetLocked();
}

void Statistics::beginRun()
{
  std::lock_guard<std::mutex> lock(mutex_);
  resetLocked();
  run_active_ = true;
}

void Statistics::endRun()
{
  std::lock_guard<std::mutex> lock(mutex_);
  if (!run_active_)
  {
    return;
  }
  const Clock::time_point now = Clock::now();
  std::vector<size_t> keys;
  keys.reserve(open_by_branch_.size());
  for (const auto& entry : open_by_branch_)
  {
    keys.push_back(entry.first);
  }
  for (const size_t key : keys)
  {
    closeOpenBranch(key, now);
  }
  std::sort(sections_.begin(), sections_.end(), sectionRowLess);
  run_active_ = false;
  closeIncrementalCsvLocked();
}

bool Statistics::isCollecting() const
{
  std::lock_guard<std::mutex> lock(mutex_);
  return run_active_;
}

bool Statistics::empty() const
{
  std::lock_guard<std::mutex> lock(mutex_);
  return sections_.empty();
}

void Statistics::closeOpenBranch(size_t timer_key, Clock::time_point now)
{
  const auto it = open_by_branch_.find(timer_key);
  if (it == open_by_branch_.end())
  {
    return;
  }
  const std::chrono::duration<double> elapsed = now - it->second.started;
  const double seconds = elapsed.count();
  SectionStats& row = ensureSection(it->second.section_id, it->second.runtime_branch_id);
  row.runtime_seconds += seconds;
  totals_.runtime_seconds += seconds;
  appendIncrementalSectionRow(row);
  open_by_branch_.erase(it);
}

void Statistics::startIncrementalCsv(const std::filesystem::path& base_path)
{
  std::lock_guard<std::mutex> lock(mutex_);
  if (!run_active_ || base_path.empty())
  {
    return;
  }
  std::filesystem::path base = base_path;
  if (base.has_filename() && base.filename() == ".")
  {
    base /= "meshing_statistics.csv";
  }
  else if (!base.has_filename())
  {
    base /= "meshing_statistics.csv";
  }
  std::string stem = base.stem().string();
  if (stem.empty())
  {
    stem = "meshing_statistics";
  }
  if (!filename_experiment_tag_.empty())
  {
    stem += "_";
    stem += filename_experiment_tag_;
  }
  stem += "_partial";
  std::string extension = base.extension().string();
  if (extension.empty())
  {
    extension = ".csv";
  }
  incremental_csv_path_ = timestampedCsvPath(base.parent_path() / (stem + extension));
  incremental_csv_file_ = std::make_unique<IncrementalCsvFile>();
  if (!incremental_csv_file_->open(incremental_csv_path_))
  {
    KINDS_WARNING("Statistics: failed to open incremental CSV " << incremental_csv_path_.generic_string());
    closeIncrementalCsvLocked();
    return;
  }
  std::ostringstream header;
  writeCsvHeader(header);
  if (!writeIncrementalCsvLocked(header.str()))
  {
    KINDS_WARNING("Statistics: failed to write incremental CSV header " << incremental_csv_path_.generic_string());
    closeIncrementalCsvLocked();
    return;
  }
  KINDS_INFO("Statistics: incremental CSV " << incremental_csv_path_.generic_string());
}

void Statistics::writeCsvHeader(std::ostream& out) const
{
  out << "section_id,runtime_branch_id,runtime_s,strand_count,branch_count,segment_count";
  for (size_t i = 0; i < kineticEventTypeCount; ++i)
  {
    out << ',' << kineticEventTypeName(static_cast<KineticEventType>(i));
  }
  out << ",alpha,triangle_count,vertex_count,failure\n";
}

void Statistics::writeCsvSectionRow(std::ostream& out, const SectionStats& row) const
{
  out << std::setprecision(std::numeric_limits<double>::max_digits10);
  out << row.section_id << ',';
  if (row.runtime_branch_id.has_value())
  {
    out << row.runtime_branch_id.value();
  }
  out << ',' << row.runtime_seconds << ',';
  if (row.strand_count.has_value())
  {
    out << row.strand_count.value();
  }
  out << ',';
  if (row.branch_count.has_value())
  {
    out << row.branch_count.value();
  }
  out << ',';
  if (row.strand_count.has_value())
  {
    out << (row.strand_count.value() + row.event_counts[static_cast<size_t>(KineticEventType::Subdivision)]);
  }
  for (size_t i = 0; i < kineticEventTypeCount; ++i)
  {
    out << ',' << row.event_counts[i];
  }
  // Section rows leave totals-only columns blank.
  out << ",,,\n";
}

void Statistics::closeIncrementalCsvLocked()
{
  incremental_csv_file_.reset();
  incremental_csv_path_.clear();
}

bool Statistics::writeIncrementalCsvLocked(const std::string& text)
{
  return incremental_csv_file_ && incremental_csv_file_->append(text);
}

void Statistics::appendIncrementalSectionRow(const SectionStats& row)
{
  if (!incremental_csv_file_ || !incremental_csv_file_->isOpen())
  {
    return;
  }
  std::ostringstream out;
  writeCsvSectionRow(out, row);
  if (!writeIncrementalCsvLocked(out.str()))
  {
    KINDS_WARNING("Statistics: failed to append incremental CSV " << incremental_csv_path_.generic_string());
  }
}

void Statistics::openBranchSection(size_t section_id, std::optional<size_t> runtime_branch_id, Clock::time_point now)
{
  ensureSection(section_id, runtime_branch_id);
  OpenSectionTimer timer;
  timer.section_id = section_id;
  timer.runtime_branch_id = runtime_branch_id;
  timer.started = now;
  open_by_branch_[branchTimerKey(runtime_branch_id)] = timer;
}

Statistics::SectionStats& Statistics::ensureSection(size_t section_id, std::optional<size_t> runtime_branch_id)
{
  for (SectionStats& row : sections_)
  {
    if (row.section_id == section_id && row.runtime_branch_id == runtime_branch_id)
    {
      return row;
    }
  }
  SectionStats row;
  row.section_id = section_id;
  row.runtime_branch_id = runtime_branch_id;
  sections_.push_back(row);
  return sections_.back();
}

void Statistics::onEventLocked(KineticEventType type, double occurrence_time, std::optional<size_t> runtime_branch_id)
{
  if (!run_active_ || type >= KineticEventType::Count || !std::isfinite(occurrence_time) || occurrence_time < 0.0)
  {
    return;
  }

  const size_t section_id = static_cast<size_t>(std::floor(occurrence_time));
  const Clock::time_point now = Clock::now();
  const size_t timer_key = branchTimerKey(runtime_branch_id);
  const auto it = open_by_branch_.find(timer_key);
  if (it == open_by_branch_.end() || it->second.section_id != section_id)
  {
    closeOpenBranch(timer_key, now);
    openBranchSection(section_id, runtime_branch_id, now);
  }

  const size_t type_index = static_cast<size_t>(type);
  ensureSection(section_id, runtime_branch_id).event_counts[type_index] += 1;
  totals_.event_counts[type_index] += 1;
}

void Statistics::onEvent(KineticEventType type, double occurrence_time, std::optional<size_t> runtime_branch_id)
{
  std::lock_guard<std::mutex> lock(mutex_);
  onEventLocked(type, occurrence_time, runtime_branch_id);
}

void Statistics::recordEvent(EventListRow row)
{
  std::lock_guard<std::mutex> lock(mutex_);
  if (!run_active_ || row.type >= KineticEventType::Count || !std::isfinite(row.occurrence_t)
    || row.occurrence_t < 0.0)
  {
    return;
  }
  onEventLocked(row.type, row.occurrence_t, row.runtime_branch_id);
  event_list_.push_back(std::move(row));
}

void Statistics::setRadiusShift(uint64_t event_id, bool shifted)
{
  std::lock_guard<std::mutex> lock(mutex_);
  if (!run_active_)
  {
    return;
  }
  for (auto it = event_list_.rbegin(); it != event_list_.rend(); ++it)
  {
    if (it->event_id == event_id)
    {
      it->radius_shift = shifted;
      return;
    }
  }
}

void Statistics::setSectionTopology(size_t section_id, std::optional<size_t> runtime_branch_id, size_t strand_count,
  size_t branch_count)
{
  std::lock_guard<std::mutex> lock(mutex_);
  if (!run_active_)
  {
    return;
  }
  SectionStats& row = ensureSection(section_id, runtime_branch_id);
  row.strand_count = strand_count;
  row.branch_count = branch_count;
}

void Statistics::setTotalsTopology(size_t strand_count, size_t branch_count)
{
  std::lock_guard<std::mutex> lock(mutex_);
  totals_.strand_count = strand_count;
  totals_.branch_count = branch_count;
}

void Statistics::setTotalsAlpha(double alpha)
{
  std::lock_guard<std::mutex> lock(mutex_);
  totals_alpha_ = alpha;
}

void Statistics::setTotalsMeshCounts(size_t triangle_count, size_t vertex_count)
{
  std::lock_guard<std::mutex> lock(mutex_);
  totals_triangle_count_ = triangle_count;
  totals_vertex_count_ = vertex_count;
}

void Statistics::setTotalsFailure(std::string message)
{
  std::lock_guard<std::mutex> lock(mutex_);
  totals_failure_ = std::move(message);
}

void Statistics::appendTotalsFailure(const std::string& message)
{
  std::lock_guard<std::mutex> lock(mutex_);
  if (message.empty())
  {
    return;
  }
  if (!totals_failure_.has_value() || totals_failure_->empty())
  {
    totals_failure_ = message;
    return;
  }
  *totals_failure_ += " | ";
  *totals_failure_ += message;
}

void Statistics::setFilenameExperimentTag(std::string tag)
{
  std::lock_guard<std::mutex> lock(mutex_);
  filename_experiment_tag_ = std::move(tag);
}

void Statistics::addWallTimeSeconds(double seconds)
{
  std::lock_guard<std::mutex> lock(mutex_);
  if (!run_active_ || !(seconds > 0.0) || !std::isfinite(seconds))
  {
    return;
  }
  totals_.runtime_seconds += seconds;
}

std::filesystem::path Statistics::timestampedCsvPath(const std::filesystem::path& path)
{
  std::filesystem::path base = path.empty() ? std::filesystem::path("meshing_statistics.csv") : path;
  if (base.has_filename() && base.filename() == ".")
  {
    base /= "meshing_statistics.csv";
  }
  else if (!base.has_filename())
  {
    base /= "meshing_statistics.csv";
  }

  const auto now = std::chrono::system_clock::now();
  const std::time_t time = std::chrono::system_clock::to_time_t(now);
  const auto milliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(now.time_since_epoch()) % 1000;
  std::tm local {};
#ifdef _WIN32
  localtime_s(&local, &time);
#else
  localtime_r(&time, &local);
#endif
  std::ostringstream stamp;
  stamp << std::put_time(&local, "%Y%m%d_%H%M%S") << '_' << std::setw(3) << std::setfill('0') << milliseconds.count();

  std::string stem = base.stem().string();
  if (stem.empty())
  {
    stem = "meshing_statistics";
  }
  std::string extension = base.extension().string();
  if (extension.empty())
  {
    extension = ".csv";
  }
  return base.parent_path() / (stem + "_" + stamp.str() + extension);
}

std::filesystem::path Statistics::eventListCsvPathBeside(const std::filesystem::path& statistics_csv_path)
{
  std::filesystem::path base
    = statistics_csv_path.empty() ? std::filesystem::path("meshing_statistics.csv") : statistics_csv_path;
  if (base.has_filename() && base.filename() == ".")
  {
    base /= "meshing_statistics.csv";
  }
  else if (!base.has_filename())
  {
    base /= "meshing_statistics.csv";
  }

  std::string stem = base.stem().string();
  if (stem.empty())
  {
    stem = "meshing_event_list";
  }
  else
  {
    const size_t pos = stem.find("statistics");
    if (pos != std::string::npos)
    {
      stem.replace(pos, std::strlen("statistics"), "event_list");
    }
    else
    {
      stem += "_event_list";
    }
  }
  std::string extension = base.extension().string();
  if (extension.empty())
  {
    extension = ".csv";
  }
  return base.parent_path() / (stem + extension);
}

std::filesystem::path Statistics::statisticsCsvStemPath(const std::filesystem::path& path) const
{
  std::filesystem::path base = path.empty() ? std::filesystem::path("meshing_statistics.csv") : path;
  if (base.has_filename() && base.filename() == ".")
  {
    base /= "meshing_statistics.csv";
  }
  else if (!base.has_filename())
  {
    base /= "meshing_statistics.csv";
  }

  std::string stem = base.stem().string();
  if (stem.empty())
  {
    stem = "meshing_statistics";
  }
  if (!filename_experiment_tag_.empty())
  {
    stem += "_";
    stem += filename_experiment_tag_;
  }
  stem += "_";
  stem += std::to_string(sections_.size());

  std::string extension = base.extension().string();
  if (extension.empty())
  {
    extension = ".csv";
  }
  return base.parent_path() / (stem + extension);
}

bool Statistics::writeCsv(const std::filesystem::path& path) const
{
  std::lock_guard<std::mutex> lock(mutex_);
  const std::filesystem::path unique_path = timestampedCsvPath(statisticsCsvStemPath(path));
  std::ofstream out(unique_path);
  if (!out)
  {
    KINDS_WARNING("Statistics: failed to open CSV " << unique_path.generic_string());
    return false;
  }

  writeCsvHeader(out);

  out << std::setprecision(std::numeric_limits<double>::max_digits10);
  auto write_optional_size = [&](const std::optional<size_t>& value)
  {
    if (value.has_value())
    {
      out << value.value();
    }
  };
  auto write_optional_double = [&](const std::optional<double>& value)
  {
    if (value.has_value())
    {
      out << value.value();
    }
  };
  auto write_csv_string = [&](const std::string& value)
  {
    out << '"';
    for (const char c : value)
    {
      if (c == '"')
      {
        out << "\"\"";
      }
      else
      {
        out << c;
      }
    }
    out << '"';
  };
  auto write_row = [&](const std::string& id, const SectionStats& row, const bool is_total)
  {
    out << id << ',';
    if (!is_total)
    {
      write_optional_size(row.runtime_branch_id);
    }
    out << ',' << row.runtime_seconds << ',';
    write_optional_size(row.strand_count);
    out << ',';
    write_optional_size(row.branch_count);
    out << ',';
    // Rod / segment count ≈ initial strands + subdivision events (one new segment per subdiv).
    if (row.strand_count.has_value())
    {
      out << (row.strand_count.value() + row.event_counts[static_cast<size_t>(KineticEventType::Subdivision)]);
    }
    for (size_t i = 0; i < kineticEventTypeCount; ++i)
    {
      out << ',' << row.event_counts[i];
    }
    out << ',';
    if (is_total)
    {
      write_optional_double(totals_alpha_);
    }
    out << ',';
    if (is_total)
    {
      write_optional_size(totals_triangle_count_);
    }
    out << ',';
    if (is_total)
    {
      write_optional_size(totals_vertex_count_);
    }
    out << ',';
    if (is_total && totals_failure_.has_value())
    {
      write_csv_string(totals_failure_.value());
    }
    out << '\n';
  };

  std::vector<SectionStats> ordered = sections_;
  std::sort(ordered.begin(), ordered.end(), sectionRowLess);
  for (const SectionStats& row : ordered)
  {
    write_row(std::to_string(row.section_id), row, /*is_total=*/false);
  }
  write_row("total", totals_, /*is_total=*/true);

  KINDS_INFO("Statistics: wrote meshing CSV to " << unique_path.generic_string() << " (" << sections_.size()
                                                 << " section/branch row(s) + total)");
  return true;
}

bool Statistics::writeEventListCsv(const std::filesystem::path& path) const
{
  std::lock_guard<std::mutex> lock(mutex_);
  // Match statistics naming (experiment tag + section count) then swap statistics → event_list.
  const std::filesystem::path unique_path
    = timestampedCsvPath(eventListCsvPathBeside(statisticsCsvStemPath(path)));
  std::ofstream out(unique_path);
  if (!out)
  {
    KINDS_WARNING("Statistics: failed to open event-list CSV " << unique_path.generic_string());
    return false;
  }

  out << "event_id,occurrence_t,occurrence_infinitesimal_t,event_type,runtime_branch_id,half_edge_id,"
         "delaunay_edge_id,voronoi_vertex_id,section_id,strand_id,parent_component_id,split_time,"
         "target_inside,radius_shift\n";
  out << std::setprecision(std::numeric_limits<double>::max_digits10);

  auto write_optional_size = [&](const std::optional<size_t>& value)
  {
    if (value.has_value())
    {
      out << value.value();
    }
  };
  auto write_optional_double = [&](const std::optional<double>& value)
  {
    if (value.has_value())
    {
      out << value.value();
    }
  };
  auto write_optional_bool01 = [&](const std::optional<bool>& value)
  {
    if (value.has_value())
    {
      out << (value.value() ? 1 : 0);
    }
  };

  for (const EventListRow& row : event_list_)
  {
    out << row.event_id << ',' << row.occurrence_t << ',' << row.occurrence_infinitesimal_t << ','
        << kineticEventTypeName(row.type) << ',';
    write_optional_size(row.runtime_branch_id);
    out << ',';
    write_optional_size(row.half_edge_id);
    out << ',';
    write_optional_size(row.delaunay_edge_id);
    out << ',';
    write_optional_size(row.voronoi_vertex_id);
    out << ',';
    write_optional_size(row.section_id);
    out << ',';
    write_optional_size(row.strand_id);
    out << ',';
    write_optional_size(row.parent_component_id);
    out << ',';
    write_optional_double(row.split_time);
    out << ',';
    write_optional_bool01(row.target_inside);
    out << ',';
    write_optional_bool01(row.radius_shift);
    out << '\n';
  }

  KINDS_INFO("Statistics: wrote event-list CSV to " << unique_path.generic_string() << " (" << event_list_.size()
                                                    << " event row(s))");
  return true;
}
} // namespace kinDS

#include "Statistics.hpp"

#include "Logger.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <limits>
#include <sstream>

namespace kinDS
{
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
  case KineticEventType::Count:
    break;
  }
  return "unknown";
}

void Statistics::reset()
{
  run_active_ = false;
  section_open_ = false;
  current_section_id_ = 0;
  sections_.clear();
  totals_ = {};
}

void Statistics::beginRun()
{
  reset();
  run_active_ = true;
}

void Statistics::endRun()
{
  if (!run_active_)
  {
    return;
  }
  closeOpenSection(Clock::now());
  run_active_ = false;
}

void Statistics::closeOpenSection(Clock::time_point now)
{
  if (!section_open_)
  {
    return;
  }
  const std::chrono::duration<double> elapsed = now - section_started_;
  const double seconds = elapsed.count();
  SectionStats& row = ensureSection(current_section_id_);
  row.runtime_seconds += seconds;
  totals_.runtime_seconds += seconds;
  section_open_ = false;
}

void Statistics::openSection(size_t section_id, Clock::time_point now)
{
  ensureSection(section_id);
  current_section_id_ = section_id;
  section_started_ = now;
  section_open_ = true;
}

Statistics::SectionStats& Statistics::ensureSection(size_t section_id)
{
  for (SectionStats& row : sections_)
  {
    if (row.section_id == section_id)
    {
      return row;
    }
  }
  SectionStats row;
  row.section_id = section_id;
  sections_.push_back(row);
  return sections_.back();
}

void Statistics::onEvent(KineticEventType type, double occurrence_time)
{
  if (!run_active_ || type >= KineticEventType::Count || !std::isfinite(occurrence_time) || occurrence_time < 0.0)
  {
    return;
  }

  const size_t section_id = static_cast<size_t>(std::floor(occurrence_time));
  const Clock::time_point now = Clock::now();
  if (!section_open_ || section_id != current_section_id_)
  {
    closeOpenSection(now);
    openSection(section_id, now);
  }

  const size_t type_index = static_cast<size_t>(type);
  ensureSection(section_id).event_counts[type_index] += 1;
  totals_.event_counts[type_index] += 1;
}

void Statistics::setSectionTopology(size_t section_id, size_t strand_count, size_t branch_count)
{
  if (!run_active_)
  {
    return;
  }
  SectionStats& row = ensureSection(section_id);
  row.strand_count = strand_count;
  row.branch_count = branch_count;
}

void Statistics::setTotalsTopology(size_t strand_count, size_t branch_count)
{
  totals_.strand_count = strand_count;
  totals_.branch_count = branch_count;
}

void Statistics::addWallTimeSeconds(double seconds)
{
  if (!run_active_ || !(seconds > 0.0) || !std::isfinite(seconds))
  {
    return;
  }
  if (!section_open_)
  {
    if (!sections_.empty())
    {
      openSection(sections_.back().section_id, Clock::now());
    }
    else
    {
      totals_.runtime_seconds += seconds;
      return;
    }
  }
  ensureSection(current_section_id_).runtime_seconds += seconds;
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

bool Statistics::writeCsv(const std::filesystem::path& path) const
{
  const std::filesystem::path unique_path = timestampedCsvPath(path);
  std::ofstream out(unique_path);
  if (!out)
  {
    KINDS_WARNING("Statistics: failed to open CSV " << unique_path.generic_string());
    return false;
  }

  out << "section_id,runtime_s,strand_count,branch_count";
  for (size_t i = 0; i < kineticEventTypeCount; ++i)
  {
    out << ',' << kineticEventTypeName(static_cast<KineticEventType>(i));
  }
  out << '\n';

  out << std::setprecision(std::numeric_limits<double>::max_digits10);
  auto write_optional_size = [&](const std::optional<size_t>& value)
  {
    if (value.has_value())
    {
      out << value.value();
    }
  };
  auto write_row = [&](const std::string& id, const SectionStats& row)
  {
    out << id << ',' << row.runtime_seconds << ',';
    write_optional_size(row.strand_count);
    out << ',';
    write_optional_size(row.branch_count);
    for (size_t i = 0; i < kineticEventTypeCount; ++i)
    {
      out << ',' << row.event_counts[i];
    }
    out << '\n';
  };

  std::vector<SectionStats> ordered = sections_;
  std::sort(ordered.begin(), ordered.end(),
    [](const SectionStats& a, const SectionStats& b) { return a.section_id < b.section_id; });
  for (const SectionStats& row : ordered)
  {
    write_row(std::to_string(row.section_id), row);
  }
  write_row("total", totals_);

  KINDS_INFO("Statistics: wrote meshing CSV to " << unique_path.generic_string() << " (" << sections_.size()
                                                 << " section row(s) + total)");
  return true;
}
} // namespace kinDS

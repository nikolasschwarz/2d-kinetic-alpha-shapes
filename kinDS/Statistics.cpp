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
#include <ostream>
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
  event_list_.clear();
  totals_alpha_.reset();
  totals_triangle_count_.reset();
  totals_vertex_count_.reset();
  totals_failure_.reset();
  filename_experiment_tag_.clear();
  incremental_csv_path_.clear();
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
  appendIncrementalSectionRow(row);
}

void Statistics::startIncrementalCsv(const std::filesystem::path& base_path)
{
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

  std::ofstream out(incremental_csv_path_);
  if (!out)
  {
    KINDS_WARNING("Statistics: failed to open incremental CSV " << incremental_csv_path_.generic_string());
    incremental_csv_path_.clear();
    return;
  }
  writeCsvHeader(out);
  out.flush();
  KINDS_INFO("Statistics: incremental CSV " << incremental_csv_path_.generic_string());
}

void Statistics::writeCsvHeader(std::ostream& out) const
{
  out << "section_id,runtime_s,strand_count,branch_count,segment_count";
  for (size_t i = 0; i < kineticEventTypeCount; ++i)
  {
    out << ',' << kineticEventTypeName(static_cast<KineticEventType>(i));
  }
  out << ",alpha,triangle_count,vertex_count,failure\n";
}

void Statistics::writeCsvSectionRow(std::ostream& out, const SectionStats& row) const
{
  out << std::setprecision(std::numeric_limits<double>::max_digits10);
  out << row.section_id << ',' << row.runtime_seconds << ',';
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

void Statistics::appendIncrementalSectionRow(const SectionStats& row)
{
  if (incremental_csv_path_.empty())
  {
    return;
  }
  std::ofstream out(incremental_csv_path_, std::ios::app);
  if (!out)
  {
    KINDS_WARNING("Statistics: failed to append incremental CSV " << incremental_csv_path_.generic_string());
    return;
  }
  writeCsvSectionRow(out, row);
  out.flush();
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

void Statistics::recordEvent(EventListRow row)
{
  if (!run_active_ || row.type >= KineticEventType::Count || !std::isfinite(row.occurrence_t)
    || row.occurrence_t < 0.0)
  {
    return;
  }
  onEvent(row.type, row.occurrence_t);
  event_list_.push_back(std::move(row));
}

void Statistics::setRadiusShift(uint64_t event_id, bool shifted)
{
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

void Statistics::setTotalsAlpha(double alpha)
{
  totals_alpha_ = alpha;
}

void Statistics::setTotalsMeshCounts(size_t triangle_count, size_t vertex_count)
{
  totals_triangle_count_ = triangle_count;
  totals_vertex_count_ = vertex_count;
}

void Statistics::setTotalsFailure(std::string message)
{
  totals_failure_ = std::move(message);
}

void Statistics::setFilenameExperimentTag(std::string tag)
{
  filename_experiment_tag_ = std::move(tag);
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
  const std::filesystem::path unique_path = timestampedCsvPath(statisticsCsvStemPath(path));
  std::ofstream out(unique_path);
  if (!out)
  {
    KINDS_WARNING("Statistics: failed to open CSV " << unique_path.generic_string());
    return false;
  }

  out << "section_id,runtime_s,strand_count,branch_count,segment_count";
  for (size_t i = 0; i < kineticEventTypeCount; ++i)
  {
    out << ',' << kineticEventTypeName(static_cast<KineticEventType>(i));
  }
  out << ",alpha,triangle_count,vertex_count,failure\n";

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
    out << id << ',' << row.runtime_seconds << ',';
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
  std::sort(ordered.begin(), ordered.end(),
    [](const SectionStats& a, const SectionStats& b) { return a.section_id < b.section_id; });
  for (const SectionStats& row : ordered)
  {
    write_row(std::to_string(row.section_id), row, /*is_total=*/false);
  }
  write_row("total", totals_, /*is_total=*/true);

  KINDS_INFO("Statistics: wrote meshing CSV to " << unique_path.generic_string() << " (" << sections_.size()
                                                 << " section row(s) + total)");
  return true;
}

bool Statistics::writeEventListCsv(const std::filesystem::path& path) const
{
  // Match statistics naming (experiment tag + section count) then swap statistics → event_list.
  const std::filesystem::path unique_path
    = timestampedCsvPath(eventListCsvPathBeside(statisticsCsvStemPath(path)));
  std::ofstream out(unique_path);
  if (!out)
  {
    KINDS_WARNING("Statistics: failed to open event-list CSV " << unique_path.generic_string());
    return false;
  }

  out << "event_id,occurrence_t,occurrence_infinitesimal_t,event_type,half_edge_id,delaunay_edge_id,"
         "voronoi_vertex_id,section_id,strand_id,parent_component_id,split_time,target_inside,radius_shift\n";
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

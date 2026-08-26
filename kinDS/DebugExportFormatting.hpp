#pragma once

#include "KineticAlgorithm.hpp"

#include <cmath>
#include <iomanip>
#include <sstream>
#include <string>

namespace kinDS
{
/// Decimal digits for kinetic times (real and non-zero infinitesimal) in debug export filenames.
/// Single source of truth so Explorer lexicographic sort stays consistent across exporters.
inline constexpr int kDebugExportTimePrecision = 17;

/// Formats @p t with @ref kDebugExportTimePrecision in fixed notation.
inline std::string formatDebugExportTime(double t)
{
  if (!std::isfinite(t))
  {
    return "unknown";
  }
  std::ostringstream o;
  o << std::fixed << std::setprecision(kDebugExportTimePrecision) << t;
  return o.str();
}

/// Infinitesimal suffix `_i0!` (zero; `!` sorts before `.` alphanumerically) or `_i{…}` with
/// @ref kDebugExportTimePrecision digits for non-zero values.
inline std::string formatDebugExportInfinitesimalSuffix(double infinitesimal_time)
{
  if (!std::isfinite(infinitesimal_time) || infinitesimal_time == 0.0)
  {
    return "_i0!";
  }
  return "_i" + formatDebugExportTime(infinitesimal_time);
}

/// Filename / path token `t{time}` or `t_unknown` (real time only; used by non-SVG exporters).
inline std::string formatDebugExportTimeToken(double t)
{
  if (!std::isfinite(t))
  {
    return "t_unknown";
  }
  return "t" + formatDebugExportTime(t);
}

/// SVG / event filename token `t{real}_i0!` or `t{real}_i{…}` (same digit count as real time).
inline std::string formatDebugExportTimeToken(EventTime t)
{
  return formatDebugExportTimeToken(t.real_time) + formatDebugExportInfinitesimalSuffix(t.infinitesimal_time);
}
} // namespace kinDS

#pragma once

#include <cmath>
#include <iomanip>
#include <sstream>
#include <string>

namespace kinDS
{
/// Decimal digits for kinetic times in debug export filenames and path tokens.
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

/// Filename / path token `t{time}` or `t_unknown`.
inline std::string formatDebugExportTimeToken(double t)
{
  if (!std::isfinite(t))
  {
    return "t_unknown";
  }
  return "t" + formatDebugExportTime(t);
}
} // namespace kinDS

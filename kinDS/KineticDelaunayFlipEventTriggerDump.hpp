#pragma once

#include "KineticDelaunay.hpp"
#include "Polynomial.hpp"
#include "Trajectory.hpp"

#include <array>
#include <filesystem>
#include <string>
#include <vector>

namespace kinDS
{
struct FlipEventTriggerDump
{
  size_t section = 0;
  double schedule_time = 0.0;
  double section_fraction = 0.0;
  size_t half_edge_id = 0;
  bool boundary_ccw = false;
  std::array<int, 4> vertex_ids { -1, -1, -1, -1 };
  size_t vertex_count = 0;
  std::vector<Trajectory<2>> site_trajectories;
  Polynomial event_trigger;
};

FlipEventTriggerDump buildFlipEventTriggerDump(const KineticDelaunay& kd, size_t half_edge_id, double schedule_time);

void writeFlipEventTriggerPolynomialDump(const KineticDelaunay& kd, const FlipEventTriggerDump& dump,
  const std::filesystem::path& filepath, double occurrence_time);

bool shouldDumpFlipPolynomialsForEvent(
  const KineticDelaunay& kd, double occurrence_time, size_t half_edge_id, double time_margin = 1e-4);

} // namespace kinDS

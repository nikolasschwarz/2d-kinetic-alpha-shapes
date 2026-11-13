#include "kinDS/Statistics.hpp"

#include <catch2/catch_test_macros.hpp>

using namespace kinDS;

TEST_CASE("Statistics rows are per section and runtime branch, sorted on endRun", "[Statistics]")
{
  Statistics stats;
  stats.beginRun();
  stats.onEvent(KineticEventType::Flip, 1.2, 2);
  stats.onEvent(KineticEventType::Flip, 0.4, 2);
  stats.onEvent(KineticEventType::Flip, 0.1, 1);
  stats.onEvent(KineticEventType::Section, 0.0);
  stats.endRun();

  const auto& rows = stats.sections();
  REQUIRE(rows.size() == 4);

  REQUIRE(rows[0].section_id == 0);
  REQUIRE_FALSE(rows[0].runtime_branch_id.has_value());
  REQUIRE(rows[0].event_counts[static_cast<size_t>(KineticEventType::Section)] == 1);

  REQUIRE(rows[1].section_id == 0);
  REQUIRE(rows[1].runtime_branch_id == 1);
  REQUIRE(rows[1].event_counts[static_cast<size_t>(KineticEventType::Flip)] == 1);

  REQUIRE(rows[2].section_id == 0);
  REQUIRE(rows[2].runtime_branch_id == 2);
  REQUIRE(rows[2].event_counts[static_cast<size_t>(KineticEventType::Flip)] == 1);

  REQUIRE(rows[3].section_id == 1);
  REQUIRE(rows[3].runtime_branch_id == 2);
  REQUIRE(rows[3].event_counts[static_cast<size_t>(KineticEventType::Flip)] == 1);
}

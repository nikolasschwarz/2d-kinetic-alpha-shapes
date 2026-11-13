#include "kinDS/ClosingMeshDebugSvg.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include <algorithm>
#include <string>

using namespace kinDS;
using Catch::Matchers::ContainsSubstring;
using Anchor = ClosingMeshWalkDebug::Anchor;

namespace
{
std::string joinedNotes(const HalfEdgeDelaunayGraphToSVG::DebugOverlay& overlay)
{
  std::string all;
  for (const std::string& note : overlay.notes)
  {
    all += note + "\n";
  }
  return all;
}
} // namespace

TEST_CASE("Closing mesh walk overlay labels strips, steps, boundary edges and the stuck point", "[ClosingMeshDebugSvg]")
{
  ClosingMeshWalkDebug debug;
  debug.strand_id = 4;
  debug.t = 35.0;
  debug.seed_index = 0;
  debug.reason = "Boundary walk reached corner site 3.";
  debug.strips = {
    { 0, { Anchor::Kind::VoronoiVertex, 23 }, { Anchor::Kind::Crossing, 13, 13 }, "VV 23", "X d13/v13", "ve 13",
      ClosingMeshWalkDebug::StripState::Walked },
    { 1, { Anchor::Kind::VoronoiVertex, 7 }, { Anchor::Kind::Crossing, 2, 9 }, "VV 7", "X d2/v9", "ve 9",
      ClosingMeshWalkDebug::StripState::Unused },
  };
  debug.walk_strip_order = { 0 };
  debug.walk_points = { { { Anchor::Kind::VoronoiVertex, 23 }, "VV 23" },
    { { Anchor::Kind::Crossing, 13, 13 }, "X d13/v13" }, { { Anchor::Kind::Site, 3 }, "site 3" } };
  debug.boundary_edges = { { 27, 5, 4 }, { 59, 4, 3 } };
  debug.stuck = ClosingMeshWalkDebug::WalkPoint { { Anchor::Kind::Site, 3 }, "corner site 3" };

  const auto overlay = buildClosingMeshWalkOverlay(debug);

  // Two boundary edges, two strips and the walk path; the walk numbers its steps.
  REQUIRE(overlay.paths.size() == 5);
  CHECK(overlay.paths[0].label == "b0 he27");
  CHECK(overlay.paths[1].label == "b1 he59");
  CHECK(overlay.paths[3].label == "s1");
  CHECK(overlay.paths[3].rgb == std::array<int, 3> { 214, 39, 40 });
  CHECK(overlay.paths[4].number_anchors);
  CHECK(overlay.paths[4].anchors.size() == 3);
  REQUIRE(overlay.markers.size() == 1);
  CHECK(overlay.markers[0].label == "STUCK: corner site 3");
  CHECK(overlay.markers[0].anchor.kind == Anchor::Kind::Site);

  const std::string notes = joinedNotes(overlay);
  CHECK_THAT(notes, ContainsSubstring("strand 4"));
  CHECK_THAT(notes, ContainsSubstring("walk order: s0"));
  CHECK_THAT(notes, ContainsSubstring("s1 [UNUSED] ve 9: VV 7 -> X d2/v9"));
  CHECK_THAT(notes, ContainsSubstring("2=site 3"));
}

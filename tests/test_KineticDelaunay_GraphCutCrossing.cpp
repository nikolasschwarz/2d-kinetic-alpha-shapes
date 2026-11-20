#include "kinDS/Logger.hpp"
#include "kinDS/StrandTree.hpp"
#include "kinDS/TreeMesher.hpp"

#include <catch2/catch_test_macros.hpp>

#include <filesystem>
#include <string>

using namespace kinDS;

namespace
{
std::filesystem::path oak10FixturePath()
{
  const char* relative = "tests/fixtures/test_oak_10_branch_241_sections_126_134.txt";
  for (const std::filesystem::path prefix : { ".", "..", "../.." })
  {
    if (std::filesystem::exists(prefix / relative))
    {
      return prefix / relative;
    }
  }
  throw std::runtime_error(std::string("test fixture not found: ") + relative);
}
} // namespace

// test_oak_10, input branch 241 at global sections 126-134. A separation cut at local t=5.292 moves a Voronoi vertex
// out of a tombstoned triangle onto the new outer face. Its crossing events used to stay keyed to the dead triangle,
// so it re-crossed the hull edge unnoticed and the crossing left its segment (later breaking the closing caps).
TEST_CASE("Graph cut reschedules crossing events of relocated Voronoi vertices", "[KineticDelaunay][CrossingData]")
{
  StrandTree tree = StrandTree::loadFromFile(oak10FixturePath());
  REQUIRE(tree.getPoints().size() == 16);

  size_t off_segment = 0;
  logger.setExternalSink(
    [&](LogLevel, const std::string& message)
    {
      if (message.find("intersection not on both segments") != std::string::npos)
      {
        ++off_segment;
      }
    });

  // Default settings (alpha 5, exact hinge); the separation that triggers the cut needs the SegmentBuilder callbacks.
  TreeMesher mesher(tree);
  TreeMesher::Settings settings;
  settings.check_delaunay_voronoi_intersections = true;
  mesher.setSettings(settings);
  REQUIRE_NOTHROW(mesher.runMeshingAlgorithm(false));
  logger.setExternalSink(nullptr);

  REQUIRE(off_segment == 0);
}

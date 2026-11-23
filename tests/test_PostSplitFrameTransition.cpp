#include "kinDS/KineticDelaunay.hpp"
#include "kinDS/StrandTree.hpp"

#include <catch2/catch_test_macros.hpp>

#include <filesystem>
#include <string>

using namespace kinDS;

namespace
{
std::filesystem::path fixturePath(const char* relative)
{
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

// A transition registered while another one hands over must start its last hold piece in the frame the previous
// piece ended in, so the section state is identical before and after the section event.
TEST_CASE("Hold pieces start in the handover frame only for the last hold piece", "[KineticDelaunay][FrameTransition]")
{
  PostSplitFrameTransition transition;
  transition.common_reference_branch = 6;
  transition.hold_end_height = 32;
  transition.start_reference_branches.emplace(52, 1);

  CHECK(holdPieceStartReference(transition, 52, 31) == 1);
  CHECK(holdPieceStartReference(transition, 52, 30) == 6);
  CHECK(holdPieceStartReference(transition, 53, 31) == 6);
}

// test_oak_10, input branch 257 at global sections 131-150 (tree top). Strands 1105/1106 leave for branch 277 at
// global height 149 and blend into its frame on the last section. At the end height the split is never cut, so the
// reference pool still says branch 272; positions there used to jump ~10 units away from where the blend piece ended,
// throwing the top cap's Voronoi vertices far outside the strands.
TEST_CASE("Positions stay continuous at the end of a post-split frame blend", "[KineticDelaunay][FrameTransition]")
{
  StrandTree tree = StrandTree::loadFromFile(fixturePath("tests/fixtures/test_oak_10_branch_257_sections_131_150.txt"));
  REQUIRE(tree.getPoints().size() == 8);

  KineticDelaunay kd(tree, 5.0, false);
  kd.init();
  REQUIRE_NOTHROW(kd.compute());

  const double end = static_cast<double>(tree.getHeight());
  for (size_t strand_id = 0; strand_id < tree.getPoints().size(); ++strand_id)
  {
    INFO("strand " << strand_id);
    CHECK(glm::distance(kd.getPointAt(strand_id, end), kd.getPointAt(strand_id, end - 1e-9)) < 1e-6);
  }
}

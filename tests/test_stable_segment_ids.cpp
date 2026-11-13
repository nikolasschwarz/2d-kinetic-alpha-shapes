#include <catch2/catch_test_macros.hpp>

#include "kinDS/SegmentBuilder.hpp"
#include "kinDS/StrandTree.hpp"

#include <vector>

using namespace kinDS;

namespace
{
StrandTree makeTinyStrandTree()
{
  // Two strands, three sections of support; strand0 has 1 subdiv -> 2 segments; strand1 has none -> 1 segment.
  std::vector<std::vector<glm::dvec2>> support = {
    { { 0.0, 0.0 }, { 0.0, 1.0 }, { 0.0, 2.0 } },
    { { 1.0, 0.0 }, { 1.0, 1.0 }, { 1.0, 2.0 } },
  };
  std::vector<std::vector<double>> subdivs = {
    { 0.5 },
    {},
  };
  std::vector<std::vector<int>> physics = {
    { 10, 11 },
    { 20 },
  };
  std::vector<std::vector<glm::dmat4>> transforms(3, std::vector<glm::dmat4>(1, glm::dmat4(1.0)));
  std::vector<std::vector<size_t>> branch_indices = {
    { 0, 0, 0 },
    { 0, 0, 0 },
  };
  std::vector<std::vector<std::vector<size_t>>> strands_by_branch = {
    { { 0, 1 } },
    { { 0, 1 } },
    { { 0, 1 } },
  };
  return StrandTree(support, subdivs, physics, transforms, branch_indices, strands_by_branch);
}
} // namespace

TEST_CASE("expectedSegmentCountForStrand uses physics map when present", "[stable_ids]")
{
  const StrandTree tree = makeTinyStrandTree();
  REQUIRE(SegmentBuilder::expectedSegmentCountForStrand(tree, 0) == 2);
  REQUIRE(SegmentBuilder::expectedSegmentCountForStrand(tree, 1) == 1);
}

TEST_CASE("expectedSegmentCountForStrand falls back to 1+subdivisions", "[stable_ids]")
{
  std::vector<std::vector<glm::dvec2>> support = {
    { { 0.0, 0.0 }, { 0.0, 1.0 } },
  };
  std::vector<std::vector<double>> subdivs = { { 0.25, 0.75 } };
  std::vector<std::vector<int>> physics = { {} }; // empty -> fallback
  std::vector<std::vector<glm::dmat4>> transforms(2, std::vector<glm::dmat4>(1, glm::dmat4(1.0)));
  std::vector<std::vector<size_t>> branch_indices = { { 0, 0 } };
  std::vector<std::vector<std::vector<size_t>>> strands_by_branch = { { { 0 } }, { { 0 } } };
  const StrandTree tree(support, subdivs, physics, transforms, branch_indices, strands_by_branch);
  REQUIRE(SegmentBuilder::expectedSegmentCountForStrand(tree, 0) == 3);
}

TEST_CASE("BranchScheduler defaults", "[branch_scheduler]")
{
  BranchScheduler scheduler;
  REQUIRE(scheduler.runtime_branch_id == BranchScheduler::no_branch);
  REQUIRE_FALSE(scheduler.independent);
  REQUIRE(scheduler.queue.empty());
  REQUIRE(scheduler.sections_advanced == 0);
}

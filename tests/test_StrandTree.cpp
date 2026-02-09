#include "kinDS/StrandTree.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <filesystem>
#include <vector>

using namespace kinDS;

TEST_CASE("StrandTree serialization", "[StrandTree]")
{
  // Create a simple test StrandTree
  std::vector<std::vector<glm::dvec2>> support_points = {
    { { 0.0, 0.0 }, { 1.0, 1.0 }, { 2.0, 0.0 } }, // strand 0
    { { 0.0, 1.0 }, { 1.0, 2.0 } }                 // strand 1
  };

  std::vector<std::vector<double>> subdivisions_by_strand = {
    { 0.0, 0.5, 1.0 }, // strand 0
    { 0.0, 1.0 }        // strand 1
  };

  std::vector<std::vector<int>> physics_strand_to_segment_indices = {
    { 0, 1 }, // strand 0
    { 0 }     // strand 1
  };

  // Create identity matrices for transforms
  std::vector<std::vector<glm::mat4>> transforms_by_height_and_branch = {
    { glm::mat4(1.0) }, // height 0, branch 0
    { glm::mat4(1.0) }  // height 1, branch 0
  };

  std::vector<std::vector<size_t>> branch_indices = {
    { 0, 0 }, // strand 0
    { 0 }     // strand 1
  };

  std::vector<std::vector<std::vector<size_t>>> strands_by_branch_id = {
    { { 0, 1 } }, // height 0: branch 0 contains strands 0, 1
    { { 0 } }     // height 1: branch 0 contains strand 0
  };

  StrandTree original_tree(support_points, subdivisions_by_strand, physics_strand_to_segment_indices,
                           transforms_by_height_and_branch, branch_indices, strands_by_branch_id);

  // Test save and load
  const std::string test_file = "test_strandtree.tmp";
  
  SECTION("Save and load preserves data")
  {
    // Save the tree
    original_tree.saveToFile(test_file);
    
    // Verify file was created
    REQUIRE(std::filesystem::exists(test_file));
    
    // Load the tree
    StrandTree loaded_tree = StrandTree::loadFromFile(test_file);
    
    // Verify basic properties
    REQUIRE(loaded_tree.getHeight() == original_tree.getHeight());
    REQUIRE(loaded_tree.getPoints().size() == original_tree.getPoints().size());
    
    // Verify support points
    const auto& original_points = original_tree.getPoints();
    const auto& loaded_points = loaded_tree.getPoints();
    REQUIRE(original_points.size() == loaded_points.size());
    
    for (size_t i = 0; i < original_points.size(); ++i)
    {
      REQUIRE(original_points[i].size() == loaded_points[i].size());
      for (size_t j = 0; j < original_points[i].size(); ++j)
      {
        REQUIRE(original_points[i][j].x == Catch::Approx(loaded_points[i][j].x));
        REQUIRE(original_points[i][j].y == Catch::Approx(loaded_points[i][j].y));
      }
    }
    
    // Verify subdivisions
    const auto& original_subdivs = original_tree.getSubdivisionsByStrand();
    const auto& loaded_subdivs = loaded_tree.getSubdivisionsByStrand();
    REQUIRE(original_subdivs.size() == loaded_subdivs.size());
    for (size_t i = 0; i < original_subdivs.size(); ++i)
    {
      REQUIRE(original_subdivs[i].size() == loaded_subdivs[i].size());
      for (size_t j = 0; j < original_subdivs[i].size(); ++j)
      {
        REQUIRE(original_subdivs[i][j] == Catch::Approx(loaded_subdivs[i][j]));
      }
    }
    
    // Clean up
    std::filesystem::remove(test_file);
  }

  SECTION("Evaluate function works correctly")
  {
    // Test evaluation at t=0.5 for strand 0 (should interpolate between points 0 and 1)
    glm::dvec2 result = original_tree.evaluate(0, 0.5);
    REQUIRE(result.x == Catch::Approx(0.5));
    REQUIRE(result.y == Catch::Approx(0.5));
    
    // Test evaluation at t=1.0 for strand 0 (should be point at index 1)
    result = original_tree.evaluate(0, 1.0);
    REQUIRE(result.x == Catch::Approx(1.0));
    REQUIRE(result.y == Catch::Approx(1.0));
  }
}

TEST_CASE("StrandTree empty tree", "[StrandTree]")
{
  // Test with empty vectors
  std::vector<std::vector<glm::dvec2>> empty_points;
  std::vector<std::vector<double>> empty_subdivs;
  std::vector<std::vector<int>> empty_physics;
  std::vector<std::vector<glm::mat4>> empty_transforms;
  std::vector<std::vector<size_t>> empty_branches;
  std::vector<std::vector<std::vector<size_t>>> empty_strands_by_branch;

  StrandTree empty_tree(empty_points, empty_subdivs, empty_physics, empty_transforms, empty_branches,
                        empty_strands_by_branch);

  REQUIRE(empty_tree.getHeight() == 0);
  REQUIRE(empty_tree.getPoints().empty());
}

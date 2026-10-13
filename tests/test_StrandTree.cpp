#include "kinDS/StrandTree.hpp"
#include "kinDS/Logger.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <filesystem>
#include <fstream>
#include <limits>
#include <string>
#include <vector>

using namespace kinDS;

static void enable_all_log_levels_for_test()
{
  logger.setLogLevelMask(
    LogLevel::Debug | LogLevel::Info | LogLevel::Warning | LogLevel::Error | LogLevel::Critical);
}

TEST_CASE("StrandTree serialization", "[StrandTree]")
{
  enable_all_log_levels_for_test();
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
  std::vector<std::vector<glm::dmat4>> transforms_by_height_and_branch = {
    { glm::dmat4(1.0) }, // height 0, branch 0
    { glm::dmat4(1.0) }  // height 1, branch 0
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

    // Verify transforms (row-major file order: m[col][row] for row, col in 0..3)
    const auto& original_transforms = original_tree.getTransformsByHeightAndBranch();
    const auto& loaded_transforms = loaded_tree.getTransformsByHeightAndBranch();
    REQUIRE(original_transforms.size() == loaded_transforms.size());
    for (size_t h = 0; h < original_transforms.size(); ++h)
    {
      REQUIRE(original_transforms[h].size() == loaded_transforms[h].size());
      for (size_t b = 0; b < original_transforms[h].size(); ++b)
      {
        for (int row = 0; row < 4; ++row)
        {
          for (int col = 0; col < 4; ++col)
          {
            REQUIRE(original_transforms[h][b][col][row] == Catch::Approx(loaded_transforms[h][b][col][row]));
          }
        }
      }
    }

    {
      std::ifstream in(test_file);
      std::string contents((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
      REQUIRE(contents.find("[Section] support_points") != std::string::npos);
      REQUIRE(contents.find("[Section] transforms_by_height_and_branch") != std::string::npos);
      REQUIRE(contents.find("# kinDS StrandTree text format") != std::string::npos);
    }

    // Clean up
    std::filesystem::remove(test_file);
  }

  SECTION("Legacy format without section prefix or comments still loads")
  {
    const std::string legacy_file = "test_strandtree_legacy.tmp";
    {
      std::ofstream out(legacy_file);
      out << "StrandTree 1\n2\n2\n";
      out << "support_points\n3\n0 0\n1 1\n2 0\n2\n0 1\n1 2\n";
      out << "subdivisions_by_strand\n3 0 0.5 1\n2 0 1\n";
      out << "physics_strand_to_segment_indices\n2 0 1\n1 0\n";
      out << "transforms_by_height_and_branch\n2\n1\n";
      out << "1 0 0 0 0 1 0 0 0 0 1 0 0 0 0 1\n";
      out << "1\n1 0 0 0 0 1 0 0 0 0 1 0 0 0 0 1\n";
      out << "branch_indices\n2 0 0\n1 0\n";
      out << "strands_by_branch_id\n2\n1\n2 0 1\n1\n1 0\n";
    }

    StrandTree loaded = StrandTree::loadFromFile(legacy_file);
    REQUIRE(loaded.getHeight() == original_tree.getHeight());
    REQUIRE(loaded.getPoints().size() == original_tree.getPoints().size());
    std::filesystem::remove(legacy_file);
  }

  SECTION("Inline comments and blank lines are ignored on load")
  {
    const std::string commented_file = "test_strandtree_comments.tmp";
    {
      std::ofstream out(commented_file);
      out << "# header comment\n";
      out << "StrandTree 1 # format tag\n";
      out << "\n2 # height\n";
      out << "2 # strands\n";
      out << "[Section] support_points\n";
      out << "3\n0 0 # base\n1 1\n2 0\n";
      out << "2\n0 1\n1 2\n";
      out << "subdivisions_by_strand\n";
      out << "3 0 0.5 1\n2 0 1\n";
      out << "physics_strand_to_segment_indices\n";
      out << "2 0 1\n1 0\n";
      out << "transforms_by_height_and_branch\n";
      out << "2\n1\n";
      out << "1 0 0 0 0 1 0 0 0 0 1 0 0 0 0 1 # identity\n";
      out << "1\n1 0 0 0 0 1 0 0 0 0 1 0 0 0 0 1\n";
      out << "branch_indices\n2 0 0\n1 0\n";
      out << "strands_by_branch_id\n2\n1\n2 0 1\n1\n1 0\n";
    }

    StrandTree loaded = StrandTree::loadFromFile(commented_file);
    REQUIRE(loaded.getPoints()[0][0].x == Catch::Approx(0.0));
    REQUIRE(loaded.getPoints()[0][1].y == Catch::Approx(1.0));
    std::filesystem::remove(commented_file);
  }

  SECTION("Transform matrices round-trip without merged tokens")
  {
    glm::dmat4 m(1.0);
    m[0][0] = 0.000999994;
    m[1][0] = 0.0;
    m[2][0] = 6.86638e-10;
    m[3][0] = 1.13621e-70;
    m[0][1] = 0.01148234;
    m[1][1] = 1.0;
    m[2][1] = 7.1555e-10;
    m[3][1] = 0.0;

    std::vector<std::vector<glm::dmat4>> transforms = { { m }, { m } };
    std::vector<std::vector<size_t>> branches = { { 0, 0, 0 }, { 0, 0 } };
    std::vector<std::vector<std::vector<size_t>>> strands_by_branch = { { { 0, 1 } }, { { 0 } } };

    StrandTree tree(support_points, subdivisions_by_strand, physics_strand_to_segment_indices, transforms,
      branches, strands_by_branch);

    tree.saveToFile(test_file);
    StrandTree loaded = StrandTree::loadFromFile(test_file);

    const glm::dmat4& loaded_m = loaded.getTransformsByHeightAndBranch()[0][0];
    for (int row = 0; row < 4; ++row)
    {
      for (int col = 0; col < 4; ++col)
      {
        REQUIRE(m[col][row] == Catch::Approx(loaded_m[col][row]));
      }
    }

    {
      std::ifstream in(test_file);
      std::string line;
      bool in_transforms = false;
      while (std::getline(in, line))
      {
        if (line.find("transforms_by_height_and_branch") != std::string::npos)
        {
          in_transforms = true;
          continue;
        }
        if (!in_transforms || line.empty() || line[0] == '#')
        {
          continue;
        }
        if (line.find('e') != std::string::npos && line.find(' ') != std::string::npos)
        {
          REQUIRE(line.find("e-700.") == std::string::npos);
          REQUIRE(line.find("e-70 1.") != std::string::npos);
          break;
        }
      }
    }

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

TEST_CASE("StrandTree evaluateTransformed before branch exists", "[StrandTree]")
{
  enable_all_log_levels_for_test();

  // Branch 1 appears only at height 2; heights 0–1 have a single branch.
  std::vector<std::vector<glm::dvec2>> support_points = {
    { { 0.0, 0.0 }, { 1.0, 0.0 }, { 2.0, 0.0 } },
  };
  std::vector<std::vector<double>> subdivisions = { { 0.0, 0.5, 1.0 } };
  std::vector<std::vector<int>> physics = { { 0 } };

  glm::dmat4 branch0 = glm::dmat4(1.0);
  glm::dmat4 branch1 = glm::translate(glm::dmat4(1.0), glm::dvec3(10.0, 0.0, 0.0));

  std::vector<std::vector<glm::dmat4>> transforms = {
    { branch0 },
    { branch0 },
    { branch0, branch1 },
  };
  std::vector<std::vector<size_t>> branch_indices = {
    { 0, 0, 1 },
  };
  std::vector<std::vector<std::vector<size_t>>> strands_by_branch = {
    { { 0 } },
    { { 0 } },
    { { 0 }, { 0 } },
  };

  StrandTree tree(support_points, subdivisions, physics, transforms, branch_indices, strands_by_branch);

  // t=1.5: lower section has only branch 0, reference branch 1 does not exist yet at that height.
  const glm::dvec2 at_lower = tree.getPointTransformed(0, 1, 1);
  REQUIRE(at_lower.x == Catch::Approx(1.0));
  REQUIRE(at_lower.y == Catch::Approx(0.0));

  const glm::dvec2 interpolated = tree.evaluateTransformed(0, 1.5, 1);
  REQUIRE(interpolated.x == Catch::Approx(1.5));
  REQUIRE(interpolated.y == Catch::Approx(0.0));

  const Trajectory<2> piece = tree.getPiecePolynomial(0, 1, 1);
  REQUIRE(interpolated.x == Catch::Approx(piece[0](0.5)));
  REQUIRE(interpolated.y == Catch::Approx(piece[1](0.5)));
}

TEST_CASE("StrandTree getPiecePolynomial matches evaluateTransformed", "[StrandTree]")
{
  enable_all_log_levels_for_test();

  glm::dmat4 branch0 = glm::dmat4(1.0);
  glm::dmat4 branch1 = glm::translate(glm::dmat4(1.0), glm::dvec3(10.0, 0.0, 0.0));

  std::vector<std::vector<glm::dvec2>> support_points = {
    { { 0.0, 0.0 }, { 1.0, 0.0 }, { 2.0, 0.0 } },
    { { 0.0, 1.0 }, { 0.0, 2.0 }, { 0.0, 3.0 } },
  };
  std::vector<std::vector<double>> subdivisions = { { 0.0, 1.0, 2.0 }, { 0.0, 1.0, 2.0 } };
  std::vector<std::vector<int>> physics = { { 0 }, { 0 } };
  std::vector<std::vector<glm::dmat4>> transforms = {
    { branch0, branch1 },
    { branch0, branch1 },
    { branch0, branch1 },
  };
  std::vector<std::vector<size_t>> branch_indices = {
    { 0, 0, 1 },
    { 1, 1, 1 },
  };
  std::vector<std::vector<std::vector<size_t>>> strands_by_branch = {
    { { 0 }, { 1 } },
    { { 0 }, { 0, 1 } },
    { { 0 }, { 0, 1 } },
  };

  StrandTree tree(support_points, subdivisions, physics, transforms, branch_indices, strands_by_branch);
  const size_t reference_branch = 1;

  for (double t : { 1.0, 1.25, 1.5, 1.75, 2.0 })
  {
    const glm::dvec2 evaluated = tree.evaluateTransformed(1, t, reference_branch);
    const size_t section = static_cast<size_t>(std::floor(t));
    const double frac = t - static_cast<double>(section);
    if (frac < std::numeric_limits<double>::epsilon())
    {
      const glm::dvec2 knot = tree.getPointTransformed(1, section, reference_branch);
      REQUIRE(evaluated.x == Catch::Approx(knot.x));
      REQUIRE(evaluated.y == Catch::Approx(knot.y));
      continue;
    }

    const Trajectory<2> piece = tree.getPiecePolynomial(1, section, reference_branch);
    REQUIRE(evaluated.x == Catch::Approx(piece[0](frac)));
    REQUIRE(evaluated.y == Catch::Approx(piece[1](frac)));
  }
}

TEST_CASE("StrandTree getBranchIndex clamps height for short strands", "[StrandTree]")
{
  enable_all_log_levels_for_test();

  // Global tree height is 3, but strand 1 only has branch entries through height 1.
  std::vector<std::vector<glm::dvec2>> support_points = {
    { { 0.0, 0.0 }, { 1.0, 0.0 }, { 2.0, 0.0 }, { 3.0, 0.0 } },
    { { 0.0, 1.0 }, { 1.0, 1.0 } },
  };
  std::vector<std::vector<double>> subdivisions = { { 0.0, 1.0, 2.0, 3.0 }, { 0.0, 1.0 } };
  std::vector<std::vector<int>> physics = { { 0 }, { 0 } };
  std::vector<std::vector<glm::dmat4>> transforms = {
    { glm::dmat4(1.0) },
    { glm::dmat4(1.0) },
    { glm::dmat4(1.0) },
    { glm::dmat4(1.0) },
  };
  std::vector<std::vector<size_t>> branch_indices = {
    { 0, 1, 2, 2 },
    { 0, 1 },
  };
  std::vector<std::vector<std::vector<size_t>>> strands_by_branch = {
    { { 0, 1 } },
    { { 0 }, { 1 } },
    { { 0 }, { 1 } },
    { { 0 }, { 1 } },
  };

  StrandTree tree(support_points, subdivisions, physics, transforms, branch_indices, strands_by_branch);

  REQUIRE(tree.getBranchIndex(1, 1) == 1);
  REQUIRE(tree.getBranchIndex(1, 2) == 1);
  REQUIRE(tree.getBranchIndex(1, 99) == 1);
}

TEST_CASE("StrandTree empty tree", "[StrandTree]")
{
  enable_all_log_levels_for_test();
  // Test with empty vectors
  std::vector<std::vector<glm::dvec2>> empty_points;
  std::vector<std::vector<double>> empty_subdivs;
  std::vector<std::vector<int>> empty_physics;
  std::vector<std::vector<glm::dmat4>> empty_transforms;
  std::vector<std::vector<size_t>> empty_branches;
  std::vector<std::vector<std::vector<size_t>>> empty_strands_by_branch;

  StrandTree empty_tree(empty_points, empty_subdivs, empty_physics, empty_transforms, empty_branches,
                        empty_strands_by_branch);

  REQUIRE(empty_tree.getHeight() == 0);
  REQUIRE(empty_tree.getPoints().empty());
}

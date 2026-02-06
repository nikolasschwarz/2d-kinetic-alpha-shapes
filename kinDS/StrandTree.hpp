#pragma once
#include "PlaneProjector.hpp"
#include "Polynomial.hpp"

#include <array>
#include <string>
#include <vector>

namespace kinDS
{

template<size_t dim> using Trajectory = std::array<Polynomial, dim>;

/**
 * This class handles the trajectories of strands according to branches, allowing to easily get points in a different
 * frame of reference as needed
 */
class StrandTree
{
 private:
  std::vector<std::vector<glm::dvec2>> support_points;
  std::vector<std::vector<double>> subdivisions_by_strand; // subdivisions for each strand [strand_id][subdivisions]
  std::vector<std::vector<int>> physics_strand_to_segment_indices; // physics strand to segment indices [strand_id]
  std::vector<std::vector<glm::mat4>>
    transforms_by_height_and_branch; // transforms for each height and branch [height][branch_id]
  std::vector<std::vector<size_t>> branch_indices; // branch indices for each strand [strand_id][height]
  std::vector<std::vector<std::vector<size_t>>>
    strands_by_branch_id; // strands by branch id [height][branch_id][strand_no]

  size_t height = 0;

 public:
  StrandTree(const std::vector<std::vector<glm::dvec2>>& support_points,
    const std::vector<std::vector<double>>& subdivisions_by_strand,
    const std::vector<std::vector<int>>& physics_strand_to_segment_indices,
    const std::vector<std::vector<glm::mat4>>& transforms_by_height_and_branch,
    const std::vector<std::vector<size_t>>& branch_indices,
    const std::vector<std::vector<std::vector<size_t>>>& strands_by_branch_id);

  const std::vector<std::vector<glm::dvec2>>& getPoints() const;

  size_t getHeight() const;

  size_t addTrajectory(const std::vector<glm::dvec2>& traj);

  glm::dvec2 evaluate(size_t strand_id, double t) const;

  glm::dvec2 getPointTransformed(size_t strand_id, size_t index, size_t reference_branch) const;

  glm::dvec2 evaluateTransformed(size_t strand_id, double t, size_t reference_branch) const;
  glm::dvec3 getPointInObjectSpace(size_t strand_id, double t) const;

  glm::dvec3 transformToObjectSpace(glm::dvec3& v_3d, size_t strand_id, double t) const;

  // TODO: also adjust to different reference frame
  std::array<Polynomial, 2> getPiecePolynomial(size_t strand_id, size_t index) const;

  // getters with named indices
  const std::vector<glm::dvec2>& getSupportPoints(size_t strand_id) const { return support_points[strand_id]; }
  const std::vector<double>& getSubdivisions(size_t strand_id) const { return subdivisions_by_strand[strand_id]; }
  const std::vector<int>& getPhysicsStrandToSegmentIndices(size_t strand_id) const
  {
    return physics_strand_to_segment_indices[strand_id];
  }
  const glm::mat4& getTransformByHeightAndBranch(size_t height, size_t branch_id) const
  {
    return transforms_by_height_and_branch[height][branch_id];
  }
  const std::vector<size_t>& getBranchIndices(size_t strand_id) const { return branch_indices[strand_id]; }
  size_t getBranchIndex(size_t strand_id, size_t height) const { return branch_indices[strand_id][height]; }
  const std::vector<std::vector<size_t>>& getStrandBranchesByHeight(size_t height) const
  {
    return strands_by_branch_id[height];
  }
  const std::vector<size_t>& getStrandsByBranch(size_t height, size_t branch_id) const
  {
    return strands_by_branch_id[height][branch_id];
  }

  // getters without indices
  const std::vector<std::vector<glm::dvec2>>& getSupportPoints() const { return support_points; }
  const std::vector<std::vector<double>>& getSubdivisionsByStrand() const { return subdivisions_by_strand; }
  const std::vector<std::vector<int>>& getPhysicsStrandToSegmentIndices() const
  {
    return physics_strand_to_segment_indices;
  }
  const std::vector<std::vector<glm::mat4>>& getTransformsByHeightAndBranch() const
  {
    return transforms_by_height_and_branch;
  }
  const std::vector<std::vector<size_t>>& getBranchIndices() const { return branch_indices; }
  const std::vector<std::vector<std::vector<size_t>>>& getStrandsByBranchId() const { return strands_by_branch_id; }
};
}; // namespace kinDS
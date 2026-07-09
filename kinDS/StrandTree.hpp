#pragma once
#include "PlaneProjector.hpp"
#include "Trajectory.hpp"

#include <array>
#include <filesystem>
#include <string>
#include <vector>

namespace kinDS
{

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
  // transforms for each height and branch [height][branch_id]
  std::vector<std::vector<glm::dmat4>> transforms_by_height_and_branch;
  // normal transforms for each height and branch [height][branch_id]
  std::vector<std::vector<glm::dmat4>> normal_transforms_by_height_and_branch;
  std::vector<std::vector<size_t>> branch_indices; // branch indices for each strand [strand_id][height]
  // strands by branch id [height][branch_id][strand_no]
  std::vector<std::vector<std::vector<size_t>>> strands_by_branch_id;

  size_t tree_height = 0;

 public:
  StrandTree(const std::vector<std::vector<glm::dvec2>>& support_points,
    const std::vector<std::vector<double>>& subdivisions_by_strand,
    const std::vector<std::vector<int>>& physics_strand_to_segment_indices,
    const std::vector<std::vector<glm::dmat4>>& transforms_by_height_and_branch,
    const std::vector<std::vector<size_t>>& branch_indices,
    const std::vector<std::vector<std::vector<size_t>>>& strands_by_branch_id);

  const std::vector<std::vector<glm::dvec2>>& getPoints() const;

  size_t getHeight() const;

  /// Number of branches at @p height (section index).
  size_t getBranchCount(size_t height) const;

  size_t addTrajectory(const std::vector<glm::dvec2>& traj);

  glm::dvec2 evaluate(size_t strand_id, double t) const;

  glm::dvec2 getPointTransformed(size_t strand_id, size_t index, size_t reference_branch) const;

  glm::dvec2 evaluateTransformed(size_t strand_id, double t, size_t reference_branch) const;
  glm::dvec3 getPointInObjectSpace(size_t strand_id, double t) const;

  glm::dvec3 transformToObjectSpace(glm::dvec3& v_3d, size_t strand_id, double t) const;
  glm::dvec3 transformToObjectSpace(
    glm::dvec3 v_3d, double t, const std::vector<size_t>& branch_indices_by_height) const;

  /**
   * Linear motion on [@p index, @p index + 1]: support points at those heights are each transformed into the
   * @p reference_branch frame at the same height, then connected by a linear segment in the section parameter.
   */
  Trajectory<2> getPiecePolynomial(size_t strand_id, size_t index, size_t reference_branch) const;

  // getters with named indices
  const std::vector<glm::dvec2>& getSupportPoints(size_t strand_id) const { return support_points[strand_id]; }
  const std::vector<double>& getSubdivisions(size_t strand_id) const { return subdivisions_by_strand[strand_id]; }
  const std::vector<int>& getPhysicsStrandToSegmentIndices(size_t strand_id) const
  {
    return physics_strand_to_segment_indices[strand_id];
  }
  const glm::dmat4& getTransformByHeightAndBranch(size_t height, size_t branch_id) const
  {
    return transforms_by_height_and_branch[height][branch_id];
  }
  const std::vector<size_t>& getBranchIndices(size_t strand_id) const { return branch_indices[strand_id]; }
  /** Input branch id at @p height; clamps @p height to the last entry when the strand ended earlier. */
  size_t getBranchIndex(size_t strand_id, size_t height) const;
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
  const std::vector<std::vector<glm::dmat4>>& getTransformsByHeightAndBranch() const
  {
    return transforms_by_height_and_branch;
  }
  const std::vector<std::vector<size_t>>& getBranchIndices() const { return branch_indices; }
  const std::vector<std::vector<std::vector<size_t>>>& getStrandsByBranchId() const { return strands_by_branch_id; }
  const std::vector<std::vector<glm::dmat4>>& getNormalTransformsByHeightAndBranch() const
  {
    return normal_transforms_by_height_and_branch;
  }

  /** Save all members to a text file. */
  void saveToFile(const std::filesystem::path& path) const;

  /** Load from a text file and return a new StrandTree. */
  static StrandTree loadFromFile(const std::filesystem::path& path);

 private:
  void computeNormalTransforms();

  /**
   * Preprocessing: snap every profile-to-world transform to a similarity (orthogonal * scalar). The source frames
   * carry a small in-plane shear and length anisotropy between the two profile axes; left unconditioned this injects
   * skew whenever points are projected between branch frames (see @ref PlaneProjector). Each transform's basis is
   * replaced by the nearest orthogonal, equal-length frame while preserving the plane, the origin, and the branch's
   * own scale.
   */
  void conditionProfileTransforms();

  glm::dvec2 getPointTransformedAtSection(
    size_t strand_id, size_t index, size_t reference_branch, size_t branch_lookup_height) const;

  /** Valid branch id at @p height; if absent, parent from a strand in @p branch_id at @p branch_lookup_height. */
  size_t resolveBranchAtHeight(size_t height, size_t branch_id, size_t branch_lookup_height) const;
};
}; // namespace kinDS
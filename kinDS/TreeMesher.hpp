#pragma once
#include <glm/glm.hpp>
#include <vector>

namespace kinDS
{

class TreeMesher
{
 public:
  void RunMeshingAlgorithm(const std::vector<std::vector<glm::dvec2>>& support_points,
    std::vector<std::vector<double>>& subdivisions_by_strand,
    std::vector<std::vector<int>>& physics_strand_to_segment_indices,
    const std::vector<std::vector<glm::mat4>>& transforms_by_height_and_branch,
    const std::vector<std::vector<size_t>>& branch_indices,
    std::vector<std::vector<std::vector<size_t>>>& strands_by_branch_id,
    std::vector<float>& boundary_distance_by_segment_id, double alpha_cutoff);
};
}
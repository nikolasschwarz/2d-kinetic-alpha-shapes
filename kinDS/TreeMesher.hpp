#pragma once
#include "StrandTree.hpp"
#include <glm/glm.hpp>
#include <vector>

namespace kinDS
{

class TreeMesher
{
 public:
  void RunMeshingAlgorithm(StrandTree& strand_tree,
    std::vector<float>& boundary_distance_by_segment_id, double alpha_cutoff);
};
}
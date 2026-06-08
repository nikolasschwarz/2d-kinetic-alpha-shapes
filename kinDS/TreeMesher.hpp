#pragma once
#include "MeshIntersection.hpp"
#include "StrandTree.hpp"
#include "VoronoiMesh.hpp"
#include <filesystem>
#include <glm/glm.hpp>
#include <memory>
#include <optional>
#include <string>
#include <vector>

namespace kinDS
{
class SegmentBuilder;
class KineticDelaunay;

class TreeMesher
{
 public:
  struct Settings
  {
    double alpha_cutoff = 10.0; // default value, can be adjusted as needed
    bool fix_missing_meshes = false; // whether to attempt to fix missing meshes by copying from neighbors

    // for debugging purposes:
    bool debug_export_meshes = false;
    size_t max_meshlet_export = size_t(-1); // maximum number of meshlets to export for debugging
    bool merge_meshlets_by_segment = true; // if false, keep one output meshlet per generated mesh pair
  };

 private:
  StrandTree& strand_tree;
  std::shared_ptr<KineticDelaunay> kinetic_delaunay;
  std::shared_ptr<SegmentBuilder> mesh_builder;
  std::vector<VoronoiMesh> segment_meshlets;
  std::vector<std::string> segment_meshlet_export_suffixes;
  std::vector<std::vector<int>> meshing_neighbor_indices; // for each mesh, the neighbor mesh index for each triangle
  std::vector<size_t> meshing_to_physics_segment_indices;
  Settings settings;
  std::function<void(size_t, std::function<void(size_t)>)> parallel_for;

 public:
  TreeMesher(StrandTree& strand_tree);
  TreeMesher(StrandTree& strand_tree, std::function<void(size_t, std::function<void(size_t)>)> parallel_for);
  /// @param visual_debug When true, SegmentBuilder callbacks export debug SVG snapshots during meshing.
  const std::vector<VoronoiMesh>& runMeshingAlgorithm(bool visual_debug = false);
  const VoronoiMesh& getBoundaryMesh() const;
  const std::vector<std::vector<int>>& getMeshingNeighborIndices() const;
  const std::vector<size_t>& getMeshingToPhysicsSegmentIndices() const;
  const std::vector<std::vector<size_t>>& getMeshingStrandToSegmentIndices() const;
  const std::vector<size_t>& getBoundaryVertexToStrandId() const;
  void transformBoundaryMesh(kinDS::VoronoiMesh& boundary_mesh, const glm::dmat4& root_transform = glm::dmat4(1));
  void transformToWorldSpace(
    VoronoiMesh& mesh, size_t strand_id, const glm::dmat4& root_transform = glm::dmat4(1)) const;
  const std::vector<VoronoiMesh>& getSegmentMeshlets() const { return segment_meshlets; }
  Settings& getSettings() { return settings; }
  const Settings& getSettings() const { return settings; }
  void setSettings(const Settings& new_settings) { settings = new_settings; }

  /// Export segment meshlets under @p export_path.
  /// @param export_path Output directory when @p separate_file_per_segment is true; output OBJ path when false.
  /// @param separate_file_per_segment If true, write one OBJ per meshlet; otherwise one OBJ with one group per meshlet.
  /// @param max_exports Maximum meshlets to export; default unlimited.
  void exportCombinedMesh(const std::filesystem::path& export_path, bool separate_file_per_segment,
    std::optional<size_t> max_exports = std::nullopt) const;
  void truncateToBoundary(const VoronoiMesh& boundary_mesh);
  void fixFailedSegments(const MeshIntersection& boundary_intersector);
  std::pair<std::vector<float>, std::vector<float>> computeTopAndBottomBoundaryDistances(
    const std::vector<float>& boundary_distance_by_segment_id);
  void runKineticDelaunay(bool visual_debug);
  void mapMeshingToPhysicsSegmentIndices();
};
}
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

enum class MeshletExportMode
{
  Combined,   ///< One OBJ; each segment meshlet is an OBJ group.
  PerSegment, ///< One OBJ per segment (meshlets merged per segment).
  Raw,        ///< One OBJ per Voronoi-edge mesh pair (no segment-level merge).
};

inline bool mergeMeshletsBySegment(MeshletExportMode mode) { return mode != MeshletExportMode::Raw; }

inline const char* meshletExportModeToString(MeshletExportMode mode)
{
  switch (mode)
  {
  case MeshletExportMode::Combined:
    return "combined";
  case MeshletExportMode::PerSegment:
    return "per_segment";
  case MeshletExportMode::Raw:
    return "raw";
  }
  return "unknown";
}

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

  /// Export meshlets under @p export_path.
  /// @param export_mode How meshlets are grouped into output file(s).
  /// @param export_path Output directory for @c PerSegment/@c Raw; output OBJ path for @c Combined.
  /// @param max_exports Maximum meshlets to export; default unlimited.
  /// @param transformed When true (default), map meshlet vertices from profile space into world space before writing.
  void exportMeshlets(MeshletExportMode export_mode, const std::filesystem::path& export_path, bool transformed = true,
    std::optional<size_t> max_exports = std::nullopt) const;
  void truncateToBoundary(const VoronoiMesh& boundary_mesh);
  void fixFailedSegments(const MeshIntersection& boundary_intersector);
  std::pair<std::vector<float>, std::vector<float>> computeTopAndBottomBoundaryDistances(
    const std::vector<float>& boundary_distance_by_segment_id);
  void runKineticDelaunay(bool visual_debug);
  void mapMeshingToPhysicsSegmentIndices();
};
}
#pragma once
#include "MeshIntersection.hpp"
#include "StrandTree.hpp"
#include "Validator.hpp"
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
    /// When true, a failed meshlet intersection keeps the uncut meshlet; when false, replaces it with an empty mesh.
    bool keep_original_on_intersection_failure = true;
    /// Seam vertices receive segment-meshlet UVs instead of clip-boundary UVs (editor intersection).
    bool intersection_prefer_meshlet_uv_on_seam = false;
    /// Clip-boundary-origin faces use interior (a,b,h) UVs; bark polar distance is treated as r=1.
    bool intersection_boundary_faces_interior_uv = false;
    /// When false, skip JSON vertex/face metadata on meshlets (material_ids still stored for OBJ export).
    /// Default off; CLI @c --store-mesh-metadata / @c --no-store-mesh-metadata.
    /// @ref validate_mesh_vertex_sources implies metadata storage when enabled.
    bool store_mesh_metadata = false;
    /// When true, run @ref Validator::validateAndReportMeshVertexSources after meshing.
    bool validate_mesh_vertex_sources = false;
    /// Output path for @ref Validator when validation is enabled (@ref Validator::defaultLogFilePath by default).
    std::string validate_mesh_vertex_sources_log_path = Validator::defaultLogFilePath();
    /// When false, skip meshlet diagnostic logging and related string assembly (@ref meshletDiagnosticLogLine).
    bool diagnostics = true;
    /// When true, mesh vertices are stored in object space during meshing (@ref SegmentBuilder::create_transformed_mesh).
    /// Export skips a second transform; polygon triangulation still uses profile-plane xy.
    /// When false (@c --untransformed), vertices stay in local profile space (no reference-branch remap, no
    /// kinetic separation), with a final Y/Z swap for viewing alongside SVG.
    bool transform_mesh_at_construction = true;
    /// When true, OBJ export alternates light/dark green and brown materials by even/odd kinetic section.
    bool alternate_section_shading = false;
    /// Inclusive first kinetic section to initialize and process (default 0).
    size_t start_section = 0;
    /// Exclusive kinetic stop / finalize time; empty means tree height (@c StrandTree::getHeight()).
    /// Section events run on `[start_section, end_section)`; events with `t >= end_section` are not processed.
    std::optional<size_t> end_section;
    /// When true, emit a closing-cap meshlet for every live strand at @ref start_section (bootstrap). Default false.
    /// Caps for input branches that finish during the run are always produced regardless of this flag.
    bool mesh_cap_at_start = false;
    /// When true, also emit closing caps for strands still continuing past a premature @ref end_section
    /// (@c --end). Default false. Natural branch endings and the tree-top finalize always get caps regardless.
    bool mesh_cap_at_end = false;
    /// When true, collect per-section runtime / event / topology statistics and write CSV after meshing.
    bool collect_meshing_statistics = false;
    /// Output path used when @ref collect_meshing_statistics is enabled.
    std::filesystem::path meshing_statistics_csv_path = "meshing_statistics.csv";

    // for debugging purposes:
    bool debug_export_meshes = false;
    /// When true, OBJ export emits one @c o object per interior/boundary contributor (iN / bN).
    bool export_separate_contributor_objects = true;
    /// When false, use @ref KineticDelaunay::ComponentSplitPolicy::InPlaceCut at section splits instead of retriangulation.
    bool retriangulate_on_component_split = false;
    size_t max_meshlet_export = size_t(-1); // maximum number of meshlets to export for debugging
    /// When true, write full visual-debug SVGs/TXTs (segmentbuilder snapshots, branch-split dumps, …).
    /// Default off; CLI @c --debug-files enables this (also enables @ref error_files dumps).
    bool visual_debug = false;
    /// When true, write failure SVGs/TXTs (ring-walk / triangulate FAIL, etc.) and the common event-style
    /// kinetic SVG with affected sites/edges/VVs highlighted, without full visual debug.
    /// Default off; CLI @c --error-files enables this. Also enabled by @ref visual_debug.
    bool error_files = false;
    /// When true with @ref visual_debug, SVG folders already split by pending child runtime branches (from the radius
    /// event that notes the split) instead of only after the graph cut. CLI @c --svg-separate-pending-splits.
    bool visual_debug_separate_pending_splits = false;
    std::optional<std::filesystem::path> visual_debug_output_root;
    /// When true, after each kinetic event verify live sites lie inside the graph convex hull.
    /// Default off; CLI @c --check-sites-in-hull.
    bool check_sites_inside_convex_hull = false;
    std::optional<double> flip_polynomial_dump_target_time;
    std::optional<size_t> flip_polynomial_dump_target_half_edge;
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
  bool mesh_vertex_source_validation_passed_ = true;

 public:
  TreeMesher(StrandTree& strand_tree);
  TreeMesher(StrandTree& strand_tree, std::function<void(size_t, std::function<void(size_t)>)> parallel_for);
  /// @param visual_debug When true, SegmentBuilder callbacks export full visual-debug SVG/TXT snapshots.
  /// Prefer setting @ref Settings::visual_debug / @ref Settings::error_files (CLI @c --debug-files / @c --error-files).
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
  std::vector<VoronoiMesh>& getSegmentMeshlets() { return segment_meshlets; }
  std::vector<std::vector<int>>& getMeshingNeighborIndices() { return meshing_neighbor_indices; }
  Settings& getSettings() { return settings; }
  const Settings& getSettings() const { return settings; }
  void setSettings(const Settings& new_settings) { settings = new_settings; }
  bool meshVertexSourceValidationPassed() const { return mesh_vertex_source_validation_passed_; }

  /// Export meshlets under @p export_path.
  /// @param export_mode How meshlets are grouped into output file(s).
  /// @param export_path Output directory for @c PerSegment/@c Raw; output OBJ path for @c Combined.
  /// @param max_exports Maximum meshlets to export; default unlimited.
  /// @param transformed When set, force export-time world/profile mapping. When unset, uses
  ///   @c !settings.transform_mesh_at_construction. When export is untransformed, OBJ filenames are
  ///   suffixed with @c _untransformed.
  void exportMeshlets(MeshletExportMode export_mode, const std::filesystem::path& export_path,
    std::optional<bool> transformed = std::nullopt, std::optional<size_t> max_exports = std::nullopt) const;
  bool meshletsTransformedAtConstruction() const { return settings.transform_mesh_at_construction; }
  /// Clip meshlets against @p boundary_mesh. Returns meshing-segment indices classified as OUTSIDE.
  std::vector<size_t> truncateToBoundary(const VoronoiMesh& boundary_mesh);
  void fixFailedSegments(const MeshIntersection& boundary_intersector);
  std::pair<std::vector<float>, std::vector<float>> computeTopAndBottomBoundaryDistances(
    const std::vector<float>& boundary_distance_by_segment_id);
  void runKineticDelaunay(bool visual_debug);
  void mapMeshingToPhysicsSegmentIndices();
};
}
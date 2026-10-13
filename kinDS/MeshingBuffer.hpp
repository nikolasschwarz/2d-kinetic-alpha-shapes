#pragma once

#include "VoronoiMesh.hpp"

#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <string>
#include <vector>

#include <glm/glm.hpp>

namespace kinDS {

/// Settings that participate in the MeshBuffers cache key (same semantics as EcoSysLab's legacy hasher).
struct MeshingBufferHashSettings {
  bool store_mesh_metadata = false;
  float spline_tension = 1.0f;
  int bark_smooth_iterations = 0;
  float bark_smooth_strength = 0.5f;
  bool bark_smooth_lock_boundary = true;
  bool bark_smooth_uvs = true;
  bool bark_subdivide = false;
  double alpha_cutoff = 10.0;
  double branch_alpha_cutoff = 10.0;
  size_t look_ahead = 0;
  bool hinge_only_profile_plane_mix = false;
};

struct MeshingBufferHashStats {
  std::string input_hash;
  std::string settings_hash;
  std::string root_hash;
  std::string support_hash;
  std::string subdiv_hash;
  std::string physics_hash;
  std::string transforms_hash;
  std::string branch_hash;
  std::string strands_by_branch_hash;
  std::string root_transform_summary;
  size_t support_strand_count = 0;
  size_t support_point_count = 0;
  size_t subdiv_strand_count = 0;
  size_t subdiv_value_count = 0;
  size_t physics_strand_count = 0;
  size_t physics_segment_count = 0;
  size_t transform_height_count = 0;
  size_t transform_matrix_count = 0;
  size_t branch_height_count = 0;
  size_t branch_id_count = 0;
  size_t strands_by_branch_outer_count = 0;
  size_t strands_by_branch_id_count = 0;
  float min_segment_length = 0.f;
  float max_segment_length = 0.f;
};

/// Canonical geometric MeshBuffers payload (KVMG). Optional GPU blobs keep EcoSysLab layout portable.
struct MeshingBufferPayload {
  glm::mat4 root_transform{1.0f};
  std::vector<VoronoiMesh> meshlets;
  std::vector<std::vector<int>> neighbors;
  std::vector<size_t> meshing_to_physics;
  std::vector<std::vector<size_t>> strand_to_segment;

  /// Optional opaque GPU vertex/triangle arrays (EcoSysLab std430 layouts). Empty ⇒ rebuild on load.
  uint32_t gpu_vertex_stride = 0;
  uint32_t gpu_triangle_stride = 0;
  std::vector<std::byte> gpu_vertices;
  std::vector<std::byte> gpu_triangles;
};

struct MeshingBufferMetadata {
  std::string input_hash;
  std::string description;
  std::string created_utc;
  double alpha_cutoff = 10.0;
  double branch_alpha_cutoff = 10.0;
  size_t look_ahead = 0;
  float spline_tension = 1.0f;
  bool store_mesh_metadata = false;
  size_t meshlet_count = 0;
  size_t gpu_vertex_count = 0;
  size_t gpu_triangle_count = 0;
  uint32_t vertex_stride = 0;
  uint32_t triangle_stride = 0;
};

[[nodiscard]] constexpr const char* meshingBufferMagic() {
  return "KVMG";
}
[[nodiscard]] constexpr uint32_t meshingBufferVersion() {
  return 1;
}

[[nodiscard]] inline std::filesystem::path meshingBufferYmlPath(const std::filesystem::path& bin_path) {
  return bin_path.parent_path() / (bin_path.stem().string() + ".yml");
}

[[nodiscard]] MeshingBufferHashStats computeMeshingInputHash(
    const std::vector<std::vector<glm::dvec2>>& support_points,
    const std::vector<std::vector<double>>& subdivisions_by_strand,
    const std::vector<std::vector<int>>& physics_strand_to_segment_indices,
    const std::vector<std::vector<glm::dmat4>>& transforms_by_height_and_branch, const glm::mat4& root_transform,
    const std::vector<std::vector<size_t>>& branch_indices,
    const std::vector<std::vector<std::vector<size_t>>>& strands_by_branch_id,
    const MeshingBufferHashSettings& settings, bool include_alpha_cutoff = true);

/// Writes `@p bin_path` (KVMG) and a sibling `.yml` next to it (legacy EcoSysLab MeshBuffers metadata format).
[[nodiscard]] bool saveMeshingBuffer(const std::filesystem::path& bin_path, const MeshingBufferPayload& payload,
                                     const MeshingBufferMetadata& metadata);

/// Loads KVMG `@p bin_path`. When `@p require_root_match` is true, fails if stored root ≠ `@p expected_root`.
[[nodiscard]] bool loadMeshingBuffer(const std::filesystem::path& bin_path, const glm::mat4& expected_root,
                                     MeshingBufferPayload& out_payload, bool require_root_match = true);

/// Loads top-level fields from a MeshBuffers `.yml` (legacy key set: `hash`, `alpha_cutoff`, ...).
[[nodiscard]] bool loadMeshingBufferMetadata(const std::filesystem::path& yml_path, MeshingBufferMetadata& out_metadata);

}  // namespace kinDS

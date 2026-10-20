#pragma once

#include "VoronoiMesh.hpp"

#include <cmath>
#include <glm/glm.hpp>
#include <mutex>
#include <string>
#include <utility>
#include <vector>

namespace kinDS
{

struct AdjustedBoundaryTriangleUvs
{
  glm::dvec3 u {};
  glm::dvec3 v {};
  glm::dvec3 w {};
};

/// Unwrap polar angles around the first corner, then scale circumference and height.
inline AdjustedBoundaryTriangleUvs adjustedBoundaryTriangleUvs(const glm::dvec2& raw_u, const glm::dvec2& raw_v,
  const glm::dvec2& raw_w, double uv_circum_factor, double uv_height_factor)
{
  AdjustedBoundaryTriangleUvs adjusted;
  adjusted.u = glm::dvec3(raw_u.x, raw_u.y, 0.0);
  adjusted.v = glm::dvec3(raw_v.x, raw_v.y, 0.0);
  adjusted.w = glm::dvec3(raw_w.x, raw_w.y, 0.0);

  const double base_angle = adjusted.u.x;
  double& angle_v = adjusted.v.x;
  angle_v -= std::round(angle_v - base_angle);
  double& angle_w = adjusted.w.x;
  angle_w -= std::round(angle_w - base_angle);

  adjusted.u.x *= uv_circum_factor;
  adjusted.v.x *= uv_circum_factor;
  adjusted.w.x *= uv_circum_factor;
  adjusted.u.y *= uv_height_factor;
  adjusted.v.y *= uv_height_factor;
  adjusted.w.y *= uv_height_factor;
  return adjusted;
}

/// The single bark/interior shell for a tree.
/// Section events on different runtime branches all append to this mesh, so every public method takes
/// @ref mutex_. The mutex is recursive: @ref withLock and methods that call each other (section completion
/// emits triangles, which appends UVs and faces) must not deadlock the thread that already holds it.
class BoundaryMesh
{
 public:
  explicit BoundaryMesh(std::vector<std::string> material_names = {});

  void setStoreMetadata(bool store);

  void resetHalfEdgeSlots(size_t slot_count);
  void ensureHalfEdgeSlots(size_t slot_count);
  size_t halfEdgeSlotCount() const;
  void clearHalfEdge(size_t he_id);

  std::pair<size_t, size_t> lastLeftAndRight(size_t he_id) const;
  void setLastLeftAndRight(size_t he_id, size_t left, size_t right);
  void setHalfEdgeBoundaryVertex(size_t he_id, size_t vertex_index);

  size_t addVertex(const glm::dvec3& position, size_t strand_id, const glm::dvec2& raw_uv, const std::string& metadata);

  size_t addTriangle(size_t u, size_t v, size_t w, const std::string& metadata, int material_id,
    double uv_circum_factor, double uv_height_factor);

  /// Emit the quads/triangles that close the previous boundary ring on @p he_id, then store @p new_left / @p new_right
  /// as that half-edge's latest corners.
  void completeSection(size_t he_id, size_t new_left, size_t new_right, const std::string& face_metadata,
    double uv_circum_factor, double uv_height_factor);

  size_t triangleCount() const;

  void applyTransform(const glm::dmat4& transform);

  /// Merge duplicate and isolated vertices, drop degenerate faces, then normals and UV layout.
  /// Remaps @ref vertexToStrandId through the same vertex map. Returns how many degenerate triangles were removed.
  size_t finalize();

  VoronoiMesh& mesh();
  const VoronoiMesh& mesh() const;
  const std::vector<size_t>& vertexToStrandId() const;

  /// Hold @ref mutex_ for the whole functor. Nested @ref BoundaryMesh calls from @p fn lock again on this thread.
  template <typename Fn>
  decltype(auto) withLock(Fn&& fn)
  {
    std::lock_guard<std::recursive_mutex> guard(mutex_);
    return std::forward<Fn>(fn)();
  }

 private:
  static constexpr size_t npos = static_cast<size_t>(-1);

  size_t addVertexUnlocked(
    const glm::dvec3& position, size_t strand_id, const glm::dvec2& raw_uv, const std::string& metadata);
  size_t addTriangleUnlocked(size_t u, size_t v, size_t w, const std::string& metadata, int material_id,
    double uv_circum_factor, double uv_height_factor);
  void completeSectionUnlocked(size_t he_id, size_t new_left, size_t new_right, const std::string& face_metadata,
    double uv_circum_factor, double uv_height_factor);

  mutable std::recursive_mutex mutex_;
  VoronoiMesh mesh_;
  std::vector<std::pair<size_t, size_t>> last_left_and_right_;
  std::vector<size_t> vertex_to_strand_id_;
  std::vector<glm::dvec2> raw_uvs_;
  /// Vertex created when a flip turned this half-edge into a boundary edge; -1 when the previous ring is still valid.
  std::vector<int> half_edge_to_boundary_vertex_index_;
};

} // namespace kinDS

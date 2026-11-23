#include "BoundaryMesh.hpp"

#include "Logger.hpp"

#include <cassert>
#include <utility>

namespace kinDS
{

BoundaryMesh::BoundaryMesh(std::vector<std::string> material_names)
  : mesh_(std::move(material_names))
{
}

void BoundaryMesh::setStoreMetadata(bool store)
{
  std::lock_guard<std::recursive_mutex> guard(mutex_);
  mesh_.setStoreMetadata(store);
}

void BoundaryMesh::resetHalfEdgeSlots(size_t slot_count)
{
  std::lock_guard<std::recursive_mutex> guard(mutex_);
  last_left_and_right_.assign(slot_count, std::make_pair(npos, npos));
  half_edge_to_boundary_vertex_index_.assign(slot_count, -1);
}

void BoundaryMesh::ensureHalfEdgeSlots(size_t slot_count)
{
  std::lock_guard<std::recursive_mutex> guard(mutex_);
  if (last_left_and_right_.size() < slot_count)
  {
    last_left_and_right_.resize(slot_count, std::make_pair(npos, npos));
  }
  if (half_edge_to_boundary_vertex_index_.size() < slot_count)
  {
    half_edge_to_boundary_vertex_index_.resize(slot_count, -1);
  }
}

size_t BoundaryMesh::halfEdgeSlotCount() const
{
  std::lock_guard<std::recursive_mutex> guard(mutex_);
  return last_left_and_right_.size();
}

std::pair<size_t, size_t> BoundaryMesh::lastLeftAndRight(size_t he_id) const
{
  std::lock_guard<std::recursive_mutex> guard(mutex_);
  return last_left_and_right_[he_id];
}

void BoundaryMesh::setLastLeftAndRight(size_t he_id, size_t left, size_t right)
{
  std::lock_guard<std::recursive_mutex> guard(mutex_);
  last_left_and_right_[he_id] = std::make_pair(left, right);
  // A pending vertex belongs to the front it was created on; a new front must not inherit it.
  half_edge_to_boundary_vertex_index_[he_id] = -1;
}

std::pair<size_t, size_t> BoundaryMesh::closeFront(
  size_t he_id, size_t apex, const std::string& metadata, double uv_circum_factor, double uv_height_factor)
{
  std::lock_guard<std::recursive_mutex> guard(mutex_);
  const auto front = last_left_and_right_[he_id];
  const int pending = half_edge_to_boundary_vertex_index_[he_id];
  if (pending < 0)
  {
    addTriangleUnlocked(front.first, front.second, apex, metadata, 0, uv_circum_factor, uv_height_factor);
  }
  else
  {
    // The front still runs through the vertex a flip or radius event left on it.
    const size_t mid = static_cast<size_t>(pending);
    addTriangleUnlocked(front.first, mid, apex, metadata, 0, uv_circum_factor, uv_height_factor);
    addTriangleUnlocked(mid, front.second, apex, metadata, 0, uv_circum_factor, uv_height_factor);
    half_edge_to_boundary_vertex_index_[he_id] = -1;
  }
  return front;
}

std::optional<size_t> BoundaryMesh::vertexOfStrandAt(size_t strand_id, double kinetic_time) const
{
  std::lock_guard<std::recursive_mutex> guard(mutex_);
  const auto it = latest_vertex_of_strand_.find(strand_id);
  return it != latest_vertex_of_strand_.end() && vertex_kinetic_time_[it->second] == kinetic_time
    ? std::optional<size_t>(it->second)
    : std::nullopt;
}

void BoundaryMesh::bridgeCorners(
  size_t from, size_t to, size_t apex, const std::string& metadata, double uv_circum_factor, double uv_height_factor)
{
  std::lock_guard<std::recursive_mutex> guard(mutex_);
  addTriangleUnlocked(from, to, apex, metadata, 0, uv_circum_factor, uv_height_factor);
  corner_links_[from].push_back(to);
  corner_links_[to].push_back(from);
}

void BoundaryMesh::closeCornerGap(
  size_t from, size_t to, size_t apex, const std::string& metadata, double uv_circum_factor, double uv_height_factor)
{
  std::lock_guard<std::recursive_mutex> guard(mutex_);
  // Breadth-first search over bridged corners: events may already have linked the two corners through others.
  std::unordered_map<size_t, size_t> parent { { from, from } };
  std::vector<size_t> queue { from };
  for (size_t i = 0; i < queue.size() && !parent.contains(to); ++i)
  {
    if (const auto it = corner_links_.find(queue[i]); it != corner_links_.end())
    {
      for (const size_t next : it->second)
      {
        if (parent.emplace(next, queue[i]).second)
        {
          queue.push_back(next);
        }
      }
    }
  }
  std::vector<size_t> chain { to };
  if (parent.contains(to))
  {
    while (chain.back() != from)
    {
      chain.push_back(parent.at(chain.back()));
    }
  }
  else
  {
    chain.push_back(from);
  }
  for (size_t i = chain.size() - 1; i > 0; --i)
  {
    addTriangleUnlocked(chain[i], chain[i - 1], apex, metadata, 0, uv_circum_factor, uv_height_factor);
  }
  // A site on the boundary twice gets one gap per occurrence; the second one must fan through the first.
  corner_links_[from].push_back(to);
  corner_links_[to].push_back(from);
}

void BoundaryMesh::setHalfEdgeBoundaryVertex(size_t he_id, size_t vertex_index)
{
  std::lock_guard<std::recursive_mutex> guard(mutex_);
  half_edge_to_boundary_vertex_index_[he_id] = static_cast<int>(vertex_index);
}

void BoundaryMesh::clearHalfEdge(size_t he_id)
{
  std::lock_guard<std::recursive_mutex> guard(mutex_);
  if (he_id < last_left_and_right_.size())
  {
    last_left_and_right_[he_id] = std::make_pair(npos, npos);
  }
  if (he_id < half_edge_to_boundary_vertex_index_.size())
  {
    half_edge_to_boundary_vertex_index_[he_id] = -1;
  }
}

size_t BoundaryMesh::addVertex(const glm::dvec3& position, size_t strand_id, double kinetic_time,
  const glm::dvec2& raw_uv, const std::string& metadata)
{
  std::lock_guard<std::recursive_mutex> guard(mutex_);
  return addVertexUnlocked(position, strand_id, kinetic_time, raw_uv, metadata);
}

size_t BoundaryMesh::addVertexUnlocked(const glm::dvec3& position, size_t strand_id, double kinetic_time,
  const glm::dvec2& raw_uv, const std::string& metadata)
{
  const size_t index = mesh_.addVertex(position, metadata);
  vertex_to_strand_id_.push_back(strand_id);
  vertex_kinetic_time_.push_back(kinetic_time);
  latest_vertex_of_strand_[strand_id] = index;
  raw_uvs_.resize(index + 1, glm::dvec2 {});
  raw_uvs_[index] = raw_uv;
  return index;
}

size_t BoundaryMesh::addTriangle(size_t u, size_t v, size_t w, const std::string& metadata, int material_id,
  double uv_circum_factor, double uv_height_factor)
{
  std::lock_guard<std::recursive_mutex> guard(mutex_);
  return addTriangleUnlocked(u, v, w, metadata, material_id, uv_circum_factor, uv_height_factor);
}

size_t BoundaryMesh::addTriangleUnlocked(size_t u, size_t v, size_t w, const std::string& metadata, int material_id,
  double uv_circum_factor, double uv_height_factor)
{
  if (u >= mesh_.getVertexCount() || v >= mesh_.getVertexCount() || w >= mesh_.getVertexCount())
  {
    KINDS_ERROR("Vertex index " + std::to_string(u) + ", " + std::to_string(v) + ", " + std::to_string(w)
      + " out of boundary mesh range.");
    return npos;
  }

  if (u >= raw_uvs_.size() || v >= raw_uvs_.size() || w >= raw_uvs_.size())
  {
    KINDS_ERROR("Vertex index " + std::to_string(u) + ", " + std::to_string(v) + ", " + std::to_string(w)
      + " out of raw uv range.");
    return npos;
  }

  const AdjustedBoundaryTriangleUvs adjusted
    = adjustedBoundaryTriangleUvs(raw_uvs_[u], raw_uvs_[v], raw_uvs_[w], uv_circum_factor, uv_height_factor);

  const size_t uv_index_u = mesh_.addUV(adjusted.u);
  const size_t uv_index_v = mesh_.addUV(adjusted.v);
  const size_t uv_index_w = mesh_.addUV(adjusted.w);
  // Reversed: profile (x, y) lifts to (x, 0, y), which mirrors the kinetic 2D winding in object space.
  return mesh_.addTriangle(u, w, v, uv_index_u, uv_index_w, uv_index_v, material_id, metadata);
}

void BoundaryMesh::completeSection(size_t he_id, size_t new_left, size_t new_right, const std::string& face_metadata,
  double uv_circum_factor, double uv_height_factor)
{
  std::lock_guard<std::recursive_mutex> guard(mutex_);
  completeSectionUnlocked(he_id, new_left, new_right, face_metadata, uv_circum_factor, uv_height_factor);
}

void BoundaryMesh::completeSectionUnlocked(size_t he_id, size_t new_left, size_t new_right,
  const std::string& face_metadata, double uv_circum_factor, double uv_height_factor)
{
  const auto& last_left_and_right = last_left_and_right_[he_id];
  if (last_left_and_right.first != npos)
  {
    if (half_edge_to_boundary_vertex_index_[he_id] == -1)
    {
      addTriangleUnlocked(
        last_left_and_right.first, new_right, new_left, face_metadata, 0, uv_circum_factor, uv_height_factor);
      if (last_left_and_right.second != npos)
      {
        addTriangleUnlocked(last_left_and_right.second, new_right, last_left_and_right.first, face_metadata, 0,
          uv_circum_factor, uv_height_factor);
      }
    }
    else
    {
      assert(last_left_and_right.second != npos);
      const size_t flipped_vertex = static_cast<size_t>(half_edge_to_boundary_vertex_index_[he_id]);
      addTriangleUnlocked(
        last_left_and_right.second, new_right, flipped_vertex, face_metadata, 0, uv_circum_factor, uv_height_factor);
      addTriangleUnlocked(
        new_left, last_left_and_right.first, flipped_vertex, face_metadata, 0, uv_circum_factor, uv_height_factor);
      addTriangleUnlocked(new_left, flipped_vertex, new_right, face_metadata, 0, uv_circum_factor, uv_height_factor);
      half_edge_to_boundary_vertex_index_[he_id] = -1;
    }
  }

  last_left_and_right_[he_id] = std::make_pair(new_left, new_right);
}

size_t BoundaryMesh::triangleCount() const
{
  std::lock_guard<std::recursive_mutex> guard(mutex_);
  return mesh_.getTriangleCount();
}

void BoundaryMesh::applyTransform(const glm::dmat4& transform)
{
  std::lock_guard<std::recursive_mutex> guard(mutex_);
  mesh_.applyTransform(transform);
}

size_t BoundaryMesh::finalize()
{
  std::lock_guard<std::recursive_mutex> guard(mutex_);
  // A site touching the boundary twice gets two coincident vertices. While the pinch lasts they share no edge and must
  // stay apart (a position merge would join the two sheets along a non-manifold edge); once events reconnect them, the
  // joining faces have a zero-length edge, which is collapsed here.
  std::vector<size_t> root(mesh_.getVertexCount());
  for (size_t i = 0; i < root.size(); ++i)
  {
    root[i] = i;
  }
  const auto find = [&](size_t v)
  {
    while (root[v] != v)
    {
      v = root[v] = root[root[v]];
    }
    return v;
  };
  auto& triangles = mesh_.getTriangles();
  const auto& vertices = mesh_.getVertices();
  for (size_t i = 0; i < triangles.size(); ++i)
  {
    const size_t a = triangles[i];
    const size_t b = triangles[i % 3 == 2 ? i - 2 : i + 1];
    if (a != b && vertices[a] == vertices[b])
    {
      root[find(a)] = find(b);
    }
  }
  for (size_t& v : triangles)
  {
    v = find(v);
  }
  const size_t degenerate_triangles_removed = mesh_.removeDegenerateTriangles();
  const std::vector<size_t> remap_isolated = mesh_.removeIsolatedVertices();
  mesh_.ensureFaceMetadataSize();
  mesh_.computeNormals(NormalMode::PerTriangleCorner);
  mesh_.validateUVLayout("finalized boundary mesh");

  std::vector<size_t> remapped_strand_ids(mesh_.getVertices().size());
  std::vector<double> remapped_kinetic_times(mesh_.getVertices().size());
  for (size_t old_index = 0; old_index < vertex_to_strand_id_.size(); ++old_index)
  {
    const size_t new_index = remap_isolated[old_index];
    if (new_index != npos)
    {
      remapped_strand_ids[new_index] = vertex_to_strand_id_[old_index];
      remapped_kinetic_times[new_index] = vertex_kinetic_time_[old_index];
    }
  }
  vertex_to_strand_id_.swap(remapped_strand_ids);
  vertex_kinetic_time_.swap(remapped_kinetic_times);
  return degenerate_triangles_removed;
}

VoronoiMesh& BoundaryMesh::mesh()
{
  return mesh_;
}

const VoronoiMesh& BoundaryMesh::mesh() const
{
  return mesh_;
}

const std::vector<size_t>& BoundaryMesh::vertexToStrandId() const
{
  return vertex_to_strand_id_;
}

const std::vector<double>& BoundaryMesh::vertexKineticTime() const { return vertex_kinetic_time_; }

} // namespace kinDS

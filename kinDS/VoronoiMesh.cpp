#include "VoronoiMesh.hpp"
#include "DebugExportFormatting.hpp"
#include "Logger.hpp"
#include "VoronoiMesh.hpp"
#include "glm/gtx/norm.hpp"
#include <glm/gtc/matrix_inverse.hpp>
#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <iomanip>
#include <limits>
#include <optional>
#include <stdexcept>
#include <sstream>
#include <string>
#include <unordered_map>
#ifdef USE_CGAL
#include <CGAL/Polygon_mesh_processing/triangulate_hole.h>
#include <CGAL/Surface_mesh.h>

#include <vector>

namespace PMP = CGAL::Polygon_mesh_processing;

typedef CGAL::Simple_cartesian<double> Kernel;
typedef Kernel::Point_3 Point_3;
typedef CGAL::Surface_mesh<Point_3> Surface_mesh;
#endif

using namespace kinDS;

namespace
{
std::string normalModeToString(NormalMode normal_mode)
{
  switch (normal_mode)
  {
  case NormalMode::NoNormals:
    return "NoNormals";
  case NormalMode::PerVertex:
    return "PerVertex";
  case NormalMode::PerTriangleCorner:
    return "PerTriangleCorner";
    default:
    return "Unknown";
  }
}

bool metadataLooksEmpty(const std::string& metadata)
{
  if (metadata.empty())
  {
    return true;
  }
  size_t i = 0;
  while (i < metadata.size() && std::isspace(static_cast<unsigned char>(metadata[i])))
  {
    ++i;
  }
  if (i >= metadata.size() || metadata[i] != '{')
  {
    return false;
  }
  ++i;
  while (i < metadata.size() && std::isspace(static_cast<unsigned char>(metadata[i])))
  {
    ++i;
  }
  return i < metadata.size() && metadata[i] == '}';
}

std::optional<std::string> metadataEventType(const std::string& metadata)
{
  const std::string needle = "\"event_type\":";
  const size_t key_pos = metadata.find(needle);
  if (key_pos == std::string::npos)
  {
    return std::nullopt;
  }
  size_t value_pos = key_pos + needle.size();
  while (value_pos < metadata.size() && std::isspace(static_cast<unsigned char>(metadata[value_pos])))
  {
    ++value_pos;
  }
  if (value_pos >= metadata.size() || metadata[value_pos] != '"')
  {
    return std::nullopt;
  }
  ++value_pos;
  std::string value;
  for (; value_pos < metadata.size(); ++value_pos)
  {
    const char ch = metadata[value_pos];
    if (ch == '"')
    {
      return value;
    }
    if (ch == '\\' && value_pos + 1 < metadata.size())
    {
      ++value_pos;
      value += metadata[value_pos];
      continue;
    }
    value += ch;
  }
  return std::nullopt;
}

/// Higher wins when coincident vertices are merged. @c glue_align is below every other event type.
int vertexMetadataMergePriority(const std::string& metadata)
{
  if (metadataLooksEmpty(metadata))
  {
    return 0;
  }
  const std::optional<std::string> event_type = metadataEventType(metadata);
  if (event_type.has_value() && event_type.value() == "glue_align")
  {
    return 1;
  }
  return 2;
}

std::string preferredMergedVertexMetadata(const std::string& kept, const std::string& incoming)
{
  if (vertexMetadataMergePriority(incoming) > vertexMetadataMergePriority(kept))
  {
    return incoming;
  }
  return kept;
}
} // namespace

std::string VoronoiMesh::creationKineticTimeFilenameSuffix() const
{
  if (!std::isfinite(creation_kinetic_time_))
  {
    return {};
  }
  return "_" + formatDebugExportTimeToken(creation_kinetic_time_);
}

std::array<double, 3> barycentricCoordinates(
  const glm::dvec3& A, const glm::dvec3& B, const glm::dvec3& C, const glm::dvec3& P)
{
  // Vectors
  const glm::dvec3 v0 = B - A;
  const glm::dvec3 v1 = C - A;
  const glm::dvec3 v2 = P - A;

  // Dot products
  const double d00 = glm::dot(v0, v0);
  const double d01 = glm::dot(v0, v1);
  const double d11 = glm::dot(v1, v1);
  const double d20 = glm::dot(v2, v0);
  const double d21 = glm::dot(v2, v1);

  // Compute barycentric coordinates
  const double denom = d00 * d11 - d01 * d01;

  // Degenerate triangle check
  if (std::abs(denom) < 1e-15)
  {
    // Fallback: put everything on A
    return { 1.0, 0.0, 0.0 };
  }

  const double v = (d11 * d20 - d01 * d21) / denom;
  const double w = (d00 * d21 - d01 * d20) / denom;
  const double u = 1.0 - v - w;

  return { u, v, w };
}

size_t VoronoiMesh::addVertex(double x, double y, double z, const std::string& metadata, const glm::dvec3& color)
{
  size_t index = vertices.size();
  vertices.emplace_back(glm::dvec3 { x, y, z });
  vertex_kinetic_times_.push_back(std::numeric_limits<double>::quiet_NaN());
  vertex_is_flexible_.push_back(false);
  if (store_metadata_)
  {
    vertex_metadata.push_back(metadata);
  }
  vertex_colors.push_back(color);
  return index;
}

size_t VoronoiMesh::addVertex(const glm::dvec3& p, const std::string& metadata, const glm::dvec3& color)
{
  size_t index = vertices.size();
  vertices.emplace_back(p);
  vertex_kinetic_times_.push_back(std::numeric_limits<double>::quiet_NaN());
  vertex_is_flexible_.push_back(false);
  if (store_metadata_)
  {
    vertex_metadata.push_back(metadata);
  }
  vertex_colors.push_back(color);
  return index;
}

void VoronoiMesh::setProfilePlanePosition(size_t vertex_index, glm::dvec2 xy)
{
  if (vertex_index >= vertices.size())
  {
    throw std::out_of_range("setProfilePlanePosition: vertex index out of range.");
  }
  if (profile_plane_xy_.size() < vertices.size())
  {
    // Unset slots must not look like valid (0,0) triangulation samples.
    profile_plane_xy_.resize(vertices.size(),
      glm::dvec2(std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN()));
  }
  profile_plane_xy_[vertex_index] = xy;
}

glm::dvec2 VoronoiMesh::triangulationPlaneXY(size_t vertex_index) const
{
  if (vertex_index >= vertices.size())
  {
    throw std::out_of_range("triangulationPlaneXY: vertex index out of range.");
  }
  if (vertex_index < profile_plane_xy_.size())
  {
    const glm::dvec2& xy = profile_plane_xy_[vertex_index];
    if (std::isfinite(xy.x) && std::isfinite(xy.y))
    {
      return xy;
    }
  }
  const glm::dvec3& v = vertices[vertex_index];
  return glm::dvec2(v.x, v.y);
}

void VoronoiMesh::setVertexKineticTime(size_t vertex_index, double t)
{
  if (vertex_index >= vertices.size())
  {
    throw std::out_of_range("setVertexKineticTime: vertex index out of range.");
  }
  if (vertex_kinetic_times_.size() < vertices.size())
  {
    vertex_kinetic_times_.resize(vertices.size(), std::numeric_limits<double>::quiet_NaN());
  }
  vertex_kinetic_times_[vertex_index] = t;
}

double VoronoiMesh::vertexKineticTime(size_t vertex_index) const
{
  if (vertex_index >= vertices.size())
  {
    throw std::out_of_range("vertexKineticTime: vertex index out of range.");
  }
  if (vertex_index < vertex_kinetic_times_.size() && std::isfinite(vertex_kinetic_times_[vertex_index]))
  {
    return vertex_kinetic_times_[vertex_index];
  }
  return vertices[vertex_index].z;
}

void VoronoiMesh::setVertexSemanticUv(size_t vertex_index, const glm::dvec3& uv)
{
  if (vertex_index >= vertices.size())
  {
    throw std::out_of_range("setVertexSemanticUv: vertex index out of range.");
  }
  if (vertex_semantic_uvs_.size() < vertices.size())
  {
    vertex_semantic_uvs_.resize(vertices.size(), glm::dvec3(std::numeric_limits<double>::quiet_NaN()));
  }
  vertex_semantic_uvs_[vertex_index] = uv;
}

std::optional<glm::dvec3> VoronoiMesh::vertexSemanticUv(size_t vertex_index) const
{
  if (vertex_index >= vertex_semantic_uvs_.size())
  {
    return std::nullopt;
  }
  const glm::dvec3& uv = vertex_semantic_uvs_[vertex_index];
  if (!std::isfinite(uv.x) || !std::isfinite(uv.y) || !std::isfinite(uv.z))
  {
    return std::nullopt;
  }
  return uv;
}

void VoronoiMesh::setVertexMetadata(size_t vertex_index, const std::string& metadata)
{
  if (vertex_index >= vertices.size())
  {
    throw std::out_of_range("setVertexMetadata: vertex index out of range.");
  }
  if (!store_metadata_)
  {
    return;
  }
  if (vertex_metadata.size() < vertices.size())
  {
    vertex_metadata.resize(vertices.size(), "{}");
  }
  vertex_metadata[vertex_index] = metadata;
}

void VoronoiMesh::setVertexFlexible(size_t vertex_index, bool is_flexible)
{
  if (vertex_index >= vertices.size())
  {
    throw std::out_of_range("setVertexFlexible: vertex index out of range.");
  }
  if (vertex_is_flexible_.size() < vertices.size())
  {
    vertex_is_flexible_.resize(vertices.size(), false);
  }
  vertex_is_flexible_[vertex_index] = is_flexible;
}

bool VoronoiMesh::isVertexFlexible(size_t vertex_index) const
{
  if (vertex_index >= vertices.size())
  {
    throw std::out_of_range("isVertexFlexible: vertex index out of range.");
  }
  return vertex_index < vertex_is_flexible_.size() && vertex_is_flexible_[vertex_index];
}

void kinDS::VoronoiMesh::replaceVertex(size_t index, const glm::dvec3& new_position) { vertices[index] = new_position; }

size_t VoronoiMesh::addTriangle(size_t v1, size_t v2, size_t v3, int material_id, const std::string& metadata)
{
  return addTriangle(v1, v2, v3, std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max(),
    std::numeric_limits<size_t>::max(), material_id, metadata);
}

size_t VoronoiMesh::addTriangle(
  size_t v1, size_t v2, size_t v3, size_t uv1, size_t uv2, size_t uv3, int material_id, const std::string& metadata)
{
  size_t index = triangles.size() / 3;

  // Check vertex indices
  if (v1 >= vertices.size() || v2 >= vertices.size() || v3 >= vertices.size())
  {
    throw std::out_of_range("Vertex index out of range when adding triangle.");
  }

  triangles.push_back(v1);
  triangles.push_back(v2);
  triangles.push_back(v3);
  uv_indices.push_back(uv1);
  uv_indices.push_back(uv2);
  uv_indices.push_back(uv3);

  material_ids.push_back(material_id);
  if (store_metadata_)
  {
    face_metadata.push_back(metadata);
  }

  return index;
}

size_t VoronoiMesh::triangleCornerIndex(size_t triangle_id, size_t vertex_id) const
{
  if (triangle_id >= triangles.size() / 3)
  {
    return static_cast<size_t>(-1);
  }
  for (size_t e = 0; e < 3; ++e)
  {
    const size_t corner = 3 * triangle_id + e;
    if (triangles[corner] == vertex_id)
    {
      return corner;
    }
  }
  return static_cast<size_t>(-1);
}

std::pair<size_t, size_t> VoronoiMesh::splitTriangle(size_t tri_vert_id0, size_t tri_vert_id1, const glm::dvec3& vertex,
  const std::string& metadata, std::optional<double> kinetic_time)
{
  constexpr std::pair<size_t, size_t> k_fail { static_cast<size_t>(-1), static_cast<size_t>(-1) };
  if (tri_vert_id0 >= triangles.size() || tri_vert_id1 >= triangles.size())
  {
    return k_fail;
  }
  // Same triangle iff integer division by 3 matches (drop remainder intentionally).
  if ((tri_vert_id0 / 3) != (tri_vert_id1 / 3) || tri_vert_id0 == tri_vert_id1)
  {
    return k_fail;
  }

  const size_t tri = tri_vert_id0 / 3;
  const size_t base = 3 * tri;

  // Orient the split edge to the triangle's existing winding: corner_a → corner_b along the cycle.
  // Callers may pass the two corners in either order.
  size_t corner_a = tri_vert_id0;
  size_t corner_b = tri_vert_id1;
  if (((corner_a % 3) + 1) % 3 != (corner_b % 3))
  {
    std::swap(corner_a, corner_b);
  }
  if (((corner_a % 3) + 1) % 3 != (corner_b % 3))
  {
    return k_fail; // not an edge of this triangle (should be unreachable for distinct corners)
  }
  const size_t corner_c = base + (((corner_b % 3) + 1) % 3);

  const size_t va = triangles[corner_a];
  const size_t vb = triangles[corner_b];
  const size_t vc = triangles[corner_c];
  if (va >= vertices.size() || vb >= vertices.size() || vc >= vertices.size())
  {
    return k_fail;
  }

  // Parameter along the winding-oriented edge va → vb.
  const glm::dvec3& pa = vertices[va];
  const glm::dvec3& pb = vertices[vb];
  const glm::dvec3 ab = pb - pa;
  const double len2 = glm::dot(ab, ab);
  double s = 0.5;
  if (len2 > 1e-24)
  {
    s = std::clamp(glm::dot(vertex - pa, ab) / len2, 0.0, 1.0);
  }

  const size_t new_vid = addVertex(vertex, metadata);

  // Interpolate per-vertex attributes from the split edge endpoints (va → vb).
  if (kinetic_time.has_value() && std::isfinite(kinetic_time.value()))
  {
    setVertexKineticTime(new_vid, kinetic_time.value());
  }
  else
  {
    const double t0 = vertexKineticTime(va);
    const double t1 = vertexKineticTime(vb);
    if (std::isfinite(t0) && std::isfinite(t1))
    {
      setVertexKineticTime(new_vid, t0 * (1.0 - s) + t1 * s);
    }
  }
  if (const auto uv0 = vertexSemanticUv(va); uv0.has_value())
  {
    if (const auto uv1 = vertexSemanticUv(vb); uv1.has_value())
    {
      setVertexSemanticUv(new_vid, uv0.value() * (1.0 - s) + uv1.value() * s);
    }
  }
  {
    const bool p0_ok = va < profile_plane_xy_.size() && std::isfinite(profile_plane_xy_[va].x)
      && std::isfinite(profile_plane_xy_[va].y);
    const bool p1_ok = vb < profile_plane_xy_.size() && std::isfinite(profile_plane_xy_[vb].x)
      && std::isfinite(profile_plane_xy_[vb].y);
    if (p0_ok && p1_ok)
    {
      setProfilePlanePosition(new_vid, profile_plane_xy_[va] * (1.0 - s) + profile_plane_xy_[vb] * s);
    }
  }
  if (isVertexFlexible(va) || isVertexFlexible(vb))
  {
    setVertexFlexible(new_vid, true);
  }
  if (va < vertex_colors.size() && vb < vertex_colors.size() && new_vid < vertex_colors.size())
  {
    vertex_colors[new_vid] = vertex_colors[va] * (1.0 - s) + vertex_colors[vb] * s;
  }
  // After addVertex, a synced PerVertex normal buffer has size == vertices.size() - 1.
  if (normal_mode == PerVertex && normals.size() + 1 == vertices.size() && va < normals.size() && vb < normals.size())
  {
    normals.push_back(normals[va] * (1.0 - s) + normals[vb] * s);
  }

  size_t uv_a_idx = std::numeric_limits<size_t>::max();
  size_t uv_b_idx = std::numeric_limits<size_t>::max();
  size_t uv_c_idx = std::numeric_limits<size_t>::max();
  size_t uv_f_idx = std::numeric_limits<size_t>::max();
  const bool has_uvs = uv_indices.size() == triangles.size();
  if (has_uvs)
  {
    uv_a_idx = uv_indices[corner_a];
    uv_b_idx = uv_indices[corner_b];
    uv_c_idx = uv_indices[corner_c];
    if (uv_a_idx < uvs.size() && uv_b_idx < uvs.size())
    {
      uv_f_idx = addUV(uvs[uv_a_idx] * (1.0 - s) + uvs[uv_b_idx] * s);
    }
  }

  glm::dvec3 n_a { 0.0 };
  glm::dvec3 n_b { 0.0 };
  glm::dvec3 n_c { 0.0 };
  glm::dvec3 n_f { 0.0 };
  const bool has_corner_normals = normal_mode == PerTriangleCorner && normals.size() == triangles.size();
  if (has_corner_normals)
  {
    n_a = normals[corner_a];
    n_b = normals[corner_b];
    n_c = normals[corner_c];
    n_f = n_a * (1.0 - s) + n_b * s;
  }

  // Keep per-face attribute buffers aligned with the existing triangle count before appending.
  const size_t tri_count_before = triangles.size() / 3;
  if (material_ids.size() < tri_count_before)
  {
    material_ids.resize(tri_count_before, -1);
  }
  if (store_metadata_ && face_metadata.size() < tri_count_before)
  {
    face_metadata.resize(tri_count_before, "{}");
  }

  const int material_id = material_ids[tri];
  std::string face_meta = "{}";
  if (store_metadata_)
  {
    face_meta = face_metadata[tri];
  }

  // Original face (…, va, vb, vc, …) → (…, va, F, vc, …): replace vb with F, preserve winding.
  triangles[corner_b] = new_vid;
  if (has_uvs)
  {
    uv_indices[corner_b] = uv_f_idx;
  }
  if (has_corner_normals)
  {
    normals[corner_b] = n_f;
  }

  // New face: F, vb, vc — same cyclic order as va → vb → vc on the original triangle.
  const size_t new_tri = addTriangle(new_vid, vb, vc, uv_f_idx, uv_b_idx, uv_c_idx, material_id, face_meta);
  if (has_corner_normals)
  {
    normals.push_back(n_f);
    normals.push_back(n_b);
    normals.push_back(n_c);
  }

  return { new_vid, new_tri };
}

size_t VoronoiMesh::addNormal(double nx, double ny, double nz) { return addNormal(glm::dvec3 { nx, ny, nz }); }

size_t VoronoiMesh::addNormal(const glm::dvec3& n)
{
  if (normal_mode == NormalMode::NoNormals)
  {
    throw std::runtime_error("Cannot add normals to a mesh with NormalMode::NoNormals.");
  }
  size_t index = normals.size();
  normals.emplace_back(n);
  return index;
}

void kinDS::VoronoiMesh::replaceNormal(size_t index, const glm::dvec3& new_normal) { normals[index] = new_normal; }

size_t VoronoiMesh::addUV(double u, double v, double w) { return addUV(glm::dvec3 { u, v, w }); }

size_t VoronoiMesh::addUV(glm::dvec3 uv)
{
  size_t index = uvs.size();
  uvs.emplace_back(uv);
  return index;
}

void VoronoiMesh::startNewGroup(const std::string& name)
{
  group_offsets.push_back(triangles.size() / 3);
  group_names.push_back(name);
}

void VoronoiMesh::setGroupOffsets(const std::vector<size_t>& offsets)
{
  group_offsets = offsets;
  if (group_names.size() > group_offsets.size())
  {
    group_names.resize(group_offsets.size());
  }
  else if (group_names.size() < group_offsets.size())
  {
    group_names.resize(group_offsets.size());
  }
}

void VoronoiMesh::setGroupNames(const std::vector<std::string>& names)
{
  group_names = names;
  if (group_names.size() < group_offsets.size())
  {
    group_names.resize(group_offsets.size());
  }
}

int VoronoiMesh::ensureMaterialName(const std::string& name)
{
  const auto it = std::find(material_names.begin(), material_names.end(), name);
  if (it != material_names.end())
  {
    return static_cast<int>(std::distance(material_names.begin(), it));
  }
  material_names.push_back(name);
  return static_cast<int>(material_names.size() - 1);
}

void VoronoiMesh::ensureMaterialIdsSize(int fill_material_id)
{
  const size_t tri_count = getTriangleCount();
  if (material_ids.size() < tri_count)
  {
    material_ids.resize(tri_count, fill_material_id);
  }
}

void VoronoiMesh::setTriangleMaterialId(size_t triangle_index, int material_id)
{
  ensureMaterialIdsSize(-1);
  if (triangle_index < material_ids.size())
  {
    material_ids[triangle_index] = material_id;
  }
}

void VoronoiMesh::setVertexColor(size_t vertex_index, const glm::dvec3& color)
{
  if (vertex_index >= vertices.size())
  {
    return;
  }
  if (vertex_colors.size() < vertices.size())
  {
    vertex_colors.resize(vertices.size(), glm::dvec3(1.0));
  }
  vertex_colors[vertex_index] = color;
}

VoronoiMesh& VoronoiMesh::operator+=(const VoronoiMesh& other)
{
  if (other.vertices.empty() && other.triangles.empty() && other.normals.empty())
  {
    return *this;
  }
  if (normal_mode != other.normal_mode)
  {
    // Maybe we should implement auto-conversion at some point, but for now just throw an error
    throw std::runtime_error("Meshes don't use the same normal mode: left=" + normalModeToString(normal_mode)
      + ", right=" + normalModeToString(other.normal_mode) + ".");
  }
  validateNormalCount("left mesh before combine");
  other.validateNormalCount("right mesh before combine");

  const size_t left_triangle_count = triangles.size() / 3;
  if (material_ids.size() < left_triangle_count)
  {
    material_ids.resize(left_triangle_count, -1);
  }
  if (triangles.empty() && material_names.empty() && !other.material_names.empty())
  {
    material_names = other.material_names;
  }

  size_t old_vertices_size = vertices.size();
  vertices.insert(vertices.end(), other.vertices.begin(), other.vertices.end());
  if (store_metadata_)
  {
    if (other.vertex_metadata.size() == other.vertices.size())
    {
      vertex_metadata.insert(vertex_metadata.end(), other.vertex_metadata.begin(), other.vertex_metadata.end());
    }
    else
    {
      vertex_metadata.insert(vertex_metadata.end(), other.vertices.size(), "{}");
    }
  }
  if (other.vertex_colors.size() == other.vertices.size())
  {
    vertex_colors.insert(vertex_colors.end(), other.vertex_colors.begin(), other.vertex_colors.end());
  }
  else
  {
    vertex_colors.insert(vertex_colors.end(), other.vertices.size(), glm::dvec3(1.0));
  }
  if (other.vertex_kinetic_times_.size() == other.vertices.size())
  {
    vertex_kinetic_times_.insert(
      vertex_kinetic_times_.end(), other.vertex_kinetic_times_.begin(), other.vertex_kinetic_times_.end());
  }
  else
  {
    vertex_kinetic_times_.insert(
      vertex_kinetic_times_.end(), other.vertices.size(), std::numeric_limits<double>::quiet_NaN());
  }
  if (other.vertex_semantic_uvs_.size() == other.vertices.size())
  {
    vertex_semantic_uvs_.insert(vertex_semantic_uvs_.end(), other.vertex_semantic_uvs_.begin(), other.vertex_semantic_uvs_.end());
  }
  else
  {
    vertex_semantic_uvs_.insert(
      vertex_semantic_uvs_.end(), other.vertices.size(), glm::dvec3(std::numeric_limits<double>::quiet_NaN()));
  }
  if (vertex_is_flexible_.size() < old_vertices_size)
  {
    vertex_is_flexible_.resize(old_vertices_size, false);
  }
  if (other.vertex_is_flexible_.size() == other.vertices.size())
  {
    vertex_is_flexible_.insert(
      vertex_is_flexible_.end(), other.vertex_is_flexible_.begin(), other.vertex_is_flexible_.end());
  }
  else
  {
    vertex_is_flexible_.insert(vertex_is_flexible_.end(), other.vertices.size(), false);
  }

  size_t old_vertex_indices_size = triangles.size();
  const size_t old_triangle_count = old_vertex_indices_size / 3;
  triangles.insert(triangles.end(), other.triangles.begin(), other.triangles.end());

  std::transform(triangles.begin() + old_vertex_indices_size, triangles.end(),
    triangles.begin() + old_vertex_indices_size, [&](size_t index) { return index + old_vertices_size; });

  normals.insert(normals.end(), other.normals.begin(), other.normals.end());

  size_t old_uvs_size = uvs.size();
  uvs.insert(uvs.end(), other.uvs.begin(), other.uvs.end());

  size_t old_uv_indices_size = uv_indices.size();
  uv_indices.insert(uv_indices.end(), other.uv_indices.begin(), other.uv_indices.end());

  std::transform(uv_indices.begin() + old_uv_indices_size, uv_indices.end(), uv_indices.begin() + old_uv_indices_size,
    [&](size_t index)
    {
      return index == std::numeric_limits<size_t>::max() ? index : index + old_uvs_size;
    });

  size_t old_group_count = group_offsets.size();
  group_offsets.insert(group_offsets.end(), other.group_offsets.begin(), other.group_offsets.end());

  std::transform(group_offsets.begin() + old_group_count, group_offsets.end(), group_offsets.begin() + old_group_count,
    [&](size_t offset) { return offset + old_triangle_count; });

  group_names.insert(group_names.end(), other.group_names.begin(), other.group_names.end());
  if (group_names.size() < group_offsets.size())
  {
    group_names.resize(group_offsets.size());
  }
  else if (group_names.size() > group_offsets.size())
  {
    group_names.resize(group_offsets.size());
  }
  if (store_metadata_)
  {
    const size_t other_tri_count = other.triangles.size() / 3;
    if (other.store_metadata_ && other.face_metadata.size() == other_tri_count)
    {
      face_metadata.insert(face_metadata.end(), other.face_metadata.begin(), other.face_metadata.end());
    }
    else if (other.store_metadata_ && !other.face_metadata.empty())
    {
      face_metadata.reserve(face_metadata.size() + other_tri_count);
      for (size_t tri_index = 0; tri_index < other_tri_count; ++tri_index)
      {
        face_metadata.push_back(
          tri_index < other.face_metadata.size() ? other.face_metadata[tri_index] : std::string("{}"));
      }
    }
    else
    {
      face_metadata.insert(face_metadata.end(), other_tri_count, "{}");
    }
    ensureFaceMetadataSize();
  }

  std::vector<int> other_material_id_remap(other.material_names.size(), -1);
  for (size_t i = 0; i < other.material_names.size(); ++i)
  {
    const auto it = std::find(material_names.begin(), material_names.end(), other.material_names[i]);
    if (it == material_names.end())
    {
      other_material_id_remap[i] = static_cast<int>(material_names.size());
      material_names.push_back(other.material_names[i]);
    }
    else
    {
      other_material_id_remap[i] = static_cast<int>(std::distance(material_names.begin(), it));
    }
  }

  const size_t other_triangle_count = other.triangles.size() / 3;
  if (other.material_ids.size() == other_triangle_count)
  {
    material_ids.reserve(material_ids.size() + other_triangle_count);
    for (int id : other.material_ids)
    {
      if (id < 0 || id >= static_cast<int>(other_material_id_remap.size()))
      {
        material_ids.push_back(-1);
      }
      else
      {
        material_ids.push_back(other_material_id_remap[static_cast<size_t>(id)]);
      }
    }
  }
  else
  {
    material_ids.insert(material_ids.end(), other_triangle_count, -1);
  }

  validateNormalCount("combined mesh after combine");
  validateUVLayout("combined mesh after combine");
  return *this;
}

void VoronoiMesh::flipOrientation()
{
  // Flip the orientation of each triangle by swapping the second and third vertex indices
  for (size_t i = 0; i < triangles.size(); i += 3)
  {
    std::swap(triangles[i + 1], triangles[i + 2]);
    if (uv_indices.size() == triangles.size())
    {
      std::swap(uv_indices[i + 1], uv_indices[i + 2]);
    }
  }

  // flip all normals
  for (size_t i = 0; i < normals.size(); i++)
  {
    normals[i] = -normals[i];
  }
}

bool VoronoiMesh::orientFacesAwayFromCentroid()
{
  if (vertices.empty() || triangles.size() < 3)
  {
    return false;
  }

  glm::dvec3 centroid(0.0);
  for (const glm::dvec3& vertex : vertices)
  {
    centroid += vertex;
  }
  centroid /= static_cast<double>(vertices.size());

  // Positive score means normals predominantly point away from the centroid (desired outward polarity).
  double outward_score = 0.0;
  for (size_t i = 0; i < triangles.size(); i += 3)
  {
    const glm::dvec3& v0 = vertices[triangles[i]];
    const glm::dvec3& v1 = vertices[triangles[i + 1]];
    const glm::dvec3& v2 = vertices[triangles[i + 2]];
    const glm::dvec3 normal = glm::cross(v1 - v0, v2 - v0);
    const glm::dvec3 face_centroid = (v0 + v1 + v2) / 3.0;
    outward_score += glm::dot(normal, face_centroid - centroid);
  }

  if (outward_score >= 0.0)
  {
    return false;
  }

  flipOrientation();
  return true;
}

glm::dmat4 VoronoiMesh::profileSpaceSwapYAndZTransform()
{
  glm::dmat4 transform(1.0);
  transform[1][1] = 0.0;
  transform[2][1] = 1.0;
  transform[1][2] = 1.0;
  transform[2][2] = 0.0;
  return transform;
}

void VoronoiMesh::applyTransform(const glm::dmat4& transform)
{
  for (glm::dvec3& vertex : vertices)
  {
    vertex = glm::dvec3(transform * glm::dvec4(vertex, 1.0));
  }

  if (normal_mode != NoNormals && !normals.empty())
  {
    const glm::dmat3 linear = glm::dmat3(transform);
    const glm::dmat3 normal_matrix = glm::transpose(glm::inverse(linear));
    for (glm::dvec3& normal : normals)
    {
      normal = glm::normalize(normal_matrix * normal);
    }
  }

  if (glm::determinant(glm::dmat3(transform)) < 0.0)
  {
    for (size_t i = 0; i < triangles.size(); i += 3)
    {
      std::swap(triangles[i + 1], triangles[i + 2]);
      std::swap(uv_indices[i + 1], uv_indices[i + 2]);
    }
  }
}

std::vector<size_t> VoronoiMesh::mergeDuplicateVertices(double epsilon)
{
  const double inv_eps = (epsilon > 0.0) ? 1.0 / epsilon : 0.0;
  std::unordered_map<glm::dvec3, size_t, VoronoiMesh::Vec3iHash> grid;
  std::vector<glm::dvec3> newVerts;
  newVerts.reserve(vertices.size());
  std::vector<std::string> new_vertex_metadata;
  new_vertex_metadata.reserve(vertices.size());
  std::vector<glm::dvec3> new_vertex_colors;
  new_vertex_colors.reserve(vertices.size());
  std::vector<double> new_vertex_kinetic_times;
  new_vertex_kinetic_times.reserve(vertices.size());
  std::vector<glm::dvec3> new_vertex_semantic_uvs;
  new_vertex_semantic_uvs.reserve(vertices.size());
  std::vector<bool> new_vertex_is_flexible;
  new_vertex_is_flexible.reserve(vertices.size());

  std::vector<size_t> remap(vertices.size(), size_t(-1));

  for (size_t i = 0; i < vertices.size(); ++i)
  {
    const auto& v = vertices[i];

    // Quantize vertex for approximate matching
    glm::dvec3 key;
    if (epsilon > 0.0)
    {
      key[0] = static_cast<int>(std::llround(v[0] * inv_eps));
      key[1] = static_cast<int>(std::llround(v[1] * inv_eps));
      key[2] = static_cast<int>(std::llround(v[2] * inv_eps));
    }
    else
    {
      key[0] = static_cast<int>(std::hash<double> {}(v[0]) & 0x7FFFFFFF);
      key[1] = static_cast<int>(std::hash<double> {}(v[1]) & 0x7FFFFFFF);
      key[2] = static_cast<int>(std::hash<double> {}(v[2]) & 0x7FFFFFFF);
    }

    const bool is_flexible = i < vertex_is_flexible_.size() && vertex_is_flexible_[i];
    auto it = grid.find(key);
    if (it == grid.end())
    {
      size_t newIndex = newVerts.size();
      grid[key] = newIndex;
      newVerts.push_back(v);
      if (i < vertex_metadata.size())
      {
        new_vertex_metadata.push_back(vertex_metadata[i]);
      }
      else
      {
        new_vertex_metadata.push_back("{}");
      }
      if (i < vertex_colors.size())
      {
        new_vertex_colors.push_back(vertex_colors[i]);
      }
      else
      {
        new_vertex_colors.push_back(glm::dvec3(1.0));
      }
      if (i < vertex_kinetic_times_.size())
      {
        new_vertex_kinetic_times.push_back(vertex_kinetic_times_[i]);
      }
      else
      {
        new_vertex_kinetic_times.push_back(std::numeric_limits<double>::quiet_NaN());
      }
      if (i < vertex_semantic_uvs_.size())
      {
        new_vertex_semantic_uvs.push_back(vertex_semantic_uvs_[i]);
      }
      else
      {
        new_vertex_semantic_uvs.push_back(glm::dvec3(std::numeric_limits<double>::quiet_NaN()));
      }
      new_vertex_is_flexible.push_back(is_flexible);
      remap[i] = newIndex;
    }
    else
    {
      remap[i] = it->second;
      if (is_flexible)
      {
        new_vertex_is_flexible[it->second] = true;
      }
      if (store_metadata_ && it->second < new_vertex_metadata.size())
      {
        const std::string incoming = (i < vertex_metadata.size()) ? vertex_metadata[i] : std::string("{}");
        new_vertex_metadata[it->second]
          = preferredMergedVertexMetadata(new_vertex_metadata[it->second], incoming);
      }
    }
  }

  // Remap triangle vertex indices only. uv_indices are per triangle corner and must keep pointing at their own UV-pool
  // entries so merged shared vertices can retain distinct corner UVs (e.g. interior vs bark on the same position).
  for (size_t& idx : triangles)
  {
    idx = remap[idx];
  }

  vertices.swap(newVerts);
  vertex_metadata.swap(new_vertex_metadata);
  vertex_colors.swap(new_vertex_colors);
  vertex_kinetic_times_.swap(new_vertex_kinetic_times);
  vertex_semantic_uvs_.swap(new_vertex_semantic_uvs);
  vertex_is_flexible_.swap(new_vertex_is_flexible);
  if (!store_metadata_)
  {
    vertex_metadata.clear();
  }

  return remap;
}

size_t VoronoiMesh::collapseDegreeTwoFlexibleVertices()
{
  const size_t n_vertices = vertices.size();
  const size_t n_triangles = triangles.size() / 3;
  if (n_vertices == 0 || n_triangles == 0)
  {
    return 0;
  }

  struct MergedTriangle
  {
    size_t v0 = 0;
    size_t v1 = 0;
    size_t v2 = 0;
    size_t uv0 = std::numeric_limits<size_t>::max();
    size_t uv1 = std::numeric_limits<size_t>::max();
    size_t uv2 = std::numeric_limits<size_t>::max();
    int material_id = -1;
    std::string face_metadata = "{}";
  };

  const auto mark_collapsed_metadata = [](const std::string& metadata) -> std::string
  {
    if (metadata.empty() || metadata == "{}")
    {
      return "{\"postprocess_merged_flexible\":true}";
    }
    if (metadata.back() == '}')
    {
      std::string out = metadata;
      out.pop_back();
      if (out.size() > 1)
      {
        out += ",";
      }
      out += "\"postprocess_merged_flexible\":true}";
      return out;
    }
    return metadata;
  };

  std::vector<bool> triangle_removed(n_triangles, false);
  std::vector<MergedTriangle> merged;
  size_t collapsed = 0;

  const bool has_uvs = !uv_indices.empty() && uv_indices.size() == triangles.size();

  for (size_t vertex_index = 0; vertex_index < n_vertices; ++vertex_index)
  {
    if (!isVertexFlexible(vertex_index))
    {
      continue;
    }

    const std::vector<size_t> corners = findTriangleCorners(vertex_index);
    if (corners.size() != 2)
    {
      continue;
    }

    const size_t tri_a = corners[0] / 3;
    const size_t tri_b = corners[1] / 3;
    if (tri_a >= n_triangles || tri_b >= n_triangles || tri_a == tri_b)
    {
      continue;
    }
    if (triangle_removed[tri_a] || triangle_removed[tri_b])
    {
      continue;
    }

    const size_t local_a = corners[0] % 3;
    const size_t local_b = corners[1] % 3;

    // Opposite directed edge of tri_a (skip flexible vertex), preserves that triangle's winding.
    const size_t u = triangles[3 * tri_a + ((local_a + 1) % 3)];
    const size_t v = triangles[3 * tri_a + ((local_a + 2) % 3)];
    if (u == v || u == vertex_index || v == vertex_index)
    {
      continue;
    }

    size_t w = std::numeric_limits<size_t>::max();
    for (size_t k = 0; k < 3; ++k)
    {
      const size_t candidate = triangles[3 * tri_b + k];
      if (candidate != vertex_index && candidate != u && candidate != v)
      {
        w = candidate;
        break;
      }
    }
    if (w == std::numeric_limits<size_t>::max() || w == u || w == v)
    {
      continue;
    }

    // Confirm the two triangles share an edge through the flexible vertex.
    size_t shared_neighbors = 0;
    for (size_t k = 0; k < 3; ++k)
    {
      const size_t candidate = triangles[3 * tri_b + k];
      if (candidate == u || candidate == v)
      {
        ++shared_neighbors;
      }
    }
    if (shared_neighbors != 1)
    {
      continue;
    }

    MergedTriangle out;
    out.v0 = u;
    out.v1 = v;
    out.v2 = w;
    out.material_id = tri_a < material_ids.size() ? material_ids[tri_a] : -1;
    if (store_metadata_ && tri_a < face_metadata.size())
    {
      out.face_metadata = mark_collapsed_metadata(face_metadata[tri_a]);
    }
    else if (store_metadata_)
    {
      out.face_metadata = mark_collapsed_metadata("{}");
    }

    if (has_uvs)
    {
      out.uv0 = uv_indices[3 * tri_a + ((local_a + 1) % 3)];
      out.uv1 = uv_indices[3 * tri_a + ((local_a + 2) % 3)];
      size_t w_local = 0;
      for (; w_local < 3; ++w_local)
      {
        if (triangles[3 * tri_b + w_local] == w)
        {
          break;
        }
      }
      out.uv2 = uv_indices[3 * tri_b + w_local];
    }

    triangle_removed[tri_a] = true;
    triangle_removed[tri_b] = true;
    merged.push_back(std::move(out));
    ++collapsed;
  }

  if (collapsed == 0)
  {
    return 0;
  }

  std::vector<size_t> new_triangles;
  std::vector<size_t> new_uv_indices;
  std::vector<int> new_material_ids;
  std::vector<std::string> new_face_metadata;
  new_triangles.reserve(triangles.size() - 3 * collapsed + 3 * merged.size());
  if (has_uvs)
  {
    new_uv_indices.reserve(new_triangles.capacity());
  }
  new_material_ids.reserve(n_triangles - collapsed + merged.size());
  if (store_metadata_)
  {
    new_face_metadata.reserve(n_triangles - collapsed + merged.size());
  }

  for (size_t t = 0; t < n_triangles; ++t)
  {
    if (triangle_removed[t])
    {
      continue;
    }
    new_triangles.push_back(triangles[3 * t + 0]);
    new_triangles.push_back(triangles[3 * t + 1]);
    new_triangles.push_back(triangles[3 * t + 2]);
    if (has_uvs)
    {
      new_uv_indices.push_back(uv_indices[3 * t + 0]);
      new_uv_indices.push_back(uv_indices[3 * t + 1]);
      new_uv_indices.push_back(uv_indices[3 * t + 2]);
    }
    new_material_ids.push_back(t < material_ids.size() ? material_ids[t] : -1);
    if (store_metadata_)
    {
      new_face_metadata.push_back(t < face_metadata.size() ? face_metadata[t] : std::string("{}"));
    }
  }

  for (const MergedTriangle& tri : merged)
  {
    new_triangles.push_back(tri.v0);
    new_triangles.push_back(tri.v1);
    new_triangles.push_back(tri.v2);
    if (has_uvs)
    {
      new_uv_indices.push_back(tri.uv0);
      new_uv_indices.push_back(tri.uv1);
      new_uv_indices.push_back(tri.uv2);
    }
    new_material_ids.push_back(tri.material_id);
    if (store_metadata_)
    {
      new_face_metadata.push_back(tri.face_metadata);
    }
  }

  triangles.swap(new_triangles);
  if (has_uvs)
  {
    uv_indices.swap(new_uv_indices);
  }
  else
  {
    uv_indices.clear();
  }
  material_ids.swap(new_material_ids);
  if (store_metadata_)
  {
    face_metadata.swap(new_face_metadata);
  }
  else
  {
    face_metadata.clear();
  }

  // Drop collapsed flexible vertices (and any other unused verts).
  removeIsolatedVertices();
  return collapsed;
}

#ifdef USE_CGAL
std::vector<Surface_mesh::Halfedge_index> get_boundary_cycle(Surface_mesh::Halfedge_index h, const Surface_mesh& mesh)
{
  std::vector<Surface_mesh::Halfedge_index> cycle;

  Surface_mesh::Halfedge_index start = h;
  Surface_mesh::Halfedge_index cur = h;

  do
  {
    cycle.push_back(cur);
    cur = mesh.next(cur);
  } while (cur != start);

  return cycle;
}
#endif

void VoronoiMesh::patchHoles(
  std::function<void(size_t)> tri_callback, std::function<void(size_t, size_t)> vertex_callback, int material_id)
{
#ifdef USE_CGAL
  Surface_mesh mesh;

  // --- 1. Copy vertices (1:1 mapping) ---
  std::vector<Surface_mesh::Vertex_index> vmap(vertices.size());

  for (size_t i = 0; i < vertices.size(); ++i)
  {
    const auto& v = vertices[i];
    vmap[i] = mesh.add_vertex(Point_3(v[0], v[1], v[2]));
  }

  // --- 2. Copy faces ---
  const size_t original_face_count = triangles.size() / 3;

  for (size_t i = 0; i < triangles.size(); i += 3)
  {
    mesh.add_face(vmap[triangles[i + 0]], vmap[triangles[i + 1]], vmap[triangles[i + 2]]);
  }

  // --- 3. Fill holes ---
  std::unordered_set<Surface_mesh::Halfedge_index> visited;
  for (auto h : mesh.halfedges())
  {
    if (!mesh.is_border(h) || visited.count(h))
      continue;

    auto cycle = get_boundary_cycle(h, mesh);

    // Mark all halfedges in this cycle as visited
    for (auto hh : cycle)
      visited.insert(hh);

    if (cycle.size() < 3)
      continue;

    PMP::triangulate_hole(mesh, h, PMP::parameters::use_delaunay_triangulation(true));
  }

  // --- 4. Append ONLY new triangles ---
  size_t face_index = 0;
  for (auto f : mesh.faces())
  {
    if (face_index++ < original_face_count)
      continue;

    auto h = mesh.halfedge(f);
    size_t corner_index = triangles.size();
    size_t tri_index = corner_index / 3;
    std::array<size_t, 3> v_indices;

    for (int i = 0; i < 3; ++i)
    {
      v_indices[i] = mesh.target(h).idx();

      if (normal_mode == PerTriangleCorner)
      {
        normals.push_back(glm::dvec3 { 0.0, 0.0, 0.0 });
      }
      h = mesh.next(h);
    }
    addTriangle(v_indices[0], v_indices[1], v_indices[2], material_id);

    tri_callback(tri_index);
    for (int i = 0; i < 3; ++i)
    {
      vertex_callback(v_indices[i], corner_index + i);
    }
  }
#endif
}

std::vector<glm::dvec3> kinDS::VoronoiMesh::computeVertexNormals()
{
  std::vector<glm::dvec3> vertex_normals(vertices.size(), glm::dvec3 { 0.0, 0.0, 0.0 });
  // Accumulate triangle normals into vertex normals
  for (size_t i = 0; i + 2 < triangles.size(); i += 3)
  {
    size_t i0 = triangles[i];
    size_t i1 = triangles[i + 1];
    size_t i2 = triangles[i + 2];

    const glm::dvec3& p0 = vertices[i0];
    const glm::dvec3& p1 = vertices[i1];
    const glm::dvec3& p2 = vertices[i2];

    glm::dvec3 e1 = p1 - p0;
    glm::dvec3 e2 = p2 - p0;

    // Unnormalized triangle normal (area-weighted)
    glm::dvec3 triNormal = glm::cross(e1, e2);

    vertex_normals[i0] += triNormal;
    vertex_normals[i1] += triNormal;
    vertex_normals[i2] += triNormal;
  }

  // Normalize the accumulated vertex normals
  for (glm::dvec3& n : vertex_normals)
  {
    if (glm::length2(n) != 0.0)
    {
      n = glm::normalize(n);
    }
  }

  return vertex_normals;
}

void kinDS::VoronoiMesh::computeNormals(NormalMode normal_mode)
{
  // Ensure normals has the correct size and is zero-initialized
  this->normal_mode = normal_mode;

  if (normal_mode == NoNormals)
  {
    normals.clear();
  }
  else if (normal_mode == PerVertex)
  {
    normals = computeVertexNormals();
  }
  else if (normal_mode == PerTriangleCorner)
  {
    computeNormalsWithCreaseAngle(30.0);
    return;
  }
  validateNormalCount("VoronoiMesh::computeNormals");
}

void kinDS::VoronoiMesh::setCornerNormals(std::vector<glm::dvec3> corner_normals)
{
  if (corner_normals.size() != triangles.size())
  {
    throw std::runtime_error("VoronoiMesh::setCornerNormals: expected " + std::to_string(triangles.size())
      + " corner normals, got " + std::to_string(corner_normals.size()));
  }
  normal_mode = PerTriangleCorner;
  normals = std::move(corner_normals);
  validateNormalCount("VoronoiMesh::setCornerNormals");
}

void kinDS::VoronoiMesh::computeNormalsWithCreaseAngle(const double crease_angle_degrees)
{
  normal_mode = PerTriangleCorner;
  const size_t triangle_count = getTriangleCount();
  normals.assign(triangles.size(), glm::dvec3(0.0));
  if (triangle_count == 0)
  {
    return;
  }

  const double cos_threshold = std::cos(glm::radians(crease_angle_degrees));

  std::vector<glm::dvec3> face_normals(triangle_count);
  for (size_t triangle_index = 0; triangle_index < triangle_count; ++triangle_index)
  {
    const size_t i0 = triangles[3 * triangle_index];
    const size_t i1 = triangles[3 * triangle_index + 1];
    const size_t i2 = triangles[3 * triangle_index + 2];
    glm::dvec3 normal = glm::cross(vertices[i1] - vertices[i0], vertices[i2] - vertices[i0]);
    face_normals[triangle_index] = glm::length2(normal) > 1e-24 ? glm::normalize(normal) : glm::dvec3(0.0, 1.0, 0.0);
  }

  std::vector<std::vector<std::pair<size_t, int>>> vertex_corners(vertices.size());
  for (size_t triangle_index = 0; triangle_index < triangle_count; ++triangle_index)
  {
    for (int corner = 0; corner < 3; ++corner)
    {
      vertex_corners[triangles[3 * triangle_index + corner]].emplace_back(triangle_index, corner);
    }
  }

  const auto edge_neighbor = [&](const size_t triangle_index, const int corner, const int step) -> size_t {
    return triangles[3 * triangle_index + (corner + step + 3) % 3];
  };

  for (size_t vertex_index = 0; vertex_index < vertices.size(); ++vertex_index)
  {
    const auto& corners = vertex_corners[vertex_index];
    const size_t corner_count = corners.size();
    if (corner_count == 0)
    {
      continue;
    }

    std::vector<std::vector<size_t>> adjacency(corner_count);
    for (size_t i = 0; i < corner_count; ++i)
    {
      const auto [triangle_i, corner_i] = corners[i];
      const size_t edge_a = edge_neighbor(triangle_i, corner_i, 1);
      const size_t edge_b = edge_neighbor(triangle_i, corner_i, 2);
      for (size_t j = i + 1; j < corner_count; ++j)
      {
        const auto [triangle_j, corner_j] = corners[j];
        if (triangle_i == triangle_j)
        {
          continue;
        }
        const size_t edge_c = edge_neighbor(triangle_j, corner_j, 1);
        const size_t edge_d = edge_neighbor(triangle_j, corner_j, 2);
        const bool shares_edge = edge_a == edge_c || edge_a == edge_d || edge_b == edge_c || edge_b == edge_d;
        if (!shares_edge)
        {
          continue;
        }
        if (glm::dot(face_normals[triangle_i], face_normals[triangle_j]) >= cos_threshold)
        {
          adjacency[i].push_back(j);
          adjacency[j].push_back(i);
        }
      }
    }

    std::vector<int> component(corner_count, -1);
    int component_count = 0;
    std::vector<int> stack;
    for (size_t i = 0; i < corner_count; ++i)
    {
      if (component[i] >= 0)
      {
        continue;
      }
      component[i] = component_count;
      stack.push_back(static_cast<int>(i));
      while (!stack.empty())
      {
        const int current = stack.back();
        stack.pop_back();
        for (const size_t neighbor : adjacency[current])
        {
          if (component[neighbor] < 0)
          {
            component[neighbor] = component_count;
            stack.push_back(static_cast<int>(neighbor));
          }
        }
      }
      ++component_count;
    }

    std::vector<glm::dvec3> averaged(component_count, glm::dvec3(0.0));
    std::vector<int> counts(component_count, 0);
    for (size_t i = 0; i < corner_count; ++i)
    {
      const int comp = component[i];
      averaged[comp] += face_normals[corners[i].first];
      ++counts[comp];
    }
    for (int comp = 0; comp < component_count; ++comp)
    {
      if (counts[comp] > 0 && glm::length2(averaged[comp]) > 1e-24)
      {
        averaged[comp] = glm::normalize(averaged[comp]);
      }
      else
      {
        averaged[comp] = glm::dvec3(0.0, 1.0, 0.0);
      }
    }

    for (size_t i = 0; i < corner_count; ++i)
    {
      const auto [triangle_index, corner] = corners[i];
      normals[3 * triangle_index + corner] = averaged[component[i]];
    }
  }

  validateNormalCount("VoronoiMesh::computeNormalsWithCreaseAngle");
}

std::array<double, 3> kinDS::VoronoiMesh::computeBarycentricCoordinates(size_t triangle_index, const glm::dvec3& point) const
{
  return barycentricCoordinates(vertices[triangles[3 * triangle_index]], vertices[triangles[3 * triangle_index + 1]],
    vertices[triangles[3 * triangle_index + 2]], point);
}

const std::vector<glm::dvec3>& VoronoiMesh::getVertices() const { return vertices; }

std::vector<glm::dvec3>& kinDS::VoronoiMesh::getVertices() { return vertices; }

const std::vector<size_t>& VoronoiMesh::getTriangles() const { return triangles; }

std::vector<size_t>& kinDS::VoronoiMesh::getTriangles() { return triangles; }

const std::vector<glm::dvec3>& VoronoiMesh::getNormals() const { return normals; }

std::vector<glm::dvec3>& VoronoiMesh::getNormals() { return normals; }

const std::vector<glm::dvec3>& VoronoiMesh::getUVs() const { return uvs; }

std::vector<glm::dvec3>& VoronoiMesh::getUVs() { return uvs; }

const std::vector<size_t>& VoronoiMesh::getUVIndices() const { return uv_indices; }

void VoronoiMesh::printStatistics() const
{
  KINDS_DEBUG("\nuv_indices.size(): " << uv_indices.size() << "\ntriangles.size(): " << triangles.size()
                                      << "\nuvs.size(): " << uvs.size() << "\nvertices.size():" << vertices.size());
}

bool kinDS::VoronoiMesh::hasValidUVIndex(size_t triangle_vertex_index) const
{
  if (triangle_vertex_index >= uv_indices.size())
  {
    KINDS_ERROR("triangle_vertex_index: " << triangle_vertex_index << "\nuv_indices.size(): " << uv_indices.size()
                                          << "\ntriangles.size(): " << triangles.size()
                                          << "\nuvs.size(): " << uvs.size() << "\nvertices.size():" << vertices.size());
    throw std::out_of_range("Triangle vertex index out of range when checking UV index.");
  }
  return uv_indices[triangle_vertex_index] < uvs.size();
}

std::vector<size_t> kinDS::VoronoiMesh::findTriangleCorners(size_t vertex_index, bool stop_early) const
{
  std::vector<size_t> corner_indices;

  for (size_t i = 0; i < triangles.size(); i++)
  {
    if (triangles[i] == vertex_index)
    {
      corner_indices.push_back(i);
      if (stop_early)
      {
      }
    }
  }

  return corner_indices;
}

std::vector<size_t>& VoronoiMesh::getUVIndices() { return uv_indices; }

std::vector<size_t> VoronoiMesh::removeIsolatedVertices()
{
  const size_t n_vertices = vertices.size();

  // 1. Mark used vertices
  std::vector<bool> used(n_vertices, false);
  for (size_t idx : triangles)
  {
    if (idx < n_vertices)
    {
      used[idx] = true;
    }
  }

  // 2. Build old -> new index map
  std::vector<size_t> remap(n_vertices, size_t(-1));
  size_t new_count = 0;

  for (size_t i = 0; i < n_vertices; ++i)
  {
    if (used[i])
    {
      remap[i] = new_count++;
    }
    else
    {
      // KINDS_DEBUG("Found unused vertex at index: " << i);
    }
  }

  // Early out: nothing to remove
  if (new_count == n_vertices)
  {
    return remap;
  }

  // 3. Compact vertex data
  std::vector<glm::dvec3> new_vertices;
  new_vertices.reserve(new_count);
  std::vector<std::string> new_vertex_metadata;
  new_vertex_metadata.reserve(new_count);
  std::vector<glm::dvec3> new_vertex_colors;
  new_vertex_colors.reserve(new_count);
  std::vector<double> new_vertex_kinetic_times;
  new_vertex_kinetic_times.reserve(new_count);
  std::vector<glm::dvec3> new_vertex_semantic_uvs;
  new_vertex_semantic_uvs.reserve(new_count);
  std::vector<bool> new_vertex_is_flexible;
  new_vertex_is_flexible.reserve(new_count);

  for (size_t i = 0; i < n_vertices; ++i)
  {
    if (used[i])
    {
      new_vertices.push_back(vertices[i]);
      if (i < vertex_metadata.size())
      {
        new_vertex_metadata.push_back(vertex_metadata[i]);
      }
      else
      {
        new_vertex_metadata.push_back("{}");
      }
      if (i < vertex_colors.size())
      {
        new_vertex_colors.push_back(vertex_colors[i]);
      }
      else
      {
        new_vertex_colors.push_back(glm::dvec3(1.0));
      }
      if (i < vertex_kinetic_times_.size())
      {
        new_vertex_kinetic_times.push_back(vertex_kinetic_times_[i]);
      }
      else
      {
        new_vertex_kinetic_times.push_back(std::numeric_limits<double>::quiet_NaN());
      }
      if (i < vertex_semantic_uvs_.size())
      {
        new_vertex_semantic_uvs.push_back(vertex_semantic_uvs_[i]);
      }
      else
      {
        new_vertex_semantic_uvs.push_back(glm::dvec3(std::numeric_limits<double>::quiet_NaN()));
      }
      new_vertex_is_flexible.push_back(i < vertex_is_flexible_.size() && vertex_is_flexible_[i]);
    }
  }
  vertices.swap(new_vertices);
  vertex_metadata.swap(new_vertex_metadata);
  vertex_colors.swap(new_vertex_colors);
  vertex_kinetic_times_.swap(new_vertex_kinetic_times);
  vertex_semantic_uvs_.swap(new_vertex_semantic_uvs);
  vertex_is_flexible_.swap(new_vertex_is_flexible);
  if (!store_metadata_)
  {
    vertex_metadata.clear();
  }

  // 4. Compact per-vertex normals if needed
  if (normal_mode == PerVertex)
  {
    std::vector<glm::dvec3> new_normals;
    new_normals.reserve(new_count);

    for (size_t i = 0; i < n_vertices; ++i)
    {
      if (used[i])
      {
        new_normals.push_back(normals[i]);
      }
    }
    normals.swap(new_normals);
  }

  // 5. Remap triangle vertex indices
  for (size_t& idx : triangles)
  {
    idx = remap[idx];
  }

  return remap;
}

size_t VoronoiMesh::removeDegenerateTriangles()
{
  const size_t n_triangles = triangles.size() / 3;
  std::vector<size_t> new_triangles;
  std::vector<size_t> new_uv_indices;
  std::vector<int> new_material_ids;
  std::vector<std::string> new_face_metadata;

  new_triangles.reserve(triangles.size());
  if (!uv_indices.empty())
  {
    new_uv_indices.reserve(uv_indices.size());
  }
  new_material_ids.reserve(material_ids.size());
  new_face_metadata.reserve(face_metadata.size());

  for (size_t t = 0; t < n_triangles; ++t)
  {
    size_t i0 = triangles[3 * t + 0];
    size_t i1 = triangles[3 * t + 1];
    size_t i2 = triangles[3 * t + 2];

    // Check for duplicate vertices
    if (i0 != i1 && i0 != i2 && i1 != i2)
    {
      new_triangles.push_back(i0);
      new_triangles.push_back(i1);
      new_triangles.push_back(i2);

      // Keep uv_indices in sync if present
      if (!uv_indices.empty())
      {
        new_uv_indices.push_back(uv_indices[3 * t + 0]);
        new_uv_indices.push_back(uv_indices[3 * t + 1]);
        new_uv_indices.push_back(uv_indices[3 * t + 2]);
      }

      const int material_id = t < material_ids.size() ? material_ids[t] : -1;
      new_material_ids.push_back(material_id);
      if (t < face_metadata.size())
      {
        new_face_metadata.push_back(face_metadata[t]);
      }
      else
      {
        new_face_metadata.push_back("{}");
      }
    }
  }

  triangles.swap(new_triangles);
  if (!uv_indices.empty())
  {
    uv_indices.swap(new_uv_indices);
  }
  material_ids.swap(new_material_ids);
  face_metadata.swap(new_face_metadata);
  if (!store_metadata_)
  {
    face_metadata.clear();
  }

  return n_triangles - (triangles.size() / 3);
}

void VoronoiMesh::ensureFaceMetadataSize(const std::string& fill_value)
{
  if (!store_metadata_)
  {
    return;
  }

  const size_t tri_count = triangles.size() / 3;
  if (face_metadata.size() < tri_count)
  {
    face_metadata.resize(tri_count, fill_value);
  }
  else if (face_metadata.size() > tri_count)
  {
    face_metadata.resize(tri_count);
  }
}

const glm::dvec3& kinDS::VoronoiMesh::getNormal(size_t triangle_vertex_index) const
{
  if (normal_mode == NoNormals)
  {
    throw std::runtime_error("Cannot access normals on a mesh with NormalMode::NoNormals.");
  }
  if (normal_mode == PerTriangleCorner)
  {
    return normals[triangle_vertex_index];
  }
  else
  {
    return normals[triangles[triangle_vertex_index]];
  }
}

const glm::dvec3& kinDS::VoronoiMesh::getUV(size_t triangle_vertex_index) const
{
  return uvs[uv_indices[triangle_vertex_index]];
}

void kinDS::VoronoiMesh::setNormal(const glm::dvec3& normal, size_t triangle_vertex_index)
{
  if (normal_mode == NoNormals)
  {
    throw std::runtime_error("Cannot set normals on a mesh with NormalMode::NoNormals.");
  }
  if (normal_mode == PerTriangleCorner)
  {
    if (triangle_vertex_index >= normals.size())
    {
      throw std::out_of_range("Triangle vertex index out of range when setting normal.");
    }
    normals[triangle_vertex_index] = normal;
  }
  else
  {
    size_t vertex_index = triangles[triangle_vertex_index];
    if (vertex_index >= normals.size())
    {
      throw std::out_of_range("Vertex index out of range when setting normal.");
    }
    normals[vertex_index] = normal;
  }
}

void kinDS::VoronoiMesh::setUV(const glm::dvec3& uv, size_t triangle_vertex_index)
{
  if (triangle_vertex_index >= uv_indices.size())
  {
    throw std::out_of_range("Triangle corner index out of range when setting UV.");
  }
  // UV references belong to corners. Always detach the edited corner so seam correction or hole patching cannot
  // accidentally mutate another corner that happens to share the same UV-pool entry.
  uv_indices[triangle_vertex_index] = addUV(uv);
}

void kinDS::VoronoiMesh::validateUVLayout(const std::string& context) const
{
  const std::string prefix = context.empty() ? std::string {} : context + ": ";
  if (uv_indices.size() != triangles.size())
  {
    throw std::runtime_error(prefix + "UV-index count " + std::to_string(uv_indices.size())
      + " does not match triangle-corner count " + std::to_string(triangles.size()) + ".");
  }
  for (size_t corner = 0; corner < uv_indices.size(); ++corner)
  {
    const size_t uv_index = uv_indices[corner];
    if (uv_index != std::numeric_limits<size_t>::max() && uv_index >= uvs.size())
    {
      throw std::runtime_error(prefix + "UV index " + std::to_string(uv_index) + " at triangle corner "
        + std::to_string(corner) + " is out of range for " + std::to_string(uvs.size()) + " UV coordinates.");
    }
  }
}

NormalMode kinDS::VoronoiMesh::getNormalMode() const { return normal_mode; }

void kinDS::VoronoiMesh::validateNormalCount(const std::string& context) const
{
  size_t expected_normal_count = 0;
  if (normal_mode == NormalMode::PerVertex)
  {
    expected_normal_count = vertices.size();
  }
  else if (normal_mode == NormalMode::PerTriangleCorner)
  {
    expected_normal_count = triangles.size();
  }
  if (normals.size() == expected_normal_count)
  {
    return;
  }

  std::ostringstream oss;
  if (!context.empty())
  {
    oss << context << ": ";
  }
  oss << "normal count mismatch for " << normalModeToString(normal_mode) << " mesh; normals.size()=" << normals.size()
      << ", expected=" << expected_normal_count << ", vertices.size()=" << vertices.size()
      << ", triangle_corner_count=" << triangles.size() << ".";
  throw std::runtime_error(oss.str());
}

void kinDS::VoronoiMesh::checkForDegenerateTriangles() const
{
  size_t degenerate_faces = 0;
  for (size_t i = 0; i < triangles.size(); i += 3)
  {
    const auto& p0 = vertices[triangles[i]];
    const auto& p1 = vertices[triangles[i + 1]];
    const auto& p2 = vertices[triangles[i + 2]];

    // compute squared area via cross product

    double area2 = glm::length2(glm::cross((p1 - p0), (p2 - p0)));

    if (area2 < 1e-20)
    {
      degenerate_faces++;
    }
  }
  if (degenerate_faces != 0)
    std::cerr << "Degenerate triangles: " << degenerate_faces << "/" << (triangles.size() / 3) << "\n";
}
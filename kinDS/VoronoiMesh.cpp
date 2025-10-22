#include "VoronoiMesh.hpp"
#include "../eigen/Eigen/Dense"
#include "Logger.hpp"
#include <unordered_map>

using namespace kinDS;

size_t VoronoiMesh::addVertex(double x, double y, double z)
{
  size_t index = vertices.size();
  vertices.emplace_back(Point<3> { x, y, z });
  return index;
}

void VoronoiMesh::addTriangle(size_t v1, size_t v2, size_t v3)
{
  vertex_indices.push_back(v1);
  vertex_indices.push_back(v2);
  vertex_indices.push_back(v3);

  // use maximum size for UV indices if not provided
  uv_indices.push_back(std::numeric_limits<size_t>::max());
  uv_indices.push_back(std::numeric_limits<size_t>::max());
  uv_indices.push_back(std::numeric_limits<size_t>::max());
}

void VoronoiMesh::addTriangle(size_t v1, size_t v2, size_t v3, size_t uv1, size_t uv2, size_t uv3)
{
  vertex_indices.push_back(v1);
  vertex_indices.push_back(v2);
  vertex_indices.push_back(v3);
  uv_indices.push_back(uv1);
  uv_indices.push_back(uv2);
  uv_indices.push_back(uv3);
}

void VoronoiMesh::addNormal(double nx, double ny, double nz)
{
  normals.emplace_back(Vector<3> { nx, ny, nz });
}

void VoronoiMesh::addUV(double u, double v)
{
  uvs.emplace_back(Vector<2> { u, v });
}

void VoronoiMesh::startNewGroup()
{
  group_offsets.push_back(vertex_indices.size()); // Store the current vertex index count as a new group offset
}

void VoronoiMesh::setGroupOffsets(const std::vector<size_t>& offsets)
{
  group_offsets = offsets;
}

VoronoiMesh& VoronoiMesh::operator+=(const VoronoiMesh& other)
{
  size_t old_vertices_size = vertices.size();
  vertices.insert(vertices.end(), other.vertices.begin(), other.vertices.end());

  size_t old_vertex_indices_size = vertex_indices.size();
  vertex_indices.insert(vertex_indices.end(), other.vertex_indices.begin(), other.vertex_indices.end());

  std::transform(vertex_indices.begin() + old_vertex_indices_size,
    vertex_indices.end(),
    vertex_indices.begin() + old_vertex_indices_size,
    [&](size_t index)
    { return index + old_vertices_size; });

  normals.insert(normals.end(), other.normals.begin(), other.normals.end());

  size_t old_uvs_size = uvs.size();
  uvs.insert(uvs.end(), other.uvs.begin(), other.uvs.end());

  size_t old_uv_indices_size = uv_indices.size();
  uv_indices.insert(uv_indices.end(), other.uv_indices.begin(), other.uv_indices.end());

  std::transform(uv_indices.begin() + old_uv_indices_size,
    uv_indices.end(),
    uv_indices.begin() + old_uv_indices_size,
    [&](size_t index)
    { return index + old_uvs_size; });

  size_t old_group_count = group_offsets.size();
  group_offsets.insert(group_offsets.end(), other.group_offsets.begin(), other.group_offsets.end());

  std::transform(group_offsets.begin() + old_group_count,
    group_offsets.end(),
    group_offsets.begin() + old_group_count,
    [&](size_t offset)
    { return offset + old_vertex_indices_size; });

  return *this;
}

void VoronoiMesh::flipOrientation()
{
  // Flip the orientation of each triangle by swapping the second and third vertex indices
  for (size_t i = 0; i < vertex_indices.size(); i += 3)
  {
    std::swap(vertex_indices[i + 1], vertex_indices[i + 2]);
    std::swap(uv_indices[i + 1], uv_indices[i + 2]);
  }
}

void VoronoiMesh::mergeDuplicateVertices(double epsilon)
{
  const double inv_eps = (epsilon > 0.0) ? 1.0 / epsilon : 0.0;
  std::unordered_map<Point<3>, size_t, VoronoiMesh::Vec3iHash> grid;
  std::vector<Point<3>> newVerts;
  newVerts.reserve(vertices.size());

  std::vector<size_t> remap(vertices.size(), size_t(-1));

  for (size_t i = 0; i < vertices.size(); ++i)
  {
    const auto& v = vertices[i];

    // Quantize vertex for approximate matching
    Point<3> key;
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

    auto it = grid.find(key);
    if (it == grid.end())
    {
      size_t newIndex = newVerts.size();
      grid[key] = newIndex;
      newVerts.push_back(v);
      remap[i] = newIndex;
    }
    else
    {
      remap[i] = it->second;
    }
  }

  // Remap triangle indices
  for (size_t& idx : vertex_indices)
  {
    idx = remap[idx];
  }

  vertices.swap(newVerts);
}

const std::vector<kinDS::Point<3>>& VoronoiMesh::getVertices() const
{
  return vertices;
}

std::vector<Point<3>>& kinDS::VoronoiMesh::getVertices()
{
  return vertices;
}

const std::vector<size_t>& VoronoiMesh::getVertexIndices() const
{
  return vertex_indices;
}

std::vector<size_t>& kinDS::VoronoiMesh::getVertexIndices()
{
  return vertex_indices;
}

const std::vector<Vector<3>>& VoronoiMesh::getNormals() const
{
  return normals;
}

const std::vector<Vector<2>>& VoronoiMesh::getUVs() const
{
  return uvs;
}

const std::vector<size_t>& VoronoiMesh::getUVIndices() const
{
  return uv_indices;
}

void kinDS::VoronoiMesh::checkForDegenerateTriangles() const
{
  std::vector<size_t> degenerate_faces;
  for (size_t i = 0; i < vertex_indices.size(); i += 3)
  {
    const auto& p0 = vertices[vertex_indices[i]];
    const auto& p1 = vertices[vertex_indices[i + 1]];
    const auto& p2 = vertices[vertex_indices[i + 2]];

    // compute squared area via cross product
    Eigen::Vector3d v0(p0[0], p0[1], p0[2]);
    Eigen::Vector3d v1(p1[0], p1[1], p1[2]);
    Eigen::Vector3d v2(p2[0], p2[1], p2[2]);
    double area2 = ((v1 - v0).cross(v2 - v0)).squaredNorm();

    if (area2 < 1e-20)
    {
      degenerate_faces.push_back(i);
    }
  }
  // if (!degenerate_faces.empty())
  std::cerr << "Degenerate triangles: " << degenerate_faces.size() << "\n";
}

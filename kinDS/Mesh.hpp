#pragma once
#include "Point.hpp"
#include <vector>

namespace kinDS
{

class Mesh
{
 private:
  std::vector<Point<3>> vertices; // Stores vertex coordinates
  std::vector<size_t> vertex_indices; // Stores indices of vertices forming triangles
  std::vector<Vector<3>> normals; // Stores normal vectors for each vertex
  std::vector<Vector<2>> uvs; // Stores texture coordinates for each vertex
  std::vector<size_t> uv_indices; // Stores indices of texture coordinates for faces
  std::vector<size_t> group_offsets; // Offsets for groups of triangles, if needed

 public:
  Mesh() = default;
  Mesh(std::vector<Point<3>> vertices,
    std::vector<size_t> vertex_indices,
    std::vector<Vector<3>> normals = {},
    std::vector<Vector<2>> uvs = {},
    std::vector<size_t> uv_indices = {})
    : vertices(std::move(vertices))
    , vertex_indices(std::move(vertex_indices))
    , normals(std::move(normals))
    , uvs(std::move(uvs))
    , uv_indices(std::move(uv_indices))
    , group_offsets(std::vector<size_t>(1, 0)) // Initialize with one group offset at 0
  {
  }

  ~Mesh() = default;

  // Add methods to manipulate the mesh, such as adding vertices, triangles, normals, and UVs
  void addVertex(double x, double y, double z);
  void addTriangle(size_t v1, size_t v2, size_t v3);
  void addTriangle(size_t v1, size_t v2, size_t v3, size_t uv1, size_t uv2, size_t uv3);
  void addNormal(double nx, double ny, double nz);
  void addUV(double u, double v);
  void startNewGroup()
  {
    group_offsets.push_back(vertex_indices.size()); // Store the current vertex index count as a new group offset
  }
  // Methods to retrieve mesh data
  const std::vector<Point<3>>& getVertices() const;
  const std::vector<size_t>& getVertexIndices() const;
  const std::vector<Vector<3>>& getNormals() const;
  const std::vector<Vector<2>>& getUVs() const;
  const std::vector<size_t>& getUVIndices() const;
  const std::vector<size_t>& getGroupOffsets() const { return group_offsets; }
};
}

#include "VoronoiMesh.hpp"

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

const std::vector<kinDS::Point<3>>& VoronoiMesh::getVertices() const
{
  return vertices;
}

const std::vector<size_t>& VoronoiMesh::getVertexIndices() const
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

#pragma once
#include <algorithm>
#include <cmath>
#include <functional>
#include <glm/glm.hpp>
#include <limits>
#include <string>
#include <vector>

namespace kinDS
{

enum NormalMode
{
  NoNormals,
  PerVertex,
  PerTriangleCorner
};

class VoronoiMesh
{
 private:
  std::vector<glm::dvec3> vertices; // Stores vertex coordinates
  std::vector<size_t> triangles; // Stores indices of vertices forming triangles
  std::vector<glm::dvec3> normals; // Stores normal vectors for each triangle or vertex depending on node
  std::vector<glm::dvec3> uvs; // Stores texture coordinates
  std::vector<glm::dvec3> vertex_colors; // Optional debug colors per vertex (RGB in [0,1])
  std::vector<size_t> uv_indices; // Stores indices of texture coordinates for face corners
  std::vector<size_t> group_offsets; // Triangle-index offsets marking the start of each face group
  std::vector<std::string> material_names; // Stores material names - needed for exporting
  std::vector<int> material_ids; // Stores material IDs per triangle
  // Optional JSON-like metadata kept separate from raw geometry/index buffers.
  std::vector<std::string> vertex_metadata;
  std::vector<std::string> face_metadata;

  /// Profile-plane (x,y) per vertex when stored positions are object-space (ear clipping).
  std::vector<glm::dvec2> profile_plane_xy_;

  NormalMode normal_mode;

  /// Kinetic time when this meshlet was first registered / last rebuilt in @ref SegmentBuilder (NaN if unset).
  double creation_kinetic_time_ = std::numeric_limits<double>::quiet_NaN();
  bool store_metadata_ = true;

  struct Vec3iHash
  {
    std::size_t operator()(const glm::dvec3& v) const noexcept
    {
      std::size_t h1 = std::hash<int> {}(v[0]);
      std::size_t h2 = std::hash<int> {}(v[1]);
      std::size_t h3 = std::hash<int> {}(v[2]);
      return h1 ^ (h2 << 1) ^ (h3 << 2);
    }
  };

 public:
  VoronoiMesh(std::vector<std::string> material_names = {}, NormalMode normal_mode = NoNormals)
    : material_names(std::move(material_names))
    , normal_mode(normal_mode) { };

  VoronoiMesh(std::vector<glm::dvec3> vertices, std::vector<size_t> triangles,
    std::vector<std::string> material_names = {}, std::vector<glm::dvec3> normals = {},
    std::vector<glm::dvec3> uvs = {}, std::vector<size_t> uv_indices = {}, NormalMode normal_mode = NoNormals)
    : vertices(std::move(vertices))
    , triangles(std::move(triangles))
    , material_names(std::move(material_names))
    , normals(std::move(normals))
    , uvs(std::move(uvs))
    , uv_indices(std::move(uv_indices))
    , group_offsets(std::vector<size_t>(1, 0))
    , // Initialize with one group offset at 0
    normal_mode(normal_mode)
  {
    if (this->uv_indices.empty())
    {
      // If no UV indices are provided, set the same length as vertex indices
      this->uv_indices.resize(this->triangles.size(), std::numeric_limits<size_t>::max());
    }
    this->vertex_metadata.resize(this->vertices.size(), "{}");
    this->face_metadata.resize(this->triangles.size() / 3, "{}");
  }

  ~VoronoiMesh() = default;

  // methods to manipulate the mesh, such as adding vertices, triangles, normals, and UVs
  size_t addVertex(
    double x, double y, double z, const std::string& metadata = "{}", const glm::dvec3& color = glm::dvec3(1.0));
  size_t addVertex(const glm::dvec3& p, const std::string& metadata = "{}", const glm::dvec3& color = glm::dvec3(1.0));
  void replaceVertex(size_t index, const glm::dvec3& new_position);
  size_t addTriangle(size_t v1, size_t v2, size_t v3, int material_id = -1, const std::string& metadata = "{}");
  size_t addTriangle(size_t v1, size_t v2, size_t v3, size_t uv1, size_t uv2, size_t uv3, int material_id = -1,
    const std::string& metadata = "{}");
  size_t addNormal(double nx, double ny, double nz);
  size_t addNormal(const glm::dvec3& n);
  void replaceNormal(size_t index, const glm::dvec3& new_normal);
  size_t addUV(double u, double v, double w);
  size_t addUV(glm::dvec3 uv);
  void startNewGroup();
  void setGroupOffsets(const std::vector<size_t>& offsets);
  void setMaterialNames(std::vector<std::string> names) { material_names = std::move(names); }
  VoronoiMesh& operator+=(const VoronoiMesh& other);
  void flipOrientation();

  // Merge duplicate vertices (within epsilon) and update triangle indices
  std::vector<size_t> mergeDuplicateVertices(double epsilon = 0.0);
  void patchHoles(std::function<void(size_t)> tri_callback, std::function<void(size_t, size_t)> vertex_callback,
    int material_id = -1);

  // compute normals
  void computeNormals(NormalMode normal_mode = PerVertex);

  std::array<double, 3> computeBarycentricCoordinates(size_t triangle_index, glm::dvec3& point) const;

  // Methods to retrieve mesh data
  const std::vector<glm::dvec3>& getVertices() const;
  std::vector<glm::dvec3>& getVertices();
  const std::vector<size_t>& getTriangles() const;
  std::vector<size_t>& getTriangles();
  const std::vector<glm::dvec3>& getNormals() const;
  std::vector<glm::dvec3>& getNormals();
  const std::vector<glm::dvec3>& getUVs() const;
  std::vector<glm::dvec3>& getUVs();
  const std::vector<size_t>& getUVIndices() const;
  void printStatistics() const;
  bool hasValidUVIndex(size_t triangle_vertex_index) const;
  std::vector<size_t> findTriangleCorners(size_t vertex_index, bool stop_early = false) const;

  std::vector<size_t>& getUVIndices();

  const std::vector<int>& getMaterialIDs() const { return material_ids; }

  const std::vector<std::string>& getMaterialNames() const { return material_names; }
  const std::vector<std::string>& getVertexMetadata() const { return vertex_metadata; }
  const std::vector<std::string>& getFaceMetadata() const { return face_metadata; }
  const std::vector<glm::dvec3>& getVertexColors() const { return vertex_colors; }

  std::vector<size_t> removeIsolatedVertices();

  void removeDegenerateTriangles();

  /**
   * Provides a mode-independent way to get the normal from a triangle vertex
   * @param triangle_vertex_index index corresponding to the indices in the triangle buffer
   * @return
   */
  const glm::dvec3& getNormal(size_t triangle_vertex_index) const;
  const glm::dvec3& getUV(size_t triangle_vertex_index) const;

  void setNormal(const glm::dvec3& normal, size_t triangle_vertex_index);
  void setUV(const glm::dvec3& uv, size_t triangle_vertex_index);

  NormalMode getNormalMode() const;
  void validateNormalCount(const std::string& context = {}) const;

  const std::vector<size_t>& getGroupOffsets() const { return group_offsets; }

  const size_t getTriangleCount() const { return triangles.size() / 3; }

  const size_t getVertexCount() const { return vertices.size(); }

  void setCreationKineticTime(double t) { creation_kinetic_time_ = t; }
  double getCreationKineticTime() const { return creation_kinetic_time_; }

  void setStoreMetadata(bool enabled) { store_metadata_ = enabled; }
  bool storeMetadata() const { return store_metadata_; }

  void setProfilePlanePosition(size_t vertex_index, glm::dvec2 xy);
  glm::dvec2 triangulationPlaneXY(size_t vertex_index) const;

  /// Empty string if unset; otherwise a stable fragment for filenames, e.g. `_t0.000000` or `_t3.500000`.
  std::string creationKineticTimeFilenameSuffix() const;

  // debug methods
  void checkForDegenerateTriangles() const;

 private:
  std::vector<glm::dvec3> computeVertexNormals();
};
} // namespace kinDS
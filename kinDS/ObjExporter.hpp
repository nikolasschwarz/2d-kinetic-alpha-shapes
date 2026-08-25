#pragma once

#include "VoronoiMesh.hpp"
#include "Logger.hpp"
#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <functional>
#include <limits>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <filesystem>

namespace kinDS
{

/// GPU-style per-vertex / per-face attributes written as a sidecar `.json` (EcoSysLab meshlet export).
/// Distinct from VoronoiMesh `#` comment metadata on vertices/faces.
struct ObjExportGpuAttributes
{
  std::vector<glm::dvec4> color;
  std::vector<double> boundary_distance;
  std::vector<glm::dvec2> profile_position;
  std::vector<glm::dvec2> profile_polar_coordinate;
  std::vector<double> HC;
  std::vector<double> HL;
  std::vector<double> RW;
  std::vector<double> RB;
  std::vector<double> moisture;
  std::vector<glm::dvec3> position0;
  std::vector<glm::dvec3> direction0;
  std::vector<double> root_distance;
  std::vector<double> uv_3;
  std::vector<bool> has_neighbor; // per triangle

  bool empty() const
  {
    return color.empty() && boundary_distance.empty() && profile_position.empty() && profile_polar_coordinate.empty()
      && HC.empty() && HL.empty() && RW.empty() && RB.empty() && moisture.empty() && position0.empty()
      && direction0.empty() && root_distance.empty() && uv_3.empty() && has_neighbor.empty();
  }

  size_t vertexCount() const
  {
    return std::max({ color.size(), boundary_distance.size(), profile_position.size(), profile_polar_coordinate.size(),
      HC.size(), HL.size(), RW.size(), RB.size(), moisture.size(), position0.size(), direction0.size(),
      root_distance.size(), uv_3.size() });
  }
};

struct ObjWriteOptions
{
  double uv_height_factor = 1.0;
  double uv_circum_factor = 1.0;
  /// Legacy simple JSON: only `boundary_distance` (used when @c gpu_attributes is empty).
  std::vector<float> boundary_distances_by_vertex = {};
  /// Full GPU-style JSON sidecar (preferred over @c boundary_distances_by_vertex when set and non-empty).
  std::optional<ObjExportGpuAttributes> gpu_attributes = std::nullopt;
  bool include_metadata = false;
  bool include_vertex_colors = false;
  bool alternate_section_shading = false;
  bool write_obj_groups = true;
  /// EcoSysLab-compatible OBJ: bark/interior face grouping, negate normal X, slim MTL.
  bool framework_compatible = false;
};

class ObjExporter
{
 private:
  static std::string sanitizeInlineComment(std::string value)
  {
    for (char& c : value)
    {
      if (c == '\n' || c == '\r')
      {
        c = ' ';
      }
    }
    return value;
  }

  [[noreturn]] static void throwObjExportRangeError(const std::string& message)
  {
    KINDS_ERROR("ObjExporter: " + message);
    throw std::runtime_error("ObjExporter: " + message);
  }

  static void validateGroupOffsets(const VoronoiMesh& mesh, const std::string& context)
  {
    const auto& offsets = mesh.getGroupOffsets();
    const size_t triangle_count = mesh.getTriangleCount();

    if (offsets.empty())
    {
      return;
    }

    for (size_t i = 0; i < offsets.size(); ++i)
    {
      if (offsets[i] > triangle_count)
      {
        std::ostringstream oss;
        oss << context << ": group_offsets[" << i << "]=" << offsets[i] << " exceeds triangle count " << triangle_count
            << " (group offsets must be triangle indices in [0, " << triangle_count << "])";
        throwObjExportRangeError(oss.str());
      }
      if (i > 0 && offsets[i] < offsets[i - 1])
      {
        std::ostringstream oss;
        oss << context << ": group_offsets[" << i << "]=" << offsets[i] << " is less than group_offsets[" << (i - 1)
            << "]=" << offsets[i - 1];
        throwObjExportRangeError(oss.str());
      }
    }
  }

  static void validateFaceWriteRange(const VoronoiMesh& mesh, size_t lb, size_t ub, const std::string& context)
  {
    const size_t triangle_count = mesh.getTriangleCount();
    const auto& indices = mesh.getTriangles();
    const auto& uv_indices = mesh.getUVIndices();
    const auto& normals = mesh.getNormals();
    const auto& material_ids = mesh.getMaterialIDs();
    const auto& face_metadata = mesh.getFaceMetadata();
    const size_t vertex_count = mesh.getVertices().size();
    const size_t uv_count = mesh.getUVs().size();

    if (lb > ub)
    {
      std::ostringstream oss;
      oss << context << ": face range lb=" << lb << " exceeds ub=" << ub;
      throwObjExportRangeError(oss.str());
    }

    if (ub > triangle_count)
    {
      std::ostringstream oss;
      oss << context << ": face range ub=" << ub << " exceeds triangle count " << triangle_count << " (lb=" << lb
          << ")";
      throwObjExportRangeError(oss.str());
    }

    if (indices.size() < 3 * ub)
    {
      std::ostringstream oss;
      oss << context << ": triangle index buffer size " << indices.size() << " is too small for face range [" << lb
          << ", " << ub << ") (needs at least " << (3 * ub) << " corner indices)";
      throwObjExportRangeError(oss.str());
    }

    if (!uv_indices.empty() && uv_indices.size() < 3 * ub)
    {
      std::ostringstream oss;
      oss << context << ": uv_indices size " << uv_indices.size() << " is too small for face range [" << lb << ", "
          << ub << ") (needs at least " << (3 * ub) << " entries)";
      throwObjExportRangeError(oss.str());
    }

    if (!material_ids.empty() && material_ids.size() < ub)
    {
      std::ostringstream oss;
      oss << context << ": material_ids size " << material_ids.size() << " is too small for face range [" << lb << ", "
          << ub << ")";
      throwObjExportRangeError(oss.str());
    }

    if (!face_metadata.empty() && face_metadata.size() < ub)
    {
      std::ostringstream oss;
      oss << context << ": face_metadata size " << face_metadata.size() << " is too small for face range [" << lb
          << ", " << ub << ")";
      throwObjExportRangeError(oss.str());
    }

    for (size_t triangle_index = lb; triangle_index < ub; ++triangle_index)
    {
      const size_t corner_base = 3 * triangle_index;
      for (size_t corner = 0; corner < 3; ++corner)
      {
        const size_t vertex_index = indices[corner_base + corner];
        if (vertex_index >= vertex_count)
        {
          std::ostringstream oss;
          oss << context << ": triangle " << triangle_index << " corner " << corner << " vertex index " << vertex_index
              << " is out of range [0, " << vertex_count << ")";
          throwObjExportRangeError(oss.str());
        }

        if (!uv_indices.empty())
        {
          const size_t uv_index = uv_indices[corner_base + corner];
          if (uv_index != std::numeric_limits<size_t>::max() && uv_index >= uv_count)
          {
            std::ostringstream oss;
            oss << context << ": triangle " << triangle_index << " corner " << corner << " uv index " << uv_index
                << " is out of range [0, " << uv_count << ")";
            throwObjExportRangeError(oss.str());
          }
        }

        size_t normal_index = std::numeric_limits<size_t>::max();
        if (mesh.getNormalMode() == NormalMode::PerTriangleCorner)
        {
          normal_index = corner_base + corner;
        }
        else if (mesh.getNormalMode() == NormalMode::PerVertex)
        {
          normal_index = vertex_index;
        }

        if (normal_index != std::numeric_limits<size_t>::max() && normal_index >= normals.size())
        {
          std::ostringstream oss;
          oss << context << ": triangle " << triangle_index << " corner " << corner << " normal index " << normal_index
              << " is out of range [0, " << normals.size() << ")";
          throwObjExportRangeError(oss.str());
        }
      }

      if (!material_ids.empty())
      {
        const int material_id = material_ids[triangle_index];
        if (material_id >= 0 && static_cast<size_t>(material_id) >= mesh.getMaterialNames().size())
        {
          std::ostringstream oss;
          oss << context << ": triangle " << triangle_index << " material id " << material_id << " is out of range [0, "
              << mesh.getMaterialNames().size() << ")";
          throwObjExportRangeError(oss.str());
        }
      }
    }
  }

  static std::optional<double> parseMetadataDoubleField(const std::string& metadata, const char* key)
  {
    std::string needle = std::string("\"") + key + "\":";
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
    if (value_pos >= metadata.size())
    {
      return std::nullopt;
    }

    const char* start = metadata.c_str() + value_pos;
    char* end = nullptr;
    const double value = std::strtod(start, &end);
    if (end == start || !std::isfinite(value))
    {
      return std::nullopt;
    }
    return value;
  }

  static std::optional<double> kineticTimeForVertex(const VoronoiMesh& mesh, size_t vertex_index)
  {
    if (vertex_index >= mesh.getVertexCount())
    {
      return std::nullopt;
    }

    const auto& vertex_metadata = mesh.getVertexMetadata();
    if (vertex_index < vertex_metadata.size())
    {
      if (const std::optional<double> metadata_time = parseMetadataDoubleField(vertex_metadata[vertex_index], "t"))
      {
        return metadata_time;
      }
    }

    const double stored_time = mesh.vertexKineticTime(vertex_index);
    if (std::isfinite(stored_time))
    {
      return stored_time;
    }

    return std::nullopt;
  }

  static std::optional<std::pair<double, double>> kineticTimeRangeForTriangle(
    const VoronoiMesh& mesh, size_t triangle_index)
  {
    const auto& triangles = mesh.getTriangles();
    const size_t corner_base = triangle_index * 3;
    if (corner_base + 2 >= triangles.size())
    {
      return std::nullopt;
    }

    double min_t = std::numeric_limits<double>::infinity();
    double max_t = -std::numeric_limits<double>::infinity();
    size_t valid_vertex_count = 0;
    for (size_t corner = 0; corner < 3; ++corner)
    {
      const size_t vertex_index = triangles[corner_base + corner];
      const std::optional<double> vertex_time = kineticTimeForVertex(mesh, vertex_index);
      if (!vertex_time.has_value())
      {
        continue;
      }
      min_t = std::min(min_t, vertex_time.value());
      max_t = std::max(max_t, vertex_time.value());
      ++valid_vertex_count;
    }

    if (valid_vertex_count == 0 || !std::isfinite(min_t) || !std::isfinite(max_t))
    {
      return std::nullopt;
    }
    return std::make_pair(min_t, max_t);
  }

  static size_t kineticSectionIndexFromMinTime(double min_t)
  {
    return static_cast<size_t>(std::floor(min_t));
  }

  static std::optional<size_t> kineticSectionForTriangle(const VoronoiMesh& mesh, size_t triangle_index)
  {
    const std::optional<std::pair<double, double>> time_range = kineticTimeRangeForTriangle(mesh, triangle_index);
    if (!time_range.has_value())
    {
      return std::nullopt;
    }
    return kineticSectionIndexFromMinTime(time_range->first);
  }

  static bool kineticSectionIsEven(size_t section_index) { return (section_index % 2) == 0; }

  static std::string sectionShadedMaterialName(const std::string& base_material_name, size_t section_index)
  {
    const bool even_section = kineticSectionIsEven(section_index);
    if (base_material_name == "green")
    {
      return even_section ? "green_light" : "green_dark";
    }
    if (base_material_name == "brown")
    {
      return even_section ? "brown_light" : "brown_dark";
    }
    return base_material_name;
  }

  static void writeFaces(std::ofstream& file, const VoronoiMesh& mesh, size_t lb, size_t ub, bool include_metadata,
    bool alternate_section_shading)
  {
    validateFaceWriteRange(mesh, lb, ub, "writeFaces");

    const auto& indices = mesh.getTriangles();
    const auto& uv_indices = mesh.getUVIndices();
    const auto& normals = mesh.getNormals();
    const auto& face_metadata = mesh.getFaceMetadata();

    std::string active_material_name;
    const auto& material_ids = mesh.getMaterialIDs();
    const auto& material_names = mesh.getMaterialNames();
    for (size_t i = 3 * lb; i < 3 * ub; i += 3)
    {
      std::string desired_material_name;
      if (!material_ids.empty() && (i / 3) < material_ids.size())
      {
        const int current_material_id = material_ids[i / 3];
        if (current_material_id >= 0 && static_cast<size_t>(current_material_id) < material_names.size())
        {
          desired_material_name = material_names[static_cast<size_t>(current_material_id)];
        }
      }
      if (alternate_section_shading && !desired_material_name.empty())
      {
        if (const std::optional<size_t> section_index = kineticSectionForTriangle(mesh, i / 3))
        {
          desired_material_name = sectionShadedMaterialName(desired_material_name, section_index.value());
        }
      }
      if (!desired_material_name.empty() && desired_material_name != active_material_name)
      {
        active_material_name = desired_material_name;
        file << "usemtl " << active_material_name << "\n";
      }

      file << "f";

      for (size_t j = 0; j < 3; j++)
      {
        file << " " << (indices[i + j] + 1); // 1-based index in OBJ

        size_t normal_index = std::numeric_limits<size_t>::max();
        if (mesh.getNormalMode() == NormalMode::PerTriangleCorner)
        {
          normal_index = i + j;
        }
        else if (mesh.getNormalMode() == NormalMode::PerVertex)
        {
          normal_index = indices[i + j];
        }

        bool has_uv = uv_indices[i + j] != std::numeric_limits<size_t>::max();
        bool has_normal = normal_index != std::numeric_limits<size_t>::max() && normals.size() > normal_index;

        if (has_uv || has_normal)
        {
          file << "/";
        }
        if (has_uv)
        {
          // writeMesh emits one vt record per triangle corner, after applying material-specific scaling. Therefore
          // OBJ faces must reference that emitted corner stream, not the mesh's internal shared UV-pool index.
          file << std::to_string(i + j + 1);
        }
        if (has_normal)
        {
          file << "/" << std::to_string(normal_index + 1);
        }
      }
      if (include_metadata)
      {
        const std::string metadata = (i / 3 < face_metadata.size()) ? face_metadata[i / 3] : "{}";
        file << " # " << sanitizeInlineComment(metadata);
      }
      file << "\n";
    }
  }

  static void writeFacesFrameworkCompatible(std::ofstream& file, const VoronoiMesh& mesh, size_t lb, size_t ub)
  {
    validateFaceWriteRange(mesh, lb, ub, "writeFacesFrameworkCompatible");

    const auto& indices = mesh.getTriangles();
    const auto& material_ids = mesh.getMaterialIDs();

    std::vector<size_t> bark_tris;
    std::vector<size_t> interior_tris;
    bark_tris.reserve(ub - lb);
    interior_tris.reserve(ub - lb);
    for (size_t tri = lb; tri < ub; ++tri)
    {
      const int material_id = (tri < material_ids.size()) ? material_ids[tri] : 1;
      if (material_id == 0)
      {
        bark_tris.push_back(tri);
      }
      else
      {
        interior_tris.push_back(tri);
      }
    }

    auto write_tri_faces = [&](const std::vector<size_t>& tris)
    {
      for (size_t tri : tris)
      {
        const size_t corner_base = 3 * tri;
        file << "f";
        for (size_t j = 0; j < 3; ++j)
        {
          const size_t corner = corner_base + j;
          // vt/vn streams are emitted once per triangle corner in triangle order.
          file << " " << (indices[corner] + 1) << "/" << (corner + 1) << "/" << (corner + 1);
        }
        file << "\n";
      }
    };

    file << "usemtl bark\n";
    write_tri_faces(bark_tris);
    file << "usemtl interior\n";
    write_tri_faces(interior_tris);
  }

  static void writeMtl(const std::filesystem::path& mtl_path)
  {
    std::ofstream file(mtl_path);
    if (!file.is_open())
    {
      throw std::runtime_error("Failed to open MTL file");
    }

    // TODO: Define proper materials or perhaps pass them as arguments

    // Brown material (boundary interval faces).
    file << "newmtl brown\n";
    file << "Ka 0.2 0.1 0.05\n";
    file << "Kd 0.4 0.25 0.1\n";
    file << "Ks 0.0 0.0 0.0\n";
    file << "d 1.0\n\n";

    // Green material (regular / interior Voronoi-edge strip meshlets).
    file << "newmtl green\n";
    file << "Ka 0.6 0.6 0.1\n";
    file << "Kd 0.9 0.85 0.2\n";
    file << "Ks 0.0 0.0 0.0\n";
    file << "d 1.0\n\n";

    // Pending-split mixed-branch radius fallback (reserved; non-shift radius cell fans use brown).
    file << "newmtl light_blue\n";
    file << "Ka 0.15 0.35 0.5\n";
    file << "Kd 0.35 0.75 1.0\n";
    file << "Ks 0.0 0.0 0.0\n";
    file << "d 1.0\n\n";

    // Section-alternating shades (even sections = light, odd sections = dark).
    file << "newmtl green_light\n";
    file << "Ka 0.6 0.6 0.1\n";
    file << "Kd 0.9 0.85 0.2\n";
    file << "Ks 0.0 0.0 0.0\n";
    file << "d 1.0\n\n";

    file << "newmtl green_dark\n";
    file << "Ka 0.35 0.35 0.05\n";
    file << "Kd 0.55 0.5 0.1\n";
    file << "Ks 0.0 0.0 0.0\n";
    file << "d 1.0\n\n";

    file << "newmtl brown_light\n";
    file << "Ka 0.2 0.1 0.05\n";
    file << "Kd 0.4 0.25 0.1\n";
    file << "Ks 0.0 0.0 0.0\n";
    file << "d 1.0\n\n";

    file << "newmtl brown_dark\n";
    file << "Ka 0.1 0.05 0.02\n";
    file << "Kd 0.22 0.12 0.05\n";
    file << "Ks 0.0 0.0 0.0\n";
    file << "d 1.0\n\n";

    // Component boundary shell (bark / interior).
    file << "newmtl bark\n";
    file << "Ka 0.3 0.2 0.1\n";
    file << "Kd 0.5 0.35 0.15\n";
    file << "Ks 0.0 0.0 0.0\n";
    file << "d 1.0\n\n";

    file << "newmtl interior\n";
    file << "Ka 0.8 0.8 0.8\n";
    file << "Kd 0.8 0.8 0.8\n";
    file << "Ks 0.0 0.0 0.0\n";
    file << "d 1.0\n\n";

    // Validation error highlight (inconsistent keyed vertex sources).
    file << "newmtl red\n";
    file << "Ka 0.3 0.0 0.0\n";
    file << "Kd 0.95 0.1 0.1\n";
    file << "Ks 0.0 0.0 0.0\n";
    file << "d 1.0\n\n";

    // Validation error highlight for triangles whose original material was green.
    file << "newmtl neon_yellow\n";
    file << "Ka 0.6 0.6 0.0\n";
    file << "Kd 1.0 1.0 0.0\n";
    file << "Ke 0.2 0.2 0.0\n";
    file << "Ks 0.0 0.0 0.0\n";
    file << "d 1.0\n";

    file.close();
  }

  /// Slim MTL matching EcoSysLab meshlet export (bark + interior only).
  static void writeFrameworkMtl(const std::filesystem::path& mtl_path)
  {
    std::ofstream file(mtl_path);
    if (!file.is_open())
    {
      throw std::runtime_error("Failed to open MTL file");
    }

    file << "newmtl bark\n";
    file << "Ka 0.2 0.1 0.05\n";
    file << "Kd 0.4 0.25 0.1\n";
    file << "Ks 0.0 0.0 0.0\n";
    file << "d 1.0\n\n";

    file << "newmtl interior\n";
    file << "Ka 0.8 0.8 0.8\n";
    file << "Kd 0.8 0.8 0.8\n";
    file << "Ks 0.0 0.0 0.0\n";
    file << "d 1.0\n";

    file.close();
  }

  static std::string jsonVec2(const glm::dvec2& v)
  {
    return "[" + std::to_string(v.x) + ", " + std::to_string(v.y) + "]";
  }

  static std::string jsonVec3(const glm::dvec3& v)
  {
    return "[" + std::to_string(v.x) + ", " + std::to_string(v.y) + ", " + std::to_string(v.z) + "]";
  }

  static std::string jsonVec4(const glm::dvec4& v)
  {
    return "[" + std::to_string(v.x) + ", " + std::to_string(v.y) + ", " + std::to_string(v.z) + ", "
      + std::to_string(v.w) + "]";
  }

 public:
  /// Legacy JSON with a single `boundary_distance` array.
  static void writeJson(const std::filesystem::path& json_path, const std::vector<float>& boundary_distances_by_vertex)
  {
    ObjExportGpuAttributes attrs;
    attrs.boundary_distance.assign(boundary_distances_by_vertex.begin(), boundary_distances_by_vertex.end());
    writeJson(json_path, attrs);
  }

  /// Framework-style GPU attribute JSON (EcoSysLab meshlet export schema).
  static void writeJson(const std::filesystem::path& json_path, const ObjExportGpuAttributes& attrs)
  {
    std::ofstream file(json_path);
    if (!file.is_open())
    {
      throw std::runtime_error("Failed to open JSON file");
    }

    struct Field
    {
      std::string key;
      size_t size = 0;
      std::function<std::string(size_t)> get_value;
    };
    std::vector<Field> fields;
    auto add_field = [&](const std::string& key, size_t size, std::function<std::string(size_t)> get_value)
    {
      if (size == 0)
      {
        return;
      }
      fields.push_back(Field { key, size, std::move(get_value) });
    };

    add_field("color", attrs.color.size(), [&](size_t i) { return jsonVec4(attrs.color[i]); });
    add_field("boundary_distance", attrs.boundary_distance.size(),
      [&](size_t i) { return std::to_string(attrs.boundary_distance[i]); });
    add_field("profile_position", attrs.profile_position.size(),
      [&](size_t i) { return jsonVec2(attrs.profile_position[i]); });
    add_field("profile_polar_coordinate", attrs.profile_polar_coordinate.size(),
      [&](size_t i) { return jsonVec2(attrs.profile_polar_coordinate[i]); });
    add_field("HC", attrs.HC.size(), [&](size_t i) { return std::to_string(attrs.HC[i]); });
    add_field("HL", attrs.HL.size(), [&](size_t i) { return std::to_string(attrs.HL[i]); });
    add_field("RW", attrs.RW.size(), [&](size_t i) { return std::to_string(attrs.RW[i]); });
    add_field("RB", attrs.RB.size(), [&](size_t i) { return std::to_string(attrs.RB[i]); });
    add_field("moisture", attrs.moisture.size(), [&](size_t i) { return std::to_string(attrs.moisture[i]); });
    add_field("position0", attrs.position0.size(), [&](size_t i) { return jsonVec3(attrs.position0[i]); });
    add_field("direction0", attrs.direction0.size(), [&](size_t i) { return jsonVec3(attrs.direction0[i]); });
    add_field("root_distance", attrs.root_distance.size(),
      [&](size_t i) { return std::to_string(attrs.root_distance[i]); });
    add_field("uv_3", attrs.uv_3.size(), [&](size_t i) { return std::to_string(attrs.uv_3[i]); });
    add_field("has_neighbor", attrs.has_neighbor.size(),
      [&](size_t i) { return attrs.has_neighbor[i] ? "true" : "false"; });

    file << "{\n";
    for (size_t field_index = 0; field_index < fields.size(); ++field_index)
    {
      const Field& field = fields[field_index];
      file << "  \"" << field.key << "\": [\n";
      for (size_t i = 0; i < field.size; ++i)
      {
        file << "    " << field.get_value(i);
        if (i + 1 < field.size)
        {
          file << ",";
        }
        file << "\n";
      }
      file << "  ]";
      if (field_index + 1 < fields.size())
      {
        file << ",";
      }
      file << "\n";
    }
    file << "}\n";
  }

  static void writeMesh(const VoronoiMesh& mesh, const std::filesystem::path& obj_path, const ObjWriteOptions& options)
  {
    mesh.validateNormalCount("ObjExporter::writeMesh(" + obj_path.string() + ")");
    mesh.validateUVLayout("ObjExporter::writeMesh(" + obj_path.string() + ")");

    std::ofstream file(obj_path);
    if (!file.is_open())
    {
      throw std::runtime_error("Failed to open file for writing: " + obj_path.filename().string());
    }
    std::filesystem::path mtl_path = obj_path;
    mtl_path.replace_extension(".mtl");
    if (options.framework_compatible)
    {
      writeFrameworkMtl(mtl_path);
    }
    else
    {
      writeMtl(mtl_path);
    }

    if (options.gpu_attributes.has_value() && !options.gpu_attributes->empty())
    {
      std::filesystem::path json_path = obj_path;
      json_path.replace_extension(".json");
      writeJson(json_path, options.gpu_attributes.value());
    }
    else if (!options.boundary_distances_by_vertex.empty())
    {
      std::filesystem::path json_path = obj_path;
      json_path.replace_extension(".json");
      writeJson(json_path, options.boundary_distances_by_vertex);
    }

    file << "# Exported by kinDS ObjExporter\n";
    file << "mtllib " << mtl_path.filename() << "\n\n";

    file << "# Vertices\n";
    const auto& vertex_metadata = mesh.getVertexMetadata();
    const auto& vertex_colors = mesh.getVertexColors();
    for (size_t i = 0; i < mesh.getVertices().size(); ++i)
    {
      const auto& vertex = mesh.getVertices()[i];
      file << "v " << vertex[0] << " " << vertex[1] << " " << vertex[2];
      if (options.include_vertex_colors && i < vertex_colors.size())
      {
        const auto& c = vertex_colors[i];
        file << " " << c[0] << " " << c[1] << " " << c[2];
      }
      if (options.include_metadata)
      {
        const std::string metadata = (i < vertex_metadata.size()) ? vertex_metadata[i] : "{}";
        file << " # " << sanitizeInlineComment(metadata);
      }
      file << "\n";
    }

    if (options.framework_compatible)
    {
      // Interleaved vt/vn per triangle corner (EcoSysLab style), with negated normal X.
      file << "# Texture coordinates and normals\n";
      for (size_t i = 0; i < mesh.getTriangleCount(); ++i)
      {
        int material = (i < mesh.getMaterialIDs().size()) ? mesh.getMaterialIDs()[i] : -1;
        for (size_t j = 0; j < 3; ++j)
        {
          const size_t corner_index = 3 * i + j;
          auto uv = mesh.hasValidUVIndex(corner_index) ? mesh.getUV(corner_index) : glm::dvec3(0.0);
          if (material == 0)
          {
            uv[0] *= options.uv_circum_factor;
            uv[1] *= options.uv_height_factor;
          }
          else if (material != -1)
          {
            uv[2] *= options.uv_height_factor;
          }

          glm::dvec3 normal(0.0);
          if (corner_index < mesh.getNormals().size())
          {
            normal = mesh.getNormals()[corner_index];
          }
          file << "vt " << uv[0] << " " << uv[1] << " " << uv[2] << "\n";
          file << "vn " << (-normal[0]) << " " << normal[1] << " " << normal[2] << "\n";
        }
      }
    }
    else
    {
      file << "# Normals\n";
      for (const auto& normal : mesh.getNormals())
      {
        file << "vn " << normal[0] << " " << normal[1] << " " << normal[2] << "\n";
      }

      file << "# UVs\n";
      for (size_t i = 0; i < mesh.getTriangleCount(); i++)
      {
        int material = -1;
        if (i < mesh.getMaterialIDs().size())
        {
          material = mesh.getMaterialIDs()[i];
        }
        for (size_t j = 0; j < 3; j++)
        {
          const size_t corner_index = 3 * i + j;
          auto uv = mesh.hasValidUVIndex(corner_index) ? mesh.getUV(corner_index) : glm::dvec3(0.0);
          if (material != -1)
          {
            if (material == 0)
            {
              uv[0] *= options.uv_circum_factor;
              uv[1] *= options.uv_height_factor;
            }
            else
            {
              uv[2] *= options.uv_height_factor;
            }
          }
          file << "vt " << uv[0] << " " << uv[1] << " " << uv[2] << "\n";
        }
      }
    }

    file << "# Faces\n";
    size_t group_count = mesh.getGroupOffsets().size();
    const std::string mesh_context = "writeMesh(" + obj_path.string() + ")";

    auto write_face_range = [&](size_t lb, size_t ub)
    {
      if (options.framework_compatible)
      {
        writeFacesFrameworkCompatible(file, mesh, lb, ub);
      }
      else
      {
        writeFaces(file, mesh, lb, ub, options.include_metadata, options.alternate_section_shading);
      }
    };

    if (options.write_obj_groups && group_count > 0)
    {
      validateGroupOffsets(mesh, mesh_context);

      for (size_t group_index = 0; group_index < group_count - 1; group_index++)
      {
        const auto& names = mesh.getGroupNames();
        const std::string group_name = (group_index < names.size() && !names[group_index].empty())
          ? names[group_index]
          : ("group_" + std::to_string(group_index));
        file << "o " << group_name << "\n";
        size_t lb = mesh.getGroupOffsets()[group_index];
        size_t ub = mesh.getGroupOffsets()[group_index + 1];
        validateFaceWriteRange(mesh, lb, ub, mesh_context + " group " + std::to_string(group_index));
        write_face_range(lb, ub);
      }

      {
        const auto& names = mesh.getGroupNames();
        const size_t last_index = group_count - 1;
        const std::string group_name = (last_index < names.size() && !names[last_index].empty())
          ? names[last_index]
          : ("group_" + std::to_string(last_index));
        file << "o " << group_name << "\n";
      }
      size_t lb = mesh.getGroupOffsets().back();
      size_t ub = mesh.getTriangles().size() / 3;
      validateFaceWriteRange(mesh, lb, ub, mesh_context + " last group");
      write_face_range(lb, ub);
    }
    else
    {
      const size_t ub = mesh.getTriangles().size() / 3;
      validateFaceWriteRange(mesh, 0, ub, mesh_context);
      write_face_range(0, ub);
    }

    file.close();
  }

  /// Backward-compatible overload used by TreeMesher / CLI.
  static void writeMesh(const VoronoiMesh& mesh, const std::filesystem::path& obj_path, double uv_height_factor = 1.0,
    double uv_circum_factor = 1.0, const std::vector<float>& boundary_distances_by_vertex = {},
    bool include_metadata = false, bool include_vertex_colors = false, bool alternate_section_shading = false,
    bool write_obj_groups = true)
  {
    ObjWriteOptions options;
    options.uv_height_factor = uv_height_factor;
    options.uv_circum_factor = uv_circum_factor;
    options.boundary_distances_by_vertex = boundary_distances_by_vertex;
    options.include_metadata = include_metadata;
    options.include_vertex_colors = include_vertex_colors;
    options.alternate_section_shading = alternate_section_shading;
    options.write_obj_groups = write_obj_groups;
    options.framework_compatible = false;
    writeMesh(mesh, obj_path, options);
  }

  static size_t parseObjFaceIndexToken(const std::string& token)
  {
    const size_t slash = token.find('/');
    const std::string index_token = slash == std::string::npos ? token : token.substr(0, slash);
    if (index_token.empty())
    {
      throw std::runtime_error("ObjExporter::readMesh: empty face index token.");
    }
    const long index = std::stol(index_token);
    if (index == 0)
    {
      throw std::runtime_error("ObjExporter::readMesh: face index 0 is invalid in OBJ (indices are 1-based).");
    }
    return static_cast<size_t>(index > 0 ? index - 1 : index);
  }

  static VoronoiMesh readMesh(const std::filesystem::path& obj_path)
  {
    std::ifstream file(obj_path);
    if (!file.is_open())
    {
      throw std::runtime_error("ObjExporter::readMesh: failed to open file: " + obj_path.string());
    }

    std::vector<glm::dvec3> positions;
    positions.reserve(1024);

    VoronoiMesh mesh({ "boundary" }, NoNormals);

    std::string line;
    while (std::getline(file, line))
    {
      if (line.empty() || line[0] == '#')
      {
        continue;
      }

      std::istringstream iss(line);
      std::string tag;
      iss >> tag;
      if (tag == "v")
      {
        double x = 0.0;
        double y = 0.0;
        double z = 0.0;
        iss >> x >> y >> z;
        positions.emplace_back(x, y, z);
      }
      else if (tag == "f")
      {
        std::vector<size_t> face_vertices;
        std::string token;
        while (iss >> token)
        {
          face_vertices.push_back(parseObjFaceIndexToken(token));
        }

        if (face_vertices.size() < 3)
        {
          continue;
        }

        for (size_t i = 1; i + 1 < face_vertices.size(); ++i)
        {
          const size_t i0 = face_vertices[0];
          const size_t i1 = face_vertices[i];
          const size_t i2 = face_vertices[i + 1];
          if (i0 >= positions.size() || i1 >= positions.size() || i2 >= positions.size())
          {
            throw std::runtime_error("ObjExporter::readMesh: face references out-of-range vertex index.");
          }
          const size_t v0 = mesh.addVertex(positions[i0]);
          const size_t v1 = mesh.addVertex(positions[i1]);
          const size_t v2 = mesh.addVertex(positions[i2]);
          mesh.addTriangle(v0, v1, v2, 0);
        }
      }
    }

    if (mesh.getTriangleCount() == 0)
    {
      throw std::runtime_error("ObjExporter::readMesh: no triangles found in " + obj_path.string());
    }

    mesh.mergeDuplicateVertices(1e-6);
    if (mesh.orientFacesAwayFromCentroid())
    {
      KINDS_INFO("Loaded mesh was flipped to orient faces away from centroid: " << obj_path.string());
    }
    mesh.computeNormals(PerTriangleCorner);
    return mesh;
  }
};
} // namespace kinDS
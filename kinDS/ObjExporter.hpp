#pragma once

#include "VoronoiMesh.hpp"
#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <limits>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>

#include <filesystem>

namespace kinDS
{
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

    // Green material (regular meshes emitted by radius events).
    file << "newmtl green\n";
    file << "Ka 0.6 0.6 0.1\n";
    file << "Kd 0.9 0.85 0.2\n";
    file << "Ks 0.0 0.0 0.0\n";
    file << "d 1.0\n\n";

    // Pending-split mixed-branch radius fallback.
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

 public:
  static void writeJson(const std::filesystem::path& json_path, const std::vector<float>& boundary_distances_by_vertex)
  {
    std::ofstream file(json_path);

    if (!file.is_open())
    {
      throw std::runtime_error("Failed to open JSON file");
    }

    auto write_values = [](std::ofstream& file, const std::string& key,
                          const std::function<std::string(size_t index)>& get_value, size_t size, bool last = false)
    {
      file << "  \"" << key << "\": [\n";
      for (size_t i = 0; i < size; ++i)
      {
        file << "    " << get_value(i);
        if (i < size - 1)
        {
          file << ",";
        }
        file << "\n";
      }

      if (!last)
      {
        file << "  ],\n";
      }
      else
      {
        file << "  ]\n";
      }
    };

    file << "{\n";

    write_values(
      file, "boundary_distance", [&](size_t index) { return std::to_string(boundary_distances_by_vertex[index]); },
      boundary_distances_by_vertex.size(), true);

    file << "}\n";
  }
  static void writeMesh(const VoronoiMesh& mesh, const std::filesystem::path& obj_path, double uv_height_factor = 1.0,
    double uv_circum_factor = 1.0, const std::vector<float>& boundary_distances_by_vertex = {},
    bool include_metadata = false, bool include_vertex_colors = false, bool alternate_section_shading = false)
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
    writeMtl(mtl_path);

    if (!boundary_distances_by_vertex.empty())
    {
      std::filesystem::path json_path = obj_path;
      json_path.replace_extension(".json");
      writeJson(json_path, boundary_distances_by_vertex);
    }

    // Write some metadata
    file << "# Exported by kinDS ObjExporter\n";
    file << "mtllib " << mtl_path.filename() << "\n\n";

    // Write vertices
    file << "# Vertices\n";
    const auto& vertex_metadata = mesh.getVertexMetadata();
    const auto& vertex_colors = mesh.getVertexColors();
    for (size_t i = 0; i < mesh.getVertices().size(); ++i)
    {
      const auto& vertex = mesh.getVertices()[i];
      file << "v " << vertex[0] << " " << vertex[1] << " " << vertex[2];
      if (include_vertex_colors && i < vertex_colors.size())
      {
        const auto& c = vertex_colors[i];
        file << " " << c[0] << " " << c[1] << " " << c[2];
      }
      if (include_metadata)
      {
        const std::string metadata = (i < vertex_metadata.size()) ? vertex_metadata[i] : "{}";
        file << " # " << sanitizeInlineComment(metadata);
      }
      file << "\n";
    }

    // Write normals
    file << "# Normals\n";
    for (const auto& normal : mesh.getNormals())
    {
      file << "vn " << normal[0] << " " << normal[1] << " " << normal[2] << "\n";
    }

    // Write UVs
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
            uv[0] *= uv_circum_factor;
            uv[1] *= uv_height_factor;
          }
          else
          {
            uv[2] *= uv_height_factor;
          }
        }
        file << "vt " << uv[0] << " " << uv[1] << " " << uv[2] << "\n";
      }
    }

    // Write faces
    file << "# Faces\n";
    size_t group_count = mesh.getGroupOffsets().size();
    const std::string mesh_context = "writeMesh(" + obj_path.string() + ")";

    if (group_count > 0)
    {
      validateGroupOffsets(mesh, mesh_context);

      for (size_t group_index = 0; group_index < group_count - 1; group_index++)
      {
        file << "o group_" << group_index << "\n";
        size_t lb = mesh.getGroupOffsets()[group_index];
        size_t ub = mesh.getGroupOffsets()[group_index + 1];
        validateFaceWriteRange(mesh, lb, ub, mesh_context + " group " + std::to_string(group_index));
        writeFaces(file, mesh, lb, ub, include_metadata, alternate_section_shading);
      }

      // Write the last group

      file << "o group_" << (group_count - 1) << "\n";
      size_t lb = mesh.getGroupOffsets().back();
      size_t ub = mesh.getTriangles().size() / 3;
      validateFaceWriteRange(mesh, lb, ub, mesh_context + " last group");
      writeFaces(file, mesh, lb, ub, include_metadata, alternate_section_shading);
    }
    else
    {
      // No groups defined, write all faces
      const size_t ub = mesh.getTriangles().size() / 3;
      validateFaceWriteRange(mesh, 0, ub, mesh_context);
      writeFaces(file, mesh, 0, ub, include_metadata, alternate_section_shading);
    }

    file.close();
  }
};
} // namespace kinDS
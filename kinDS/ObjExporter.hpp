#pragma once

#include "VoronoiMesh.hpp"
#include <fstream>
#include <limits>
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

  static void writeFaces(std::ofstream& file, const VoronoiMesh& mesh, size_t lb, size_t ub, bool include_metadata)
  {
    const auto& indices = mesh.getTriangles();
    const auto& uv_indices = mesh.getUVIndices();
    const auto& normals = mesh.getNormals();
    const auto& face_metadata = mesh.getFaceMetadata();

    std::string active_material_name;
    const auto& material_ids = mesh.getMaterialIDs();
    auto materialFromMetadata = [](const std::string& metadata) -> std::string
    {
      // Boundary-interval meshes carry segment_action in metadata.
      if (metadata.find("\"segment_action\"") != std::string::npos)
      {
        return "brown";
      }
      // Radius-event regular meshes: metadata uses lowercase event tag from @ref boundaryEventTypeToString.
      if ((metadata.find("\"event_type\":\"radius_event\"") != std::string::npos
            || metadata.find("\"event_type\":\"Radius\"") != std::string::npos)
        && metadata.find("\"mesh_type\":\"regular\"") != std::string::npos)
      {
        return "yellow";
      }
      return {};
    };
    for (size_t i = 3 * lb; i < 3 * ub; i += 3)
    {
      std::string desired_material_name;
      if (!material_ids.empty() && (i / 3) < material_ids.size())
      {
        const int current_material_id = material_ids[i / 3];
        if (current_material_id >= 0 && static_cast<size_t>(current_material_id) < mesh.getMaterialNames().size())
        {
          desired_material_name = mesh.getMaterialNames()[static_cast<size_t>(current_material_id)];
        }
      }
      if (desired_material_name.empty() && include_metadata)
      {
        const std::string metadata = (i / 3 < face_metadata.size()) ? face_metadata[i / 3] : "{}";
        desired_material_name = materialFromMetadata(metadata);
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

        size_t normal_index;
        if (mesh.getNormalMode() == NormalMode::PerTriangleCorner)
        {
          normal_index = i + j;
        }
        else
        {
          normal_index = indices[i + j];
        }

        bool has_uv = uv_indices[i + j] != std::numeric_limits<size_t>::max();
        bool has_normal = normals.size() > normal_index;

        if (has_uv || has_normal)
        {
          file << "/";
        }
        if (has_uv)
        {
          file << std::to_string(uv_indices[i + j] + 1);
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

    // Yellow material (regular meshes emitted by radius events).
    file << "newmtl yellow\n";
    file << "Ka 0.6 0.6 0.1\n";
    file << "Kd 0.9 0.85 0.2\n";
    file << "Ks 0.0 0.0 0.0\n";
    file << "d 1.0\n\n";

    // Legacy material kept for existing exports.
    file << "newmtl interior\n";
    file << "Ka 0.8 0.8 0.8\n";
    file << "Kd 0.8 0.8 0.8\n";
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
    bool include_metadata = false, bool include_vertex_colors = false)
  {
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
        auto uv = mesh.getUV(3 * i + j);

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

    if (group_count > 0)
    {
      for (size_t group_index = 0; group_index < group_count - 1; group_index++)
      {
        file << "o group_" << group_index << "\n";
        size_t lb = mesh.getGroupOffsets()[group_index];
        size_t ub = mesh.getGroupOffsets()[group_index + 1];
        writeFaces(file, mesh, lb, ub, include_metadata);
      }

      // Write the last group

      file << "o group_" << (group_count - 1) << "\n";
      size_t lb = mesh.getGroupOffsets().back();
      size_t ub = mesh.getTriangles().size() / 3;
      writeFaces(file, mesh, lb, ub, include_metadata);
    }
    else
    {
      // No groups defined, write all faces
      writeFaces(file, mesh, 0, mesh.getTriangles().size() / 3, include_metadata);
    }

    file.close();
  }
};
} // namespace kinDS
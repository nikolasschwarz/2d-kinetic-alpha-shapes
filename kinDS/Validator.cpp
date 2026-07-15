#include "Validator.hpp"

#include <algorithm>
#include <cmath>
#include <cctype>
#include <cstring>
#include <ctime>
#include <fstream>
#include <glm/gtx/norm.hpp>
#include <iomanip>
#include <iostream>
#include <limits>
#include <optional>
#include <sstream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

using namespace kinDS;

namespace
{
std::string validation_log_path = Validator::defaultLogFilePath();
std::ofstream validation_log_file;

void ensureValidationLogOpen()
{
  if (validation_log_file.is_open())
  {
    return;
  }

  validation_log_file.open(validation_log_path, std::ios::app);
  if (!validation_log_file.is_open())
  {
    std::cerr << "Warning: Could not open validation log file: " << validation_log_path << std::endl;
  }
}

void writeValidationLog(const char* level, const std::string& message)
{
  ensureValidationLogOpen();
  if (!validation_log_file.is_open())
  {
    std::cerr << "[validation][" << level << "] " << message << std::endl;
    return;
  }

  const std::time_t now = std::time(nullptr);
  std::tm timeinfo {};
#if defined(_WIN32)
  localtime_s(&timeinfo, &now);
#else
  localtime_r(&now, &timeinfo);
#endif
  char timestamp[20];
  std::strftime(timestamp, sizeof(timestamp), "%Y-%m-%d %H:%M:%S", &timeinfo);

  validation_log_file << '[' << timestamp << "] [" << level << "] " << message << '\n';
  validation_log_file.flush();
}

void beginValidationSession(const char* scope)
{
  writeValidationLog("INFO", std::string("=== Mesh vertex source validation: ") + scope + " ===");
}

std::string jsonStringLiteral(const std::string& value)
{
  std::string out = "\"";
  for (char ch : value)
  {
    switch (ch)
    {
    case '\\':
      out += "\\\\";
      break;
    case '"':
      out += "\\\"";
      break;
    case '\n':
      out += "\\n";
      break;
    case '\r':
      out += "\\r";
      break;
    case '\t':
      out += "\\t";
      break;
    default:
      out += ch;
      break;
    }
  }
  out += '"';
  return out;
}

std::string numberLiteral(double value)
{
  std::ostringstream o;
  o << std::setprecision(17) << std::showpoint << value;
  return o.str();
}

std::optional<std::string> metadataStringField(const std::string& metadata, const char* key)
{
  const std::string needle = jsonStringLiteral(key) + ":";
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

std::optional<size_t> metadataSizeField(const std::string& metadata, const char* key)
{
  const std::string needle = jsonStringLiteral(key) + ":";
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
  if (value_pos >= metadata.size() || !std::isdigit(static_cast<unsigned char>(metadata[value_pos])))
  {
    return std::nullopt;
  }
  size_t value = 0;
  for (; value_pos < metadata.size() && std::isdigit(static_cast<unsigned char>(metadata[value_pos])); ++value_pos)
  {
    value = value * 10 + static_cast<size_t>(metadata[value_pos] - '0');
  }
  return value;
}

std::optional<double> metadataDoubleField(const std::string& metadata, const char* key)
{
  const std::string needle = jsonStringLiteral(key) + ":";
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

enum class MeshVertexSourceKind : uint8_t
{
  Intersection = 0,
  VoronoiVertex = 1,
  Site = 2,
};

struct MeshVertexSourceKey
{
  MeshVertexSourceKind kind = MeshVertexSourceKind::Intersection;
  size_t primary_id = static_cast<size_t>(-1);
  size_t secondary_id = static_cast<size_t>(-1);
  uint64_t kinetic_time_bits = 0;

  bool operator==(const MeshVertexSourceKey& other) const noexcept
  {
    return kind == other.kind && primary_id == other.primary_id && secondary_id == other.secondary_id
      && kinetic_time_bits == other.kinetic_time_bits;
  }
};

struct MeshVertexSourceKeyHash
{
  size_t operator()(const MeshVertexSourceKey& key) const noexcept
  {
    size_t h = std::hash<uint8_t> {}(static_cast<uint8_t>(key.kind));
    h ^= std::hash<size_t> {}(key.primary_id) + 0x9e3779b97f4a7c15ull + (h << 6) + (h >> 2);
    h ^= std::hash<size_t> {}(key.secondary_id) + 0x9e3779b97f4a7c15ull + (h << 6) + (h >> 2);
    h ^= std::hash<uint64_t> {}(key.kinetic_time_bits) + 0x9e3779b97f4a7c15ull + (h << 6) + (h >> 2);
    return h;
  }
};

struct MeshVertexSourceRecord
{
  size_t meshlet_index = static_cast<size_t>(-1);
  std::string meshlet_label;
  size_t vertex_index = static_cast<size_t>(-1);
  glm::dvec3 world_position {};
};

uint64_t kineticTimeKeyBits(double t)
{
  uint64_t bits = 0;
  static_assert(sizeof(bits) == sizeof(t));
  std::memcpy(&bits, &t, sizeof(bits));
  return bits;
}

size_t profileXYKeyBits(double x, double y)
{
  uint64_t x_bits = 0;
  uint64_t y_bits = 0;
  std::memcpy(&x_bits, &x, sizeof(x));
  std::memcpy(&y_bits, &y, sizeof(y));
  uint64_t combined = x_bits;
  combined ^= y_bits + 0x9e3779b97f4a7c15ull + (combined << 6) + (combined >> 2);
  return static_cast<size_t>(combined);
}

std::optional<MeshVertexSourceKey> meshVertexSourceKeyFromMetadata(const std::string& metadata)
{
  if (metadata.empty() || metadata == "{}")
  {
    return std::nullopt;
  }
  if (metadata.find("\"intersection_flexible_placeholder\":true") != std::string::npos)
  {
    return std::nullopt;
  }

  const std::optional<double> kinetic_time = metadataDoubleField(metadata, "t");
  if (!kinetic_time.has_value())
  {
    return std::nullopt;
  }

  const std::optional<std::string> source = metadataStringField(metadata, "source");
  const std::string source_value = source.value_or("");

  if (source_value == "intersection"
    || metadataSizeField(metadata, "delaunay_edge_id").has_value()
    || metadataSizeField(metadata, "crossing_delaunay_edge_id").has_value())
  {
    std::optional<size_t> delaunay_edge_id = metadataSizeField(metadata, "conceptual_delaunay_edge_id");
    if (!delaunay_edge_id.has_value())
    {
      delaunay_edge_id = metadataSizeField(metadata, "delaunay_edge_id");
    }
    if (!delaunay_edge_id.has_value())
    {
      delaunay_edge_id = metadataSizeField(metadata, "crossing_delaunay_edge_id");
    }

    std::optional<size_t> voronoi_edge_id = metadataSizeField(metadata, "conceptual_voronoi_edge_id");
    if (!voronoi_edge_id.has_value())
    {
      voronoi_edge_id = metadataSizeField(metadata, "voronoi_edge_id");
    }
    if (!voronoi_edge_id.has_value())
    {
      voronoi_edge_id = metadataSizeField(metadata, "crossing_voronoi_edge_id");
    }

    if (!delaunay_edge_id.has_value() || !voronoi_edge_id.has_value())
    {
      return std::nullopt;
    }

    MeshVertexSourceKey key;
    key.kind = MeshVertexSourceKind::Intersection;
    key.primary_id = delaunay_edge_id.value();
    key.secondary_id = voronoi_edge_id.value();
    key.kinetic_time_bits = kineticTimeKeyBits(kinetic_time.value());
    return key;
  }

  if (source_value == "Voronoi vertex")
  {
    const std::optional<size_t> voronoi_vertex_id = metadataSizeField(metadata, "voronoi_vertex_id");
    if (!voronoi_vertex_id.has_value())
    {
      return std::nullopt;
    }

    MeshVertexSourceKey key;
    key.kind = MeshVertexSourceKind::VoronoiVertex;
    key.primary_id = voronoi_vertex_id.value();
    key.secondary_id = static_cast<size_t>(-1);
    key.kinetic_time_bits = kineticTimeKeyBits(kinetic_time.value());
    return key;
  }

  if (source_value == "site")
  {
    std::optional<size_t> strand_id = metadataSizeField(metadata, "strand_id");
    if (!strand_id.has_value())
    {
      strand_id = metadataSizeField(metadata, "strand_cell_id");
    }
    if (!strand_id.has_value())
    {
      return std::nullopt;
    }

    const std::optional<double> profile_x = metadataDoubleField(metadata, "x");
    const std::optional<double> profile_y = metadataDoubleField(metadata, "y");
    if (!profile_x.has_value() || !profile_y.has_value())
    {
      return std::nullopt;
    }

    MeshVertexSourceKey key;
    key.kind = MeshVertexSourceKind::Site;
    key.primary_id = strand_id.value();
    key.secondary_id = profileXYKeyBits(profile_x.value(), profile_y.value());
    key.kinetic_time_bits = kineticTimeKeyBits(kinetic_time.value());
    return key;
  }

  return std::nullopt;
}

std::string describeMeshVertexSourceKey(const MeshVertexSourceKey& key)
{
  double t = 0.0;
  std::memcpy(&t, &key.kinetic_time_bits, sizeof(t));
  std::ostringstream oss;
  if (key.kind == MeshVertexSourceKind::Intersection)
  {
    oss << "intersection d=" << key.primary_id << " v=" << key.secondary_id << " t=" << numberLiteral(t);
  }
  else if (key.kind == MeshVertexSourceKind::VoronoiVertex)
  {
    oss << "voronoi_vertex id=" << key.primary_id << " t=" << numberLiteral(t);
  }
  else
  {
    oss << "site strand=" << key.primary_id << " profile_key=" << key.secondary_id << " t=" << numberLiteral(t);
  }
  return oss.str();
}

MeshVertexSourceValidationResult validateMeshVertexSourceRecords(
  const std::unordered_map<MeshVertexSourceKey, std::vector<MeshVertexSourceRecord>, MeshVertexSourceKeyHash>& buckets,
  double position_tolerance)
{
  MeshVertexSourceValidationResult result;
  result.unique_source_count = buckets.size();

  for (const auto& [key, records] : buckets)
  {
    result.keyed_vertex_count += records.size();
    if (records.size() < 2)
    {
      continue;
    }

    const glm::dvec3& reference = records.front().world_position;
    if (!std::isfinite(reference.x) || !std::isfinite(reference.y) || !std::isfinite(reference.z))
    {
      ++result.inconsistent_source_count;
      writeValidationLog("WARNING",
        "non-finite reference position for " + describeMeshVertexSourceKey(key) + " at "
          + records.front().meshlet_label + " v" + std::to_string(records.front().vertex_index));
      for (const MeshVertexSourceRecord& record : records)
      {
        result.discrepancies.push_back(MeshVertexSourceDiscrepancy { record.meshlet_index, record.vertex_index });
      }
      continue;
    }

    double max_distance = 0.0;
    for (size_t i = 1; i < records.size(); ++i)
    {
      const glm::dvec3& candidate = records[i].world_position;
      if (!std::isfinite(candidate.x) || !std::isfinite(candidate.y) || !std::isfinite(candidate.z))
      {
        max_distance = std::numeric_limits<double>::infinity();
        break;
      }
      max_distance = std::max(max_distance, glm::length(candidate - reference));
    }

    if (max_distance <= position_tolerance)
    {
      continue;
    }

    ++result.inconsistent_source_count;
    std::ostringstream oss;
    oss << "inconsistent world positions for " << describeMeshVertexSourceKey(key) << " (max_distance="
        << numberLiteral(max_distance) << ", tolerance=" << numberLiteral(position_tolerance) << "):";
    for (const MeshVertexSourceRecord& record : records)
    {
      result.discrepancies.push_back(MeshVertexSourceDiscrepancy { record.meshlet_index, record.vertex_index });
      oss << " [" << record.meshlet_label << " v" << record.vertex_index << " pos=(" << numberLiteral(record.world_position.x)
          << "," << numberLiteral(record.world_position.y) << "," << numberLiteral(record.world_position.z) << ")]";
    }
    writeValidationLog("WARNING", oss.str());
  }

  return result;
}
} // namespace

void Validator::setLogFile(const std::string& path)
{
  if (validation_log_file.is_open())
  {
    validation_log_file.close();
  }
  validation_log_path = path.empty() ? defaultLogFilePath() : path;
}

const std::string& Validator::logFilePath() { return validation_log_path; }

MeshVertexSourceValidationResult Validator::validateMeshVertexSourcesHaveConsistentWorldPositions(
  const std::vector<VoronoiMesh*>& meshlets, const std::vector<std::string>* meshlet_labels, double position_tolerance)
{
  std::unordered_map<MeshVertexSourceKey, std::vector<MeshVertexSourceRecord>, MeshVertexSourceKeyHash> buckets;

  for (size_t meshlet_index = 0; meshlet_index < meshlets.size(); ++meshlet_index)
  {
    if (meshlets[meshlet_index] == nullptr)
    {
      continue;
    }
    const VoronoiMesh& mesh = *meshlets[meshlet_index];
    if (!mesh.storeMetadata())
    {
      continue;
    }

    const auto& vertices = mesh.getVertices();
    const auto& vertex_metadata = mesh.getVertexMetadata();
    const size_t metadata_count = std::min(vertices.size(), vertex_metadata.size());

    std::string label = "meshlet_" + std::to_string(meshlet_index);
    if (meshlet_labels != nullptr && meshlet_index < meshlet_labels->size() && !(*meshlet_labels)[meshlet_index].empty())
    {
      label = (*meshlet_labels)[meshlet_index];
    }

    for (size_t vertex_index = 0; vertex_index < metadata_count; ++vertex_index)
    {
      const std::optional<MeshVertexSourceKey> key = meshVertexSourceKeyFromMetadata(vertex_metadata[vertex_index]);
      if (!key.has_value())
      {
        continue;
      }

      buckets[key.value()].push_back(
        MeshVertexSourceRecord { meshlet_index, label, vertex_index, vertices[vertex_index] });
    }
  }

  return validateMeshVertexSourceRecords(buckets, position_tolerance);
}

size_t Validator::markMeshVertexSourceDiscrepancies(const std::vector<VoronoiMesh*>& meshlets,
  const std::vector<MeshVertexSourceDiscrepancy>& discrepancies)
{
  if (meshlets.empty() || discrepancies.empty())
  {
    return 0;
  }

  static const glm::dvec3 kRed(1.0, 0.0, 0.0);
  size_t marked_triangle_count = 0;

  std::unordered_map<size_t, std::vector<size_t>> vertices_by_meshlet;
  vertices_by_meshlet.reserve(discrepancies.size());
  for (const MeshVertexSourceDiscrepancy& discrepancy : discrepancies)
  {
    if (discrepancy.meshlet_index >= meshlets.size() || discrepancy.vertex_index == static_cast<size_t>(-1))
    {
      continue;
    }
    vertices_by_meshlet[discrepancy.meshlet_index].push_back(discrepancy.vertex_index);
  }

  for (const auto& [meshlet_index, vertex_indices] : vertices_by_meshlet)
  {
    VoronoiMesh* mesh = meshlets[meshlet_index];
    if (mesh == nullptr)
    {
      continue;
    }

    const int red_material_id = mesh->ensureMaterialName(validationErrorMaterialName());
    std::unordered_set<size_t> triangle_indices;
    triangle_indices.reserve(vertex_indices.size() * 2);

    for (size_t vertex_index : vertex_indices)
    {
      mesh->setVertexColor(vertex_index, kRed);
      const std::vector<size_t> corner_indices = mesh->findTriangleCorners(vertex_index, false);
      for (size_t corner_index : corner_indices)
      {
        triangle_indices.insert(corner_index / 3);
      }
    }

    for (size_t triangle_index : triangle_indices)
    {
      mesh->setTriangleMaterialId(triangle_index, red_material_id);
      ++marked_triangle_count;
    }
  }

  return marked_triangle_count;
}

bool Validator::meshUsesValidationErrorMaterial(const VoronoiMesh& mesh)
{
  const auto& material_names = mesh.getMaterialNames();
  const auto material_it = std::find(material_names.begin(), material_names.end(), validationErrorMaterialName());
  if (material_it == material_names.end())
  {
    return false;
  }

  const int red_material_id = static_cast<int>(std::distance(material_names.begin(), material_it));
  for (int material_id : mesh.getMaterialIDs())
  {
    if (material_id == red_material_id)
    {
      return true;
    }
  }
  return false;
}

void Validator::logMeshVertexSourceValidationResult(const MeshVertexSourceValidationResult& result, const char* scope)
{
  if (result.inconsistent_source_count > 0)
  {
    std::ostringstream oss;
    oss << scope << ": " << result.inconsistent_source_count << " inconsistent source group(s) among "
        << result.keyed_vertex_count << " keyed vertices (" << result.unique_source_count << " unique sources)";
    if (result.marked_triangle_count > 0)
    {
      oss << "; marked " << result.marked_triangle_count << " triangle(s) with material '" << validationErrorMaterialName()
          << "'";
    }
    writeValidationLog("WARNING", oss.str());
    return;
  }

  if (result.keyed_vertex_count > 0)
  {
    std::ostringstream oss;
    oss << scope << ": all " << result.unique_source_count << " keyed source groups agree in world position ("
        << result.keyed_vertex_count << " vertices)";
    writeValidationLog("INFO", oss.str());
    return;
  }

  writeValidationLog("INFO", std::string(scope) + ": no keyed vertices found.");
}

bool Validator::validateAndReportMeshVertexSources(const std::vector<VoronoiMesh*>& meshlets, const char* scope,
  const std::vector<std::string>* meshlet_labels, double position_tolerance, bool mark_discrepancies)
{
  beginValidationSession(scope);
  writeValidationLog("INFO", std::string("log_file=") + validation_log_path);
  MeshVertexSourceValidationResult result
    = validateMeshVertexSourcesHaveConsistentWorldPositions(meshlets, meshlet_labels, position_tolerance);
  if (mark_discrepancies && !result.discrepancies.empty())
  {
    result.marked_triangle_count = markMeshVertexSourceDiscrepancies(meshlets, result.discrepancies);
  }
  logMeshVertexSourceValidationResult(result, scope);
  return result.passed();
}

bool Validator::validateAndReportMeshVertexSources(std::vector<VoronoiMesh>& meshlets, const char* scope,
  const std::vector<std::string>* meshlet_labels, double position_tolerance, bool mark_discrepancies)
{
  std::vector<VoronoiMesh*> meshlet_ptrs;
  meshlet_ptrs.reserve(meshlets.size());
  for (VoronoiMesh& mesh : meshlets)
  {
    meshlet_ptrs.push_back(&mesh);
  }
  return validateAndReportMeshVertexSources(
    meshlet_ptrs, scope, meshlet_labels, position_tolerance, mark_discrepancies);
}

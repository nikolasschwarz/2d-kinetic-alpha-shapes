#include "MeshingBuffer.hpp"

#include "Logger.hpp"

#include <chrono>
#include <cstddef>
#include <cstring>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <sstream>

namespace kinDS {
namespace {

constexpr uint64_t kMaxCount = 500000000ull;
constexpr uint64_t kMaxString = 16ull * 1024ull * 1024ull;
constexpr uint64_t kFnvOffset = 14695981039346656037ull;
constexpr uint64_t kFnvPrime = 1099511628211ull;

struct Fnv64 {
  uint64_t value = kFnvOffset;
  void MixBytes(const void* data, size_t size) {
    const auto* bytes = static_cast<const uint8_t*>(data);
    for (size_t i = 0; i < size; ++i) {
      value ^= bytes[i];
      value *= kFnvPrime;
    }
  }
  void MixCString(const char* text) {
    MixBytes(text, std::strlen(text));
  }
  template <typename T>
  void MixPod(const T& value_pod) {
    MixBytes(&value_pod, sizeof(T));
  }
  template <typename T>
  void MixVec(const std::vector<T>& values) {
    MixPod(static_cast<uint64_t>(values.size()));
    if (!values.empty()) {
      MixBytes(values.data(), values.size() * sizeof(T));
    }
  }
  template <typename T>
  void MixNested(const std::vector<std::vector<T>>& values) {
    MixPod(static_cast<uint64_t>(values.size()));
    for (const auto& inner : values) {
      MixVec(inner);
    }
  }
};

std::string HashToHex(uint64_t hash) {
  std::ostringstream stream;
  stream << std::hex << std::setw(16) << std::setfill('0') << hash;
  return stream.str();
}

template <typename Nested>
size_t CountNestedElements(const Nested& nested) {
  size_t count = 0;
  for (const auto& inner : nested) {
    count += inner.size();
  }
  return count;
}

size_t CountTripleNestedElements(const std::vector<std::vector<std::vector<size_t>>>& nested) {
  size_t count = 0;
  for (const auto& by_height : nested) {
    count += CountNestedElements(by_height);
  }
  return count;
}

std::string FormatRootTransformSummary(const glm::mat4& m) {
  const glm::vec3 t = m[3];
  std::ostringstream oss;
  oss << std::fixed << std::setprecision(6) << "t=(" << t.x << "," << t.y << "," << t.z << ")"
      << " m00=" << m[0][0] << " m11=" << m[1][1] << " m22=" << m[2][2];
  return oss.str();
}

class BinaryWriter {
 public:
  explicit BinaryWriter(const std::filesystem::path& path) : out_(path, std::ios::binary | std::ios::trunc) {
  }
  bool Good() const {
    return static_cast<bool>(out_);
  }
  template <typename T>
  void WritePod(const T& value) {
    out_.write(reinterpret_cast<const char*>(&value), static_cast<std::streamsize>(sizeof(T)));
  }
  void WriteBytes(const void* data, size_t size) {
    if (size == 0) {
      return;
    }
    out_.write(reinterpret_cast<const char*>(data), static_cast<std::streamsize>(size));
  }
  template <typename T>
  void WriteVec(const std::vector<T>& values) {
    WritePod(static_cast<uint64_t>(values.size()));
    if (!values.empty()) {
      WriteBytes(values.data(), values.size() * sizeof(T));
    }
  }
  void WriteByteBlob(uint64_t element_count, const std::vector<std::byte>& bytes) {
    WritePod(element_count);
    WriteBytes(bytes.data(), bytes.size());
  }
  void WriteSizeTVec(const std::vector<size_t>& values) {
    WritePod(static_cast<uint64_t>(values.size()));
    if constexpr (sizeof(size_t) == 8) {
      if (!values.empty()) {
        WriteBytes(values.data(), values.size() * sizeof(size_t));
      }
    } else {
      for (const size_t value : values) {
        WritePod(static_cast<uint64_t>(value));
      }
    }
  }
  void WriteString(const std::string& value) {
    WritePod(static_cast<uint64_t>(value.size()));
    WriteBytes(value.data(), value.size());
  }
  void WriteStrings(const std::vector<std::string>& values) {
    WritePod(static_cast<uint64_t>(values.size()));
    for (const auto& value : values) {
      WriteString(value);
    }
  }
  void WriteNestedInt(const std::vector<std::vector<int>>& values) {
    WritePod(static_cast<uint64_t>(values.size()));
    for (const auto& inner : values) {
      WriteVec(inner);
    }
  }
  void WriteNestedSizeT(const std::vector<std::vector<size_t>>& values) {
    WritePod(static_cast<uint64_t>(values.size()));
    for (const auto& inner : values) {
      WriteSizeTVec(inner);
    }
  }

 private:
  std::ofstream out_;
};

class BinaryReader {
 public:
  explicit BinaryReader(const std::filesystem::path& path) : in_(path, std::ios::binary) {
  }
  bool Good() const {
    return !failed_ && static_cast<bool>(in_);
  }
  void Fail() {
    failed_ = true;
  }
  template <typename T>
  T ReadPod() {
    T value{};
    in_.read(reinterpret_cast<char*>(&value), static_cast<std::streamsize>(sizeof(T)));
    if (!in_) {
      failed_ = true;
    }
    return value;
  }
  template <typename T>
  std::vector<T> ReadVec() {
    const uint64_t count = ReadPod<uint64_t>();
    if (failed_ || count > kMaxCount) {
      failed_ = true;
      return {};
    }
    std::vector<T> values(static_cast<size_t>(count));
    if (count > 0) {
      in_.read(reinterpret_cast<char*>(values.data()), static_cast<std::streamsize>(count * sizeof(T)));
      if (!in_) {
        failed_ = true;
        return {};
      }
    }
    return values;
  }
  std::vector<std::byte> ReadByteBlob(uint32_t stride) {
    const uint64_t count = ReadPod<uint64_t>();
    if (failed_ || count > kMaxCount) {
      failed_ = true;
      return {};
    }
    if (stride == 0) {
      if (count != 0) {
        failed_ = true;
      }
      return {};
    }
    const uint64_t byte_count = count * static_cast<uint64_t>(stride);
    if (byte_count / stride != count) {
      failed_ = true;
      return {};
    }
    std::vector<std::byte> values(static_cast<size_t>(byte_count));
    if (byte_count > 0) {
      in_.read(reinterpret_cast<char*>(values.data()), static_cast<std::streamsize>(byte_count));
      if (!in_) {
        failed_ = true;
        return {};
      }
    }
    return values;
  }
  std::vector<size_t> ReadSizeTVec() {
    const uint64_t count = ReadPod<uint64_t>();
    if (failed_ || count > kMaxCount) {
      failed_ = true;
      return {};
    }
    std::vector<size_t> values(static_cast<size_t>(count));
    if constexpr (sizeof(size_t) == 8) {
      if (count > 0) {
        in_.read(reinterpret_cast<char*>(values.data()), static_cast<std::streamsize>(count * sizeof(size_t)));
        if (!in_) {
          failed_ = true;
          return {};
        }
      }
    } else {
      for (uint64_t i = 0; i < count; ++i) {
        values[static_cast<size_t>(i)] = static_cast<size_t>(ReadPod<uint64_t>());
      }
    }
    return values;
  }
  std::string ReadString() {
    const uint64_t count = ReadPod<uint64_t>();
    if (failed_ || count > kMaxString) {
      failed_ = true;
      return {};
    }
    std::string value(static_cast<size_t>(count), '\0');
    if (count > 0) {
      in_.read(value.data(), static_cast<std::streamsize>(count));
      if (!in_) {
        failed_ = true;
        return {};
      }
    }
    return value;
  }
  std::vector<std::string> ReadStrings() {
    const uint64_t count = ReadPod<uint64_t>();
    if (failed_ || count > kMaxCount) {
      failed_ = true;
      return {};
    }
    std::vector<std::string> values;
    values.reserve(static_cast<size_t>(count));
    for (uint64_t i = 0; i < count; ++i) {
      values.push_back(ReadString());
      if (failed_) {
        return {};
      }
    }
    return values;
  }
  std::vector<std::vector<int>> ReadNestedInt() {
    const uint64_t count = ReadPod<uint64_t>();
    if (failed_ || count > kMaxCount) {
      failed_ = true;
      return {};
    }
    std::vector<std::vector<int>> values(static_cast<size_t>(count));
    for (auto& inner : values) {
      inner = ReadVec<int>();
      if (failed_) {
        return {};
      }
    }
    return values;
  }
  std::vector<std::vector<size_t>> ReadNestedSizeT() {
    const uint64_t count = ReadPod<uint64_t>();
    if (failed_ || count > kMaxCount) {
      failed_ = true;
      return {};
    }
    std::vector<std::vector<size_t>> values(static_cast<size_t>(count));
    for (auto& inner : values) {
      inner = ReadSizeTVec();
      if (failed_) {
        return {};
      }
    }
    return values;
  }

 private:
  std::ifstream in_;
  bool failed_ = false;
};

void WriteVoronoiMesh(BinaryWriter& writer, const VoronoiMesh& mesh) {
  writer.WritePod(static_cast<int32_t>(mesh.getNormalMode()));
  writer.WritePod(static_cast<uint8_t>(mesh.storeMetadata() ? 1 : 0));
  writer.WritePod(mesh.getCreationKineticTime());
  writer.WriteStrings(mesh.getMaterialNames());
  writer.WriteVec(mesh.getVertices());
  writer.WriteSizeTVec(mesh.getTriangles());
  writer.WriteVec(mesh.getNormals());
  writer.WriteVec(mesh.getUVs());
  writer.WriteSizeTVec(mesh.getUVIndices());
  writer.WriteVec(mesh.getMaterialIDs());
  writer.WriteSizeTVec(mesh.getGroupOffsets());
  writer.WriteStrings(mesh.getGroupNames());
  writer.WriteVec(mesh.getVertexColors());
  writer.WriteStrings(mesh.getVertexMetadata());
  writer.WriteStrings(mesh.getFaceMetadata());
  writer.WriteVec(mesh.getProfilePlaneXY());
  writer.WriteVec(mesh.getVertexKineticTimes());
  writer.WriteVec(mesh.getVertexSemanticUvs());
  writer.WritePod(static_cast<uint64_t>(mesh.getVertexCount()));
  for (size_t i = 0; i < mesh.getVertexCount(); ++i) {
    writer.WritePod(static_cast<uint8_t>(mesh.isVertexFlexible(i) ? 1 : 0));
  }
}

VoronoiMesh ReadVoronoiMesh(BinaryReader& reader) {
  const auto normal_mode = static_cast<NormalMode>(reader.ReadPod<int32_t>());
  const bool store_metadata = reader.ReadPod<uint8_t>() != 0;
  const double creation_time = reader.ReadPod<double>();
  auto material_names = reader.ReadStrings();
  VoronoiMesh mesh(std::move(material_names), normal_mode);
  mesh.setStoreMetadata(store_metadata);
  mesh.setCreationKineticTime(creation_time);
  mesh.getVertices() = reader.ReadVec<glm::dvec3>();
  mesh.getTriangles() = reader.ReadSizeTVec();
  mesh.getNormals() = reader.ReadVec<glm::dvec3>();
  mesh.getUVs() = reader.ReadVec<glm::dvec3>();
  mesh.getUVIndices() = reader.ReadSizeTVec();
  mesh.getMaterialIDs() = reader.ReadVec<int>();
  mesh.setGroupOffsets(reader.ReadSizeTVec());
  mesh.setGroupNames(reader.ReadStrings());
  mesh.getVertexColors() = reader.ReadVec<glm::dvec3>();
  mesh.getVertexMetadata() = reader.ReadStrings();
  mesh.getFaceMetadata() = reader.ReadStrings();
  mesh.getProfilePlaneXY() = reader.ReadVec<glm::dvec2>();
  mesh.getVertexKineticTimes() = reader.ReadVec<double>();
  mesh.getVertexSemanticUvs() = reader.ReadVec<glm::dvec3>();
  const uint64_t flexible_count = reader.ReadPod<uint64_t>();
  if (!reader.Good() || flexible_count > kMaxCount) {
    reader.Fail();
    return mesh;
  }
  for (uint64_t i = 0; i < flexible_count; ++i) {
    const uint8_t flag = reader.ReadPod<uint8_t>();
    if (flag && i < mesh.getVertexCount()) {
      mesh.setVertexFlexible(static_cast<size_t>(i), true);
    }
  }
  return mesh;
}

std::string CurrentUtcTimestamp() {
  const std::time_t time = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
  std::tm utc{};
#if defined(_WIN32)
  gmtime_s(&utc, &time);
#else
  gmtime_r(&time, &utc);
#endif
  std::ostringstream timestamp;
  timestamp << std::put_time(&utc, "%Y-%m-%dT%H:%M:%SZ");
  return timestamp.str();
}

std::string YamlQuoteIfNeeded(const std::string& value) {
  bool needs_quotes = value.empty();
  if (!needs_quotes) {
    for (const char c : value) {
      if (c == ':' || c == '#' || c == '"' || c == '\'' || c == '\\' || c == '\n' || c == '\r' || c == '\t' ||
          c == '{' || c == '}' || c == '[' || c == ']' || c == ',' || c == '&' || c == '*' || c == '!' || c == '|' ||
          c == '>' || c == '%' || c == '@' || c == '`') {
        needs_quotes = true;
        break;
      }
    }
    if (!needs_quotes && (value.front() == ' ' || value.back() == ' ')) {
      needs_quotes = true;
    }
  }
  if (!needs_quotes) {
    return value;
  }
  std::string out = "\"";
  for (const char c : value) {
    if (c == '\\' || c == '"') {
      out += '\\';
    }
    if (c == '\n') {
      out += "\\n";
      continue;
    }
    if (c == '\r') {
      continue;
    }
    out += c;
  }
  out += '"';
  return out;
}

/// Legacy EcoSysLab MeshBuffers sidecar layout (top-level keys only; statistics filled by EcoSysLab when available).
bool WriteMeshingBufferYml(const std::filesystem::path& yml_path, const MeshingBufferMetadata& metadata) {
  std::ofstream out(yml_path);
  if (!out) {
    return false;
  }
  const std::string created =
      metadata.created_utc.empty() ? CurrentUtcTimestamp() : metadata.created_utc;
  out << "format: " << meshingBufferMagic() << "\n";
  out << "version: " << meshingBufferVersion() << "\n";
  out << "hash: " << YamlQuoteIfNeeded(metadata.input_hash) << "\n";
  out << "created_utc: " << YamlQuoteIfNeeded(created) << "\n";
  out << "gpu_vertex_count: " << metadata.gpu_vertex_count << "\n";
  out << "gpu_triangle_count: " << metadata.gpu_triangle_count << "\n";
  out << "meshlet_count: " << metadata.meshlet_count << "\n";
  out << "vertex_stride: " << metadata.vertex_stride << "\n";
  out << "triangle_stride: " << metadata.triangle_stride << "\n";
  out << "spline_tension: " << metadata.spline_tension << "\n";
  out << "alpha_cutoff: " << metadata.alpha_cutoff << "\n";
  out << "branch_alpha_cutoff: " << metadata.branch_alpha_cutoff << "\n";
  out << "look_ahead: " << metadata.look_ahead << "\n";
  out << "store_mesh_metadata: " << (metadata.store_mesh_metadata ? "true" : "false") << "\n";
  out << "mesh_cap_at_start: true\n";
  out << "transform_mesh_at_construction: true\n";
  out << "description: " << YamlQuoteIfNeeded(metadata.description) << "\n";
  return static_cast<bool>(out);
}

std::string StripYamlScalar(std::string value) {
  while (!value.empty() && (value.front() == ' ' || value.front() == '\t')) {
    value.erase(value.begin());
  }
  while (!value.empty() && (value.back() == ' ' || value.back() == '\t' || value.back() == '\r')) {
    value.pop_back();
  }
  if (value.size() >= 2 && ((value.front() == '"' && value.back() == '"') || (value.front() == '\'' && value.back() == '\''))) {
    value = value.substr(1, value.size() - 2);
  }
  return value;
}

}  // namespace

MeshingBufferHashStats computeMeshingInputHash(
    const std::vector<std::vector<glm::dvec2>>& support_points,
    const std::vector<std::vector<double>>& subdivisions_by_strand,
    const std::vector<std::vector<int>>& physics_strand_to_segment_indices,
    const std::vector<std::vector<glm::dmat4>>& transforms_by_height_and_branch, const glm::mat4& root_transform,
    const std::vector<std::vector<size_t>>& branch_indices,
    const std::vector<std::vector<std::vector<size_t>>>& strands_by_branch_id,
    const MeshingBufferHashSettings& settings, bool include_alpha_cutoff) {
  const auto mix_settings = [&settings, include_alpha_cutoff](Fnv64& hash) {
    hash.MixCString("DsKineticVoronoiMeshing.v1");
    hash.MixPod(meshingBufferVersion());
    hash.MixPod(static_cast<uint8_t>(1));  // mesh_cap_at_start
    hash.MixPod(static_cast<uint8_t>(1));  // transform_mesh_at_construction
    hash.MixPod(static_cast<uint8_t>(settings.store_mesh_metadata ? 1 : 0));
    hash.MixPod(settings.spline_tension);
    if (settings.bark_smooth_iterations != 0) {
      hash.MixPod(settings.bark_smooth_iterations);
      hash.MixPod(settings.bark_smooth_strength);
      if (!settings.bark_smooth_lock_boundary) {
        hash.MixPod(static_cast<uint8_t>(0));
      }
      if (!settings.bark_smooth_uvs) {
        hash.MixPod(static_cast<uint8_t>(0));
      }
    }
    if (settings.bark_subdivide) {
      hash.MixPod(static_cast<uint8_t>(1));
    }
    if (include_alpha_cutoff) {
      hash.MixPod(settings.alpha_cutoff);
      if (settings.branch_alpha_cutoff != settings.alpha_cutoff) {
        hash.MixPod(settings.branch_alpha_cutoff);
      }
      if (settings.look_ahead != 0) {
        hash.MixPod(settings.look_ahead);
      }
      if (settings.hinge_only_profile_plane_mix) {
        hash.MixPod(static_cast<uint8_t>(1));
      }
    }
  };

  Fnv64 settings_hash;
  mix_settings(settings_hash);

  Fnv64 root_hash;
  root_hash.MixPod(root_transform);

  Fnv64 support_hash;
  support_hash.MixNested(support_points);

  Fnv64 subdiv_hash;
  subdiv_hash.MixNested(subdivisions_by_strand);

  Fnv64 physics_hash;
  physics_hash.MixNested(physics_strand_to_segment_indices);

  Fnv64 transforms_hash;
  transforms_hash.MixNested(transforms_by_height_and_branch);

  Fnv64 branch_hash;
  branch_hash.MixNested(branch_indices);

  Fnv64 strands_by_branch_hash;
  strands_by_branch_hash.MixPod(static_cast<uint64_t>(strands_by_branch_id.size()));
  for (const auto& by_height : strands_by_branch_id) {
    strands_by_branch_hash.MixNested(by_height);
  }

  Fnv64 hash;
  mix_settings(hash);
  hash.MixPod(root_transform);
  hash.MixNested(support_points);
  hash.MixNested(subdivisions_by_strand);
  hash.MixNested(physics_strand_to_segment_indices);
  hash.MixNested(transforms_by_height_and_branch);
  hash.MixNested(branch_indices);
  hash.MixPod(static_cast<uint64_t>(strands_by_branch_id.size()));
  for (const auto& by_height : strands_by_branch_id) {
    hash.MixNested(by_height);
  }

  MeshingBufferHashStats stats;
  stats.input_hash = HashToHex(hash.value);
  stats.settings_hash = HashToHex(settings_hash.value);
  stats.root_hash = HashToHex(root_hash.value);
  stats.support_hash = HashToHex(support_hash.value);
  stats.subdiv_hash = HashToHex(subdiv_hash.value);
  stats.physics_hash = HashToHex(physics_hash.value);
  stats.transforms_hash = HashToHex(transforms_hash.value);
  stats.branch_hash = HashToHex(branch_hash.value);
  stats.strands_by_branch_hash = HashToHex(strands_by_branch_hash.value);
  stats.root_transform_summary = FormatRootTransformSummary(root_transform);
  stats.support_strand_count = support_points.size();
  stats.support_point_count = CountNestedElements(support_points);
  stats.subdiv_strand_count = subdivisions_by_strand.size();
  stats.subdiv_value_count = CountNestedElements(subdivisions_by_strand);
  stats.physics_strand_count = physics_strand_to_segment_indices.size();
  stats.physics_segment_count = CountNestedElements(physics_strand_to_segment_indices);
  stats.transform_height_count = transforms_by_height_and_branch.size();
  stats.transform_matrix_count = CountNestedElements(transforms_by_height_and_branch);
  stats.branch_height_count = branch_indices.size();
  stats.branch_id_count = CountNestedElements(branch_indices);
  stats.strands_by_branch_outer_count = strands_by_branch_id.size();
  stats.strands_by_branch_id_count = CountTripleNestedElements(strands_by_branch_id);
  return stats;
}

bool saveMeshingBuffer(const std::filesystem::path& bin_path, const MeshingBufferPayload& payload,
                       const MeshingBufferMetadata& metadata) {
  std::error_code error;
  std::filesystem::create_directories(bin_path.parent_path(), error);
  if (error) {
    KINDS_ERROR("Failed to create meshing buffer directory: " << error.message());
    return false;
  }

  const uint32_t vertex_stride = payload.gpu_vertex_stride;
  const uint32_t triangle_stride = payload.gpu_triangle_stride;
  const uint64_t gpu_vertex_count =
      vertex_stride == 0 ? 0ull : static_cast<uint64_t>(payload.gpu_vertices.size() / vertex_stride);
  const uint64_t gpu_triangle_count =
      triangle_stride == 0 ? 0ull : static_cast<uint64_t>(payload.gpu_triangles.size() / triangle_stride);
  if ((vertex_stride == 0 && !payload.gpu_vertices.empty()) ||
      (triangle_stride == 0 && !payload.gpu_triangles.empty()) ||
      (vertex_stride != 0 && payload.gpu_vertices.size() % vertex_stride != 0) ||
      (triangle_stride != 0 && payload.gpu_triangles.size() % triangle_stride != 0)) {
    KINDS_ERROR("Meshing buffer GPU blob size does not match stride");
    return false;
  }

  BinaryWriter writer(bin_path);
  writer.WriteBytes(meshingBufferMagic(), 4);
  writer.WritePod(meshingBufferVersion());
  writer.WritePod(vertex_stride);
  writer.WritePod(triangle_stride);
  writer.WritePod(payload.root_transform);
  writer.WriteByteBlob(gpu_vertex_count, payload.gpu_vertices);
  writer.WriteByteBlob(gpu_triangle_count, payload.gpu_triangles);
  writer.WritePod(static_cast<uint64_t>(payload.meshlets.size()));
  for (const auto& meshlet : payload.meshlets) {
    WriteVoronoiMesh(writer, meshlet);
  }
  writer.WriteNestedInt(payload.neighbors);
  writer.WriteSizeTVec(payload.meshing_to_physics);
  writer.WriteNestedSizeT(payload.strand_to_segment);
  if (!writer.Good()) {
    KINDS_ERROR("Failed to write meshing buffer " << bin_path.string());
    return false;
  }

  MeshingBufferMetadata meta = metadata;
  meta.meshlet_count = payload.meshlets.size();
  meta.gpu_vertex_count = static_cast<size_t>(gpu_vertex_count);
  meta.gpu_triangle_count = static_cast<size_t>(gpu_triangle_count);
  meta.vertex_stride = vertex_stride;
  meta.triangle_stride = triangle_stride;
  const std::filesystem::path yml_path = meshingBufferYmlPath(bin_path);
  if (!WriteMeshingBufferYml(yml_path, meta)) {
    KINDS_ERROR("Wrote meshing buffer binary but failed to write metadata " << yml_path.string());
  }
  KINDS_INFO("Saved meshing buffer " << bin_path.string() << " (" << meta.meshlet_count << " meshlets, "
                                     << meta.gpu_vertex_count << " gpu verts)");
  return true;
}

bool loadMeshingBuffer(const std::filesystem::path& bin_path, const glm::mat4& expected_root,
                       MeshingBufferPayload& out_payload, bool require_root_match) {
  BinaryReader reader(bin_path);
  char magic[4]{};
  magic[0] = reader.ReadPod<char>();
  magic[1] = reader.ReadPod<char>();
  magic[2] = reader.ReadPod<char>();
  magic[3] = reader.ReadPod<char>();
  if (!reader.Good() || std::memcmp(magic, meshingBufferMagic(), 4) != 0) {
    return false;
  }
  if (reader.ReadPod<uint32_t>() != meshingBufferVersion()) {
    return false;
  }
  out_payload = {};
  out_payload.gpu_vertex_stride = reader.ReadPod<uint32_t>();
  out_payload.gpu_triangle_stride = reader.ReadPod<uint32_t>();
  out_payload.root_transform = reader.ReadPod<glm::mat4>();
  if (require_root_match && out_payload.root_transform != expected_root) {
    return false;
  }
  out_payload.gpu_vertices = reader.ReadByteBlob(out_payload.gpu_vertex_stride);
  out_payload.gpu_triangles = reader.ReadByteBlob(out_payload.gpu_triangle_stride);
  const uint64_t meshlet_count = reader.ReadPod<uint64_t>();
  if (!reader.Good() || meshlet_count > kMaxCount) {
    return false;
  }
  out_payload.meshlets.clear();
  out_payload.meshlets.reserve(static_cast<size_t>(meshlet_count));
  for (uint64_t i = 0; i < meshlet_count; ++i) {
    out_payload.meshlets.push_back(ReadVoronoiMesh(reader));
    if (!reader.Good()) {
      return false;
    }
  }
  out_payload.neighbors = reader.ReadNestedInt();
  out_payload.meshing_to_physics = reader.ReadSizeTVec();
  out_payload.strand_to_segment = reader.ReadNestedSizeT();
  return reader.Good();
}

bool loadMeshingBufferMetadata(const std::filesystem::path& yml_path, MeshingBufferMetadata& out_metadata) {
  std::ifstream in(yml_path);
  if (!in) {
    return false;
  }
  out_metadata = {};
  std::string line;
  while (std::getline(in, line)) {
    if (line.empty() || line[0] == '#' || line[0] == ' ' || line[0] == '\t') {
      continue;  // skip blanks, comments, and nested map entries
    }
    const auto colon = line.find(':');
    if (colon == std::string::npos) {
      continue;
    }
    const std::string key = line.substr(0, colon);
    const std::string value = StripYamlScalar(line.substr(colon + 1));
    if (key == "hash" || key == "input_hash") {
      out_metadata.input_hash = value;
    } else if (key == "description") {
      out_metadata.description = value;
    } else if (key == "created_utc") {
      out_metadata.created_utc = value;
    } else if (key == "alpha_cutoff") {
      out_metadata.alpha_cutoff = std::stod(value);
    } else if (key == "branch_alpha_cutoff") {
      out_metadata.branch_alpha_cutoff = std::stod(value);
    } else if (key == "look_ahead") {
      out_metadata.look_ahead = static_cast<size_t>(std::stoull(value));
    } else if (key == "spline_tension") {
      out_metadata.spline_tension = std::stof(value);
    } else if (key == "store_mesh_metadata") {
      out_metadata.store_mesh_metadata = (value == "true" || value == "1");
    } else if (key == "meshlet_count") {
      out_metadata.meshlet_count = static_cast<size_t>(std::stoull(value));
    } else if (key == "gpu_vertex_count") {
      out_metadata.gpu_vertex_count = static_cast<size_t>(std::stoull(value));
    } else if (key == "gpu_triangle_count") {
      out_metadata.gpu_triangle_count = static_cast<size_t>(std::stoull(value));
    } else if (key == "vertex_stride") {
      out_metadata.vertex_stride = static_cast<uint32_t>(std::stoul(value));
    } else if (key == "triangle_stride") {
      out_metadata.triangle_stride = static_cast<uint32_t>(std::stoul(value));
    }
  }
  return !out_metadata.input_hash.empty() || out_metadata.meshlet_count > 0;
}

}  // namespace kinDS

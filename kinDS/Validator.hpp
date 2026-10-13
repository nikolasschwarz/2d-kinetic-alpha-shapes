#pragma once

#include "VoronoiMesh.hpp"

#include <cstddef>
#include <string>
#include <vector>

namespace kinDS
{
struct MeshVertexSourceDiscrepancy
{
  size_t meshlet_index = static_cast<size_t>(-1);
  size_t vertex_index = static_cast<size_t>(-1);
};

struct MeshVertexSourceValidationResult
{
  size_t keyed_vertex_count = 0;
  size_t unique_source_count = 0;
  size_t inconsistent_source_count = 0;
  size_t marked_triangle_count = 0;
  std::vector<MeshVertexSourceDiscrepancy> discrepancies;

  bool passed() const { return inconsistent_source_count == 0; }
  bool hasMarkedGeometry() const { return marked_triangle_count > 0; }
};

class Validator
{
 public:
  static constexpr const char* defaultLogFilePath() { return "mesh_vertex_validation.log"; }
  /// OBJ/MTL material used for errors on non-green source triangles.
  static constexpr const char* validationErrorMaterialName() { return "red"; }
  /// OBJ/MTL material used for errors on green source triangles.
  static constexpr const char* validationGreenErrorMaterialName() { return "neon_yellow"; }

  /// Route validation output to @p path (append mode). Reopens the file when the path changes.
  static void setLogFile(const std::string& path);

  /// Current validation log path (explicitly set or @ref defaultLogFilePath).
  static const std::string& logFilePath();

  /// Hash-bucket check: vertices sharing the same metadata source identity must agree in stored
  /// world/object-space coordinates. Ignores vertices without a comparable key.
  static MeshVertexSourceValidationResult validateMeshVertexSourcesHaveConsistentWorldPositions(
    const std::vector<VoronoiMesh*>& meshlets, const std::vector<std::string>* meshlet_labels = nullptr,
    double position_tolerance = 1e-5);

  /// Mark triangles incident to discrepancies red, or neon yellow when their original material was green.
  static size_t markMeshVertexSourceDiscrepancies(const std::vector<VoronoiMesh*>& meshlets,
    const std::vector<MeshVertexSourceDiscrepancy>& discrepancies);

  /// @return true when @p mesh uses either validation-error material on at least one triangle.
  static bool meshUsesValidationErrorMaterial(const VoronoiMesh& mesh);

  /// Write a summary of @p result to the validation log file.
  static void logMeshVertexSourceValidationResult(const MeshVertexSourceValidationResult& result, const char* scope);

  /// Run validation, optionally mark discrepancies on the meshlets, and write the report to the validation log file.
  /// @return @c true when no inconsistent source groups were found.
  static bool validateAndReportMeshVertexSources(const std::vector<VoronoiMesh*>& meshlets, const char* scope,
    const std::vector<std::string>* meshlet_labels = nullptr, double position_tolerance = 1e-5,
    bool mark_discrepancies = true);

  static bool validateAndReportMeshVertexSources(std::vector<VoronoiMesh>& meshlets, const char* scope,
    const std::vector<std::string>* meshlet_labels = nullptr, double position_tolerance = 1e-5,
    bool mark_discrepancies = true);

  /// Check green (interior) triangle corners only: uv.z == vertex kinetic time * @p uv_height_factor.
  static bool validateAndReportInteriorMeshletUvHeights(const std::vector<VoronoiMesh*>& meshlets, const char* scope,
    double uv_height_factor, const std::vector<std::string>* meshlet_labels = nullptr, double tolerance = 1e-9);

  /// Check bark/boundary triangle corners only: uv.y == vertex kinetic time * @p uv_height_factor.
  static bool validateAndReportBarkMeshletUvHeights(const std::vector<VoronoiMesh*>& meshlets, const char* scope,
    double uv_height_factor, const std::vector<std::string>* meshlet_labels = nullptr, double tolerance = 1e-9);
};
} // namespace kinDS

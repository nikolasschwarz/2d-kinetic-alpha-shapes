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
  /// OBJ/MTL material name used for triangles incident to inconsistent keyed vertices.
  static constexpr const char* validationErrorMaterialName() { return "red"; }

  /// Route validation output to @p path (append mode). Reopens the file when the path changes.
  static void setLogFile(const std::string& path);

  /// Current validation log path (explicitly set or @ref defaultLogFilePath).
  static const std::string& logFilePath();

  /// Hash-bucket check: vertices sharing the same metadata source identity must agree in stored
  /// world/object-space coordinates. Ignores vertices without a comparable key.
  static MeshVertexSourceValidationResult validateMeshVertexSourcesHaveConsistentWorldPositions(
    const std::vector<VoronoiMesh*>& meshlets, const std::vector<std::string>* meshlet_labels = nullptr,
    double position_tolerance = 1e-5);

  /// Apply @ref validationErrorMaterialName and red vertex colors to triangles incident to @p discrepancies.
  static size_t markMeshVertexSourceDiscrepancies(const std::vector<VoronoiMesh*>& meshlets,
    const std::vector<MeshVertexSourceDiscrepancy>& discrepancies);

  /// @return true when @p mesh uses @ref validationErrorMaterialName on at least one triangle.
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
};
} // namespace kinDS

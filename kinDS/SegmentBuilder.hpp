#pragma once
#include "KineticDelaunay.hpp"
#include "MeshStructure.hpp"
#include "VoronoiMesh.hpp"
#include <functional>
#include <list>
#include <map>
#include <memory>
#include <unordered_map>
#include <vector>

namespace kinDS
{
class SegmentBuilderSectionCallback;
class SegmentBuilderFlipCallback;
class SegmentBuilderRadiusCallback;
class SegmentBuilderCrossingCallback;
class SegmentBuilderSubdivisionCallback;

class SegmentBuilder : public KineticDelaunay::CallbackManager
{
  friend class SegmentBuilderSectionCallback;
  friend class SegmentBuilderFlipCallback;
  friend class SegmentBuilderRadiusCallback;
  friend class SegmentBuilderCrossingCallback;
  friend class SegmentBuilderSubdivisionCallback;

 private:
  // Maps strand IDs to their corresponding segment indices in correct order
  std::vector<std::vector<size_t>> strand_to_segment_indices;
  std::vector<MeshStructure::SegmentProperties> segment_properties; // Properties for each segment mesh
  // Pairs of segments and their corresponding mesh data
  std::vector<MeshStructure::SegmentMeshPair> segment_mesh_pairs;
  std::vector<size_t>
    half_edge_index_to_segment_mesh_pair_index; // Maps edge indices to their corresponding segment mesh pair indices
  std::vector<VoronoiMesh> meshes; // List of all generated meshes

  struct MeshingData
  {
    int mesh_start_vertex_id;
    int mesh_end_vertex_id;
    int start_half_edge_id;
    int end_half_edge_id;
    const void* start_intersection_ref = nullptr;
    const void* end_intersection_ref = nullptr;
    /// Cyclic index around the strand when this meshlet is a closing cap (see `createClosingMesh`).
    int closing_incident_edge_index = -1;
    /// Undirected dual Delaunay edge id (= CrossingData `voronoi_edge_id`) for closing-cap segments; `(size_t)-1` if unset.
    size_t closing_voronoi_edge_id = static_cast<size_t>(-1);
    /// `CrossingData` / `delaunay_edge_intersections` order follows the even directed Delaunay half-edge; when the strand
    /// lies on the odd Voronoi half-edge of this dual edge, walk those lists in the opposite direction along the boundary.
    bool closing_strand_at_voronoi_even_he = true;
  };

  /**
   * @brief Lookup structures built from complete raw segments for closing-cap polygon tracing.
   *
   * @details @p ordered_segments points into @c closing_segments list elements that have both mesh endpoints set.
   * @p start_ref_to_segment maps each crossing `start_intersection_ref` pointer to its index in @p ordered_segments.
   */
  struct ClosingMeshOrderedIndex
  {
    std::vector<MeshingData*> ordered_segments;              ///< Pointers into list storage; valid while list is alive.
    std::unordered_map<const void*, size_t> start_ref_to_segment; ///< CrossingData intersection address -> ordered index.
  };

  /**
   * @brief Output of @ref closingMeshTraceCapPolygons: closed boundary loops and bookkeeping.
   *
   * @details @p polygons are vertex id rings (into the cap mesh) ready to fan-triangulate. @p segment_used records which
   * ordered segments were consumed by the walk. @p mesh_vertex_ids is the id list after extraction plus any vertices
   * added along Delaunay boundary crossings during the trace.
   */
  struct ClosingMeshPolygonsTraceResult
  {
    std::vector<std::vector<size_t>> polygons; ///< One simple polygon per closed walk (mesh vertex indices).
    std::vector<bool> segment_used;             ///< Parallel to ordered_segments: true if that segment was visited.
    std::vector<size_t> mesh_vertex_ids;        ///< Extended cap vertex id trail including boundary intersection verts.
  };

  std::vector<std::list<MeshingData>> segment_mesh_pair_last_left_and_right_vertex;
  // Maps corner indices (correspoding to outgoing half-edge inside the cell) to the index of the cutoff mesh, -1 if no
  // cutoff mesh exists
  std::vector<int> corner_to_cutoff_mesh_indices;
  bool create_transformed_mesh;

  // We no longer use these two factors, they are instead adjusted dynamically in the shader at runtime
  double uv_height_factor = 1.0;
  // 0.1;
  double uv_circum_factor = 1.0;
  // 2.0;
  double texture_diameter = 0.9;

  std::vector<VoronoiMesh> boundary_meshes;

  // for the boundary
  VoronoiMesh boundary_mesh; // Mesh for the boundary cuts
  std::vector<std::pair<size_t, size_t>> boundary_mesh_last_left_and_right_vertex;
  std::vector<size_t> boundary_vertex_to_strand_id;

  // UVs must be adjusted to avoid seams, these are the raw UVs before adjustment
  std::vector<glm::dvec2> boundary_mesh_raw_uvs;

  // Map half-edges to a vertex index in the boundary mesh if a flip created a new boundary edge
  std::vector<int> half_edge_to_boundary_vertex_index;

  KineticDelaunay& kin_del;
  std::unique_ptr<SegmentBuilderSectionCallback> section_callback_;
  std::unique_ptr<SegmentBuilderFlipCallback> flip_callback_;
  std::unique_ptr<SegmentBuilderRadiusCallback> radius_callback_;
  std::unique_ptr<SegmentBuilderCrossingCallback> crossing_callback_;
  std::unique_ptr<SegmentBuilderSubdivisionCallback> subdivision_callback_;
  bool finalized = false; // Flag to indicate if the mesh has been finalized
  bool visual_debug = true; // Always-on visual debug for now (SVG exports)

  glm::dvec3 computeVoronoiVertex(size_t half_edge_id, double t) const;

  void finishMesh(size_t half_edge_id, double t, const std::vector<BoundaryPoint>& boundary_points);

  void startNewMesh(size_t half_edge_id, double t);

  void completeBoundaryMeshSection(size_t he_id, size_t new_left, size_t new_right);

  /**
   * Add a triangle to the boundary mesh. Automatically takes the raw UVs and adjusts them to avoid seams.
   * @param u first vertex index
   * @param v second vertex index
   * @param w third vertex index
   */
  size_t addBoundaryTriangle(size_t u, size_t v, size_t w);

  size_t addBoundaryVertex(glm::dvec3 vertex, glm::dvec2 centroid, size_t strand_id, double t);

  size_t addMeshletTriangle(VoronoiMesh& mesh, size_t u, size_t v, size_t w);

  size_t addMeshletVertex(VoronoiMesh& mesh, const std::vector<BoundaryPoint>& boundary_polygon,
    const glm::dvec2& centroid, glm::dvec3 vertex, size_t strand_id, double t);

  void addVoronoiTriangulationToBoundaryMesh(double t, bool invert_orientation, double offset);

  std::vector<BoundaryPoint> traceConvexHull(double t) const;

  void advanceBoundaryMesh(double t, const std::vector<BoundaryPoint>& boundary_points, const glm::dvec2& centroid);

  void updateBoundary(double t, std::vector<bool>& visited, size_t component_index);

  void updateBoundaries(double t);

  void advanceBoundaryMeshes(double t);

  /**
   * @brief Counts Delaunay edges incident to a strand vertex (graph iterator pass).
   * @param strand_id Delaunay vertex id of the strand.
   * @return Number of incident undirected edges; defines cyclic order for adjacent closing-cap strips.
   */
  size_t closingMeshCountStrandIncidentEdges(size_t strand_id) const;

  /**
   * @brief Appends one vertex to a closing-cap Voronoi meshlet.
   * @param mesh Target cap mesh (modified).
   * @param boundary_polygon Component boundary polygon used for inside tests / UV context.
   * @param centroid Component centroid passed to meshlet vertex creation.
   * @param strand_id Strand owning this cap.
   * @param t Time coordinate for kinetic positions.
   * @param position World-space position of the new vertex.
   * @return New vertex index inside @p mesh (same convention as @ref addMeshletVertex).
   */
  int closingMeshAppendVertex(VoronoiMesh& mesh, const std::vector<BoundaryPoint>& boundary_polygon,
    const glm::dvec2& centroid, size_t strand_id, double t, const glm::dvec3& position);

  /**
   * @brief Finds the CrossingData intersection object for a Voronoi/Delaunay edge pair.
   * @param voronoi_edge_id Undirected dual edge id (half-edge index / 2).
   * @param crossed_delaunay_he_id Directed Delaunay half-edge id whose undirected edge carries the crossing.
   * @return Address of the stored intersection (for `start_intersection_ref` / `end_intersection_ref`), or nullptr.
   */
  const void* closingMeshFindVoronoiEdgeIntersectionPtr(size_t voronoi_edge_id, size_t crossed_delaunay_he_id) const;

  /**
   * @brief 2D position of a Voronoi–Delaunay edge crossing at time @p t (lifted to z = t).
   * @param t Time.
   * @param voronoi_edge_id Undirected Voronoi (dual) edge id.
   * @param delaunay_edge_id Undirected Delaunay edge id.
   * @return Intersection in the xy-plane with z = t, or NaNs if endpoints are degenerate / missing.
   */
  glm::dvec3 closingMeshVoronoiDelaunayCrossingPosition(
    double t, size_t voronoi_edge_id, size_t delaunay_edge_id) const;

  /**
   * @brief Raw inside Voronoi polylines on one dual edge, oriented by @p reverse (no strand lookup).
   *
   * @details When @p reverse is false, the walk starts at the even Voronoi half-edge circumcenter and ends at the odd;
   * when true, odd→even. @ref closingMeshExtractRawSegmentsForVoronoiEdge sets @p reverse from strand incidence.
   * Cap vertices are created only via @p track_vertex.
   */
  std::vector<MeshingData> extractSegmentsForVoronoiEdge(double t, int incident_edge_index, size_t voronoi_edge_id,
    const std::function<int(const glm::dvec3&)>& track_vertex, bool reverse = false);

  /**
   * @brief Extracts raw inside Voronoi polylines for one dual edge incident to the strand (one Voronoi edge id).
   *
   * @details Validates strand incidence, derives orientation, then @ref extractSegmentsForVoronoiEdge.
   *
   * @param incident_edge_index Cyclic index around the strand (matches @ref closingMeshCountStrandIncidentEdges order).
   * @param incident_he Directed Delaunay half-edge id for this incident undirected edge (iterator value).
   */
  std::vector<MeshingData> closingMeshExtractRawSegmentsForVoronoiEdge(size_t strand_id, double t, VoronoiMesh& mesh,
    const std::vector<BoundaryPoint>& boundary_polygon, const glm::dvec2& centroid, int incident_edge_index,
    size_t incident_he, const std::function<int(const glm::dvec3&)>& track_vertex);

  /**
   * @brief Builds ordered segment pointers and a start-ref map from a raw segment list.
   * @param closing_segments In/out list; must stay alive while @p ordered_segments pointers are used.
   * @return Complete segments only, plus map from `start_intersection_ref` to ordered index.
   */
  ClosingMeshOrderedIndex closingMeshBuildOrderedIndex(std::list<MeshingData>& closing_segments);

  /**
   * @brief Logs raw segment count, skipped incomplete entries, and each ordered segment (debug).
   * @param strand_id Strand id (for log prefix).
   * @param t Time (for log prefix).
   * @param closing_segments All raw segments including incomplete.
   * @param ordered_segments Pointers to complete segments only.
   */
  void closingMeshLogDiscoveredSegments(size_t strand_id, double t, const std::list<MeshingData>& closing_segments,
    const std::vector<MeshingData*>& ordered_segments);

  /**
   * @brief Validates ordered segments against the dual edge and CrossingData refs (errors on mismatch).
   * @param t Time used for geometric checks.
   * @param mesh Cap mesh whose vertex positions are compared to refs.
   * @param ordered_segments Segments to validate.
   */
  void closingMeshValidateOrderedSegmentGeometry(double t, const VoronoiMesh& mesh,
    const std::vector<MeshingData*>& ordered_segments);

  /**
   * @brief Walks Voronoi inside legs and Delaunay component boundary to assemble closing polygons.
   *
   * @details For each unused ordered segment as seed, walks along mesh_start→mesh_end on the Voronoi strip, then
   * follows the Delaunay boundary consuming every crossing in list order (with direction adjusted per segment by
   * @ref MeshingData::closing_strand_at_voronoi_even_he), hands off at strand-incident crossings matching
   * @p start_ref_to_segment, and closes loops back to the seed. Emits one vertex ring per closed walk.
   *
   * @param strand_id Strand vertex id.
   * @param t Time.
   * @param num_incident_edges Value from @ref closingMeshCountStrandIncidentEdges (cyclic next-edge checks).
   * @param mesh Cap mesh (modified: new vertices at boundary crossings).
   * @param boundary_polygon Boundary context for new vertices.
   * @param centroid Component centroid.
   * @param mesh_vertex_ids Vertex ids from closing-cap raw extraction; extended with new ids during the walk.
   * @param ordered_segments Complete segments to trace.
   * @param start_ref_to_segment Map for boundary handoff to the next strip.
   * @return Polygon rings, per-segment used flags, and updated vertex-id list.
   */
  ClosingMeshPolygonsTraceResult closingMeshTraceCapPolygons(size_t strand_id, double t, size_t num_incident_edges,
    VoronoiMesh& mesh, const std::vector<BoundaryPoint>& boundary_polygon, const glm::dvec2& centroid,
    std::vector<size_t> mesh_vertex_ids, const std::vector<MeshingData*>& ordered_segments,
    const std::unordered_map<const void*, size_t>& start_ref_to_segment);

  /**
   * @brief Warns about ordered segments that were never visited by the boundary trace.
   * @param strand_id Strand id (log prefix).
   * @param t Time (log prefix).
   * @param ordered_segments All ordered segments.
   * @param segment_used Parallel flags from @ref closingMeshTraceCapPolygons.
   */
  void closingMeshLogUnmatchedOrderedSegments(size_t strand_id, double t,
    const std::vector<MeshingData*>& ordered_segments, const std::vector<bool>& segment_used);

  /**
   * @brief Fan-triangulates each polygon ring into the cap mesh (in-place).
   * @param mesh Cap mesh (triangles appended).
   * @param polygons Vertex index rings; first vertex is the fan pivot for each polygon.
   */
  void closingMeshTriangulatePolygonsFan(VoronoiMesh& mesh, const std::vector<std::vector<size_t>>& polygons);

  size_t createClosingMesh(
    size_t strand_id, double t, const std::vector<BoundaryPoint>& boundary_polygon, const glm::dvec2& centroid);

  void accumulateSegmentProperties();

 public:
  SegmentBuilder(
    KineticDelaunay& kin_del, std::vector<std::pair<size_t, double>> subdivisions, bool create_transformed_mesh);

  SegmentBuilder(KineticDelaunay& kin_del, bool create_transformed_mesh);

  ~SegmentBuilder() override;

  void init() override;

  void finalize(double t) override;

  std::vector<VoronoiMesh> extractMeshes() const;

  std::pair<std::vector<VoronoiMesh>, std::vector<std::vector<int>>> extractSegmentMeshlets(
    bool merge_by_segment = true) const;

  const VoronoiMesh& getBoundaryMesh() const;

  const std::vector<size_t>& getBoundaryVertexToStrandId() const;

  const std::vector<std::vector<size_t>>& getStrandToSegmentIndices() const;

  std::vector<glm::dvec3> computeClampedVoronoiVertices(
    size_t strand_id, double t, const std::vector<BoundaryPoint>& boundary_polygon, const glm::dvec2& centroid);

  void splitComponent(size_t component_id, const std::vector<std::vector<size_t>>& new_components, double t);
};
} // namespace kinDS
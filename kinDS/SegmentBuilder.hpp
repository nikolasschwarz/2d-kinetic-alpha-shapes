#pragma once
#include "KineticDelaunay.hpp"
#include "KineticDelaunayCrossingEvent.hpp"
#include "MeshStructure.hpp"
#include "VoronoiMesh.hpp"
#include <array>
#include <functional>
#include <list>
#include <limits>
#include <map>
#include <memory>
#include <optional>
#include <tuple>
#include <unordered_map>
#include <vector>

namespace kinDS
{
/// Source/target Delaunay edges for one radius 2↔1 boundary transition (passed into intersection meshing only).
struct RadiusBoundaryTransitionShiftContext
{
  bool roles_valid = false;
  std::array<size_t, 2> source_delaunay_edges {};
  size_t target_delaunay_edge = 0;
};

class SegmentBuilderSectionCallback;
class SegmentBuilderFlipCallback;
class SegmentBuilderRadiusCallback;
class SegmentBuilderCrossingCallback;
class SegmentBuilderSubdivisionCallback;

class SegmentBuilder : public KineticDelaunay::CallbackManager
{
 public:
  enum class BoundaryEventType
  {
    Init,
    Section,
    Radius,
    Crossing,
    Subdivision
  };

  enum class BoundarySegmentAction
  {
    NewSegment,
    SegmentCompleted,
    SegmentRemoved,
    SegmentRemapped
  };

  static std::string composeBoundaryMetadata(BoundaryEventType event_type, BoundarySegmentAction segment_action);

  /// JSON for Voronoi-edge strip meshlet vertices (OBJ inline comments when exporting with metadata).
  static std::string composeRegularStripVertexMetadata(double kinetic_time, size_t voronoi_edge_id,
    size_t even_half_edge_id, int strand_even_origin, int strand_odd_origin, BoundaryEventType event_type,
    BoundarySegmentAction segment_action,
    const std::optional<KineticDelaunay::CrossingData::EdgeIntersectionRef>& crossing, const char* pos,
    const char* op = nullptr);

  /// JSON for Voronoi-edge strip meshlet faces (quads emitted as two triangles in @ref finishMesh, etc.).
  static std::string composeRegularStripFaceMetadata(double kinetic_time, size_t voronoi_edge_id,
    size_t even_half_edge_id, int strand_even_origin, int strand_odd_origin, BoundaryEventType event_type,
    BoundarySegmentAction segment_action, const char* op);

  /// When true, radius 2↔1 transitions snap intersection-mesh crossing vertices along the internal Voronoi edge (XY only).
  bool radius_boundary_transition_shift_enabled = false;

  /// Material table for meshlet OBJ export (`material_ids` index into this list).
  static constexpr int RegularMeshletMaterialId = 0;
  static constexpr int BoundaryIntervalMeshletMaterialId = 1;
  static inline const std::vector<std::string> MeshletExportMaterialNames = { "yellow", "brown" };

 private:
  friend class SegmentBuilderSectionCallback;
  friend class SegmentBuilderFlipCallback;
  friend class SegmentBuilderRadiusCallback;
  friend class SegmentBuilderCrossingCallback;
  friend class SegmentBuilderSubdivisionCallback;

  // Maps strand IDs to their corresponding segment indices in correct order
  std::vector<std::vector<size_t>> strand_to_segment_indices;
  std::vector<MeshStructure::SegmentProperties> segment_properties; // Properties for each segment mesh
  // Pairs of segments and their corresponding mesh data
  std::vector<MeshStructure::SegmentMeshPair> segment_mesh_pairs;
  std::vector<MeshStructure::SegmentMeshPair> intersection_segment_mesh_pairs;
  std::vector<MeshStructure::IntersectionMeshPairMetadata> intersection_mesh_pair_metadata;
  std::vector<size_t>
    half_edge_index_to_segment_mesh_pair_index; // Maps edge indices to their corresponding segment mesh pair indices
  std::vector<VoronoiMesh> meshes; // List of all generated meshes
  std::vector<std::string> meshlet_export_suffixes; // Parallel to `meshes`, e.g. "_strandX" / "_voronoiX".

  struct MeshingData
  {
    int mesh_start_vertex_id;
    int mesh_end_vertex_id;
    int start_half_edge_id;
    int end_half_edge_id;
    /// Iterator into `KineticDelaunay::CrossingData::edge_intersections` (via per-edge lists), when this endpoint is a
    /// stored Voronoi–Delaunay crossing.
    std::optional<KineticDelaunay::CrossingData::EdgeIntersectionRef> start_crossing;
    std::optional<KineticDelaunay::CrossingData::EdgeIntersectionRef> end_crossing;
    /// Cyclic index around the strand when this meshlet is a closing cap (see `createClosingMesh`).
    int closing_incident_edge_index = -1;
    /// Undirected dual Delaunay edge id (= CrossingData `voronoi_edge_id`) for closing-cap segments; `(size_t)-1` if unset.
    size_t closing_voronoi_edge_id = static_cast<size_t>(-1);
    /// `CrossingData` / `delaunay_edge_intersections` order follows the even directed Delaunay half-edge; when the strand
    /// lies on the odd Voronoi half-edge of this dual edge, walk those lists in the opposite direction along the boundary.
    bool closing_strand_at_voronoi_even_he = true;
    /// Placeholders on interval start side ("left"); XY filled when the next fixed left endpoint is inserted.
    std::vector<int> flexible_left_vertex_ids;
    /// Placeholders on interval end side ("right"); XY filled when the next fixed right endpoint is inserted.
    std::vector<int> flexible_right_vertex_ids;
  };

  /// One inside-alpha strip on a Voronoi edge: each end is either an open circumcenter (@c *_open_voronoi_half_edge_id)
  /// or a stored boundary crossing (@c *_crossing plus the inside directed Delaunay half-edge id).
  struct RegularMeshStripIntervalEndpoints
  {
    std::optional<KineticDelaunay::CrossingData::EdgeIntersectionRef> start_crossing;
    std::optional<size_t> start_open_voronoi_half_edge_id;
    int start_crossed_inside_half_edge_id = -1;

    std::optional<KineticDelaunay::CrossingData::EdgeIntersectionRef> end_crossing;
    std::optional<size_t> end_open_voronoi_half_edge_id;
    int end_crossed_inside_half_edge_id = -1;
  };

  /**
   * @brief Lookup structures built from complete raw segments for closing-cap polygon tracing.
   *
   * @details @p ordered_segments points into @c closing_segments list elements with both mesh endpoints. Strips may
   * omit boundary `end_crossing` when the strip ends at a circumcenter; the cap walk then joins the next strip whose
   * strand-side endpoint is the same dual Voronoi vertex (circumcenter id: `half_edge.face` on the closing Voronoi
   * edge), not by comparing mesh vertex indices. @p start_crossing_to_segment maps each segment's `start_crossing` iterator
   * to its index in @p ordered_segments (duplicate start crossings are rejected).
   */
  struct ClosingMeshCrossingIteratorHash
  {
    size_t operator()(KineticDelaunay::CrossingData::EdgeIntersectionRef it) const noexcept
    {
      return std::hash<const void*> {}(static_cast<const void*>(std::addressof(*it)));
    }
  };
  struct ClosingMeshCrossingIteratorEq
  {
    bool operator()(KineticDelaunay::CrossingData::EdgeIntersectionRef a,
      KineticDelaunay::CrossingData::EdgeIntersectionRef b) const noexcept
    {
      return a == b;
    }
  };

  struct ClosingMeshOrderedIndex
  {
    std::vector<MeshingData*> ordered_segments; ///< Pointers into list storage; valid while list is alive.
    std::unordered_map<KineticDelaunay::CrossingData::EdgeIntersectionRef, size_t, ClosingMeshCrossingIteratorHash,
      ClosingMeshCrossingIteratorEq>
      start_crossing_to_segment;
  };

  /** Delaunay-boundary strip in canonical interval order (even-half-edge list / increasing parameter). */
  struct BoundaryIntersectionInterval
  {
    size_t voronoi_cell_id = static_cast<size_t>(-1);
    std::optional<KineticDelaunay::CrossingData::EdgeIntersectionRef> start_intersection;
    std::optional<KineticDelaunay::CrossingData::EdgeIntersectionRef> end_intersection;
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
    std::vector<BoundaryIntersectionInterval> traced_boundary_intervals; ///< Boundary intervals in canonical order (see @ref BoundaryIntersectionInterval).
  };

  std::vector<std::list<MeshingData>> segment_mesh_pair_last_left_and_right_vertex;
  // Boundary-interval meshes (built from crossing intersections) are stored separately from regular Voronoi-edge strips.
  std::vector<VoronoiMesh> intersection_meshes;
  std::vector<std::string> intersection_meshlet_export_suffixes;
  std::vector<std::list<MeshingData>> intersection_mesh_pair_last_left_and_right_vertex;
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
  /// When true, one-sided intersection-strip updates append a flexible placeholder on the opposite side (full scheme).
  /// When false (default), ablation: same triangles/endpoints as before flex vectors existed; `MeshingData` flex lists stay empty.
  bool intersection_strip_flexible_vertices_enabled = true;

  glm::dvec3 computeVoronoiVertex(size_t half_edge_id, double t) const;

  std::vector<RegularMeshStripIntervalEndpoints> collectRegularMeshStripIntervalsOnVoronoiEdge(size_t even_half_edge_id,
    size_t voronoi_edge_id, size_t left_containing_tri_id) const;

  MeshingData meshRegularStripInterval(VoronoiMesh& mesh, const std::vector<BoundaryPoint>& boundary_polygon,
    const glm::dvec2& centroid, size_t even_half_edge_id, size_t voronoi_edge_id, double t, int strand_even_origin_i,
    int strand_odd_origin_i, BoundaryEventType event_type, BoundarySegmentAction segment_action,
    const RegularMeshStripIntervalEndpoints& interval,
    const RadiusBoundaryTransitionShiftContext* boundary_transition_shift = nullptr);

  static RegularMeshStripIntervalEndpoints regularMeshStripIntervalFromMeshingData(const MeshingData& segment,
    size_t even_half_edge_id, size_t odd_half_edge_id);

  glm::dvec3 regularMeshStripIntervalEndpointPositionAt(const RegularMeshStripIntervalEndpoints& interval, bool at_start,
    size_t even_half_edge_id, size_t odd_half_edge_id, size_t voronoi_edge_id, double t,
    const RadiusBoundaryTransitionShiftContext* boundary_transition_shift = nullptr) const;

  /// @return `{ new_start_vertex_index, new_end_vertex_index }` — use structured binding, e.g. `auto [left, right] = ...`.
  std::tuple<size_t, size_t> finishRegularMeshStripInterval(VoronoiMesh& mesh, const std::vector<BoundaryPoint>& boundary_polygon,
    const glm::dvec2& centroid, size_t even_half_edge_id, size_t voronoi_edge_id, double t, size_t strand_vertex_id,
    int strand_even_origin_i, int strand_odd_origin_i, BoundaryEventType event_type, BoundarySegmentAction segment_action,
    const RegularMeshStripIntervalEndpoints& interval, size_t last_start_vertex_index, size_t last_end_vertex_index,
    const std::string& finish_face_metadata,
    const RadiusBoundaryTransitionShiftContext* boundary_transition_shift);

  /// True when some crossing on this strip's Voronoi-edge list between @c start_crossing and @c end_crossing (inclusive;
  /// open ends use the list head/tail) has @c voronoi_edge_id equal to @p voronoi_edge_id. Used when @c only_adjacent_segment is set.
  bool regularMeshStripCrossingTouchesVoronoiEdge(const MeshingData& segment, size_t voronoi_edge_id) const;

  /// Extends strips on one Voronoi edge to @p t (quads via @ref finishRegularMeshStripInterval). See implementation for strip state.
  /// @return Copies of @ref MeshingData strips that were extended (empty on early exit or when every strip was skipped).
  std::vector<MeshingData> finishMesh(size_t half_edge_id, double t, const std::vector<BoundaryPoint>& boundary_points,
    BoundaryEventType event_type = BoundaryEventType::Init,
    BoundarySegmentAction segment_action = BoundarySegmentAction::SegmentCompleted,
    const RadiusBoundaryTransitionShiftContext* boundary_transition_shift = nullptr, bool only_adjacent_segment = false);

  /// Seeds strip corner vertices on one Voronoi edge at @p t (no quads). See implementation for strip state.
  /// @return Copies of @ref MeshingData strips that were seeded (empty on early exit or when every interval was skipped).
  std::vector<MeshingData> startNewMesh(size_t half_edge_id, double t, bool reuse_existing_pair_and_mesh = false,
    BoundaryEventType event_type = BoundaryEventType::Init,
    BoundarySegmentAction segment_action = BoundarySegmentAction::NewSegment,
    const RadiusBoundaryTransitionShiftContext* boundary_transition_shift = nullptr, bool only_adjacent_segment = false);
  size_t startNewMeshFromIntersections(size_t voronoi_cell_id, double t,
    std::optional<KineticDelaunay::CrossingData::EdgeIntersectionRef> start_intersection,
    std::optional<KineticDelaunay::CrossingData::EdgeIntersectionRef> end_intersection,
    bool reuse_existing_pair_and_mesh = false, BoundaryEventType event_type = BoundaryEventType::Init,
    BoundarySegmentAction segment_action = BoundarySegmentAction::NewSegment, bool force_single_seed_vertex = false,
    const RadiusBoundaryTransitionShiftContext* boundary_transition_shift = nullptr);
  size_t resolveIntersectionMeshPairIndex(size_t voronoi_cell_id,
    std::optional<KineticDelaunay::CrossingData::EdgeIntersectionRef> start_intersection,
    std::optional<KineticDelaunay::CrossingData::EdgeIntersectionRef> end_intersection,
    double event_time = std::numeric_limits<double>::quiet_NaN()) const;
  void finishMeshFromIntersections(size_t voronoi_cell_id, double t,
    std::optional<KineticDelaunay::CrossingData::EdgeIntersectionRef> start_intersection,
    std::optional<KineticDelaunay::CrossingData::EdgeIntersectionRef> end_intersection,
    BoundaryEventType event_type = BoundaryEventType::Section,
    BoundarySegmentAction segment_action = BoundarySegmentAction::SegmentCompleted,
    const RadiusBoundaryTransitionShiftContext* boundary_transition_shift = nullptr);

  /// If @p crossing_ref is on a radius transition source Delaunay edge, returns an adjacent crossing on the same
  /// Voronoi-edge list whose Delaunay edge is the transition target (vertex position only).
  std::optional<KineticDelaunay::CrossingData::EdgeIntersectionRef> neighborIntersectionOnTargetAlongVoronoiEdge(
    KineticDelaunay::CrossingData::EdgeIntersectionRef crossing_ref, size_t voronoi_edge_id,
    const RadiusBoundaryTransitionShiftContext* boundary_transition_shift) const;

  glm::dvec3 crossingPositionWithRadiusBoundaryTransitionShift(double t,
    KineticDelaunay::CrossingData::EdgeIntersectionRef orig_ref,
    const RadiusBoundaryTransitionShiftContext* boundary_transition_shift) const;

  void logRadiusBoundaryTransitionVertexShift(const char* context, double t,
    KineticDelaunay::CrossingData::EdgeIntersectionRef from_ref, KineticDelaunay::CrossingData::EdgeIntersectionRef to_ref,
    const glm::dvec3& old_pos, const glm::dvec3& new_pos) const;

  static bool delaunayUndirectedEdgeHasVertex(
    const HalfEdgeDelaunayGraph& graph, size_t delaunay_edge_id, size_t vertex_id);

  static std::optional<size_t> oppositeFiniteDelaunayVertexOnUndirectedEdge(
    const HalfEdgeDelaunayGraph& graph, size_t delaunay_edge_id, size_t site_vertex_id);

  /// When @p site_vertex_id is the junction of the two transition source edges, blends shifted corner crossings
  /// (inverse-distance weights from the unshifted site in XY). Otherwise returns @c nullopt.
  std::optional<glm::dvec3> radiusTransitionInterpolatedSitePosition(double t, size_t site_vertex_id,
    size_t strip_delaunay_edge_id, const RadiusBoundaryTransitionShiftContext* boundary_transition_shift) const;

  /**
   * @brief Returns Delaunay-edge intersections in component-boundary traversal order.
   *
   * @details `CrossingData::delaunay_edge_intersections[e]` is stored in canonical even-half-edge order.
   * This helper reorders (possibly reverses) that list so callers can build `[null,first]`, `[k,k+1]`,
   * `[last,null]` intervals consistently in boundary-walk direction.
   */
  std::vector<KineticDelaunay::CrossingData::EdgeIntersectionRef> getBoundaryIntersectionsInBoundaryOrder(
    size_t delaunay_edge_id) const;
  /**
   * @brief Writes one-null interval links by identifying which endpoint vertex is null.
   *
   * @details The null endpoint corresponds to one Delaunay edge origin (even or odd). If it maps to
   * the even-half-edge origin, we write `prev`; otherwise we write `next`.
   *
   * @param intersection_pair_index Pair id to write into crossing link fields.
   * @param null_vertex_id Delaunay vertex id corresponding to the null endpoint.
   * @param ref Non-null crossing reference of the one-null interval.
   * @param interval_is_ref_to_null True for `[ref,null]`, false for `[null,ref]` (debug context).
   * @return true if the pair index was stored in @c prev_segment_mesh_pair_index, false if in @c next_segment_mesh_pair_index.
   */
  bool writeOneNullIntersectionPairLinkByNullVertex(size_t intersection_pair_index, size_t null_vertex_id,
    KineticDelaunay::CrossingData::EdgeIntersectionRef ref, bool interval_is_ref_to_null);
  /**
   * @brief Centralized writer for crossing `prev`/`next` mesh-pair links.
   *
   * @details Handles all interval forms:
   * - `[ref0,ref1]`: writes by actual crossing list order on the Delaunay edge.
   * - one-null: maps `prev`/`next` from whether `voronoi_cell_id` matches the even-half-edge origin on the edge.
   *
   * @param intersection_pair_index Pair id to write into crossing link fields.
   * @param voronoi_cell_id Current interval Voronoi cell id (identifies the open site for one-null intervals).
   * @param start_intersection Interval start crossing (or null endpoint).
   * @param end_intersection Interval end crossing (or null endpoint).
   */
  void writeIntersectionPairLinks(size_t intersection_pair_index, size_t voronoi_cell_id,
    std::optional<KineticDelaunay::CrossingData::EdgeIntersectionRef> start_intersection,
    std::optional<KineticDelaunay::CrossingData::EdgeIntersectionRef> end_intersection);
  size_t determineVoronoiCellForBoundaryIntersectionInterval(size_t delaunay_edge_id,
    std::optional<KineticDelaunay::CrossingData::EdgeIntersectionRef> start_intersection,
    std::optional<KineticDelaunay::CrossingData::EdgeIntersectionRef> end_intersection) const;

  void completeBoundaryMeshSection(size_t he_id, size_t new_left, size_t new_right);

  /**
   * Add a triangle to the boundary mesh. Automatically takes the raw UVs and adjusts them to avoid seams.
   * @param u first vertex index
   * @param v second vertex index
   * @param w third vertex index
   */
  size_t addBoundaryTriangle(size_t u, size_t v, size_t w);

  size_t addBoundaryVertex(glm::dvec3 vertex, glm::dvec2 centroid, size_t strand_id, double t);

  size_t addMeshletTriangle(VoronoiMesh& mesh, size_t u, size_t v, size_t w, const std::string& metadata = "{}",
    int material_id = RegularMeshletMaterialId);

  /**
   * Same as @ref addMeshletTriangle after orienting `(u,v,w)` from combinatorics only: @p inside_boundary_he_id is the
   * *inside* boundary Delaunay half-edge on this strip; its twin is the *outside* boundary half-edge (`index ^ 1`).
   * Apex order relative to strip `(u,v)` is chosen from that outside half-edge index alone (no vertex cross products).
   */
  size_t addBoundaryIntervalTriangleOriented(VoronoiMesh& mesh, size_t u, size_t v, size_t w, int inside_boundary_he_id,
    double t, const std::string& metadata = "{}");

  size_t addMeshletVertex(VoronoiMesh& mesh, const std::vector<BoundaryPoint>& boundary_polygon,
    const glm::dvec2& centroid, glm::dvec3 vertex, size_t strand_id, double t,
    std::optional<size_t> meshlet_voronoi_vertex_for_alpha_check = std::nullopt, const std::string& metadata = "{}",
    const std::optional<glm::dvec3>& debug_color = std::nullopt);

  /// Effective strip corner for triangles: latest flexible on @p left_side, else fixed @c mesh_start / @c mesh_end.
  size_t intersectionStripEffectiveVertexIndex(const MeshingData& seg, bool left_side) const;

  /// When @ref intersection_strip_flexible_vertices_enabled: appends a placeholder (xy = @p centroid, z=@p t, metadata
  /// `pos` left/right) and a wedge triangle `(eff_left, eff_right, flex)` per @ref intersectionStripEffectiveVertexIndex
  /// (snapshot before the new flex is appended). Otherwise no-op.
  void addFlexibleVertexToIntersectionMesh(VoronoiMesh& mesh, MeshingData& seg, bool flexible_on_left_side,
    const std::vector<BoundaryPoint>& boundary_polygon, const glm::dvec2& centroid, size_t strand_id, double t,
    const std::string& metadata = "{}");

  /**
   * @brief Wedge-extend a boundary-interval intersection mesh at a shared Voronoi–Delaunay crossing.
   *
   * @details Uses @p neighbor_pair_idx from the crossing's @c prev_segment_mesh_pair_index (update strip start) or
   * @c next_segment_mesh_pair_index (update strip end). No-op when the pair index is invalid or @p skip_pair_idx matches.
   * When @p append_flexible_placeholder is false, endpoint assignment still runs but no flex placeholder wedge is added.
   */
  void extendIntersectionMeshAtSharedCrossing(size_t neighbor_pair_idx,
    KineticDelaunay::CrossingData::EdgeIntersectionRef shared_ref, bool update_start_on_neighbor, double t,
    BoundaryEventType event_type, BoundarySegmentAction segment_action,
    const RadiusBoundaryTransitionShiftContext* boundary_transition_shift = nullptr,
    bool append_flexible_placeholder = true, std::optional<size_t> skip_pair_idx = std::nullopt);

  /**
   * After a one-sided @c addMeshletVertex + triangle: lerp flexibles on the updated side, assign endpoint,
   * optionally append one flexible on the opposite side. If @p keep_strip_alive is false, only interpolation runs
   * (e.g. strip is about to be cleared).
   */
  void applyIntersectionStripOneSidedFixedVertex(VoronoiMesh& mesh, MeshingData& seg, bool fixed_start_side,
    size_t new_fixed_vertex_index, int inside_half_edge_id,
    const std::optional<KineticDelaunay::CrossingData::EdgeIntersectionRef>& new_crossing_for_updated_side,
    const std::vector<BoundaryPoint>& boundary_polygon, const glm::dvec2& centroid, size_t strand_id, double t,
    bool keep_strip_alive, bool append_flexible_placeholder = true);

  /**
   * For a merge closure vertex (`pos` uniform): both flex chains lerp from their fixed corners to the same
   * @p closure_vertex_index, then lists are cleared. No endpoint reassignment and no new flex.
   */
  void applyIntersectionStripUniformClosureVertex(VoronoiMesh& mesh, MeshingData& seg, size_t closure_vertex_index);

  /// If the containing Delaunay triangle for @p voronoi_vertex_id is not inside the alpha-shape, log a warning with @p position.
  void warnIfVoronoiVertexOutsideAlphaShape(
    const char* context, size_t voronoi_vertex_id, const glm::dvec3& position) const;

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
    const glm::dvec2& centroid, size_t strand_id, double t, const glm::dvec3& position,
    std::optional<size_t> voronoi_vertex_for_alpha_check = std::nullopt);

  /**
   * @brief Finds the CrossingData intersection record for a Voronoi/Delaunay edge pair.
   * @param voronoi_edge_id Undirected dual edge id (half-edge index / 2).
   * @param crossed_delaunay_he_id Directed Delaunay half-edge id whose undirected edge carries the crossing.
   * @return Iterator into `CrossingData::edge_intersections`, or empty if not found.
   */
  std::optional<KineticDelaunay::CrossingData::EdgeIntersectionRef> closingMeshFindVoronoiEdgeIntersection(
    size_t voronoi_edge_id, size_t crossed_delaunay_he_id) const;

  /**
   * Re-resolve `start_crossing` / `end_crossing` from boundary half-edges. Call after `CrossingData` mutates lists
   * (erase/insert), since stored iterators into `edge_intersections` become invalid while `std::optional` may still hold
   * a value.
   */
  void refreshMeshingDataCrossingRefs(MeshingData& seg, size_t voronoi_edge_id) const;

  /// Refresh every segment strip whose dual Voronoi edge is incident to this Voronoi vertex (after a crossing event).
  void refreshStripCrossingRefsIncidentToVoronoiVertex(size_t voronoi_vertex_id);

  /// Refresh all segment strips (after flips / broad `CrossingData` updates).
  void refreshCrossingRefsForAllStrips();

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
   * Cap vertices are created only via @p track_vertex. Pass the Voronoi vertex id (dual circumcenter, `half_edge.face`)
   * as the second argument for strand/end circumcenters so alpha-shape warnings can be logged; use @c std::nullopt for
   * boundary crossing points.
   */
  std::vector<MeshingData> extractSegmentsForVoronoiEdge(double t, int incident_edge_index, size_t voronoi_edge_id,
    const std::function<int(const glm::dvec3&, std::optional<size_t>)>& track_vertex, bool reverse = false);

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
    size_t incident_he, const std::function<int(const glm::dvec3&, std::optional<size_t>)>& track_vertex);

  /**
   * @brief Builds ordered segment pointers and a start-ref map from a raw segment list.
   * @param closing_segments In/out list; must stay alive while @p ordered_segments pointers are used.
   * @return Complete segments only, plus map from `start_crossing` iterator to ordered index.
   */
  ClosingMeshOrderedIndex closingMeshBuildOrderedIndex(std::list<MeshingData>& closing_segments);

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
   * @p start_crossing_to_segment, and closes loops back to the seed. Emits one vertex ring per closed walk.
   *
   * @param strand_id Strand vertex id.
   * @param t Time.
   * @param num_incident_edges Value from @ref closingMeshCountStrandIncidentEdges (cyclic next-edge checks).
   * @param mesh Cap mesh (modified: new vertices at boundary crossings).
   * @param boundary_polygon Boundary context for new vertices.
   * @param centroid Component centroid.
   * @param mesh_vertex_ids Vertex ids from closing-cap raw extraction; extended with new ids during the walk.
   * @param ordered_segments Complete segments to trace.
   * @param start_crossing_to_segment Map for boundary handoff to the next strip.
   * @return Polygon rings, per-segment used flags, and updated vertex-id list.
   */
  ClosingMeshPolygonsTraceResult closingMeshTraceCapPolygons(size_t strand_id, double t, size_t num_incident_edges,
    VoronoiMesh& mesh, const std::vector<BoundaryPoint>& boundary_polygon, const glm::dvec2& centroid,
    std::vector<size_t> mesh_vertex_ids, const std::vector<MeshingData*>& ordered_segments,
    const std::unordered_map<KineticDelaunay::CrossingData::EdgeIntersectionRef, size_t, ClosingMeshCrossingIteratorHash,
      ClosingMeshCrossingIteratorEq>& start_crossing_to_segment);

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

  size_t createClosingMesh(size_t strand_id, double t, const std::vector<BoundaryPoint>& boundary_polygon,
    const glm::dvec2& centroid, std::vector<BoundaryIntersectionInterval>* traced_boundary_intervals = nullptr);

  void accumulateSegmentProperties();

  /// Registers a newly created meshlet and stores export suffix metadata.
  /// @param creation_kinetic_time If finite, stored on the mesh for export/debug (see @ref VoronoiMesh::setCreationKineticTime).
  size_t registerMeshletWithSuffix(VoronoiMesh&& mesh, std::string suffix,
    double creation_kinetic_time = std::numeric_limits<double>::quiet_NaN());

   public:
  /// One-line Voronoi meshlet diagnostics: dual edge id, pair slot, verts, tris, strip counts (@p extra_note optional).
  void meshletDiagnosticLogLine(const char* tag, size_t half_edge_id, double t, const char* extra_note = "") const;

  /// After @ref startNewMesh strip build: warn if topology/metadata suggests a non-empty mesh but vertices are missing.
  void meshletDiagnosticWarnIfUnexpectedEmptyAfterStartNewMesh(
    size_t half_edge_even, double t, bool initial_left_inside, const VoronoiMesh& mesh, const std::list<MeshingData>& strips) const;


  SegmentBuilder(
    KineticDelaunay& kin_del, std::vector<std::pair<size_t, double>> subdivisions, bool create_transformed_mesh);

  SegmentBuilder(KineticDelaunay& kin_del, bool create_transformed_mesh);

  ~SegmentBuilder() override;

  void init() override;

  void finalize(double t) override;

  std::vector<VoronoiMesh> extractMeshes() const;

  std::pair<std::vector<VoronoiMesh>, std::vector<std::vector<int>>> extractSegmentMeshlets(
    bool merge_by_segment = true) const;

  std::vector<std::string> extractSegmentMeshletExportSuffixes(bool merge_by_segment = true) const;
  std::vector<VoronoiMesh> extractBoundaryIntervalMeshlets() const;
  std::vector<std::string> extractBoundaryIntervalMeshletExportSuffixes() const;


  const VoronoiMesh& getBoundaryMesh() const;

  const std::vector<size_t>& getBoundaryVertexToStrandId() const;

  const std::vector<std::vector<size_t>>& getStrandToSegmentIndices() const;

  std::vector<glm::dvec3> computeClampedVoronoiVertices(
    size_t strand_id, double t, const std::vector<BoundaryPoint>& boundary_polygon, const glm::dvec2& centroid);

  void splitComponent(size_t component_id, const std::vector<std::vector<size_t>>& new_components, double t);
};
} // namespace kinDS
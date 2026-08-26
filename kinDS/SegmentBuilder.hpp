#pragma once
#include "KineticDelaunay.hpp"
#include "KineticDelaunayCrossingEvent.hpp"
#include "MeshStructure.hpp"
#include "VoronoiMesh.hpp"
#include <array>
#include <functional>
#include <limits>
#include <list>
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
  /// When set, at least one registered Voronoi vertex lies in the affected triangle. Shift must not run; use the
  /// no-shift traced cell / triangle-cap path (same outcome as the pending branch-split fallback).
  std::optional<size_t> interior_voronoi_vertex_id {};
};

/// True when radius boundary-transition vertex shift should drive meshing (2↔1 roles valid, no interior VV).
inline bool radiusBoundaryTransitionShiftApplicable(
  bool shift_enabled, const RadiusBoundaryTransitionShiftContext& ctx)
{
  return shift_enabled && ctx.roles_valid && !ctx.interior_voronoi_vertex_id.has_value();
}

/// Mesh-space segment used to place an intersection vertex:
/// @c endpoint0 + @c param * (@c endpoint1 - @c endpoint0), where each endpoint is placed like @c source @c site.
struct IntersectionInterpolationDebug
{
  glm::dvec3 endpoint0 {};
  glm::dvec3 endpoint1 {};
  double param = std::numeric_limits<double>::quiet_NaN();

  bool valid() const
  {
    return std::isfinite(param) && std::isfinite(endpoint0.x) && std::isfinite(endpoint0.y)
      && std::isfinite(endpoint0.z) && std::isfinite(endpoint1.x) && std::isfinite(endpoint1.y)
      && std::isfinite(endpoint1.z);
  }
};

/// One endpoint of a radius-transition projection segment. Exactly one source should be set.
struct RadiusTransitionProjectionEndpoint
{
  std::optional<KineticDelaunay::CrossingData::EdgeIntersectionRef> intersection {};
  std::optional<size_t> site_vertex_id {};
  std::optional<size_t> voronoi_vertex_id {};
};

/// Semantic segment plus parameter used to reconstruct the same projected point in profile or mesh space.
struct RadiusTransitionProjection
{
  std::array<RadiusTransitionProjectionEndpoint, 2> endpoints {};
  double param = std::numeric_limits<double>::quiet_NaN();
};

/// Conceptual crossing (topology / mesh-pair links) vs profile position source during a radius boundary transition.
struct RadiusBoundaryTransitionCrossingPlacement
{
  KineticDelaunay::CrossingData::EdgeIntersectionRef conceptual_intersection {};
  KineticDelaunay::CrossingData::EdgeIntersectionRef position_intersection {};
  /// Legacy fallback for synthetic shifts that cannot be represented semantically.
  std::optional<glm::dvec3> explicit_profile_position {};
  /// When set, mesh placement should use Voronoi barycentric transfer for this vertex (unshifted mesh sites).
  std::optional<size_t> snap_voronoi_vertex_id {};
  /// Semantic projection used to reconstruct synthetic shifts in either profile or mesh space.
  std::optional<RadiusTransitionProjection> projection {};

  bool positionDiffersFromConceptual() const
  {
    return explicit_profile_position.has_value() || projection.has_value()
      || position_intersection != conceptual_intersection;
  }
};

class SegmentBuilderSectionCallback;
class SegmentBuilderFlipCallback;
class SegmentBuilderRadiusCallback;
class SegmentBuilderCrossingCallback;
class SegmentBuilderSubdivisionCallback;
class SegmentBuilderSeparationCallback;

class SegmentBuilder : public KineticDelaunay::CallbackManager
{
 public:
  class MetadataBuilder
  {
   public:
    static MetadataBuilder fromObject(const std::string& metadata);

    MetadataBuilder& addString(const char* key, const char* value);
    MetadataBuilder& addString(const char* key, const std::string& value);
    MetadataBuilder& addSize(const char* key, size_t value);
    MetadataBuilder& addInt(const char* key, int value);
    MetadataBuilder& addDouble(const char* key, double value);
    MetadataBuilder& addBool(const char* key, bool value);
    MetadataBuilder& ensureSize(const char* key, size_t value);
    MetadataBuilder& ensureDouble(const char* key, double value);
    MetadataBuilder& addRaw(const char* key, const std::string& raw_json_value);
    MetadataBuilder& ensureString(const char* key, const char* value);
    MetadataBuilder& ensureBool(const char* key, bool value);
    std::string build() const;

   private:
    bool hasKey(const char* key) const;
    std::string raw_body_;
    std::vector<std::pair<std::string, std::string>> fields_;
  };

  class ScopedMetadataCallbackPhase
  {
   public:
    ScopedMetadataCallbackPhase(SegmentBuilder& segment_builder, const char* callback_phase);
    ~ScopedMetadataCallbackPhase();

   private:
    SegmentBuilder& segment_builder_;
    std::string previous_phase_;
  };

  /// Clears the crossing-event vertex buffer when @c afterEvent returns (including early exits).
  class ScopedEndActiveCrossingEvent
  {
   public:
    explicit ScopedEndActiveCrossingEvent(SegmentBuilder& segment_builder);
    ~ScopedEndActiveCrossingEvent();

   private:
    SegmentBuilder& segment_builder_;
  };

  enum class BoundaryEventType
  {
    Init,
    Section,
    Radius,
    Crossing,
    Subdivision,
    Separation
  };

  enum class BoundarySegmentAction
  {
    NewSegment,
    SegmentCompleted,
    SegmentRemoved,
    SegmentRemapped
  };

  std::string composeBoundaryMetadata(BoundaryEventType event_type, BoundarySegmentAction segment_action) const;

  /// JSON for Voronoi-edge strip meshlet vertices (OBJ inline comments when exporting with metadata).
  std::string composeRegularStripVertexMetadata(double kinetic_time, size_t voronoi_edge_id,
    size_t even_half_edge_id, int strand_even_origin, int strand_odd_origin, BoundaryEventType event_type,
    BoundarySegmentAction segment_action,
    const std::optional<KineticDelaunay::CrossingData::EdgeIntersectionRef>& crossing, const char* pos,
    const char* op = nullptr, const char* source = nullptr,
    const std::optional<KineticDelaunay::CrossingData::EdgeIntersectionRef>& position_crossing = std::nullopt,
    bool radius_shift_explicit_profile_position = false) const;

  static std::string intersectionCrossingVertexMetadata(const std::string& base_metadata,
    KineticDelaunay::CrossingData::EdgeIntersectionRef conceptual_ref,
    KineticDelaunay::CrossingData::EdgeIntersectionRef position_ref, const char* pos_label,
    bool radius_shift_explicit_profile_position = false);

  /// JSON for Voronoi-edge strip meshlet faces (quads emitted as two triangles in @ref finishMesh, etc.).
  std::string composeRegularStripFaceMetadata(double kinetic_time, size_t voronoi_edge_id, size_t even_half_edge_id,
    int strand_even_origin, int strand_odd_origin, BoundaryEventType event_type, BoundarySegmentAction segment_action,
    const char* op) const;

  /// JSON for boundary-mesh (bark/interior) triangles.
  std::string composeBoundaryMeshFaceMetadata(double kinetic_time, const char* mesh_type,
    size_t half_edge_id = static_cast<size_t>(-1), size_t delaunay_face_id = static_cast<size_t>(-1),
    size_t input_branch_id = static_cast<size_t>(-1)) const;

  /// JSON for closing-cap meshlet triangles.
  std::string composeClosingMeshFaceMetadata(double kinetic_time, size_t strand_id) const;

  void configureMeshletStorage(VoronoiMesh& mesh) const;

  /// Index into @c intersection_meshes for a boundary-interval meshlet, if @p mesh is one of them.
  std::optional<size_t> intersectionMeshletIndexForMesh(const VoronoiMesh& mesh) const;

  /// When true, radius 2↔1 transitions use boundary-transition vertex shift when a 2↔1 boundary-edge correspondence
  /// exists and no interior Voronoi vertex lies in the triangle; otherwise traced Voronoi-cell meshlets are used.
  bool radius_boundary_transition_shift_enabled = true;
  /// When true, reject radius shifting and use the triangle-closing fallback if a pending split's triangle spans
  /// multiple future runtime branches.
  bool radius_pending_split_triangle_fallback_enabled = true;
  /// When false, skip building/storing per-vertex and per-face JSON metadata on meshlets.
  bool store_mesh_metadata = false;
  /// When true (e.g. via @c --validate), run mesh vertex source consistency checks in @ref finalize.
  bool validate_mesh_vertex_sources = false;
  /// When false, skip meshlet diagnostic logging and related string assembly.
  bool diagnostics = false;
  /// When true, finalize collapses degree-2 flexible vertices on intersection meshlets.
  bool collapse_degree_two_flexible_vertices_postprocess_enabled = true;
  /// When true, segment meshlet assembly aligns glue seams across contributing mesh copies
  /// (multi-match boundary seams; propagate unmatched interiors onto partners before combine).
  bool align_flexible_glue_edges_postprocess_enabled = true;
  /// Emit closing-cap meshlets at kinetic @c start_section / premature @c end_section.
  /// Ending input branches (including tree top) always get caps at finalize; @c mesh_cap_at_end only
  /// seals strands truncated by @c --end.
  bool mesh_cap_at_start = false;
  bool mesh_cap_at_end = false;
  /// Failure SVG/TXT dumps (ring-walk / triangulate FAIL, …) plus common highlighted kinetic SVG.
  /// CLI --error-files; also implied by visual debug.
  void setErrorFiles(bool enabled) { error_files = enabled; }
  bool errorFilesEnabled() const { return error_files; }
  bool shouldDumpErrorFiles() const { return error_files || visual_debug; }

  /// Material table for meshlet OBJ export (`material_ids` index into this list).
  static constexpr int RegularMeshletMaterialId = 0;
  static constexpr int BoundaryIntervalMeshletMaterialId = 1;
  static constexpr int PendingSplitFallbackMeshletMaterialId = 2;
  static inline const std::vector<std::string> MeshletExportMaterialNames = { "green", "brown", "light_blue" };

 private:
  friend class SegmentBuilderSectionCallback;
  friend class SegmentBuilderFlipCallback;
  friend class SegmentBuilderRadiusCallback;
  friend class SegmentBuilderCrossingCallback;
  friend class SegmentBuilderSubdivisionCallback;
  friend class SegmentBuilderSeparationCallback;

  // Maps strand IDs to their corresponding segment indices in correct order
  std::vector<std::vector<size_t>> strand_to_segment_indices;
  std::vector<MeshStructure::SegmentProperties> segment_properties; // Properties for each segment mesh
  // Pairs of segments and their corresponding mesh data
  std::vector<MeshStructure::SegmentMeshPair> segment_mesh_pairs;
  std::vector<MeshStructure::SegmentMeshPair> intersection_segment_mesh_pairs;
  std::vector<MeshStructure::IntersectionMeshPairMetadata> intersection_mesh_pair_metadata;
  std::string metadata_callback_phase_;
  std::optional<size_t> active_crossing_voronoi_vertex_id_;
  std::optional<size_t> active_crossing_delaunay_edge_id_;
  std::array<size_t, 3> active_crossing_voronoi_edge_ids_ {
    static_cast<size_t>(-1), static_cast<size_t>(-1), static_cast<size_t>(-1)
  };
  std::optional<glm::dvec3> active_crossing_mesh_position_;
  std::optional<glm::dvec2> active_crossing_delaunay_xy_;
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
    /// Undirected dual Delaunay edge id (= CrossingData `voronoi_edge_id`) for closing-cap segments; `(size_t)-1` if
    /// unset.
    size_t closing_voronoi_edge_id = static_cast<size_t>(-1);
    /// `CrossingData` / `delaunay_edge_intersections` order follows the even directed Delaunay half-edge; when the
    /// strand lies on the odd Voronoi half-edge of this dual edge, walk those lists in the opposite direction along the
    /// boundary.
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

  /// One oriented boundary crossing along a Voronoi edge from `CrossingData::voronoi_edge_intersections`.
  struct DirectedVoronoiEdgeCrossing
  {
    size_t crossed_half_edge_id;
    KineticDelaunay::CrossingData::EdgeIntersectionRef ref;
  };

  std::vector<DirectedVoronoiEdgeCrossing> orientCrossingsAlongVoronoiEdge(
    size_t voronoi_edge_id, size_t start_containing_tri_id, bool reverse_traversal,
    size_t even_half_edge_id_for_diagnostics = static_cast<size_t>(-1)) const;

  /**
   * @brief Lookup structures built from complete raw segments for closing-cap polygon tracing.
   *
   * @details @p ordered_segments points into @c closing_segments list elements with both mesh endpoints. Strips may
   * omit boundary `end_crossing` when the strip ends at a circumcenter; the cap walk then joins the next strip whose
   * strand-side endpoint is the same dual Voronoi vertex (circumcenter id: `half_edge.face` on the closing Voronoi
   * edge), not by comparing mesh vertex indices. @p start_crossing_to_segment maps each segment's `start_crossing`
   * iterator to its index in @p ordered_segments (duplicate start crossings are rejected).
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
   * @details @p polygons are vertex id rings (into the cap mesh) ready to triangulate. @p segment_used records
   * which ordered segments were consumed by the walk. @p mesh_vertex_ids is the id list after extraction plus any
   * vertices added along Delaunay boundary crossings during the trace.
   */
  struct ClosingMeshPolygonsTraceResult
  {
    std::vector<std::vector<size_t>> polygons; ///< One simple polygon per closed walk (mesh vertex indices).
    std::vector<bool> segment_used; ///< Parallel to ordered_segments: true if that segment was visited.
    std::vector<size_t> mesh_vertex_ids; ///< Extended cap vertex id trail including boundary intersection verts.
    std::vector<BoundaryIntersectionInterval>
      traced_boundary_intervals; ///< Boundary intervals in canonical order (see @ref BoundaryIntersectionInterval).
  };

  std::vector<std::list<MeshingData>> segment_mesh_pair_last_left_and_right_vertex;
  // Boundary-interval meshes (built from crossing intersections) are stored separately from regular Voronoi-edge
  // strips.
  std::vector<VoronoiMesh> intersection_meshes;
  std::vector<std::string> intersection_meshlet_export_suffixes;
  std::vector<std::list<MeshingData>> intersection_mesh_pair_last_left_and_right_vertex;
  /// Parallel to @c meshes (interior / regular Voronoi-edge strips) and @c intersection_meshes
  /// (boundary-interval strips). Set when a meshlet is retired (currently: subdivision of its cell).
  std::vector<bool> interior_meshlet_completed_;
  std::vector<bool> boundary_meshlet_completed_;

  std::optional<size_t> regularMeshletIndexForMesh(const VoronoiMesh& mesh) const;
  /// Warns when @p mesh is a retired interior or boundary meshlet. Returns true so callers can skip the mutation.
  bool warnAndSkipIfMeshletCompleted(const VoronoiMesh& mesh, const char* operation, double t) const;
  bool isBoundaryMeshletCompleted(size_t meshlet_index) const;
  /// If @p pair_idx names a meshlet retired by subdivision, return the newest live pair with the same
  /// interval metadata (cell + start/end Delaunay edges). Otherwise return @p pair_idx.
  size_t preferLiveBoundaryMeshPair(size_t pair_idx) const;
  void markInteriorMeshletCompleted(size_t meshlet_index);
  void markBoundaryMeshletCompleted(size_t meshlet_index);
  /// Retire every interior strip on Voronoi edges of @p voronoi_cell_id (Delaunay site / strand).
  void markInteriorMeshletsCompletedForVoronoiCell(size_t voronoi_cell_id);
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

  // for the boundary
  VoronoiMesh boundary_mesh; // Mesh for the boundary cuts
  std::vector<std::pair<size_t, size_t>> boundary_mesh_last_left_and_right_vertex;
  std::vector<size_t> boundary_vertex_to_strand_id;

  // UVs must be adjusted to avoid seams, these are the raw UVs before adjustment
  std::vector<glm::dvec2> boundary_mesh_raw_uvs;
  /// Parallel to @ref intersection_meshes: unscaled polar raw UV per boundary-interval vertex (angle/2pi, t).
  std::vector<std::vector<glm::dvec2>> intersection_mesh_raw_uvs;

  // Map half-edges to a vertex index in the boundary mesh if a flip created a new boundary edge
  std::vector<int> half_edge_to_boundary_vertex_index;

  KineticDelaunay& kin_del;
  std::function<void(size_t, std::function<void(size_t)>)> parallel_for;
  std::unique_ptr<SegmentBuilderSectionCallback> section_callback_;
  std::unique_ptr<SegmentBuilderFlipCallback> flip_callback_;
  std::unique_ptr<SegmentBuilderRadiusCallback> radius_callback_;
  std::unique_ptr<SegmentBuilderCrossingCallback> crossing_callback_;
  std::unique_ptr<SegmentBuilderSubdivisionCallback> subdivision_callback_;
  std::unique_ptr<SegmentBuilderSeparationCallback> separation_callback_;
  bool finalized = false; // Flag to indicate if the mesh has been finalized
  bool visual_debug = false; // Full visual-debug SVG/TXT exports (CLI --debug-files)
  /// Failure SVG/TXT dumps; toggled via @ref setErrorFiles / implied by @ref visual_debug.
  bool error_files = false;
  /// When true, one-sided intersection-strip updates append a flexible placeholder on the opposite side (full scheme).
  /// When false (default), ablation: same triangles/endpoints as before flex vectors existed; `MeshingData` flex lists
  /// stay empty.
  bool intersection_strip_flexible_vertices_enabled = true;

  glm::dvec3 transformFromInputBranchToObjectSpace(glm::dvec3 vertex, size_t strand_id, double t) const;
  /// Site position for meshing: never includes kinetic separation.
  /// When @c create_transformed_mesh is false (CLI @c --untransformed), local profile coordinates (no
  /// reference-branch remapping); when true, object space. Prefer this over raw
  /// @ref KineticDelaunay::getPointAt for mesh vertices. Kinetic SVG / event geometry still uses
  /// @ref KineticDelaunay::getPointInDelaunaySpace (reference-branch frame + separation).
  glm::dvec3 getPointInMeshSpace(size_t strand_id, double t) const;
  /// Site → stored mesh coordinate, matching @c addMeshletVertex for @c source @c site (@ref getPointInMeshSpace).
  glm::dvec3 computeMeshSiteVertexPosition(glm::dvec3 profile_vertex, size_t strand_id, double t) const;
  /// Crossing position for meshing: never includes kinetic separation; refreshes stale @c delaunay_edge_param at @p t.
  glm::dvec3 getCrossingCoordsInMeshSpace(KineticDelaunay::CrossingData::EdgeIntersectionRef intersection,
    double t) const;
  struct MeshletVertexRuntimeInfo
  {
    bool is_flexible_placeholder = false;
    bool radius_shift_explicit_profile_position = false;
    std::optional<KineticDelaunay::CrossingData::EdgeIntersectionRef> position_intersection;
    std::optional<KineticDelaunay::CrossingData::EdgeIntersectionRef> conceptual_intersection;
    std::optional<RadiusTransitionProjection> radius_transition_projection;
    /// When set, stored as @ref VoronoiMesh::setProfilePlanePosition for ear-clip triangulation.
    /// Use the same XY the caller used to build the polygon ring (typically Delaunay space).
    std::optional<glm::dvec2> triangulation_plane_xy {};
    /// When set, skip recomputing mesh placement and use these coordinates (flip coincidence buffer, or
    /// radius-shifted corner reuse via @c RadiusShiftedSiteCacheEntry::explicit_mesh_position).
    /// Crossing-event coincidence uses @c SegmentBuilder's event buffer instead of this field.
    std::optional<glm::dvec3> explicit_mesh_position {};
    std::optional<glm::dvec2> explicit_delaunay_xy {};

    bool isIntersectionVertex() const { return position_intersection.has_value(); }
  };
  struct MeshIntersectionObjectSpaceResult
  {
    glm::dvec3 position {};
    std::optional<IntersectionInterpolationDebug> mesh_interpolation {};
  };

  MeshIntersectionObjectSpaceResult computeMeshIntersectionObjectSpace(
    const MeshletVertexRuntimeInfo& runtime_info, glm::dvec3 fallback_profile_vertex, size_t fallback_strand_id,
    double t) const;
  struct MeshVoronoiVertexObjectSpaceResult
  {
    glm::dvec3 position {};
    /// Mesh-space positions of the three *containing* triangle sites (not necessarily the dual triangle).
    std::array<std::pair<size_t, glm::dvec3>, 3> site_mesh_positions {};
    size_t site_count = 0;
    std::optional<size_t> containing_tri_id {};
    /// True when the containing Delaunay triangle has an infinite vertex (barycentric transfer is ill-defined).
    std::optional<bool> containing_tri_infinite {};
    std::optional<std::array<double, 3>> barycentric {};
    /// True when the active crossing-event buffer supplied this position (not barycentric transfer).
    bool from_crossing_event_buffer = false;
  };
  /// Place a Voronoi vertex by barycentric transfer: compute it and the containing triangle in Delaunay space
  /// *with* kinetic separation shifts, then interpolate the same sites in mesh space *without* those shifts.
  /// If the containing triangle is infinite, the vertex lies on the finite edge and is interpolated from those
  /// two sites (same as an intersection), with barycentric weight 0 on the infinite vertex.
  /// During a crossing event, the event Voronoi vertex reuses the position buffered at @c beforeEvent.
  MeshVoronoiVertexObjectSpaceResult computeMeshVoronoiVertexObjectSpace(size_t voronoi_vertex_id,
    glm::dvec3 fallback_profile_vertex, size_t fallback_strand_id, double t) const;
  /// Intersection of an incident Voronoi edge of @p voronoi_vertex_id with @p delaunay_edge_id, if any.
  std::optional<KineticDelaunay::CrossingData::EdgeIntersectionRef> findIncidentIntersectionOnDelaunayEdge(
    size_t voronoi_vertex_id, size_t delaunay_edge_id) const;
  /// Compute one mesh/Delaunay position for this crossing (prefer a pre-event incident intersection) and
  /// reuse it for the event Voronoi vertex and for added/removed coincidence intersections.
  void beginActiveCrossingEvent(size_t voronoi_vertex_id, size_t delaunay_edge_id, double t, size_t strand_id);
  void endActiveCrossingEvent();
  bool activeCrossingEventBufferApplies(std::optional<size_t> voronoi_vertex_id,
    const std::optional<KineticDelaunay::CrossingData::EdgeIntersectionRef>& intersection) const;
  glm::dvec2 meshVirtualShiftForStrand(size_t strand_id, double t) const;
  void applyMeshVirtualShiftToProfileVertex(
    glm::dvec3& vertex, glm::dvec2& profile_xy, size_t strand_id, double t, bool& includes_virtual_shift) const;
  /// Y/Z swap so profile-space meshes view with the same XY layout as SVG.
  void applyUntransformedMeshViewTransform();

  glm::dvec3 computeVoronoiVertex(size_t half_edge_id, double t) const;

  std::vector<RegularMeshStripIntervalEndpoints> collectRegularMeshStripIntervalsOnVoronoiEdge(
    size_t even_half_edge_id, size_t voronoi_edge_id, size_t left_containing_tri_id) const;

  MeshingData meshRegularStripInterval(VoronoiMesh& mesh, const std::vector<BoundaryPoint>& boundary_polygon,
    const glm::dvec2& centroid, size_t even_half_edge_id, size_t voronoi_edge_id, double t, int strand_even_origin_i,
    int strand_odd_origin_i, BoundaryEventType event_type, BoundarySegmentAction segment_action,
    const RegularMeshStripIntervalEndpoints& interval,
    const RadiusBoundaryTransitionShiftContext* boundary_transition_shift = nullptr);

  static RegularMeshStripIntervalEndpoints regularMeshStripIntervalFromMeshingData(
    const MeshingData& segment, size_t even_half_edge_id, size_t odd_half_edge_id);

  /// @return `{ new_start_vertex_index, new_end_vertex_index }` — use structured binding, e.g. `auto [left, right] =
  /// ...`.
  std::tuple<size_t, size_t> finishRegularMeshStripInterval(VoronoiMesh& mesh,
    const std::vector<BoundaryPoint>& boundary_polygon, const glm::dvec2& centroid, size_t even_half_edge_id,
    size_t voronoi_edge_id, double t, size_t strand_vertex_id, int strand_even_origin_i, int strand_odd_origin_i,
    BoundaryEventType event_type, BoundarySegmentAction segment_action,
    const RegularMeshStripIntervalEndpoints& interval, size_t last_start_vertex_index, size_t last_end_vertex_index,
    const std::string& finish_face_metadata, const RadiusBoundaryTransitionShiftContext* boundary_transition_shift);

  /// True when some crossing on this strip's Voronoi-edge list between @c start_crossing and @c end_crossing
  /// (inclusive; open ends use the list head/tail) has @c voronoi_edge_id equal to @p voronoi_edge_id. Used when @c
  /// only_adjacent_segment is set.
  bool regularMeshStripCrossingTouchesVoronoiEdge(const MeshingData& segment, size_t voronoi_edge_id) const;

  /// Extends strips on one Voronoi edge to @p t (quads via @ref finishRegularMeshStripInterval). See implementation for
  /// strip state.
  /// @return Copies of @ref MeshingData strips that were extended (empty on early exit or when every strip was
  /// skipped).
  std::vector<MeshingData> finishMesh(size_t half_edge_id, double t, const std::vector<BoundaryPoint>& boundary_points,
    BoundaryEventType event_type = BoundaryEventType::Init,
    BoundarySegmentAction segment_action = BoundarySegmentAction::SegmentCompleted,
    const RadiusBoundaryTransitionShiftContext* boundary_transition_shift = nullptr,
    bool only_adjacent_segment = false);

  /// Seeds strip corner vertices on one Voronoi edge at @p t (no quads). See implementation for strip state.
  /// @return Copies of @ref MeshingData strips that were seeded (empty on early exit or when every interval was
  /// skipped).
  std::vector<MeshingData> startNewMesh(size_t half_edge_id, double t, bool reuse_existing_pair_and_mesh = false,
    BoundaryEventType event_type = BoundaryEventType::Init,
    BoundarySegmentAction segment_action = BoundarySegmentAction::NewSegment,
    const RadiusBoundaryTransitionShiftContext* boundary_transition_shift = nullptr,
    bool only_adjacent_segment = false);
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
  /// Same as the resolve overload, but uses a caller-known @p intersection_pair_index (no re-resolve).
  void finishMeshFromIntersections(size_t voronoi_cell_id, double t,
    std::optional<KineticDelaunay::CrossingData::EdgeIntersectionRef> start_intersection,
    std::optional<KineticDelaunay::CrossingData::EdgeIntersectionRef> end_intersection, BoundaryEventType event_type,
    BoundarySegmentAction segment_action, const RadiusBoundaryTransitionShiftContext* boundary_transition_shift,
    size_t intersection_pair_index);

  /// Adjacent crossing on the same Voronoi-edge list whose Delaunay edge is @p target_delaunay_edge.
  /// Uses @p crossing_ref's cached @c voronoi_ref; caller must only invoke for source-edge intersections.
  std::optional<KineticDelaunay::CrossingData::EdgeIntersectionRef> neighborIntersectionOnTargetAlongVoronoiEdge(
    KineticDelaunay::CrossingData::EdgeIntersectionRef crossing_ref, size_t target_delaunay_edge) const;

  RadiusBoundaryTransitionCrossingPlacement resolveRadiusBoundaryTransitionCrossingPlacement(double t,
    KineticDelaunay::CrossingData::EdgeIntersectionRef conceptual_ref,
    const RadiusBoundaryTransitionShiftContext* boundary_transition_shift) const;

  glm::dvec3 crossingProfilePosition(double t,
    KineticDelaunay::CrossingData::EdgeIntersectionRef intersection_ref) const;

  glm::dvec3 crossingProfilePositionFromPlacement(
    double t, const RadiusBoundaryTransitionCrossingPlacement& placement) const;

  void logRadiusBoundaryTransitionVertexShift(const char* context, double t,
    KineticDelaunay::CrossingData::EdgeIntersectionRef from_ref,
    KineticDelaunay::CrossingData::EdgeIntersectionRef to_ref, const glm::dvec3& old_pos,
    const glm::dvec3& new_pos) const;

  static bool delaunayUndirectedEdgeHasVertex(
    const HalfEdgeDelaunayGraph& graph, size_t delaunay_edge_id, size_t vertex_id);

  static std::optional<size_t> oppositeFiniteDelaunayVertexOnUndirectedEdge(
    const HalfEdgeDelaunayGraph& graph, size_t delaunay_edge_id, size_t site_vertex_id);

  /// When @p site_vertex_id is the junction of the two transition source edges, either snaps to the interior
  /// Voronoi vertex or projects the site onto the segment between the two unshifted anchors. The returned semantic
  /// endpoints and parameter reconstruct the shifted point independently in profile and mesh space.
  struct RadiusTransitionSitePlacement
  {
    glm::dvec3 position {};
    std::optional<size_t> snap_voronoi_vertex_id {};
    std::optional<RadiusTransitionProjection> projection {};
  };
  /// Event-local cache: one canonical projection per (site, unordered source-edge pair); mesh XYZ filled on first insert.
  struct RadiusShiftedSiteCacheKey
  {
    size_t site_vertex_id = 0;
    size_t source_edge_lo = 0;
    size_t source_edge_hi = 0;

    bool operator==(const RadiusShiftedSiteCacheKey& other) const
    {
      return site_vertex_id == other.site_vertex_id && source_edge_lo == other.source_edge_lo
        && source_edge_hi == other.source_edge_hi;
    }
  };
  struct RadiusShiftedSiteCacheKeyHash
  {
    size_t operator()(const RadiusShiftedSiteCacheKey& key) const noexcept
    {
      size_t h = std::hash<size_t> {}(key.site_vertex_id);
      h ^= std::hash<size_t> {}(key.source_edge_lo) + 0x9e3779b9 + (h << 6) + (h >> 2);
      h ^= std::hash<size_t> {}(key.source_edge_hi) + 0x9e3779b9 + (h << 6) + (h >> 2);
      return h;
    }
  };
  struct RadiusShiftedSiteCacheEntry
  {
    RadiusTransitionSitePlacement placement {};
    std::optional<glm::dvec3> explicit_mesh_position {};
  };

  /// Deferred complementary-strip insert: recorded when the target meshlet/strip is known
  /// (often freshly started with verts 0/1 and no triangles yet), applied in @ref finalize.
  struct PendingRadiusComplementarySplit
  {
    size_t intersection_pair_index = 0;
    size_t edge_vertex_a = 0;
    size_t edge_vertex_b = 1;
    glm::dvec3 mesh_position {};
    double t = 0.0;
    size_t site_vertex_id = 0;
    size_t target_delaunay_edge = 0;
    /// True when insert came from radius-shift projection; false for non-shift VV/site insert.
    bool from_shift = true;
    /// True: split the closing edge of a just-finished strip (last triangle).
    /// False: split the seed edge of a newly started strip (first triangle, verts 0/1).
    bool split_last_triangle = false;
    /// True: @c intersection_pair_index refers to @c meshes / interior Voronoi-edge strip.
    /// False: boundary-interval meshlet in @c intersection_meshes.
    bool interior_strip = false;
    std::optional<size_t> snap_voronoi_vertex_id {};
    std::optional<size_t> insert_voronoi_edge_id {};

    /// True only for a logic-error duplicate (same pair, time, site, insert, and related fields).
    /// Delaunay-edge identity alone is not unique: many splits can share a target edge.
    bool isExactDuplicateOf(const PendingRadiusComplementarySplit& other) const
    {
      return intersection_pair_index == other.intersection_pair_index && interior_strip == other.interior_strip
        && t == other.t && site_vertex_id == other.site_vertex_id
        && target_delaunay_edge == other.target_delaunay_edge && insert_voronoi_edge_id == other.insert_voronoi_edge_id
        && split_last_triangle == other.split_last_triangle && from_shift == other.from_shift
        && snap_voronoi_vertex_id == other.snap_voronoi_vertex_id && mesh_position.x == other.mesh_position.x
        && mesh_position.y == other.mesh_position.y && mesh_position.z == other.mesh_position.z;
    }
  };

  void clearRadiusShiftedSiteCache();
  static RadiusShiftedSiteCacheKey makeRadiusShiftedSiteCacheKey(
    size_t site_vertex_id, const RadiusBoundaryTransitionShiftContext& ctx);
  /// Fill-once cache lookup; canonical strip = min(source edges). Returns nullptr when shift is inapplicable.
  RadiusShiftedSiteCacheEntry* getOrFillRadiusShiftedSiteCache(double t, size_t site_vertex_id,
    const RadiusBoundaryTransitionShiftContext* boundary_transition_shift);
  void noteRadiusShiftedSiteMeshPosition(
    RadiusShiftedSiteCacheEntry& entry, const glm::dvec3& mesh_position, const char* consumer_context);

  /**
   * @brief True when @p start/@p end are the two remapped shift-anchor crossings on @c ctx.target_delaunay_edge.
   *
   * @details A Delaunay boundary edge with crossings @c refs ordered along the boundary is partitioned into
   * intervals @c [null,refs.front()], @c [refs[k],refs[k+1]] for each consecutive pair, and
   * @c [refs.back(),null]. The consecutive @c [refs[k],refs[k+1]] pieces are the mid-intervals (as opposed to the
   * one-null end intervals). For radius boundary-transition shift, the complementary mid-interval is the one whose
   * endpoints are exactly the two projected anchor crossings on the surviving / target Delaunay edge.
   */
  bool midIntervalMatchesRadiusShiftAnchors(KineticDelaunay::CrossingData::EdgeIntersectionRef start,
    KineticDelaunay::CrossingData::EdgeIntersectionRef end, size_t site_vertex_id,
    const RadiusBoundaryTransitionShiftContext* boundary_transition_shift, double t);

  /// Queue a deferred split of edge (@p edge_a,@p edge_b) on @p intersection_pair_index (typically 0,1).
  /// @p split_last_triangle: finished/closing strip (last face) vs newly started strip (first face).
  /// @p interior_strip: target @c meshes (Voronoi-edge strip) instead of @c intersection_meshes.
  void queuePendingRadiusComplementarySplit(size_t intersection_pair_index, size_t edge_a, size_t edge_b,
    const glm::dvec3& mesh_position, double t, size_t site_vertex_id, size_t target_delaunay_edge,
    bool from_shift = true, std::optional<size_t> snap_voronoi_vertex_id = std::nullopt,
    bool split_last_triangle = false, bool interior_strip = false,
    std::optional<size_t> insert_voronoi_edge_id = std::nullopt);

  /**
   * @brief Queue a complementary insert on a mid-interval strip that was just started (2→1 radius shift).
   *
   * @details After reseeding @c target_delaunay_edge, @c startNewMeshFromIntersections creates a fresh boundary
   * meshlet for each mid-interval @c [crossing,crossing]. When that interval is the shift complementary mid
   * (see @ref midIntervalMatchesRadiusShiftAnchors), this queues a split of its seed edge (mesh verts 0–1) that
   * inserts the shifted shared-site corner. The strip usually has no triangles yet, so the split is deferred until
   * @ref flushPendingRadiusComplementarySplits / finalize.
   *
   * @param intersection_pair_index New @c intersection_meshes pair from @c startNewMeshFromIntersections.
   * @param start,@p end Mid-interval endpoints (ordered along the boundary).
   * @param site_vertex_id Shared Delaunay site of the two source boundary edges (shift junction).
   */
  void maybeQueueRadiusComplementarySplitForStartedMid(size_t intersection_pair_index,
    KineticDelaunay::CrossingData::EdgeIntersectionRef start, KineticDelaunay::CrossingData::EdgeIntersectionRef end,
    double t, size_t site_vertex_id, const RadiusBoundaryTransitionShiftContext* boundary_transition_shift);

  /**
   * @brief Split (or queue) the complementary mid-interval strip that already existed before a 1→2 radius shift.
   *
   * @details On 1→2, @c beforeEvent finished the mid-interval meshlet on the old single boundary edge
   * (@c target_delaunay_edge). After the event that edge leaves the alpha boundary, so CrossingData
   * @c prev/@c next links are cleared before this runs — do not rely on link lookup. Instead pass the
   * pair indices captured while finishing that edge (@p finished_one_edge_pair_indices); this picks the
   * complementary mid among them (both endpoints on @c target_d matching the remapped shift anchors),
   * then applies @ref applyPendingRadiusComplementarySplit immediately (@c split_last_triangle), or
   * queues it for finalize if the immediate apply fails.
   *
   * Contrasts with @ref maybeQueueRadiusComplementarySplitForStartedMid (2→1 / newly started seed edge).
   *
   * @param site_vertex_id Shared Delaunay site of the two source boundary edges (shift junction).
   * @param finished_one_edge_pair_indices Meshlets finished on the pre-flip single boundary edge.
   */
  void maybeQueueRadiusComplementarySplitForExistingMid(double t, size_t site_vertex_id,
    const RadiusBoundaryTransitionShiftContext* boundary_transition_shift,
    const std::vector<size_t>& finished_one_edge_pair_indices);

  /// Interior Voronoi-edge strip whose current seed/top edge contains @p crossing as an interior point.
  struct InteriorStripSplitTarget
  {
    size_t segment_mesh_pair_index = static_cast<size_t>(-1);
    size_t edge_vertex_a = 0;
    size_t edge_vertex_b = 1;
  };
  std::optional<InteriorStripSplitTarget> findInteriorStripSplitTarget(size_t voronoi_edge_id,
    KineticDelaunay::CrossingData::EdgeIntersectionRef crossing) const;
  /// Non-shift: insert a crossing that left or joined the alpha boundary onto the interior strip of its Voronoi edge.
  /// @p apply_immediately: II→IO after @c finishMesh (top triangle exists). Otherwise queue until @ref finalize
  /// (IO→II seed edge; the first triangle exists only after later events).
  void maybeQueueOrApplyRadiusNonShiftInteriorCrossingSplit(
    KineticDelaunay::CrossingData::EdgeIntersectionRef crossing, double t, bool apply_immediately);
  /// Non-shift IO↔II inserts, grouped by Voronoi-edge strip. Neighboring crossings that share one strip edge
  /// (triangle between them flipped) are applied in order along that edge so the second split hits the correct
  /// remaining sub-triangle. @p apply_immediately: II→IO after @c finishMesh; otherwise queue the seed-edge case
  /// until @ref finalize.
  void applyOrQueueRadiusNonShiftInteriorCrossingSplits(
    const std::vector<KineticDelaunay::CrossingData::EdgeIntersectionRef>& crossings, double t, bool apply_immediately);

  /// Apply buffered complementary splits (before other intersection-mesh postprocess in @ref finalize).
  void flushPendingRadiusComplementarySplits();
  /// Apply any pending complementary split targeting @p intersection_pair_index (e.g. just before it is retired).
  void flushPendingRadiusComplementarySplitsForPair(size_t intersection_pair_index);
  /// @p allow_completed: apply even if the pair was just retired (last-triangle insert on the finished strip).
  bool applyPendingRadiusComplementarySplit(const PendingRadiusComplementarySplit& pending, bool allow_completed);

  /// Cleared at the start of each radius beforeEvent / afterEvent phase.
  std::unordered_map<RadiusShiftedSiteCacheKey, RadiusShiftedSiteCacheEntry, RadiusShiftedSiteCacheKeyHash>
    radius_shifted_site_cache_;
  std::vector<PendingRadiusComplementarySplit> pending_radius_complementary_splits_;

  glm::dvec3 radiusTransitionProjectionPosition(
    const RadiusTransitionProjection& projection, double t, bool mesh_space) const;
  std::optional<RadiusTransitionSitePlacement> radiusTransitionInterpolatedSitePosition(double t, size_t site_vertex_id,
    size_t strip_delaunay_edge_id, const RadiusBoundaryTransitionShiftContext* boundary_transition_shift) const;

  static void appendIntersectionInterpolationDebugToMetadata(
    MetadataBuilder& builder, const IntersectionInterpolationDebug& debug);

  /**
   * @brief Returns Delaunay-edge intersections in component-boundary traversal order.
   *
   * @details `CrossingData::delaunay_edge_intersections[e]` is stored in canonical even-half-edge order.
   * This helper reorders (possibly reverses) that list so callers can build `[null,first]`, `[k,k+1]`,
   * `[last,null]` intervals consistently in boundary-walk direction.
   */
  std::vector<KineticDelaunay::CrossingData::EdgeIntersectionRef> getBoundaryIntersectionsInBoundaryOrder(
    size_t delaunay_edge_id) const;
  /// Clears @c prev_segment_mesh_pair_index / @c next_segment_mesh_pair_index on all crossings along @p delaunay_edge_id.
  void clearIntersectionMeshPairLinksOnDelaunayEdge(size_t delaunay_edge_id);
  /**
   * When @p even_half_edge_id is no longer on the alpha-shape component boundary, drop crossing mesh-pair links and
   * boundary-mesh strip anchors for that Delaunay edge only (does not modify open intersection-strip vertex state).
   */
  void invalidateStaleAlphaBoundaryMeshLinksOnEvenHalfEdge(size_t even_half_edge_id);
  /**
   * When a Delaunay edge of @p triangle_half_edges (post-event triangle) is no longer on the alpha boundary, clear
   * @c prev_segment_mesh_pair_index / @c next_segment_mesh_pair_index on all crossings along that edge.
   */
  void invalidateStaleAlphaBoundaryMeshLinksOnTriangleEdgesLeftBoundary(
    const std::array<size_t, 3>& triangle_half_edges);
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
   * @return true if the pair index was stored in @c prev_segment_mesh_pair_index, false if in @c
   * next_segment_mesh_pair_index.
   */
  bool writeOneNullIntersectionPairLinkByNullVertex(size_t intersection_pair_index, size_t null_vertex_id,
    KineticDelaunay::CrossingData::EdgeIntersectionRef ref, bool interval_is_ref_to_null,
    double t = std::numeric_limits<double>::quiet_NaN());
  /**
   * Assign @p new_pair_index to @c prev_segment_mesh_pair_index or @c next_segment_mesh_pair_index on @p ref.
   * When @ref diagnostics is enabled and @p new_pair_index equals @ref kDiagnosticsMonitoredMeshPairId, emits
   * @c KINDS_INFO with @p context, crossing coordinates, and optional @p t.
   */
  void assignIntersectionMeshPairLink(KineticDelaunay::CrossingData::EdgeIntersectionRef ref, bool is_prev_link,
    size_t new_pair_index, const char* context, double t = std::numeric_limits<double>::quiet_NaN());
  /**
   * @brief Verifies that adjacent crossings on @p delaunay_edge_id share mesh-pair links:
   * `crossing[i].next_segment_mesh_pair_index == crossing[i+1].prev_segment_mesh_pair_index`.
   *
   * @throws std::runtime_error with every crossing's prev/next indices when the invariant fails.
   */
  void validateDelaunayEdgeIntersectionMeshPairLinks(
    size_t delaunay_edge_id, double event_time, const char* context) const;
  /// Mesh-pair indices as @c "(prev0, next0, prev1, next1, ...)" (@c -1 for unset).
  static std::string formatCrossingMeshPairLinkSequence(
    const std::vector<std::pair<size_t, size_t>>& prev_next_pairs);
  /// Same as @ref formatCrossingMeshPairLinkSequence from current CrossingData on @p delaunay_edge_id.
  std::string formatDelaunayEdgeCrossingMeshPairLinkSequence(size_t delaunay_edge_id) const;
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
    std::optional<KineticDelaunay::CrossingData::EdgeIntersectionRef> end_intersection,
    double t = std::numeric_limits<double>::quiet_NaN());
  size_t determineVoronoiCellForBoundaryIntersectionInterval(size_t delaunay_edge_id,
    std::optional<KineticDelaunay::CrossingData::EdgeIntersectionRef> start_intersection,
    std::optional<KineticDelaunay::CrossingData::EdgeIntersectionRef> end_intersection) const;

  void completeBoundaryMeshSection(size_t he_id, size_t new_left, size_t new_right, double t);

  /**
   * Add a triangle to the boundary mesh. Automatically takes the raw UVs and adjusts them to avoid seams.
   * @param u first vertex index
   * @param v second vertex index
   * @param w third vertex index
   */
  size_t addBoundaryTriangle(
    size_t u, size_t v, size_t w, const std::string& metadata = "{}", int material_id = 0);
  /// Per-corner bark UV for boundary-interval meshlets (same pipeline as @ref addBoundaryTriangle).
  size_t addBoundaryIntervalTriangle(
    VoronoiMesh& mesh, size_t u, size_t v, size_t w, const std::string& metadata = "{}", int material_id = 1);
  std::vector<glm::dvec2>& boundaryIntervalRawUvs(VoronoiMesh& mesh);
  const std::vector<glm::dvec2>& boundaryIntervalRawUvs(const VoronoiMesh& mesh) const;
  std::optional<glm::dvec2> boundaryIntervalRawUvAtVertex(const VoronoiMesh& mesh, size_t vertex_index) const;
  void setBoundaryIntervalRawUv(VoronoiMesh& mesh, size_t vertex_index, const glm::dvec2& raw_uv);
  void refreshBoundaryIntervalTrianglesIncidentToVertex(VoronoiMesh& mesh, size_t vertex_index);

  /// Unscaled boundary UV in Delaunay space: normalized polar angle and kinetic height.
  static glm::dvec2 boundaryRawUv(const glm::dvec2& delaunay_xy, const glm::dvec2& centroid, double t);
  /// Interior UV from Delaunay space: radial distance normalized by the component boundary, converted back to Cartesian
  /// disk coordinates, plus kinetic height.
  glm::dvec3 interiorMeshUv(const std::vector<BoundaryPoint>& boundary_polygon, const glm::dvec2& centroid,
    const glm::dvec2& delaunay_xy, double t) const;
  size_t addBoundaryVertex(
    glm::dvec3 vertex, glm::dvec2 centroid, size_t strand_id, double t, bool includes_virtual_shift);

  size_t addMeshletTriangle(VoronoiMesh& mesh, size_t u, size_t v, size_t w, const std::string& metadata = "{}",
    int material_id = RegularMeshletMaterialId);

  /**
   * Same as @ref addMeshletTriangle after orienting `(u,v,w)` from combinatorics only: @p inside_boundary_he_id is the
   * *inside* boundary Delaunay half-edge on this strip; its twin is the *outside* boundary half-edge (`index ^ 1`).
   * Apex order relative to strip `(u,v)` is chosen from that outside half-edge index alone (no vertex cross products).
   */
  size_t addBoundaryIntervalTriangleOriented(VoronoiMesh& mesh, size_t u, size_t v, size_t w, int inside_boundary_he_id,
    double t, const std::string& metadata = "{}", std::optional<size_t> boundary_meshlet_id = std::nullopt);

  /// @p strand_id Delaunay site id for transform frame lookup; may be a Voronoi vertex id when it matches
  /// @p meshlet_voronoi_vertex_for_alpha_check (resolved to an adjacent site automatically).
  size_t addMeshletVertex(VoronoiMesh& mesh, const std::vector<BoundaryPoint>& boundary_polygon,
    const glm::dvec2& centroid, glm::dvec3 vertex, size_t strand_id, double t, bool includes_virtual_shift,
    std::optional<size_t> meshlet_voronoi_vertex_for_alpha_check = std::nullopt, const std::string& metadata = "{}",
    const std::optional<glm::dvec3>& debug_color = std::nullopt);
  size_t addMeshletVertex(VoronoiMesh& mesh, const std::vector<BoundaryPoint>& boundary_polygon,
    const glm::dvec2& centroid, glm::dvec3 vertex, size_t strand_id, double t, bool includes_virtual_shift,
    std::optional<size_t> meshlet_voronoi_vertex_for_alpha_check, const std::string& metadata,
    const std::optional<glm::dvec3>& debug_color, const MeshletVertexRuntimeInfo& runtime_info);

  /// Effective strip corner for triangles: latest flexible on @p left_side, else fixed @c mesh_start / @c mesh_end.
  size_t intersectionStripEffectiveVertexIndex(const MeshingData& seg, bool left_side) const;

  /// When @ref intersection_strip_flexible_vertices_enabled: appends a placeholder (xy = @p centroid, z=@p t, metadata
  /// `pos` left/right) and a wedge triangle `(eff_left, eff_right, flex)` per @ref
  /// intersectionStripEffectiveVertexIndex (snapshot before the new flex is appended). Otherwise no-op.
  void addFlexibleVertexToIntersectionMesh(VoronoiMesh& mesh, MeshingData& seg, bool flexible_on_left_side,
    const std::vector<BoundaryPoint>& boundary_polygon, const glm::dvec2& centroid, size_t strand_id, double t,
    const std::string& metadata = "{}");

  /**
   * @brief Wedge-extend a boundary-interval intersection mesh at a shared Voronoi–Delaunay crossing.
   *
   * @details Uses @p neighbor_pair_idx from the crossing's @c prev_segment_mesh_pair_index (update strip start) or
   * @c next_segment_mesh_pair_index (update strip end). No-op when the pair index is invalid or @p skip_pair_idx
   * matches. When @p append_flexible_placeholder is false, endpoint assignment still runs but no flex placeholder wedge
   * is added.
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
    bool keep_strip_alive, bool append_flexible_placeholder = true, const std::string& flexible_base_metadata = "{}");

  /**
   * For a merge closure vertex (`pos` uniform): both flex chains lerp from their fixed corners to the same
   * @p closure_vertex_index, then lists are cleared. No endpoint reassignment and no new flex.
   */
  void applyIntersectionStripUniformClosureVertex(VoronoiMesh& mesh, MeshingData& seg, size_t closure_vertex_index);
  void resolveRemainingFlexibleVertices(VoronoiMesh& mesh, MeshingData& seg, const char* context);
  bool interpolateFlexibleVerticesAlongEdge(VoronoiMesh& mesh, std::vector<int>& flex, size_t anchor_old_vertex,
    size_t anchor_new_vertex);
  void snapFlexibleVerticesToAnchor(VoronoiMesh& mesh, const std::vector<int>& flex, size_t anchor_vertex);
  void resolveAllIntersectionFlexibleVertices(const char* context);
  /// Collapse flexible vertices that only bind two triangles (see @ref VoronoiMesh::collapseDegreeTwoFlexibleVertices).
  void collapseDegreeTwoFlexibleVerticesOnIntersectionMeshes();

  /// Glue-edge unmatched-vertex alignment runs during @ref extractSegmentMeshlets on contributing copies.
  /// Reserved for any remaining segment-meshlet post-steps after combine/merge.
  void postProcessFlexibleVerticesOnSegmentMeshlet(VoronoiMesh& segment_mesh) const;

  /// If the containing Delaunay triangle for @p voronoi_vertex_id is not inside the alpha-shape, log a warning with @p
  /// position and strand/branch ids at @p t.
  void warnIfVoronoiVertexOutsideAlphaShape(const char* context, size_t voronoi_vertex_id, const glm::dvec3& position,
    size_t strand_id, double t) const;

  /// Compact strand / runtime-branch / input-branch suffix for warning logs.
  std::string formatStrandBranchLogInfo(size_t strand_id, double t) const;

  /// Validates both circumcenters (Delaunay face ids) incident to a Voronoi edge.
  void requireLiveRegisteredVoronoiEdgeEndpoints(size_t voronoi_edge_id, const char* context) const;

  void addDelaunayTriangulationToBoundaryMesh(
    double t, size_t input_branch_id, bool invert_orientation, double offset);

  std::vector<BoundaryPoint> traceConvexHull(double t) const;

  void advanceBoundaryMesh(double t, const std::vector<BoundaryPoint>& boundary_points, const glm::dvec2& centroid);

  void updateBoundary(double t, std::vector<bool>& visited, size_t component_index);

  void updateBoundaries(double t, const std::vector<size_t>& component_indices);

  bool isComponentLive(size_t component_index) const;

  std::vector<size_t> collectLiveComponentIndices() const;

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
    bool includes_virtual_shift, std::optional<size_t> voronoi_vertex_for_alpha_check = std::nullopt,
    const std::string& metadata = "{}");
  int closingMeshAppendVertex(VoronoiMesh& mesh, const std::vector<BoundaryPoint>& boundary_polygon,
    const glm::dvec2& centroid, size_t strand_id, double t, const glm::dvec3& position, bool includes_virtual_shift,
    std::optional<size_t> voronoi_vertex_for_alpha_check, const std::string& metadata,
    const MeshletVertexRuntimeInfo& runtime_info);

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
   * (erase/insert), since stored iterators into `edge_intersections` become invalid while `std::optional` may still
   * hold a value.
   */
  void refreshMeshingDataCrossingRefs(MeshingData& seg, size_t voronoi_edge_id) const;

  /// Refresh every segment strip whose dual Voronoi edge is incident to this Voronoi vertex (after a crossing event).
  void refreshStripCrossingRefsIncidentToVoronoiVertex(size_t voronoi_vertex_id);

  /// Refresh all segment strips (after flips / broad `CrossingData` updates).
  void refreshCrossingRefsForAllStrips();

  void growGraphSlotArrays();
  void clearDeadHalfEdgeState();
  void initializeNewHalfEdgesAfterGraphUpdate(double t, size_t first_new_he_slot);
  void refreshCrossingRefsForAllIntersectionStrips();

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
    const std::function<int(
      const glm::dvec3&, std::optional<size_t>, const std::string&, const MeshletVertexRuntimeInfo&)>& track_vertex,
    bool reverse = false);

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
    size_t incident_he,
    const std::function<int(
      const glm::dvec3&, std::optional<size_t>, const std::string&, const MeshletVertexRuntimeInfo&)>& track_vertex);

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
  void closingMeshValidateOrderedSegmentGeometry(
    double t, const VoronoiMesh& mesh, const std::vector<MeshingData*>& ordered_segments) const;

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
    const std::unordered_map<KineticDelaunay::CrossingData::EdgeIntersectionRef, size_t,
      ClosingMeshCrossingIteratorHash, ClosingMeshCrossingIteratorEq>& start_crossing_to_segment);

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
   * @brief Triangulates one simple polygon ring into @p mesh using ear clipping.
   * Rings that revisit the same Delaunay-plane XY position are split into sub-rings first.
   * Uses @ref VoronoiMesh::triangulationPlaneXY (Delaunay space when meshes are object-space transformed).
   * @param polygon Vertex index ring into @p mesh.
   * @param orient_upwards If true, emits CCW XY triangles; otherwise emits CW triangles.
   * @param occurrence_time Optional kinetic @ref EventTime for failure/split debug dumps (path + TXT header).
   * @param runtime_branch_id Optional runtime branch for failure/split debug dump folder naming.
   */
  void triangulateSimplePolygon(VoronoiMesh& mesh, const std::vector<size_t>& polygon,
    const std::string& metadata = "{}", int material_id = RegularMeshletMaterialId, bool orient_upwards = true,
    std::optional<EventTime> occurrence_time = std::nullopt, std::optional<size_t> runtime_branch_id = std::nullopt);

  /**
   * @brief Fan-triangulates a convex polygon ring (no geometric tests).
   * Used by radius traced cell / triangle-cap meshlets. Does not use plane XY or ear clipping.
   * @param polygon Vertex index ring into @p mesh (assumed convex, correctly ordered).
   * @param orient_upwards If true, emits (0,i,i+1); otherwise emits (0,i+1,i).
   */
  void fanTriangulateConvexPolygon(VoronoiMesh& mesh, const std::vector<size_t>& polygon,
    const std::string& metadata = "{}", int material_id = RegularMeshletMaterialId, bool orient_upwards = true);

  /**
   * @brief Triangulates each traced closing polygon ring into the cap mesh (in-place).
   * @param mesh Cap mesh (triangles appended).
   * @param polygons Vertex index rings.
   */
  void closingMeshTriangulatePolygons(VoronoiMesh& mesh, const std::vector<std::vector<size_t>>& polygons, double t,
    size_t strand_id);

  size_t createClosingMesh(size_t strand_id, double t, const std::vector<BoundaryPoint>& boundary_polygon,
    const glm::dvec2& centroid, std::vector<BoundaryIntersectionInterval>* traced_boundary_intervals = nullptr);

  /// Top closing cap for one live strand at @p t; wires the new mesh pair to the strand's latest segment.
  void createClosingCapForStrand(size_t strand_id, double t);

  /// Top closing caps for every live strand in an input branch whose last section is @p t.
  void createClosingCapsForInputBranchFinishingAtSection(double t, size_t input_branch_id);

  /// Top closing caps for all input branches finishing at kinetic section @p t.
  void createClosingCapsForInputBranchesFinishingAtSection(double t);

  /// Finishes Voronoi strip meshlets on all half-edges incident to @p strand_id at @p t.
  void finishIncidentStripMeshesForStrandAtSection(size_t strand_id, double t);

  void accumulateSegmentProperties();

  /// Registers a newly created meshlet and stores export suffix metadata.
  /// @param creation_kinetic_time If finite, stored on the mesh for export/debug (see @ref
  /// VoronoiMesh::setCreationKineticTime).
  size_t registerMeshletWithSuffix(
    VoronoiMesh&& mesh, std::string suffix, double creation_kinetic_time = std::numeric_limits<double>::quiet_NaN());

 public:
  double getUvHeightFactor() const { return uv_height_factor; }
  double getUvCircumFactor() const { return uv_circum_factor; }
  double getTextureDiameter() const { return texture_diameter; }

  /// One-line Voronoi meshlet diagnostics: dual edge id, pair slot, verts, tris, strip counts (@p extra_note optional).
  void meshletDiagnosticLogLine(const char* tag, size_t half_edge_id, double t, const char* extra_note = "") const;

  /// When @ref diagnostics is enabled, log @ref KineticDelaunay::kDiagnosticsMonitoredFaceId inside state at @p t.
  void logDiagnosticsMonitoredFaceInsideState(double t, const char* event_context) const;

  /// Debug target: Delaunay edge whose crossing @c prev/@c next mesh-pair links are traced.
  /// Set to @ref KineticDelaunay::kDiagnosticsMonitorDisabledId to disable.
  static constexpr size_t kDiagnosticsMonitoredDelaunayEdgeId = KineticDelaunay::kDiagnosticsMonitorDisabledId;
  /// Suspected missed crossing time for @ref kDiagnosticsMonitoredDelaunayEdgeId.
  static constexpr double kDiagnosticsMonitoredCrossingTime = 10.0;
  /// Debug target: boundary-interval mesh pair id suspected of incorrect wiring.
  /// Set to @ref KineticDelaunay::kDiagnosticsMonitorDisabledId to disable (must not collide with cleared-link sentinel).
  static constexpr size_t kDiagnosticsMonitoredMeshPairId = KineticDelaunay::kDiagnosticsMonitorDisabledId;
  /// Suspected incorrect flip event (log monitored edge state in [floor(t), floor(t)+1)).
  /// Keep in sync with @ref KineticDelaunay::kDiagnosticsMonitoredFlipTime.
  static constexpr double kDiagnosticsMonitoredFlipTime = 20.0;
  static constexpr double kDiagnosticsMonitoredTimeEpsilon = 0.05;

  /// Full snapshot of crossings on a Delaunay edge (default: @ref kDiagnosticsMonitoredDelaunayEdgeId)
  /// and mesh-pair metadata they reference. No-op when @p delaunay_edge_id is disabled (@ref KineticDelaunay::kDiagnosticsMonitorDisabledId).
  void logDiagnosticsMonitoredDelaunayEdgeState(double t, const char* event_context,
    size_t delaunay_edge_id = kDiagnosticsMonitoredDelaunayEdgeId) const;

  /// Log monitored-edge snapshot when @p delaunay_edge_id / @p mesh_pair_index / @p t matches debug targets.
  void maybeLogDiagnosticsMonitoredDelaunayEdgeTrigger(double t, const char* event_context,
    std::optional<size_t> delaunay_edge_id = std::nullopt, std::optional<size_t> mesh_pair_index = std::nullopt) const;

  /// After @ref startNewMesh strip build: warn if topology/metadata suggests a non-empty mesh but vertices are missing.
  void meshletDiagnosticWarnIfUnexpectedEmptyAfterStartNewMesh(size_t half_edge_even, double t,
    bool initial_left_inside, const VoronoiMesh& mesh, const std::list<MeshingData>& strips) const;

  /// Initial-cap / strand-segment wiring diagnostics (@p tag is breakpoint anchor, e.g. `init_cap_begin`).
  void strandInitDiagnosticLogLine(const char* tag, size_t strand_id, double t, const char* extra_note = "") const;

  /// Post-init snapshot for strand 0 (and general wiring) when @ref diagnostics is enabled.
  void logStrandInitDiagnosticsSummary(double t) const;

  SegmentBuilder(KineticDelaunay& kin_del, std::vector<std::pair<size_t, double>> subdivisions,
    bool create_transformed_mesh, bool visual_debug = false,
    std::function<void(size_t, std::function<void(size_t)>)> parallel_for = {});

  SegmentBuilder(KineticDelaunay& kin_del, bool create_transformed_mesh, bool visual_debug = false,
    std::function<void(size_t, std::function<void(size_t)>)> parallel_for = {});

  ~SegmentBuilder() override;

  void init() override;

  void onGraphRetriangulated(double t, size_t prev_face_slots, size_t prev_he_slots) override;
  void onGraphCutApplied(double t, size_t prev_face_slots, size_t prev_he_slots) override;
  void onBeforeComponentGraphSplit(double t) override;

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

  /// Strand whose profile frame maps @p segment_id mesh vertices to world space.
  size_t strandIdForSegment(size_t segment_id) const;
  /// Strand used to transform a raw (unmerged) meshlet at @p meshlet_index.
  size_t strandIdForRawMeshlet(size_t meshlet_index) const;

  std::vector<glm::dvec3> computeClampedVoronoiVertices(
    size_t strand_id, double t, const std::vector<BoundaryPoint>& boundary_polygon, const glm::dvec2& centroid);

  void splitComponent(size_t component_id, const std::vector<std::vector<size_t>>& new_components, double t);
};
} // namespace kinDS
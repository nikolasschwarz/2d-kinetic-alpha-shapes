#pragma once
#include "HalfEdgeDelaunayGraph.hpp"
#include "HalfEdgeDelaunayGraphToSVG.hpp"
#include "KineticAlgorithm.hpp"
#include "Logger.hpp"
#include "Polynomial.hpp"
#include "ProgressBar.hpp"
#include "StrandTree.hpp"
#include <filesystem>
#include <format>
#include <glm/gtx/exterior_product.hpp>
#include <glm/gtx/string_cast.hpp>
#include <map>
#include <memory>
#include <unordered_map>
#include <unordered_set>
#include <optional>
#include <queue>
#include <string>
#include <utility>
#include <vector>

namespace kinDS
{
struct BoundaryPoint
{
  size_t vertex_id;
  size_t he_id;
  glm::dvec2 p;
};

/// Data recorded when a component split awaits graph retriangulation.
struct PendingBranchSplit
{
  size_t parent_component_id = static_cast<size_t>(-1);
  /// Reference strand from the parent component used for the frozen frame.
  size_t reference_vertex = static_cast<size_t>(-1);
  /// Full pre-split parent-component strand list.
  std::vector<size_t> frozen_parent_strands;
  /// Component ids from the split: retained parent first, then split-off components.
  std::vector<size_t> split_component_ids;
  glm::dvec2 old_branch_centroid { 0.0, 0.0 };
  glm::dvec2 new_branch_centroid { 0.0, 0.0 };
  glm::dvec2 separation_direction { 0.0, 0.0 };
  double separation_t0 = 0.0;
  double separation_te = 0.0;
  size_t separation_iteration = 0;
  bool separation_trajectory_active = false;
};

/// Registry of pending branch splits and per-strand lookup into the frozen parent frame.
struct PendingBranchSplitState
{
  static constexpr size_t no_parent_component = static_cast<size_t>(-1);

  void clear();
  void resetStrandLookup(size_t strand_count);
  PendingBranchSplit& getOrCreate(size_t parent_component_id);
  const PendingBranchSplit* findByParent(size_t parent_component_id) const;
  void registerStrandsForSplit(
    size_t parent_component_id, const std::vector<std::vector<size_t>>& new_components);
  const std::vector<size_t>* frozenStrandsForStrand(size_t strand_id) const;

  std::unordered_map<size_t, PendingBranchSplit> by_parent_;
  std::vector<size_t> strand_parent_component_;
};

static glm::dvec2 polygonCentroid(const std::vector<BoundaryPoint>& polygon)
{
  double A = 0.0;
  glm::dvec2 C { 0.0, 0.0 };

  const size_t n = polygon.size();
  for (size_t i = 0; i < n; ++i)
  {
    const glm::dvec2& p = polygon[i].p;
    const glm::dvec2& q = polygon[(i + 1) % n].p;

    double cross = glm::cross(p, q);
    A += cross;
    C += (p + q) * cross;
  }

  A *= 0.5;

  if (std::abs(A) < 1e-12)
    return C; // degenerate polygon

  return C / (6.0 * A);
}

static std::vector<size_t> buildComponentMap(const std::vector<std::vector<size_t>>& components, size_t vertex_count)
{
  std::vector<size_t> component_map(vertex_count);

  for (size_t i = 0; i < components.size(); i++)
  {
    for (const auto v : components[i])
    {
      component_map[v] = i;
    }
  }

  return component_map;
}

/**
 * \brief Class for computing the Delaunay triangulation of a set of cubic Hermite splines.
 *
 * This follows Guibas, L.J., Mitchell, J.S.B., Roos, T. (1992).
 * Voronoi diagrams of moving points in the plane. In:
 * Schmidt, G., Berghammer, R. (eds) Graph-Theoretic Concepts in Computer Science. WG 1991.
 * Lecture Notes in Computer Science, vol 570. Springer, Berlin, Heidelberg.
 * https://doi.org/10.1007/3-540-55121-2_11
 */
class KineticDelaunay
{
 public:
  using Event = KineticAlgorithm::Event;
  using EventCallback = KineticAlgorithm::EventCallback;
  using EventManager = KineticAlgorithm::EventManager;

  enum class ComponentSplitPolicy
  {
    Retriangulate,
    InPlaceCut,
  };

  class FlipEvent;
  class RadiusEvent;
  class CrossingEvent;

  class CallbackManager
  {
   public:
    virtual ~CallbackManager() = default;

    virtual void init() { }
    virtual void finalize(double t) { }
    virtual void onGraphRetriangulated(double t, size_t prev_face_slots, size_t prev_he_slots) { (void)t; (void)prev_face_slots; (void)prev_he_slots; }
    virtual void onGraphCutApplied(double t, size_t prev_face_slots, size_t prev_he_slots) { (void)t; (void)prev_face_slots; (void)prev_he_slots; }
  };

  class FlipEvent;
  class RadiusEvent;
  class CrossingEvent;
  class SectionEvent;

  class FlipEventManager;
  class RadiusEventManager;
  class CrossingEventManager;
  class SectionEventManager;
  class SubdivisionEvent;
  class SubdivisionEventManager;
  class SeparationEvent;
  class SeparationEventManager;

  /// Forward declaration; full definition in `KineticDelaunayCrossingEvent.hpp`.
  struct CrossingData;

  struct ComponentData
  {
    std::vector<std::vector<size_t>> components;
    std::vector<size_t> component_map;
    // [component_index][boundary_no][point_no] - the first boundary is the outer one, any additional ones are holes in
    // the polygon
    std::vector<std::vector<std::vector<BoundaryPoint>>> component_boundaries;
    std::vector<glm::dvec2> component_centroids;
    std::vector<double> component_last_updated;
  };

  ComponentData component_data;

  std::pair<glm::dvec2, glm::dvec2> computeAngularBisector(size_t he_id, double t) const;

  std::pair<double, double> delaunayVoronoiEdgeIntersection(
    size_t delaunay_edge_id, size_t voronoi_edge_id, double t) const;

  /**
   * \brief Compute the half-edges crossed by the Voronoi edge between the given start point and destination, starting
   * from the given face.
   *
   * Currently assumes that the start face is finite and that the destination lies outside of only one edge of the start
   * triangle.
   */
  std::pair<std::vector<size_t>, std::vector<double>> computeCrossedHalfEdges(
    size_t start_face_id, const glm::dvec2& destination, const glm::dvec2& start_point, double t) const;

 private:
  friend class SectionEventManager;
  friend class SubdivisionEvent;
  friend class SeparationEvent;

  StrandTree branch_trajs;
  HalfEdgeDelaunayGraph graph;
  std::unique_ptr<KineticAlgorithm> kinetic_algorithm_;
  // Reused managers: one per event type.
  // Kept as pointers to avoid forcing complete manager types in this header.
  std::unique_ptr<FlipEventManager> flip_event_manager_;
  std::unique_ptr<RadiusEventManager> radius_event_manager_;
  std::unique_ptr<CrossingEventManager> crossing_event_manager_;
  std::unique_ptr<SectionEventManager> section_event_manager_;
  std::unique_ptr<SubdivisionEventManager> subdivision_event_manager_;
  std::unique_ptr<SeparationEventManager> separation_event_manager_;
  /// Strand mesh subdivision schedule; consumed on the first @ref enqueueScheduledSubdivisionEvents in @ref compute.
  std::vector<std::pair<size_t, double>> subdivision_schedule_;
  CallbackManager* callback_manager_ = nullptr;
  size_t sections_advanced = 0; // Counter for the number of sections advanced
  double cutoff; // Cutoff radius for boundary events
  std::vector<bool> face_inside; // Tracks whether faces are inside or outside the boundary

  std::vector<std::vector<size_t>> branches; // track which vertices/splines belong to which branch
  std::vector<glm::dvec2> dummy_boundary;
  bool add_dummy_boundary;
  size_t prev_component_count = 1;
  PendingBranchSplitState pending_branch_splits_;
  /// Per-strand runtime branch id (see @ref getRuntimeBranchIdForStrand). Differs from @ref ComponentData::component_map
  /// (inside-face connectivity) and from @ref StrandTree::getBranchIndex (input branch). Updated only when Delaunay
  /// graph connectivity is severed (@ref onGraphCutApplied / @ref onGraphRetriangulated), using live-face adjacency.
  std::vector<size_t> runtime_branch_map_;
  /// Monotonic runtime-branch id allocator; ids are never reused.
  size_t next_runtime_branch_id_ = 0;
  /// Liveness per runtime branch id (index = id).
  std::vector<bool> runtime_branch_alive_;
  ComponentSplitPolicy component_split_policy_ = ComponentSplitPolicy::InPlaceCut;
  double separation_offset_scale_ = 100.0;
  std::vector<double> quadrilateral_last_updated;
  std::vector<double> face_last_updated;
  bool on_the_fly_boundary = true;
  std::optional<std::filesystem::path> visual_debug_output_root_;
  std::optional<double> flip_polynomial_dump_target_time_;
  std::optional<size_t> flip_polynomial_dump_target_half_edge_;

  // crossing-related data (see public `CrossingData` forward declaration above).

  // Owned by `CrossingEventManager`. This is only an alias reference so the existing
  // code can keep using the name `crossing_data`.
  CrossingData& crossing_data;

  glm::dvec3 computeVoronoiVertexHomogenous(size_t voronoi_vertex_id, double t) const;

  void reassignVoronoiVerticesOnBoundary(size_t he_id, double t);

  void reassignVoronoiVerticesInQuadrilateral(
    size_t quad_index, double t, const std::map<size_t, size_t>& pre_flip_quad_faces);

  void precomputeStep(double t);

  /**
   * Piecewise-linear site motion on [@p section, @p section + 1] using @ref StrandTree::getPiecePolynomial with the
   * same per-strand reference branch as @ref getPointAt for times in that interval.
   */
  Trajectory<2> getSitePiecePolynomial(size_t strand_id, size_t section, double schedule_time) const;

  void growGraphSlotArrays();
  size_t findContainingTriForVoronoiVertex(size_t voronoi_vertex_id, double t) const;
  void initializeFaceState(size_t face_index, double t);
  void initializeNewFacesAfterGraphUpdate(double t, size_t first_new_face_slot);
  void clearPendingBranchSplits();
  bool pendingSplitSeamsAreConvex(size_t parent_component_id, double t) const;
  void maybeScheduleSeparationOrApplyPendingSplit(size_t parent_component_id, double split_time);
  void handleSeparationEventAtTime(size_t parent_component_id, double t);
  void applyPendingComponentGraphSplit(double t);
  void startSeparationSchedule(size_t parent_component_id, double segment_start_time);
  void continueSeparationSchedule(size_t parent_component_id, double segment_start_time);
  void recomputeEventsAfterSeparationTrajectory(size_t parent_component_id, double t);
  void collectSeparationRecomputeTargets(size_t parent_component_id, std::unordered_set<size_t>& affected_quads,
    std::unordered_set<size_t>& affected_faces) const;
  const PendingBranchSplit* activeSeparationForStrand(size_t strand_id) const;
  glm::dvec2 separationOffsetAt(size_t strand_id, double t) const;
  Trajectory<2> addSeparationOffsetToPiecePolynomial(
    const Trajectory<2>& base, size_t strand_id, size_t section) const;
  /// Strand ids whose min input branch defines the motion frame for @p strand_id.
  void collectReferenceBranchStrandPool(size_t strand_id, std::vector<size_t>& pool) const;
  void updateRuntimeBranchMapFromInputBranches(double t);
  void updateRuntimeBranchMapFromLiveGraph(double t, const char* trigger);
  size_t allocateRuntimeBranch();
  void markRuntimeBranchDead(size_t runtime_branch_id, double t, const char* reason,
    std::optional<size_t> input_branch_id = std::nullopt);
  void logRuntimeBranchEvent(double t, const std::string& event_line);
  bool runtime_branch_log_header_written_ = false;
  void retireFinishedInputBranches(double t);
  void validateFinishedInputBranchMatchesRuntime(size_t section, size_t input_branch_id) const;
  std::vector<std::vector<size_t>> extractGraphConnectedComponents() const;
  std::vector<size_t> extractGraphConnectedComponent(size_t u, std::vector<bool>& visited) const;
  size_t getRuntimeBranchIdForFace(size_t face_id) const;
  void onGraphRetriangulated(double t, size_t prev_face_slots, size_t prev_he_slots);
  void onGraphCutApplied(double t, size_t prev_face_slots, size_t prev_he_slots, bool update_runtime_branch_map = true);
  bool isStrandLiveInGraph(size_t strand_id) const;

  void handleEvents();

  /// Enqueues one @ref SubdivisionEvent per entry in @ref subdivision_schedule_ (called once from @ref compute).
  void enqueueScheduledSubdivisionEvents();

  size_t getBranchIndex(size_t strand_id, size_t t) const;

  const std::vector<std::vector<size_t>>& getBranches(size_t t) const;

  const std::vector<size_t>& getBranchStrands(size_t t, size_t branch_id);

  std::vector<double> findEvents(
    Polynomial& event_trigger, double min_fraction, bool only_positive_to_negative = false);

 public:
  KineticDelaunay(const StrandTree& branch_trajs, double cutoff, bool add_dummy_splines);
  ~KineticDelaunay();

  bool isDummyBoundary(size_t v) const;

  bool computeBoundaryOnTheFly() const;

  glm::dvec2 getPointAt(size_t v, double t) const;

  glm::dvec2 getPointAt(double t, size_t v) const;

  /** Branch frame used by getPointAt for sites in the same component as @p strand_id. */
  size_t getReferenceBranch(size_t strand_id, double t) const;

  /// Smallest live strand id sharing @p strand_id's runtime branch; used for transform frame lookup.
  size_t representativeStrandIdForRuntimeBranch(size_t strand_id) const;

  /// One shared reference branch for all @p strand_ids (e.g. flip quadrilateral vertices).
  size_t getSharedReferenceBranchForStrands(const std::vector<size_t>& strand_ids, double branch_lookup_time) const;

  glm::dvec2 getPointAtWithReferenceBranch(size_t v, double t, size_t reference_branch) const;

  Trajectory<2> getSitePiecePolynomialWithReferenceBranch(
    size_t strand_id, size_t section, size_t reference_branch) const;

  std::vector<glm::dvec2> getPointsAt(double t) const;

  glm::dvec3 getPointInObjectSpace(size_t v, double t) const;

  const StrandTree& getStrandTree() const;

  /// Input branch ids whose last real-strand section is @p t (same criterion as @ref retireFinishedInputBranches).
  std::vector<size_t> inputBranchesFinishingAtSection(double t) const;

  void setVisualDebugOutputRoot(const std::filesystem::path& root);
  const std::optional<std::filesystem::path>& getVisualDebugOutputRoot() const;
  /// Path to @c branches.txt under @ref getVisualDebugOutputRoot when visual debug is enabled.
  std::optional<std::filesystem::path> getRuntimeBranchLogPath() const;

  void setFlipPolynomialDumpTargetTime(std::optional<double> target_time);
  const std::optional<double>& getFlipPolynomialDumpTargetTime() const;

  void setFlipPolynomialDumpTargetHalfEdge(std::optional<size_t> half_edge_id);
  const std::optional<size_t>& getFlipPolynomialDumpTargetHalfEdge() const;

  void computeComponentData(double t);

  /// True once @ref component_data reflects the current @ref HalfEdgeDelaunayGraph after the latest section
  /// retriangulation. While false, connected components may already be split in component data but the mesh
  /// is still the pre-split triangulation (see @ref prev_component_count).
  bool isGraphRetriangulatedForComponents() const;

  void setComponentSplitPolicy(ComponentSplitPolicy policy) { component_split_policy_ = policy; }
  ComponentSplitPolicy getComponentSplitPolicy() const { return component_split_policy_; }
  void setSeparationOffsetScale(double scale) { separation_offset_scale_ = scale; }
  double getSeparationOffsetScale() const { return separation_offset_scale_; }

  /// Record a pending branch split until the next graph cut/retriangulation.
  void notePendingBranchSplit(size_t parent_component_id, double split_time,
    const std::vector<size_t>& pre_split_parent_strands, const std::vector<std::vector<size_t>>& new_components,
    const std::vector<size_t>& split_component_ids);

  /// Pending split data keyed by the parent component id, if recorded.
  std::optional<PendingBranchSplit> getPendingBranchSplit(size_t parent_component_id) const;

  /**
   * Seam outline between @p component_id and its pending-split partner components.
   * Walks live cross-component Delaunay edges that will be severed by @ref HalfEdgeDelaunayGraph::applyRuntimeBranchSplit.
   */
  std::vector<BoundaryPoint> collectFutureRuntimeBranchOutline(size_t component_id, double t) const;

  /// One outer outline per future runtime branch listed in the pending split metadata.
  std::vector<std::vector<BoundaryPoint>> collectPendingSplitBranchOutlines(
    size_t parent_component_id, double t) const;

  /// Delaunay simplices whose flip/radius/crossing events are recomputed during separation scheduling.
  VisualDebugHighlight buildSeparationRecomputeHighlight(size_t parent_component_id) const;

  /// Base-to-offset segments for separated strands in an active pending split.
  std::vector<HalfEdgeDelaunayGraphToSVG::SeparationOffsetSegment> collectSeparationOffsetSegments(
    size_t parent_component_id, double t) const;

  /// Smallest @ref StrandTree::getBranchIndex among live strands in @p component_id at @p branch_lookup_height.
  size_t minInputBranchForComponent(size_t component_id, size_t branch_lookup_height) const;

  /// Smallest @ref StrandTree::getBranchIndex among @p strand_ids at @p branch_lookup_height.
  size_t minInputBranchForStrands(const std::vector<size_t>& strand_ids, size_t branch_lookup_height) const;

  CrossingData& getCrossingDataMutable();
  const CrossingData& getCrossingData() const;

  /** Per-edge intersection list consistency; intended after @ref EventCallback::afterEvent (e.g. debug SVG export). */
  void validateCrossingIntersectionInvariants(const char* context, double t) const;

  /** Cached Voronoi-vertex list iterators must reference their own id in the containing triangle list. */
  void validateVoronoiVertexIteratorInvariants(const char* context, double t) const;

  std::vector<HalfEdgeDelaunayGraphToSVG::IntersectionDebugInfo> getCrossingIntersectionDebugData() const;

  const HalfEdgeDelaunayGraph& init(CallbackManager* callback_manager = nullptr);
  void registerSectionEventCallback(EventCallback* callback);
  void registerFlipEventCallback(EventCallback* callback);
  void registerRadiusEventCallback(EventCallback* callback);
  void registerCrossingEventCallback(EventCallback* callback);
  void registerSubdivisionEventCallback(EventCallback* callback);
  void registerSeparationEventCallback(EventCallback* callback);
  /// Pairs `(strand_id, t)` sorted by non-decreasing `t`; enqueued once when @ref compute builds the initial queue.
  void setSubdivisionSchedule(std::vector<std::pair<size_t, double>> schedule);
  void registerEventCallbacks(EventCallback* section_callback, EventCallback* flip_callback,
    EventCallback* radius_callback, EventCallback* crossing_callback);

  const HalfEdgeDelaunayGraph& getGraph() const;

  size_t getSectionCount() const;

  // Computes the Delaunay triangulation of the given splines
  void compute();

  std::vector<size_t> extractConnectedComponent(size_t u, std::vector<bool>& visited) const;

  const std::vector<glm::dvec2>& getDummyBoundary() const;

  std::vector<std::vector<size_t>> checkForSplit(const std::array<int, 3>& tri_vertices) const;

  std::vector<std::vector<size_t>> extractConnectedComponents() const;

  std::vector<BoundaryPoint> traverseBoundary(size_t start_he_id, double t) const;

  std::vector<std::vector<BoundaryPoint>> extractComponentBoundaries(
    const std::vector<size_t>& component, double t, std::vector<bool>& he_visited) const;

  std::vector<BoundaryPoint> extractComponentBoundary(const std::vector<size_t>& component, double t) const;

  const std::vector<bool>& getFacesInside() const;

  bool getFaceInside(size_t face_index) const;

  void setFaceInside(size_t face_index, bool value, double t);

  /** All three vertices share one input branch with exactly three strands at @p t. */
  bool isMinimalInputBranchTriangle(const std::array<int, 3>& vertices, double t) const;

  /**
   * Runtime branch of @p strand_id (see @ref runtime_branch_map_).
   * Not the inside-face component id and not the input @ref StrandTree branch index.
   */
  size_t getRuntimeBranchIdForStrand(size_t strand_id) const;

  /** Inside-face connected component id for @p strand_id. */
  size_t getInsideFaceComponentId(size_t strand_id) const;

  /** Lowest strand id in inside-face component @p component_id (BFS seed order). */
  size_t getComponentLowestStrandId(size_t component_id) const;

  /**
   * Runtime branch of a live half-edge, resolved from its incident face(s).
   * Returns 0 when the half-edge is dead or the branch cannot be resolved.
   */
  size_t getRuntimeBranchIdForHalfEdge(size_t half_edge_id) const;

  /** True iff the live Delaunay graph contains exactly one finite triangle in @p runtime_branch_id. */
  bool runtimeBranchHasSingleFiniteTriangle(size_t runtime_branch_id) const;

  const std::vector<size_t>& getRuntimeBranchMap() const { return runtime_branch_map_; }

  bool mustRemainInside(size_t face_index, double t) const;

  bool isOnComponentBoundary(size_t he_id) const;

  bool isOnComponentBoundaryOutside(size_t he_id) const;

  bool isOnFutureBranchSeamForComponent(
    size_t he_id, size_t component_id, const std::unordered_set<size_t>& partner_component_ids) const;
  size_t nextOnFutureBranchSeamId(
    size_t he_id, size_t component_id, const std::unordered_set<size_t>& partner_component_ids) const;
  std::vector<BoundaryPoint> traverseFutureBranchSeam(size_t start_he_id, size_t component_id,
    const std::unordered_set<size_t>& partner_component_ids, double t) const;
  std::optional<std::unordered_set<size_t>> seamPartnerComponentsFor(size_t component_id) const;

  size_t nextOnComponentBoundaryId(size_t he_id) const;

  // CrossingData accessors. Parameter `voronoi_vertex_id` is a Delaunay face index (dual circumcenter),
  // NOT a Voronoi cell / site id (Delaunay vertex index).
  bool isCrossingDataVoronoiVertexRegistered(size_t voronoi_vertex_id) const;
  void requireLiveRegisteredVoronoiVertex(size_t voronoi_vertex_id, const char* context) const;
  size_t getCrossingDataContainingTriId(size_t voronoi_vertex_id) const;
  std::vector<size_t> getCrossingDataVoronoiVerticesInTri(size_t tri_id) const;
  glm::dvec3 getVoronoiVertexHomogeneous(size_t voronoi_vertex_id, double t) const;

  /**
   * \brief Compute the (possibly clamped) Voronoi vertex position for the Delaunay edge represented by half_edge_id.
   *
   * For finite triangles this is the circumcenter; for infinite / hull cases, this returns a finite point obtained
   * by moving a neighboring circumcenter along a perpendicular direction so it can be used for meshing and
   * intersection computations.
   */
  glm::dvec3 computeVoronoiVertexClampedInfinity(size_t half_edge_id, double t) const;
  glm::dvec3 computeVoronoiVertexClampedInfinityWithReferenceBranch(
    size_t half_edge_id, double t, size_t reference_branch) const;
};
} // namespace kinDS
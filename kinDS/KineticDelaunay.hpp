#pragma once
#include "HalfEdgeDelaunayGraph.hpp"
#include "HalfEdgeDelaunayGraphToSVG.hpp"
#include "KineticAlgorithm.hpp"
#include "Logger.hpp"
#include "Polynomial.hpp"
#include "ProgressBar.hpp"
#include "StrandTree.hpp"
#include <format>
#include <glm/gtx/exterior_product.hpp>
#include <glm/gtx/string_cast.hpp>
#include <map>
#include <memory>
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

  class FlipEvent;
  class RadiusEvent;
  class CrossingEvent;

  class CallbackManager
  {
   public:
    virtual ~CallbackManager() = default;

    virtual void init() { }
    virtual void finalize(double t) { }
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
  std::vector<double> quadrilateral_last_updated;
  std::vector<double> face_last_updated;
  bool on_the_fly_boundary = true;

  // crossing-related data (see public `CrossingData` forward declaration above).

  // Owned by `CrossingEventManager`. This is only an alias reference so the existing
  // code can keep using the name `crossing_data`.
  CrossingData& crossing_data;

  glm::dvec3 computeVoronoiVertexHomogenous(size_t voronoi_vertex_id, double t) const;

  void reassignVoronoiVerticesOnBoundary(size_t he_id, double t);

  void reassignVoronoiVerticesInQuadrilateral(
    size_t quad_index, double t, const std::map<size_t, size_t>& pre_flip_quad_faces);

  void precomputeStep(double t);

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

  bool isDummyBoundary(size_t v);

  bool computeBoundaryOnTheFly() const;

  glm::dvec2 getPointAt(size_t v, double t) const;

  glm::dvec2 getPointAt(double t, size_t v) const;

  std::vector<glm::dvec2> getPointsAt(double t) const;

  glm::dvec3 getPointInObjectSpace(size_t v, double t) const;

  const StrandTree& getStrandTree() const;

  void computeComponentData(double t);

  const CrossingData& getCrossingData() const;
  /// Rebuild all Voronoi–Delaunay edge intersection lists at time @p t (invalidates existing `EdgeIntersectionRef`s).
  void recomputeEdgeIntersections(double t);
  std::vector<std::array<size_t, 4>> getCrossingIntersectionDebugData() const;

  const HalfEdgeDelaunayGraph& init(CallbackManager* callback_manager = nullptr);
  void registerSectionEventCallback(EventCallback* callback);
  void registerFlipEventCallback(EventCallback* callback);
  void registerRadiusEventCallback(EventCallback* callback);
  void registerCrossingEventCallback(EventCallback* callback);
  void registerSubdivisionEventCallback(EventCallback* callback);
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

  void setFaceInside(size_t face_index, bool value);

  bool isOnComponentBoundary(size_t he_id) const;

  bool isOnComponentBoundaryOutside(size_t he_id) const;

  size_t nextOnComponentBoundaryId(size_t he_id) const;

  // Getters for CrossingData (for testing/validation)
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
};
} // namespace kinDS
#pragma once
#include "HalfEdgeDelaunayGraph.hpp"
#include "HalfEdgeDelaunayGraphToSVG.hpp"
#include "Logger.hpp"
#include "Polynomial.hpp"
#include "ProgressBar.hpp"
#include "StrandTree.hpp"
#include <format>
#include <glm/gtx/exterior_product.hpp>
#include <glm/gtx/string_cast.hpp>
#include <map>
#include <queue>
#include <string>

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
  class Event
  {
   public:
    double time; // Time of the event
    size_t half_edge_id; // Half-edge index associated with the event
    double creation_time; // Time when the event was created, used do check validity after a quadrilateral is updated
    glm::dvec2 position; // Position of the event
    size_t voronoi_vertex_id; // Voronoi vertex index associated with the event, only used for crossing events

    enum Type
    {
      FLIP,
      RADIUS,
      CROSSING
    } type;

    Event(double t, size_t he_id, double creation_time, glm::dvec2 position, Type type)
      : time(t)
      , half_edge_id(he_id)
      , creation_time(creation_time)
      , position(position)
      , type(type)
    {
      if (type == Event::CROSSING)
      {
        throw std::invalid_argument("Crossing events must be created with a voronoi vertex id");
      }
    }

    Event(double t, size_t he_id, double creation_time, glm::dvec2 position, size_t voronoi_vertex_id, Type type)
      : time(t)
      , half_edge_id(he_id)
      , creation_time(creation_time)
      , position(position)
      , type(type)
      , voronoi_vertex_id(voronoi_vertex_id)
    {
    }

    bool operator<(const Event& other) const
    {
      return time > other.time; // For priority queue, we want the earliest event first
    }
  };

  // EventHandler class, inherit from this class to handle events in the KineticDelaunay algorithm
  class EventHandler
  {
   public:
    virtual ~EventHandler() = default;
    /**
     * \brief Handle a FLIP event before it is processed, i.e. before any edges are swapped
     *
     * @param e The event to handle.
     */
    virtual void beforeFlipEvent(Event& e) { }

    /**
     * \brief Handle a FLIP event after it is processed, i.e. after edges are swapped
     *
     * @param e The event to handle.
     */
    virtual void afterFlipEvent(Event& e) { }

    /**
     * \brief Handle a RADIUS event before it is processed
     *
     * @param e The event to handle
     */
    virtual void beforeRadiusEvent(Event& e) { }

    /**
     * \brief Handle a RADIUS event after it is processed
     *
     * @param e The event to handle
     */
    virtual void afterRadiusEvent(Event& e) { }

    virtual void beforeCrossingEvent(Event& e) { }

    virtual void afterCrossingEvent(Event& e) { }

    virtual void betweenSections(size_t index) { }

    /**
     * \brief initialize event handler.
     */
    virtual void init() { }

    /**
     * \brief Finalize after all events have been handled
     */
    virtual void finalize(double t) { }
  };

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

  std::pair<double, double> delaunayVoronoiEdgeIntersection(size_t delaunay_edge_id, size_t voronoi_edge_id, double t) const;

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
  typedef std::priority_queue<Event> EventQueue;

  StrandTree branch_trajs;
  HalfEdgeDelaunayGraph graph;
  EventQueue events;
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

  // crossing-related data
  struct CrossingData
  {
    struct VoronoiDelaunayEdgeIntersection;
    typedef std::list<VoronoiDelaunayEdgeIntersection>::iterator EdgeIntersectionRef;

   private:
    std::vector<size_t> voronoi_vertex_to_containing_tri_id;
    std::vector<std::list<size_t>> tri_id_to_voronoi_vertices;
    std::vector<std::list<size_t>::iterator> voronoi_vertex_to_iterator;

   public:
    std::list<VoronoiDelaunayEdgeIntersection> edge_intersections;
    std::vector<std::list<EdgeIntersectionRef>> voronoi_edge_intersections;
    std::vector<std::list<EdgeIntersectionRef>> delaunay_edge_intersections;


    std::vector<double> last_crossing;

    struct VoronoiDelaunayEdgeIntersection {
      std::list<EdgeIntersectionRef>::iterator delaunay_ref;
      std::list<EdgeIntersectionRef>::iterator voronoi_ref;
      size_t delaunay_edge_id;
      size_t voronoi_edge_id;
      double delaunay_edge_param;
    };

    void init(size_t face_count)
    {
      voronoi_vertex_to_containing_tri_id.clear();
      voronoi_vertex_to_containing_tri_id.resize(face_count, -1);

      tri_id_to_voronoi_vertices.clear();
      tri_id_to_voronoi_vertices.resize(face_count);

      voronoi_vertex_to_iterator.clear();
      voronoi_vertex_to_iterator.resize(face_count);

      last_crossing.clear();
      last_crossing.resize(face_count, 0.0);
    }

    void setVoronoiVertexTriId(size_t voronoi_vertex_id, size_t tri_id)
    {
      voronoi_vertex_to_containing_tri_id[voronoi_vertex_id] = tri_id;
      voronoi_vertex_to_iterator[voronoi_vertex_id]
        = tri_id_to_voronoi_vertices[tri_id].emplace(tri_id_to_voronoi_vertices[tri_id].end(), voronoi_vertex_id);
    }

    void moveVertex(size_t voronoi_vertex_id, size_t target_tri_id, double t)
    {
      std::list<size_t>::iterator v_it = voronoi_vertex_to_iterator[voronoi_vertex_id];
      size_t current_tri_id = voronoi_vertex_to_containing_tri_id[voronoi_vertex_id];

      if (current_tri_id == target_tri_id)
      {
        return;
      }

      tri_id_to_voronoi_vertices[current_tri_id].erase(v_it);

      setVoronoiVertexTriId(voronoi_vertex_id, target_tri_id);

      KINDS_DEBUG("Voronoi vertex " << voronoi_vertex_id << " moved from triangle " << current_tri_id << " to "
                                    << target_tri_id << " at t = " << t);
    }

    size_t getContainingTriId(size_t voronoi_vertex_id) const
    {
      return voronoi_vertex_to_containing_tri_id[voronoi_vertex_id];
    }

    // Note: we copy this into a vector because we need this for reassigning Voronoi vertices in a quadrilateral and we
    // will be modifying the underlying list while iterating
    std::vector<size_t> getVoronoiVerticesInTri(size_t tri_id) const
    {
      std::vector<size_t> result(tri_id_to_voronoi_vertices[tri_id].begin(), tri_id_to_voronoi_vertices[tri_id].end());
      return result;
    }

    void computeEdgeIntersections(const KineticDelaunay& kd, double t);

    // Update Voronoi–Delaunay edge intersections after a single crossing event.
    void updateAfterCrossingEvent(const KineticDelaunay& kd, const Event& e);

    // Remove a single intersection from all three data structures (global list,
    // per-Voronoi-edge list, and per-Delaunay-edge list).
    void removeIntersection(EdgeIntersectionRef intersection_ref);

  } crossing_data;

  glm::dvec3 computeVoronoiVertexHomogenous(size_t voronoi_vertex_id, double t) const;

  void computeCrossingEvents(double t, size_t voronoi_vertex_id);

  void reassignVoronoiVerticesOnBoundary(size_t he_id, double t);

  void reassignVoronoiVerticesInQuadrilateral(
    size_t quad_index, double t, const std::map<size_t, size_t>& pre_flip_quad_faces);

  void computeRadiusEvents(double t, size_t he_id);

  void computeFlipEvents(double t, size_t quad_id);

  void precomputeStep(double t);

  void handleFlipEvent(EventHandler& event_handler, Event& event);

  void handleRadiusEvent(EventHandler& event_handler, Event& event);

  void handleCrossingEvent(EventHandler& event_handler, Event& event);

  void handleEvents(EventHandler& event_handler);

  size_t getBranchIndex(size_t strand_id, size_t t) const;

  const std::vector<std::vector<size_t>>& getBranches(size_t t) const;

  const std::vector<size_t>& getBranchStrands(size_t t, size_t branch_id);

  std::vector<double> findEvents(Polynomial& event_trigger, double min_fraction, bool only_positive_to_negative = false);

 public:
  KineticDelaunay(const StrandTree& branch_trajs, double cutoff, bool add_dummy_splines);

  bool isDummyBoundary(size_t v);

  bool computeBoundaryOnTheFly() const;

  glm::dvec2 getPointAt(size_t v, double t) const;

  glm::dvec2 getPointAt(double t, size_t v) const;

  std::vector<glm::dvec2> getPointsAt(double t) const;

  glm::dvec3 getPointInObjectSpace(size_t v, double t) const;

  const StrandTree& getStrandTree() const;

  void computeComponentData(double t);

  const HalfEdgeDelaunayGraph& init();

  const HalfEdgeDelaunayGraph& advanceOneSection(EventHandler& event_handler);

  const HalfEdgeDelaunayGraph& getGraph() const;

  size_t getSectionCount() const;

  // Computes the Delaunay triangulation of the given splines
  void compute(EventHandler& event_handler);

  std::vector<size_t> extractConnectedComponent(size_t u, std::vector<bool>& visited) const;

  const std::vector<glm::dvec2>& getDummyBoundary() const;

  std::vector<std::vector<size_t>> checkForSplit(const std::array<int, 3>& tri_vertices) const;

  std::vector<std::vector<size_t>> extractConnectedComponents() const;

  std::vector<BoundaryPoint> traverseBoundary(size_t start_he_id, double t) const;

  std::vector<std::vector<BoundaryPoint>> extractComponentBoundaries(
    const std::vector<size_t>& component, double t, std::vector<bool>& he_visited) const;

  std::vector<BoundaryPoint> extractComponentBoundary(const std::vector<size_t>& component, double t) const;

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
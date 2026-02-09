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

    enum Type
    {
      SWAP,
      BOUNDARY,
      RIGHT_ANGLED
    } type;

    Event(double t, size_t he_id, double creation_time, glm::dvec2 position, Type type)
      : time(t)
      , half_edge_id(he_id)
      , creation_time(creation_time)
      , position(position)
      , type(type)
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
     * \brief Handle a SWAP event before it is processed, i.e. before any edges are swapped
     *
     * @param e The event to handle.
     */
    virtual void beforeEvent(Event& e) { }

    /**
     * \brief Handle a SWAP event after it is processed, i.e. after edges are swapped
     *
     * @param e The event to handle.
     */
    virtual void afterEvent(Event& e) { }

    /**
     * \brief Handle a BOUNDARY event before it is processed
     *
     * @param e The event to handle
     */
    virtual void beforeBoundaryEvent(Event& e) { }

    /**
     * \brief Handle a BOUNDARY event after it is processed
     *
     * @param e The event to handle
     */
    virtual void afterBoundaryEvent(Event& e) { }

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

  void computeRightAngleEvents(double t, size_t he_id);

  void computeBoundaryEvents(double t, size_t he_id);

  void computeSwapEvents(double t, size_t quad_id);

  void precomputeStep(double t);

  void handleSwapEvent(EventHandler& event_handler, Event& event);

  void handleBoundaryEvent(EventHandler& event_handler, Event& event);

  void handleEvents(EventHandler& event_handler);

  size_t getBranchIndex(size_t strand_id, size_t t) const;

  const std::vector<std::vector<size_t>>& getBranches(size_t t) const;

  const std::vector<size_t>& getBranchStrands(size_t t, size_t branch_id);

 public:
  KineticDelaunay(const StrandTree& branch_trajs, double cutoff, bool add_dummy_splines);

  bool isDummyBoundary(size_t v);

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
};
} // namespace kinDS
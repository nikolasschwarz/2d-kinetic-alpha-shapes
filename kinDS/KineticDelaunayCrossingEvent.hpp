#pragma once

#include "KineticDelaunay.hpp"

#include <glm/glm.hpp>
#include <optional>
#include <stdexcept>
#include <string>

namespace kinDS
{
struct KineticDelaunay::CrossingData
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

  struct VoronoiDelaunayEdgeIntersection
  {
    std::list<EdgeIntersectionRef>::iterator delaunay_ref;
    std::list<EdgeIntersectionRef>::iterator voronoi_ref;
    size_t delaunay_edge_id;
    size_t voronoi_edge_id;
    double delaunay_edge_param;
    // SegmentBuilder boundary-interval mesh linkage along one Delaunay edge.
    size_t prev_segment_mesh_pair_index = static_cast<size_t>(-1);
    size_t next_segment_mesh_pair_index = static_cast<size_t>(-1);
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

  const std::vector<size_t>& getContainingTriIds() const { return voronoi_vertex_to_containing_tri_id; }

  // Note: we copy this into a vector because we need this for reassigning Voronoi vertices in a quadrilateral and we
  // will be modifying the underlying list while iterating
  std::vector<size_t> getVoronoiVerticesInTri(size_t tri_id) const
  {
    std::vector<size_t> result(tri_id_to_voronoi_vertices[tri_id].begin(), tri_id_to_voronoi_vertices[tri_id].end());
    return result;
  }

  void computeEdgeIntersections(const KineticDelaunay& kd, double t);

  // Update Voronoi–Delaunay edge intersections after a single crossing event.
  void updateAfterCrossingEvent(const KineticDelaunay& kd, const KineticDelaunay::CrossingEvent& e);

  // Remove a single intersection from all three data structures (global list,
  // per-Voronoi-edge list, and per-Delaunay-edge list).
  void removeIntersection(EdgeIntersectionRef intersection_ref);
};

class KineticDelaunay::CrossingEvent final : public KineticDelaunay::Event
{
 public:
  KineticDelaunay* kd_;
  size_t half_edge_id;
  glm::dvec2 position;
  size_t voronoi_vertex_id;

  CrossingEvent(
    KineticDelaunay* kd, double t, size_t he_id, double creation_time, glm::dvec2 position, size_t voronoi_vertex_id)
    : KineticDelaunay::Event(t, creation_time, 40u)
    , kd_(kd)
    , half_edge_id(he_id)
    , position(position)
    , voronoi_vertex_id(voronoi_vertex_id)
  {
  }

  void handleEvent() override;
};

class KineticDelaunay::CrossingEventManager final : public KineticDelaunay::EventManager
{
 public:
  explicit CrossingEventManager(KineticDelaunay* kd)
    : kd_(kd)
  {
  }

  void computeEvents(double t, size_t event_id) override;

  // Expose crossing data to the KineticDelaunay orchestration code and to events.
  CrossingData& getCrossingDataMutable() { return crossing_data_; }
  const CrossingData& getCrossingData() const { return crossing_data_; }

 private:
  KineticDelaunay* kd_;
  CrossingData crossing_data_;
};

inline void KineticDelaunay::CrossingEvent::handleEvent()
{
  auto* kd = kd_;
  if (!kd)
  {
    throw std::runtime_error("CrossingEvent has no KineticDelaunay pointer");
  }

  auto& graph = kd->graph;

  // Check if the event is still valid
  if (creation_time < kd->crossing_data.last_crossing[voronoi_vertex_id])
  {
    return;
  }

  size_t containing_tri_id = kd->crossing_data.getContainingTriId(voronoi_vertex_id);

  // The event is also outdated if the face has been updated in a flip event
  if (creation_time < kd->face_last_updated[containing_tri_id])
  {
    return;
  }

  auto* event_handler = kd->crossing_event_manager_->getCallback();
  if (event_handler)
  {
    event_handler->beforeEvent(*this);
  }

  kd->crossing_data.last_crossing[voronoi_vertex_id] = occurrence_time;

  KINDS_DEBUG("Processing crossing event at time " << occurrence_time << " for Voronoi vertex ID " << voronoi_vertex_id
                                                   << " crossing half-edge ID " << half_edge_id);

  // move to neighboring triangle
  KINDS_DEBUG("Moving Voronoi vertex " << voronoi_vertex_id << " from triangle " << containing_tri_id << " to triangle "
                                       << graph.getHalfEdges()[half_edge_id ^ 1].face);
  kd->crossing_data.moveVertex(voronoi_vertex_id, graph.getHalfEdges()[half_edge_id ^ 1].face, occurrence_time);

  // Update Voronoi–Delaunay edge intersections stored in crossing_data in response to this crossing.
  kd->crossing_data.updateAfterCrossingEvent(*kd, *this);

  if (event_handler)
  {
    event_handler->afterEvent(*this);
  }

  // Re-compute crossing events for this Voronoi vertex
  kd->crossing_event_manager_->computeEvents(occurrence_time, voronoi_vertex_id);
}

/** Empty optional formats as "null". Log form is V{idx}xD{idx} for list positions along voronoi/delaunay edge lists. */
std::string formatCrossingIntersectionForLog(const KineticDelaunay& kd,
  std::optional<KineticDelaunay::CrossingData::EdgeIntersectionRef> intersection);

bool tryComputeCrossingIntersectionPosition2D(const KineticDelaunay& kd,
  std::optional<KineticDelaunay::CrossingData::EdgeIntersectionRef> intersection, double t, glm::dvec2& out_xy);

/** No-op if @p intersection is empty; otherwise log KINDS_ERROR on mismatch with expected dual edge / half-edge. */
void validateClosingCapCrossingRef(const KineticDelaunay& kd, const char* context_msg,
  std::optional<KineticDelaunay::CrossingData::EdgeIntersectionRef> intersection, size_t expected_voronoi_edge_id,
  int delaunay_half_edge_id);

} // namespace kinDS

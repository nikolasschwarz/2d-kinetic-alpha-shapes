#pragma once

#include "KineticDelaunay.hpp"
#include "KineticDelaunayRadiusEvent.hpp"

#include <limits>
#include <map>

namespace kinDS
{
class KineticDelaunay::FlipEvent final : public KineticDelaunay::Event
{
 public:
  KineticDelaunay* kd_;
  size_t half_edge_id;
  glm::dvec2 position;

  FlipEvent(KineticDelaunay* kd, double t, size_t he_id, double creation_time, glm::dvec2 position)
    : KineticDelaunay::Event(t, creation_time, 20u)
    , kd_(kd)
    , half_edge_id(he_id)
    , position(position)
  {
  }

  void handleEvent() override;
};

class KineticDelaunay::FlipEventManager final : public KineticDelaunay::EventManager
{
 public:
  explicit FlipEventManager(KineticDelaunay* kd)
    : kd_(kd)
  {
  }

  void computeEvents(double t, size_t event_id) override;

 private:
  KineticDelaunay* kd_;
};

inline void KineticDelaunay::FlipEvent::handleEvent()
{
  auto* kd = kd_;
  if (!kd)
  {
    throw std::runtime_error("FlipEvent has no KineticDelaunay pointer");
  }

  auto& graph = kd->graph;

  // Check if the event is still valid
  if (creation_time < kd->quadrilateral_last_updated[half_edge_id / 2])
  {
    return;
  }

  // Before modifying the topology, store the face id for each half-edge in the quadrilateral
  // (three per triangle) so we can reason about pre-flip topology if needed.
  std::map<size_t, size_t> pre_flip_quad_faces;
  {
    size_t he0 = half_edge_id;
    size_t he1 = graph.getHalfEdges()[he0].next;
    size_t he2 = graph.getHalfEdges()[he1].next;
    size_t he3 = he0 ^ 1;
    size_t he4 = graph.getHalfEdges()[he3].next;
    size_t he5 = graph.getHalfEdges()[he4].next;

    pre_flip_quad_faces[he0] = graph.getHalfEdges()[he0].face;
    pre_flip_quad_faces[he1] = graph.getHalfEdges()[he1].face;
    pre_flip_quad_faces[he2] = graph.getHalfEdges()[he2].face;
    pre_flip_quad_faces[he3] = graph.getHalfEdges()[he3].face;
    pre_flip_quad_faces[he4] = graph.getHalfEdges()[he4].face;
    pre_flip_quad_faces[he5] = graph.getHalfEdges()[he5].face;
  }

  // Process the event at the given time
  size_t face_id = graph.getHalfEdges()[half_edge_id].face;
  size_t twin_face_id = graph.getHalfEdges()[half_edge_id ^ 1].face;
  KINDS_DEBUG("Processing flip event at time " << occurrence_time << " for half-edge ID " << half_edge_id
                                               << ". Faces inside " << kd->face_inside[face_id] << " | "
                                               << kd->face_inside[twin_face_id]);

  auto* event_handler = kd->flip_event_manager_->getCallback();
  if (event_handler)
  {
    event_handler->beforeEvent(*this);
  }

  // Faces swapped to the inside start out with an infinite circumradius, therefore their state depends on the cutoff
  if (graph.getHalfEdges()[half_edge_id].origin == -1)
  {
    kd->face_inside[twin_face_id] = (kd->cutoff == std::numeric_limits<double>::infinity());
  }

  if (graph.getHalfEdges()[half_edge_id ^ 1].origin == -1)
  {
    kd->face_inside[face_id] = (kd->cutoff == std::numeric_limits<double>::infinity());
  }

  graph.flipEdge(half_edge_id);

  // one of the triangles might have been swapped outside
  auto tri_verts1 = graph.adjacentTriangleVertices(half_edge_id);
  for (auto& v : tri_verts1)
  {
    if (v == -1)
    {
      size_t swapped_face_id = graph.getHalfEdges()[half_edge_id].face;
      kd->setFaceInside(swapped_face_id, false);
    }
  }

  auto tri_verts2 = graph.adjacentTriangleVertices(half_edge_id ^ 1);
  for (auto& v : tri_verts2)
  {
    if (v == -1)
    {
      size_t swapped_face_id = graph.getHalfEdges()[half_edge_id ^ 1].face;
      kd->setFaceInside(swapped_face_id, false);
    }
  }

  // After flipping the edge, we need to recompute the events for all surrounding half-edges
  size_t next1 = graph.getHalfEdges()[half_edge_id].next;
  size_t next2 = graph.getHalfEdges()[next1].next;

  size_t twin_next1 = graph.getHalfEdges()[half_edge_id ^ 1].next;
  size_t twin_next2 = graph.getHalfEdges()[twin_next1].next;

  kd->flip_event_manager_->computeEvents(occurrence_time, next1 / 2);
  kd->quadrilateral_last_updated[next1 / 2] = occurrence_time;

  kd->flip_event_manager_->computeEvents(occurrence_time, next2 / 2);
  kd->quadrilateral_last_updated[next2 / 2] = occurrence_time;

  kd->flip_event_manager_->computeEvents(occurrence_time, twin_next1 / 2);
  kd->quadrilateral_last_updated[twin_next1 / 2] = occurrence_time;

  kd->flip_event_manager_->computeEvents(occurrence_time, twin_next2 / 2);
  kd->quadrilateral_last_updated[twin_next2 / 2] = occurrence_time;

  // re-compute radius events for both triangles
  kd->radius_event_manager_->computeEvents(occurrence_time, half_edge_id);
  kd->face_last_updated[face_id] = occurrence_time;

  kd->radius_event_manager_->computeEvents(occurrence_time, half_edge_id ^ 1);
  kd->face_last_updated[twin_face_id] = occurrence_time;

  // trigger re-assignment of voronoi vertices needed for crossing events
  if (!graph.isOnConvexBoundary(half_edge_id))
  {
    kd->reassignVoronoiVerticesInQuadrilateral(half_edge_id / 2, occurrence_time, pre_flip_quad_faces);
  }
  else
  {
    kd->reassignVoronoiVerticesOnBoundary(half_edge_id, occurrence_time);
  }

  if (event_handler)
  {
    event_handler->afterEvent(*this);
  }
}

} // namespace kinDS

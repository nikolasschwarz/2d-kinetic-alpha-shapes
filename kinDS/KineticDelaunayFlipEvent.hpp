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
    size_t he1 = graph.halfEdge(he0).next;
    size_t he2 = graph.halfEdge(he1).next;
    size_t he3 = he0 ^ 1;
    size_t he4 = graph.halfEdge(he3).next;
    size_t he5 = graph.halfEdge(he4).next;

    pre_flip_quad_faces[he0] = graph.halfEdge(he0).face;
    pre_flip_quad_faces[he1] = graph.halfEdge(he1).face;
    pre_flip_quad_faces[he2] = graph.halfEdge(he2).face;
    pre_flip_quad_faces[he3] = graph.halfEdge(he3).face;
    pre_flip_quad_faces[he4] = graph.halfEdge(he4).face;
    pre_flip_quad_faces[he5] = graph.halfEdge(he5).face;
  }

  // Process the event at the given time
  size_t face_id = graph.halfEdge(half_edge_id).face;
  size_t twin_face_id = graph.halfEdge(half_edge_id ^ 1).face;
  KINDS_DEBUG("Processing flip event at time " << occurrence_time << " for half-edge ID " << half_edge_id
                                               << ". Faces inside " << kd->face_inside[face_id] << " | "
                                               << kd->face_inside[twin_face_id]);

  auto* event_handler = kd->flip_event_manager_->getCallback();
  if (event_handler)
  {
    event_handler->beforeEvent(*this);
  }

  // Faces swapped to the inside start out with an infinite circumradius, therefore their state depends on the cutoff
  if (graph.halfEdge(half_edge_id).origin == -1)
  {
    kd->face_inside[twin_face_id] = (kd->cutoff == std::numeric_limits<double>::infinity());
  }

  if (graph.halfEdge(half_edge_id ^ 1).origin == -1)
  {
    kd->face_inside[face_id] = (kd->cutoff == std::numeric_limits<double>::infinity());
  }

  // Special case if there is only one triangle
  const size_t branch = kd->getRuntimeBranchIdForHalfEdge(half_edge_id);

  bool is_single_triangle = kd->runtimeBranchHasSingleFiniteTriangle(branch);
  
  if(is_single_triangle){
    // First determine which edge is inside the triangle and which is outside.
    size_t inside_edge_id;

    if(kd->isOnComponentBoundaryOutside(half_edge_id)){
      inside_edge_id = half_edge_id ^ 1;
    } else if(kd->isOnComponentBoundaryOutside(half_edge_id ^ 1)){
      inside_edge_id = half_edge_id;
    } else {
      throw std::runtime_error("Single triangle flip event: neither edge is on the component boundary!");
    }

    size_t opposite_vertex_id = graph.triangleOppositeVertex(inside_edge_id);
    if(opposite_vertex_id == -1){
      throw std::runtime_error("Single triangle flip event: opposite vertex is infinite!"); 
    }

    // Now find an infinite outgoing half-edge from the opposite vertex
    size_t other_flip_edge_id = -1;
    for(auto incident_he_id = graph.incidentEdgesBegin(opposite_vertex_id); incident_he_id != graph.incidentEdgesEnd(opposite_vertex_id); ++incident_he_id){
      if(graph.destination(*incident_he_id) == -1){
        other_flip_edge_id = *incident_he_id;
        break;
      }
    }

    if (other_flip_edge_id == static_cast<size_t>(-1))
    {
      throw std::runtime_error("Single triangle flip event: no infinite outgoing half-edge at opposite vertex");
    }

    // order shouldn't matter, so we do this edge flip first, then the other one
    graph.flipEdge(other_flip_edge_id);

    const auto is_finite_live_face = [&](size_t flipped_face_id) -> bool
    {
      if (!graph.isLiveFace(flipped_face_id))
      {
        return false;
      }
      const auto vertices = graph.getTriangleVertexIndices(flipped_face_id);
      return vertices[0] != -1 && vertices[1] != -1 && vertices[2] != -1;
    };

    const size_t flipped_face0 = static_cast<size_t>(graph.halfEdge(other_flip_edge_id).face);
    const size_t flipped_face1 = static_cast<size_t>(graph.halfEdge(other_flip_edge_id ^ 1).face);

    if (is_finite_live_face(flipped_face0))
    {
      kd->setFaceInside(flipped_face0, true, occurrence_time);
    }
    else if (is_finite_live_face(flipped_face1))
    {
      kd->setFaceInside(flipped_face1, true, occurrence_time);
    }
    else
    {
      throw std::runtime_error("Single triangle flip event: auxiliary flip did not produce a finite triangle");
    }
  }

  graph.flipEdge(half_edge_id);

  // one of the triangles might have been swapped outside
  auto tri_verts1 = graph.adjacentTriangleVertices(half_edge_id);
  for (auto& v : tri_verts1)
  {
    if (v == -1)
    {
      size_t swapped_face_id = graph.halfEdge(half_edge_id).face;
      kd->setFaceInside(swapped_face_id, false, occurrence_time);
    }
  }

  auto tri_verts2 = graph.adjacentTriangleVertices(half_edge_id ^ 1);
  for (auto& v : tri_verts2)
  {
    if (v == -1)
    {
      size_t swapped_face_id = graph.halfEdge(half_edge_id ^ 1).face;
      kd->setFaceInside(swapped_face_id, false, occurrence_time);
    }
  }

  // After flipping the edge, we need to recompute the events for all surrounding half-edges
  size_t next1 = graph.halfEdge(half_edge_id).next;
  size_t next2 = graph.halfEdge(next1).next;

  size_t twin_next1 = graph.halfEdge(half_edge_id ^ 1).next;
  size_t twin_next2 = graph.halfEdge(twin_next1).next;

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

  // After callbacks (e.g. debug SVG export); intersection lists must be consistent.
  kd->validateVoronoiVertexIteratorInvariants("FlipEvent:afterEvent", occurrence_time);
  kd->validateCrossingIntersectionInvariants("FlipEvent:afterEvent", occurrence_time);
}

} // namespace kinDS

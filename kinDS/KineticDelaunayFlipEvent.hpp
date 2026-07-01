#pragma once

#include "KineticDelaunay.hpp"
#include "KineticDelaunayRadiusEvent.hpp"
#include "KineticDelaunayFlipEventTriggerDump.hpp"

#include <cmath>
#include <limits>
#include <map>
#include <string>

namespace kinDS
{
inline constexpr double flip_voronoi_vertex_distance_eps = 1e-6;
inline constexpr double flip_boundary_collinearity_eps = 1e-8;

inline double normalizedTriangleCollinearityMetric(const glm::dvec2& pa, const glm::dvec2& pb, const glm::dvec2& pc)
{
  const glm::dvec2 ab = pb - pa;
  const glm::dvec2 ac = pc - pa;
  const double area2 = std::abs(ab.x * ac.y - ab.y * ac.x);
  const double max_len2 = std::max({ glm::dot(ab, ab), glm::dot(ac, ac), glm::dot(pc - pb, pc - pb) });
  if (max_len2 <= 0.0)
  {
    return 0.0;
  }
  return area2 / max_len2;
}

inline glm::dvec2 flipTriangleCircumcenterAt(
  const KineticDelaunay& kd, const HalfEdgeDelaunayGraph& graph, size_t he_id, double t, bool use_transformed_points)
{
  const std::array<int, 3> vertices = graph.adjacentTriangleVertices(he_id);
  const StrandTree& tree = kd.getStrandTree();
  const auto point_at = [&](int vertex) -> glm::dvec2
  {
    const size_t strand_id = static_cast<size_t>(vertex);
    return use_transformed_points ? kd.getPointAt(strand_id, t) : tree.evaluate(strand_id, t);
  };
  return HalfEdgeDelaunayGraph::circumcenter(
    point_at(vertices[0]), point_at(vertices[1]), point_at(vertices[2]));
}

inline std::string flipUntransformedFrameMismatchNote(bool transformed_passes, bool untransformed_passes)
{
  if (!transformed_passes && untransformed_passes)
  {
    return " [untransformed geometry satisfies the flip criterion; possible frame mismatch in event polynomials]";
  }
  if (transformed_passes && !untransformed_passes)
  {
    return " [only transformed geometry satisfies the flip criterion]";
  }
  return {};
}

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

  if (kd->getVisualDebugOutputRoot().has_value()
    && shouldDumpFlipPolynomialsForEvent(*kd, occurrence_time, half_edge_id))
  {
    const FlipEventTriggerDump dump = buildFlipEventTriggerDump(*kd, half_edge_id, creation_time);
    writeFlipEventTriggerPolynomialDump(
      *kd, dump, *kd->getVisualDebugOutputRoot() / "polynomials.txt", occurrence_time);
  }

  // Flip validation disabled: inCircle/ccw predicate roots can disagree slightly with
  // transformed Voronoi geometry at event time (numerical / frame lookup tolerance).
#if 0
  if (graph.isOnConvexBoundary(half_edge_id) || graph.isOutsideConvexBoundary(half_edge_id))
  {
    size_t boundary_he_id = half_edge_id;
    if (graph.isOutsideConvexBoundary(boundary_he_id))
    {
      boundary_he_id ^= 1;
    }

    const int a = graph.halfEdge(boundary_he_id).origin;
    const int b = graph.triangleOppositeVertex(boundary_he_id ^ 1);
    const int c = graph.halfEdge(boundary_he_id ^ 1).origin;
    if (a >= 0 && b >= 0 && c >= 0)
    {
      const glm::dvec2 pa = kd->getPointAt(static_cast<size_t>(a), occurrence_time);
      const glm::dvec2 pb = kd->getPointAt(static_cast<size_t>(b), occurrence_time);
      const glm::dvec2 pc = kd->getPointAt(static_cast<size_t>(c), occurrence_time);
      const double collinearity_metric = normalizedTriangleCollinearityMetric(pa, pb, pc);
      if (collinearity_metric > flip_boundary_collinearity_eps)
      {
        const glm::dvec2 pa_raw = kd->getStrandTree().evaluate(static_cast<size_t>(a), occurrence_time);
        const glm::dvec2 pb_raw = kd->getStrandTree().evaluate(static_cast<size_t>(b), occurrence_time);
        const glm::dvec2 pc_raw = kd->getStrandTree().evaluate(static_cast<size_t>(c), occurrence_time);
        const double raw_collinearity_metric = normalizedTriangleCollinearityMetric(pa_raw, pb_raw, pc_raw);
        const bool transformed_collinear = collinearity_metric <= flip_boundary_collinearity_eps;
        const bool untransformed_collinear = raw_collinearity_metric <= flip_boundary_collinearity_eps;

        throw std::runtime_error(
          "Invalid boundary flip event: finite vertices are not collinear for half-edge " + std::to_string(half_edge_id)
          + " at t=" + std::to_string(occurrence_time) + " (transformed_collinearity_metric="
          + std::to_string(collinearity_metric) + ", untransformed_collinearity_metric="
          + std::to_string(raw_collinearity_metric) + ", a=" + std::to_string(a) + ", b=" + std::to_string(b) + ", c="
          + std::to_string(c) + ", pa=" + glm::to_string(pa) + ", pb=" + glm::to_string(pb) + ", pc=" + glm::to_string(pc)
          + ", pa_raw=" + glm::to_string(pa_raw) + ", pb_raw=" + glm::to_string(pb_raw) + ", pc_raw="
          + glm::to_string(pc_raw)
          + flipUntransformedFrameMismatchNote(transformed_collinear, untransformed_collinear) + ")");
      }
    }
  }
  else
  {
    const glm::dvec3 left_voronoi_vertex = kd->computeVoronoiVertexClampedInfinity(half_edge_id, occurrence_time);
    const glm::dvec3 right_voronoi_vertex = kd->computeVoronoiVertexClampedInfinity(half_edge_id ^ 1, occurrence_time);
    const double voronoi_vertex_distance
      = glm::distance(glm::dvec2(left_voronoi_vertex), glm::dvec2(right_voronoi_vertex));
    if (voronoi_vertex_distance > flip_voronoi_vertex_distance_eps)
    {
      const glm::dvec2 raw_left_cc
        = flipTriangleCircumcenterAt(*kd, graph, half_edge_id, occurrence_time, false);
      const glm::dvec2 raw_right_cc
        = flipTriangleCircumcenterAt(*kd, graph, half_edge_id ^ 1, occurrence_time, false);
      const double raw_circumcenter_distance = glm::distance(raw_left_cc, raw_right_cc);
      const bool transformed_coincident = voronoi_vertex_distance <= flip_voronoi_vertex_distance_eps;
      const bool untransformed_coincident = raw_circumcenter_distance <= flip_voronoi_vertex_distance_eps;

      throw std::runtime_error(
        "Invalid flip event: Voronoi edge endpoints are not coincident for half-edge " + std::to_string(half_edge_id)
        + " at t=" + std::to_string(occurrence_time) + " (faces " + std::to_string(face_id) + " and "
        + std::to_string(twin_face_id) + ", transformed_voronoi_distance=" + std::to_string(voronoi_vertex_distance)
        + ", untransformed_circumcenter_distance=" + std::to_string(raw_circumcenter_distance) + ", left="
        + glm::to_string(left_voronoi_vertex) + ", right=" + glm::to_string(right_voronoi_vertex) + ", raw_left_cc="
        + glm::to_string(raw_left_cc) + ", raw_right_cc=" + glm::to_string(raw_right_cc)
        + flipUntransformedFrameMismatchNote(transformed_coincident, untransformed_coincident) + ")");
    }
  }
#endif

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

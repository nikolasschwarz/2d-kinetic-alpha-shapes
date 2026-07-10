#include "KineticDelaunayCrossingEvent.hpp"

#include <algorithm>
#include <limits>

#include "KineticDelaunayEventPredicates.hpp"

using namespace kinDS;

void KineticDelaunay::CrossingEventManager::computeEvents(double t, size_t voronoi_vertex_id)
{
  auto* kd = kd_;
  auto& graph = kd->graph;
  auto& crossing_data = kd->crossing_data;
  auto& branch_trajs = kd->branch_trajs;

  if (!kd->on_the_fly_boundary)
  {
    return;
  }
  if (voronoi_vertex_id >= graph.faceSlotCount() || !graph.isLiveFace(voronoi_vertex_id)
    || !crossing_data.isVoronoiVertexRegistered(voronoi_vertex_id))
  {
    return;
  }
  const size_t section = static_cast<size_t>(t);
  const float fraction = t - section;

  auto& dual_triangle = graph.face(voronoi_vertex_id);
  auto& containing_triangle = graph.face(crossing_data.getContainingTriId(voronoi_vertex_id));

  std::vector<size_t> event_strand_ids;
  event_strand_ids.reserve(6);
  const auto append_strand = [&](int vertex) {
    if (vertex >= 0)
    {
      const size_t strand_id = static_cast<size_t>(vertex);
      if (std::find(event_strand_ids.begin(), event_strand_ids.end(), strand_id) == event_strand_ids.end())
      {
        event_strand_ids.push_back(strand_id);
      }
    }
  };
  for (size_t he_id : dual_triangle.half_edges)
  {
    append_strand(static_cast<int>(graph.halfEdge(he_id).origin));
    append_strand(static_cast<int>(graph.destination(he_id)));
  }
  for (size_t he_id : containing_triangle.half_edges)
  {
    append_strand(static_cast<int>(graph.halfEdge(he_id).origin));
    append_strand(static_cast<int>(graph.destination(he_id)));
  }

  const auto piece_poly = [&](size_t strand_id) {
    return kd->getSitePiecePolynomialForEventStrands(strand_id, section, t, event_strand_ids);
  };

  // compute polynomials of two bisectors in homogeneous coordinates
  size_t v_i = graph.halfEdge(dual_triangle.half_edges[0]).origin;
  size_t v_j = graph.halfEdge(dual_triangle.half_edges[1]).origin;
  size_t v_k = graph.halfEdge(dual_triangle.half_edges[2]).origin;

  // If a vertex is infinite, so is the Voronoi vertex and it cannot cross any edge, so we can skip this event.
  if (v_i == -1 || v_j == -1 || v_k == -1)
  {
    return;
  }

  // Check a special case: the containing triangle is infinite and adjacent to the dual triangle. In this case, we need
  // a different predicate
  bool adjacent = false;
  size_t adjacent_edge_index = -1;
  size_t finite_he_id = -1;
  if (graph.isInfinite(containing_triangle.half_edges[0]) || graph.isInfinite(containing_triangle.half_edges[1])
    || graph.isInfinite(containing_triangle.half_edges[2]))
  {
    // One edge must be finite, find it

    size_t finite_edge_index = -1;
    for (size_t edge_index = 0; edge_index < 3; edge_index++)
    {
      if (!graph.isInfinite(containing_triangle.half_edges[edge_index]))
      {
        finite_he_id = containing_triangle.half_edges[edge_index];
        finite_edge_index = edge_index;
        break;
      }
    }

    // Now iterate over the edges of the dual triangle and check if any of them is the twin of the finite edge.
    for (size_t edge_index = 0; edge_index < 3; edge_index++)
    {
      if (graph.twin(dual_triangle.half_edges[edge_index]) == finite_he_id)
      {
        adjacent = true;
        adjacent_edge_index = edge_index;
        break;
      }
    }
  }

  if (adjacent)
  {

    // re-assign vertices
    v_i = graph.triangleOppositeVertex(dual_triangle.half_edges[adjacent_edge_index]);
    v_j = graph.halfEdge(dual_triangle.half_edges[adjacent_edge_index]).origin;
    v_k = graph.halfEdge(dual_triangle.half_edges[adjacent_edge_index] ^ 1).origin;

    Trajectory<2> traj_i = piece_poly(v_i);
    Trajectory<2> traj_j = piece_poly(v_j);
    Trajectory<2> traj_k = piece_poly(v_k);

    Trajectory<2> vector_ij;
    vector_ij[0] = traj_j[0] - traj_i[0];
    vector_ij[1] = traj_j[1] - traj_i[1];
    Trajectory<2> vector_ik;
    vector_ik[0] = traj_k[0] - traj_i[0];
    vector_ik[1] = traj_k[1] - traj_i[1];

    Polynomial event_trigger = -(vector_ij[0] * vector_ik[0] + vector_ij[1] * vector_ik[1]);
    auto fractional_event_times = kd->findEvents(event_trigger, fraction, true);

    // Only need the first event as any following events will be invalidated by the first crossing event. TODO: The
    // exception is the edge being crossed, but that would make this more complex. We can optimize this later if needed.
    if (!fractional_event_times.empty())
    {
      double fractional_event_time = fractional_event_times.front();
      double event_time = fractional_event_time + section;
      // Position must be the midpoint of the two vertices
      glm::dvec2 position = glm::vec2((traj_j[0](fractional_event_time) + traj_k[0](fractional_event_time)) / 2.0,
        (traj_j[1](fractional_event_time) + traj_k[1](fractional_event_time)) / 2.0);

      KINDS_DEBUG("Crossing (right angle) Event queued at time "
        << event_time << " for Voronoi vertex ID " << voronoi_vertex_id << " crossing half-edge ID " << finite_he_id
        << " at position " << glm::to_string(position));
      kd->kinetic_algorithm_->enqueueEvent(
        std::make_shared<CrossingEvent>(kd, event_time, finite_he_id, t, position, voronoi_vertex_id));
    }
  }
  else
  {
    Trajectory<2> traj_i = piece_poly(v_i);
    Trajectory<2> traj_j = piece_poly(v_j);
    Trajectory<2> traj_k = piece_poly(v_k);
    Trajectory<3> bisector_ij;

    bisector_ij[0] = 2 * (traj_j[0] - traj_i[0]);
    bisector_ij[1] = 2 * (traj_j[1] - traj_i[1]);
    bisector_ij[2] = (traj_i[0] * traj_i[0] + traj_i[1] * traj_i[1]) - (traj_j[0] * traj_j[0] + traj_j[1] * traj_j[1]);

    Trajectory<3> bisector_ik;

    bisector_ik[0] = 2 * (traj_k[0] - traj_i[0]);
    bisector_ik[1] = 2 * (traj_k[1] - traj_i[1]);
    bisector_ik[2] = (traj_i[0] * traj_i[0] + traj_i[1] * traj_i[1]) - (traj_k[0] * traj_k[0] + traj_k[1] * traj_k[1]);

    // We only need the first event as any following events will be invalidated by the first crossing event.
    // TODO: The exception is the edge being crossed, but that would make this more complex. We can optimize this later
    // if needed.
    double event_time = std::numeric_limits<double>::infinity();
    size_t event_he_id = -1;
    Polynomial event_trigger;
    // Construct polynomial predicates for each of the three edges of the containing triangle
    for (size_t edge_index = 0; edge_index < 3; edge_index++)
    {
      size_t he_id = containing_triangle.half_edges[edge_index];
      size_t a = graph.halfEdge(he_id).origin;
      size_t b = graph.halfEdge(he_id ^ 1).origin;

      Trajectory<3> line_ab;
      if (a != -1 && b != -1)
      {

        Trajectory<2> traj_a = piece_poly(static_cast<size_t>(a));
        Trajectory<2> traj_b = piece_poly(static_cast<size_t>(b));

        // line through a and b in homogeneous coordinates

        line_ab[0] = traj_a[1] - traj_b[1];
        line_ab[1] = traj_b[0] - traj_a[0];
        line_ab[2] = traj_a[0] * traj_b[1] - traj_a[1] * traj_b[0];

        // now compute the determinant of the matrix with bisector_ij, bisector_ik and line_ab as columns
        event_trigger = bisector_ij[0] * bisector_ik[1] * line_ab[2] + bisector_ij[1] * bisector_ik[2] * line_ab[0]
          + bisector_ij[2] * bisector_ik[0] * line_ab[1] - bisector_ij[2] * bisector_ik[1] * line_ab[0]
          - bisector_ij[1] * bisector_ik[0] * line_ab[2] - bisector_ij[0] * bisector_ik[2] * line_ab[1];
      }
      else
      {
        size_t finite_vertex = (a != -1) ? a : b;

        if (a == finite_vertex)
        {
          size_t prev_he_id = graph.prev(he_id);
          size_t next_he_id = graph.halfEdge(he_id ^ 1).next;

          size_t c = graph.halfEdge(prev_he_id).origin;
          size_t c_prime = graph.halfEdge(next_he_id).origin;

          Trajectory<2> traj_a = piece_poly(static_cast<size_t>(a));
          Trajectory<2> traj_c = piece_poly(static_cast<size_t>(c));
          Trajectory<2> traj_c_prime = piece_poly(static_cast<size_t>(c_prime));
          Trajectory<3> voronoi_homogeneous = Trajectory<3>::cross(bisector_ij, bisector_ik);

          event_trigger = angularBisector(traj_a, traj_c, traj_c_prime, voronoi_homogeneous);
        }
        else
        {
          size_t prev_he_id = graph.prev(he_id ^ 1);
          size_t next_he_id = graph.halfEdge(he_id).next;

          size_t c_prime = graph.halfEdge(prev_he_id).origin;
          size_t c = graph.halfEdge(next_he_id).origin;

          Trajectory<2> traj_b = piece_poly(static_cast<size_t>(b));
          Trajectory<2> traj_c = piece_poly(static_cast<size_t>(c));
          Trajectory<2> traj_c_prime = piece_poly(static_cast<size_t>(c_prime));
          Trajectory<3> voronoi_homogeneous = Trajectory<3>::cross(bisector_ij, bisector_ik);

          event_trigger = angularBisector(traj_b, traj_c, traj_c_prime, voronoi_homogeneous);
        }
      }

      auto fractional_event_times = kd->findEvents(event_trigger, fraction, true);
      if (!fractional_event_times.empty())
      {
        double fractional_event_time = fractional_event_times.front();
        double candidate_event_time = fractional_event_time + section;
        if (candidate_event_time < event_time)
        {
          event_time = candidate_event_time;
          event_he_id = he_id;
        }
      }
    }

    if (event_time != std::numeric_limits<double>::infinity())
    {
      const double local_event_t = event_time - static_cast<double>(section);
      glm::dvec3 position_homogeneous;

      position_homogeneous[0] = bisector_ij[1](local_event_t) * bisector_ik[2](local_event_t)
        - bisector_ij[2](local_event_t) * bisector_ik[1](local_event_t);
      position_homogeneous[1] = bisector_ij[2](local_event_t) * bisector_ik[0](local_event_t)
        - bisector_ij[0](local_event_t) * bisector_ik[2](local_event_t);
      position_homogeneous[2] = bisector_ij[0](local_event_t) * bisector_ik[1](local_event_t)
        - bisector_ij[1](local_event_t) * bisector_ik[0](local_event_t);

      glm::dvec2 position(
        position_homogeneous.x / position_homogeneous.z, position_homogeneous.y / position_homogeneous.z);
      KINDS_DEBUG("Crossing Event scheduled at time " << event_time << " for Voronoi vertex ID " << voronoi_vertex_id
                                            << " crossing half-edge ID " << event_he_id << " at position "
                                            << glm::to_string(position));
      kd->kinetic_algorithm_->enqueueEvent(
        std::make_shared<CrossingEvent>(kd, event_time, event_he_id, t, position, voronoi_vertex_id));
    }
  }
}
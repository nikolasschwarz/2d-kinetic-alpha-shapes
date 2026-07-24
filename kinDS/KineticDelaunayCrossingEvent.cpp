#include "KineticDelaunayCrossingEvent.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
#include <optional>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "KineticDelaunayEventPredicates.hpp"

using namespace kinDS;

namespace
{
std::string formatSignLabel(double sign)
{
  return sign > 0.0 ? "+" : "-";
}

std::string formatSignChange(double sign_before, double sign_after)
{
  return "(" + formatSignLabel(sign_before) + " => " + formatSignLabel(sign_after) + ")";
}

std::pair<double, double> signChangeAtRoot(const Polynomial& event_trigger, double root,
  const std::vector<double>& sorted_zeros)
{
  const auto root_it = std::lower_bound(sorted_zeros.begin(), sorted_zeros.end(), root);
  const size_t root_index = static_cast<size_t>(root_it - sorted_zeros.begin());

  const double left_bound = root_index == 0 ? 0.0 : sorted_zeros[root_index - 1];
  const double right_bound
    = root_index + 1 < sorted_zeros.size() ? sorted_zeros[root_index + 1] : 1.0;

  const double before_test = (left_bound + root) * 0.5;
  const double after_test = (root + right_bound) * 0.5;
  const double sign_before = event_trigger(before_test) > 0.0 ? 1.0 : -1.0;
  const double sign_after = event_trigger(after_test) > 0.0 ? 1.0 : -1.0;
  return { sign_before, sign_after };
}

bool isTriangleAdjacentToMonitoredCrossingEdge(const KineticDelaunay& kd, size_t tri_id)
{
  if (!KineticDelaunay::isDiagnosticsMonitorIdEnabled(KineticDelaunay::kDiagnosticsMonitoredCrossingDelaunayEdgeId))
  {
    return false;
  }

  const size_t monitored_even_he_id = 2 * KineticDelaunay::kDiagnosticsMonitoredCrossingDelaunayEdgeId;
  const size_t monitored_odd_he_id = monitored_even_he_id + 1;

  for (size_t monitored_he_id : { monitored_even_he_id, monitored_odd_he_id })
  {
    if (!kd.isDiagnosticsHalfEdgeIdValid(monitored_he_id))
    {
      continue;
    }
    if (kd.getGraph().halfEdge(monitored_he_id).face == tri_id)
    {
      return true;
    }
  }
  return false;
}

bool shouldLogCrossingDiagnostics(const KineticDelaunay& kd, size_t voronoi_vertex_id, double schedule_t)
{
  if (!kd.diagnosticsEnabled()
    || !KineticDelaunay::matchesDiagnosticsMonitorId(
         voronoi_vertex_id, KineticDelaunay::kDiagnosticsMonitoredCrossingVoronoiVertexId)
    || !kd.isDiagnosticsFaceIdValid(voronoi_vertex_id)
    || schedule_t < std::floor(KineticDelaunay::kDiagnosticsMonitoredCrossingTime)
    || schedule_t >= std::floor(KineticDelaunay::kDiagnosticsMonitoredCrossingTime) + 1.0
    || !kd.getCrossingData().isVoronoiVertexRegistered(voronoi_vertex_id))
  {
    return false;
  }

  const size_t containing_tri_id = kd.getCrossingData().getContainingTriId(voronoi_vertex_id);
  return kd.diagnosticsEnabled()
    && kd.isDiagnosticsFaceIdValid(containing_tri_id)
    && isTriangleAdjacentToMonitoredCrossingEdge(kd, containing_tri_id);
}

/// Handle-side filter: match monitored vv/edge/time window without requiring registration, so skip reasons stay visible.
bool shouldLogCrossingHandleDiagnostics(
  const KineticDelaunay& kd, size_t voronoi_vertex_id, size_t half_edge_id, double occurrence_t)
{
  return kd.diagnosticsEnabled()
    && KineticDelaunay::matchesDiagnosticsMonitorId(
         voronoi_vertex_id, KineticDelaunay::kDiagnosticsMonitoredCrossingVoronoiVertexId)
    && KineticDelaunay::matchesDiagnosticsMonitorId(
         half_edge_id / 2, KineticDelaunay::kDiagnosticsMonitoredCrossingDelaunayEdgeId)
    && occurrence_t >= std::floor(KineticDelaunay::kDiagnosticsMonitoredCrossingTime)
    && occurrence_t < std::floor(KineticDelaunay::kDiagnosticsMonitoredCrossingTime) + 1.0;
}

void logCrossingComputeContext(const KineticDelaunay& kd, size_t voronoi_vertex_id, double schedule_t,
  size_t section, double min_fraction, size_t containing_tri_id, bool adjacent, size_t finite_he_id)
{
  const double event_interval_upper_bound = eventIntervalUpperBound(schedule_t);
  std::ostringstream header;
  header << "Crossing computeEvents monitor (voronoi_vertex=" << voronoi_vertex_id << ", schedule_t=" << schedule_t
         << ", section=" << section << ", min_fraction=" << min_fraction
         << ", event_interval_upper_bound=" << event_interval_upper_bound
         << ", containing_tri=" << containing_tri_id << ", adjacent_case=" << (adjacent ? "true" : "false");
  if (adjacent)
  {
    header << ", finite_he_id=" << finite_he_id << ", delaunay_edge=" << (finite_he_id / 2);
  }
  header << ", monitored_crossing_t=" << KineticDelaunay::kDiagnosticsMonitoredCrossingTime << " +/- "
         << KineticDelaunay::kDiagnosticsMonitoredCrossingTimeEpsilon
         << ", monitored_delaunay_edge=" << KineticDelaunay::kDiagnosticsMonitoredCrossingDelaunayEdgeId
         << ", monitored_half_edges=["
         << (2 * KineticDelaunay::kDiagnosticsMonitoredCrossingDelaunayEdgeId) << ","
         << (2 * KineticDelaunay::kDiagnosticsMonitoredCrossingDelaunayEdgeId + 1) << "])";
  KINDS_DEBUG(header.str());

  KINDS_DEBUG("  voronoi_dual_tri=" << voronoi_vertex_id << " containing_delaunay_tri=" << containing_tri_id
                                   << " registered="
                                   << (kd.getCrossingData().isVoronoiVertexRegistered(voronoi_vertex_id) ? "true"
                                                                                                       : "false"));
}

void logCrossingTriggerRoots(const KineticDelaunay& kd, size_t voronoi_vertex_id, size_t he_id, size_t edge_index,
  double schedule_t, double min_fraction, const Polynomial& event_trigger_in, bool only_positive_to_negative)
{
  Polynomial event_trigger = event_trigger_in;
  const size_t section = static_cast<size_t>(schedule_t);
  const size_t delaunay_edge_id = he_id / 2;
  const bool monitored_edge = KineticDelaunay::matchesDiagnosticsMonitorId(
    delaunay_edge_id, KineticDelaunay::kDiagnosticsMonitoredCrossingDelaunayEdgeId);
  std::ostringstream header;
  header << "  crossing trigger roots";
  if (monitored_edge)
  {
    header << " **MONITORED_EDGE**";
  }
  header << " (voronoi_vertex=" << voronoi_vertex_id << ", edge_index=" << edge_index << ", he_id=" << he_id
         << ", delaunay_edge=" << delaunay_edge_id << ", schedule_t=" << schedule_t << ", section=" << section
         << ", min_fraction=" << min_fraction << ", only_positive_to_negative=" << (only_positive_to_negative ? "true"
                                                                                                              : "false")
         << ", trigger_degree=" << event_trigger.degree() << ")";
  KINDS_DEBUG(header.str());

  if (event_trigger.degree() == -1)
  {
    KINDS_DEBUG("    trigger empty (degree -1)");
    return;
  }

  event_trigger.trim();
  const auto zeros = event_trigger.realRoots();
  if (zeros.empty())
  {
    KINDS_DEBUG("    no real roots");
    return;
  }

  std::vector<double> sorted_zeros;
  sorted_zeros.reserve(zeros.size());
  for (double root : zeros)
  {
    if (!std::isnan(root))
    {
      sorted_zeros.push_back(root);
    }
  }
  std::sort(sorted_zeros.begin(), sorted_zeros.end());

  std::optional<double> first_enqueued_fraction;
  for (size_t root_index = 0; root_index < zeros.size(); ++root_index)
  {
    const double root = zeros[root_index];
    std::ostringstream line;
    line << std::setprecision(17) << "    root[" << root_index << "] fraction=" << root << " absolute_t="
         << (root + static_cast<double>(section));

    if (std::isnan(root))
    {
      line << " discarded (nan)";
      KINDS_DEBUG(line.str());
      continue;
    }

    const auto [sign_before, sign_after] = signChangeAtRoot(event_trigger, root, sorted_zeros);
    line << " sign_change=" << formatSignChange(sign_before, sign_after);

    if (root <= min_fraction)
    {
      line << " discarded (fraction <= min_fraction)";
      KINDS_DEBUG(line.str());
      continue;
    }
    if (root > kEventIntervalFractionUpperBound)
    {
      line << " discarded (fraction > " << kEventIntervalFractionUpperBound << ")";
      KINDS_DEBUG(line.str());
      continue;
    }

    if (sign_before == sign_after)
    {
      line << " discarded (no_sign_change)";
      KINDS_DEBUG(line.str());
      continue;
    }

    if (only_positive_to_negative && sign_before < 0.0)
    {
      line << " discarded (negative_to_positive, only_positive_to_negative=true)";
      KINDS_DEBUG(line.str());
      continue;
    }

    const bool positive_to_negative = sign_before > 0.0 && sign_after < 0.0;
    line << " queued_by_findEvents positive_to_negative=" << (positive_to_negative ? "true" : "false");
    if (!first_enqueued_fraction.has_value())
    {
      first_enqueued_fraction = root;
      line << " **first_enqueued_root**";
    }
    KINDS_DEBUG(line.str());
  }

  if (!first_enqueued_fraction.has_value())
  {
    KINDS_DEBUG("    findEvents would return empty for this trigger");
  }
  else
  {
    KINDS_DEBUG("    findEvents would return first root fraction=" << *first_enqueued_fraction << " absolute_t="
                                                                  << (*first_enqueued_fraction
                                                                    + static_cast<double>(section)));
  }
}

struct CrossingEdgeCandidateSummary
{
  size_t he_id = size_t(-1);
  size_t delaunay_edge_id = size_t(-1);
  bool monitored_edge = false;
  bool had_enqueued_root = false;
  double first_enqueued_fraction = std::numeric_limits<double>::infinity();
  double first_enqueued_absolute_t = std::numeric_limits<double>::infinity();
};

void logCrossingCandidateSelection(const KineticDelaunay& kd, size_t voronoi_vertex_id, double schedule_t,
  const std::vector<CrossingEdgeCandidateSummary>& candidates, size_t selected_he_id, double selected_event_time)
{
  KINDS_DEBUG("  crossing candidate selection summary (voronoi_vertex=" << voronoi_vertex_id
                                                                         << ", schedule_t=" << schedule_t << ")");
  for (size_t edge_index = 0; edge_index < candidates.size(); ++edge_index)
  {
    const auto& candidate = candidates[edge_index];
    std::ostringstream line;
    line << "    edge_index=" << edge_index << " he_id=" << candidate.he_id
         << " delaunay_edge=" << candidate.delaunay_edge_id;
    if (candidate.monitored_edge)
    {
      line << " **MONITORED_EDGE**";
    }
    if (!candidate.had_enqueued_root)
    {
      line << " discarded (no queued root from findEvents)";
    }
    else
    {
      line << std::setprecision(17) << " first_queued_fraction=" << candidate.first_enqueued_fraction
           << " first_queued_absolute_t=" << candidate.first_enqueued_absolute_t;
      if (selected_he_id != size_t(-1) && candidate.he_id == selected_he_id)
      {
        line << " **SELECTED**";
      }
      else if (candidate.first_enqueued_absolute_t > selected_event_time)
      {
        line << " discarded (later than selected candidate at t=" << selected_event_time << ")";
      }
      else if (selected_he_id == size_t(-1))
      {
        line << " discarded (no overall winner)";
      }
      else
      {
        line << " discarded (not earliest among candidates)";
      }
    }
    KINDS_DEBUG(line.str());
  }

  if (selected_he_id == size_t(-1))
  {
    KINDS_DEBUG("  crossing event NOT queued (no candidate roots)");
  }
  else
  {
    KINDS_DEBUG("  crossing event queued he_id=" << selected_he_id << " delaunay_edge=" << (selected_he_id / 2)
                                                << " event_time=" << std::setprecision(17) << selected_event_time);
  }
}
} // namespace

bool KineticDelaunay::isDiagnosticsMonitoredCrossingValid() const
{
  return isDiagnosticsMonitorIdEnabled(kDiagnosticsMonitoredCrossingVoronoiVertexId)
    && isDiagnosticsFaceIdValid(kDiagnosticsMonitoredCrossingVoronoiVertexId);
}

void KineticDelaunay::logCrossingEventTriggerRoots(size_t voronoi_vertex_id, size_t he_id, size_t edge_index,
  double t, double min_fraction, const Polynomial& event_trigger, bool only_positive_to_negative) const
{
  logCrossingTriggerRoots(
    *this, voronoi_vertex_id, he_id, edge_index, t, min_fraction, event_trigger, only_positive_to_negative);
}

void KineticDelaunay::CrossingEventManager::computeEvents(double t, size_t voronoi_vertex_id)
{
  auto* kd = kd_;
  auto& graph = kd->graph;
  auto& crossing_data = kd->crossing_data;
  auto& branch_trajs = kd->branch_trajs;

  const bool log_crossing_diag = shouldLogCrossingDiagnostics(*kd, voronoi_vertex_id, t);
  const auto log_early_exit = [&](const char* reason)
  {
    if (log_crossing_diag)
    {
      KINDS_DEBUG("Crossing computeEvents early exit (voronoi_vertex=" << voronoi_vertex_id << ", schedule_t=" << t
                                                                       << "): " << reason);
    }
  };

  if (!kd->on_the_fly_boundary)
  {
    log_early_exit("on_the_fly_boundary=false");
    return;
  }
  if (voronoi_vertex_id >= graph.faceSlotCount() || !graph.isLiveFace(voronoi_vertex_id)
    || !crossing_data.isVoronoiVertexRegistered(voronoi_vertex_id))
  {
    log_early_exit("voronoi vertex not live or not registered");
    return;
  }
  const size_t section = static_cast<size_t>(t);
  const float fraction = t - section;

  auto& dual_triangle = graph.face(voronoi_vertex_id);
  const size_t containing_tri_id = crossing_data.getContainingTriId(voronoi_vertex_id);
  auto& containing_triangle = graph.face(containing_tri_id);

  std::vector<size_t> event_strand_ids;
  event_strand_ids.reserve(6);
  const auto append_strand = [&](int vertex)
  {
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

  const auto piece_poly = [&](size_t strand_id)
  {
    return kd->getSitePiecePolynomialForEventStrands(strand_id, section, t, event_strand_ids);
  };

  // compute polynomials of two bisectors in homogeneous coordinates
  size_t v_i = graph.halfEdge(dual_triangle.half_edges[0]).origin;
  size_t v_j = graph.halfEdge(dual_triangle.half_edges[1]).origin;
  size_t v_k = graph.halfEdge(dual_triangle.half_edges[2]).origin;

  // If a vertex is infinite, so is the Voronoi vertex and it cannot cross any edge, so we can skip this event.
  if (v_i == -1 || v_j == -1 || v_k == -1)
  {
    log_early_exit("infinite Voronoi vertex (dual triangle has infinite site)");
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

  if (log_crossing_diag)
  {
    logCrossingComputeContext(*kd, voronoi_vertex_id, t, section, fraction, containing_tri_id, adjacent, finite_he_id);
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
    if (log_crossing_diag)
    {
      kd->logCrossingEventTriggerRoots(
        voronoi_vertex_id, finite_he_id, adjacent_edge_index, t, fraction, event_trigger, true);
    }
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

      if (log_crossing_diag)
      {
        const CrossingEdgeCandidateSummary candidate { finite_he_id, finite_he_id / 2,
          KineticDelaunay::matchesDiagnosticsMonitorId(
            finite_he_id / 2, KineticDelaunay::kDiagnosticsMonitoredCrossingDelaunayEdgeId),
          true, fractional_event_time, event_time };
        logCrossingCandidateSelection(*kd, voronoi_vertex_id, t, { candidate }, finite_he_id, event_time);
      }

      KINDS_DEBUG("Crossing (right angle) Event queued at time "
        << event_time << " for Voronoi vertex ID " << voronoi_vertex_id << " crossing half-edge ID " << finite_he_id
        << " at position " << glm::to_string(position));
      kd->kinetic_algorithm_->enqueueEvent(
        std::make_shared<CrossingEvent>(kd, event_time, finite_he_id, t, position, voronoi_vertex_id));
    }
    else if (log_crossing_diag)
    {
      const CrossingEdgeCandidateSummary candidate { finite_he_id, finite_he_id / 2,
        KineticDelaunay::matchesDiagnosticsMonitorId(
          finite_he_id / 2, KineticDelaunay::kDiagnosticsMonitoredCrossingDelaunayEdgeId),
        false, std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity() };
      logCrossingCandidateSelection(*kd, voronoi_vertex_id, t, { candidate }, size_t(-1),
        std::numeric_limits<double>::infinity());
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
    std::vector<CrossingEdgeCandidateSummary> candidate_summaries;
    candidate_summaries.reserve(3);
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

      if (log_crossing_diag)
      {
        kd->logCrossingEventTriggerRoots(voronoi_vertex_id, he_id, edge_index, t, fraction, event_trigger, true);
      }

      auto fractional_event_times = kd->findEvents(event_trigger, fraction, true);
      CrossingEdgeCandidateSummary candidate_summary;
      candidate_summary.he_id = he_id;
      candidate_summary.delaunay_edge_id = he_id / 2;
      candidate_summary.monitored_edge = KineticDelaunay::matchesDiagnosticsMonitorId(
        candidate_summary.delaunay_edge_id, KineticDelaunay::kDiagnosticsMonitoredCrossingDelaunayEdgeId);
      if (!fractional_event_times.empty())
      {
        candidate_summary.had_enqueued_root = true;
        candidate_summary.first_enqueued_fraction = fractional_event_times.front();
        candidate_summary.first_enqueued_absolute_t = fractional_event_times.front() + section;
      }
      candidate_summaries.push_back(candidate_summary);

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

    if (log_crossing_diag)
    {
      logCrossingCandidateSelection(
        *kd, voronoi_vertex_id, t, candidate_summaries, event_he_id, event_time);
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

void KineticDelaunay::CrossingEvent::handleEvent()
{
  auto* kd = kd_;
  if (!kd)
  {
    throw std::runtime_error("CrossingEvent has no KineticDelaunay pointer");
  }

  auto& graph = kd->graph;
  const bool log_crossing_diag
    = shouldLogCrossingHandleDiagnostics(*kd, voronoi_vertex_id, half_edge_id, occurrence_time);
  const auto log_skip = [&](const char* reason)
  {
    if (log_crossing_diag)
    {
      KINDS_DEBUG("Crossing handleEvent SKIP (voronoi_vertex=" << voronoi_vertex_id << ", he_id=" << half_edge_id
                                                              << ", delaunay_edge=" << (half_edge_id / 2)
                                                              << ", occurrence_t=" << std::setprecision(17)
                                                              << occurrence_time << ", creation_t=" << creation_time
                                                              << "): " << reason);
    }
  };

  if (log_crossing_diag)
  {
    const bool registered = kd->crossing_data.isVoronoiVertexRegistered(voronoi_vertex_id);
    const std::optional<size_t> containing_tri
      = registered ? std::optional<size_t>(kd->crossing_data.getContainingTriId(voronoi_vertex_id)) : std::nullopt;
    const double last_crossing
      = (voronoi_vertex_id < kd->crossing_data.last_crossing.size()) ? kd->crossing_data.last_crossing[voronoi_vertex_id]
                                                                     : std::numeric_limits<double>::quiet_NaN();
    const double containing_last_updated
      = (containing_tri.has_value() && *containing_tri < kd->face_last_updated.size())
      ? kd->face_last_updated[*containing_tri]
      : std::numeric_limits<double>::quiet_NaN();

    KINDS_DEBUG("Crossing handleEvent ENTER (voronoi_vertex=" << voronoi_vertex_id << ", he_id=" << half_edge_id
                                                             << ", delaunay_edge=" << (half_edge_id / 2)
                                                             << ", occurrence_t=" << std::setprecision(17)
                                                             << occurrence_time << ", creation_t=" << creation_time
                                                             << ", registered=" << (registered ? "true" : "false")
                                                             << ", containing_tri="
                                                             << (containing_tri.has_value()
                                                                   ? std::to_string(*containing_tri)
                                                                   : std::string("n/a"))
                                                             << ", last_crossing=" << last_crossing
                                                             << ", containing_face_last_updated="
                                                             << containing_last_updated
                                                             << ", he_live="
                                                             << (graph.isLiveHalfEdge(half_edge_id) ? "true" : "false")
                                                             << ")");
  }

  // Check if the event is still valid
  if (creation_time < kd->crossing_data.last_crossing[voronoi_vertex_id])
  {
    log_skip("creation_time < last_crossing (stale vs later crossing/recompute stamp)");
    return;
  }

  // Outdated if the crossed Delaunay edge or the Delaunay triangle dual to this Voronoi vertex
  // was tombstoned (e.g. after a branch split).
  if (!graph.isLiveHalfEdge(half_edge_id))
  {
    log_skip("half_edge not live");
    return;
  }
  if (voronoi_vertex_id >= graph.faceSlotCount() || !graph.isLiveFace(voronoi_vertex_id)
    || !kd->crossing_data.isVoronoiVertexRegistered(voronoi_vertex_id))
  {
    log_skip("voronoi vertex not live or not registered");
    return;
  }

  size_t containing_tri_id = kd->crossing_data.getContainingTriId(voronoi_vertex_id);

  // Outdated if the containing Delaunay triangle was updated after this event was scheduled.
  if (creation_time < kd->face_last_updated[containing_tri_id])
  {
    log_skip("creation_time < face_last_updated[containing_tri]");
    return;
  }

  if (log_crossing_diag)
  {
    KINDS_DEBUG("Crossing handleEvent PROCEED (voronoi_vertex=" << voronoi_vertex_id << ", he_id=" << half_edge_id
                                                               << ", containing_tri=" << containing_tri_id
                                                               << ", occurrence_t=" << std::setprecision(17)
                                                               << occurrence_time << ")");
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
                                       << graph.halfEdge(half_edge_id ^ 1).face);
  kd->crossing_data.moveVertex(voronoi_vertex_id, graph.halfEdge(half_edge_id ^ 1).face, occurrence_time);

  // Update Voronoi–Delaunay edge intersections stored in crossing_data in response to this crossing.
  kd->crossing_data.updateAfterCrossingEvent(*kd, *this);

  // Recompute params on all three Delaunay edges of this Voronoi vertex (Delaunay face) before callbacks mesh.
  kd->refreshTriangleDelaunayEdgeIntersectionParams(voronoi_vertex_id, occurrence_time);

  if (event_handler)
  {
    event_handler->afterEvent(*this);
  }

  // After callbacks (e.g. debug SVG export); intersection lists must be consistent.
  kd->validateVoronoiVertexIteratorInvariants("CrossingEvent:afterEvent", occurrence_time);
  kd->validateCrossingIntersectionInvariants("CrossingEvent:afterEvent", occurrence_time);

  // Re-compute crossing events for this Voronoi vertex
  kd->crossing_event_manager_->computeEvents(occurrence_time, voronoi_vertex_id);

  if (kd->diagnosticsEnabled() && kd->isDiagnosticsMonitoredFaceValid())
  {
    kd->logFaceInsideStateAtTime(kDiagnosticsMonitoredFaceId, occurrence_time, "crossing_event");
  }
}

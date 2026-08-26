#include "KineticDelaunayFlipEvent.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
#include <optional>
#include <sstream>
#include <utility>
#include <vector>

#include "KineticDelaunayEventPredicates.hpp"
#include "Logger.hpp"

using namespace kinDS;

namespace
{
double referenceBranchLookupTimeForSection(size_t section, double schedule_time)
{
  const double section_start = static_cast<double>(section);
  double lookup_time = schedule_time;
  if (lookup_time <= section_start + std::numeric_limits<double>::epsilon())
  {
    lookup_time = section_start + std::numeric_limits<double>::epsilon();
  }
  return lookup_time;
}

std::string formatSignLabel(double sign)
{
  return sign > 0.0 ? "+" : "-";
}

std::string formatSignChange(double sign_before, double sign_after)
{
  return "(" + formatSignLabel(sign_before) + " => " + formatSignLabel(sign_after) + ")";
}

std::pair<double, double> signChangeAtRoot(const Polynomial& event_trigger, double root,
  const std::vector<double>& sorted_zeros, double left_default = 0.0, bool virtual_unbounded = false)
{
  const auto root_it = std::lower_bound(sorted_zeros.begin(), sorted_zeros.end(), root);
  const size_t root_index = static_cast<size_t>(root_it - sorted_zeros.begin());

  const double left_bound = root_index == 0 ? left_default : sorted_zeros[root_index - 1];
  const double right_bound = root_index + 1 < sorted_zeros.size()
    ? sorted_zeros[root_index + 1]
    : (virtual_unbounded ? root + 1.0 : 1.0);

  const double before_test = (left_bound + root) * 0.5;
  const double after_test = (root + right_bound) * 0.5;
  const double sign_before = event_trigger(before_test) > 0.0 ? 1.0 : -1.0;
  const double sign_after = event_trigger(after_test) > 0.0 ? 1.0 : -1.0;
  return { sign_before, sign_after };
}

bool shouldLogFlipDiagnostics(const KineticDelaunay& kd, size_t he_id, double schedule_t)
{
  return kd.diagnosticsEnabled()
    && KineticDelaunay::matchesDiagnosticsMonitorId(he_id / 2, KineticDelaunay::kDiagnosticsMonitoredFlipDelaunayEdgeId)
    && schedule_t >= std::floor(KineticDelaunay::kDiagnosticsMonitoredFlipTime)
    && schedule_t < std::floor(KineticDelaunay::kDiagnosticsMonitoredFlipTime) + 1.0;
}

void logFlipTriggerRoots(const KineticDelaunay& kd, size_t he_id, double schedule_t, double min_fraction,
  const Polynomial& event_trigger_in, const char* trigger_pass, const char* trigger_predicate, bool virtual_mode,
  double frozen_real_t)
{
  Polynomial event_trigger = event_trigger_in;
  const size_t section = static_cast<size_t>(schedule_t);
  const size_t delaunay_edge_id = he_id / 2;
  const double root_min = min_fraction;

  const EventTime schedule_event_time(schedule_t, virtual_mode ? root_min : 0.0);
  std::ostringstream header;
  header << "  flip trigger roots **MONITORED_EDGE** (he_id=" << he_id << ", delaunay_edge=" << delaunay_edge_id
         << ", schedule_t=" << schedule_event_time << ", section=" << section
         << ", " << (virtual_mode ? "min_infinitesimal_t=" : "min_fraction=") << root_min
         << ", pass=" << trigger_pass << ", predicate=" << trigger_predicate
         << ", trigger_degree=" << event_trigger.degree() << ")";
  if (virtual_mode)
  {
    header << ", virtual=true, frozen_real_t=" << frozen_real_t;
  }
  KINDS_MONITOR(header.str());

  if (event_trigger.degree() == -1)
  {
    KINDS_MONITOR("    trigger empty (degree -1) predicate=" << trigger_predicate);
    return;
  }

  event_trigger.trim();
  const auto zeros = event_trigger.realRoots();
  if (zeros.empty())
  {
    KINDS_MONITOR("    no real roots");
    return;
  }

  std::vector<double> sorted_zeros;
  sorted_zeros.reserve(zeros.size());
  for (double root : zeros)
  {
    if (!std::isnan(root) && std::isfinite(root))
    {
      sorted_zeros.push_back(root);
    }
  }
  std::sort(sorted_zeros.begin(), sorted_zeros.end());

  size_t queued_count = 0;
  for (size_t root_index = 0; root_index < zeros.size(); ++root_index)
  {
    const double root = zeros[root_index];
    std::ostringstream line;
    line << std::setprecision(17) << "    root[" << root_index << "] ";
    if (virtual_mode)
    {
      // Polynomial parameter is virtual / infinitesimal time, not a section fraction.
      line << "infinitesimal_t=" << root << " frozen_real_t=" << frozen_real_t;
    }
    else
    {
      line << "fraction=" << root << " absolute_t=" << (root + static_cast<double>(section));
    }

    if (std::isnan(root) || !std::isfinite(root))
    {
      line << " discarded (nan/non-finite)";
      KINDS_MONITOR(line.str());
      continue;
    }

    const auto [sign_before, sign_after]
      = signChangeAtRoot(event_trigger, root, sorted_zeros, root_min, virtual_mode);
    line << " sign_change=" << formatSignChange(sign_before, sign_after);

    if (root <= root_min)
    {
      line << " discarded (" << (virtual_mode ? "infinitesimal_t" : "fraction") << " <= min)";
      KINDS_MONITOR(line.str());
      continue;
    }
    if (!virtual_mode && root > kEventIntervalFractionUpperBound)
    {
      line << " discarded (fraction > " << kEventIntervalFractionUpperBound << ")";
      KINDS_MONITOR(line.str());
      continue;
    }

    if (sign_before == sign_after)
    {
      line << " discarded (no_sign_change)";
      KINDS_MONITOR(line.str());
      continue;
    }

    line << (virtual_mode ? " queued_by_findVirtualEvents" : " queued_by_findEvents");
    ++queued_count;
    KINDS_MONITOR(line.str());
  }

  if (queued_count == 0)
  {
    KINDS_MONITOR("    " << (virtual_mode ? "findVirtualEvents" : "findEvents")
                         << " would return empty for this trigger");
  }
  else
  {
    KINDS_MONITOR("    " << (virtual_mode ? "findVirtualEvents" : "findEvents") << " would queue "
                         << queued_count << " root(s) for this trigger");
  }
  (void)kd;
}
} // namespace

void KineticDelaunay::logFlipEventTriggerRoots(size_t he_id, double t, double min_fraction,
  const Polynomial& event_trigger, const char* trigger_pass, const char* trigger_predicate) const
{
  if (!shouldLogFlipDiagnostics(*this, he_id, t))
  {
    return;
  }
  logFlipTriggerRoots(*this, he_id, t, min_fraction, event_trigger, trigger_pass, trigger_predicate,
    computing_infinitesimal_events_, infinitesimal_schedule_t_);
}

void KineticDelaunay::FlipEventManager::computeEvents(double t, size_t quad_id,
  std::optional<InfinitesimalComputeContext> infinitesimal)
{
  auto* kd = kd_;
  auto& graph = kd->graph;

  std::optional<KineticDelaunay::ScopedInfinitesimalEventCompute> scope;
  if (infinitesimal.has_value())
  {
    scope.emplace(*kd, infinitesimal->parent_component_id, t, infinitesimal->min_infinitesimal_t);
    if (!scope->active())
    {
      return;
    }
  }

  const size_t section = static_cast<size_t>(t);
  const float fraction = t - section;

  size_t he_id = quad_id * 2;
  const bool log_flip_diag = shouldLogFlipDiagnostics(*kd, he_id, t);
  const std::vector<size_t> quad_strand_ids = collectFlipQuadrilateralStrandIds(graph, he_id);

  const bool virtual_mode = infinitesimal.has_value();
  const double root_min = virtual_mode ? kd->infinitesimal_recompute_min_x_ : static_cast<double>(fraction);
  if (log_flip_diag)
  {
    const EventTime schedule_event_time(t, virtual_mode ? root_min : 0.0);
    std::ostringstream header;
    header << "Flip computeEvents monitor (he_id=" << he_id << "/" << (he_id ^ 1) << ", delaunay_edge=" << quad_id
           << ", schedule_t=" << schedule_event_time << ", section=" << section
           << ", " << (virtual_mode ? "min_infinitesimal_t=" : "min_fraction=") << root_min
           << ", event_interval_upper_bound=" << eventIntervalUpperBound(t)
           << ", pass=" << (virtual_mode ? "infinitesimal" : "primary")
           << ", monitored_flip_t=" << KineticDelaunay::kDiagnosticsMonitoredFlipTime
           << ", he_live=" << (graph.isLiveHalfEdge(he_id) ? "true" : "false")
           << ", on_convex_boundary=" << (graph.isOnConvexBoundary(he_id) ? "true" : "false")
           << ", outside_convex_boundary=" << (graph.isOutsideConvexBoundary(he_id) ? "true" : "false")
           << ", quad_strands=[";
    for (size_t i = 0; i < quad_strand_ids.size(); ++i)
    {
      if (i > 0)
      {
        header << ',';
      }
      header << quad_strand_ids[i];
    }
    header << "], quadrilateral_last_updated="
           << (quad_id < kd->quadrilateral_last_updated.size() ? kd->quadrilateral_last_updated[quad_id]
                                                               : EventTime(std::numeric_limits<double>::quiet_NaN()))
           << ")";
    KINDS_MONITOR(header.str());
  }

  const auto build_trigger = [&](size_t active_he_id, double schedule_time, Polynomial& event_trigger_out,
                                 std::vector<Trajectory<2>>& trajs_out, const char*& predicate_out) {
    trajs_out.clear();
    std::vector<size_t> trigger_strand_ids;
    const auto append_trigger_strand = [&](int vertex) {
      if (vertex < 0)
      {
        return;
      }
      const size_t strand_id = static_cast<size_t>(vertex);
      if (std::find(trigger_strand_ids.begin(), trigger_strand_ids.end(), strand_id) == trigger_strand_ids.end())
      {
        trigger_strand_ids.push_back(strand_id);
      }
    };
    const auto piece_for_trigger = [&](size_t strand_id, double schedule_time) {
      return kd->getSitePiecePolynomialForEventStrands(strand_id, section, schedule_time, trigger_strand_ids);
    };

    if (graph.isOnConvexBoundary(active_he_id) || graph.isOutsideConvexBoundary(active_he_id))
    {
      predicate_out = "ccw";
      if (graph.isOutsideConvexBoundary(active_he_id))
      {
        active_he_id = active_he_id ^ 1;
      }

      int indices[4];
      indices[0] = graph.halfEdge(active_he_id).origin;
      indices[1] = graph.triangleOppositeVertex(active_he_id ^ 1);
      indices[2] = graph.halfEdge(active_he_id ^ 1).origin;
      indices[3] = graph.triangleOppositeVertex(active_he_id);

      std::vector<int> filtered_indices;
      std::copy_if(indices, indices + 4, std::back_inserter(filtered_indices), [](int index) { return index != -1; });
      if (filtered_indices.size() < 3)
      {
        event_trigger_out = Polynomial();
        return;
      }

      for (int vertex : filtered_indices)
      {
        append_trigger_strand(vertex);
      }

      trajs_out.push_back(piece_for_trigger(static_cast<size_t>(filtered_indices[0]), schedule_time));
      trajs_out.push_back(piece_for_trigger(static_cast<size_t>(filtered_indices[1]), schedule_time));
      trajs_out.push_back(piece_for_trigger(static_cast<size_t>(filtered_indices[2]), schedule_time));
      event_trigger_out = ccw(trajs_out[0][0], trajs_out[0][1], trajs_out[1][0], trajs_out[1][1], trajs_out[2][0],
        trajs_out[2][1]);
      return;
    }

    predicate_out = "inCircle";
    const int a = graph.halfEdge(active_he_id).origin;
    const int b = graph.triangleOppositeVertex(active_he_id ^ 1);
    const int c = graph.halfEdge(active_he_id ^ 1).origin;
    const int d = graph.triangleOppositeVertex(active_he_id);
    for (int vertex : { a, b, c, d })
    {
      append_trigger_strand(vertex);
    }
    trajs_out.push_back(piece_for_trigger(static_cast<size_t>(a), schedule_time));
    trajs_out.push_back(piece_for_trigger(static_cast<size_t>(b), schedule_time));
    trajs_out.push_back(piece_for_trigger(static_cast<size_t>(c), schedule_time));
    trajs_out.push_back(piece_for_trigger(static_cast<size_t>(d), schedule_time));
    event_trigger_out = inCircle(trajs_out[0][0], trajs_out[0][1], trajs_out[1][0], trajs_out[1][1], trajs_out[2][0],
      trajs_out[2][1], trajs_out[3][0], trajs_out[3][1]);
  };

  const auto enqueue_flip_roots = [&](const std::vector<Trajectory<2>>& trajs, Polynomial& event_trigger,
                                    double min_fraction, size_t enqueue_he_id, double creation_time,
                                    const char* trigger_pass, const char* trigger_predicate) {
    const bool virtual_mode = kd->computing_infinitesimal_events_;
    const double root_min = virtual_mode ? kd->infinitesimal_recompute_min_x_ : min_fraction;
    if (log_flip_diag)
    {
      kd->logFlipEventTriggerRoots(enqueue_he_id, t, root_min, event_trigger, trigger_pass, trigger_predicate);
    }

    if (event_trigger.degree() < 0)
    {
      if (log_flip_diag)
      {
        KINDS_MONITOR("  flip enqueue skipped (empty trigger) pass=" << trigger_pass
                                                                     << " predicate=" << trigger_predicate);
      }
      return;
    }

    auto event_times = virtual_mode ? kd->findVirtualEvents(event_trigger, root_min)
                                    : kd->findEvents(event_trigger, root_min);
    if (log_flip_diag && event_times.empty())
    {
      KINDS_MONITOR("  flip enqueue: " << (virtual_mode ? "findVirtualEvents" : "findEvents")
                                       << " returned empty pass=" << trigger_pass
                                       << " predicate=" << trigger_predicate);
    }
    for (const auto& event_time : event_times)
    {
      glm::dvec2 center {};
      for (const auto& traj : trajs)
      {
        center[0] += traj[0](event_time);
        center[1] += traj[1](event_time);
      }
      center[0] /= static_cast<double>(trajs.size());
      center[1] /= static_cast<double>(trajs.size());

      const EventTime occurrence = virtual_mode ? EventTime(kd->infinitesimal_schedule_t_, event_time)
                                                : EventTime(event_time + section);
      const EventTime creation = virtual_mode
        ? EventTime(kd->infinitesimal_schedule_t_, kd->infinitesimal_recompute_min_x_)
        : EventTime(creation_time);

      if (log_flip_diag)
      {
        if (virtual_mode)
        {
          KINDS_MONITOR("  flip event QUEUED pass=" << trigger_pass << " predicate=" << trigger_predicate
                                                    << " frozen_real_t=" << std::setprecision(17)
                                                    << occurrence.real_time
                                                    << " infinitesimal_t=" << occurrence.infinitesimal_time
                                                    << " creation_inf=" << creation.infinitesimal_time
                                                    << " he_id=" << enqueue_he_id << " center=(" << center[0]
                                                    << "," << center[1] << ")");
        }
        else
        {
          KINDS_MONITOR("  flip event QUEUED absolute_t=" << std::setprecision(17) << occurrence.real_time
                                                         << " fraction=" << event_time << " he_id=" << enqueue_he_id
                                                         << " creation_t=" << creation.real_time << " pass="
                                                         << trigger_pass << " predicate=" << trigger_predicate
                                                         << " center=(" << center[0] << "," << center[1] << ")");
        }
      }

      KINDS_DEBUG("Scheduled flip event at time " << occurrence.real_time << " infinitesimal_t="
                                                << occurrence.infinitesimal_time << " for half-edge ID "
                                                << enqueue_he_id << " at center position "
                                                << glm::to_string(center));

      auto flip_event = std::make_shared<FlipEvent>(kd, occurrence.real_time, enqueue_he_id, creation.real_time, center);
      flip_event->occurrence_time = occurrence;
      flip_event->creation_time = creation;
      if (virtual_mode)
      {
        flip_event->infinitesimal_epoch_ = kd->infinitesimal_schedule_epoch_;
      }
      kd->kinetic_algorithm_->enqueueEvent(flip_event);
    }
  };

  Polynomial event_trigger;
  std::vector<Trajectory<2>> trajs;
  const char* trigger_predicate = "inCircle";
  build_trigger(he_id, t, event_trigger, trajs, trigger_predicate);
  const char* trigger_pass = kd->computing_infinitesimal_events_ ? "infinitesimal" : "primary";
  enqueue_flip_roots(trajs, event_trigger, fraction, he_id, t, trigger_pass, trigger_predicate);
}

void KineticDelaunay::FlipEvent::handleEvent()
{
  auto* kd = kd_;
  if (!kd)
  {
    throw std::runtime_error("FlipEvent has no KineticDelaunay pointer");
  }

  auto& graph = kd->graph;
  const double t = occurrence_time.real_time;
  const double infinitesimal_t = occurrence_time.infinitesimal_time;
  const bool is_infinitesimal = infinitesimal_t > 0.0;
  const bool log_flip_diag = shouldLogFlipDiagnostics(*kd, half_edge_id, t);
  const auto log_skip = [&](const char* reason)
  {
    if (log_flip_diag)
    {
      KINDS_MONITOR("Flip handleEvent SKIP (he_id=" << half_edge_id << ", delaunay_edge=" << (half_edge_id / 2)
                                                 << ", occurrence_t=" << std::setprecision(17) << t
                                                 << ", infinitesimal_t=" << infinitesimal_t
                                                 << ", creation_t=" << creation_time.real_time
                                                 << "+inf=" << creation_time.infinitesimal_time << "): " << reason);
    }
  };

  size_t parent_component_id = static_cast<size_t>(-1);
  std::optional<InfinitesimalComputeContext> infinitesimal_compute;
  if (is_infinitesimal)
  {
    bool epoch_ok = false;
    for (const auto& entry : kd->pending_branch_splits_.by_parent_)
    {
      if (entry.second.infinitesimal_active && entry.second.infinitesimal_epoch == infinitesimal_epoch_)
      {
        epoch_ok = true;
        parent_component_id = entry.first;
        break;
      }
    }
    if (!epoch_ok)
    {
      log_skip("stale infinitesimal epoch");
      return;
    }
    kd->current_infinitesimal_t_ = infinitesimal_t;
    infinitesimal_compute = InfinitesimalComputeContext { infinitesimal_t, parent_component_id };
  }

  if (log_flip_diag)
  {
    const size_t quad_id = half_edge_id / 2;
    const EventTime quad_last = (quad_id < kd->quadrilateral_last_updated.size())
      ? kd->quadrilateral_last_updated[quad_id]
      : EventTime(std::numeric_limits<double>::quiet_NaN());
    KINDS_MONITOR("Flip handleEvent ENTER (he_id=" << half_edge_id << "/" << (half_edge_id ^ 1)
                                                << ", delaunay_edge=" << quad_id << ", occurrence_t="
                                                << std::setprecision(17) << t
                                                << ", creation_t=" << creation_time.real_time
                                                << "+inf=" << creation_time.infinitesimal_time << ", he_live="
                                                << (graph.isLiveHalfEdge(half_edge_id) ? "true" : "false")
                                                << ", quadrilateral_last_updated=" << quad_last.real_time
                                                << "+inf=" << quad_last.infinitesimal_time << ")");
  }

  // Outdated if the flip edge was tombstoned (e.g. after a branch split).
  if (!graph.isLiveHalfEdge(half_edge_id))
  {
    log_skip("half_edge not live");
    if (is_infinitesimal)
    {
      kd->current_infinitesimal_t_ = 0.0;
    }
    return;
  }

  // Check if the event is still valid
  if (creation_time < kd->quadrilateral_last_updated[half_edge_id / 2])
  {
    log_skip("creation_time < quadrilateral_last_updated");
    if (is_infinitesimal)
    {
      kd->current_infinitesimal_t_ = 0.0;
    }
    return;
  }

  if (log_flip_diag)
  {
    KINDS_MONITOR("Flip handleEvent PROCEED (he_id=" << half_edge_id << ", occurrence_t=" << std::setprecision(17)
                                                  << t << ")");
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
  KINDS_DEBUG("Processing flip event at time " << t << " for half-edge ID " << half_edge_id
                                               << ". Faces inside " << kd->face_inside[face_id] << " | "
                                               << kd->face_inside[twin_face_id]);

  auto* event_handler = kd->flip_event_manager_->getCallback();
  if (event_handler)
  {
    event_handler->beforeEvent(*this);
  }

  if (kd->isVisualDebugEnabled() && kd->getVisualDebugOutputRoot().has_value()
    && shouldDumpFlipPolynomialsForEvent(*kd, t, half_edge_id))
  {
    const FlipEventTriggerDump dump = buildFlipEventTriggerDump(*kd, half_edge_id, creation_time.real_time);
    writeFlipEventTriggerPolynomialDump(
      *kd, dump, *kd->getVisualDebugOutputRoot() / "polynomials.txt", t);
  }

  // Sanity-check Voronoi coincidence / boundary collinearity on every flip (log only; never throws).
  // FAIL → WARNING (ungated). OK → MONITOR only under flip diagnostic guards.
  // Intentionally compares both flip-edge Voronoi vertices; do not use
  // @ref canonicalFlipEdgeVoronoiVertexIdForMeshing here.
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
      const glm::dvec2 pa = kd->getPointAt(static_cast<size_t>(a), t);
      const glm::dvec2 pb = kd->getPointAt(static_cast<size_t>(b), t);
      const glm::dvec2 pc = kd->getPointAt(static_cast<size_t>(c), t);
      const double collinearity_metric = normalizedTriangleCollinearityMetric(pa, pb, pc);
      const bool transformed_collinear = collinearity_metric <= flip_boundary_collinearity_eps;

      if (!transformed_collinear)
      {
        const glm::dvec2 pa_raw = kd->getStrandTree().evaluate(static_cast<size_t>(a), t);
        const glm::dvec2 pb_raw = kd->getStrandTree().evaluate(static_cast<size_t>(b), t);
        const glm::dvec2 pc_raw = kd->getStrandTree().evaluate(static_cast<size_t>(c), t);
        const double raw_collinearity_metric = normalizedTriangleCollinearityMetric(pa_raw, pb_raw, pc_raw);
        const bool untransformed_collinear = raw_collinearity_metric <= flip_boundary_collinearity_eps;

        KINDS_WARNING("Flip sanity FAIL boundary collinearity (he_id="
          << half_edge_id << ", occurrence_t=" << std::setprecision(17) << t
          << ", creation_t=" << creation_time
          << ", transformed_collinearity_metric=" << collinearity_metric
          << ", untransformed_collinearity_metric=" << raw_collinearity_metric << ", eps="
          << flip_boundary_collinearity_eps << ", a=" << a << ", b=" << b << ", c=" << c << ", pa="
          << glm::to_string(pa) << ", pb=" << glm::to_string(pb) << ", pc=" << glm::to_string(pc)
          << ", pa_raw=" << glm::to_string(pa_raw) << ", pb_raw=" << glm::to_string(pb_raw)
          << ", pc_raw=" << glm::to_string(pc_raw)
          << flipUntransformedFrameMismatchNote(transformed_collinear, untransformed_collinear) << ")");
      }
      else if (log_flip_diag)
      {
        const glm::dvec2 pa_raw = kd->getStrandTree().evaluate(static_cast<size_t>(a), t);
        const glm::dvec2 pb_raw = kd->getStrandTree().evaluate(static_cast<size_t>(b), t);
        const glm::dvec2 pc_raw = kd->getStrandTree().evaluate(static_cast<size_t>(c), t);
        const double raw_collinearity_metric = normalizedTriangleCollinearityMetric(pa_raw, pb_raw, pc_raw);

        KINDS_MONITOR("Flip sanity OK boundary collinearity (he_id="
          << half_edge_id << ", occurrence_t=" << std::setprecision(17) << t
          << ", creation_t=" << creation_time
          << ", transformed_collinearity_metric=" << collinearity_metric
          << ", untransformed_collinearity_metric=" << raw_collinearity_metric << ")");
      }
    }
  }
  else
  {
    const std::vector<size_t> quad_strand_ids = collectFlipQuadrilateralStrandIds(graph, half_edge_id);
    const size_t shared_reference_branch
      = kd->getSharedReferenceBranchForStrands(quad_strand_ids, t);
    const glm::dvec3 left_voronoi_vertex = kd->computeVoronoiVertexClampedInfinityWithReferenceBranch(
      half_edge_id, t, shared_reference_branch);
    const glm::dvec3 right_voronoi_vertex = kd->computeVoronoiVertexClampedInfinityWithReferenceBranch(
      half_edge_id ^ 1, t, shared_reference_branch);
    const double voronoi_vertex_distance
      = glm::distance(glm::dvec2(left_voronoi_vertex), glm::dvec2(right_voronoi_vertex));
    const bool transformed_coincident = voronoi_vertex_distance <= flip_voronoi_vertex_distance_eps;

    if (!transformed_coincident)
    {
      const glm::dvec2 raw_left_cc
        = flipTriangleCircumcenterAt(*kd, graph, half_edge_id, t, false);
      const glm::dvec2 raw_right_cc
        = flipTriangleCircumcenterAt(*kd, graph, half_edge_id ^ 1, t, false);
      const glm::dvec2 transformed_left_cc = flipTriangleCircumcenterAt(
        *kd, graph, half_edge_id, t, true, shared_reference_branch);
      const glm::dvec2 transformed_right_cc = flipTriangleCircumcenterAt(
        *kd, graph, half_edge_id ^ 1, t, true, shared_reference_branch);
      const double raw_circumcenter_distance = glm::distance(raw_left_cc, raw_right_cc);
      const double shared_frame_circumcenter_distance
        = glm::distance(transformed_left_cc, transformed_right_cc);
      const bool untransformed_coincident = raw_circumcenter_distance <= flip_voronoi_vertex_distance_eps;
      const bool shared_frame_coincident
        = shared_frame_circumcenter_distance <= flip_voronoi_vertex_distance_eps;

      KINDS_WARNING("Flip sanity FAIL Voronoi coincidence (he_id="
        << half_edge_id << ", occurrence_t=" << std::setprecision(17) << t
        << ", creation_t=" << creation_time << ", faces " << face_id << " and " << twin_face_id
        << ", transformed_voronoi_distance=" << voronoi_vertex_distance
        << ", untransformed_circumcenter_distance=" << raw_circumcenter_distance
        << ", shared_frame_circumcenter_distance=" << shared_frame_circumcenter_distance << ", eps="
        << flip_voronoi_vertex_distance_eps << ", shared_reference_branch=" << shared_reference_branch
        << ", left=" << glm::to_string(left_voronoi_vertex)
        << ", right=" << glm::to_string(right_voronoi_vertex) << ", raw_left_cc="
        << glm::to_string(raw_left_cc) << ", raw_right_cc=" << glm::to_string(raw_right_cc)
        << flipUntransformedFrameMismatchNote(transformed_coincident, untransformed_coincident)
        << (shared_frame_coincident && !transformed_coincident
              ? " [shared-frame circumcenters coincide; per-vertex getPointAt frame mismatch]"
              : "")
        << ")");
    }
    else if (log_flip_diag)
    {
      const glm::dvec2 raw_left_cc
        = flipTriangleCircumcenterAt(*kd, graph, half_edge_id, t, false);
      const glm::dvec2 raw_right_cc
        = flipTriangleCircumcenterAt(*kd, graph, half_edge_id ^ 1, t, false);
      const glm::dvec2 transformed_left_cc = flipTriangleCircumcenterAt(
        *kd, graph, half_edge_id, t, true, shared_reference_branch);
      const glm::dvec2 transformed_right_cc = flipTriangleCircumcenterAt(
        *kd, graph, half_edge_id ^ 1, t, true, shared_reference_branch);
      const double raw_circumcenter_distance = glm::distance(raw_left_cc, raw_right_cc);
      const double shared_frame_circumcenter_distance
        = glm::distance(transformed_left_cc, transformed_right_cc);

      KINDS_MONITOR("Flip sanity OK Voronoi coincidence (he_id="
        << half_edge_id << ", occurrence_t=" << std::setprecision(17) << t
        << ", creation_t=" << creation_time
        << ", transformed_voronoi_distance=" << voronoi_vertex_distance
        << ", untransformed_circumcenter_distance=" << raw_circumcenter_distance
        << ", shared_frame_circumcenter_distance=" << shared_frame_circumcenter_distance
        << ", shared_reference_branch=" << shared_reference_branch << ")");
    }
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
      kd->setFaceInside(flipped_face0, true, t);
    }
    else if (is_finite_live_face(flipped_face1))
    {
      kd->setFaceInside(flipped_face1, true, t);
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
      kd->setFaceInside(swapped_face_id, false, t);
    }
  }

  auto tri_verts2 = graph.adjacentTriangleVertices(half_edge_id ^ 1);
  for (auto& v : tri_verts2)
  {
    if (v == -1)
    {
      size_t swapped_face_id = graph.halfEdge(half_edge_id ^ 1).face;
      kd->setFaceInside(swapped_face_id, false, t);
    }
  }

  // After flipping the edge, reassign Voronoi vertices, run callbacks (SVG/mesh), then recompute.
  // Keep current_infinitesimal_t_ set through afterEvent so debug exports / meshing see the virtual shift.
  if (!graph.isOnConvexBoundary(half_edge_id))
  {
    kd->reassignVoronoiVerticesInQuadrilateral(half_edge_id / 2, t, pre_flip_quad_faces, infinitesimal_compute);
  }
  else
  {
    kd->reassignVoronoiVerticesOnBoundary(half_edge_id, t, infinitesimal_compute);
  }

  if (event_handler)
  {
    event_handler->afterEvent(*this);
  }

  // After callbacks (e.g. debug SVG export); intersection lists must be consistent.
  kd->validateVoronoiVertexIteratorInvariants("FlipEvent:afterEvent", t);
  kd->validateCrossingIntersectionInvariants("FlipEvent:afterEvent", t);
  kd->validateSitesInsideConvexHull("FlipEvent:afterEvent", t);

  if (is_infinitesimal)
  {
    // Finalize first: if the cut applies, epoch bump invalidates queued virtual events — no recompute.
    if (kd->maybeFinalizeInfinitesimalSeparation(parent_component_id, t))
    {
      kd->current_infinitesimal_t_ = 0.0;
      return;
    }
  }

  // Local neighbor recompute (same paradigm for kinetic and infinitesimal).
  {
    size_t next1 = graph.halfEdge(half_edge_id).next;
    size_t next2 = graph.halfEdge(next1).next;

    size_t twin_next1 = graph.halfEdge(half_edge_id ^ 1).next;
    size_t twin_next2 = graph.halfEdge(twin_next1).next;

    kd->flip_event_manager_->computeEvents(t, next1 / 2, infinitesimal_compute);
    kd->quadrilateral_last_updated[next1 / 2] = occurrence_time;

    kd->flip_event_manager_->computeEvents(t, next2 / 2, infinitesimal_compute);
    kd->quadrilateral_last_updated[next2 / 2] = occurrence_time;

    kd->flip_event_manager_->computeEvents(t, twin_next1 / 2, infinitesimal_compute);
    kd->quadrilateral_last_updated[twin_next1 / 2] = occurrence_time;

    kd->flip_event_manager_->computeEvents(t, twin_next2 / 2, infinitesimal_compute);
    kd->quadrilateral_last_updated[twin_next2 / 2] = occurrence_time;

    // re-compute radius events for both triangles
    kd->radius_event_manager_->computeEvents(t, half_edge_id, infinitesimal_compute);
    kd->face_last_updated[face_id] = occurrence_time;

    kd->radius_event_manager_->computeEvents(t, half_edge_id ^ 1, infinitesimal_compute);
    kd->face_last_updated[twin_face_id] = occurrence_time;
  }

  if (is_infinitesimal)
  {
    kd->current_infinitesimal_t_ = 0.0;
  }
}

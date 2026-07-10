#include "KineticDelaunayFlipEvent.hpp"

#include <algorithm>
#include <limits>

using namespace kinDS;

#include "KineticDelaunayEventPredicates.hpp"

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

std::optional<double> separationRampEndFractionForQuad(
  const KineticDelaunay& kd, const std::vector<size_t>& quad_strand_ids, size_t section)
{
  for (size_t strand_id : quad_strand_ids)
  {
    if (const std::optional<double> ramp_end = kd.separationRampEndSectionFraction(strand_id, section);
      ramp_end.has_value())
    {
      return ramp_end;
    }
  }
  return std::nullopt;
}
} // namespace

void KineticDelaunay::FlipEventManager::computeEvents(double t, size_t quad_id)
{
  auto* kd = kd_;
  auto& graph = kd->graph;

  const size_t section = static_cast<size_t>(t);
  const float fraction = t - section;

  size_t he_id = quad_id * 2;
  const std::vector<size_t> quad_strand_ids = collectFlipQuadrilateralStrandIds(graph, he_id);

  const auto build_trigger = [&](size_t active_he_id, double schedule_time, Polynomial& event_trigger_out,
                                 std::vector<Trajectory<2>>& trajs_out) {
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
                                    double min_fraction, size_t enqueue_he_id, double creation_time) {
    if (event_trigger.degree() < 0)
    {
      return;
    }

    auto event_times = kd->findEvents(event_trigger, min_fraction);
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

      KINDS_DEBUG("Scheduled flip event at time " << event_time + section << " for half-edge ID " << enqueue_he_id
                                                << " at center position " << glm::to_string(center));

      kd->kinetic_algorithm_->enqueueEvent(
        std::make_shared<FlipEvent>(kd, event_time + section, enqueue_he_id, creation_time, center));
    }
  };

  Polynomial event_trigger;
  std::vector<Trajectory<2>> trajs;
  build_trigger(he_id, t, event_trigger, trajs);
  enqueue_flip_roots(trajs, event_trigger, fraction, he_id, t);

  if (const std::optional<double> ramp_end_fraction = separationRampEndFractionForQuad(*kd, quad_strand_ids, section);
    ramp_end_fraction.has_value() && fraction < *ramp_end_fraction && *ramp_end_fraction < kEventIntervalFractionUpperBound)
  {
    const double post_ramp_schedule_time = static_cast<double>(section) + *ramp_end_fraction;
    Polynomial post_ramp_trigger;
    std::vector<Trajectory<2>> post_ramp_trajs;
    build_trigger(he_id, post_ramp_schedule_time, post_ramp_trigger, post_ramp_trajs);
    enqueue_flip_roots(post_ramp_trajs, post_ramp_trigger, *ramp_end_fraction, he_id, t);
  }
}

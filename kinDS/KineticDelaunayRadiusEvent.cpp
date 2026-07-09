#include "KineticDelaunayRadiusEvent.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

using namespace kinDS;

#include "KineticDelaunayEventPredicates.hpp"

namespace
{
struct SignedRadiusRoot
{
  double fraction = 0.0;
  bool target_inside = false;
};

std::vector<SignedRadiusRoot> findSignedRadiusRoots(Polynomial& event_trigger, double min_fraction)
{
  if (event_trigger.degree() == -1)
  {
    return {};
  }

  event_trigger.trim();
  auto zeros = event_trigger.realRoots();
  std::vector<double> filtered_sorted_zeros;
  for (double root : zeros)
  {
    if (!std::isnan(root) && root > min_fraction && root <= 1.0)
    {
      filtered_sorted_zeros.push_back(root);
    }
  }

  if (filtered_sorted_zeros.empty())
  {
    return {};
  }

  std::sort(filtered_sorted_zeros.begin(), filtered_sorted_zeros.end());

  std::vector<double> interval_signs(filtered_sorted_zeros.size() + 1);
  double test_point = (min_fraction + filtered_sorted_zeros[0]) / 2.0;
  interval_signs[0] = event_trigger(test_point) > 0.0 ? 1.0 : -1.0;
  for (size_t i = 0; i < filtered_sorted_zeros.size(); ++i)
  {
    test_point
      = (filtered_sorted_zeros[i] + (i + 1 < filtered_sorted_zeros.size() ? filtered_sorted_zeros[i + 1] : 1.0)) / 2.0;
    interval_signs[i + 1] = event_trigger(test_point) > 0.0 ? 1.0 : -1.0;
  }

  std::vector<SignedRadiusRoot> roots;
  for (size_t i = 0; i < filtered_sorted_zeros.size(); ++i)
  {
    if (interval_signs[i] == interval_signs[i + 1])
    {
      KINDS_DEBUG("Radius root has no sign change at fraction " << filtered_sorted_zeros[i] << ", skipping.");
      continue;
    }

    // Positive -> negative means the circumradius predicate moved into the inside range: add the triangle.
    const bool target_inside = interval_signs[i] > 0.0 && interval_signs[i + 1] < 0.0;
    roots.push_back(SignedRadiusRoot { filtered_sorted_zeros[i], target_inside });
    KINDS_DEBUG("Radius event will be queued at fraction " << filtered_sorted_zeros[i] << " with sign change from "
                                                          << interval_signs[i] << " to " << interval_signs[i + 1]
                                                          << " target_inside=" << target_inside);
  }
  return roots;
}
} // namespace

void KineticDelaunay::RadiusEventManager::computeEvents(double t, size_t he_id)
{
  auto* kd = kd_;
  auto& graph = kd->graph;
  auto& branch_trajs = kd->branch_trajs;

  if (kd->cutoff == std::numeric_limits<double>::infinity())
  {
    // no radius events wanted
    return;
  }

  const size_t section = static_cast<size_t>(t);
  const float fraction = t - section;

  size_t face_id = graph.halfEdge(he_id).face;
  size_t u = graph.halfEdge(he_id).origin;
  size_t v = graph.destination(he_id);
  size_t w = graph.triangleOppositeVertex(he_id);

  if (u == -1 || v == -1 || w == -1)
  {
    // one of the vertices is at infinity, no event possible
    return;
  }

  if (kd->mustRemainInside(face_id, t))
  {
    return;
  }

  const auto piece_poly = [&](size_t strand_id) { return kd->getSitePiecePolynomial(strand_id, section, t); };

  std::vector<Trajectory<2>> trajs;

  trajs.push_back(piece_poly(static_cast<size_t>(u)));
  trajs.push_back(piece_poly(static_cast<size_t>(v)));
  trajs.push_back(piece_poly(static_cast<size_t>(w)));

  Polynomial event_trigger
    = circumradiusEquals(
      trajs[0][0], trajs[0][1], trajs[1][0], trajs[1][1], trajs[2][0], trajs[2][1], kd->cutoff);

  auto event_times = findSignedRadiusRoots(event_trigger, fraction);
  for (const auto& event_time : event_times)
  {
    glm::dvec2 center {};

    for (const auto& traj : trajs)
    {
      center[0] += traj[0](event_time.fraction);
      center[1] += traj[1](event_time.fraction);
    }
    center[0] /= trajs.size();
    center[1] /= trajs.size();
    // KINDS_DEBUG("Boundary Event at time " << event_time + section << " for half-edge ID " << he_id << " at center
    // position"
    //                                       << glm::to_string(center));

    kd->kinetic_algorithm_->enqueueEvent(
      std::make_shared<RadiusEvent>(kd, event_time.fraction + section, he_id, t, center, event_time.target_inside));
  }
}
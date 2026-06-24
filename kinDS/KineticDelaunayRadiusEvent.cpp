#include "KineticDelaunayRadiusEvent.hpp"

#include <limits>

using namespace kinDS;

#include "KineticDelaunayEventPredicates.hpp"

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

  const size_t reference_branch = kd->getReferenceBranch(static_cast<size_t>(u), t);
  const auto piece_poly = [&](size_t strand_id)
  { return branch_trajs.getPiecePolynomial(strand_id, section, reference_branch); };

  std::vector<Trajectory<2>> trajs;

  trajs.push_back(piece_poly(static_cast<size_t>(u)));
  trajs.push_back(piece_poly(static_cast<size_t>(v)));
  trajs.push_back(piece_poly(static_cast<size_t>(w)));

  Polynomial event_trigger
    = circumradiusEquals(
      trajs[0][0], trajs[0][1], trajs[1][0], trajs[1][1], trajs[2][0], trajs[2][1], kd->cutoff);

  auto event_times = kd->findEvents(event_trigger, fraction);
  for (const auto& event_time : event_times)
  {
    glm::dvec2 center {};

    for (const auto& traj : trajs)
    {
      center[0] += traj[0](event_time);
      center[1] += traj[1](event_time);
    }
    center[0] /= trajs.size();
    center[1] /= trajs.size();
    // KINDS_DEBUG("Boundary Event at time " << event_time + section << " for half-edge ID " << he_id << " at center
    // position"
    //                                       << glm::to_string(center));

    kd->kinetic_algorithm_->enqueueEvent(std::make_shared<RadiusEvent>(kd, event_time + section, he_id, t, center));
  }
}
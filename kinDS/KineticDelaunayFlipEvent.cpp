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
} // namespace

void KineticDelaunay::FlipEventManager::computeEvents(double t, size_t quad_id)
{
  auto* kd = kd_;
  auto& graph = kd->graph;
  auto& branch_trajs = kd->branch_trajs;

  const size_t section = static_cast<size_t>(t);
  const float fraction = t - section;

  size_t he_id = quad_id * 2;
  Polynomial event_trigger;

  const double branch_lookup_time = referenceBranchLookupTimeForSection(section, t);
  const std::vector<size_t> quad_strand_ids = collectFlipQuadrilateralStrandIds(graph, he_id);
  const size_t shared_reference_branch = kd->getSharedReferenceBranchForStrands(quad_strand_ids, branch_lookup_time);

  std::vector<Trajectory<2>> trajs;
  const auto piece_poly = [&](size_t strand_id) {
    return kd->getSitePiecePolynomialWithReferenceBranch(strand_id, section, shared_reference_branch);
  };

  if (graph.isOnConvexBoundary(he_id) || graph.isOutsideConvexBoundary(he_id))
  {
    // boundary edges must be treated separately using ccw

    // need to get the inner half-edge so we have access to the triangle
    // in case both are outside, this swap does not matter, so we just let it happen
    if (graph.isOutsideConvexBoundary(he_id))
    {
      he_id = he_id ^ 1; // use the twin half-edge if the current one is on the boundary
    }

    // Depending on the half-edge, the infinite vertex could be in different places, so we just collect all and filter
    // it out
    int indices[4];
    indices[0] = graph.halfEdge(he_id).origin; // First vertex
    indices[1] = graph.triangleOppositeVertex(he_id ^ 1); // Second vertex
    indices[2] = graph.halfEdge(he_id ^ 1).origin; // Third vertex
    indices[3] = graph.triangleOppositeVertex(he_id); // Fourth vertex

    std::vector<int> filtered_indices;

    std::copy_if(indices, indices + 4, std::back_inserter(filtered_indices), [](int index) { return index != -1; });

    int& a = filtered_indices[0]; // First vertex
    int& b = filtered_indices[1]; // Second vertex
    int& c = filtered_indices[2]; // Third vertex

    // print the triangle vertices:
    // std::cout << "Triangle vertices: " << a << ", " << b << ", " << c << std::endl;

    trajs.push_back(piece_poly(static_cast<size_t>(a)));
    trajs.push_back(piece_poly(static_cast<size_t>(b)));
    trajs.push_back(piece_poly(static_cast<size_t>(c)));

    event_trigger = ccw(trajs[0][0], trajs[0][1], trajs[1][0], trajs[1][1], trajs[2][0], trajs[2][1]);
  }
  else
  {
    int a = graph.halfEdge(he_id).origin; // First vertex
    int b = graph.triangleOppositeVertex(he_id ^ 1); // Second vertex
    int c = graph.halfEdge(he_id ^ 1).origin; // Third vertex
    int d = graph.triangleOppositeVertex(he_id); // Fourth vertex

    // print the quadrilateral vertices:
    // std::cout << "Quadrilateral vertices: " << a << ", " << b << ", " << c << ", " << d << std::endl;

    trajs.push_back(piece_poly(static_cast<size_t>(a)));
    trajs.push_back(piece_poly(static_cast<size_t>(b)));
    trajs.push_back(piece_poly(static_cast<size_t>(c)));
    trajs.push_back(piece_poly(static_cast<size_t>(d)));

    event_trigger = inCircle(
      trajs[0][0], trajs[0][1], trajs[1][0], trajs[1][1], trajs[2][0], trajs[2][1], trajs[3][0], trajs[3][1]);
  }

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

    KINDS_DEBUG("Scheduled flip event at time " << event_time + section << " for half-edge ID " << he_id << " at center position "
                                      << glm::to_string(center));

    kd->kinetic_algorithm_->enqueueEvent(std::make_shared<FlipEvent>(kd, event_time + section, he_id, t, center));
  }
}
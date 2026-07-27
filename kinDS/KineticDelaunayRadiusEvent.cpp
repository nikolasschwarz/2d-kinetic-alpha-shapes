#include "KineticDelaunayRadiusEvent.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <limits>
#include <sstream>
#include <utility>
#include <vector>

#include "KineticDelaunayEventPredicates.hpp"

using namespace kinDS;

namespace
{
struct SignedRadiusRoot
{
  double fraction = 0.0;
  bool target_inside = false;
};

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
    if (!std::isnan(root) && root > min_fraction && root <= kEventIntervalFractionUpperBound)
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

void logMonitoredFaceRadiusTrajectories(const KineticDelaunay& kd, size_t section, double schedule_t,
  const std::array<size_t, 3>& strand_ids, const std::array<Trajectory<2>, 3>& trajectories)
{
  for (size_t vertex_index = 0; vertex_index < 3; ++vertex_index)
  {
    if (!kd.isDiagnosticsStrandIdValid(strand_ids[vertex_index]))
    {
      return;
    }
  }

  const StrandTree& tree = kd.getStrandTree();
  for (size_t vertex_index = 0; vertex_index < 3; ++vertex_index)
  {
    const size_t strand_id = strand_ids[vertex_index];
    const auto& support_points = tree.getSupportPoints(strand_id);
    if (section >= support_points.size())
    {
      return;
    }
  }

  const std::vector<size_t> event_strand_ids = { strand_ids[0], strand_ids[1], strand_ids[2] };
  const double event_interval_upper_bound = eventIntervalUpperBound(schedule_t);
  const size_t branch_section_index = kd.inputBranchSectionIndexAtIntervalUpperBound(event_interval_upper_bound);
  const std::vector<size_t> distinct_input_branches
    = kd.collectDistinctInputBranchesForEventTrigger(event_strand_ids, event_interval_upper_bound);
  const bool use_shared_transformed_frame
    = kd.eventTriggerUsesSharedTransformedFrame(event_strand_ids, event_interval_upper_bound);
  std::ostringstream branch_summary;
  branch_summary << "  event_trigger_strands=" << event_strand_ids[0] << "," << event_strand_ids[1] << ","
                 << event_strand_ids[2] << " schedule_t=" << schedule_t
                 << " event_interval_upper_bound=" << event_interval_upper_bound
                 << " branch_section_index=" << branch_section_index << " distinct_input_branches=[";
  for (size_t i = 0; i < distinct_input_branches.size(); ++i)
  {
    if (i > 0)
    {
      branch_summary << ',';
    }
    branch_summary << distinct_input_branches[i];
  }
  branch_summary << "] uses_shared_transformed_frame=" << (use_shared_transformed_frame ? "true" : "false");
  if (use_shared_transformed_frame)
  {
    branch_summary << " shared_reference_branch="
                   << kd.sharedReferenceBranchForEventTrigger(event_strand_ids, event_interval_upper_bound);
  }
  branch_summary << " frame_policy=" << (use_shared_transformed_frame ? "shared_transformed" : "local_support");
  KINDS_DEBUG(branch_summary.str());

  for (size_t vertex_index = 0; vertex_index < 3; ++vertex_index)
  {
    const size_t strand_id = strand_ids[vertex_index];
    const auto& support_points = tree.getSupportPoints(strand_id);
    const size_t input_branch_at_section = tree.getBranchIndex(strand_id, section);
    const glm::dvec2 raw_section_start = support_points[section];
    const glm::dvec2 raw_section_end
      = (section + 1 < support_points.size()) ? support_points[section + 1] : support_points[section];
    KINDS_DEBUG("  trajectory[" << vertex_index << "] strand=" << strand_id << " input_branch_at_section="
                               << input_branch_at_section << " raw_support_at_section=(" << raw_section_start.x
                               << "," << raw_section_start.y << ") raw_support_at_section_end=(" << raw_section_end.x
                               << "," << raw_section_end.y << ")");

    for (double eval_t : { static_cast<double>(section), static_cast<double>(section) + 1.0 })
    {
      if (eval_t > static_cast<double>(section) && section + 1 >= support_points.size())
      {
        continue;
      }
      const double fraction = eval_t - static_cast<double>(section);
      const double x = trajectories[vertex_index][0](fraction);
      const double y = trajectories[vertex_index][1](fraction);
      std::ostringstream line;
      line << std::setprecision(17) << "  trajectory[" << vertex_index << "] strand=" << strand_id << " at t="
           << eval_t << " fraction=" << fraction << " pos=(" << x << "," << y << ")";
      KINDS_DEBUG(line.str());
    }
  }
}

void logRadiusTriggerRootsForMonitoredFace(const KineticDelaunay& kd, size_t face_id, size_t he_id, double t,
  double min_fraction, Polynomial event_trigger, const std::array<size_t, 3>& strand_ids,
  const std::array<Trajectory<2>, 3>& trajectories)
{
  if (!kd.diagnosticsEnabled()
    || !KineticDelaunay::matchesDiagnosticsMonitorId(face_id, KineticDelaunay::kDiagnosticsMonitoredFaceId))
  {
    return;
  }
  if (!kd.isDiagnosticsFaceIdValid(face_id))
  {
    return;
  }

  const size_t section = static_cast<size_t>(t);
  const double event_interval_upper_bound = eventIntervalUpperBound(t);
  std::ostringstream header;
  header << "Radius trigger roots monitor (face " << face_id << ", he_id=" << he_id << ", schedule_t=" << t
         << ", event_interval_upper_bound=" << event_interval_upper_bound << ", section=" << section
         << ", min_fraction=" << min_fraction << ", cutoff=" << kd.getCutoff()
         << ", trigger_degree=" << event_trigger.degree() << ")";
  KINDS_DEBUG(header.str());
  logMonitoredFaceRadiusTrajectories(kd, section, t, strand_ids, trajectories);

  if (event_trigger.degree() == -1)
  {
    KINDS_DEBUG("  trigger empty (degree -1)");
    return;
  }

  event_trigger.trim();
  const auto zeros = event_trigger.realRoots();
  if (zeros.empty())
  {
    KINDS_DEBUG("  no real roots");
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

  for (size_t root_index = 0; root_index < zeros.size(); ++root_index)
  {
    const double root = zeros[root_index];
    std::ostringstream line;
    line << std::setprecision(17) << "  root[" << root_index << "] fraction=" << root << " absolute_t="
         << (root + static_cast<double>(section));

    if (std::isnan(root))
    {
      line << " filtered (nan)";
      KINDS_DEBUG(line.str());
      continue;
    }

    const auto [sign_before, sign_after] = signChangeAtRoot(event_trigger, root, sorted_zeros);
    line << " sign_change=" << formatSignChange(sign_before, sign_after);

    if (root <= min_fraction)
    {
      line << " filtered (fraction <= min_fraction)";
      KINDS_DEBUG(line.str());
      continue;
    }
    if (root > kEventIntervalFractionUpperBound)
    {
      line << " filtered (fraction > " << kEventIntervalFractionUpperBound << ")";
      KINDS_DEBUG(line.str());
      continue;
    }

    if (sign_before == sign_after)
    {
      line << " filtered (no_sign_change)";
      KINDS_DEBUG(line.str());
      continue;
    }

    const bool target_inside = sign_before > 0.0 && sign_after < 0.0;
    line << " enqueued target_inside=" << (target_inside ? "true" : "false");
    KINDS_DEBUG(line.str());
  }
}
} // namespace

void KineticDelaunay::setDiagnosticsEnabled(bool enabled) { diagnostics_enabled_ = enabled; }

bool KineticDelaunay::diagnosticsEnabled() const { return diagnostics_enabled_; }

void KineticDelaunay::setSitesInsideConvexHullCheckEnabled(bool enabled)
{
  sites_inside_convex_hull_check_enabled_ = enabled;
}

bool KineticDelaunay::sitesInsideConvexHullCheckEnabled() const
{
  return sites_inside_convex_hull_check_enabled_;
}

bool KineticDelaunay::isDiagnosticsStrandIdValid(size_t strand_id) const
{
  return isDiagnosticsMonitorIdEnabled(strand_id) && strand_id < graph.getVertexCount();
}

bool KineticDelaunay::isDiagnosticsFaceIdValid(size_t face_id) const
{
  return isDiagnosticsMonitorIdEnabled(face_id) && face_id < graph.faceSlotCount();
}

bool KineticDelaunay::isDiagnosticsHalfEdgeIdValid(size_t half_edge_id) const
{
  return isDiagnosticsMonitorIdEnabled(half_edge_id) && graph.isLiveHalfEdge(half_edge_id);
}

bool KineticDelaunay::isDiagnosticsMonitoredFaceValid() const
{
  return isDiagnosticsMonitorIdEnabled(kDiagnosticsMonitoredFaceId)
    && isDiagnosticsFaceIdValid(kDiagnosticsMonitoredFaceId);
}

void KineticDelaunay::logRadiusEventTriggerRoots(size_t face_id, size_t he_id, double t, double min_fraction,
  Polynomial event_trigger, const std::array<size_t, 3>& strand_ids,
  const std::array<Trajectory<2>, 3>& trajectories) const
{
  logRadiusTriggerRootsForMonitoredFace(
    *this, face_id, he_id, t, min_fraction, std::move(event_trigger), strand_ids, trajectories);
}

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

  const std::vector<size_t> event_strand_ids
    = { static_cast<size_t>(u), static_cast<size_t>(v), static_cast<size_t>(w) };
  const auto piece_poly = [&](size_t strand_id) {
    return kd->getSitePiecePolynomialForEventStrands(strand_id, section, t, event_strand_ids);
  };

  std::vector<Trajectory<2>> trajs;

  trajs.push_back(piece_poly(static_cast<size_t>(u)));
  trajs.push_back(piece_poly(static_cast<size_t>(v)));
  trajs.push_back(piece_poly(static_cast<size_t>(w)));

  Polynomial event_trigger
    = circumradiusEquals(
      trajs[0][0], trajs[0][1], trajs[1][0], trajs[1][1], trajs[2][0], trajs[2][1], kd->cutoff);

  if (kd->diagnosticsEnabled()
    && KineticDelaunay::matchesDiagnosticsMonitorId(face_id, KineticDelaunay::kDiagnosticsMonitoredFaceId)
    && kd->isDiagnosticsFaceIdValid(face_id))
  {
    kd->logRadiusEventTriggerRoots(face_id, he_id, t, fraction, event_trigger,
      { static_cast<size_t>(u), static_cast<size_t>(v), static_cast<size_t>(w) },
      { trajs[0], trajs[1], trajs[2] });
  }

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

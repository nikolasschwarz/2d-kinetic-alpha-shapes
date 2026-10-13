#include "KineticDelaunayFlipEventTriggerDump.hpp"

#include "KineticDelaunayEventPredicates.hpp"
#include "KineticDelaunayFlipEvent.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <sstream>

namespace kinDS
{
namespace
{
std::string formatCoefficient(double value)
{
  std::ostringstream oss;
  oss << std::setprecision(17) << value;
  return oss.str();
}

std::string formatPolynomialHumanReadable(const Polynomial& polynomial, const char* variable = "X")
{
  const Eigen::VectorXd& coeffs = polynomial.getCoefficients();
  if (coeffs.size() == 0)
  {
    return "0";
  }

  std::ostringstream oss;
  bool first_term = true;
  for (int exponent = static_cast<int>(coeffs.size()) - 1; exponent >= 0; --exponent)
  {
    const double coefficient = coeffs[exponent];
    if (coefficient == 0.0)
    {
      continue;
    }

    if (!first_term)
    {
      oss << (coefficient > 0.0 ? " + " : " - ");
    }
    else if (coefficient < 0.0)
    {
      oss << '-';
    }

    const double magnitude = std::abs(coefficient);
    if (exponent == 0)
    {
      oss << formatCoefficient(magnitude);
    }
    else if (magnitude == 1.0)
    {
      oss << variable;
      if (exponent > 1)
      {
        oss << '^' << exponent;
      }
    }
    else
    {
      oss << formatCoefficient(magnitude) << '*' << variable;
      if (exponent > 1)
      {
        oss << '^' << exponent;
      }
    }

    first_term = false;
  }

  const std::string result = oss.str();
  return result.empty() ? "0" : result;
}

std::string formatTrajectoryHumanReadable(const Trajectory<2>& trajectory, const char* variable = "X")
{
  return "x(" + std::string(variable) + ") = " + formatPolynomialHumanReadable(trajectory[0], variable) + "\n  y("
    + std::string(variable) + ") = " + formatPolynomialHumanReadable(trajectory[1], variable);
}

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

Trajectory<2> sitePiecePolynomialAtScheduleTime(const KineticDelaunay& kd, size_t strand_id, size_t section,
  double schedule_time, const std::vector<size_t>& event_strand_ids)
{
  return kd.getSitePiecePolynomialForEventStrands(strand_id, section, schedule_time, event_strand_ids);
}

void appendSiteTrajectorySection(std::ostringstream& oss, const KineticDelaunay& kd, int vertex_id, size_t section,
  double schedule_time, const Trajectory<2>& trajectory, char label, const std::vector<size_t>& event_strand_ids)
{
  if (vertex_id < 0)
  {
    return;
  }

  const size_t strand_id = static_cast<size_t>(vertex_id);
  const double event_interval_upper_bound = eventIntervalUpperBound(schedule_time);
  const size_t branch_section_index = kd.inputBranchSectionIndexAtIntervalUpperBound(event_interval_upper_bound);
  const size_t input_branch_at_section = kd.getStrandTree().getBranchIndex(strand_id, branch_section_index);
  const std::vector<size_t> distinct_input_branches
    = kd.collectDistinctInputBranchesForEventTrigger(event_strand_ids, event_interval_upper_bound);
  const bool use_shared_transformed_frame
    = kd.eventTriggerUsesSharedTransformedFrame(event_strand_ids, event_interval_upper_bound);
  const size_t runtime_branch = kd.getRuntimeBranchIdForStrand(strand_id);
  const size_t component_id = kd.getInsideFaceComponentId(strand_id);

  oss << "\nSite " << label << " (strand " << vertex_id << ")\n";
  oss << "  component_id=" << component_id << '\n';
  oss << "  runtime_branch=" << runtime_branch << '\n';
  oss << "  input_branch_at_section=" << input_branch_at_section << " (branch_section_index=" << branch_section_index
      << ")\n";
  oss << "  distinct_input_branches=[";
  for (size_t i = 0; i < distinct_input_branches.size(); ++i)
  {
    if (i > 0)
    {
      oss << ',';
    }
    oss << distinct_input_branches[i];
  }
  oss << "]\n";
  oss << "  uses_shared_transformed_frame=" << (use_shared_transformed_frame ? "true" : "false") << '\n';
  if (use_shared_transformed_frame)
  {
    oss << "  shared_reference_branch_for_predicate="
        << kd.sharedReferenceBranchForEventTrigger(event_strand_ids, event_interval_upper_bound) << '\n';
  }
  oss << "  reference_branch_for_motion=" << kd.getReferenceBranch(strand_id, schedule_time) << '\n';
  oss << "  schedule_time=" << std::setprecision(17) << schedule_time << '\n';
  oss << "  event_interval_upper_bound=" << event_interval_upper_bound << '\n';
  oss << "  section=" << section << '\n';
  oss << "  " << formatTrajectoryHumanReadable(trajectory) << '\n';
}

void appendInCircleComposition(std::ostringstream& oss)
{
  oss << "\nComposition (Guibas-Stolfi inCircle test; X is the section fraction parameter)\n";
  oss << "  Let a,b,c be three sites and p the fourth site.\n";
  oss << "  dx = ax - px,  dy = ay - py\n";
  oss << "  ex = bx - px,  ey = by - py\n";
  oss << "  fx = cx - px,  fy = cy - py\n";
  oss << "  ap = dx^2 + dy^2,  bp = ex^2 + ey^2,  cp = fx^2 + fy^2\n";
  oss << "  event_trigger(X) = dx*(ey*cp - bp*fy) - dy*(ex*cp - bp*fx) + ap*(ex*fy - ey*fx)\n";
  oss << "  Flip roots are real roots of event_trigger(X) in (section_fraction, "
      << kEventIntervalFractionUpperBound << "].\n";
}

void appendCcwComposition(std::ostringstream& oss)
{
  oss << "\nComposition (ccw orientation test; X is the section fraction parameter)\n";
  oss << "  event_trigger(X) = ax*by + bx*cy + cx*ay - ay*bx - by*cx - cy*ax\n";
  oss << "  Flip roots are real roots of event_trigger(X) in (section_fraction, "
      << kEventIntervalFractionUpperBound << "].\n";
}
} // namespace

FlipEventTriggerDump buildFlipEventTriggerDump(const KineticDelaunay& kd, size_t half_edge_id, double schedule_time)
{
  const auto& graph = kd.getGraph();
  const size_t section = static_cast<size_t>(schedule_time);
  const double section_fraction = schedule_time - static_cast<double>(section);
  const std::vector<size_t> quad_strand_ids = collectFlipQuadrilateralStrandIds(graph, half_edge_id);

  FlipEventTriggerDump dump;
  dump.section = section;
  dump.schedule_time = schedule_time;
  dump.section_fraction = section_fraction;
  dump.half_edge_id = half_edge_id;

  size_t he_id = half_edge_id;
  if (graph.isOnConvexBoundary(he_id) || graph.isOutsideConvexBoundary(he_id))
  {
    dump.boundary_ccw = true;
    if (graph.isOutsideConvexBoundary(he_id))
    {
      he_id ^= 1;
    }

    int indices[4] = {
      graph.halfEdge(he_id).origin,
      graph.triangleOppositeVertex(he_id ^ 1),
      graph.halfEdge(he_id ^ 1).origin,
      graph.triangleOppositeVertex(he_id),
    };

    std::vector<int> filtered_indices;
    std::copy_if(indices, indices + 4, std::back_inserter(filtered_indices), [](int index) { return index != -1; });

    dump.vertex_count = filtered_indices.size();
    for (size_t i = 0; i < filtered_indices.size() && i < dump.vertex_ids.size(); ++i)
    {
      dump.vertex_ids[i] = filtered_indices[i];
      dump.site_trajectories.push_back(sitePiecePolynomialAtScheduleTime(kd,
        static_cast<size_t>(filtered_indices[i]), section, schedule_time, quad_strand_ids));
    }

    if (filtered_indices.size() >= 3)
    {
      int& a = filtered_indices[0];
      int& b = filtered_indices[1];
      int& c = filtered_indices[2];
      dump.event_trigger = ccw(dump.site_trajectories[0][0], dump.site_trajectories[0][1], dump.site_trajectories[1][0],
        dump.site_trajectories[1][1], dump.site_trajectories[2][0], dump.site_trajectories[2][1]);
      (void)a;
      (void)b;
      (void)c;
    }
  }
  else
  {
    dump.boundary_ccw = false;
    dump.vertex_ids[0] = graph.halfEdge(he_id).origin;
    dump.vertex_ids[1] = graph.triangleOppositeVertex(he_id ^ 1);
    dump.vertex_ids[2] = graph.halfEdge(he_id ^ 1).origin;
    dump.vertex_ids[3] = graph.triangleOppositeVertex(he_id);
    dump.vertex_count = 4;

    for (size_t i = 0; i < dump.vertex_count; ++i)
    {
      dump.site_trajectories.push_back(sitePiecePolynomialAtScheduleTime(kd,
        static_cast<size_t>(dump.vertex_ids[i]), section, schedule_time, quad_strand_ids));
    }

    dump.event_trigger = inCircle(dump.site_trajectories[0][0], dump.site_trajectories[0][1], dump.site_trajectories[1][0],
      dump.site_trajectories[1][1], dump.site_trajectories[2][0], dump.site_trajectories[2][1],
      dump.site_trajectories[3][0], dump.site_trajectories[3][1]);
  }

  return dump;
}

void writeFlipEventTriggerPolynomialDump(
  const KineticDelaunay& kd, const FlipEventTriggerDump& dump, const std::filesystem::path& filepath,
  double occurrence_time)
{
  const double event_fraction = occurrence_time - static_cast<double>(dump.section);

  std::ostringstream oss;
  oss << std::setprecision(17);
  oss << "# Flip event trigger polynomial dump\n";
  oss << "half_edge_id: " << dump.half_edge_id << '\n';
  oss << "occurrence_time: " << occurrence_time << '\n';
  oss << "schedule_time_used_for_polynomials: " << dump.schedule_time << '\n';
  const double event_interval_upper_bound = eventIntervalUpperBound(dump.schedule_time);
  oss << "event_interval_upper_bound: " << event_interval_upper_bound << '\n';
  const double branch_lookup_time = referenceBranchLookupTimeForSection(dump.section, dump.schedule_time);
  const std::vector<size_t> quad_strand_ids
    = collectFlipQuadrilateralStrandIds(kd.getGraph(), dump.half_edge_id);
  oss << "event_trigger_strands: ";
  for (size_t i = 0; i < quad_strand_ids.size(); ++i)
  {
    if (i > 0)
    {
      oss << ',';
    }
    oss << quad_strand_ids[i];
  }
  oss << '\n';
  const std::vector<size_t> distinct_input_branches
    = kd.collectDistinctInputBranchesForEventTrigger(quad_strand_ids, event_interval_upper_bound);
  const size_t branch_section_index = kd.inputBranchSectionIndexAtIntervalUpperBound(event_interval_upper_bound);
  oss << "distinct_input_branches: [";
  for (size_t i = 0; i < distinct_input_branches.size(); ++i)
  {
    if (i > 0)
    {
      oss << ',';
    }
    oss << distinct_input_branches[i];
  }
  oss << "]\n";
  oss << "branch_section_index: " << branch_section_index << '\n';
  oss << "uses_shared_transformed_frame: "
      << (kd.eventTriggerUsesSharedTransformedFrame(quad_strand_ids, event_interval_upper_bound) ? "true" : "false")
      << '\n';
  if (kd.eventTriggerUsesSharedTransformedFrame(quad_strand_ids, event_interval_upper_bound))
  {
    oss << "shared_reference_branch_for_predicate: "
        << kd.sharedReferenceBranchForEventTrigger(quad_strand_ids, event_interval_upper_bound) << '\n';
  }
  oss << "section: " << dump.section << '\n';
  oss << "section_fraction_at_schedule_time: " << dump.section_fraction << '\n';
  oss << "section_fraction_at_occurrence_time: " << event_fraction << '\n';
  oss << "predicate: " << (dump.boundary_ccw ? "ccw (boundary)" : "inCircle (interior)") << '\n';

  if (dump.boundary_ccw)
  {
    oss << "triangle_vertices: ";
    for (size_t i = 0; i < dump.vertex_count; ++i)
    {
      if (i > 0)
      {
        oss << ", ";
      }
      oss << dump.vertex_ids[i];
    }
    oss << '\n';
    appendCcwComposition(oss);
    const char labels[] = { 'a', 'b', 'c' };
    for (size_t i = 0; i < dump.site_trajectories.size() && i < 3; ++i)
    {
      appendSiteTrajectorySection(
        oss, kd, dump.vertex_ids[i], dump.section, dump.schedule_time, dump.site_trajectories[i], labels[i],
        quad_strand_ids);
    }
  }
  else
  {
    oss << "quadrilateral_vertices: a=" << dump.vertex_ids[0] << ", b=" << dump.vertex_ids[1] << ", c="
        << dump.vertex_ids[2] << ", d=" << dump.vertex_ids[3] << " (inCircle uses a,b,c,p with p=d)\n";
    appendInCircleComposition(oss);
    const char labels[] = { 'a', 'b', 'c', 'd' };
    for (size_t i = 0; i < dump.site_trajectories.size() && i < 4; ++i)
    {
      appendSiteTrajectorySection(
        oss, kd, dump.vertex_ids[i], dump.section, dump.schedule_time, dump.site_trajectories[i], labels[i],
        quad_strand_ids);
    }
  }

  oss << "\nevent_trigger(X) = " << formatPolynomialHumanReadable(dump.event_trigger) << '\n';
  oss << "event_trigger(" << event_fraction << ") = " << dump.event_trigger(event_fraction) << '\n';
  oss << "event_trigger(0) = " << dump.event_trigger(0.0) << '\n';

  if (filepath.has_parent_path())
  {
    std::filesystem::create_directories(filepath.parent_path());
  }
  std::ofstream out(filepath);
  out << oss.str();
}

bool shouldDumpFlipPolynomialsForEvent(
  const KineticDelaunay& kd, double occurrence_time, size_t half_edge_id, double time_margin)
{
  const auto& target_time = kd.getFlipPolynomialDumpTargetTime();
  const auto& target_half_edge = kd.getFlipPolynomialDumpTargetHalfEdge();
  if (!target_time.has_value())
  {
    return false;
  }

  if (std::abs(occurrence_time - target_time.value()) > time_margin)
  {
    return false;
  }

  if (target_half_edge.has_value())
  {
    if (!kd.isDiagnosticsHalfEdgeIdValid(target_half_edge.value()))
    {
      return false;
    }
    if (half_edge_id != target_half_edge.value())
    {
      return false;
    }
  }

  return true;
}

} // namespace kinDS

#include "KineticDelaunay.hpp"
#include "HalfEdgeDelaunayGraphToSVG.hpp"
#include "KineticDelaunayCrossingEvent.hpp"
#include "KineticDelaunayFlipEvent.hpp"
#include "KineticDelaunayHelpers.hpp"
#include "KineticDelaunayRadiusEvent.hpp"
#include "KineticDelaunaySectionEvent.hpp"
#include "KineticDelaunaySubdivisionEvent.hpp"
#include "KineticDelaunaySeparationEvent.hpp"
#include "KineticDelaunayFlipEvent.hpp"
#include "KineticDelaunayFlipEventTriggerDump.hpp"
#include "KineticDelaunayEventPredicates.hpp"
#include "Polygon2D.hpp"
#include "DebugExportFormatting.hpp"
#include "VisualDebugHighlight.hpp"
#include "SegmentBuilderVisualDebug.hpp"
#include "PlaneProjector.hpp"
#include <algorithm>
#include <array>
#include <cctype>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <functional>
#include <iomanip>
#include <glm/geometric.hpp>
#include <limits>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

using namespace kinDS;

void PendingBranchSplitState::clear()
{
  by_parent_.clear();
  strand_parent_component_.clear();
}

void PendingBranchSplitState::resetStrandLookup(size_t strand_count)
{
  strand_parent_component_.assign(strand_count, no_parent_component);
}

PendingBranchSplit& PendingBranchSplitState::getOrCreate(size_t parent_component_id)
{
  const auto [it, inserted] = by_parent_.emplace(parent_component_id, PendingBranchSplit{});
  if (inserted)
  {
    it->second.parent_component_id = parent_component_id;
  }
  return it->second;
}

const PendingBranchSplit* PendingBranchSplitState::findByParent(size_t parent_component_id) const
{
  const auto it = by_parent_.find(parent_component_id);
  if (it == by_parent_.end())
  {
    return nullptr;
  }
  return &it->second;
}

void PendingBranchSplitState::registerStrandsForSplit(
  size_t parent_component_id, const std::vector<std::vector<size_t>>& new_components)
{
  for (const auto& component : new_components)
  {
    for (size_t strand_id : component)
    {
      if (strand_id >= strand_parent_component_.size())
      {
        strand_parent_component_.resize(strand_id + 1, no_parent_component);
      }
      strand_parent_component_[strand_id] = parent_component_id;
    }
  }
}

const std::vector<size_t>* PendingBranchSplitState::frozenStrandsForStrand(size_t strand_id) const
{
  if (strand_id < strand_parent_component_.size()
    && strand_parent_component_[strand_id] != no_parent_component)
  {
    const PendingBranchSplit* split = findByParent(strand_parent_component_[strand_id]);
    if (split != nullptr)
    {
      return &split->frozen_parent_strands;
    }
  }

  for (const auto& entry : by_parent_)
  {
    const std::vector<size_t>& frozen_strands = entry.second.frozen_parent_strands;
    if (std::find(frozen_strands.begin(), frozen_strands.end(), strand_id) != frozen_strands.end())
    {
      return &frozen_strands;
    }
  }

  return nullptr;
}

void PostSplitFrameTransitionState::clear()
{
  transitions_.clear();
  strand_transition_index_.clear();
}

void PostSplitFrameTransitionState::resetStrandLookup(size_t strand_count)
{
  strand_transition_index_.assign(strand_count, no_transition);
}

PostSplitFrameTransition& PostSplitFrameTransitionState::add(PostSplitFrameTransition transition)
{
  const size_t index = transitions_.size();
  transitions_.push_back(std::move(transition));
  PostSplitFrameTransition& stored = transitions_.back();
  for (size_t strand_id : stored.strand_ids)
  {
    if (strand_id >= strand_transition_index_.size())
    {
      strand_transition_index_.resize(strand_id + 1, no_transition);
    }
    strand_transition_index_[strand_id] = index;
  }
  return stored;
}

const PostSplitFrameTransition* PostSplitFrameTransitionState::findForStrand(size_t strand_id) const
{
  if (strand_id >= strand_transition_index_.size())
  {
    return nullptr;
  }
  const size_t index = strand_transition_index_[strand_id];
  if (index == no_transition || index >= transitions_.size())
  {
    return nullptr;
  }
  return &transitions_[index];
}

void PostSplitFrameTransitionState::expireBeforeHeight(size_t height)
{
  // Drop transitions whose blend section has fully finished (native from hold_end_height+1 onward).
  // A transition with hold_end_height = S+1 / blend_section = S+1 is done once height >= S+2.
  std::vector<PostSplitFrameTransition> kept;
  kept.reserve(transitions_.size());
  for (PostSplitFrameTransition& transition : transitions_)
  {
    if (height < transition.hold_end_height + 1)
    {
      kept.push_back(std::move(transition));
    }
  }
  transitions_ = std::move(kept);
  strand_transition_index_.assign(strand_transition_index_.size(), no_transition);
  for (size_t index = 0; index < transitions_.size(); ++index)
  {
    for (size_t strand_id : transitions_[index].strand_ids)
    {
      if (strand_id >= strand_transition_index_.size())
      {
        strand_transition_index_.resize(strand_id + 1, no_transition);
      }
      strand_transition_index_[strand_id] = index;
    }
  }
}

namespace
{
constexpr double voronoi_degenerate_circumcenter_eps = 1e-12;
constexpr double voronoi_infinity_clamp_length = 1000.0;

struct VoronoiVertexPosition2D
{
  glm::dvec2 position_or_direction;
  bool infinite = false;
};

glm::dvec2 normalizeOrFallback(const glm::dvec2& v)
{
  const double len2 = glm::dot(v, v);
  if (len2 <= 1e-24)
  {
    return glm::dvec2(1.0, 0.0);
  }
  return v / std::sqrt(len2);
}

VoronoiVertexPosition2D computeVoronoiVertexPosition2D(const KineticDelaunay& kd,
  const HalfEdgeDelaunayGraph& graph, size_t face_id, double t,
  const std::function<glm::dvec2(int)>& vertex_at);

VoronoiVertexPosition2D computeVoronoiVertexPosition2D(
  const KineticDelaunay& kd, const HalfEdgeDelaunayGraph& graph, size_t face_id, double t)
{
  const auto vertex_at = [&](int vertex_index) { return kd.getPointAt(static_cast<size_t>(vertex_index), t); };
  return computeVoronoiVertexPosition2D(kd, graph, face_id, t, vertex_at);
}

VoronoiVertexPosition2D computeVoronoiVertexPosition2D(const KineticDelaunay& kd,
  const HalfEdgeDelaunayGraph& graph, size_t face_id, double t,
  const std::function<glm::dvec2(int)>& vertex_at)
{
  if (const std::optional<glm::dvec2> infinite_dir = graph.infiniteVoronoiRayDirection(face_id, vertex_at))
  {
    return { normalizeOrFallback(*infinite_dir), true };
  }

  const auto& face = graph.face(face_id);
  const std::array<int, 3> cycle_vertices = graph.adjacentTriangleVertices(face.half_edges[0]);
  const glm::dvec2 v0 = vertex_at(cycle_vertices[0]);
  const glm::dvec2 v1 = vertex_at(cycle_vertices[1]);
  const glm::dvec2 v2 = vertex_at(cycle_vertices[2]);

  const glm::dvec2 e01 = v1 - v0;
  const glm::dvec2 e12 = v2 - v1;
  const glm::dvec2 e20 = v0 - v2;
  const double len01_2 = glm::dot(e01, e01);
  const double len12_2 = glm::dot(e12, e12);
  const double len20_2 = glm::dot(e20, e20);
  const double max_edge_len2 = std::max({ len01_2, len12_2, len20_2 });
  const double area2 = std::abs(e01.x * (v2.y - v0.y) - e01.y * (v2.x - v0.x));

  if (area2 <= voronoi_degenerate_circumcenter_eps * max_edge_len2)
  {
    glm::dvec2 longest_edge_vec = e01;
    if (len12_2 >= len01_2 && len12_2 >= len20_2)
    {
      longest_edge_vec = e12;
    }
    else if (len20_2 >= len01_2 && len20_2 >= len12_2)
    {
      longest_edge_vec = e20;
    }

    glm::dvec2 circumcenter_dir { longest_edge_vec[1], -longest_edge_vec[0] };
    return { normalizeOrFallback(-circumcenter_dir), true };
  }

  return { HalfEdgeDelaunayGraph::circumcenter(v0, v1, v2), false };
}

glm::dvec2 finiteFaceAnchor(const KineticDelaunay& kd, const HalfEdgeDelaunayGraph& graph, size_t face_id, double t)
{
  glm::dvec2 anchor(0.0);
  size_t count = 0;
  for (size_t he_id : graph.face(face_id).half_edges)
  {
    const int origin = graph.halfEdge(he_id).origin;
    if (origin >= 0)
    {
      anchor += kd.getPointAt(static_cast<size_t>(origin), t);
      ++count;
    }
  }
  return count == 0 ? anchor : anchor / static_cast<double>(count);
}

/// Delaunay-edge axis used by @ref KineticDelaunay::delaunayVoronoiEdgeIntersection (even half-edge / dParam).
struct DelaunayEdgeAxis2D
{
  glm::dvec2 origin {};
  glm::dvec2 direction {}; // unnormalized; parameter t satisfies origin + t * direction
};

DelaunayEdgeAxis2D delaunayEdgeAxisForDParam(const KineticDelaunay& kd, const HalfEdgeDelaunayGraph& graph,
  size_t delaunay_edge_id, double t)
{
  const size_t d_he = 2 * delaunay_edge_id;
  if (graph.isInfinite(d_he))
  {
    const auto ray = kd.computeAngularBisector(d_he, t);
    return { ray.first, ray.second };
  }

  const glm::dvec2 edge_start = kd.getPointAt(static_cast<size_t>(graph.halfEdge(d_he).origin), t);
  const glm::dvec2 edge_end = kd.getPointAt(static_cast<size_t>(graph.destination(d_he)), t);
  return { edge_start, edge_end - edge_start };
}

/// Orthogonal projection parameter of @p point onto the Delaunay edge axis (same t as dParam).
double projectPointOntoDelaunayEdgeAxis(const DelaunayEdgeAxis2D& axis, const glm::dvec2& point)
{
  const double len2 = glm::dot(axis.direction, axis.direction);
  if (len2 < 1e-24)
  {
    return 0.0;
  }
  return glm::dot(point - axis.origin, axis.direction) / len2;
}

/// Other Voronoi-edge endpoint away from the crossing vertex. Infinite edges: a finite sample in the ray direction.
glm::dvec2 voronoiEdgeOtherEndpointForDelaunayOrder(const KineticDelaunay& kd, const HalfEdgeDelaunayGraph& graph,
  size_t voronoi_edge_id, size_t crossing_vv, const glm::dvec2& crossing_xy, double t)
{
  const size_t he_even = 2 * voronoi_edge_id;
  const int face_even = graph.halfEdge(he_even).face;
  const int face_odd = graph.halfEdge(he_even + 1).face;

  size_t other_face_id = static_cast<size_t>(-1);
  if (face_even >= 0 && static_cast<size_t>(face_even) == crossing_vv)
  {
    if (face_odd < 0)
    {
      return glm::dvec2(std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN());
    }
    other_face_id = static_cast<size_t>(face_odd);
  }
  else if (face_odd >= 0 && static_cast<size_t>(face_odd) == crossing_vv)
  {
    if (face_even < 0)
    {
      return glm::dvec2(std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN());
    }
    other_face_id = static_cast<size_t>(face_even);
  }
  else
  {
    std::ostringstream oss;
    oss << "voronoiEdgeOtherEndpointForDelaunayOrder: crossing VV " << crossing_vv
        << " is not an endpoint of Voronoi edge " << voronoi_edge_id << " (faces " << face_even << ", " << face_odd
        << ")";
    throw std::runtime_error(oss.str());
  }

  const auto vertex_at = [&](int vertex_index) { return kd.getPointAt(static_cast<size_t>(vertex_index), t); };
  const VoronoiVertexPosition2D other_vv
    = computeVoronoiVertexPosition2D(kd, graph, other_face_id, t, vertex_at);
  if (!other_vv.infinite)
  {
    return other_vv.position_or_direction;
  }

  // Infinite Voronoi edge: sample a point along the infinite direction from the crossing endpoint.
  return crossing_xy + normalizeOrFallback(other_vv.position_or_direction) * voronoi_infinity_clamp_length;
}

double referenceBranchLookupTimeForSection(size_t section, double schedule_time)
{
  const double section_start = static_cast<double>(section);
  double lookup_time = schedule_time;
  // getReferenceBranch uses floor(t) at whole numbers and floor(t)+1 for fractional t.
  // Flip/radius roots in (section, section+1] follow the fractional rule even when scheduled at section_start.
  if (lookup_time <= section_start + std::numeric_limits<double>::epsilon())
  {
    lookup_time = section_start + std::numeric_limits<double>::epsilon();
  }
  return lookup_time;
}

size_t branchLookupHeightForTime(const StrandTree& tree, double t)
{
  const size_t tree_height = tree.getHeight();
  const size_t lower_index = static_cast<size_t>(std::floor(t));
  const double frac = t - static_cast<double>(lower_index);
  size_t branch_lookup_height = lower_index;
  if (frac > std::numeric_limits<double>::epsilon())
  {
    branch_lookup_height = lower_index + 1;
  }
  if (tree_height > 0 && branch_lookup_height >= tree_height)
  {
    branch_lookup_height = tree_height - 1;
  }
  return branch_lookup_height;
}
} // namespace

glm::dvec3 kinDS::KineticDelaunay::computeVoronoiVertexHomogenous(
  size_t voronoi_vertex_id, double t, bool apply_reference_transform, bool include_virtual_offset) const
{
  requireLiveRegisteredVoronoiVertex(voronoi_vertex_id, "computeVoronoiVertexHomogenous");
  const auto vertex_at = [&](int vertex_index) {
    return getPointAt(static_cast<size_t>(vertex_index), t, apply_reference_transform, include_virtual_offset);
  };
  if (const std::optional<glm::dvec2> infinite_dir = graph.infiniteVoronoiRayDirection(voronoi_vertex_id, vertex_at))
  {
    return glm::dvec3(normalizeOrFallback(*infinite_dir), 0.0);
  }

  auto face = graph.face(voronoi_vertex_id);
  auto vertices = graph.adjacentTriangleVertices(face.half_edges[0]);
  const glm::dvec2 p0 = vertex_at(vertices[0]);
  const glm::dvec2 p1 = vertex_at(vertices[1]);
  const glm::dvec2 p2 = vertex_at(vertices[2]);
  auto circumcenter = graph.circumcenter(p0, p1, p2);
  return glm::dvec3(circumcenter, 1.0);
}

glm::dvec3 KineticDelaunay::computeVoronoiVertexClampedInfinity(
  size_t half_edge_id, double t, bool apply_reference_transform, bool include_virtual_offset) const
{
  const auto& graph = getGraph();
  const int face_id = graph.halfEdge(half_edge_id).face;
  if (face_id < 0)
  {
    return glm::dvec3(0.0, 0.0, t);
  }

  size_t reference_strand = 0;
  for (size_t he_id : graph.face(static_cast<size_t>(face_id)).half_edges)
  {
    const int origin = graph.halfEdge(he_id).origin;
    if (origin >= 0)
    {
      reference_strand = static_cast<size_t>(origin);
      break;
    }
  }

  return computeVoronoiVertexClampedInfinityWithReferenceBranch(
    half_edge_id, t, getReferenceBranch(reference_strand, t), apply_reference_transform, include_virtual_offset);
}

glm::dvec3 KineticDelaunay::computeVoronoiVertexClampedInfinityWithReferenceBranch(
  size_t half_edge_id, double t, size_t reference_branch, bool apply_reference_transform,
  bool include_virtual_offset) const
{
  const auto& graph = getGraph();
  const int face_id = graph.halfEdge(half_edge_id).face;
  if (face_id < 0)
  {
    return glm::dvec3(0.0, 0.0, t);
  }

  requireLiveRegisteredVoronoiVertex(static_cast<size_t>(face_id), "computeVoronoiVertexClampedInfinity");

  const auto vertex_at = [&](int vertex_index) {
    return getPointAtWithReferenceBranch(
      static_cast<size_t>(vertex_index), t, reference_branch, apply_reference_transform, include_virtual_offset);
  };

  const VoronoiVertexPosition2D endpoint
    = computeVoronoiVertexPosition2D(*this, graph, static_cast<size_t>(face_id), t, vertex_at);
  if (!endpoint.infinite)
  {
    return glm::dvec3(endpoint.position_or_direction, t);
  }

  const int twin_face_id = graph.halfEdge(half_edge_id ^ 1).face;
  glm::dvec2 anchor(0.0);
  size_t anchor_count = 0;
  for (size_t he_id : graph.face(static_cast<size_t>(face_id)).half_edges)
  {
    const int origin = graph.halfEdge(he_id).origin;
    if (origin >= 0)
    {
      anchor += vertex_at(origin);
      ++anchor_count;
    }
  }
  if (anchor_count > 0)
  {
    anchor /= static_cast<double>(anchor_count);
  }

  if (twin_face_id >= 0)
  {
    const VoronoiVertexPosition2D twin_endpoint
      = computeVoronoiVertexPosition2D(*this, graph, static_cast<size_t>(twin_face_id), t, vertex_at);
    if (!twin_endpoint.infinite)
    {
      anchor = twin_endpoint.position_or_direction;
    }
  }

  const glm::dvec2 clamped = anchor + voronoi_infinity_clamp_length * endpoint.position_or_direction;
  return glm::dvec3(clamped, t);
}

size_t KineticDelaunay::getCrossingDataContainingTriId(size_t voronoi_vertex_id) const
{
  requireLiveRegisteredVoronoiVertex(voronoi_vertex_id, "getCrossingDataContainingTriId");
  return crossing_data.getContainingTriId(voronoi_vertex_id);
}

bool KineticDelaunay::isCrossingDataVoronoiVertexRegistered(size_t voronoi_vertex_id) const
{
  return crossing_data.isVoronoiVertexRegistered(voronoi_vertex_id);
}

void KineticDelaunay::requireLiveRegisteredVoronoiVertex(size_t voronoi_vertex_id, const char* context) const
{
  const std::string ctx
    = (context != nullptr && context[0] != '\0') ? context : "requireLiveRegisteredVoronoiVertex";

  if (voronoi_vertex_id >= graph.faceSlotCount())
  {
    throw std::runtime_error(ctx + ": Voronoi vertex " + std::to_string(voronoi_vertex_id) + " is out of range (faceSlotCount="
      + std::to_string(graph.faceSlotCount()) + ")");
  }
  if (!graph.isLiveFace(voronoi_vertex_id))
  {
    throw std::runtime_error(ctx + ": Voronoi vertex " + std::to_string(voronoi_vertex_id)
      + " is dead (face slot is no longer live after graph update)");
  }
  if (!crossing_data.isVoronoiVertexRegistered(voronoi_vertex_id))
  {
    throw std::runtime_error(ctx + ": Voronoi vertex " + std::to_string(voronoi_vertex_id)
      + " is not registered in CrossingData");
  }
}

std::vector<size_t> KineticDelaunay::getCrossingDataVoronoiVerticesInTri(size_t tri_id) const
{
  return crossing_data.getVoronoiVerticesInTri(tri_id);
}

glm::dvec3 KineticDelaunay::getVoronoiVertexHomogeneous(
  size_t voronoi_vertex_id, double t, bool apply_reference_transform, bool include_virtual_offset) const
{
  return computeVoronoiVertexHomogenous(voronoi_vertex_id, t, apply_reference_transform, include_virtual_offset);
}

std::vector<double> KineticDelaunay::findEvents(
  Polynomial& event_trigger, double min_fraction, bool only_positive_to_negative)
{
  if (event_trigger.degree() == -1)
  {
    // No events possible, return empty vector
    return {};
  }

  event_trigger.trim();
  auto zeros = event_trigger.realRoots();
  std::vector<double> filtered_sorted_zeros;

  for (const auto& root : zeros)
  {
    if (isnan(root))
    {
      continue; // Skip NaN roots
    }

    if (root > min_fraction && root <= kEventIntervalFractionUpperBound)
    {
      double event_time = root;
      filtered_sorted_zeros.emplace_back(event_time);
    }
  }

  if (filtered_sorted_zeros.empty())
  {
    // No valid events found, return empty vector
    return {};
  }

  // Sort events ascending by time
  std::sort(filtered_sorted_zeros.begin(), filtered_sorted_zeros.end());

  // Determine sign changes for each root
  std::vector<double> interval_signs(filtered_sorted_zeros.size() + 1);
  double test_point = (min_fraction + filtered_sorted_zeros[0]) / 2.0; // Start with a test point before the first root
  interval_signs[0] = event_trigger(test_point) > 0 ? 1 : -1;

  for (size_t i = 0; i < filtered_sorted_zeros.size(); ++i)
  {
    test_point
      = (filtered_sorted_zeros[i] + (i + 1 < filtered_sorted_zeros.size() ? filtered_sorted_zeros[i + 1] : 1.0)) / 2.0;
    interval_signs[i + 1] = event_trigger(test_point) > 0 ? 1 : -1;
  }

  std::vector<double> found_event_times;

  for (size_t i = 0; i < filtered_sorted_zeros.size(); ++i)
  {
    if (interval_signs[i] != interval_signs[i + 1])
    {
      if (only_positive_to_negative && interval_signs[i] < 0)
      {
        KINDS_DEBUG("Sign change from negative to positive at root "
          << filtered_sorted_zeros[i] << ", skipping event creation due to only_positive_to_negative flag.");
        continue; // Skip if we only want positive to negative sign changes
      }
      else
      {
        // Sign change detected, create an event
        double event_time = filtered_sorted_zeros[i];
        found_event_times.push_back(event_time);
        KINDS_DEBUG("Event will be queued at time " << event_time << " with sign change from " << interval_signs[i]
                                                    << " to " << interval_signs[i + 1]);
      }
    }
    else
    {
      KINDS_DEBUG("No sign change at root " << filtered_sorted_zeros[i] << ", skipping event creation.");
    }
  }

  return found_event_times;
}

std::vector<double> KineticDelaunay::findVirtualEvents(
  Polynomial& event_trigger, double min_x, bool only_positive_to_negative)
{
  if (event_trigger.degree() == -1)
  {
    return {};
  }

  event_trigger.trim();
  auto zeros = event_trigger.realRoots();
  std::vector<double> filtered_sorted_zeros;

  for (const auto& root : zeros)
  {
    if (std::isnan(root) || !std::isfinite(root))
    {
      continue;
    }
    if (root > min_x)
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
  double test_point = (min_x + filtered_sorted_zeros[0]) * 0.5;
  if (!(test_point > min_x))
  {
    test_point = filtered_sorted_zeros[0] - std::max(1e-9, std::abs(filtered_sorted_zeros[0]) * 1e-12);
  }
  interval_signs[0] = event_trigger(test_point) > 0.0 ? 1.0 : -1.0;

  for (size_t i = 0; i < filtered_sorted_zeros.size(); ++i)
  {
    const double next_bound = (i + 1 < filtered_sorted_zeros.size())
      ? filtered_sorted_zeros[i + 1]
      : (filtered_sorted_zeros[i] + 1.0);
    test_point = (filtered_sorted_zeros[i] + next_bound) * 0.5;
    interval_signs[i + 1] = event_trigger(test_point) > 0.0 ? 1.0 : -1.0;
  }

  std::vector<double> found_event_times;
  for (size_t i = 0; i < filtered_sorted_zeros.size(); ++i)
  {
    if (interval_signs[i] == interval_signs[i + 1])
    {
      continue;
    }
    if (only_positive_to_negative && interval_signs[i] < 0.0)
    {
      continue;
    }
    found_event_times.push_back(filtered_sorted_zeros[i]);
  }
  return found_event_times;
}

void KineticDelaunay::reassignVoronoiVerticesOnBoundary(size_t he_id, double t,
  std::optional<InfinitesimalComputeContext> infinitesimal)
{
  // Find which side of he_id is on the outside of the convex hull
  size_t face0 = graph.halfEdge(he_id).face;
  size_t face1 = graph.halfEdge(he_id ^ 1).face;

  // If either face has an infinite vertex, that face is the outside
  bool face0_is_outside = false;
  bool face1_is_outside = false;

  // A face is outside if any of its triangle vertices is -1
  auto isFaceOutside = [&](size_t face_id) -> bool
  {
    const auto& face = graph.face(face_id);
    for (size_t i = 0; i < 3; ++i)
    {
      if (graph.halfEdge(face.half_edges[i]).origin == size_t(-1))
      {
        return true;
      }
    }
    return false;
  };

  face0_is_outside = isFaceOutside(face0);
  face1_is_outside = isFaceOutside(face1);

  // Determine the outside face id
  size_t outside_face_id = face0_is_outside ? face0 : face1;

  // For each voronoi vertex in face0, reassign to outside face if not already there
  auto vertices0 = crossing_data.getVoronoiVerticesInTri(face0);
  for (size_t voronoi_vertex : vertices0)
  {
    if (crossing_data.getContainingTriId(voronoi_vertex) != outside_face_id)
    {
      crossing_data.moveVertex(voronoi_vertex, outside_face_id, t);
    }
  }

  // For each voronoi vertex in face1, reassign to outside face if not already there
  auto vertices1 = crossing_data.getVoronoiVerticesInTri(face1);
  for (size_t voronoi_vertex : vertices1)
  {
    if (crossing_data.getContainingTriId(voronoi_vertex) != outside_face_id)
    {
      crossing_data.moveVertex(voronoi_vertex, outside_face_id, t);
    }
  }

  // Secondly, we need to add all intersections to the new Delaunay edge (he_id/2).
  // To do this, find the "inner" edge (the one whose face is not outside) and get its triangle's two other half-edges.
  // Those two edges will have Voronoi–Delaunay intersections; we copy those intersections over to the new (now
  // boundary) edge, adjusting appropriately.

  auto& boundary_edge_d_intersections = crossing_data.delaunay_edge_intersections[he_id / 2];
  auto it = boundary_edge_d_intersections.begin();
  while (it != boundary_edge_d_intersections.end())
  {
    // Capture next BEFORE erasing, because removeIntersection() erases from this list.
    auto next = std::next(it);
    crossing_data.removeIntersection(*it);
    it = next;
  }

  size_t inner_face_id = face0_is_outside ? face1 : face0;
  size_t inner_he_id = (graph.halfEdge(he_id).face == inner_face_id) ? he_id : (he_id ^ 1);

  // Get the three half-edges for this triangle, and identify the other two that aren't inner_he_id
  const auto& inner_face = graph.face(inner_face_id);
  std::vector<size_t> triangle_half_edges;
  size_t next_he_id = graph.halfEdge(inner_he_id).next;
  triangle_half_edges.push_back(next_he_id);
  next_he_id = graph.halfEdge(next_he_id).next;
  triangle_half_edges.push_back(next_he_id);

  // Now, for each intersection on these two edges, copy the intersection (with changes) to he_id/2 boundary edge
  // for (size_t tri_he : triangle_half_edges) {
  auto copy_intersections = [&](size_t tri_he, bool backwards)
  {
    auto& d_intersections = crossing_data.delaunay_edge_intersections[tri_he / 2];
    bool even = tri_he % 2 == 0;

    auto process_intersection = [&](CrossingData::VoronoiDelaunayEdgeIntersection intersection)
    {
      intersection.delaunay_edge_id = inner_he_id / 2;
      assignCrossingIntersectionDelaunayParam(this, intersection,
        delaunayVoronoiEdgeIntersection(intersection.delaunay_edge_id, intersection.voronoi_edge_id, t).first, t,
        "copyBoundaryIntersections");
      // Copied intersections move to a new boundary Delaunay edge; mesh-pair links from the source edge are invalid.
      intersection.prev_segment_mesh_pair_index = static_cast<size_t>(-1);
      intersection.next_segment_mesh_pair_index = static_cast<size_t>(-1);

      auto& v_intersections = crossing_data.voronoi_edge_intersections[intersection.voronoi_edge_id];

      crossing_data.edge_intersections.emplace_back(intersection);
      auto intersection_it = std::prev(crossing_data.edge_intersections.end());

      // check if the start or end voronoi vertex lies in the outer triangle to determine on which side to insert
      size_t start_voronoi_vertex_id = graph.halfEdge(2 * intersection.voronoi_edge_id).face;
      size_t end_voronoi_vertex_id = graph.halfEdge(2 * intersection.voronoi_edge_id + 1).face;
      size_t start_containing_triangle_id = crossing_data.getContainingTriId(start_voronoi_vertex_id);
      size_t end_containing_triangle_id = crossing_data.getContainingTriId(end_voronoi_vertex_id);
      std::list<CrossingData::EdgeIntersectionRef>::iterator v_ref;
      if (start_containing_triangle_id == outside_face_id)
      {
        v_ref = v_intersections.insert(v_intersections.begin(), intersection_it);
      }
      else if (end_containing_triangle_id == outside_face_id)
      {
        v_ref = v_intersections.insert(v_intersections.end(), intersection_it);
      }
      else
      {
        KINDS_DEBUG("Skipping intersection copy during flip: Voronoi edge " << intersection.voronoi_edge_id
          << " is not on the boundary of outside face " << outside_face_id);
        crossing_data.edge_intersections.pop_back();
        return;
      }

      auto d_ref = boundary_edge_d_intersections.insert(boundary_edge_d_intersections.end(), intersection_it);
      intersection_it->voronoi_ref = v_ref;
      intersection_it->delaunay_ref = d_ref;
    };

    if (even != backwards)
    {
      for (auto iter : d_intersections)
      {
        process_intersection(*iter);
      }
    }
    else
    {
      for (auto iter = d_intersections.rbegin(); iter != d_intersections.rend(); iter++)
      {
        process_intersection(**iter);
      }
    }
  };

  if (inner_he_id % 2 == 0)
  {
    copy_intersections(triangle_half_edges[1], true);
    copy_intersections(triangle_half_edges[0], true);
  }
  else
  {
    copy_intersections(triangle_half_edges[0], false);
    copy_intersections(triangle_half_edges[1], false);
  }

  vertices0 = crossing_data.getVoronoiVerticesInTri(face0);
  vertices1 = crossing_data.getVoronoiVerticesInTri(face1);
  std::unordered_set<size_t> recompute_voronoi_vertices;
  for (size_t voronoi_vertex : vertices0)
  {
    recompute_voronoi_vertices.insert(voronoi_vertex);
  }
  for (size_t voronoi_vertex : vertices1)
  {
    recompute_voronoi_vertices.insert(voronoi_vertex);
  }
  recompute_voronoi_vertices.insert(face0);
  recompute_voronoi_vertices.insert(face1);
  for (size_t voronoi_vertex : recompute_voronoi_vertices)
  {
    crossing_event_manager_->computeEvents(t, voronoi_vertex, infinitesimal);
    if (voronoi_vertex < crossing_data.last_crossing.size())
    {
      crossing_data.last_crossing[voronoi_vertex]
        = infinitesimal.has_value() ? EventTime(t, infinitesimal->min_infinitesimal_t) : EventTime(t);
    }
  }
}

void KineticDelaunay::reassignVoronoiVerticesInQuadrilateral(size_t quad_index, double t,
  const std::map<size_t, size_t>& pre_flip_quad_faces, std::optional<InfinitesimalComputeContext> infinitesimal)
{
  size_t he_id = quad_index * 2;
  size_t face_id0 = graph.halfEdge(he_id).face;
  size_t face_id1 = graph.halfEdge(he_id ^ 1).face;

  // After a CCW (or similar) hull flip the diagonal can become an infinite convex-hull edge. Crossings on that
  // Delaunay edge are meaningless; discard them and do not recompute params onto it.
  const bool flip_edge_infinite = graph.isInfinite(he_id);

  size_t u = graph.halfEdge(he_id).origin;
  size_t v = graph.destination(he_id);

  glm::dvec2 edge_vector;
  glm::dvec2 pu, pv;
  // if a vertex is infinite, we need to compute a placeholder vector. Just as for the predicates we choose the vector
  // perpendicular to the line through the two neighboring vertices on the convex hull.
  if (u == -1 || v == -1)
  {
    size_t opposite0 = graph.triangleOppositeVertex(he_id);
    size_t opposite1 = graph.triangleOppositeVertex(he_id ^ 1);

    glm::vec2 p0 = getPointAt(opposite0, t);
    glm::vec2 p1 = getPointAt(opposite1, t);
    edge_vector = glm::normalize(glm::vec2(p1.y - p0.y, p0.x - p1.x));

    if (u == -1)
    {
      pv = getPointAt(v, t);
      pu = pv - edge_vector; // placeholder position for the infinite vertex
    }
    else
    {
      pu = getPointAt(u, t);
      pv = pu + edge_vector; // placeholder position for the infinite vertex
    }
  }
  else
  {
    pu = getPointAt(u, t);
    pv = getPointAt(v, t);

    edge_vector = pv - pu;
  }

  auto vertices0 = crossing_data.getVoronoiVerticesInTri(face_id0);
  auto vertices1 = crossing_data.getVoronoiVerticesInTri(face_id1);

  // TODO: what if the edge function evaluates to 0. Perhaps we should look at the derivative.
  auto reassign_vertices = [&](const std::vector<size_t>& vertices, size_t target_face_id, bool target_positive_side)
  {
    for (size_t voronoi_vertex : vertices)
    {
      // Compute the position of the Voronoi vertex at time t and check on which side of the edge it is.
      glm::dvec3 voronoi_pos = computeVoronoiVertexHomogenous(voronoi_vertex, t);
      if (voronoi_pos.z != 0)
      {
        glm::dvec2 event_pos = glm::dvec2(voronoi_pos.x / voronoi_pos.z, voronoi_pos.y / voronoi_pos.z);
        const double side = glm::cross(edge_vector, event_pos - pu);
        const bool is_on_target_side = target_positive_side ? (side > 0.0) : (side < 0.0);
        if (is_on_target_side)
        {
          crossing_data.moveVertex(voronoi_vertex, target_face_id, t);
        }
      }
      else
      {
        // Must belong to the dual triangle
        crossing_data.moveVertex(voronoi_vertex, voronoi_vertex, t);
      }
    }
  };

  reassign_vertices(vertices0, face_id1, true);
  reassign_vertices(vertices1, face_id0, false);

  // get updated vertices
  vertices0 = crossing_data.getVoronoiVerticesInTri(face_id0);
  vertices1 = crossing_data.getVoronoiVerticesInTri(face_id1);

  if (flip_edge_infinite)
  {
    KINDS_DEBUG("reassignVoronoiVerticesInQuadrilateral: flip edge de=" << quad_index
                                                                        << " is infinite after flip at t=" << t
                                                                        << "; discarding all intersections on that edge");
    auto& d_list = crossing_data.delaunay_edge_intersections[quad_index];
    while (!d_list.empty())
    {
      crossing_data.removeIntersection(d_list.front());
    }
  }

  auto update_intersections = [&](const std::vector<size_t>& vertices)
  {
    // Find and handle all voronoi edges that lie entirely within the quadrilateral
    for (size_t voronoi_vertex : vertices)
    {
      size_t containing_tri_id = crossing_data.getContainingTriId(voronoi_vertex);
      // Detect for all Voronoi edges incident to the vertex if they cross the newly flipped edge
      auto adjacent_voronoi_edges = graph.face(voronoi_vertex).half_edges;
      for (size_t adj_he_id : adjacent_voronoi_edges)
      {
        size_t other_voronoi_vertex = graph.halfEdge(adj_he_id ^ 1).face;

        size_t other_containing_tri_id = crossing_data.getContainingTriId(other_voronoi_vertex);
        if (other_containing_tri_id == face_id0 || other_containing_tri_id == face_id1)
        {
          // fully within the quadrilateral, might have one intersection removed or added with the flipped edge
          auto v_intersections = crossing_data.voronoi_edge_intersections[adj_he_id / 2];
          if (containing_tri_id == other_containing_tri_id)
          {
            // no intersection exists, loop through intersections and remove them

            for (auto intersection : v_intersections)
            {
              crossing_data.removeIntersection(intersection);
            }
          }
          else if (flip_edge_infinite)
          {
            // Infinite flip edge: drop any leftover crossings on this Voronoi edge that pointed at the flip edge.
            for (auto intersection : v_intersections)
            {
              if (intersection->delaunay_edge_id == quad_index)
              {
                crossing_data.removeIntersection(intersection);
              }
            }
          }
          else
          {
            // there is an intersection, either update the existing intersection or add a new one
            // Compute the parameter along the delaunay edge

            const size_t voronoi_edge_id = adj_he_id / 2;
            if (v_intersections.empty())
            {
              // New intersection with the flipped (quad) Delaunay edge; sole entry on this Voronoi edge.
              crossing_data.edge_intersections.emplace_back();
              const CrossingData::EdgeIntersectionRef intersection = std::prev(crossing_data.edge_intersections.end());
              intersection->delaunay_edge_id = quad_index;
              intersection->voronoi_edge_id = voronoi_edge_id;
              intersection->prev_segment_mesh_pair_index = static_cast<size_t>(-1);
              intersection->next_segment_mesh_pair_index = static_cast<size_t>(-1);
              assignCrossingIntersectionDelaunayParam(this, *intersection,
                delaunayVoronoiEdgeIntersection(quad_index, voronoi_edge_id, t).first, t,
                "flipUpdate:newVoronoiIntersection", voronoi_vertex, other_voronoi_vertex);

              auto& v_list = crossing_data.voronoi_edge_intersections[voronoi_edge_id];
              intersection->voronoi_ref = v_list.emplace(v_list.end(), intersection);

              auto& d_list = crossing_data.delaunay_edge_intersections[quad_index];
              intersection->delaunay_ref = d_list.emplace(d_list.end(), intersection);
            }
            else
            {
              // Update the existing sole intersection (parameter along the flip edge).
              const CrossingData::EdgeIntersectionRef intersection = v_intersections.front();
              const size_t old_delaunay_edge_id = intersection->delaunay_edge_id;
              intersection->delaunay_edge_id = quad_index;
              assignCrossingIntersectionDelaunayParam(this, *intersection,
                delaunayVoronoiEdgeIntersection(quad_index, voronoi_edge_id, t).first, t,
                "flipUpdate:existingVoronoiIntersection", voronoi_vertex, other_voronoi_vertex);

              if (old_delaunay_edge_id != quad_index)
              {
                // I don't think this case can occurr, but the AI suggested it, so I'm leaving it in.
                KINDS_DEBUG("Old Delaunay edge id: " << old_delaunay_edge_id
                                                     << " is not the same as the new quad index: " << quad_index);
                if (old_delaunay_edge_id < crossing_data.delaunay_edge_intersections.size())
                {
                  if (intersection->delaunay_ref.has_value())
                  {
                    crossing_data.delaunay_edge_intersections[old_delaunay_edge_id].erase(*intersection->delaunay_ref);
                  }
                }
                auto& d_list = crossing_data.delaunay_edge_intersections[quad_index];
                intersection->delaunay_ref = d_list.emplace(d_list.end(), intersection);
              }
            }
          }
        }
      }
    }
  };

  update_intersections(vertices0);
  update_intersections(vertices1);

  // Now handle all voronoi edges that are partially outside the quadrilateral
  auto quad_he_ids = graph.getQuadBoundaryHalfEdgeIndices(quad_index);

  for (size_t he_id : quad_he_ids)
  {
    KINDS_DEBUG("Processing he_id: " << he_id);
    auto& d_edge_intersections = crossing_data.delaunay_edge_intersections[he_id / 2];

    KINDS_DEBUG("Found " << d_edge_intersections.size() << " intersections on Delaunay edge.");

    for (auto d_it = d_edge_intersections.begin(); d_it != d_edge_intersections.end(); )
    {
      CrossingData::EdgeIntersectionRef intersection = *d_it;
      bool stale_intersection = false;
      KINDS_DEBUG("Processing intersection with Voronoi edge: "
        << intersection->voronoi_edge_id << " and Delaunay edge: " << intersection->delaunay_edge_id);
      // get the next and previous intersections to check if they match with the flipped edge
      auto& v_edge_intersections = crossing_data.voronoi_edge_intersections[intersection->voronoi_edge_id];
      auto v_ref = intersection->voronoi_ref;
      if (!v_ref.has_value())
      {
        throw std::runtime_error("CrossingData update: intersection has unset voronoi_ref");
      }
      std::list<CrossingData::EdgeIntersectionRef>::iterator v_next;
      std::list<CrossingData::EdgeIntersectionRef>::iterator v_prev;

      // check if next or previous matches the face
      size_t face_inside_old = pre_flip_quad_faces.at(he_id);
      size_t face_inside_new = graph.halfEdge(he_id).face;
      KINDS_DEBUG("face_inside_old: " << face_inside_old << ", face_inside_new: " << face_inside_new);

      v_next = std::next(*v_ref);
      if (*v_ref == v_edge_intersections.begin())
      {
        v_prev = v_edge_intersections.end();
      }
      else
      {
        v_prev = std::prev(*v_ref);
      }

      bool use_prev = false;
      bool use_next = false;
      bool at_end = false;
      bool intersected_before = false;
      std::list<CrossingData::EdgeIntersectionRef>::iterator v_intersection = v_edge_intersections.end();

      if (v_prev != v_edge_intersections.end())
      {
        KINDS_DEBUG("Found previous intersection with Voronoi edge: "
          << (*v_prev)->voronoi_edge_id << " and Delaunay edge: " << (*v_prev)->delaunay_edge_id);
        size_t d_edge_id = (*v_prev)->delaunay_edge_id;
        size_t prev_face0 = graph.halfEdge(2 * d_edge_id).face;
        size_t prev_face1 = graph.halfEdge(2 * d_edge_id + 1).face;
        KINDS_DEBUG("Faces of Delaunay edge at previous intersection: " << prev_face0 << ", " << prev_face1);

        // need to compare to both face ids of the quad because the flip might have caused a mismatch between two
        // neighboring intersections
        use_prev = (face_id0 == prev_face0) || (face_id0 == prev_face1) || (face_id1 == prev_face0)
          || (face_id1 == prev_face1);

        if (use_prev)
        {
          // check if it intersected before
          if ((prev_face0 == face_id0 && prev_face1 == face_id1) || (prev_face0 == face_id1 && prev_face1 == face_id0))
          {
            intersected_before = true;
            v_intersection = v_prev;

            if (v_prev == v_edge_intersections.begin())
            {
              v_prev = v_edge_intersections.end();
            }
            else
            {
              v_prev = std::prev(v_prev);
            }

            KINDS_DEBUG("Intersection existed before, moving v_prev to "
              << (v_prev != v_edge_intersections.end() ? "previous" : "end"));
          }
        }
      }
      else
      {

        size_t start_voronoi_vertex_id = graph.halfEdge(2 * intersection->voronoi_edge_id).face;
        size_t containing_triangle_id = crossing_data.getContainingTriId(start_voronoi_vertex_id);
        KINDS_DEBUG("start_voronoi_vertex_id: " << start_voronoi_vertex_id
                                               << ", containing_triangle_id: " << containing_triangle_id
                                               << ", face_id0: " << face_id0 << ", face_id1: " << face_id1);

        use_prev = containing_triangle_id == face_id0 || containing_triangle_id == face_id1;
        if (use_prev)
        {
          KINDS_DEBUG("No previous intersection found for Voronoi edge: "
            << intersection->voronoi_edge_id << " and Delaunay edge: " << intersection->delaunay_edge_id
            << ", now using location of endpoint");
          at_end = true;
        }
        else
        {
          KINDS_DEBUG("No previous intersection found for Voronoi edge: "
            << intersection->voronoi_edge_id << " and Delaunay edge: " << intersection->delaunay_edge_id
            << ", and endpoint is not in the quadrilateral, attempting to use next intersection");
        }
      }

      // now do the same for next
      if (!use_prev)
      {
        if (v_next != v_edge_intersections.end())
        {
          KINDS_DEBUG("Found next intersection with Voronoi edge: "
            << (*v_next)->voronoi_edge_id << " and Delaunay edge: " << (*v_next)->delaunay_edge_id);
          size_t d_edge_id = (*v_next)->delaunay_edge_id;
          size_t next_face0 = graph.halfEdge(2 * d_edge_id).face;
          size_t next_face1 = graph.halfEdge(2 * d_edge_id + 1).face;

          use_next = (face_id0 == next_face0) || (face_id0 == next_face1) || (face_id1 == next_face0)
            || (face_id1 == next_face1);
          if (use_next)
          {
            // check if it intersected before
            if ((next_face0 == face_id0 && next_face1 == face_id1)
              || (next_face0 == face_id1 && next_face1 == face_id0))
            {
              intersected_before = true;
              v_intersection = v_next;
              v_next = std::next(v_next);
              KINDS_DEBUG("Intersection existed before, moving v_next to "
                << (v_next != v_edge_intersections.end() ? "next" : "end"));
            }
          }
        }
        else
        {
          size_t end_voronoi_vertex_id = graph.halfEdge(2 * intersection->voronoi_edge_id + 1).face;
          size_t containing_triangle_id = crossing_data.getContainingTriId(end_voronoi_vertex_id);
          // output values
          KINDS_DEBUG("end_voronoi_vertex_id: " << end_voronoi_vertex_id
                                               << ", containing_triangle_id: " << containing_triangle_id
                                               << ", face_id0: " << face_id0 << ", face_id1: " << face_id1);
          use_next = containing_triangle_id == face_id0 || containing_triangle_id == face_id1;
          if (use_next)
          {
            KINDS_DEBUG("No next intersection found for Voronoi edge: "
              << intersection->voronoi_edge_id << " and Delaunay edge: " << intersection->delaunay_edge_id
              << ", now using location of endpoint");
            at_end = true;
          }
          else
          {
            KINDS_DEBUG("Stale flip intersection: Voronoi edge " << intersection->voronoi_edge_id << " Delaunay edge "
              << intersection->delaunay_edge_id
              << " has no next Voronoi-list neighbor and its far endpoint is outside the flip quadrilateral; removing");
            stale_intersection = true;
          }
        }
      }

      if (!stale_intersection && use_next == use_prev)
      {
        if (!use_next)
        {
          KINDS_DEBUG("Stale flip intersection: Voronoi edge " << intersection->voronoi_edge_id << " Delaunay edge "
            << intersection->delaunay_edge_id
            << " cannot be oriented relative to the flip quadrilateral; removing");
          stale_intersection = true;
        }
        else
        {
          KINDS_ERROR("use_next: " << use_next << ", use_prev: " << use_prev);
          throw std::runtime_error("Logic error: use_next and use_prev cannot be equal!");
        }
      }

      if (stale_intersection)
      {
        auto next_d_it = std::next(d_it);
        crossing_data.removeIntersection(intersection);
        d_it = next_d_it;
        continue;
      }

      KINDS_DEBUG("Using " << (use_prev ? "previous" : "next") << " intersection for Voronoi edge: "
                          << intersection->voronoi_edge_id << " and Delaunay edge: " << intersection->delaunay_edge_id);

      if (use_prev)
      {
        v_next = v_prev;
      }

      auto add_new_quad_intersection
        = [&](size_t delaunay_edge_id, size_t voronoi_edge_id, double delaunay_edge_param, bool insert_after,
            std::list<CrossingData::EdgeIntersectionRef>::iterator insert_position)
      {
        crossing_data.edge_intersections.emplace_back();

        CrossingData::EdgeIntersectionRef new_intersection = std::prev(crossing_data.edge_intersections.end());
        new_intersection->delaunay_edge_id = delaunay_edge_id;
        new_intersection->voronoi_edge_id = voronoi_edge_id;
        assignCrossingIntersectionDelaunayParam(
          this, *new_intersection, delaunay_edge_param, t, "flipUpdate:addNewQuadIntersection");

        if (!insert_after)
        {
          new_intersection->voronoi_ref = v_edge_intersections.insert(insert_position, new_intersection);
        }
        else
        {
          new_intersection->voronoi_ref = v_edge_intersections.insert(std::next(insert_position), new_intersection);
        }
        crossing_data.delaunay_edge_intersections[delaunay_edge_id].push_back(new_intersection); // sort later
        // In a list, an iterator will not be invalidated by sorting, so we can safely set the delaunay_ref now and
        // sort later.
        new_intersection->delaunay_ref = std::prev(crossing_data.delaunay_edge_intersections[delaunay_edge_id].end());
      };

      if (v_next != v_edge_intersections.end())
      { // Voronoi edge extends beyond quadrilateral
        auto quad_he_id_it = std::find(quad_he_ids.begin(), quad_he_ids.end(), (*v_next)->delaunay_edge_id * 2);
        if (quad_he_id_it == quad_he_ids.end())
        {
          // try again with odd id
          quad_he_id_it = std::find(quad_he_ids.begin(), quad_he_ids.end(), (*v_next)->delaunay_edge_id * 2 + 1);
        }

        if (quad_he_id_it == quad_he_ids.end())
        {
          throw std::runtime_error("Logic error: could not find matching half-edge in quadrilateral for intersection");
        }

        size_t other_he_id = *quad_he_id_it;

        // now compare the face ids to determine if there is an intersection
        size_t face_id = graph.halfEdge(he_id).face;
        size_t other_face_id = graph.halfEdge(other_he_id).face;

        if (face_id != other_face_id)
        {
          // new intersection
          if (!intersected_before)
          {
            if (flip_edge_infinite)
            {
              // Infinite flip edge carries no crossings.
            }
            else
            {
              const double delaunay_edge_param
                = delaunayVoronoiEdgeIntersection(quad_index, intersection->voronoi_edge_id, t).first;
              add_new_quad_intersection(quad_index, intersection->voronoi_edge_id, delaunay_edge_param, use_prev, v_next);
            }
          }
          else if (flip_edge_infinite)
          {
            crossing_data.removeIntersection(*v_intersection);
          }
          else
          {
            assignCrossingIntersectionDelaunayParam(this, **v_intersection,
              delaunayVoronoiEdgeIntersection(quad_index, (*v_intersection)->voronoi_edge_id, t).first, t,
              "flipUpdate:extendBeyondQuad");
          }
        }
        else
        {
          if (intersected_before)
          {
            crossing_data.removeIntersection(*v_intersection);
          }
        }
      }
      else // Voronoi edge is partially in quadrilateral, check which face its end belongs to.
      {
        size_t containing_triangle_id;
        if (use_next)
        {
          size_t end_voronoi_vertex_id = graph.halfEdge(2 * intersection->voronoi_edge_id + 1).face;
          containing_triangle_id = crossing_data.getContainingTriId(end_voronoi_vertex_id);
        }
        else
        {
          size_t start_voronoi_vertex_id = graph.halfEdge(2 * intersection->voronoi_edge_id).face;
          containing_triangle_id = crossing_data.getContainingTriId(start_voronoi_vertex_id);
        }

        if (containing_triangle_id == face_inside_new && intersected_before)
        {
          // remove intersection
          crossing_data.removeIntersection(*v_intersection);
        }
        else if (containing_triangle_id != face_inside_new && !intersected_before)
        {
          if (!flip_edge_infinite)
          {
            // add new intersection
            const double delaunay_edge_param
              = delaunayVoronoiEdgeIntersection(quad_index, intersection->voronoi_edge_id, t).first;

            if (!use_prev)
            {

              add_new_quad_intersection(quad_index, intersection->voronoi_edge_id, delaunay_edge_param, false, v_next);
            }
            else
            {
              add_new_quad_intersection(
                quad_index, intersection->voronoi_edge_id, delaunay_edge_param, false, v_edge_intersections.begin());
            }
          }
        }
        else if (containing_triangle_id != face_inside_new && intersected_before)
        {
          if (flip_edge_infinite)
          {
            crossing_data.removeIntersection(*v_intersection);
          }
          else
          {
            assignCrossingIntersectionDelaunayParam(this, **v_intersection,
              delaunayVoronoiEdgeIntersection(quad_index, (*v_intersection)->voronoi_edge_id, t).first, t,
              "flipUpdate:containedInNewFace");
          }
        }
      }
      ++d_it;
    }
  }

  {
    auto& quad_d_intersections = crossing_data.delaunay_edge_intersections[quad_index];
    // TODO: Sorting by recomputed dParam is no longer needed — list order is maintained by insert/erase.
    // Re-sorting can scramble coincident crossings whose params differ only by noise.
    // quad_d_intersections.sort(
    //   [&](const CrossingData::EdgeIntersectionRef& a, const CrossingData::EdgeIntersectionRef& b)
    //   { return a->delaunay_edge_param < b->delaunay_edge_param; });
    for (auto it = quad_d_intersections.begin(); it != quad_d_intersections.end(); ++it)
    {
      (*it)->delaunay_ref = it;
    }
  }

  // Recompute crossing events for Voronoi vertices in both triangles, including dual vertices
  // (face_id0/face_id1) that may not appear in getVoronoiVerticesInTri after reassignment.
  std::unordered_set<size_t> recompute_voronoi_vertices;
  for (size_t voronoi_vertex : vertices0)
  {
    recompute_voronoi_vertices.insert(voronoi_vertex);
  }
  for (size_t voronoi_vertex : vertices1)
  {
    recompute_voronoi_vertices.insert(voronoi_vertex);
  }
  recompute_voronoi_vertices.insert(face_id0);
  recompute_voronoi_vertices.insert(face_id1);
  for (size_t voronoi_vertex : recompute_voronoi_vertices)
  {
    crossing_event_manager_->computeEvents(t, voronoi_vertex, infinitesimal);
    if (voronoi_vertex < crossing_data.last_crossing.size())
    {
      crossing_data.last_crossing[voronoi_vertex]
        = infinitesimal.has_value() ? EventTime(t, infinitesimal->min_infinitesimal_t) : EventTime(t);
    }
  }
}

void KineticDelaunay::growGraphSlotArrays()
{
  const size_t face_slots = graph.faceSlotCount();
  if (face_inside.size() < face_slots)
  {
    face_inside.resize(face_slots, false);
  }
  if (face_last_updated.size() < face_slots)
  {
    face_last_updated.resize(face_slots, EventTime{});
  }
  crossing_data.growTo(face_slots);
  crossing_data.growEdgeSlotsTo(graph.halfEdgeSlotCount() / 2);

  const size_t quad_slots = graph.halfEdgeSlotCount() / 2;
  if (quadrilateral_last_updated.size() < quad_slots)
  {
    quadrilateral_last_updated.resize(quad_slots, EventTime{});
  }
}

size_t KineticDelaunay::findContainingTriForVoronoiVertex(size_t voronoi_vertex_id, double t) const
{
  const HalfEdgeDelaunayGraph::Triangle& tri = graph.face(voronoi_vertex_id);

  auto vertices = graph.adjacentTriangleVertices(tri.half_edges[0]);
  std::vector<glm::dvec2> points;
  bool outer_face = false;
  for (const auto& v : vertices)
  {
    if (v == -1)
    {
      outer_face = true;
      break;
    }
    points.push_back(getPointAt(static_cast<size_t>(v), t));
  }

  if (!outer_face)
  {
    double r = circumradius(points[0], points[1], points[2]);
    if (r < cutoff || isMinimalInputBranchTriangle(vertices, t))
    {
      // For existing Voronoi vertices we only recompute containment; face-inside state is managed elsewhere.
    }
  }

  if (outer_face)
  {
    return voronoi_vertex_id;
  }

  glm::dvec2 circumcenter = HalfEdgeDelaunayGraph::circumcenter(points[0], points[1], points[2]);

  auto edge_function = [](const glm::dvec2& a, const glm::dvec2& b, const glm::dvec2& c)
  { return -((b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x)); };

  auto edge_function_01 = edge_function(points[0], points[1], circumcenter);
  auto edge_function_12 = edge_function(points[1], points[2], circumcenter);
  auto edge_function_20 = edge_function(points[2], points[0], circumcenter);

  bool inside_triangle = false;
  glm::dvec2 start_point;

  if (edge_function_01 < 0)
  {
    start_point = (points[0] + points[1]) / 2.0;
  }
  else if (edge_function_12 < 0)
  {
    start_point = (points[1] + points[2]) / 2.0;
  }
  else if (edge_function_20 < 0)
  {
    start_point = (points[2] + points[0]) / 2.0;
  }
  else
  {
    inside_triangle = true;
  }

  if (inside_triangle)
  {
    return voronoi_vertex_id;
  }

  auto crossed_half_edges = computeCrossedHalfEdges(voronoi_vertex_id, circumcenter, start_point, t).first;
  if (crossed_half_edges.empty())
  {
    return voronoi_vertex_id;
  }

  const size_t last_crossed_edge_id = crossed_half_edges.back();
  return graph.halfEdge(last_crossed_edge_id ^ 1).face;
}

void KineticDelaunay::initializeFaceState(size_t face_index, double t)
{
  const HalfEdgeDelaunayGraph::Triangle& tri = graph.face(face_index);

  auto vertices = graph.adjacentTriangleVertices(tri.half_edges[0]);
  std::vector<glm::dvec2> points;
  bool outer_face = false;
  for (const auto& v : vertices)
  {
    if (v == -1)
    {
      outer_face = true;
      break;
    }
    points.push_back(getPointAt(static_cast<size_t>(v), t));
  }

  if (!outer_face)
  {
    double r = circumradius(points[0], points[1], points[2]);
    if (r < cutoff || isMinimalInputBranchTriangle(vertices, t))
    {
      setFaceInside(face_index, true, t);
    }
  }

  crossing_data.setVoronoiVertexTriId(face_index, findContainingTriForVoronoiVertex(face_index, t));
}

void KineticDelaunay::initializeNewFacesAfterGraphUpdate(double t, size_t first_new_face_slot)
{
  for (size_t face_index : graph.liveFaces())
  {
    if (face_index < first_new_face_slot)
    {
      continue;
    }
    initializeFaceState(face_index, t);
  }
}

void KineticDelaunay::onGraphRetriangulated(double t, size_t prev_face_slots, size_t prev_he_slots)
{
  growGraphSlotArrays();
  initializeNewFacesAfterGraphUpdate(t, prev_face_slots);
  // Runtime branch ids for the targeted split are already assigned at notePendingBranchSplit; do not remap all
  // live-graph components (that would finalize / reassign unrelated pending splits).
  crossing_data.computeEdgeIntersections(*this, t);
  validateVoronoiVertexIteratorInvariants("onGraphRetriangulated", t);
  if (callback_manager_)
  {
    callback_manager_->onGraphRetriangulated(t, prev_face_slots, prev_he_slots);
  }
}

void KineticDelaunay::onGraphCutApplied(double t, size_t prev_face_slots, size_t prev_he_slots,
  bool update_runtime_branch_map, const HalfEdgeDelaunayGraph::RuntimeBranchSplitResult* split_result)
{
  growGraphSlotArrays();
  initializeNewFacesAfterGraphUpdate(t, prev_face_slots);
  for (size_t face_id = 0; face_id < prev_face_slots; ++face_id)
  {
    if (face_id >= graph.faceSlotCount() || graph.isLiveFace(face_id))
    {
      continue;
    }
    const std::vector<size_t> contained_vertices = crossing_data.getVoronoiVerticesInTri(face_id);
    for (size_t voronoi_vertex_id : contained_vertices)
    {
      if (voronoi_vertex_id >= graph.faceSlotCount() || !graph.isLiveFace(voronoi_vertex_id))
      {
        continue;
      }
      const size_t target_tri_id = findContainingTriForVoronoiVertex(voronoi_vertex_id, t);
      if (crossing_data.isVoronoiVertexRegistered(voronoi_vertex_id)
        && crossing_data.getContainingTriId(voronoi_vertex_id) == target_tri_id)
      {
        continue;
      }
      crossing_data.moveVertex(voronoi_vertex_id, target_tri_id, t);
    }
  }

  for (size_t voronoi_vertex_id = 0; voronoi_vertex_id < graph.faceSlotCount(); ++voronoi_vertex_id)
  {
    if (!graph.isLiveFace(voronoi_vertex_id))
    {
      crossing_data.unsetVoronoiVertex(voronoi_vertex_id);
    }
  }

  crossing_data.removeIntersectionsOnDeadDelaunayEdges(graph);

  if (split_result != nullptr)
  {
    refreshEventsAfterGraphCut(t, *split_result);
  }

  if (update_runtime_branch_map)
  {
    updateRuntimeBranchMapFromLiveGraph(t, "onGraphCutApplied");
  }
  validateVoronoiVertexIteratorInvariants("onGraphCutApplied", t);
  validateCrossingIntersectionInvariants("onGraphCutApplied", t);
  if (callback_manager_)
  {
    callback_manager_->onGraphCutApplied(t, prev_face_slots, prev_he_slots);
  }
}

void KineticDelaunay::refreshEventsAfterGraphCut(
  double t, const HalfEdgeDelaunayGraph::RuntimeBranchSplitResult& split_result)
{
  // Flip schedules are keyed by undirected edge id (half_edge / 2).
  // 1) Infinite edges that bordered one live + one tombstoned face: recompute that edge's own quad.
  // 2) Finite capped outer edges (interior→hull): recompute those quads separately.
  std::unordered_set<size_t> flip_quads;
  std::vector<size_t> infinite_edge_ids = split_result.infinite_half_edges_bordering_tombstone;

  for (size_t he_id : infinite_edge_ids)
  {
    flip_quads.insert(he_id / 2);
  }

  for (size_t he_id : split_result.capped_outer_half_edges)
  {
    if (!graph.isLiveHalfEdge(he_id))
    {
      continue;
    }
    if (graph.halfEdge(he_id).origin >= 0 && graph.destination(he_id) >= 0)
    {
      flip_quads.insert(he_id / 2);
    }
  }

  std::vector<size_t> recomputed_quads;
  recomputed_quads.reserve(flip_quads.size());
  for (size_t quad_id : flip_quads)
  {
    if (!graph.isLiveHalfEdge(quad_id * 2) && !graph.isLiveHalfEdge(quad_id * 2 + 1))
    {
      continue;
    }
    flip_event_manager_->computeEvents(t, quad_id);
    if (quad_id < quadrilateral_last_updated.size())
    {
      quadrilateral_last_updated[quad_id] = t;
    }
    recomputed_quads.push_back(quad_id);
  }

  {
    std::sort(infinite_edge_ids.begin(), infinite_edge_ids.end());
    std::sort(recomputed_quads.begin(), recomputed_quads.end());

    std::ostringstream infinite_list;
    for (size_t i = 0; i < infinite_edge_ids.size(); ++i)
    {
      if (i > 0)
      {
        infinite_list << ',';
      }
      infinite_list << infinite_edge_ids[i];
    }
    std::ostringstream quad_list;
    for (size_t i = 0; i < recomputed_quads.size(); ++i)
    {
      if (i > 0)
      {
        quad_list << ',';
      }
      quad_list << recomputed_quads[i];
    }
    KINDS_INFO("Graph-cut flip refresh t=" << t
                                           << " infinite_half_edges_bordering_tombstone=[" << infinite_list.str()
                                           << "] flip_quad_ids=[" << quad_list.str() << "]");
  }
}

void KineticDelaunay::precomputeStep(double t)
{
  // TODO: Make sure it works where no change of sign occurs in the polynomial, i.e., roots that do not lead to a
  // change in the triangulation.
  for (size_t he_id : graph.liveDelaunayEdges())
  {
    flip_event_manager_->computeEvents(t, he_id / 2);
  }

  for (size_t face_id : graph.liveFaces())
  {
    const size_t he_id = graph.face(face_id).half_edges[0];
    radius_event_manager_->computeEvents(t, he_id);
  }

  for (size_t tri_id : graph.liveFaces())
  {
    crossing_event_manager_->computeEvents(t, tri_id);
  }
}

void KineticDelaunay::handleEvents()
{
  kinetic_algorithm_->processEvents(
    static_cast<double>(getEndSection()), collect_statistics_ ? &statistics_ : nullptr);
}

size_t KineticDelaunay::getBranchIndex(size_t strand_id, size_t t) const
{
  return branch_trajs.getBranchIndex(strand_id, t);
}

const std::vector<std::vector<size_t>>& KineticDelaunay::getBranches(size_t t) const
{
  return branch_trajs.getStrandsByBranchId()[t];
}

const std::vector<size_t>& KineticDelaunay::getBranchStrands(size_t t, size_t branch_id)
{
  return branch_trajs.getStrandsByBranchId()[t][branch_id];
}

KineticDelaunay::KineticDelaunay(const StrandTree& branch_trajs, double cutoff, bool add_dummy_splines)
  : branch_trajs(branch_trajs)
  , kinetic_algorithm_(std::make_unique<KineticAlgorithm>())
  , cutoff(cutoff)
  , add_dummy_boundary(add_dummy_splines)
  , flip_event_manager_(std::make_unique<FlipEventManager>(this))
  , radius_event_manager_(std::make_unique<RadiusEventManager>(this))
  , crossing_event_manager_(std::make_unique<CrossingEventManager>(this))
  , section_event_manager_(std::make_unique<SectionEventManager>(this))
  , subdivision_event_manager_(std::make_unique<SubdivisionEventManager>(this))
  , separation_event_manager_(std::make_unique<SeparationEventManager>(this))
  , crossing_data(crossing_event_manager_->getCrossingDataMutable())
{
  if (add_dummy_splines)
  {
    // first compute a bounding box:
    glm::dvec2 p_min { std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity() };
    glm::dvec2 p_max { -std::numeric_limits<double>::infinity(), -std::numeric_limits<double>::infinity() };

    for (const auto& points : branch_trajs.getPoints())
    {
      for (const auto& p : points)
      {
        for (int dim = 0; dim < 2; dim++)
        {
          if (p[dim] < p_min[dim])
          {
            p_min[dim] = p[dim];
          }
          if (p[dim] > p_max[dim])
          {
            p_max[dim] = p[dim];
          }
        }
      }
    }

    // We will need dummy points such that no voronoi vertices can slip outside
    double range = std::max(p_max[0] - p_min[0], p_max[1] - p_min[1]);
    double dist_from_bb = std::max(range, 2 * cutoff);

    dummy_boundary = {
      { p_min[0] - 0.75 * dist_from_bb, p_max[1] + 0.75 * dist_from_bb }, // corner_tl
      { p_min[0], p_max[1] + dist_from_bb }, // top_left
      { p_max[0], p_max[1] + dist_from_bb }, // top_right
      { p_max[0] + 0.75 * dist_from_bb, p_max[1] + 0.75 * dist_from_bb }, // corner_tr
      { p_max[0] + dist_from_bb, p_max[1] }, // right_top
      { p_max[0] + dist_from_bb, p_min[1] }, // right_bottom
      { p_max[0] + 0.75 * dist_from_bb, p_min[1] - 0.75 * dist_from_bb }, // corner_br
      { p_max[0], p_min[1] - dist_from_bb }, // bottom_right
      { p_min[0], p_min[1] - dist_from_bb }, // bottom_left
      { p_min[0] - 0.75 * dist_from_bb, p_min[1] - 0.75 * dist_from_bb }, // corner_bl
      { p_min[0] - dist_from_bb, p_min[1] }, // left_bottom
      { p_min[0] - dist_from_bb, p_max[1] } // left_top
    };

    size_t length = branch_trajs.getHeight() + 1;

    for (const auto& p : dummy_boundary)
    {
      std::vector<glm::dvec2> new_spline;
      for (size_t i = 0; i < length; i++)
      {
        new_spline.push_back(p);
      }
      this->branch_trajs.addTrajectory(new_spline);
    }
  }
}

KineticDelaunay::~KineticDelaunay() = default;

bool KineticDelaunay::isDummyBoundary(size_t v) const
{
  if (add_dummy_boundary)
  {
    return v >= branch_trajs.getPoints().size() - 12;
  }
  return false;
}

void KineticDelaunay::collectReferenceBranchStrandPool(size_t strand_id, std::vector<size_t>& pool) const
{
  pool.clear();

  if (!isGraphRetriangulatedForComponents())
  {
    if (const std::vector<size_t>* frozen_strands = pending_branch_splits_.frozenStrandsForStrand(strand_id);
      frozen_strands != nullptr)
    {
      pool = *frozen_strands;
      return;
    }
  }

  if (strand_id < component_data.component_map.size())
  {
    const size_t component_id = component_data.component_map[strand_id];
    if (component_id < component_data.components.size())
    {
      pool = component_data.components[component_id];
    }
  }
}

size_t KineticDelaunay::getSharedReferenceBranchForStrands(
  const std::vector<size_t>& strand_ids, double branch_lookup_time) const
{
  const size_t branch_lookup_height = branchLookupHeightForTime(branch_trajs, branch_lookup_time);

  std::vector<size_t> pool;
  pool.reserve(strand_ids.size() * 4);
  for (size_t strand_id : strand_ids)
  {
    std::vector<size_t> per_strand_pool;
    collectReferenceBranchStrandPool(strand_id, per_strand_pool);
    for (size_t pooled_strand_id : per_strand_pool)
    {
      if (std::find(pool.begin(), pool.end(), pooled_strand_id) == pool.end())
      {
        pool.push_back(pooled_strand_id);
      }
    }
  }

  if (pool.empty())
  {
    if (strand_ids.empty())
    {
      return 0;
    }
    return branch_trajs.getBranchIndex(strand_ids.front(), branch_lookup_height);
  }

  return minInputBranchForStrands(pool, branch_lookup_height);
}

size_t KineticDelaunay::getReferenceBranch(size_t strand_id, double t) const
{
  const size_t section = static_cast<size_t>(std::floor(std::max(0.0, t)));
  if (const PostSplitFrameTransition* transition = postSplitFrameTransitionForStrand(strand_id);
    transition != nullptr && section < transition->hold_end_height)
  {
    // During the hold window, motion stays in the shared pre-split frame.
    return transition->common_reference_branch;
  }

  std::vector<size_t> pool;
  collectReferenceBranchStrandPool(strand_id, pool);
  const size_t branch_lookup_height = branchLookupHeightForTime(branch_trajs, t);
  if (pool.empty())
  {
    return branch_trajs.getBranchIndex(strand_id, branch_lookup_height);
  }

  return minInputBranchForStrands(pool, branch_lookup_height);
}

size_t KineticDelaunay::inputBranchSectionIndexAtIntervalUpperBound(double event_interval_upper_bound) const
{
  const size_t tree_height = branch_trajs.getHeight();
  size_t section_index = static_cast<size_t>(event_interval_upper_bound);
  if (tree_height > 0 && section_index >= tree_height)
  {
    section_index = tree_height - 1;
  }
  return section_index;
}

std::vector<size_t> KineticDelaunay::collectDistinctInputBranchesForEventTrigger(
  const std::vector<size_t>& event_strand_ids, double event_interval_upper_bound) const
{
  const size_t branch_section_index = inputBranchSectionIndexAtIntervalUpperBound(event_interval_upper_bound);
  std::vector<size_t> distinct_branches;
  for (size_t strand_id : event_strand_ids)
  {
    if (isDummyBoundary(strand_id))
    {
      continue;
    }

    const size_t input_branch = branch_trajs.getBranchIndex(strand_id, branch_section_index);
    if (std::find(distinct_branches.begin(), distinct_branches.end(), input_branch) == distinct_branches.end())
    {
      distinct_branches.push_back(input_branch);
    }
  }

  std::sort(distinct_branches.begin(), distinct_branches.end());
  return distinct_branches;
}

bool KineticDelaunay::eventTriggerUsesSharedTransformedFrame(
  const std::vector<size_t>& event_strand_ids, double event_interval_upper_bound) const
{
  return collectDistinctInputBranchesForEventTrigger(event_strand_ids, event_interval_upper_bound).size() > 1;
}

size_t KineticDelaunay::sharedReferenceBranchForEventTrigger(
  const std::vector<size_t>& event_strand_ids, double event_interval_upper_bound) const
{
  const size_t branch_section_index = inputBranchSectionIndexAtIntervalUpperBound(event_interval_upper_bound);
  return minInputBranchForStrands(event_strand_ids, branch_section_index);
}

Trajectory<2> KineticDelaunay::getSitePiecePolynomialForEventStrands(size_t strand_id, size_t section,
  double schedule_time, const std::vector<size_t>& event_strand_ids) const
{
  if (computing_infinitesimal_events_)
  {
    return buildInfinitesimalSiteTrajectory(strand_id, infinitesimal_schedule_t_, event_strand_ids);
  }

  const double event_interval_upper_bound = eventIntervalUpperBound(schedule_time);
  Trajectory<2> base;

  const PostSplitFrameTransition* transition = postSplitFrameTransitionForStrand(strand_id);
  const bool transition_active
    = transition != nullptr && section < transition->hold_end_height + 1;

  if (transition_active && section < transition->hold_end_height)
  {
    // Hold: both piece endpoints in the shared pre-split frame.
    base = branch_trajs.getPiecePolynomial(strand_id, section, transition->common_reference_branch);
  }
  else if (transition_active && section == transition->blend_section)
  {
    // Blend section [S+1, S+2]: common frame at S+1 → native branch frame at S+2.
    const size_t native_ref = branch_trajs.getBranchIndex(strand_id, section + 1);
    base = branch_trajs.getPiecePolynomialBlendingReference(
      strand_id, section, transition->common_reference_branch, native_ref);
  }
  else if (eventTriggerUsesSharedTransformedFrame(event_strand_ids, event_interval_upper_bound))
  {
    const size_t reference_branch
      = sharedReferenceBranchForEventTrigger(event_strand_ids, event_interval_upper_bound);
    base = branch_trajs.getPiecePolynomial(strand_id, section, reference_branch);
  }
  else
  {
    base = branch_trajs.getLocalPiecePolynomial(strand_id, section);
  }

  // Kinetic separation ramp removed: infinitesimal virtual offset is handled via buildInfinitesimalSiteTrajectory.
  (void)schedule_time;
  return base;
}

glm::dvec2 KineticDelaunay::getPointAtWithReferenceBranch(size_t v, double t, size_t reference_branch,
  bool apply_reference_transform, bool include_virtual_offset) const
{
  const std::vector<glm::dvec2>& support = branch_trajs.getSupportPoints(v);
  if (support.empty())
  {
    throw std::out_of_range("getPointAt: strand " + std::to_string(v) + " has no support points");
  }

  const size_t last_index = support.size() - 1;
  double query_t = t;
  if (static_cast<size_t>(std::floor(query_t)) >= last_index)
  {
    query_t = static_cast<double>(last_index);
  }

  const size_t section = static_cast<size_t>(std::floor(query_t));
  const double frac = query_t - static_cast<double>(section);

  glm::dvec2 position;
  if (!apply_reference_transform)
  {
    position = branch_trajs.evaluate(v, query_t);
  }
  else
  {
    const PostSplitFrameTransition* transition = postSplitFrameTransitionForStrand(v);
    const bool transition_active
      = transition != nullptr && section < transition->hold_end_height + 1;

    if (transition_active && section < transition->hold_end_height)
    {
      if (frac < std::numeric_limits<double>::epsilon())
      {
        position = branch_trajs.getPointTransformed(v, section, transition->common_reference_branch);
      }
      else
      {
        const Trajectory<2> piece
          = branch_trajs.getPiecePolynomial(v, section, transition->common_reference_branch);
        position = glm::dvec2(piece[0](frac), piece[1](frac));
      }
    }
    else if (transition_active && section == transition->blend_section)
    {
      const size_t native_ref = branch_trajs.getBranchIndex(v, section + 1);
      if (frac < std::numeric_limits<double>::epsilon())
      {
        position = branch_trajs.getPointTransformedAtSection(
          v, section, transition->common_reference_branch, section);
      }
      else
      {
        const Trajectory<2> piece = branch_trajs.getPiecePolynomialBlendingReference(
          v, section, transition->common_reference_branch, native_ref);
        position = glm::dvec2(piece[0](frac), piece[1](frac));
      }
    }
    else if (frac < std::numeric_limits<double>::epsilon())
    {
      position = branch_trajs.getPointTransformed(v, section, reference_branch);
    }
    else
    {
      const Trajectory<2> piece = branch_trajs.getPiecePolynomial(v, section, reference_branch);
      position = glm::dvec2(piece[0](frac), piece[1](frac));
    }
  }

  if (include_virtual_offset)
  {
    position += separationOffsetAt(v, t);
  }
  return position;
}

glm::dvec2 KineticDelaunay::getPointAt(
  size_t v, double t, bool apply_reference_transform, bool include_virtual_offset) const
{
  return getPointAtWithReferenceBranch(
    v, t, getReferenceBranch(v, t), apply_reference_transform, include_virtual_offset);
}

glm::dvec2 KineticDelaunay::getPointAt(
  double t, size_t v, bool apply_reference_transform, bool include_virtual_offset) const
{
  return getPointAt(v, t, apply_reference_transform, include_virtual_offset);
}

glm::dvec2 KineticDelaunay::getPointInDelaunaySpace(size_t v, double t) const
{
  return getPointAt(v, t, /*apply_reference_transform=*/true, /*include_virtual_offset=*/true);
}

glm::dvec2 KineticDelaunay::getPointInDelaunaySpace(size_t v, double t, size_t reference_branch) const
{
  return getPointAtWithReferenceBranch(
    v, t, reference_branch, /*apply_reference_transform=*/true, /*include_virtual_offset=*/true);
}

Trajectory<2> KineticDelaunay::getSitePiecePolynomialWithReferenceBranch(
  size_t strand_id, size_t section, size_t reference_branch, double schedule_time) const
{
  const Trajectory<2> base = branch_trajs.getPiecePolynomial(strand_id, section, reference_branch);
  return addSeparationOffsetToPiecePolynomial(base, strand_id, section, schedule_time);
}

Trajectory<2> KineticDelaunay::getSitePiecePolynomial(size_t strand_id, size_t section, double schedule_time) const
{
  const PostSplitFrameTransition* transition = postSplitFrameTransitionForStrand(strand_id);
  const bool transition_active
    = transition != nullptr && section < transition->hold_end_height + 1;

  Trajectory<2> base;
  if (transition_active && section < transition->hold_end_height)
  {
    base = branch_trajs.getPiecePolynomial(strand_id, section, transition->common_reference_branch);
  }
  else if (transition_active && section == transition->blend_section)
  {
    const size_t native_ref = branch_trajs.getBranchIndex(strand_id, section + 1);
    base = branch_trajs.getPiecePolynomialBlendingReference(
      strand_id, section, transition->common_reference_branch, native_ref);
  }
  else
  {
    const double branch_lookup_time = referenceBranchLookupTimeForSection(section, schedule_time);
    const size_t reference_branch = getReferenceBranch(strand_id, branch_lookup_time);
    base = branch_trajs.getPiecePolynomial(strand_id, section, reference_branch);
  }

  return addSeparationOffsetToPiecePolynomial(base, strand_id, section, schedule_time);
}

std::vector<glm::dvec2> KineticDelaunay::getPointsAt(
  double t, bool apply_reference_transform, bool include_virtual_offset) const
{
  size_t vertex_count = graph.getVertexCount();
  std::vector<glm::dvec2> points;
  points.reserve(vertex_count);
  for (size_t v = 0; v < vertex_count; v++)
  {
    points.push_back(getPointAt(v, t, apply_reference_transform, include_virtual_offset));
  }
  return points;
}

glm::dvec3 KineticDelaunay::getPointInObjectSpace(size_t v, double t) const
{
  return branch_trajs.getPointInObjectSpace(v, t);
}

const StrandTree& KineticDelaunay::getStrandTree() const { return branch_trajs; }

void KineticDelaunay::setVisualDebugOutputRoot(const std::filesystem::path& root)
{
  visual_debug_output_root_ = root;
  runtime_branch_log_header_written_ = false;

  // Truncate any previous run's log so appends in logRuntimeBranchEvent start fresh.
  std::error_code ec;
  std::filesystem::create_directories(root, ec);
  std::ofstream clear_log(root / "branches.txt", std::ios::trunc);
}

void KineticDelaunay::setVisualDebugEnabled(bool enabled)
{
  visual_debug_enabled_ = enabled;
}

bool KineticDelaunay::isVisualDebugEnabled() const
{
  return visual_debug_enabled_;
}

void KineticDelaunay::setErrorFilesEnabled(bool enabled)
{
  error_files_enabled_ = enabled;
}

bool KineticDelaunay::isErrorFilesEnabled() const
{
  return error_files_enabled_;
}

bool KineticDelaunay::shouldDumpErrorFiles() const
{
  return error_files_enabled_ || visual_debug_enabled_;
}

void KineticDelaunay::setVisualDebugSeparatePendingSplits(bool enabled)
{
  visual_debug_separate_pending_splits_ = enabled;
}

bool KineticDelaunay::visualDebugSeparatePendingSplits() const
{
  return visual_debug_separate_pending_splits_;
}

std::optional<std::filesystem::path> KineticDelaunay::getRuntimeBranchLogPath() const
{
  if (!visual_debug_enabled_ || !visual_debug_output_root_.has_value())
  {
    return std::nullopt;
  }
  return *visual_debug_output_root_ / "branches.txt";
}

void KineticDelaunay::logRuntimeBranchEvent(double t, const std::string& event_line)
{
  const auto log_path = getRuntimeBranchLogPath();
  if (!log_path.has_value())
  {
    return;
  }

  std::ofstream out(*log_path, std::ios::app);
  if (!out)
  {
    return;
  }

  if (!runtime_branch_log_header_written_)
  {
    out << "# Runtime branch event log. Ids are monotonic and never reused.\n"
        << "# event=create  new runtime branch id allocated\n"
        << "# event=assign  existing runtime branch id assigned to a live-graph component\n"
        << "# event=retire  runtime branch id marked dead (no longer active)\n"
        << "#\n";
    runtime_branch_log_header_written_ = true;
  }

  out << "t=" << std::fixed << std::setprecision(6) << t << ' ' << event_line << '\n';
  out.flush();
}

const std::optional<std::filesystem::path>& KineticDelaunay::getVisualDebugOutputRoot() const
{
  return visual_debug_output_root_;
}

void KineticDelaunay::setFlipPolynomialDumpTargetTime(std::optional<double> target_time)
{
  flip_polynomial_dump_target_time_ = target_time;
}

const std::optional<double>& KineticDelaunay::getFlipPolynomialDumpTargetTime() const
{
  return flip_polynomial_dump_target_time_;
}

void KineticDelaunay::setFlipPolynomialDumpTargetHalfEdge(std::optional<size_t> half_edge_id)
{
  flip_polynomial_dump_target_half_edge_ = half_edge_id;
}

const std::optional<size_t>& KineticDelaunay::getFlipPolynomialDumpTargetHalfEdge() const
{
  return flip_polynomial_dump_target_half_edge_;
}

void KineticDelaunay::computeComponentData(double t)
{
  auto& graph = getGraph();
  component_data.components = extractConnectedComponents();
  // KINDS_DEBUG("Extracted " << component_data.components.size() << " components.");
  component_data.component_map = buildComponentMap(component_data.components, graph.getVertexCount());
  component_data.component_boundaries.resize(component_data.components.size());

  std::vector<bool> he_visited(graph.halfEdgeSlotCount(), false);

  for (size_t component_index = 0; component_index < component_data.components.size(); component_index++)
  {
    component_data.component_boundaries[component_index]
      = extractComponentBoundaries(component_data.components[component_index], t, he_visited);
  }

  component_data.component_centroids.resize(component_data.components.size());
  for (size_t component_index = 0; component_index < component_data.components.size(); component_index++)
  {
    if (!component_data.component_boundaries[component_index].empty())
    {
      component_data.component_centroids[component_index]
        = polygonCentroid(component_data.component_boundaries[component_index][0]);
    }
    else
    {
      // compute centroid from points in the component
      glm::dvec2 centroid { 0.0, 0.0 };
      for (auto& v : component_data.components[component_index])
      {
        glm::dvec2 p = getPointAt(t, v);
        centroid += p;
      }
      component_data.component_centroids[component_index]
        = centroid / double(component_data.components[component_index].size());
    }
  }

  component_data.component_last_updated.resize(component_data.components.size(), t);
}

void KineticDelaunay::clearComponentSupportData(size_t component_id)
{
  if (component_id >= component_data.components.size())
  {
    return;
  }

  component_data.components[component_id].clear();
  if (component_id < component_data.component_boundaries.size())
  {
    component_data.component_boundaries[component_id].clear();
  }
  if (component_id < component_data.component_centroids.size())
  {
    component_data.component_centroids[component_id] = glm::dvec2(0.0, 0.0);
  }
  if (component_id < component_data.component_last_updated.size())
  {
    component_data.component_last_updated[component_id] = std::numeric_limits<double>::quiet_NaN();
  }
}

void KineticDelaunay::consolidateComponentsAtTriangle(const std::array<int, 3>& triangle_vertices, double t)
{
  if (component_data.components.empty() || component_data.component_map.empty())
  {
    return;
  }

  std::vector<size_t> triangle_component_ids;
  triangle_component_ids.reserve(3);
  for (int vertex : triangle_vertices)
  {
    if (vertex < 0)
    {
      continue;
    }
    const size_t strand_id = static_cast<size_t>(vertex);
    if (isDummyBoundary(strand_id) || strand_id >= component_data.component_map.size())
    {
      continue;
    }
    triangle_component_ids.push_back(component_data.component_map[strand_id]);
  }
  if (triangle_component_ids.size() < 2)
  {
    return;
  }

  std::sort(triangle_component_ids.begin(), triangle_component_ids.end());
  triangle_component_ids.erase(
    std::unique(triangle_component_ids.begin(), triangle_component_ids.end()), triangle_component_ids.end());
  if (triangle_component_ids.size() <= 1)
  {
    return;
  }

  // Smallest component index survives; every other id on this triangle is absorbed into it.
  const size_t survivor = triangle_component_ids.front();
  if (survivor >= component_data.components.size())
  {
    return;
  }

  // Merging pending-split pieces makes finalize ill-formed; keep separation running but put those splits on hiatus.
  putPendingSplitsOnHiatusTouchingComponents(triangle_component_ids, t);

  for (size_t i = 1; i < triangle_component_ids.size(); ++i)
  {
    const size_t absorbed_id = triangle_component_ids[i];
    if (absorbed_id >= component_data.components.size())
    {
      continue;
    }

    for (size_t strand_id : component_data.components[absorbed_id])
    {
      if (strand_id < component_data.component_map.size())
      {
        component_data.component_map[strand_id] = survivor;
      }
      component_data.components[survivor].push_back(strand_id);
    }
    clearComponentSupportData(absorbed_id);
  }

  component_data.component_boundaries.resize(component_data.components.size());
  component_data.component_centroids.resize(component_data.components.size());
  component_data.component_last_updated.resize(component_data.components.size());

  std::vector<bool> he_visited(graph.halfEdgeSlotCount(), false);
  component_data.component_boundaries[survivor]
    = extractComponentBoundaries(component_data.components[survivor], t, he_visited, false);
  if (!component_data.component_boundaries[survivor].empty()
    && !component_data.component_boundaries[survivor][0].empty())
  {
    component_data.component_centroids[survivor]
      = polygonCentroid(component_data.component_boundaries[survivor][0]);
  }
  else
  {
    glm::dvec2 centroid { 0.0, 0.0 };
    for (size_t strand_id : component_data.components[survivor])
    {
      centroid += getPointAt(t, strand_id);
    }
    const double n = static_cast<double>(component_data.components[survivor].size());
    component_data.component_centroids[survivor] = n > 0.0 ? centroid / n : glm::dvec2(0.0, 0.0);
  }
  component_data.component_last_updated[survivor] = t;

  KINDS_DEBUG("Component consolidate at t=" << t << " survivor=" << survivor << " absorbed_count="
                                           << (triangle_component_ids.size() - 1));
}

void KineticDelaunay::putPendingSplitsOnHiatusTouchingComponents(const std::vector<size_t>& component_ids, double t)
{
  if (component_ids.empty() || pending_branch_splits_.by_parent_.empty())
  {
    return;
  }

  std::unordered_set<size_t> touched(component_ids.begin(), component_ids.end());
  for (auto& entry : pending_branch_splits_.by_parent_)
  {
    PendingBranchSplit& split = entry.second;
    if (split.on_hiatus || split.split_component_ids.size() < 2)
    {
      continue;
    }

    bool overlaps = false;
    for (size_t piece_id : split.split_component_ids)
    {
      if (touched.count(piece_id) > 0)
      {
        overlaps = true;
        break;
      }
    }
    if (!overlaps)
    {
      continue;
    }

    split.on_hiatus = true;
    KINDS_DEBUG("Pending split parent_component_id=" << entry.first << " parent_runtime="
                                                     << split.parent_runtime_branch << " put on hiatus at t=" << t
                                                     << " (component consolidate)");
  }
}

bool KineticDelaunay::isPendingRuntimeBranchOnHiatus(size_t parent_runtime_branch_id) const
{
  if (parent_runtime_branch_id == RuntimeBranchData::no_branch)
  {
    return false;
  }
  for (const auto& entry : pending_branch_splits_.by_parent_)
  {
    const PendingBranchSplit& split = entry.second;
    if (split.parent_runtime_branch == parent_runtime_branch_id && split.on_hiatus)
    {
      return true;
    }
  }
  return false;
}

bool KineticDelaunay::isGraphRetriangulatedForComponents() const
{
  return component_data.components.size() <= prev_component_count;
}

void KineticDelaunay::notePendingBranchSplit(size_t parent_component_id, double split_time,
  const std::vector<size_t>& pre_split_parent_strands, const std::vector<std::vector<size_t>>& new_components,
  const std::vector<size_t>& split_component_ids)
{
  if (parent_component_id >= component_data.components.size() || new_components.empty()
    || pre_split_parent_strands.empty() || split_component_ids.size() < 2)
  {
    return;
  }

  const size_t branch_section = pendingSplitBranchSection(split_time);
  const std::vector<size_t> input_branch_ids
    = collectInputBranchIdsForStrandGroups(new_components, branch_section);

  // Same input-branch family on hiatus: reactivate that pending split. Separation is already driving it.
  if (PendingBranchSplit* hiatus_split = findHiatusPendingSplitWithInputBranches(input_branch_ids))
  {
    const size_t hiatus_parent = hiatus_split->parent_component_id;
    hiatus_split->on_hiatus = false;
    hiatus_split->split_component_ids = split_component_ids;
    hiatus_split->input_branch_ids = input_branch_ids;
    pending_branch_splits_.registerStrandsForSplit(hiatus_parent, new_components);
    assignPendingSplitChildRuntimeBranches(
      *hiatus_split, hiatus_parent, split_component_ids, branch_section, split_time, /*allow_allocate=*/false);
    refreshPendingSplitSeparationCentroids(hiatus_parent, split_time);
    KINDS_DEBUG("Pending split parent_component_id=" << hiatus_parent << " parent_runtime="
                                                     << hiatus_split->parent_runtime_branch
                                                     << " reactivated from hiatus at t=" << split_time
                                                     << " (noted_parent=" << parent_component_id << ")");
    // Clear active so activate re-enters virtual separation / cut at this time.
    hiatus_split->infinitesimal_active = false;
    activateInfinitesimalSeparationOrApplyCut(hiatus_parent, split_time);
    return;
  }

  PendingBranchSplit& split = pending_branch_splits_.getOrCreate(parent_component_id);
  const bool first_note_for_parent = split.frozen_parent_strands.empty();
  if (first_note_for_parent)
  {
    split.reference_vertex = pre_split_parent_strands.front();
    split.frozen_parent_strands = pre_split_parent_strands;
  }

  pending_branch_splits_.registerStrandsForSplit(parent_component_id, new_components);
  split.split_component_ids = split_component_ids;
  split.input_branch_ids = input_branch_ids;
  split.on_hiatus = false;

  // Mirror the pending split into the runtime-branch structure and write branches.txt before post-split
  // frame registration / separation scheduling (those paths can throw).
  const size_t parent_branch = runtime_branch_data_.branchForStrand(split.reference_vertex);
  split.parent_runtime_branch = parent_branch;
  assignPendingSplitChildRuntimeBranches(
    split, parent_component_id, split_component_ids, branch_section, split_time, /*allow_allocate=*/true);

  if (first_note_for_parent)
  {
    registerPostSplitFrameTransition(parent_component_id, split_time, pre_split_parent_strands);
  }

  refreshPendingSplitSeparationCentroids(parent_component_id, split_time);
}

size_t KineticDelaunay::pendingSplitBranchSection(double split_time) const
{
  const size_t height = branch_trajs.getHeight();
  if (height == 0)
  {
    return 0;
  }
  size_t section_index = static_cast<size_t>(std::ceil(split_time));
  if (section_index >= height)
  {
    section_index = height - 1;
  }
  return section_index;
}

std::vector<size_t> KineticDelaunay::collectInputBranchIdsForStrandGroups(
  const std::vector<std::vector<size_t>>& strand_groups, size_t branch_section) const
{
  std::vector<size_t> input_branch_ids;
  for (const std::vector<size_t>& strands : strand_groups)
  {
    for (size_t strand_id : strands)
    {
      if (isDummyBoundary(strand_id))
      {
        continue;
      }
      const size_t input_branch = branch_trajs.getBranchIndex(strand_id, branch_section);
      input_branch_ids.push_back(input_branch);
    }
  }
  std::sort(input_branch_ids.begin(), input_branch_ids.end());
  input_branch_ids.erase(std::unique(input_branch_ids.begin(), input_branch_ids.end()), input_branch_ids.end());
  return input_branch_ids;
}

PendingBranchSplit* KineticDelaunay::findHiatusPendingSplitWithInputBranches(
  const std::vector<size_t>& input_branch_ids)
{
  if (input_branch_ids.empty())
  {
    return nullptr;
  }
  for (auto& entry : pending_branch_splits_.by_parent_)
  {
    PendingBranchSplit& split = entry.second;
    if (!split.on_hiatus || split.input_branch_ids != input_branch_ids)
    {
      continue;
    }
    return &split;
  }
  return nullptr;
}

void KineticDelaunay::assignPendingSplitChildRuntimeBranches(PendingBranchSplit& split, size_t parent_component_id,
  const std::vector<size_t>& split_component_ids, size_t branch_section, double split_time, bool allow_allocate)
{
  const size_t parent_branch = split.parent_runtime_branch;
  if (parent_branch == RuntimeBranchData::no_branch || split_component_ids.size() < 2)
  {
    return;
  }

  const std::vector<size_t>* existing_children = runtime_branch_data_.pendingChildBranches(parent_branch);
  if (!allow_allocate && existing_children == nullptr)
  {
    return;
  }
  if (allow_allocate && existing_children != nullptr)
  {
    // Already mirrored for this parent runtime; still refresh strand → child assignments below.
  }
  else if (allow_allocate)
  {
    std::ostringstream line;
    line << "event=pending_split parent_component=" << parent_component_id << " parent_runtime=" << parent_branch
         << " piece_count=" << split_component_ids.size() << " section=" << branch_section;
    logRuntimeBranchEvent(split_time, line.str());
  }

  // Prefer reconstructing input_branch → child runtime from strands already on pending children.
  std::unordered_map<size_t, size_t> input_branch_to_child_runtime;
  if (existing_children != nullptr)
  {
    for (size_t child_branch : *existing_children)
    {
      for (size_t strand_id = 0; strand_id < runtime_branch_data_.branch_map.size(); ++strand_id)
      {
        if (runtime_branch_data_.branch_map[strand_id] != child_branch || isDummyBoundary(strand_id))
        {
          continue;
        }
        input_branch_to_child_runtime.emplace(
          branch_trajs.getBranchIndex(strand_id, branch_section), child_branch);
      }
    }
  }

  std::vector<size_t> newly_allocated_children;
  bool assigned_any_strand = false;
  for (size_t i = 1; i < split_component_ids.size(); ++i)
  {
    const size_t split_off_component = split_component_ids[i];
    if (split_off_component >= component_data.components.size())
    {
      continue;
    }

    std::optional<size_t> child_input_branch;
    for (size_t strand_id : component_data.components[split_off_component])
    {
      if (isDummyBoundary(strand_id))
      {
        continue;
      }
      child_input_branch = branch_trajs.getBranchIndex(strand_id, branch_section);
      break;
    }
    if (!child_input_branch.has_value())
    {
      continue;
    }

    size_t child_branch = RuntimeBranchData::no_branch;
    const auto existing = input_branch_to_child_runtime.find(child_input_branch.value());
    bool allocated_new_child = false;
    if (existing != input_branch_to_child_runtime.end())
    {
      child_branch = existing->second;
    }
    else if (allow_allocate && existing_children == nullptr)
    {
      child_branch = runtime_branch_data_.allocate(parent_branch);
      input_branch_to_child_runtime.emplace(child_input_branch.value(), child_branch);
      newly_allocated_children.push_back(child_branch);
      allocated_new_child = true;
    }
    else
    {
      continue;
    }

    // Reassign the split-off strands to the child runtime branch now, so runtime branch ids distinguish the future
    // sub-branches while the graph is still connected. Motion frames use component data (not the runtime branch
    // map), so this does not perturb geometry. When this parent's cut succeeds, completePendingRuntimeBranchSplit
    // keeps these child ids as final without remapping unrelated pending splits.
    size_t assigned_strand_count = 0;
    for (size_t strand_id : component_data.components[split_off_component])
    {
      if (isDummyBoundary(strand_id) || strand_id >= runtime_branch_data_.branch_map.size())
      {
        continue;
      }
      runtime_branch_data_.branch_map[strand_id] = child_branch;
      ++assigned_strand_count;
    }
    assigned_any_strand = assigned_any_strand || assigned_strand_count > 0;

    std::ostringstream line;
    line << "event=pending_split_child id=" << child_branch << " parent=" << parent_branch
         << " component=" << split_off_component << " strand_count=" << assigned_strand_count
         << " input_branch=" << child_input_branch.value() << " section=" << branch_section;
    if (!allocated_new_child)
    {
      line << " reused=1";
    }
    if (!allow_allocate)
    {
      line << " reactivate=1";
    }
    logRuntimeBranchEvent(split_time, line.str());
  }

  bool did_note_split = false;
  if (!newly_allocated_children.empty())
  {
    runtime_branch_data_.noteSplit(parent_branch, std::move(newly_allocated_children));
    did_note_split = true;
  }
  if (assigned_any_strand || existing_children != nullptr || did_note_split)
  {
    runtime_branch_data_.rebuildBranchStrandLists();
  }
}

void KineticDelaunay::refreshPendingSplitSeparationCentroids(size_t parent_component_id, double t)
{
  auto it = pending_branch_splits_.by_parent_.find(parent_component_id);
  if (it == pending_branch_splits_.by_parent_.end() || it->second.split_component_ids.size() < 2)
  {
    return;
  }

  PendingBranchSplit& split = it->second;
  // Store direction in native-transformed Delaunay ambient space (matches getPointAt / SVG offsets).
  const glm::dvec2 dir
    = computeSeparationDirection(split, t, /*apply_reference_transform=*/true, /*shared_reference_branch=*/std::nullopt);
  if (glm::dot(dir, dir) <= 0.0)
  {
    return;
  }

  // Keep centroids for diagnostics in the same frame as the stored direction.
  const auto piece_centroid = [&](size_t begin, size_t end) -> std::optional<glm::dvec2>
  {
    glm::dvec2 sum(0.0);
    double weight = 0.0;
    for (size_t i = begin; i < end; ++i)
    {
      const size_t component_id = split.split_component_ids[i];
      if (component_id >= component_data.components.size())
      {
        continue;
      }
      for (size_t strand_id : component_data.components[component_id])
      {
        if (isDummyBoundary(strand_id))
        {
          continue;
        }
        sum += getPointAt(strand_id, t, /*apply_reference_transform=*/true, /*include_virtual_offset=*/false);
        weight += 1.0;
      }
    }
    if (weight <= 0.0)
    {
      return std::nullopt;
    }
    return sum / weight;
  };

  const std::optional<glm::dvec2> old_centroid = piece_centroid(0, 1);
  const std::optional<glm::dvec2> new_centroid = piece_centroid(1, split.split_component_ids.size());
  if (old_centroid.has_value() && new_centroid.has_value())
  {
    split.old_branch_centroid = old_centroid.value();
    split.new_branch_centroid = new_centroid.value();
    split.separation_direction = new_centroid.value() - old_centroid.value();
  }
}

glm::dvec2 KineticDelaunay::computeSeparationDirection(const PendingBranchSplit& split, double t,
  bool apply_reference_transform, std::optional<size_t> shared_reference_branch) const
{
  const auto eval = [&](size_t strand_id) -> glm::dvec2
  {
    if (shared_reference_branch.has_value())
    {
      return getPointAtWithReferenceBranch(
        strand_id, t, *shared_reference_branch, apply_reference_transform, /*include_virtual_offset=*/false);
    }
    return getPointAt(strand_id, t, apply_reference_transform, /*include_virtual_offset=*/false);
  };

  const auto piece_centroid = [&](size_t begin, size_t end) -> std::optional<glm::dvec2>
  {
    glm::dvec2 sum(0.0);
    double weight = 0.0;
    for (size_t i = begin; i < end; ++i)
    {
      const size_t component_id = split.split_component_ids[i];
      if (component_id >= component_data.components.size())
      {
        continue;
      }
      for (size_t strand_id : component_data.components[component_id])
      {
        if (isDummyBoundary(strand_id))
        {
          continue;
        }
        sum += eval(strand_id);
        weight += 1.0;
      }
    }
    if (weight <= 0.0)
    {
      return std::nullopt;
    }
    return sum / weight;
  };

  if (split.split_component_ids.size() < 2)
  {
    return split.separation_direction;
  }

  const std::optional<glm::dvec2> old_centroid = piece_centroid(0, 1);
  const std::optional<glm::dvec2> new_centroid = piece_centroid(1, split.split_component_ids.size());
  if (!old_centroid.has_value() || !new_centroid.has_value())
  {
    return split.separation_direction;
  }
  return new_centroid.value() - old_centroid.value();
}

// Activate infinitesimal virtual separation (or cut immediately). Formerly enqueued a SeparationEvent ramp.
void KineticDelaunay::maybeScheduleSeparationOrApplyPendingSplit(
  size_t parent_component_id, double split_time)
{
  activateInfinitesimalSeparationOrApplyCut(parent_component_id, split_time);
}

void KineticDelaunay::activateInfinitesimalSeparationOrApplyCut(size_t parent_component_id, double t)
{
  PendingBranchSplit* split = nullptr;
  {
    auto it = pending_branch_splits_.by_parent_.find(parent_component_id);
    if (it == pending_branch_splits_.by_parent_.end())
    {
      return;
    }
    split = &it->second;
  }
  if (split->split_component_ids.size() < 2)
  {
    return;
  }
  if (split->on_hiatus)
  {
    KINDS_DEBUG("activateInfinitesimalSeparation: parent_component_id=" << parent_component_id
                                                                        << " on hiatus at t=" << t << "; skipping");
    return;
  }
  if (split->infinitesimal_active)
  {
    // Already driving virtual separation; only attempt finalize (do not reset the epoch schedule).
    maybeFinalizeInfinitesimalSeparation(parent_component_id, t);
    return;
  }

  refreshPendingSplitSeparationCentroids(parent_component_id, t);

  if (pendingSplitSeamsAreConvex(parent_component_id, t))
  {
    KINDS_DEBUG("activateInfinitesimalSeparation: parent_component_id=" << parent_component_id
                                                                        << " already convex at t=" << t
                                                                        << "; applying graph split");
    ++split->infinitesimal_epoch;
    split->infinitesimal_active = false;
    applyPendingRuntimeBranchSplit(t, split->parent_runtime_branch);
    return;
  }

  ++split->infinitesimal_epoch;
  split->infinitesimal_active = true;
  split->infinitesimal_t_event = t;
  KINDS_DEBUG("activateInfinitesimalSeparation: parent_component_id=" << parent_component_id << " t=" << t
                                                                      << " epoch=" << split->infinitesimal_epoch
                                                                      << " dir=" << glm::to_string(split->separation_direction));
  // Seed the virtual event queue for the seam region (same local post-event paradigm applies afterward).
  recomputeEventsAfterInfinitesimalSeparation(parent_component_id, t, /*min_virtual_x=*/0.0);
}

bool KineticDelaunay::maybeFinalizeInfinitesimalSeparation(size_t parent_component_id, double t)
{
  auto it = pending_branch_splits_.by_parent_.find(parent_component_id);
  if (it == pending_branch_splits_.by_parent_.end())
  {
    return false;
  }
  PendingBranchSplit& split = it->second;
  if (split.on_hiatus || !split.infinitesimal_active || split.split_component_ids.size() < 2)
  {
    return false;
  }
  if (!pendingSplitSeamsAreConvex(parent_component_id, t))
  {
    return false;
  }

  KINDS_DEBUG("maybeFinalizeInfinitesimalSeparation: parent_component_id=" << parent_component_id << " t=" << t
                                                                           << " epoch=" << split.infinitesimal_epoch
                                                                           << "; applying graph split");
  ++split.infinitesimal_epoch;
  split.infinitesimal_active = false;
  applyPendingRuntimeBranchSplit(t, split.parent_runtime_branch);
  return true;
}

KineticDelaunay::ScopedInfinitesimalEventCompute::ScopedInfinitesimalEventCompute(
  KineticDelaunay& kd, size_t parent_component_id, double t, double min_virtual_x)
  : kd_(kd)
  , previous_computing_(kd.computing_infinitesimal_events_)
  , previous_min_x_(kd.infinitesimal_recompute_min_x_)
  , previous_schedule_t_(kd.infinitesimal_schedule_t_)
  , previous_parent_(kd.infinitesimal_schedule_parent_)
  , previous_epoch_(kd.infinitesimal_schedule_epoch_)
{
  const PendingBranchSplit* split = kd_.pending_branch_splits_.findByParent(parent_component_id);
  if (split == nullptr || !split->infinitesimal_active || split->on_hiatus)
  {
    return;
  }
  active_ = true;
  kd_.computing_infinitesimal_events_ = true;
  kd_.infinitesimal_recompute_min_x_ = min_virtual_x;
  kd_.infinitesimal_schedule_t_ = t;
  kd_.infinitesimal_schedule_parent_ = parent_component_id;
  kd_.infinitesimal_schedule_epoch_ = split->infinitesimal_epoch;
}

KineticDelaunay::ScopedInfinitesimalEventCompute::~ScopedInfinitesimalEventCompute()
{
  if (!active_)
  {
    return;
  }
  kd_.computing_infinitesimal_events_ = previous_computing_;
  kd_.infinitesimal_recompute_min_x_ = previous_min_x_;
  kd_.infinitesimal_schedule_t_ = previous_schedule_t_;
  kd_.infinitesimal_schedule_parent_ = previous_parent_;
  kd_.infinitesimal_schedule_epoch_ = previous_epoch_;
}

void KineticDelaunay::recomputeEventsAfterInfinitesimalSeparation(
  size_t parent_component_id, double t, double min_virtual_x)
{
  const PendingBranchSplit* split = pending_branch_splits_.findByParent(parent_component_id);
  if (split == nullptr || !split->infinitesimal_active || split->on_hiatus)
  {
    return;
  }

  std::unordered_set<size_t> affected_quads;
  std::unordered_set<size_t> affected_faces;
  collectSeparationRecomputeTargets(parent_component_id, affected_quads, affected_faces);
  growGraphSlotArrays();

  const InfinitesimalComputeContext infinitesimal { min_virtual_x, parent_component_id };
  const EventTime stamp(t, min_virtual_x);

  for (size_t quad_id : affected_quads)
  {
    if (quad_id < quadrilateral_last_updated.size())
    {
      quadrilateral_last_updated[quad_id] = stamp;
    }
    flip_event_manager_->computeEvents(t, quad_id, infinitesimal);
  }

  for (size_t face_id : affected_faces)
  {
    if (face_id < face_last_updated.size())
    {
      face_last_updated[face_id] = stamp;
    }
    const size_t he_id = graph.face(face_id).half_edges[0];
    radius_event_manager_->computeEvents(t, he_id, infinitesimal);

    const auto stamp_voronoi_crossing_invalidation = [&](size_t voronoi_vertex_id)
    {
      if (voronoi_vertex_id < crossing_data.last_crossing.size())
      {
        crossing_data.last_crossing[voronoi_vertex_id] = stamp;
      }
    };

    const auto recompute_crossing_for_voronoi_vertex = [&](size_t voronoi_vertex_id)
    {
      if (!crossing_data.isVoronoiVertexRegistered(voronoi_vertex_id))
      {
        return;
      }
      stamp_voronoi_crossing_invalidation(voronoi_vertex_id);
      crossing_event_manager_->computeEvents(t, voronoi_vertex_id, infinitesimal);
    };

    recompute_crossing_for_voronoi_vertex(face_id);
    for (size_t voronoi_vertex_id : crossing_data.getVoronoiVerticesInTri(face_id))
    {
      if (voronoi_vertex_id == face_id)
      {
        stamp_voronoi_crossing_invalidation(voronoi_vertex_id);
        continue;
      }
      recompute_crossing_for_voronoi_vertex(voronoi_vertex_id);
    }
  }
}

namespace
{
/// Sibling components sharing the seam with @p component_id: all of @p split_component_ids except itself.
std::unordered_set<size_t> partnerComponentsExcluding(
  const std::vector<size_t>& split_component_ids, size_t component_id)
{
  std::unordered_set<size_t> partners;
  for (size_t partner : split_component_ids)
  {
    if (partner != component_id)
    {
      partners.insert(partner);
    }
  }
  return partners;
}
} // namespace

bool KineticDelaunay::pendingSplitSeamsAreConvex(size_t parent_component_id, double t) const
{
  const PendingBranchSplit* split = pending_branch_splits_.findByParent(parent_component_id);
  if (split == nullptr || split->on_hiatus || split->split_component_ids.size() < 2)
  {
    return false;
  }

  for (size_t component_id : split->split_component_ids)
  {
    const std::unordered_set<size_t> partner_components
      = partnerComponentsExcluding(split->split_component_ids, component_id);
    const std::vector<size_t> start_edges = collectFutureBranchSeamStartEdges(component_id, partner_components);
    if (!isConvexPolygonOutline(collectFutureRuntimeBranchOutline(start_edges, t)))
    {
      return false;
    }
  }
  return true;
}

namespace
{
constexpr size_t kSeparationDebugEvenHalfEdge = 22;
constexpr size_t kSeparationDebugSiteA = 2;
constexpr size_t kSeparationDebugSiteB = 18;

double inCircleNumeric(const glm::dvec2& a, const glm::dvec2& b, const glm::dvec2& c, const glm::dvec2& p)
{
  const glm::dvec2 da = a - p;
  const glm::dvec2 db = b - p;
  const glm::dvec2 dc = c - p;
  const double apa = glm::dot(da, da);
  const double apb = glm::dot(db, db);
  const double apc = glm::dot(dc, dc);
  return da.x * (db.y * apc - apb * dc.y) - da.y * (db.x * apc - apb * dc.x) + apa * (db.x * dc.y - db.y * dc.x);
}

bool edgeConnectsSites(const HalfEdgeDelaunayGraph& graph, size_t even_he, size_t site_a, size_t site_b)
{
  if (!graph.isLiveHalfEdge(even_he))
  {
    return false;
  }
  const int origin = graph.halfEdge(even_he).origin;
  const int destination = graph.destination(even_he);
  return (origin == static_cast<int>(site_a) && destination == static_cast<int>(site_b))
    || (origin == static_cast<int>(site_b) && destination == static_cast<int>(site_a));
}
} // namespace // namespace

void KineticDelaunay::debugSeparationTrackedFlipProbe(
  size_t parent_component_id, double t, const char* phase, size_t even_half_edge_id) const
{
  if (!isDiagnosticsHalfEdgeIdValid(even_half_edge_id) || !isDiagnosticsHalfEdgeIdValid(even_half_edge_id ^ 1)
    || !isDiagnosticsStrandIdValid(kSeparationDebugSiteA) || !isDiagnosticsStrandIdValid(kSeparationDebugSiteB))
  {
    return;
  }

  const size_t he_id = even_half_edge_id;
  const size_t quad_id = he_id / 2;
  const size_t section = static_cast<size_t>(t);
  const double fraction = t - static_cast<double>(section);

  std::unordered_set<size_t> affected_quads;
  std::unordered_set<size_t> affected_faces;
  collectSeparationRecomputeTargets(parent_component_id, affected_quads, affected_faces);

  KINDS_DEBUG("Separation flip probe [" << phase << "] parent=" << parent_component_id << " t=" << t << " he="
                                        << he_id << "/" << (he_id ^ 1) << " quad_in_recompute_set="
                                        << (affected_quads.count(quad_id) > 0));

  if (!graph.isLiveHalfEdge(he_id))
  {
    KINDS_DEBUG("Separation flip probe: half-edge " << he_id << " is not live");
    return;
  }

  if (!edgeConnectsSites(graph, he_id, kSeparationDebugSiteA, kSeparationDebugSiteB))
  {
  const int origin = graph.halfEdge(he_id).origin;
  const int destination = graph.destination(he_id);
    KINDS_DEBUG("Separation flip probe: half-edge " << he_id << " connects sites " << origin << " and "
                                                    << destination << " (expected " << kSeparationDebugSiteA << " and "
                                                    << kSeparationDebugSiteB << ")");
  }

  const PendingBranchSplit* split = pending_branch_splits_.findByParent(parent_component_id);
  if (split != nullptr)
  {
    KINDS_DEBUG("Separation flip probe: infinitesimal_active=" << split->infinitesimal_active
                                                    << " epoch=" << split->infinitesimal_epoch
                                                    << " t_event=" << split->infinitesimal_t_event << " dir="
                                                    << glm::to_string(split->separation_direction));
  }

  for (size_t site_id : { kSeparationDebugSiteA, kSeparationDebugSiteB })
  {
    const glm::dvec2 offset = separationOffsetAt(site_id, t);
    const glm::dvec2 shifted = getPointAt(site_id, t);
    const glm::dvec2 base = shifted - offset;
    const PendingBranchSplit* site_split = activeSeparationForStrand(site_id);
    KINDS_DEBUG("Separation flip probe site " << site_id << " base=" << glm::to_string(base) << " offset="
                                             << glm::to_string(offset) << " shifted=" << glm::to_string(shifted)
                                             << " active_split=" << (site_split != nullptr) << " component="
                                             << (site_id < component_data.component_map.size()
                                                  ? component_data.component_map[site_id]
                                                  : static_cast<size_t>(-1)));
  }

  if (graph.isOnConvexBoundary(he_id) || graph.isOutsideConvexBoundary(he_id))
  {
    KINDS_DEBUG("Separation flip probe: half-edge " << he_id << " is on convex boundary; skipping inCircle checks");
    return;
  }

  const int a = graph.halfEdge(he_id).origin;
  const int b = graph.triangleOppositeVertex(he_id ^ 1);
  const int c = graph.halfEdge(he_id ^ 1).origin;
  const int d = graph.triangleOppositeVertex(he_id);
  if (a < 0 || b < 0 || c < 0 || d < 0)
  {
    KINDS_DEBUG("Separation flip probe: non-finite quad for he " << he_id);
    return;
  }

  const glm::dvec2 pa = getPointAt(static_cast<size_t>(a), t);
  const glm::dvec2 pb = getPointAt(static_cast<size_t>(b), t);
  const glm::dvec2 pc = getPointAt(static_cast<size_t>(c), t);
  const glm::dvec2 pd = getPointAt(static_cast<size_t>(d), t);
  const glm::dvec2 pa0 = getPointAt(static_cast<size_t>(a), t) - separationOffsetAt(static_cast<size_t>(a), t);
  const glm::dvec2 pb0 = getPointAt(static_cast<size_t>(b), t) - separationOffsetAt(static_cast<size_t>(b), t);
  const glm::dvec2 pc0 = getPointAt(static_cast<size_t>(c), t) - separationOffsetAt(static_cast<size_t>(c), t);
  const glm::dvec2 pd0 = getPointAt(static_cast<size_t>(d), t) - separationOffsetAt(static_cast<size_t>(d), t);

  const double in_circle_shifted = inCircleNumeric(pa, pb, pc, pd);
  const double in_circle_base = inCircleNumeric(pa0, pb0, pc0, pd0);
  const double in_circle_opposite_shifted = inCircleNumeric(pc, pb, pa, pd);
  const double in_circle_opposite_base = inCircleNumeric(pc0, pb0, pa0, pd0);

  KINDS_DEBUG("Separation flip probe quad a=" << a << " b=" << b << " c=" << c << " d=" << d
                                              << " inCircle(abc,d) shifted=" << in_circle_shifted << " base="
                                              << in_circle_base << " opposite_shifted=" << in_circle_opposite_shifted
                                              << " opposite_base=" << in_circle_opposite_base
                                              << " (positive => Delaunay violated / flip needed)");

  const FlipEventTriggerDump dump = buildFlipEventTriggerDump(*this, he_id, t);
  const double trigger_at_fraction = dump.event_trigger.degree() >= 0 ? dump.event_trigger(fraction) : 0.0;
  Polynomial trigger_copy = dump.event_trigger;
  const std::vector<double> trigger_roots
    = trigger_copy.degree() >= 0 ? const_cast<KineticDelaunay*>(this)->findEvents(trigger_copy, fraction)
                                   : std::vector<double> {};

  KINDS_DEBUG("Separation flip probe polynomial section=" << section << " fraction=" << fraction
                                                        << " trigger(fraction)=" << trigger_at_fraction
                                                        << " trigger_roots_after_fraction="
                                                        << trigger_roots.size());

  if (!trigger_roots.empty())
  {
    std::ostringstream roots;
    roots.setf(std::ios::fixed);
    roots.precision(6);
    for (size_t i = 0; i < trigger_roots.size(); ++i)
    {
      if (i > 0)
      {
        roots << ", ";
      }
      roots << (trigger_roots[i] + static_cast<double>(section));
    }
    KINDS_DEBUG("Separation flip probe scheduled flip roots at t: " << roots.str());
  }

  if (dump.site_trajectories.size() == 4)
  {
    for (size_t i = 0; i < dump.site_trajectories.size(); ++i)
    {
      const size_t site_id = static_cast<size_t>(dump.vertex_ids[i]);
      const glm::dvec2 poly_pos(dump.site_trajectories[i][0](fraction), dump.site_trajectories[i][1](fraction));
      const glm::dvec2 runtime_pos = getPointAt(site_id, t);
      const glm::dvec2 delta = poly_pos - runtime_pos;
      KINDS_DEBUG("Separation flip probe site " << site_id << " poly@fraction=" << glm::to_string(poly_pos)
                                               << " runtime@t=" << glm::to_string(runtime_pos) << " delta="
                                               << glm::to_string(delta) << " |delta|=" << glm::length(delta));
    }
  }
}

void KineticDelaunay::logSplitTransformOrthonormalityDiagnostic(size_t parent_component_id, double t) const
{
  const PendingBranchSplit* split = pending_branch_splits_.findByParent(parent_component_id);
  if (split == nullptr)
  {
    return;
  }

  const size_t tree_height = branch_trajs.getHeight();
  const size_t lower_height = static_cast<size_t>(std::floor(t));

  const auto logForHeight = [&](size_t height) {
    if (tree_height > 0 && height > tree_height)
    {
      return;
    }

    // Input branches touched by the split at this height.
    std::unordered_set<size_t> input_branches;
    for (size_t component_id : split->split_component_ids)
    {
      if (component_id >= component_data.components.size())
      {
        continue;
      }
      for (size_t strand_id : component_data.components[component_id])
      {
        if (!isDummyBoundary(strand_id))
        {
          input_branches.insert(branch_trajs.getBranchIndex(strand_id, height));
        }
      }
    }

    for (size_t input_branch : input_branches)
    {
      if (input_branch >= branch_trajs.getBranchCount(height))
      {
        continue;
      }

      const glm::dmat4& transform = branch_trajs.getTransformByHeightAndBranch(height, input_branch);
      const glm::dvec3 col0(transform[0]);
      const glm::dvec3 col1(transform[1]); // profile-plane normal axis
      const glm::dvec3 col2(transform[2]);
      const double len0 = glm::length(col0);
      const double len1 = glm::length(col1);
      const double len2 = glm::length(col2);
      const double cos02 = (len0 > 0.0 && len2 > 0.0) ? glm::dot(col0, col2) / (len0 * len2) : 0.0;
      const double cos01 = (len0 > 0.0 && len1 > 0.0) ? glm::dot(col0, col1) / (len0 * len1) : 0.0;
      const double cos12 = (len1 > 0.0 && len2 > 0.0) ? glm::dot(col1, col2) / (len1 * len2) : 0.0;
      const double det3 = glm::dot(glm::cross(col0, col1), col2); // >0 right-handed, <0 reflected
      const double inplane_len_ratio = (std::max(len0, len2) > 0.0) ? std::min(len0, len2) / std::max(len0, len2) : 0.0;

      KINDS_DEBUG("split_transform_diag t=" << t << " height=" << height << " input_branch=" << input_branch
        << " |col0|=" << len0 << " |col1|=" << len1 << " |col2|=" << len2
        << " cos(col0,col2)=" << cos02 << " cos(col0,col1)=" << cos01 << " cos(col1,col2)=" << cos12
        << " det3=" << det3 << " inplane_len_ratio=" << inplane_len_ratio
        << " inplane_orthogonal=" << (std::abs(cos02) <= 1e-9));
    }
  };

  logForHeight(lower_height);
  logForHeight(lower_height + 1);
}

void KineticDelaunay::handleSeparationEventAtTime(size_t parent_component_id, double t)
{
  // Legacy SeparationEvent entry point: redirect to instantaneous infinitesimal separation.
  activateInfinitesimalSeparationOrApplyCut(parent_component_id, t);
}

void KineticDelaunay::applyPendingRuntimeBranchSplit(double t, size_t parent_runtime_branch_id)
{
  if (parent_runtime_branch_id == RuntimeBranchData::no_branch)
  {
    return;
  }
  if (isPendingRuntimeBranchOnHiatus(parent_runtime_branch_id))
  {
    KINDS_DEBUG("applyPendingRuntimeBranchSplit: parent_runtime=" << parent_runtime_branch_id
                                                                  << " on hiatus at t=" << t << "; skipping finalize");
    return;
  }
  if (runtime_branch_data_.pendingChildBranches(parent_runtime_branch_id) == nullptr)
  {
    KINDS_DEBUG("applyPendingRuntimeBranchSplit: parent_runtime=" << parent_runtime_branch_id
                                                                  << " has no pending children; skipping");
    return;
  }

  if (callback_manager_)
  {
    callback_manager_->onBeforeComponentGraphSplit(t);
  }

  const size_t prev_face_slots = face_inside.size();
  const size_t prev_he_slots = graph.halfEdgeSlotCount();
  if (getComponentSplitPolicy() == ComponentSplitPolicy::Retriangulate)
  {
    graph.update(graph.getVertexCount(), component_data.components,
      [this, t](size_t v) { return getPointAt(v, t); });
    onGraphRetriangulated(t, prev_face_slots, prev_he_slots);
  }
  else
  {
    // Pass occurrence EventTime only for full visual-debug dumps (not --error-files alone).
    const std::optional<EventTime> branch_split_debug_time
      = isVisualDebugEnabled() ? std::optional<EventTime>(eventTimeAt(t)) : std::nullopt;
    // Cut only this parent runtime branch's pending children. Other pending-split child ids are collapsed in the
    // cut map so untargeted branches are not severed.
    const std::vector<size_t> cut_map = buildRuntimeBranchCutMapForParent(parent_runtime_branch_id);
    const HalfEdgeDelaunayGraph::RuntimeBranchSplitResult split_result = graph.applyRuntimeBranchSplit(
      cut_map, [this, t](size_t v) { return getPointAt(v, t); }, branch_split_debug_time);
    onGraphCutApplied(t, prev_face_slots, prev_he_slots, false, &split_result);
  }

  completePendingRuntimeBranchSplit(parent_runtime_branch_id, t);
  syncPrevComponentCountWithPendingSplits();
}

std::vector<size_t> KineticDelaunay::buildRuntimeBranchCutMapForParent(size_t target_parent_runtime_branch) const
{
  std::vector<size_t> cut_map = runtime_branch_data_.branch_map;
  for (size_t strand_id = 0; strand_id < cut_map.size(); ++strand_id)
  {
    const size_t branch_id = cut_map[strand_id];
    if (branch_id == RuntimeBranchData::no_branch)
    {
      continue;
    }
    if (!runtime_branch_data_.isPendingSplitChild(branch_id))
    {
      continue;
    }
    if (branch_id >= runtime_branch_data_.parent.size())
    {
      continue;
    }
    const size_t parent_branch = runtime_branch_data_.parent[branch_id];
    if (parent_branch != target_parent_runtime_branch)
    {
      // Leave other pending splits connected in the mesh until they pass their own cut check.
      cut_map[strand_id] = runtime_branch_data_.unsplitBranchId(branch_id);
    }
  }
  return cut_map;
}

void KineticDelaunay::syncPrevComponentCountWithPendingSplits()
{
  size_t uncut_extra = 0;
  for (const auto& entry : pending_branch_splits_.by_parent_)
  {
    const size_t piece_count = entry.second.split_component_ids.size();
    if (piece_count >= 2)
    {
      uncut_extra += piece_count - 1;
    }
  }
  if (component_data.components.size() >= uncut_extra)
  {
    prev_component_count = component_data.components.size() - uncut_extra;
  }
  else
  {
    prev_component_count = component_data.components.size();
  }
}

void KineticDelaunay::completePendingRuntimeBranchSplit(size_t parent_runtime_branch_id, double t)
{
  if (parent_runtime_branch_id == RuntimeBranchData::no_branch)
  {
    return;
  }

  std::vector<size_t> finalized_children;
  if (const std::vector<size_t>* children = runtime_branch_data_.pendingChildBranches(parent_runtime_branch_id))
  {
    finalized_children = *children;
  }
  runtime_branch_data_.completeSplit(parent_runtime_branch_id);

  // A runtime branch may own multiple kinetic PendingBranchSplit rows (one per split component parent).
  std::vector<size_t> kinetic_parents_to_erase;
  for (const auto& entry : pending_branch_splits_.by_parent_)
  {
    if (entry.second.parent_runtime_branch == parent_runtime_branch_id)
    {
      kinetic_parents_to_erase.push_back(entry.first);
    }
  }
  for (size_t parent_component_id : kinetic_parents_to_erase)
  {
    for (size_t strand_id = 0; strand_id < pending_branch_splits_.strand_parent_component_.size(); ++strand_id)
    {
      if (pending_branch_splits_.strand_parent_component_[strand_id] == parent_component_id)
      {
        pending_branch_splits_.strand_parent_component_[strand_id] = PendingBranchSplitState::no_parent_component;
      }
    }
    pending_branch_splits_.by_parent_.erase(parent_component_id);
  }

  // Child runtime ids were assigned at notePendingBranchSplit; they are now final. Only rebuild strand lists.
  runtime_branch_data_.rebuildBranchStrandLists();

  std::ostringstream line;
  line << "event=complete_split parent_runtime=" << parent_runtime_branch_id << " child_runtime_ids=";
  for (size_t i = 0; i < finalized_children.size(); ++i)
  {
    if (i > 0)
    {
      line << ',';
    }
    line << finalized_children[i];
  }
  logRuntimeBranchEvent(t, line.str());
}

const PendingBranchSplit* KineticDelaunay::activeSeparationForStrand(size_t strand_id) const
{
  if (strand_id >= component_data.component_map.size())
  {
    return nullptr;
  }

  const size_t component_id = component_data.component_map[strand_id];
  for (const auto& entry : pending_branch_splits_.by_parent_)
  {
    const PendingBranchSplit& split = entry.second;
    if (!split.infinitesimal_active || split.split_component_ids.size() < 2)
    {
      continue;
    }

    for (size_t i = 1; i < split.split_component_ids.size(); ++i)
    {
      if (split.split_component_ids[i] == component_id)
      {
        return &split;
      }
    }
  }

  return nullptr;
}

const PendingBranchSplit* KineticDelaunay::infinitesimalSplitForStrand(size_t strand_id) const
{
  if (strand_id >= component_data.component_map.size())
  {
    return nullptr;
  }

  const size_t component_id = component_data.component_map[strand_id];
  for (const auto& entry : pending_branch_splits_.by_parent_)
  {
    const PendingBranchSplit& split = entry.second;
    if (!split.infinitesimal_active || split.split_component_ids.size() < 2)
    {
      continue;
    }
    for (size_t piece_id : split.split_component_ids)
    {
      if (piece_id == component_id)
      {
        return &split;
      }
    }
  }
  return nullptr;
}

glm::dvec2 KineticDelaunay::separationOffsetAt(size_t strand_id, double t) const
{
  const PendingBranchSplit* split = activeSeparationForStrand(strand_id);
  if (split == nullptr || !split->infinitesimal_active)
  {
    return glm::dvec2(0.0, 0.0);
  }
  // Direction must live in the same ambient space as getPointAt(..., transform=true).
  const glm::dvec2 dir
    = computeSeparationDirection(*split, t, /*apply_reference_transform=*/true, /*shared_reference_branch=*/std::nullopt);
  return current_infinitesimal_t_ * dir;
}

Trajectory<2> KineticDelaunay::buildInfinitesimalSiteTrajectory(
  size_t strand_id, double t_event, const std::vector<size_t>& event_strand_ids) const
{
  // Match kinetic event frame policy: shared transformed frame when the trigger spans multiple input
  // branches, otherwise local support coordinates. Separation direction is computed in that same frame.
  const double event_interval_upper_bound = eventIntervalUpperBound(t_event);
  const bool use_shared_frame = !event_strand_ids.empty()
    && eventTriggerUsesSharedTransformedFrame(event_strand_ids, event_interval_upper_bound);

  std::optional<size_t> shared_reference_branch;
  bool apply_reference_transform = false;
  if (use_shared_frame)
  {
    shared_reference_branch
      = sharedReferenceBranchForEventTrigger(event_strand_ids, event_interval_upper_bound);
    apply_reference_transform = true;
  }

  const auto eval_frozen = [&](size_t sid) -> glm::dvec2
  {
    if (shared_reference_branch.has_value())
    {
      return getPointAtWithReferenceBranch(
        sid, t_event, *shared_reference_branch, apply_reference_transform, /*include_virtual_offset=*/false);
    }
    return getPointAt(sid, t_event, apply_reference_transform, /*include_virtual_offset=*/false);
  };

  const glm::dvec2 p = eval_frozen(strand_id);

  Trajectory<2> result;
  const PendingBranchSplit* child_split = activeSeparationForStrand(strand_id);
  if (child_split != nullptr)
  {
    const glm::dvec2 dir = computeSeparationDirection(
      *child_split, t_event, apply_reference_transform, shared_reference_branch);
    for (int dim = 0; dim < 2; ++dim)
    {
      result[dim] = POLYNOMIAL(p[dim] + dir[dim] * x);
    }
  }
  else
  {
    for (int dim = 0; dim < 2; ++dim)
    {
      result[dim] = POLYNOMIAL(p[dim]);
    }
  }
  return result;
}

Trajectory<2> KineticDelaunay::addSeparationOffsetToPiecePolynomial(
  const Trajectory<2>& base, size_t strand_id, size_t section, double schedule_time) const
{
  // Kinetic-time ramp removed; virtual offset is applied only via buildInfinitesimalSiteTrajectory.
  (void)strand_id;
  (void)section;
  (void)schedule_time;
  return base;
}

void KineticDelaunay::collectSeparationRecomputeTargets(size_t parent_component_id,
  std::unordered_set<size_t>& affected_quads, std::unordered_set<size_t>& affected_faces) const
{
  affected_quads.clear();
  affected_faces.clear();

  const PendingBranchSplit* split = pending_branch_splits_.findByParent(parent_component_id);
  if (split == nullptr || split->split_component_ids.empty())
  {
    return;
  }

  std::unordered_set<size_t> retained_strands;
  std::unordered_set<size_t> separated_strands;
  for (size_t strand_id : component_data.components[split->split_component_ids[0]])
  {
    retained_strands.insert(strand_id);
  }
  for (size_t i = 1; i < split->split_component_ids.size(); ++i)
  {
    for (size_t strand_id : component_data.components[split->split_component_ids[i]])
    {
      separated_strands.insert(strand_id);
    }
  }

  const auto touches_differently_shifted_pending_branches = [&](const std::vector<size_t>& strand_ids) {
    // Retained piece has no separation offset; separated pieces share the pending-split shift.
    // Only mixed sets need recomputation — a common shift leaves predicate roots invariant.
    bool has_retained = false;
    bool has_separated = false;
    for (size_t strand_id : strand_ids)
    {
      if (retained_strands.count(strand_id) > 0)
      {
        has_retained = true;
      }
      if (separated_strands.count(strand_id) > 0)
      {
        has_separated = true;
      }
    }
    return has_retained && has_separated;
  };

  for (size_t he_id : graph.liveDelaunayEdges())
  {
    if (touches_differently_shifted_pending_branches(collectFlipQuadrilateralStrandIds(graph, he_id)))
    {
      affected_quads.insert(he_id / 2);
    }
  }

  for (size_t face_id : graph.liveFaces())
  {
    const size_t he_id = graph.face(face_id).half_edges[0];
    const int a = graph.halfEdge(he_id).origin;
    const int b = graph.triangleOppositeVertex(he_id ^ 1);
    const int c = graph.halfEdge(he_id ^ 1).origin;
    std::vector<size_t> triangle_strands;
    triangle_strands.reserve(3);
    for (int vertex : { a, b, c })
    {
      if (vertex >= 0)
      {
        triangle_strands.push_back(static_cast<size_t>(vertex));
      }
    }
    if (touches_differently_shifted_pending_branches(triangle_strands))
    {
      affected_faces.insert(face_id);
    }
  }

  // Match VisualDebugHighlight::forSeparationRecompute: also include both triangles
  // adjacent to each affected flip quadrilateral.
  for (size_t quad_id : affected_quads)
  {
    const size_t flip_half_edge_id = quad_id * 2;
    if (flip_half_edge_id + 1 >= graph.halfEdgeSlotCount() || !graph.isLiveHalfEdge(flip_half_edge_id))
    {
      continue;
    }

    const auto& flip_he = graph.halfEdge(flip_half_edge_id);
    const auto& flip_twin_he = graph.halfEdge(flip_half_edge_id ^ 1);
    if (flip_he.face != static_cast<int>(-1))
    {
      affected_faces.insert(static_cast<size_t>(flip_he.face));
    }
    if (flip_twin_he.face != static_cast<int>(-1))
    {
      affected_faces.insert(static_cast<size_t>(flip_twin_he.face));
    }
  }
}

VisualDebugHighlight KineticDelaunay::buildSeparationRecomputeHighlight(size_t parent_component_id) const
{
  std::unordered_set<size_t> affected_quads;
  std::unordered_set<size_t> affected_faces;
  collectSeparationRecomputeTargets(parent_component_id, affected_quads, affected_faces);
  return VisualDebugHighlight::forSeparationRecompute(graph, affected_quads, affected_faces);
}

std::vector<HalfEdgeDelaunayGraphToSVG::SeparationOffsetSegment> KineticDelaunay::collectSeparationOffsetSegments(
  size_t parent_component_id, double t) const
{
  std::vector<HalfEdgeDelaunayGraphToSVG::SeparationOffsetSegment> segments;
  const PendingBranchSplit* split = pending_branch_splits_.findByParent(parent_component_id);
  if (split == nullptr)
  {
    return segments;
  }

  constexpr double min_offset_len2 = 1e-18;
  for (size_t i = 1; i < split->split_component_ids.size(); ++i)
  {
    const size_t component_id = split->split_component_ids[i];
    if (component_id >= component_data.components.size())
    {
      continue;
    }

    for (size_t strand_id : component_data.components[component_id])
    {
      if (isDummyBoundary(strand_id))
      {
        continue;
      }

      const glm::dvec2 offset = separationOffsetAt(strand_id, t);
      if (glm::dot(offset, offset) < min_offset_len2)
      {
        continue;
      }

      const glm::dvec2 end = getPointAt(strand_id, t);
      segments.push_back(HalfEdgeDelaunayGraphToSVG::SeparationOffsetSegment { end - offset, end });
    }
  }

  return segments;
}

std::optional<double> KineticDelaunay::separationRampEndSectionFraction(size_t strand_id, size_t section) const
{
  // Kinetic-time separation ramp removed; no mid-section ramp endpoint.
  (void)strand_id;
  (void)section;
  return std::nullopt;
}

std::optional<PendingBranchSplit> KineticDelaunay::getPendingBranchSplit(size_t parent_component_id) const
{
  const PendingBranchSplit* split = pending_branch_splits_.findByParent(parent_component_id);
  if (split == nullptr)
  {
    return std::nullopt;
  }
  return *split;
}

void KineticDelaunay::visitPendingBranchSplits(
  const std::function<void(size_t parent_component_id, const PendingBranchSplit& split)>& visitor) const
{
  for (const auto& entry : pending_branch_splits_.by_parent_)
  {
    visitor(entry.first, entry.second);
  }
}

std::vector<BoundaryPoint> KineticDelaunay::collectFutureRuntimeBranchOutline(
  const std::vector<size_t>& start_edges, double t) const
{
  if (start_edges.empty())
  {
    return {};
  }

  std::vector<bool> he_visited(graph.halfEdgeSlotCount(), false);
  std::vector<std::vector<BoundaryPoint>> seam_loops;

  for (size_t start_he_id : start_edges)
  {
    if (start_he_id >= he_visited.size() || he_visited[start_he_id])
    {
      continue;
    }
    he_visited[start_he_id] = true;

    std::vector<BoundaryPoint> seam_loop = traverseFutureBranchSeam(start_he_id, t);
    for (const BoundaryPoint& boundary_point : seam_loop)
    {
      he_visited[boundary_point.he_id] = true;
      he_visited[graph.twin(boundary_point.he_id)] = true;
    }
    if (seam_loop.size() >= 3)
    {
      seam_loops.push_back(std::move(seam_loop));
    }
  }

  if (seam_loops.empty())
  {
    return {};
  }

  const auto signed_area2 = [](const std::vector<BoundaryPoint>& loop) {
    double area2 = 0.0;
    for (size_t i = 0; i < loop.size(); ++i)
    {
      const glm::dvec2& p0 = loop[i].p;
      const glm::dvec2& p1 = loop[(i + 1) % loop.size()].p;
      area2 += p0.x * p1.y - p1.x * p0.y;
    }
    return area2;
  };

  auto best_it = seam_loops.begin();
  double best_area = std::abs(signed_area2(*best_it));
  for (auto it = std::next(seam_loops.begin()); it != seam_loops.end(); ++it)
  {
    const double area = std::abs(signed_area2(*it));
    if (area > best_area)
    {
      best_area = area;
      best_it = it;
    }
  }
  return *best_it;
}

std::vector<std::vector<BoundaryPoint>> KineticDelaunay::collectPendingSplitBranchOutlines(
  size_t parent_component_id, double t) const
{
  const std::optional<PendingBranchSplit> split = getPendingBranchSplit(parent_component_id);
  if (!split.has_value())
  {
    return {};
  }

  std::vector<std::vector<BoundaryPoint>> outlines;
  outlines.reserve(split->split_component_ids.size());
  for (size_t component_id : split->split_component_ids)
  {
    const std::unordered_set<size_t> partner_components
      = partnerComponentsExcluding(split->split_component_ids, component_id);
    const std::vector<size_t> start_edges = collectFutureBranchSeamStartEdges(component_id, partner_components);
    outlines.push_back(collectFutureRuntimeBranchOutline(start_edges, t));
  }
  return outlines;
}

void KineticDelaunay::clearPendingBranchSplits()
{
  pending_branch_splits_.clear();
  // Drop all pending-split lookup entries. Prefer @ref completePendingRuntimeBranchSplit for a single executed parent.
  // Post-split frame hold/blend state is intentionally retained until after the blend section.
  runtime_branch_data_.pending_splits.clear();
}

void KineticDelaunay::registerUpcomingPostSplitFrameTransitions(size_t section_index)
{
  const size_t tree_height = branch_trajs.getHeight();
  const size_t hold_end_height = section_index + 1;
  if (tree_height == 0 || hold_end_height >= tree_height)
  {
    return;
  }

  for (size_t component_id = 0; component_id < component_data.components.size(); ++component_id)
  {
    std::vector<size_t> live_strands;
    live_strands.reserve(component_data.components[component_id].size());
    for (size_t strand_id : component_data.components[component_id])
    {
      if (!isDummyBoundary(strand_id))
      {
        live_strands.push_back(strand_id);
      }
    }
    if (live_strands.size() < 2)
    {
      continue;
    }

    std::vector<size_t> distinct_branches_at_end;
    distinct_branches_at_end.reserve(live_strands.size());
    for (size_t strand_id : live_strands)
    {
      const size_t input_branch = branch_trajs.getBranchIndex(strand_id, hold_end_height);
      if (std::find(distinct_branches_at_end.begin(), distinct_branches_at_end.end(), input_branch)
        == distinct_branches_at_end.end())
      {
        distinct_branches_at_end.push_back(input_branch);
      }
    }
    if (distinct_branches_at_end.size() <= 1)
    {
      continue;
    }

    bool already_registered = false;
    for (size_t strand_id : live_strands)
    {
      if (const PostSplitFrameTransition* existing = postSplitFrameTransitionForStrand(strand_id);
        existing != nullptr && existing->hold_end_height == hold_end_height)
      {
        already_registered = true;
        break;
      }
    }
    if (already_registered)
    {
      continue;
    }

    // Parent/common frame at the section start; endpoints at hold_end_height map into this frame.
    const size_t common_reference_branch = minInputBranchForStrands(live_strands, section_index);

    PostSplitFrameTransition transition;
    transition.parent_component_id = component_id;
    transition.common_reference_branch = common_reference_branch;
    transition.hold_end_height = hold_end_height;
    transition.blend_section = hold_end_height;
    transition.strand_ids = std::move(live_strands);
    PostSplitFrameTransition& stored = post_split_frame_transitions_.add(std::move(transition));
    maybeWarnPostSplitFrameBlendRotation(stored);
  }
}

void KineticDelaunay::registerPostSplitFrameTransition(size_t parent_component_id, double split_time,
  const std::vector<size_t>& affected_strands)
{
  if (affected_strands.empty() || !std::isfinite(split_time) || split_time < 0.0)
  {
    return;
  }

  const size_t split_section = static_cast<size_t>(std::floor(split_time));
  const size_t hold_end_height = split_section + 1;
  const size_t tree_height = branch_trajs.getHeight();
  if (tree_height == 0 || hold_end_height >= tree_height)
  {
    // Need a full blend section [hold_end, hold_end+1] inside the tree.
    return;
  }

  for (size_t strand_id : affected_strands)
  {
    if (const PostSplitFrameTransition* existing = postSplitFrameTransitionForStrand(strand_id);
      existing != nullptr && existing->hold_end_height == hold_end_height)
    {
      // Section lookahead already registered the hold for this split window.
      return;
    }
  }

  const size_t common_reference_branch
    = minInputBranchForStrands(affected_strands, branchLookupHeightForTime(branch_trajs, split_time));

  PostSplitFrameTransition transition;
  transition.parent_component_id = parent_component_id;
  transition.common_reference_branch = common_reference_branch;
  transition.hold_end_height = hold_end_height;
  transition.blend_section = hold_end_height;
  transition.strand_ids = affected_strands;
  PostSplitFrameTransition& stored = post_split_frame_transitions_.add(std::move(transition));
  maybeWarnPostSplitFrameBlendRotation(stored);
}

const PostSplitFrameTransition* KineticDelaunay::postSplitFrameTransitionForStrand(size_t strand_id) const
{
  return post_split_frame_transitions_.findForStrand(strand_id);
}

void KineticDelaunay::maybeWarnPostSplitFrameBlendRotation(PostSplitFrameTransition& transition) const
{
  if (transition.rotation_warning_emitted)
  {
    return;
  }

  const size_t seam_height = transition.hold_end_height;
  const size_t native_height = transition.hold_end_height + 1;
  if (seam_height >= branch_trajs.getHeight() || native_height > branch_trajs.getHeight())
  {
    return;
  }
  if (transition.common_reference_branch >= branch_trajs.getBranchCount(seam_height))
  {
    return;
  }

  const glm::dmat4& common_transform
    = branch_trajs.getTransformByHeightAndBranch(seam_height, transition.common_reference_branch);

  std::unordered_set<size_t> native_branches;
  for (size_t strand_id : transition.strand_ids)
  {
    if (isDummyBoundary(strand_id))
    {
      continue;
    }
    native_branches.insert(branch_trajs.getBranchIndex(strand_id, native_height));
  }

  double max_abs_degrees = 0.0;
  size_t max_native_branch = transition.common_reference_branch;
  for (size_t native_branch : native_branches)
  {
    if (native_branch == transition.common_reference_branch)
    {
      continue;
    }
    if (native_branch >= branch_trajs.getBranchCount(native_height))
    {
      continue;
    }

    // Relative in-plane rotation: how the common-frame +x axis maps into the native frame at the blend end height.
    const glm::dmat4& native_transform
      = branch_trajs.getTransformByHeightAndBranch(native_height, native_branch);
    const PlaneProjector projector(common_transform, native_transform);
    const glm::dvec2 origin_native = projector.project(glm::dvec2(0.0, 0.0));
    const glm::dvec2 x_axis_native = projector.project(glm::dvec2(1.0, 0.0)) - origin_native;
    const double len2 = glm::dot(x_axis_native, x_axis_native);
    if (len2 <= 1e-24)
    {
      continue;
    }
    const double angle_rad = std::atan2(x_axis_native.y, x_axis_native.x);
    const double angle_deg = angle_rad * (180.0 / 3.14159265358979323846);
    if (std::abs(angle_deg) > std::abs(max_abs_degrees))
    {
      max_abs_degrees = angle_deg;
      max_native_branch = native_branch;
    }
  }

  transition.rotation_warning_emitted = true;
  if (std::abs(max_abs_degrees) > kPostSplitFrameBlendRotationWarnDegrees)
  {
    KINDS_WARNING("Post-split frame blend rotation "
      << max_abs_degrees << " deg exceeds threshold " << kPostSplitFrameBlendRotationWarnDegrees
      << " deg (common_branch=" << transition.common_reference_branch << ", native_branch=" << max_native_branch
      << ", hold_end_height=" << transition.hold_end_height << ", blend_section=" << transition.blend_section
      << ", parent_component=" << transition.parent_component_id << ")");
  }
}

size_t KineticDelaunay::minInputBranchForComponent(size_t component_id, size_t branch_lookup_height) const
{
  if (component_id >= component_data.components.size())
  {
    return 0;
  }

  return minInputBranchForStrands(component_data.components[component_id], branch_lookup_height);
}

size_t KineticDelaunay::minInputBranchForStrands(
  const std::vector<size_t>& strand_ids, size_t branch_lookup_height) const
{
  size_t min_branch = std::numeric_limits<size_t>::max();
  for (size_t strand_id : strand_ids)
  {
    if (isDummyBoundary(strand_id))
    {
      continue;
    }
    const size_t input_branch = branch_trajs.getBranchIndex(strand_id, branch_lookup_height);
    min_branch = std::min(min_branch, input_branch);
  }

  if (min_branch == std::numeric_limits<size_t>::max())
  {
    return 0;
  }
  return min_branch;
}

size_t KineticDelaunay::representativeStrandIdForRuntimeBranch(size_t strand_id) const
{
  if (strand_id == static_cast<size_t>(-1))
  {
    throw std::runtime_error("representativeStrandIdForRuntimeBranch: invalid strand id.");
  }

  const size_t runtime_branch_id = getRuntimeBranchIdForStrand(strand_id);
  if (runtime_branch_id == static_cast<size_t>(-1))
  {
    return strand_id;
  }

  std::optional<size_t> representative;
  for (size_t candidate = 0; candidate < runtime_branch_data_.branch_map.size(); ++candidate)
  {
    if (isDummyBoundary(candidate))
    {
      continue;
    }
    if (runtime_branch_data_.branch_map[candidate] != runtime_branch_id)
    {
      continue;
    }
    if (!representative.has_value() || candidate < representative.value())
    {
      representative = candidate;
    }
  }

  return representative.value_or(strand_id);
}

namespace
{
size_t runtimeBranchSectionIndex(const KineticDelaunay& kd, double t)
{
  const size_t height = kd.getStrandTree().getHeight();
  if (height == 0)
  {
    return 0;
  }

  size_t section_index = static_cast<size_t>(std::ceil(t));
  if (section_index >= height)
  {
    section_index = height - 1;
  }
  return section_index;
}
} // namespace

size_t KineticDelaunay::allocateRuntimeBranch() { return runtime_branch_data_.allocate(); }

void KineticDelaunay::markRuntimeBranchDead(size_t runtime_branch_id, double t, const char* reason,
  std::optional<size_t> input_branch_id)
{
  if (runtime_branch_id < runtime_branch_data_.alive.size())
  {
    runtime_branch_data_.alive[runtime_branch_id] = false;
  }

  std::ostringstream line;
  line << "event=retire id=" << runtime_branch_id << " reason=" << reason;
  if (input_branch_id.has_value())
  {
    line << " input_branch=" << input_branch_id.value();
  }
  logRuntimeBranchEvent(t, line.str());
}

void KineticDelaunay::updateRuntimeBranchMapFromInputBranches(double t)
{
  const size_t vertex_count = graph.getVertexCount();
  runtime_branch_data_.branch_map.resize(vertex_count);

  const size_t section_index = runtimeBranchSectionIndex(*this, t);
  std::vector<size_t> input_branch_runtime;
  constexpr size_t no_runtime_branch = static_cast<size_t>(-1);
  std::unordered_map<size_t, size_t> input_branch_strand_count;

  for (size_t strand_id = 0; strand_id < vertex_count; ++strand_id)
  {
    if (isDummyBoundary(strand_id))
    {
      continue;
    }

    const size_t input_branch = branch_trajs.getBranchIndex(strand_id, section_index);
    if (input_branch >= input_branch_runtime.size())
    {
      input_branch_runtime.resize(input_branch + 1, no_runtime_branch);
    }
    if (input_branch_runtime[input_branch] == no_runtime_branch)
    {
      input_branch_runtime[input_branch] = allocateRuntimeBranch();
      std::ostringstream line;
      line << "event=create id=" << input_branch_runtime[input_branch] << " reason=init_input_branch"
           << " input_branch=" << input_branch << " section=" << section_index;
      logRuntimeBranchEvent(t, line.str());
    }
    runtime_branch_data_.branch_map[strand_id] = input_branch_runtime[input_branch];
    ++input_branch_strand_count[input_branch];
  }

  for (const auto& [input_branch, strand_count] : input_branch_strand_count)
  {
    if (input_branch >= input_branch_runtime.size() || input_branch_runtime[input_branch] == no_runtime_branch)
    {
      continue;
    }
    std::ostringstream line;
    line << "event=assign id=" << input_branch_runtime[input_branch] << " reason=init_input_branch"
         << " input_branch=" << input_branch << " section=" << section_index << " strand_count=" << strand_count;
    logRuntimeBranchEvent(t, line.str());
  }

  runtime_branch_data_.rebuildBranchStrandLists();
}

namespace
{
std::string formatIdList(const std::vector<size_t>& ids)
{
  std::ostringstream out;
  for (size_t i = 0; i < ids.size(); ++i)
  {
    if (i > 0)
    {
      out << ',';
    }
    out << ids[i];
  }
  return out.str();
}

std::string formatStrandSample(const std::vector<size_t>& strands, size_t max_count = 8)
{
  std::ostringstream out;
  const size_t show = std::min(strands.size(), max_count);
  for (size_t i = 0; i < show; ++i)
  {
    if (i > 0)
    {
      out << ',';
    }
    out << strands[i];
  }
  if (strands.size() > show)
  {
    out << ",...";
  }
  return out.str();
}

std::vector<size_t> distinctSortedIds(const std::vector<size_t>& values)
{
  std::vector<size_t> ids = values;
  std::sort(ids.begin(), ids.end());
  ids.erase(std::unique(ids.begin(), ids.end()), ids.end());
  return ids;
}
} // namespace

void KineticDelaunay::updateRuntimeBranchMapFromLiveGraph(double t, const char* trigger)
{
  const size_t vertex_count = graph.getVertexCount();
  runtime_branch_data_.branch_map.resize(vertex_count);

  const size_t section_index = runtimeBranchSectionIndex(*this, t);
  const std::vector<std::vector<size_t>> graph_components = extractGraphConnectedComponents();

  for (size_t component_index = 0; component_index < graph_components.size(); ++component_index)
  {
    const auto& component = graph_components[component_index];
    if (component.empty())
    {
      continue;
    }

    std::optional<size_t> component_runtime_id;
    std::vector<size_t> prior_runtime_ids;
    std::vector<size_t> input_branch_ids;
    size_t real_strand_count = 0;

    for (size_t strand_id : component)
    {
      if (isDummyBoundary(strand_id))
      {
        continue;
      }

      ++real_strand_count;
      if (strand_id < runtime_branch_data_.branch_map.size())
      {
        prior_runtime_ids.push_back(runtime_branch_data_.branch_map[strand_id]);
      }
      input_branch_ids.push_back(branch_trajs.getBranchIndex(strand_id, section_index));

      if (strand_id >= runtime_branch_data_.branch_map.size())
      {
        continue;
      }

      const size_t runtime_id = runtime_branch_data_.branch_map[strand_id];
      if (!component_runtime_id.has_value())
      {
        component_runtime_id = runtime_id;
      }
      else if (component_runtime_id.value() != runtime_id)
      {
        component_runtime_id = std::nullopt;
        break;
      }
    }

    size_t runtime_branch = 0;
    bool reused_existing = false;
    // If every strand in this graph component already shares one designated runtime branch id (e.g. a pending-split
    // child that spanned multiple kinetic pieces), keep that id across disconnected pieces. Do not allocate a fresh
    // runtime branch per graph component.
    if (component_runtime_id.has_value() && component_runtime_id.value() < runtime_branch_data_.alive.size()
      && runtime_branch_data_.alive[component_runtime_id.value()])
    {
      runtime_branch = component_runtime_id.value();
      reused_existing = true;
    }
    else
    {
      runtime_branch = allocateRuntimeBranch();
    }

    for (size_t strand_id : component)
    {
      if (!isDummyBoundary(strand_id))
      {
        runtime_branch_data_.branch_map[strand_id] = runtime_branch;
      }
    }

    const auto sorted_prior_runtime_ids = distinctSortedIds(prior_runtime_ids);
    const auto sorted_input_branch_ids = distinctSortedIds(input_branch_ids);

    std::ostringstream line;
    if (!reused_existing)
    {
      line << "event=create id=" << runtime_branch << " reason=graph_split trigger=" << trigger
           << " graph_component=" << component_index << " strand_count=" << real_strand_count
           << " prior_runtime_ids=" << formatIdList(sorted_prior_runtime_ids)
           << " input_branches=" << formatIdList(sorted_input_branch_ids)
           << " strands=" << formatStrandSample(component);
    }
    else
    {
      line << "event=assign id=" << runtime_branch << " reason=graph_split_reuse trigger=" << trigger
           << " graph_component=" << component_index << " strand_count=" << real_strand_count
           << " prior_runtime_ids=" << formatIdList(sorted_prior_runtime_ids)
           << " input_branches=" << formatIdList(sorted_input_branch_ids)
           << " strands=" << formatStrandSample(component);
    }
    logRuntimeBranchEvent(t, line.str());
  }

  runtime_branch_data_.rebuildBranchStrandLists();
}

void KineticDelaunay::validateFinishedInputBranchMatchesRuntime(size_t section, size_t input_branch_id) const
{
  const auto& branch_strands = branch_trajs.getStrandsByBranch(section, input_branch_id);

  std::optional<size_t> runtime_branch_id;
  for (size_t strand_id : branch_strands)
  {
    if (isDummyBoundary(strand_id))
    {
      continue;
    }

    if (strand_id >= runtime_branch_data_.branch_map.size())
    {
      throw std::runtime_error("Finished input branch validation: runtime branch map is out of date");
    }

    if (!runtime_branch_id.has_value())
    {
      runtime_branch_id = runtime_branch_data_.branch_map[strand_id];
    }
    else if (runtime_branch_id.value() != runtime_branch_data_.branch_map[strand_id])
    {
      throw std::runtime_error("Finished input branch " + std::to_string(input_branch_id)
        + " spans multiple runtime branches at section " + std::to_string(section));
    }
  }

  if (!runtime_branch_id.has_value())
  {
    return;
  }

  for (size_t face_id : graph.liveFaces())
  {
    const std::array<int, 3> tri_vertices = graph.getTriangleVertexIndices(face_id);
    bool touches_finished_branch = false;

    for (int vertex : tri_vertices)
    {
      if (vertex < 0 || isDummyBoundary(static_cast<size_t>(vertex)))
      {
        continue;
      }
      if (branch_trajs.getBranchIndex(static_cast<size_t>(vertex), section) == input_branch_id)
      {
        touches_finished_branch = true;
        break;
      }
    }

    if (!touches_finished_branch)
    {
      continue;
    }

    if (getRuntimeBranchIdForFace(face_id) != runtime_branch_id.value())
    {
      throw std::runtime_error("Finished input branch " + std::to_string(input_branch_id)
        + " does not match a single runtime branch at section " + std::to_string(section));
    }

    for (int vertex : tri_vertices)
    {
      if (vertex < 0 || isDummyBoundary(static_cast<size_t>(vertex)))
      {
        continue;
      }
      if (branch_trajs.getBranchIndex(static_cast<size_t>(vertex), section) != input_branch_id)
      {
        throw std::runtime_error("Finished input branch " + std::to_string(input_branch_id)
          + " still shares a live triangle with input branch "
          + std::to_string(branch_trajs.getBranchIndex(static_cast<size_t>(vertex), section)) + " at section "
          + std::to_string(section));
      }
    }
  }

  for (size_t vertex = 0; vertex < graph.getVertexCount(); ++vertex)
  {
    if (isDummyBoundary(vertex) || vertex >= runtime_branch_data_.branch_map.size())
    {
      continue;
    }

    if (runtime_branch_data_.branch_map[vertex] != runtime_branch_id.value())
    {
      continue;
    }

    if (graph.incidentEdgesBegin(vertex) == graph.incidentEdgesEnd(vertex))
    {
      continue;
    }

    if (branch_trajs.getBranchIndex(vertex, section) != input_branch_id)
    {
      throw std::runtime_error("Finished input branch " + std::to_string(input_branch_id)
        + " does not match runtime branch " + std::to_string(runtime_branch_id.value()) + " at section "
        + std::to_string(section));
    }
  }
}

std::vector<size_t> KineticDelaunay::inputBranchesFinishingAtSection(double t) const
{
  const size_t section = static_cast<size_t>(t);
  std::vector<size_t> finished;
  if (section >= branch_trajs.getStrandsByBranchId().size())
  {
    return finished;
  }

  const auto& branches_at_section = branch_trajs.getStrandBranchesByHeight(section);
  for (size_t input_branch_id = 0; input_branch_id < branches_at_section.size(); ++input_branch_id)
  {
    const auto& branch_strands = branches_at_section[input_branch_id];
    if (branch_strands.empty())
    {
      continue;
    }

    bool any_real_strand = false;
    bool all_have_next_piece = true;
    bool none_have_next_piece = true;

    for (size_t strand_id : branch_strands)
    {
      if (isDummyBoundary(strand_id))
      {
        continue;
      }

      any_real_strand = true;
      const bool has_next_piece = branch_trajs.getSupportPoints(strand_id).size() > section + 1;
      all_have_next_piece = all_have_next_piece && has_next_piece;
      none_have_next_piece = none_have_next_piece && !has_next_piece;
    }

    if (!any_real_strand || all_have_next_piece)
    {
      continue;
    }

    if (!none_have_next_piece)
    {
      throw std::runtime_error("Input branch " + std::to_string(input_branch_id)
        + " has inconsistent strand heights at section " + std::to_string(section));
    }

    finished.push_back(input_branch_id);
  }

  return finished;
}

void KineticDelaunay::retireFinishedInputBranches(double t)
{
  const size_t section = static_cast<size_t>(t);
  if (section >= branch_trajs.getStrandsByBranchId().size())
  {
    return;
  }

  std::vector<bool> dead_vertex_mask(graph.getVertexCount(), false);
  bool has_finished_branch = false;

  for (size_t input_branch_id : inputBranchesFinishingAtSection(t))
  {
    validateFinishedInputBranchMatchesRuntime(section, input_branch_id);
    has_finished_branch = true;

    const auto& branch_strands = branch_trajs.getStrandsByBranch(section, input_branch_id);

    for (size_t strand_id : branch_strands)
    {
      if (!isDummyBoundary(strand_id) && strand_id < runtime_branch_data_.branch_map.size())
      {
        markRuntimeBranchDead(runtime_branch_data_.branch_map[strand_id], t, "input_branch_finished", input_branch_id);
        break;
      }
    }

    for (size_t strand_id : branch_strands)
    {
      if (!isDummyBoundary(strand_id))
      {
        dead_vertex_mask[strand_id] = true;
        if (strand_id < runtime_branch_data_.branch_map.size())
        {
          runtime_branch_data_.branch_map[strand_id] = static_cast<size_t>(-1);
        }
      }
    }
  }

  if (!has_finished_branch)
  {
    return;
  }

  runtime_branch_data_.rebuildBranchStrandLists();

  const size_t prev_face_slots = face_inside.size();
  const size_t prev_he_slots = graph.halfEdgeSlotCount();
  graph.killVertices(dead_vertex_mask);
  logRuntimeBranchEvent(t, "event=skip_remap reason=input_branch_vertex_kill");
  onGraphCutApplied(t, prev_face_slots, prev_he_slots, false);
}

KineticDelaunay::CrossingData& kinDS::KineticDelaunay::getCrossingDataMutable() { return crossing_data; }

const KineticDelaunay::CrossingData& kinDS::KineticDelaunay::getCrossingData() const { return crossing_data; }

void KineticDelaunay::validateCrossingIntersectionInvariants(const char* context, double t) const
{
  (void)context;
  (void)t;
}

void KineticDelaunay::validateVoronoiVertexIteratorInvariants(const char* context, double t) const
{
  crossing_data.validateVoronoiVertexIteratorInvariants(context, graph, this, t);
}

void KineticDelaunay::validateSitesInsideConvexHull(const char* context, double t) const
{
  if (!sites_inside_convex_hull_check_enabled_)
  {
    return;
  }

  constexpr double kHullContainmentEps = 1e-9;
  const auto cross2 = [](const glm::dvec2& u, const glm::dvec2& v) { return u.x * v.y - u.y * v.x; };
  const auto signed_area2 = [&](const std::vector<glm::dvec2>& poly)
  {
    double area2 = 0.0;
    for (size_t i = 0; i < poly.size(); ++i)
    {
      const glm::dvec2& p0 = poly[i];
      const glm::dvec2& p1 = poly[(i + 1) % poly.size()];
      area2 += cross2(p0, p1);
    }
    return area2;
  };
  // Convex hull edges from the graph are oriented so the infinite face is to the left of the outside half-edge
  // walk; for containment we want the finite component on the left, so reverse if needed.
  const auto point_inside_or_on_convex = [&](const glm::dvec2& p, const std::vector<glm::dvec2>& poly_ccw) -> bool
  {
    for (size_t i = 0; i < poly_ccw.size(); ++i)
    {
      const glm::dvec2& a = poly_ccw[i];
      const glm::dvec2& b = poly_ccw[(i + 1) % poly_ccw.size()];
      if (cross2(b - a, p - a) < -kHullContainmentEps)
      {
        return false;
      }
    }
    return true;
  };

  const std::vector<std::vector<size_t>> components = extractGraphConnectedComponents();
  for (size_t component_index = 0; component_index < components.size(); ++component_index)
  {
    const std::vector<size_t>& component = components[component_index];
    if (component.size() < 3)
    {
      continue;
    }

    size_t start_he = static_cast<size_t>(-1);
    for (size_t vertex_id : component)
    {
      for (auto it = graph.incidentEdgesBegin(vertex_id); it != graph.incidentEdgesEnd(vertex_id); ++it)
      {
        const size_t he_id = *it;
        if (graph.isOnConvexBoundaryOutside(he_id))
        {
          start_he = he_id;
          break;
        }
      }
      if (start_he != static_cast<size_t>(-1))
      {
        break;
      }
    }
    if (start_he == static_cast<size_t>(-1))
    {
      KINDS_WARNING("Sites-in-hull check (" << (context ? context : "?") << "): component " << component_index
                                           << " at t=" << t << " has no outside convex-hull half-edge");
      continue;
    }

    std::vector<size_t> hull_vertex_ids;
    size_t he_id = start_he;
    do
    {
      const int origin = graph.halfEdge(he_id).origin;
      if (origin < 0)
      {
        KINDS_WARNING("Sites-in-hull check (" << (context ? context : "?") << "): infinite origin on hull he="
                                             << he_id << " at t=" << t);
        hull_vertex_ids.clear();
        break;
      }
      hull_vertex_ids.push_back(static_cast<size_t>(origin));
      he_id = graph.nextOnConvexBoundaryId(he_id);
      if (!graph.isLiveHalfEdge(he_id) || !graph.isOnConvexBoundaryOutside(he_id))
      {
        KINDS_WARNING("Sites-in-hull check (" << (context ? context : "?") << "): broken hull walk from he="
                                             << start_he << " at t=" << t);
        hull_vertex_ids.clear();
        break;
      }
    } while (he_id != start_he);

    if (hull_vertex_ids.size() < 3)
    {
      continue;
    }

    const size_t reference_branch = getSharedReferenceBranchForStrands(component, t);
    std::vector<glm::dvec2> hull_poly;
    hull_poly.reserve(hull_vertex_ids.size());
    for (size_t vertex_id : hull_vertex_ids)
    {
      hull_poly.push_back(getPointInDelaunaySpace(vertex_id, t, reference_branch));
    }
    if (signed_area2(hull_poly) < 0.0)
    {
      std::reverse(hull_poly.begin(), hull_poly.end());
    }

    for (size_t site_id : component)
    {
      const glm::dvec2 site_p = getPointInDelaunaySpace(site_id, t, reference_branch);
      if (point_inside_or_on_convex(site_p, hull_poly))
      {
        continue;
      }

      double worst_cross = 0.0;
      for (size_t i = 0; i < hull_poly.size(); ++i)
      {
        const glm::dvec2& a = hull_poly[i];
        const glm::dvec2& b = hull_poly[(i + 1) % hull_poly.size()];
        worst_cross = std::min(worst_cross, cross2(b - a, site_p - a));
      }

      KINDS_WARNING("Sites-in-hull FAIL (" << (context ? context : "?") << "): site " << site_id
                                          << " outside graph convex hull of component " << component_index
                                          << " at t=" << std::setprecision(17) << t << " pos=(" << site_p.x << ","
                                          << site_p.y << ") worst_edge_cross=" << worst_cross
                                          << " (possible incorrect CCW flip)");
    }
  }
}

namespace
{
glm::dvec2 intersect_segments_2d_closing(
  const glm::dvec2& p, const glm::dvec2& p2, const glm::dvec2& q, const glm::dvec2& q2)
{
  glm::dvec2 r = p2 - p;
  glm::dvec2 s = q2 - q;
  const double rxs = r.x * s.y - r.y * s.x;
  const double qpxr = (q - p).x * r.y - (q - p).y * r.x;
  if (rxs == 0.0 && qpxr == 0.0)
  {
    return glm::dvec2(std::numeric_limits<double>::quiet_NaN());
  }
  if (rxs == 0.0 && qpxr != 0.0)
  {
    return glm::dvec2(std::numeric_limits<double>::infinity());
  }
  const double tt = ((q - p).x * s.y - (q - p).y * s.x) / rxs;
  return p + tt * r;
}
} // namespace

std::string kinDS::formatCrossingIntersectionForLog(
  const KineticDelaunay& kd, std::optional<KineticDelaunay::CrossingData::EdgeIntersectionRef> intersection)
{
  if (!intersection.has_value())
  {
    return "null";
  }
  const auto& crossing_data = kd.getCrossingData();
  const auto& ir = *intersection.value();
  const size_t d_edge = ir.delaunay_edge_id;
  const size_t v_edge = ir.voronoi_edge_id;

  auto fmt_fallback = [&]() -> std::string
  {
    return std::string("V?xD? v_edge=") + std::to_string(v_edge) + " d_edge=" + std::to_string(d_edge)
      + " tD=" + std::to_string(ir.delaunay_edge_param);
  };

  if (d_edge >= crossing_data.delaunay_edge_intersections.size()
    || v_edge >= crossing_data.voronoi_edge_intersections.size())
  {
    return fmt_fallback();
  }

  size_t d_list_idx = 0;
  bool found_on_d = false;
  for (auto d_it = crossing_data.delaunay_edge_intersections[d_edge].begin();
    d_it != crossing_data.delaunay_edge_intersections[d_edge].end(); ++d_it, ++d_list_idx)
  {
    const KineticDelaunay::CrossingData::EdgeIntersectionRef ref = *d_it;
    if (&(*ref) == &ir)
    {
      found_on_d = true;
      break;
    }
  }
  if (!found_on_d)
  {
    return fmt_fallback();
  }

  size_t v_list_idx = 0;
  bool found_on_v = false;
  for (auto v_it = crossing_data.voronoi_edge_intersections[v_edge].begin();
    v_it != crossing_data.voronoi_edge_intersections[v_edge].end(); ++v_it, ++v_list_idx)
  {
    const KineticDelaunay::CrossingData::EdgeIntersectionRef ref = *v_it;
    if (&(*ref) == &ir)
    {
      found_on_v = true;
      break;
    }
  }
  if (!found_on_v)
  {
    return fmt_fallback();
  }

  return std::string("V") + std::to_string(v_list_idx) + "xD" + std::to_string(d_list_idx) + " v_edge="
    + std::to_string(v_edge) + " d_edge=" + std::to_string(d_edge) + " tD=" + std::to_string(ir.delaunay_edge_param);
}

bool kinDS::tryComputeCrossingIntersectionPosition2D(const KineticDelaunay& kd,
  std::optional<KineticDelaunay::CrossingData::EdgeIntersectionRef> intersection, double t, glm::dvec2& out_xy,
  bool apply_reference_transform, bool include_virtual_offset)
{
  (void)apply_reference_transform;
  (void)include_virtual_offset;
  if (!intersection.has_value())
  {
    return false;
  }
  const glm::dvec3 pos = getCrossingCoordsInDelaunaySpace(const_cast<KineticDelaunay&>(kd), intersection.value(), t);
  if (!std::isfinite(pos.x) || !std::isfinite(pos.y))
  {
    return false;
  }
  out_xy = glm::dvec2(pos.x, pos.y);
  return true;
}

void kinDS::validateClosingCapCrossingRef([[maybe_unused]] const KineticDelaunay& kd, const char* context_msg,
  std::optional<KineticDelaunay::CrossingData::EdgeIntersectionRef> intersection, size_t expected_voronoi_edge_id,
  int delaunay_half_edge_id)
{
  if (!intersection.has_value())
  {
    return;
  }
  const auto& ir = *intersection.value();
  if (expected_voronoi_edge_id != static_cast<size_t>(-1) && ir.voronoi_edge_id != expected_voronoi_edge_id)
  {
    KINDS_ERROR("validateClosingCapCrossingRef ("
      << context_msg << "): crossing v_edge=" << ir.voronoi_edge_id << " but expected closing_voronoi_edge_id="
      << expected_voronoi_edge_id << " (d_edge=" << ir.delaunay_edge_id << ")");
  }
  if (delaunay_half_edge_id >= 0)
  {
    const size_t ud = static_cast<size_t>(delaunay_half_edge_id) / 2;
    if (ir.delaunay_edge_id != ud)
    {
      KINDS_ERROR("validateClosingCapCrossingRef (" << context_msg << "): crossing d_edge=" << ir.delaunay_edge_id
                                                    << " does not match delaunay_half_edge_id " << delaunay_half_edge_id
                                                    << " (undirected " << ud << ")");
    }
  }
}

std::vector<HalfEdgeDelaunayGraphToSVG::IntersectionDebugInfo> KineticDelaunay::getCrossingIntersectionDebugData() const
{
  std::vector<HalfEdgeDelaunayGraphToSVG::IntersectionDebugInfo> result;
  result.reserve(crossing_data.edge_intersections.size());

  for (size_t d_edge_id = 0; d_edge_id < crossing_data.delaunay_edge_intersections.size(); ++d_edge_id)
  {
    const auto& d_list = crossing_data.delaunay_edge_intersections[d_edge_id];
    size_t d_index = 0;
    for (auto d_it = d_list.begin(); d_it != d_list.end(); ++d_it, ++d_index)
    {
      CrossingData::EdgeIntersectionRef ref = *d_it;
      const auto& v_list = crossing_data.voronoi_edge_intersections[ref->voronoi_edge_id];
      size_t v_index = 0;
      for (auto v_it = v_list.begin(); v_it != v_list.end(); ++v_it, ++v_index)
      {
        if (*v_it == ref)
        {
          HalfEdgeDelaunayGraphToSVG::IntersectionDebugInfo info;
          info.delaunay_edge_id = ref->delaunay_edge_id;
          info.voronoi_edge_id = ref->voronoi_edge_id;
          info.delaunay_list_index = d_index;
          info.voronoi_list_index = v_index;
          info.prev_segment_mesh_pair_index = ref->prev_segment_mesh_pair_index;
          info.next_segment_mesh_pair_index = ref->next_segment_mesh_pair_index;
          info.delaunay_edge_param = ref->delaunay_edge_param;
          result.push_back(info);
          break;
        }
      }
    }
  }

  return result;
}

bool KineticDelaunay::computeBoundaryOnTheFly() const { return on_the_fly_boundary; }

glm::dvec3 lineToHomegeneous(const glm::dvec2& p, const glm::dvec2& dir)
{
  glm::dvec2 normal = glm::dvec2(-dir.y, dir.x);
  return glm::dvec3(normal, -glm::dot(normal, p));
}

std::pair<double, double> segmentIntersectionParameters(
  const glm::dvec2& p0, const glm::dvec2& p1, const glm::dvec2& q0, const glm::dvec2& q1)
{
  auto cross2D = [](const glm::dvec2& a, const glm::dvec2& b) { return a.x * b.y - a.y * b.x; };

  const glm::dvec2 p = p0;
  const glm::dvec2 r = p1 - p0;
  const glm::dvec2 q = q0;
  const glm::dvec2 s = q1 - q0;

  const double rxs = glm::cross(r, s);

  // Parallel (or nearly so): treat as no proper intersection
  if (std::abs(rxs) < 1e-12)
  {
    return std::pair<double, double>(std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity());
  }

  const double t = cross2D(q - p, s) / rxs;
  const double u = cross2D(q - p, r) / rxs;

  return std::pair<double, double>(t, u);
}

bool segmentRayIntersection(
  const glm::dvec2& segment_p0, const glm::dvec2& segment_p1, const glm::dvec2& ray_origin, const glm::dvec2& ray_dir)
{
  auto [t, u] = segmentIntersectionParameters(segment_p0, segment_p1, ray_origin, ray_origin + ray_dir);
  return t >= 0.0 && u >= 0.0 && u <= 1.0;
}

std::pair<glm::dvec2, glm::dvec2> KineticDelaunay::computeAngularBisector(size_t he_id, double t) const
{
  int a = graph.halfEdge(he_id).origin;
  int b = graph.destination(he_id);
  size_t finite_he_id = graph.halfEdge(he_id).next;
  size_t c, c_prime;

  if (a == -1)
  {
    a = b;
    c = graph.destination(finite_he_id);
    assert(c != size_t(-1));
    size_t prev_he_id = graph.prevOnConvexBoundaryId(finite_he_id);
    c_prime = graph.halfEdge(prev_he_id).origin;
    assert(c_prime != size_t(-1));
  }
  else
  {
    finite_he_id = graph.halfEdge(finite_he_id).next;
    c = graph.halfEdge(finite_he_id).origin;
    assert(c != size_t(-1));
    size_t next_he_id = graph.nextOnConvexBoundaryId(finite_he_id);
    c_prime = graph.destination(next_he_id);
    assert(c_prime != size_t(-1));
  }

  glm::dvec2 p_a = getPointAt(a, t);
  glm::dvec2 p_c = getPointAt(c, t);
  glm::dvec2 p_c_prime = getPointAt(c_prime, t);

  glm::dvec2 angular_bisector_direction = glm::normalize(p_c_prime - p_a) + glm::normalize(p_c - p_a);

  // Fallback: If the direction is zero, we can just use the normal of the segment.
  if (glm::length(angular_bisector_direction) < 1e-12)
  {
    glm::dvec2 tangent = glm::normalize(p_c_prime - p_c);
    angular_bisector_direction = glm::dvec2(-tangent.y, tangent.x);
    // TODO: I'm not sure if the sign is correct here.
  }

  // Since we are on the convex hull, the direction points inward, so we need to negate it.
  return std::pair(p_a, -angular_bisector_direction);
};

std::pair<double, double> KineticDelaunay::delaunayVoronoiEdgeIntersection(
  size_t delaunay_edge_id, size_t voronoi_edge_id, double t) const
{

  glm::dvec2 start_point = computeVoronoiVertexClampedInfinity(voronoi_edge_id * 2, t);
  glm::dvec2 destination = computeVoronoiVertexClampedInfinity(voronoi_edge_id * 2 + 1, t);

  if (graph.isInfinite(2 * delaunay_edge_id))
  {
    auto ray = computeAngularBisector(2 * delaunay_edge_id, t);

    return segmentIntersectionParameters(ray.first, ray.first + ray.second, start_point, destination);
  }
  else
  {
    glm::dvec2 edge_start = getPointAt(static_cast<size_t>(graph.halfEdge(2 * delaunay_edge_id).origin), t);
    glm::dvec2 edge_end = getPointAt(static_cast<size_t>(graph.destination(2 * delaunay_edge_id)), t);
    return segmentIntersectionParameters(edge_start, edge_end, start_point, destination);
  }
}

double KineticDelaunay::delaunayVoronoiEdgeIntersectionParameter(
  size_t delaunay_edge_id, size_t voronoi_edge_id, double t) const
{
  return delaunayVoronoiEdgeIntersection(delaunay_edge_id, voronoi_edge_id, t).first;
}

namespace
{
glm::dvec3 crossingCoordsFromDelaunayEdgeInterpolation(const KineticDelaunay& kd,
  const KineticDelaunay::CrossingData::VoronoiDelaunayEdgeIntersection& intersection, double t,
  const std::function<glm::dvec2(size_t, double)>& site_xy_at)
{
  const size_t d_he0 = 2 * intersection.delaunay_edge_id;
  const size_t d_he1 = d_he0 + 1;
  const auto& graph = kd.getGraph();
  if (d_he1 >= graph.halfEdgeSlotCount())
  {
    return glm::dvec3(std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN(), t);
  }

  const int a = graph.halfEdge(d_he0).origin;
  const int b = graph.halfEdge(d_he1).origin;
  if (a < 0 || b < 0)
  {
    return glm::dvec3(std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN(), t);
  }

  const double param = intersection.delaunay_edge_param;
  if (!std::isfinite(param))
  {
    return glm::dvec3(std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN(), t);
  }

  const glm::dvec2 pa = site_xy_at(static_cast<size_t>(a), t);
  const glm::dvec2 pb = site_xy_at(static_cast<size_t>(b), t);
  const glm::dvec2 xy = pa * (1.0 - param) + pb * param;
  if (!std::isfinite(xy.x) || !std::isfinite(xy.y))
  {
    return glm::dvec3(std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN(), t);
  }
  return glm::dvec3(xy.x, xy.y, t);
}

/// Warning-only range check. Finite edges: [0,1]. Infinite boundary rays: [0,+inf).
bool isDelaunayEdgeParamInExpectedWarningRange(
  const KineticDelaunay* kd, size_t delaunay_edge_id, double param, double eps)
{
  if (!std::isfinite(param))
  {
    return false;
  }

  const bool infinite_edge
    = kd != nullptr && kd->getGraph().isInfinite(2 * delaunay_edge_id);
  if (infinite_edge)
  {
    return param >= -eps;
  }
  return param >= -eps && param <= 1.0 + eps;
}

const char* delaunayEdgeParamExpectedRangeDescription(const KineticDelaunay* kd, size_t delaunay_edge_id)
{
  if (kd != nullptr && kd->getGraph().isInfinite(2 * delaunay_edge_id))
  {
    return "[0,+inf)";
  }
  return "[0,1]";
}

std::string sanitizeCrossingParamErrorContextForFilename(const char* context)
{
  std::string safe = (context != nullptr && context[0] != '\0') ? context : "unknown";
  for (char& c : safe)
  {
    if (!(std::isalnum(static_cast<unsigned char>(c)) || c == '-' || c == '_'))
    {
      c = '_';
    }
  }
  return safe;
}
} // namespace

void kinDS::assignCrossingIntersectionDelaunayParam(const KineticDelaunay* kd,
  KineticDelaunay::CrossingData::VoronoiDelaunayEdgeIntersection& intersection, double param, double t,
  const char* context, std::optional<size_t> near_voronoi_vertex, std::optional<size_t> opposite_voronoi_vertex)
{
  intersection.delaunay_edge_param = param;
  intersection.param_last_updated = t;

  constexpr double param_range_eps = 1e-9;
  if (!isDelaunayEdgeParamInExpectedWarningRange(kd, intersection.delaunay_edge_id, param, param_range_eps))
  {
    KINDS_WARNING((context != nullptr ? context : "assignCrossingIntersectionDelaunayParam")
      << ": delaunay_edge_param outside expected range "
      << delaunayEdgeParamExpectedRangeDescription(kd, intersection.delaunay_edge_id) << " param=" << param << " de="
      << intersection.delaunay_edge_id << " ve=" << intersection.voronoi_edge_id << " t=" << t
      << (near_voronoi_vertex.has_value() ? (" near_vv=" + std::to_string(near_voronoi_vertex.value())) : std::string())
      << (opposite_voronoi_vertex.has_value()
            ? (" opposite_vv=" + std::to_string(opposite_voronoi_vertex.value()))
            : std::string()));

    if (kd != nullptr && kd->shouldDumpErrorFiles())
    {
      // SVG writer needs a non-const KineticDelaunay (branch labels / path helpers); dump is read-mostly.
      KineticDelaunay& kd_mut = const_cast<KineticDelaunay&>(*kd);
      const auto& graph = kd_mut.getGraph();
      const size_t delaunay_edge_id = intersection.delaunay_edge_id;
      const size_t voronoi_edge_id = intersection.voronoi_edge_id;

      VisualDebugHighlight highlight;
      // Delaunay edge @c delaunay_edge_id: primal only (directed half-edges + site endpoints).
      const size_t de_even = 2 * delaunay_edge_id;
      if (de_even + 1 < graph.halfEdgeSlotCount())
      {
        highlight.directed_half_edges.insert(de_even);
        highlight.directed_half_edges.insert(de_even ^ 1);
        const int origin0 = graph.halfEdge(de_even).origin;
        const int origin1 = graph.halfEdge(de_even ^ 1).origin;
        if (origin0 >= 0)
        {
          highlight.delaunay_vertices.insert(static_cast<size_t>(origin0));
        }
        if (origin1 >= 0)
        {
          highlight.delaunay_vertices.insert(static_cast<size_t>(origin1));
        }
      }
      // Voronoi edge @c voronoi_edge_id: dual only (do not also mark directed half-edges / primary for this id).
      highlight.voronoi_edges.insert(voronoi_edge_id);
      highlight.label_crossings_on_delaunay_edges.insert(delaunay_edge_id);
      highlight.label_crossings_on_voronoi_edges.insert(voronoi_edge_id);
      highlight.crossing_intersection_keys.insert(
        (static_cast<uint64_t>(delaunay_edge_id) << 32) | static_cast<uint64_t>(voronoi_edge_id));

      if (near_voronoi_vertex.has_value())
      {
        highlight.voronoi_vertices.insert(near_voronoi_vertex.value());
      }
      if (opposite_voronoi_vertex.has_value())
      {
        highlight.voronoi_vertices.insert(opposite_voronoi_vertex.value());
      }
      // If callers did not pass ends, fall back to the Voronoi edge's two face endpoints.
      if (!near_voronoi_vertex.has_value() || !opposite_voronoi_vertex.has_value())
      {
        const size_t he_even = 2 * voronoi_edge_id;
        if (he_even + 1 < graph.halfEdgeSlotCount())
        {
          const int face0 = graph.halfEdge(he_even).face;
          const int face1 = graph.halfEdge(he_even ^ 1).face;
          if (face0 >= 0)
          {
            highlight.voronoi_vertices.insert(static_cast<size_t>(face0));
          }
          if (face1 >= 0)
          {
            highlight.voronoi_vertices.insert(static_cast<size_t>(face1));
          }
        }
      }

      std::ostringstream descriptor;
      descriptor << sanitizeCrossingParamErrorContextForFilename(context) << "_param_FAIL_de" << delaunay_edge_id
                 << "_ve" << voronoi_edge_id << "_dparam";
      {
        std::ostringstream param_token;
        param_token.setf(std::ios::fixed);
        param_token.precision(kDebugExportTimePrecision);
        param_token << std::showpoint << param;
        std::string param_str = param_token.str();
        for (char& ch : param_str)
        {
          if (ch == '-')
          {
            ch = 'm';
          }
          else if (ch == '+')
          {
            ch = 'p';
          }
          else if (ch == '.')
          {
            ch = '_';
          }
          else if (!std::isalnum(static_cast<unsigned char>(ch)))
          {
            ch = '_';
          }
        }
        descriptor << param_str;
      }
      if (opposite_voronoi_vertex.has_value())
      {
        descriptor << "_oppvv" << opposite_voronoi_vertex.value();
      }

      // Route by the Delaunay edge's two site endpoints. Do not rely on highlight inference: this dump also
      // marks Voronoi faces/edges whose vertices can span multiple runtime branches, which used to fan the
      // SVG out to every active branch folder. With --svg-separate-pending-splits, endpoints on different
      // pending children get a duplicate in each folder; otherwise both collapse to the unsplit parent.
      std::vector<size_t> endpoint_runtime_branch_ids;
      const bool separate_pending_splits = kd_mut.visualDebugSeparatePendingSplits();
      auto add_endpoint_branch = [&](int vertex)
      {
        if (vertex < 0 || kd_mut.isDummyBoundary(static_cast<size_t>(vertex)))
        {
          return;
        }
        const size_t strand_id = static_cast<size_t>(vertex);
        const size_t raw_branch_id = kd_mut.getRuntimeBranchData().branchForStrand(strand_id);
        if (raw_branch_id == KineticDelaunay::RuntimeBranchData::no_branch)
        {
          return;
        }
        const size_t folder_branch_id
          = separate_pending_splits ? raw_branch_id : kd_mut.unsplitRuntimeBranchId(raw_branch_id);
        if (std::find(endpoint_runtime_branch_ids.begin(), endpoint_runtime_branch_ids.end(), folder_branch_id)
          == endpoint_runtime_branch_ids.end())
        {
          endpoint_runtime_branch_ids.push_back(folder_branch_id);
        }
      };
      if (de_even + 1 < graph.halfEdgeSlotCount())
      {
        add_endpoint_branch(graph.halfEdge(de_even).origin);
        add_endpoint_branch(graph.halfEdge(de_even ^ 1).origin);
      }

      const std::optional<size_t> preferred_branch = !endpoint_runtime_branch_ids.empty()
        ? std::optional<size_t>(endpoint_runtime_branch_ids.front())
        : std::nullopt;
      const std::vector<size_t>* explicit_branches
        = endpoint_runtime_branch_ids.empty() ? nullptr : &endpoint_runtime_branch_ids;
      writeSegmentBuilderVisualDebugSvg(true, kd_mut, graph, kd_mut.eventTimeAt(t), "error", descriptor.str(),
        highlight, preferred_branch, /*separation_offset_segments=*/nullptr, /*seam_outlines=*/nullptr,
        explicit_branches);
    }
  }
}

void kinDS::updateCrossingIntersectionParam(KineticDelaunay& kd,
  KineticDelaunay::CrossingData::EdgeIntersectionRef intersection, double t)
{
  assignCrossingIntersectionDelaunayParam(&kd, *intersection,
    kd.delaunayVoronoiEdgeIntersection(intersection->delaunay_edge_id, intersection->voronoi_edge_id, t).first, t,
    "updateCrossingIntersectionParam");
}

void kinDS::ensureCrossingIntersectionParamUpToDate(KineticDelaunay& kd,
  KineticDelaunay::CrossingData::EdgeIntersectionRef intersection, double t)
{
  if (intersection->param_last_updated != t)
  {
    updateCrossingIntersectionParam(kd, intersection, t);
  }
}

glm::dvec3 kinDS::getCrossingCoordsInDelaunaySpace(KineticDelaunay& kd,
  KineticDelaunay::CrossingData::EdgeIntersectionRef intersection, double t)
{
  ensureCrossingIntersectionParamUpToDate(kd, intersection, t);
  return crossingCoordsFromDelaunayEdgeInterpolation(
    kd, *intersection, t, [&kd](size_t strand_id, double query_t) { return kd.getPointInDelaunaySpace(strand_id, query_t); });
}

void KineticDelaunay::refreshDelaunayEdgeIntersectionParams(size_t delaunay_edge_id, double t)
{
  if (delaunay_edge_id >= crossing_data.delaunay_edge_intersections.size())
  {
    return;
  }

  auto& d_list = crossing_data.delaunay_edge_intersections[delaunay_edge_id];
  for (auto ref : d_list)
  {
    updateCrossingIntersectionParam(*this, ref, t);
  }

  for (auto it = d_list.begin(); it != d_list.end(); ++it)
  {
    (*it)->delaunay_ref = it;
  }
}

void KineticDelaunay::refreshAndSortDelaunayEdgeIntersectionParams(size_t delaunay_edge_id, double t)
{
  if (delaunay_edge_id >= crossing_data.delaunay_edge_intersections.size())
  {
    return;
  }

  auto& d_list = crossing_data.delaunay_edge_intersections[delaunay_edge_id];
  for (auto ref : d_list)
  {
    updateCrossingIntersectionParam(*this, ref, t);
  }

  d_list.sort([&](const CrossingData::EdgeIntersectionRef& a, const CrossingData::EdgeIntersectionRef& b)
    { return a->delaunay_edge_param < b->delaunay_edge_param; });
  for (auto it = d_list.begin(); it != d_list.end(); ++it)
  {
    (*it)->delaunay_ref = it;
  }
}

void KineticDelaunay::refreshTriangleDelaunayEdgeIntersectionParams(size_t face_id, double t)
{
  const auto& face_half_edges = graph.face(face_id).half_edges;
  std::unordered_set<size_t> visited_delaunay_edges;
  for (size_t he_id : face_half_edges)
  {
    const size_t d_edge_id = he_id / 2;
    if (visited_delaunay_edges.insert(d_edge_id).second)
    {
      refreshDelaunayEdgeIntersectionParams(d_edge_id, t);
    }
  }
}

std::pair<std::vector<size_t>, std::vector<double>> KineticDelaunay::computeCrossedHalfEdges(
  size_t start_face_id, const glm::dvec2& destination, const glm::dvec2& start_point, double t) const
{
  std::vector<size_t> crossed_half_edge_ids;
  std::vector<double> crossed_half_edge_params;
  auto edge_function = [](const glm::dvec2& a, const glm::dvec2& b, const glm::dvec2& c)
  { return -((b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x)); };

  auto compute_angular_bisector_homogeneous = [&](size_t he_id)
  {
    auto ray = computeAngularBisector(he_id, t);
    return lineToHomegeneous(ray.first, ray.second);
  };

  auto compute_he_graph_edge_function = [&](size_t he_id, const glm::dvec2& query_point)
  {
    size_t origin = graph.halfEdge(he_id).origin;
    size_t dest = graph.destination(he_id);
    // KINDS_DEBUG("Computing edge function for half-edge " << he_id << " with origin " << origin << " and destination "
    // << dest);
    if (origin != static_cast<size_t>(-1) && dest != static_cast<size_t>(-1))
    {
      glm::dvec2 p0 = getPointAt(origin, t);
      glm::dvec2 p1 = getPointAt(dest, t);
      return -((p1.x - p0.x) * (query_point.y - p0.y) - (p1.y - p0.y) * (query_point.x - p0.x));
    }
    else
    {
      glm::dvec3 edge = compute_angular_bisector_homogeneous(he_id);

      if (dest == -1)
      {
        return -glm::dot(edge, glm::dvec3(query_point, 1.0));
      }
      else
      {
        return glm::dot(edge, glm::dvec3(query_point, 1.0));
      }
    }
  };

  int next_crossed_edge_id = -1;
  bool inside_triangle = false;
  size_t next_face_id = start_face_id;
  auto next_vertices = graph.getTriangleVertexIndices(start_face_id);
  auto next_tri_half_edges = graph.face(start_face_id).half_edges;

  // KINDS_DEBUG("Following line from " << glm::to_string(start_point) << " to " << glm::to_string(destination));

  while (!inside_triangle)
  {
    // KINDS_DEBUG("Next face ID: " << next_face_id);
    //  check each edge for an intersection with the line we need to follow.
    inside_triangle = true;

    // first test if we are inside
    for (int edge_index = 0; edge_index < 3; edge_index++)
    {
      size_t he_id = next_tri_half_edges[edge_index];
      double edge_function = compute_he_graph_edge_function(he_id, destination);
      if (edge_function < 0)
      {
        inside_triangle = false;
        break;
      }
    }

    // All edges are on the correct side, so we are inside the triangle and can stop.
    if (inside_triangle)
    {
      break;
    }

    // Now determine direction of walk
    size_t max_s_index = -1;
    double max_s = -1.0;
    double crossed_edge_param;

    for (int edge_index = 0; edge_index < 3; edge_index++)
    {
      if (edge_index == 0 && next_face_id != start_face_id)
      {
        // This is where we came from, so we don't need to check it as it will make us go backwards
        continue;
      }
      size_t he_id = next_tri_half_edges[edge_index];

      if (graph.isInfinite(he_id))
      {
        auto ray = computeAngularBisector(he_id, t);
        /*if(segmentRayIntersection(start_point, destination, ray.first, ray.second)){
          next_crossed_edge_id = he_id;
          inside_triangle = false;
          break;
        }*/
        auto [r, s] = segmentIntersectionParameters(ray.first, ray.first + ray.second, start_point, destination);
        if (r >= 0.0 && s <= 1.0)
        {
          if (s > max_s)
          {
            max_s = s;
            max_s_index = edge_index;
            crossed_edge_param = r;
          }
        }
      }
      else
      {
        glm::dvec2 edge_start = getPointAt(static_cast<size_t>(graph.halfEdge(he_id).origin), t);
        glm::dvec2 edge_end = getPointAt(static_cast<size_t>(graph.destination(he_id)), t);
        auto [r, s] = segmentIntersectionParameters(edge_start, edge_end, start_point, destination);
        if (r >= 0.0 && r <= 1.0 && s <= 1.0)
        {
          if (s > max_s)
          {
            max_s = s;
            max_s_index = edge_index;
            crossed_edge_param = r;
          }
        }
      }
    }

    if(max_s_index == static_cast<size_t>(-1))
    {
      KINDS_ERROR("computeCrossedHalfEdges: no valid edge intersection found for line from "
                  << glm::to_string(start_point) << " to " << glm::to_string(destination) << " in triangle "
                  << next_face_id << ", start face: " << start_face_id << ". Crossed edges: " << crossed_half_edge_ids.size());
      break;
    }
    next_crossed_edge_id = next_tri_half_edges[max_s_index];
    next_face_id = graph.halfEdge(next_crossed_edge_id ^ 1).face;
    // Record the edge we are about to cross
    crossed_half_edge_ids.push_back(next_crossed_edge_id);
    crossed_half_edge_params.push_back(crossed_edge_param);
    next_vertices = graph.adjacentTriangleVertices(next_crossed_edge_id ^ 1);
    next_tri_half_edges = graph.getTriangleHalfEdgeIndices(next_crossed_edge_id ^ 1);
  }

  return std::make_pair(crossed_half_edge_ids, crossed_half_edge_params);
}

const HalfEdgeDelaunayGraph& KineticDelaunay::init(CallbackManager* callback_manager)
{
  callback_manager_ = callback_manager;
  const size_t section_count = getSectionCount();
  if (section_count == 0)
  {
    throw std::runtime_error("KineticDelaunay::init: StrandTree has no sections");
  }
  if (start_section_ >= section_count)
  {
    throw std::runtime_error("KineticDelaunay::init: start_section "
      + std::to_string(start_section_) + " is out of range [0, " + std::to_string(section_count) + ")");
  }
  if (end_section_.has_value() && end_section_.value() < start_section_)
  {
    throw std::runtime_error("KineticDelaunay::init: end_section < start_section");
  }
  if (end_section_.has_value() && end_section_.value() > section_count)
  {
    throw std::runtime_error("KineticDelaunay::init: end_section "
      + std::to_string(end_section_.value()) + " is out of range [0, " + std::to_string(section_count) + "]");
  }

  const size_t bootstrap_section = start_section_;
  const double bootstrap_t = static_cast<double>(bootstrap_section);

  const size_t vertex_count = branch_trajs.getPoints().size();
  // Bootstrap with one HalfEdgeDelaunayGraph per input branch at the start height, then combine into a
  // single graph whose vertex indices remain global strand ids.
  const auto& branches_at_bootstrap = branch_trajs.getStrandBranchesByHeight(bootstrap_section);
  std::vector<HalfEdgeDelaunayGraph> branch_graphs;
  std::vector<std::vector<size_t>> local_to_global;
  branch_graphs.reserve(branches_at_bootstrap.size() + (add_dummy_boundary ? 1 : 0));
  local_to_global.reserve(branches_at_bootstrap.size() + (add_dummy_boundary ? 1 : 0));

  const auto append_branch_graph = [this, bootstrap_t, bootstrap_section, &branch_graphs, &local_to_global](
                                     std::vector<size_t> strand_ids)
  {
    if (strand_ids.size() < 3)
    {
      return;
    }

    std::vector<glm::dvec2> sites;
    sites.reserve(strand_ids.size());
    for (size_t strand_id : strand_ids)
    {
      const size_t reference_branch
        = isDummyBoundary(strand_id) ? 0 : branch_trajs.getBranchIndex(strand_id, bootstrap_section);
      sites.push_back(branch_trajs.evaluateTransformed(strand_id, bootstrap_t, reference_branch));
    }

    HalfEdgeDelaunayGraph part;
    part.init(sites);
    branch_graphs.push_back(std::move(part));
    local_to_global.push_back(std::move(strand_ids));
  };

  for (const auto& branch_strands : branches_at_bootstrap)
  {
    std::vector<size_t> component;
    component.reserve(branch_strands.size());
    for (size_t strand_id : branch_strands)
    {
      if (!isDummyBoundary(strand_id))
      {
        component.push_back(strand_id);
      }
    }
    append_branch_graph(std::move(component));
  }
  if (add_dummy_boundary)
  {
    std::vector<size_t> dummy_component;
    dummy_component.reserve(12);
    for (size_t v = 0; v < vertex_count; ++v)
    {
      if (isDummyBoundary(v))
      {
        dummy_component.push_back(v);
      }
    }
    append_branch_graph(std::move(dummy_component));
  }
  if (branch_graphs.empty())
  {
    throw std::runtime_error("KineticDelaunay::init: no triangulable branch components (>= 3 strands) at start_section "
      + std::to_string(bootstrap_section));
  }

  graph.combine(vertex_count, branch_graphs, local_to_global);
  sections_advanced = bootstrap_section;
  pending_branch_splits_.clear();
  pending_branch_splits_.resetStrandLookup(vertex_count);

  /*section_event_manager_->setCallback(section_callback);
  flip_event_manager_->setCallback(flip_callback);
  radius_event_manager_->setCallback(radius_callback);
  crossing_event_manager_->setCallback(crossing_callback);*/

  // graph.printDebug();

  face_inside.clear();
  quadrilateral_last_updated.clear();
  face_last_updated.clear();

  crossing_data.init(graph.faceSlotCount());

  face_inside.resize(graph.faceSlotCount(), false);
  quadrilateral_last_updated.resize(graph.halfEdgeSlotCount() / 2, EventTime{});
  face_last_updated.resize(graph.faceSlotCount(), EventTime{});

  // Bootstrap component_map from the per-branch triangulation so getPointAt works during face initialization.
  computeComponentData(bootstrap_t);
  // Branches already start as separate graph components; do not treat this as a pending split.
  prev_component_count = component_data.components.size();

  for (size_t face_index : graph.liveFaces())
  {
    initializeFaceState(face_index, bootstrap_t);
  }

  validateVoronoiVertexIteratorInvariants("init", bootstrap_t);

  // Input branches at the bootstrap height become runtime branches directly.
  computeComponentData(bootstrap_t);
  prev_component_count = component_data.components.size();
  updateRuntimeBranchMapFromInputBranches(bootstrap_t);

  // Precompute Voronoi–Delaunay edge intersections at the bootstrap time and store them in crossing_data.
  crossing_data.computeEdgeIntersections(*this, bootstrap_t);

  if (callback_manager)
  {
    callback_manager->init();
  }

  return graph;
}

void KineticDelaunay::registerSectionEventCallback(EventCallback* callback)
{
  section_event_manager_->setCallback(callback);
}

void KineticDelaunay::registerFlipEventCallback(EventCallback* callback) { flip_event_manager_->setCallback(callback); }

void KineticDelaunay::registerRadiusEventCallback(EventCallback* callback)
{
  radius_event_manager_->setCallback(callback);
}

void KineticDelaunay::registerCrossingEventCallback(EventCallback* callback)
{
  crossing_event_manager_->setCallback(callback);
}

void KineticDelaunay::registerSubdivisionEventCallback(EventCallback* callback)
{
  subdivision_event_manager_->setCallback(callback);
}

void KineticDelaunay::registerSeparationEventCallback(EventCallback* callback)
{
  separation_event_manager_->setCallback(callback);
}

void KineticDelaunay::setSubdivisionSchedule(std::vector<std::pair<size_t, double>> schedule)
{
  subdivision_schedule_ = std::move(schedule);
}

void KineticDelaunay::enqueueScheduledSubdivisionEvents()
{
  const double min_t = static_cast<double>(getStartSection());
  const double max_t = static_cast<double>(getEndSection());
  for (const auto& strand_and_t : subdivision_schedule_)
  {
    if (strand_and_t.second < min_t || !(strand_and_t.second < max_t))
    {
      continue;
    }
    kinetic_algorithm_->enqueueEvent(
      std::make_shared<SubdivisionEvent>(this, strand_and_t.second, strand_and_t.first, strand_and_t.second));
  }
  subdivision_schedule_.clear();
}

void KineticDelaunay::registerEventCallbacks(EventCallback* section_callback, EventCallback* flip_callback,
  EventCallback* radius_callback, EventCallback* crossing_callback)
{
  registerSectionEventCallback(section_callback);
  registerFlipEventCallback(flip_callback);
  registerRadiusEventCallback(radius_callback);
  registerCrossingEventCallback(crossing_callback);
}

void KineticDelaunay::CrossingData::computeEdgeIntersections(const KineticDelaunay& kd, double t)
{
  // Clear any previous intersections
  voronoi_edge_intersections.clearAllLists();
  edge_intersections.clear();

  size_t num_edges = kd.getGraph().halfEdgeSlotCount() / 2;
  voronoi_edge_intersections.resize(num_edges);
  delaunay_edge_intersections.resize(num_edges);

  const auto& graph = kd.getGraph();

  // Fill intersections for each Voronoi edge (dual to a Delaunay edge)
  for (size_t voronoi_edge_id = 0; voronoi_edge_id < num_edges; ++voronoi_edge_id)
  {
    size_t he_id0 = voronoi_edge_id * 2;
    size_t he_id1 = he_id0 + 1;

    if (graph.isInfinite(he_id0))
    {
      continue;
    }

    glm::dvec3 left_pos = kd.computeVoronoiVertexClampedInfinity(he_id0, t);
    glm::dvec3 right_pos = kd.computeVoronoiVertexClampedInfinity(he_id1, t);

    glm::dvec2 left_vertex(left_pos.x, left_pos.y);
    glm::dvec2 right_vertex(right_pos.x, right_pos.y);

    size_t left_voronoi_vertex_id = graph.halfEdge(he_id0).face;
    size_t left_containing_tri_id = kd.getCrossingDataContainingTriId(left_voronoi_vertex_id);

    auto crossed_half_edges_params = kd.computeCrossedHalfEdges(left_containing_tri_id, right_vertex, left_vertex, t);

    for (size_t i = 0; i < crossed_half_edges_params.first.size(); i++)
    {
      size_t delaunay_he_id = crossed_half_edges_params.first[i];
      double param = crossed_half_edges_params.second[i];

      edge_intersections.emplace_back();
      auto edge_itr = std::prev(edge_intersections.end());

      edge_itr->delaunay_edge_id = delaunay_he_id / 2;
      edge_itr->voronoi_edge_id = voronoi_edge_id;
      edge_itr->prev_segment_mesh_pair_index = static_cast<size_t>(-1);
      edge_itr->next_segment_mesh_pair_index = static_cast<size_t>(-1);

      if (delaunay_he_id % 2 == 0)
      {
        assignCrossingIntersectionDelaunayParam(
          &kd, *edge_itr, param, t, "initializeVoronoiEdgeIntersections:even");
      }
      else
      {
        assignCrossingIntersectionDelaunayParam(
          &kd, *edge_itr, 1.0 - param, t, "initializeVoronoiEdgeIntersections:odd");
      }

      auto voronoi_ref = voronoi_edge_intersections[voronoi_edge_id].emplace(
        voronoi_edge_intersections[voronoi_edge_id].end(), edge_itr);
      edge_itr->voronoi_ref = voronoi_ref;
    }
  }

  // Populate delaunay_edge_intersections and sort by parameter along the Delaunay edge.
  for (auto edge_itr = edge_intersections.begin(); edge_itr != edge_intersections.end(); ++edge_itr)
  {
    size_t delaunay_edge_id = edge_itr->delaunay_edge_id;
    if (delaunay_edge_id >= delaunay_edge_intersections.size())
    {
      KINDS_ERROR(
        "Delaunay edge id out of bounds: " << delaunay_edge_id << " >= " << delaunay_edge_intersections.size());
      continue;
    }
    delaunay_edge_intersections[delaunay_edge_id].push_back(edge_itr);
  }

  for (size_t delaunay_edge_id = 0; delaunay_edge_id < delaunay_edge_intersections.size(); ++delaunay_edge_id)
  {
    auto& edge_list = delaunay_edge_intersections[delaunay_edge_id];
    edge_list.sort([&](const EdgeIntersectionRef& a, const EdgeIntersectionRef& b)
      { return a->delaunay_edge_param < b->delaunay_edge_param; });

    for (auto it = edge_list.begin(); it != edge_list.end(); ++it)
    {
      (*it)->delaunay_ref = it;
    }
  }

  validateIntersectionInvariants("computeEdgeIntersections", &kd, t);
}

namespace
{
struct IntersectionLogEntry
{
  size_t delaunay_edge_id;
  size_t voronoi_edge_id;
  double delaunay_edge_param;
};

struct InvariantViolationScope
{
  std::optional<size_t> primary_dual_edge;
  std::unordered_set<size_t> auxiliary_dual_edges;
  std::unordered_set<uint64_t> crossing_keys;
  std::vector<IntersectionLogEntry> intersection_log;

  void noteIntersection(size_t delaunay_edge_id, size_t voronoi_edge_id, double delaunay_edge_param)
  {
    crossing_keys.insert((static_cast<uint64_t>(delaunay_edge_id) << 32) | voronoi_edge_id);
    intersection_log.push_back({ delaunay_edge_id, voronoi_edge_id, delaunay_edge_param });
  }

  void setPrimaryDualEdge(size_t edge_id) { primary_dual_edge = edge_id; }

  void noteAuxiliaryDualEdge(size_t edge_id) { auxiliary_dual_edges.insert(edge_id); }
};

std::string formatIntersectionLogEntry(const IntersectionLogEntry& entry)
{
  std::ostringstream oss;
  oss.setf(std::ios::fixed);
  oss.precision(15);
  oss << "(d=" << entry.delaunay_edge_id << ",v=" << entry.voronoi_edge_id << ",param=" << entry.delaunay_edge_param
      << ")";
  return oss.str();
}

double computeDelaunayEdgeParamForInvariant(
  const KineticDelaunay& kd, KineticDelaunay::CrossingData::EdgeIntersectionRef ref, double t)
{
  return kd.delaunayVoronoiEdgeIntersection(ref->delaunay_edge_id, ref->voronoi_edge_id, t).first;
}

std::string formatDelaunayEdgeIntersectionList(
  const KineticDelaunay::CrossingData& crossing_data, const KineticDelaunay& kd, size_t delaunay_edge_id, double t)
{
  std::ostringstream oss;
  oss << "delaunay_edge_intersections[" << delaunay_edge_id << "]={";
  if (delaunay_edge_id < crossing_data.delaunay_edge_intersections.size())
  {
    const auto& d_list = crossing_data.delaunay_edge_intersections[delaunay_edge_id];
    size_t list_index = 0;
    for (auto d_it = d_list.begin(); d_it != d_list.end(); ++d_it, ++list_index)
    {
      if (list_index > 0)
      {
        oss << ", ";
      }
      const auto& ir = **d_it;
      const double d_param = computeDelaunayEdgeParamForInvariant(kd, *d_it, t);
      oss << "[" << list_index << "]"
          << formatIntersectionLogEntry({ ir.delaunay_edge_id, ir.voronoi_edge_id, ir.delaunay_edge_param })
          << ",recomputedDParam=" << d_param;
    }
  }
  oss << "}";
  return oss.str();
}

double computeVoronoiEdgeParamForInvariant(
  const KineticDelaunay& kd, KineticDelaunay::CrossingData::EdgeIntersectionRef ref, double t)
{
  return kd.delaunayVoronoiEdgeIntersection(ref->delaunay_edge_id, ref->voronoi_edge_id, t).second;
}

std::string formatVoronoiEdgeIntersectionList(
  const KineticDelaunay::CrossingData& crossing_data, const KineticDelaunay& kd, size_t voronoi_edge_id, double t)
{
  std::ostringstream oss;
  oss << "voronoi_edge_intersections[" << voronoi_edge_id << "]={";
  if (voronoi_edge_id < crossing_data.voronoi_edge_intersections.size())
  {
    const auto& v_list = crossing_data.voronoi_edge_intersections[voronoi_edge_id];
    size_t list_index = 0;
    for (auto v_it = v_list.begin(); v_it != v_list.end(); ++v_it, ++list_index)
    {
      if (list_index > 0)
      {
        oss << ", ";
      }
      const auto& ir = **v_it;
      const double v_param = computeVoronoiEdgeParamForInvariant(kd, *v_it, t);
      oss << "[" << list_index << "]"
          << formatIntersectionLogEntry({ ir.delaunay_edge_id, ir.voronoi_edge_id, ir.delaunay_edge_param })
          << ",vParam=" << v_param;
    }
  }
  oss << "}";
  return oss.str();
}

void logInvariantViolationScope(const InvariantViolationScope& scope, double t)
{
  if (scope.primary_dual_edge.has_value())
  {
    KINDS_ERROR("  t=" << t << " primary dual edge (magenta): " << *scope.primary_dual_edge);
  }
  if (!scope.auxiliary_dual_edges.empty())
  {
    KINDS_ERROR("  t=" << t << " auxiliary dual edges:");
    for (size_t edge_id : scope.auxiliary_dual_edges)
    {
      KINDS_ERROR("    t=" << t << " " << edge_id);
    }
  }
  if (!scope.intersection_log.empty())
  {
    KINDS_ERROR("  t=" << t << " highlighted intersections:");
    for (const IntersectionLogEntry& entry : scope.intersection_log)
    {
      KINDS_ERROR("    t=" << t << " " << formatIntersectionLogEntry(entry));
    }
  }
}

std::string sanitizeContextForFilename(const char* context)
{
  std::string safe = (context != nullptr && context[0] != '\0') ? context : "unknown";
  for (char& c : safe)
  {
    if (!(std::isalnum(static_cast<unsigned char>(c)) || c == '-' || c == '_'))
    {
      c = '_';
    }
  }
  return safe;
}

void exportCrossingInvariantFailureDebugSvg(
  const KineticDelaunay& kd, const char* context, double t, const InvariantViolationScope& scope)
{
  (void)kd;
  (void)context;
  (void)t;
  (void)scope;
  return;

  const VisualDebugHighlight highlight = VisualDebugHighlight::forInvariantViolation(
    kd.getGraph(), scope.primary_dual_edge, scope.auxiliary_dual_edges, scope.crossing_keys);

  const std::vector<glm::dvec2> points = kd.getPointsAt(t);
  const std::string filename
    = formatDebugExportTimeToken(kd.eventTimeAt(t)) + "_crossing_invariant_FAIL_"
    + sanitizeContextForFilename(context) + ".svg";
  const auto& containing_tri_ids = kd.getCrossingData().getContainingTriIds();
  const auto intersection_debug_data = kd.getCrossingIntersectionDebugData();

  KINDS_ERROR("Writing crossing invariant failure debug SVG at t=" << kd.eventTimeAt(t) << ": " << filename);
  HalfEdgeDelaunayGraphToSVG::write(points, kd.getGraph(), filename, 0.1, &kd.getFacesInside(), true,
    &containing_tri_ids, &intersection_debug_data, &highlight);
}

[[noreturn]] void failCrossingIntersectionInvariant(const char* context, const std::string& detail,
  const InvariantViolationScope& scope, const KineticDelaunay* kd, double t)
{
  if (kd != nullptr)
  {
    exportCrossingInvariantFailureDebugSvg(*kd, context, t, scope);
  }
  KINDS_ERROR("validateIntersectionInvariants(" << context << ", t=" << t << "): " << detail);
  logInvariantViolationScope(scope, t);
  throw std::runtime_error(
    std::string("validateIntersectionInvariants(") + context + ", t=" + std::to_string(t) + "): " + detail);
}

bool listContainsIntersectionRef(const std::list<KineticDelaunay::CrossingData::EdgeIntersectionRef>& list,
  std::list<KineticDelaunay::CrossingData::VoronoiDelaunayEdgeIntersection>::const_iterator target)
{
  for (auto it = list.begin(); it != list.end(); ++it)
  {
    if (*it == target)
    {
      return true;
    }
  }
  return false;
}

bool eraseIntersectionRefFromList(std::list<KineticDelaunay::CrossingData::EdgeIntersectionRef>& list,
  std::list<KineticDelaunay::CrossingData::VoronoiDelaunayEdgeIntersection>::const_iterator target)
{
  for (auto it = list.begin(); it != list.end(); ++it)
  {
    if (*it == target)
    {
      list.erase(it);
      return true;
    }
  }
  return false;
}

void noteAllIntersections(const KineticDelaunay::CrossingData& crossing_data, InvariantViolationScope& scope)
{
  for (auto ref = crossing_data.edge_intersections.begin(); ref != crossing_data.edge_intersections.end(); ++ref)
  {
    scope.noteIntersection(ref->delaunay_edge_id, ref->voronoi_edge_id, ref->delaunay_edge_param);
  }
}
} // namespace

void KineticDelaunay::CrossingData::validateVoronoiVertexIteratorInvariants(
  const char* context, const HalfEdgeDelaunayGraph& graph, const KineticDelaunay* kd, double t) const
{
  (void)kd;
  (void)t;
  const char* ctx = (context != nullptr && context[0] != '\0') ? context : "unknown";

  for (size_t voronoi_vertex_id = 0; voronoi_vertex_id < voronoi_vertex_to_containing_tri_id.size(); ++voronoi_vertex_id)
  {
    const size_t containing_tri_id = voronoi_vertex_to_containing_tri_id[voronoi_vertex_id];
    const std::optional<std::list<size_t>::iterator>& cached_it = voronoi_vertex_to_iterator[voronoi_vertex_id];
    const bool registered = containing_tri_id != invalid_containing_tri_id;

    if (!registered)
    {
      if (cached_it.has_value())
      {
        throw std::runtime_error(std::string(ctx) + ": Voronoi vertex " + std::to_string(voronoi_vertex_id)
          + " has no containing triangle but a cached list iterator.");
      }
      continue;
    }

    if (voronoi_vertex_id < graph.faceSlotCount() && !graph.isLiveFace(voronoi_vertex_id))
    {
      throw std::runtime_error(std::string(ctx) + ": dead dual Voronoi vertex " + std::to_string(voronoi_vertex_id)
        + " is still registered in triangle " + std::to_string(containing_tri_id) + ".");
    }

    if (containing_tri_id >= tri_id_to_voronoi_vertices.size() || !tri_id_to_voronoi_vertices[containing_tri_id])
    {
      throw std::runtime_error(std::string(ctx) + ": Voronoi vertex " + std::to_string(voronoi_vertex_id)
        + " has out-of-range containing triangle " + std::to_string(containing_tri_id) + ".");
    }

    if (!cached_it.has_value())
    {
      throw std::runtime_error(std::string(ctx) + ": Voronoi vertex " + std::to_string(voronoi_vertex_id)
        + " is registered in triangle " + std::to_string(containing_tri_id) + " without a cached list iterator.");
    }

    if (**cached_it != voronoi_vertex_id)
    {
      throw std::runtime_error(std::string(ctx) + ": Voronoi vertex " + std::to_string(voronoi_vertex_id)
        + " cached iterator refers to id " + std::to_string(**cached_it) + ".");
    }

    bool iterator_in_containing_list = false;
    for (auto list_it = triList(containing_tri_id).begin(); list_it != triList(containing_tri_id).end(); ++list_it)
    {
      if (list_it == *cached_it)
      {
        iterator_in_containing_list = true;
        break;
      }
    }
    if (!iterator_in_containing_list)
    {
      throw std::runtime_error(std::string(ctx) + ": Voronoi vertex " + std::to_string(voronoi_vertex_id)
        + " cached iterator is not an element of tri_id_to_voronoi_vertices[" + std::to_string(containing_tri_id)
        + "].");
    }
  }

  for (size_t tri_id = 0; tri_id < tri_id_to_voronoi_vertices.size(); ++tri_id)
  {
    if (!tri_id_to_voronoi_vertices[tri_id])
    {
      continue;
    }

    for (auto list_it = triList(tri_id).begin(); list_it != triList(tri_id).end(); ++list_it)
    {
      const size_t voronoi_vertex_id = *list_it;
      if (voronoi_vertex_id >= voronoi_vertex_to_containing_tri_id.size())
      {
        throw std::runtime_error(std::string(ctx) + ": triangle " + std::to_string(tri_id) + " lists unknown Voronoi vertex "
          + std::to_string(voronoi_vertex_id) + ".");
      }
      if (voronoi_vertex_to_containing_tri_id[voronoi_vertex_id] != tri_id)
      {
        throw std::runtime_error(std::string(ctx) + ": triangle " + std::to_string(tri_id) + " lists Voronoi vertex "
          + std::to_string(voronoi_vertex_id) + " but its containing tri is "
          + std::to_string(voronoi_vertex_to_containing_tri_id[voronoi_vertex_id]) + ".");
      }
      if (!voronoi_vertex_to_iterator[voronoi_vertex_id].has_value()
        || *voronoi_vertex_to_iterator[voronoi_vertex_id] != list_it)
      {
        throw std::runtime_error(std::string(ctx) + ": triangle " + std::to_string(tri_id)
          + " list iterator for Voronoi vertex " + std::to_string(voronoi_vertex_id)
          + " does not match voronoi_vertex_to_iterator.");
      }
    }
  }
}

void KineticDelaunay::CrossingData::validateIntersectionInvariants(
  const char* context, const KineticDelaunay* kd, double t) const
{
  (void)context;
  (void)kd;
  (void)t;
  return;

  const char* ctx = (context != nullptr && context[0] != '\0') ? context : "unknown";

  size_t delaunay_list_entries = 0;
  for (size_t d_id = 0; d_id < delaunay_edge_intersections.size(); ++d_id)
  {
    delaunay_list_entries += delaunay_edge_intersections[d_id].size();
  }

  size_t voronoi_list_entries = 0;
  for (size_t v_id = 0; v_id < voronoi_edge_intersections.size(); ++v_id)
  {
    voronoi_list_entries += voronoi_edge_intersections[v_id].size();
  }

  const size_t global_count = edge_intersections.size();
  if (delaunay_list_entries != global_count)
  {
    InvariantViolationScope scope;
    noteAllIntersections(*this, scope);
    for (size_t d_id = 0; d_id < delaunay_edge_intersections.size(); ++d_id)
    {
      if (!delaunay_edge_intersections[d_id].empty())
      {
        scope.noteAuxiliaryDualEdge(d_id);
      }
    }
    failCrossingIntersectionInvariant(ctx,
      "delaunay_edge_intersections has " + std::to_string(delaunay_list_entries)
        + " entries but edge_intersections has " + std::to_string(global_count),
      scope, kd, t);
  }
  if (voronoi_list_entries != global_count)
  {
    InvariantViolationScope scope;
    noteAllIntersections(*this, scope);
    for (size_t v_id = 0; v_id < voronoi_edge_intersections.size(); ++v_id)
    {
      if (!voronoi_edge_intersections[v_id].empty())
      {
        scope.noteAuxiliaryDualEdge(v_id);
      }
    }
    failCrossingIntersectionInvariant(ctx,
      "voronoi_edge_intersections has " + std::to_string(voronoi_list_entries) + " entries but edge_intersections has "
        + std::to_string(global_count),
      scope, kd, t);
  }

  for (auto ref = edge_intersections.begin(); ref != edge_intersections.end(); ++ref)
  {
    const size_t d_id = ref->delaunay_edge_id;
    const size_t v_id = ref->voronoi_edge_id;

    if (d_id >= delaunay_edge_intersections.size())
    {
      InvariantViolationScope scope;
      scope.setPrimaryDualEdge(d_id);
      scope.noteIntersection(d_id, v_id, ref->delaunay_edge_param);
      failCrossingIntersectionInvariant(ctx,
        "intersection has delaunay_edge_id=" + std::to_string(d_id)
          + " >= delaunay_edge_intersections.size()=" + std::to_string(delaunay_edge_intersections.size()),
        scope, kd, t);
    }
    if (v_id >= voronoi_edge_intersections.size())
    {
      InvariantViolationScope scope;
      scope.setPrimaryDualEdge(v_id);
      scope.noteIntersection(d_id, v_id, ref->delaunay_edge_param);
      failCrossingIntersectionInvariant(ctx,
        "intersection has voronoi_edge_id=" + std::to_string(v_id)
          + " >= voronoi_edge_intersections.size()=" + std::to_string(voronoi_edge_intersections.size()),
        scope, kd, t);
    }

    if (!ref->delaunay_ref.has_value())
    {
      InvariantViolationScope scope;
      scope.setPrimaryDualEdge(d_id);
      scope.noteIntersection(d_id, v_id, ref->delaunay_edge_param);
      failCrossingIntersectionInvariant(
        ctx, "intersection has unset delaunay_ref (d=" + std::to_string(d_id) + ", v=" + std::to_string(v_id) + ")",
        scope, kd, t);
    }
    if (!ref->voronoi_ref.has_value())
    {
      InvariantViolationScope scope;
      scope.setPrimaryDualEdge(d_id);
      scope.noteIntersection(d_id, v_id, ref->delaunay_edge_param);
      failCrossingIntersectionInvariant(
        ctx, "intersection has unset voronoi_ref (d=" + std::to_string(d_id) + ", v=" + std::to_string(v_id) + ")",
        scope, kd, t);
    }

    const auto& d_list = delaunay_edge_intersections[d_id];
    const auto& v_list = voronoi_edge_intersections[v_id];

    if (!listContainsIntersectionRef(d_list, ref))
    {
      InvariantViolationScope scope;
      scope.setPrimaryDualEdge(d_id);
      scope.noteIntersection(d_id, v_id, ref->delaunay_edge_param);
      failCrossingIntersectionInvariant(ctx,
        "intersection is not present in delaunay_edge_intersections[" + std::to_string(d_id)
          + "] (v=" + std::to_string(v_id) + ", param=" + std::to_string(ref->delaunay_edge_param) + ")",
        scope, kd, t);
    }
    if (!listContainsIntersectionRef(v_list, ref))
    {
      InvariantViolationScope scope;
      scope.setPrimaryDualEdge(d_id);
      scope.noteIntersection(d_id, v_id, ref->delaunay_edge_param);
      failCrossingIntersectionInvariant(ctx,
        "intersection is not present in voronoi_edge_intersections[" + std::to_string(v_id)
          + "] (d=" + std::to_string(d_id) + ", param=" + std::to_string(ref->delaunay_edge_param) + ")",
        scope, kd, t);
    }
  }

  for (size_t d_id = 0; d_id < delaunay_edge_intersections.size(); ++d_id)
  {
    const auto& d_list = delaunay_edge_intersections[d_id];
    double prev_recomputed_d_param = -std::numeric_limits<double>::infinity();
    std::optional<EdgeIntersectionRef> prev_ref;
    double prev_ref_recomputed_d_param = -std::numeric_limits<double>::infinity();
    size_t list_index = 0;
    for (auto d_it = d_list.begin(); d_it != d_list.end(); ++d_it, ++list_index)
    {
      const EdgeIntersectionRef ref = *d_it;
      if (ref->delaunay_edge_id != d_id)
      {
        InvariantViolationScope scope;
        scope.setPrimaryDualEdge(d_id);
        scope.noteIntersection(ref->delaunay_edge_id, ref->voronoi_edge_id, ref->delaunay_edge_param);
        if (ref->delaunay_edge_id != d_id)
        {
          scope.noteAuxiliaryDualEdge(ref->delaunay_edge_id);
        }
        failCrossingIntersectionInvariant(ctx,
          "delaunay_edge_intersections[" + std::to_string(d_id) + "] entry at list index " + std::to_string(list_index)
            + " has delaunay_edge_id=" + std::to_string(ref->delaunay_edge_id) + " but expected " + std::to_string(d_id)
            + "; entry="
            + formatIntersectionLogEntry({ ref->delaunay_edge_id, ref->voronoi_edge_id, ref->delaunay_edge_param }),
          scope, kd, t);
      }
      if (kd != nullptr)
      {
        const double recomputed_d_param = computeDelaunayEdgeParamForInvariant(*kd, ref, t);
        if (!std::isfinite(recomputed_d_param))
        {
          InvariantViolationScope scope;
          scope.setPrimaryDualEdge(d_id);
          scope.noteIntersection(ref->delaunay_edge_id, ref->voronoi_edge_id, ref->delaunay_edge_param);
          scope.noteAuxiliaryDualEdge(ref->voronoi_edge_id);
          failCrossingIntersectionInvariant(ctx,
            "could not compute finite Delaunay-edge parameter for delaunay_edge_intersections[" + std::to_string(d_id)
              + "] list index " + std::to_string(list_index) + "; entry="
              + formatIntersectionLogEntry({ ref->delaunay_edge_id, ref->voronoi_edge_id, ref->delaunay_edge_param })
              + ", recomputedDParam=" + std::to_string(recomputed_d_param),
            scope, kd, t);
        }

        constexpr double param_order_eps = 1e-12;
        constexpr double param_range_eps = 1e-9;
        if (!isDelaunayEdgeParamInExpectedWarningRange(kd, ref->delaunay_edge_id, recomputed_d_param, param_range_eps))
        {
          KINDS_WARNING("validateIntersectionInvariants(" << ctx << ", t=" << t
            << "): recomputed Delaunay-edge parameter is outside expected range "
            << delaunayEdgeParamExpectedRangeDescription(kd, ref->delaunay_edge_id)
            << " for delaunay_edge_intersections[" << d_id << "] list index " << list_index << "; entry="
            << formatIntersectionLogEntry({ ref->delaunay_edge_id, ref->voronoi_edge_id, ref->delaunay_edge_param })
            << ", recomputedDParam=" << recomputed_d_param << "; "
            << formatDelaunayEdgeIntersectionList(*this, *kd, d_id, t));
        }
        if (recomputed_d_param + param_order_eps < prev_recomputed_d_param)
        {
          InvariantViolationScope scope;
          scope.setPrimaryDualEdge(d_id);
          if (prev_ref.has_value())
          {
            const EdgeIntersectionRef& prev = *prev_ref;
            scope.noteIntersection(prev->delaunay_edge_id, prev->voronoi_edge_id, prev->delaunay_edge_param);
            scope.noteAuxiliaryDualEdge(prev->voronoi_edge_id);
          }
          scope.noteIntersection(ref->delaunay_edge_id, ref->voronoi_edge_id, ref->delaunay_edge_param);
          scope.noteAuxiliaryDualEdge(ref->voronoi_edge_id);

          std::string detail = "delaunay_edge_intersections[" + std::to_string(d_id)
            + "] not sorted by recomputed Delaunay-edge parameter at list index " + std::to_string(list_index);
          detail += "; " + formatDelaunayEdgeIntersectionList(*this, *kd, d_id, t);
          if (prev_ref.has_value())
          {
            const EdgeIntersectionRef& prev = *prev_ref;
            detail += "; out-of-order pair: "
              + formatIntersectionLogEntry({ prev->delaunay_edge_id, prev->voronoi_edge_id, prev->delaunay_edge_param })
              + ",recomputedDParam=" + std::to_string(prev_ref_recomputed_d_param) + " then "
              + formatIntersectionLogEntry({ ref->delaunay_edge_id, ref->voronoi_edge_id, ref->delaunay_edge_param })
              + ",recomputedDParam=" + std::to_string(recomputed_d_param);
            detail += " (refs " + formatCrossingIntersectionForLog(*kd, prev) + " then "
              + formatCrossingIntersectionForLog(*kd, ref) + ")";
          }
          failCrossingIntersectionInvariant(ctx, detail, scope, kd, t);
        }
        prev_recomputed_d_param = recomputed_d_param;
        prev_ref_recomputed_d_param = recomputed_d_param;
        prev_ref = ref;
      }
    }
  }

  for (size_t v_id = 0; v_id < voronoi_edge_intersections.size(); ++v_id)
  {
    const auto& v_list = voronoi_edge_intersections[v_id];
    double prev_v_param = -std::numeric_limits<double>::infinity();
    std::optional<EdgeIntersectionRef> prev_ref;
    double prev_ref_v_param = -std::numeric_limits<double>::infinity();
    size_t list_index = 0;
    for (auto v_it = v_list.begin(); v_it != v_list.end(); ++v_it, ++list_index)
    {
      const EdgeIntersectionRef ref = *v_it;
      if (ref->voronoi_edge_id != v_id)
      {
        InvariantViolationScope scope;
        scope.setPrimaryDualEdge(v_id);
        scope.noteIntersection(ref->delaunay_edge_id, ref->voronoi_edge_id, ref->delaunay_edge_param);
        if (ref->voronoi_edge_id != v_id)
        {
          scope.noteAuxiliaryDualEdge(ref->voronoi_edge_id);
        }
        failCrossingIntersectionInvariant(ctx,
          "voronoi_edge_intersections[" + std::to_string(v_id) + "] entry has voronoi_edge_id="
            + std::to_string(ref->voronoi_edge_id) + " but expected " + std::to_string(v_id) + "; entry="
            + formatIntersectionLogEntry({ ref->delaunay_edge_id, ref->voronoi_edge_id, ref->delaunay_edge_param }),
          scope, kd, t);
      }
      if (kd != nullptr)
      {
        const double v_param = computeVoronoiEdgeParamForInvariant(*kd, ref, t);
        if (!std::isfinite(v_param))
        {
          InvariantViolationScope scope;
          scope.setPrimaryDualEdge(v_id);
          scope.noteIntersection(ref->delaunay_edge_id, ref->voronoi_edge_id, ref->delaunay_edge_param);
          scope.noteAuxiliaryDualEdge(ref->delaunay_edge_id);
          failCrossingIntersectionInvariant(ctx,
            "could not compute finite Voronoi-edge parameter for voronoi_edge_intersections[" + std::to_string(v_id)
              + "] list index " + std::to_string(list_index) + "; entry="
              + formatIntersectionLogEntry({ ref->delaunay_edge_id, ref->voronoi_edge_id, ref->delaunay_edge_param })
              + ", computed vParam=" + std::to_string(v_param),
            scope, kd, t);
        }

        constexpr double param_order_eps = 1e-12;
        constexpr double param_range_eps = 1e-9;
        if (v_param < -param_range_eps || v_param > 1.0 + param_range_eps)
        {
          KINDS_WARNING("validateIntersectionInvariants(" << ctx << ", t=" << t
            << "): recomputed Voronoi-edge parameter is outside [0,1] for voronoi_edge_intersections[" << v_id
            << "] list index " << list_index << "; entry="
            << formatIntersectionLogEntry({ ref->delaunay_edge_id, ref->voronoi_edge_id, ref->delaunay_edge_param })
            << ", vParam=" << v_param << "; " << formatVoronoiEdgeIntersectionList(*this, *kd, v_id, t));
        }
        if (v_param + param_order_eps < prev_v_param)
        {
          InvariantViolationScope scope;
          scope.setPrimaryDualEdge(v_id);
          if (prev_ref.has_value())
          {
            const EdgeIntersectionRef& prev = *prev_ref;
            scope.noteIntersection(prev->delaunay_edge_id, prev->voronoi_edge_id, prev->delaunay_edge_param);
            scope.noteAuxiliaryDualEdge(prev->delaunay_edge_id);
          }
          scope.noteIntersection(ref->delaunay_edge_id, ref->voronoi_edge_id, ref->delaunay_edge_param);
          scope.noteAuxiliaryDualEdge(ref->delaunay_edge_id);

          std::string detail = "voronoi_edge_intersections[" + std::to_string(v_id)
            + "] not sorted by recomputed Voronoi-edge parameter at list index " + std::to_string(list_index);
          detail += "; " + formatVoronoiEdgeIntersectionList(*this, *kd, v_id, t);
          if (prev_ref.has_value())
          {
            const EdgeIntersectionRef& prev = *prev_ref;
            detail += "; out-of-order pair: "
              + formatIntersectionLogEntry({ prev->delaunay_edge_id, prev->voronoi_edge_id, prev->delaunay_edge_param })
              + ",vParam=" + std::to_string(prev_ref_v_param) + " then "
              + formatIntersectionLogEntry({ ref->delaunay_edge_id, ref->voronoi_edge_id, ref->delaunay_edge_param })
              + ",vParam=" + std::to_string(v_param);
            detail += " (refs " + formatCrossingIntersectionForLog(*kd, prev) + " then "
              + formatCrossingIntersectionForLog(*kd, ref) + ")";
          }
          failCrossingIntersectionInvariant(ctx, detail, scope, kd, t);
        }
        prev_v_param = v_param;
        prev_ref_v_param = v_param;
        prev_ref = ref;
      }
    }
  }
}

void KineticDelaunay::CrossingData::removeIntersection(EdgeIntersectionRef intersection_ref)
{
  // Remove from Delaunay list if the cached iterator is still initialized and points into the expected list.
  auto& d_list = delaunay_edge_intersections[intersection_ref->delaunay_edge_id];
  eraseIntersectionRefFromList(d_list, intersection_ref);
  intersection_ref->delaunay_ref.reset();

  // Remove from Voronoi list if the cached iterator is still initialized and points into the expected list.
  auto& v_list = voronoi_edge_intersections[intersection_ref->voronoi_edge_id];
  eraseIntersectionRefFromList(v_list, intersection_ref);
  intersection_ref->voronoi_ref.reset();

  // Finally remove from the global list.
  edge_intersections.erase(intersection_ref);
}

void KineticDelaunay::CrossingData::removeIntersectionsOnDeadDelaunayEdges(const HalfEdgeDelaunayGraph& graph)
{
  for (size_t delaunay_edge_id = 0; delaunay_edge_id < delaunay_edge_intersections.size(); ++delaunay_edge_id)
  {
    const size_t he_id = 2 * delaunay_edge_id;
    if (he_id < graph.halfEdgeSlotCount() && graph.isLiveHalfEdge(he_id))
    {
      continue;
    }

    auto doomed = std::vector<EdgeIntersectionRef>(
      delaunay_edge_intersections[delaunay_edge_id].begin(), delaunay_edge_intersections[delaunay_edge_id].end());
    for (const EdgeIntersectionRef ref : doomed)
    {
      removeIntersection(ref);
    }
  }
}

void KineticDelaunay::CrossingData::updateAfterCrossingEvent(
  const KineticDelaunay& kd, const KineticDelaunay::CrossingEvent& e)
{
  auto& graph = kd.getGraph();
  size_t voronoi_vertex_id = e.voronoi_vertex_id;
  size_t crossed_delaunay_edge_id = e.half_edge_id / 2;
  auto& d_intersections = delaunay_edge_intersections[crossed_delaunay_edge_id];

  glm::dvec3 voronoi_vertex_position = glm::dvec3(e.position, e.occurrence_time);
  auto half_edges = graph.face(voronoi_vertex_id).half_edges;

  bool erased[3] = { false, false, false };
  std::list<EdgeIntersectionRef>::iterator next_after_deletion = d_intersections.end();

  // First remove any intersection entries that involved this Voronoi vertex and the crossed Delaunay edge.
  for (size_t i = 0; i < half_edges.size(); i++)
  {
    size_t voronoi_he_id = half_edges[i];
    size_t voronoi_edge_id = voronoi_he_id / 2;
    auto& v_intersections
      = voronoi_edge_intersections[voronoi_edge_id]; // wrong as of now, first edge not deleted in preceding flip event

    auto is_matching = [&](EdgeIntersectionRef ref)
    { return ref->delaunay_edge_id == crossed_delaunay_edge_id && ref->voronoi_edge_id == voronoi_edge_id; };

    if (!v_intersections.empty() && is_matching(v_intersections.front()))
    {
      auto main_ref = v_intersections.front();
      if (!main_ref->delaunay_ref.has_value() || !main_ref->voronoi_ref.has_value())
      {
        throw std::runtime_error("updateAfterCrossingEvent: matched front intersection has unset cached refs");
      }
      next_after_deletion = std::next(*main_ref->delaunay_ref);
      d_intersections.erase(*main_ref->delaunay_ref);
      v_intersections.erase(*main_ref->voronoi_ref);
      edge_intersections.erase(main_ref);
      erased[i] = true;
    }
    else if (!v_intersections.empty() && is_matching(v_intersections.back()))
    {
      auto main_ref = v_intersections.back();
      if (!main_ref->delaunay_ref.has_value() || !main_ref->voronoi_ref.has_value())
      {
        throw std::runtime_error("updateAfterCrossingEvent: matched back intersection has unset cached refs");
      }
      next_after_deletion = std::next(*main_ref->delaunay_ref);
      d_intersections.erase(*main_ref->delaunay_ref);
      v_intersections.erase(*main_ref->voronoi_ref);
      edge_intersections.erase(main_ref);
      erased[i] = true;
    }
  }

  // Now insert new intersection entries for the updated configuration.
  // When two crossings are inserted, order them by orthogonal projection of each Voronoi edge's *other*
  // endpoint onto the crossed Delaunay edge (same axis as dParam). Skip that for a single insert —
  // ordering is unnecessary and the projection can be ill-defined if the Voronoi edge is orthogonal.
  struct PendingCrossingInsert
  {
    size_t voronoi_edge_id = 0;
    double other_endpoint_dparam = 0.0;
    VoronoiDelaunayEdgeIntersection intersection {};
  };

  const double t = e.occurrence_time;
  std::vector<PendingCrossingInsert> pending;
  pending.reserve(2);

  for (size_t i = 0; i < half_edges.size(); i++)
  {
    if (erased[i])
    {
      continue;
    }

    const size_t voronoi_edge_id = half_edges[i] / 2;

    PendingCrossingInsert pending_insert;
    pending_insert.voronoi_edge_id = voronoi_edge_id;
    pending_insert.intersection.delaunay_edge_id = crossed_delaunay_edge_id;
    pending_insert.intersection.voronoi_edge_id = voronoi_edge_id;
    assignCrossingIntersectionDelaunayParam(&kd, pending_insert.intersection,
      kd.delaunayVoronoiEdgeIntersection(crossed_delaunay_edge_id, voronoi_edge_id, t).first, t,
      "updateAfterCrossingEvent:pendingInsert");
    pending.push_back(pending_insert);
  }

  if (pending.size() == 2)
  {
    const glm::dvec2 crossing_xy(e.position);
    const DelaunayEdgeAxis2D d_axis = delaunayEdgeAxisForDParam(kd, graph, crossed_delaunay_edge_id, t);
    const double crossing_dparam = projectPointOntoDelaunayEdgeAxis(d_axis, crossing_xy);

    for (PendingCrossingInsert& pending_insert : pending)
    {
      const glm::dvec2 other_xy = voronoiEdgeOtherEndpointForDelaunayOrder(
        kd, graph, pending_insert.voronoi_edge_id, voronoi_vertex_id, crossing_xy, t);
      pending_insert.other_endpoint_dparam = projectPointOntoDelaunayEdgeAxis(d_axis, other_xy);

      // Relative to the crossing: using the crossing endpoint by mistake yields ~0.
      const double relative_dparam = pending_insert.other_endpoint_dparam - crossing_dparam;
      if (!std::isfinite(relative_dparam) || std::abs(relative_dparam) < 1e-12)
      {
        std::ostringstream oss;
        oss << "updateAfterCrossingEvent: other Voronoi-edge endpoint projects to the crossing on de="
            << crossed_delaunay_edge_id << " ve=" << pending_insert.voronoi_edge_id << " t=" << t
            << " (crossing_dparam=" << crossing_dparam
            << ", other_dparam=" << pending_insert.other_endpoint_dparam
            << ", relative=" << relative_dparam
            << "). Likely used the crossing endpoint instead of the opposite Voronoi vertex.";
        throw std::runtime_error(oss.str());
      }
    }

    std::sort(pending.begin(), pending.end(),
      [](const PendingCrossingInsert& a, const PendingCrossingInsert& b)
      { return a.other_endpoint_dparam < b.other_endpoint_dparam; });

    // Both crossings coincide on the Delaunay edge at the event; force identical dParams so a later
    // param refresh + sort cannot reorder them due to tiny numeric differences.
    const double averaged_dparam
      = 0.5 * (pending[0].intersection.delaunay_edge_param + pending[1].intersection.delaunay_edge_param);
    assignCrossingIntersectionDelaunayParam(
      &kd, pending[0].intersection, averaged_dparam, t, "updateAfterCrossingEvent:averagedDParam");
    assignCrossingIntersectionDelaunayParam(
      &kd, pending[1].intersection, averaged_dparam, t, "updateAfterCrossingEvent:averagedDParam");

    KINDS_DEBUG("updateAfterCrossingEvent split order on de="
      << crossed_delaunay_edge_id << " vv=" << voronoi_vertex_id << " t=" << t
      << " crossing_dparam=" << crossing_dparam << " ve0=" << pending[0].voronoi_edge_id
      << " other_dparam0=" << pending[0].other_endpoint_dparam
      << " relative0=" << (pending[0].other_endpoint_dparam - crossing_dparam)
      << " ve1=" << pending[1].voronoi_edge_id << " other_dparam1=" << pending[1].other_endpoint_dparam
      << " relative1=" << (pending[1].other_endpoint_dparam - crossing_dparam)
      << " averaged_dparam=" << averaged_dparam);
  }

  auto insert_pos = next_after_deletion;
  for (const PendingCrossingInsert& pending_insert : pending)
  {
    const size_t voronoi_edge_id = pending_insert.voronoi_edge_id;
    auto& v_intersections = voronoi_edge_intersections[voronoi_edge_id];

    edge_intersections.push_back(pending_insert.intersection);
    auto main_iter = std::prev(edge_intersections.end());

    auto inserted_it = d_intersections.insert(insert_pos, main_iter);
    main_iter->delaunay_ref = inserted_it;
    insert_pos = std::next(inserted_it);

    if (e.voronoi_vertex_id == graph.halfEdge(2 * voronoi_edge_id).face)
    {
      v_intersections.emplace_front(main_iter);
      main_iter->voronoi_ref = v_intersections.begin();
    }
    else
    {
      v_intersections.emplace_back(main_iter);
      main_iter->voronoi_ref = std::prev(v_intersections.end());
    }
  }
}

const HalfEdgeDelaunayGraph& KineticDelaunay::getGraph() const { return graph; }

size_t KineticDelaunay::getSectionCount() const { return branch_trajs.getHeight(); }

void KineticDelaunay::setSectionRange(size_t start_section, std::optional<size_t> end_section)
{
  start_section_ = start_section;
  end_section_ = end_section;
}

size_t KineticDelaunay::getEndSection() const
{
  const size_t section_count = getSectionCount();
  // Default stop/finalize time is the tree top (height). Section events run on [start, end).
  const size_t default_end = section_count;
  if (!end_section_.has_value())
  {
    return default_end;
  }
  return std::min(end_section_.value(), default_end);
}

// Computes the Delaunay triangulation of the given splines
void KineticDelaunay::compute()
{
  const size_t start_section = getStartSection();
  const size_t end_section = getEndSection();
  if (collect_statistics_)
  {
    statistics_.beginRun();
    size_t total_strands = 0;
    const size_t point_count = branch_trajs.getPoints().size();
    for (size_t strand_id = 0; strand_id < point_count; ++strand_id)
    {
      if (!isDummyBoundary(strand_id))
      {
        ++total_strands;
      }
    }
    size_t total_branches = 0;
    {
      std::unordered_set<size_t> unique_input_branches;
      const auto& branch_indices = branch_trajs.getBranchIndices();
      for (size_t strand_id = 0; strand_id < branch_indices.size(); ++strand_id)
      {
        if (isDummyBoundary(strand_id))
        {
          continue;
        }
        for (size_t branch_id : branch_indices[strand_id])
        {
          unique_input_branches.insert(branch_id);
        }
      }
      total_branches = unique_input_branches.size();
    }
    statistics_.setTotalsTopology(total_strands, total_branches);
  }
  section_event_manager_->computeEvents(static_cast<double>(start_section), static_cast<size_t>(-1));
  enqueueScheduledSubdivisionEvents();
  handleEvents();
  section_event_manager_->finishProgress();

  // Finalize at the exclusive stop time (no section event / kinetic events at this time).
  const double end_time = static_cast<double>(end_section);

  if (callback_manager_)
  {
    if (collect_statistics_)
    {
      using Clock = std::chrono::steady_clock;
      const Clock::time_point finalize_started = Clock::now();
      callback_manager_->finalize(end_time);
      const std::chrono::duration<double> finalize_elapsed = Clock::now() - finalize_started;
      statistics_.addWallTimeSeconds(finalize_elapsed.count());
    }
    else
    {
      callback_manager_->finalize(end_time);
    }
  }
  if (collect_statistics_)
  {
    statistics_.endRun();
  }
}

std::vector<size_t> KineticDelaunay::extractConnectedComponent(size_t u, std::vector<bool>& visited) const
{
  std::vector<size_t> component;

  // Perform an iterative DFS with edges induced by inside faces

  std::vector<size_t> stack;
  stack.push_back(u);

  while (!stack.empty())
  {
    size_t v = stack.back();
    stack.pop_back();

    if (visited[v])
      continue;

    visited[v] = true;
    component.push_back(v);

    const auto nbrs = graph.inducedNeighbors(v, face_inside);

    // Push neighbors in reverse order, the same order as recursive DFS
    for (auto it = nbrs.rbegin(); it != nbrs.rend(); ++it)
    {
      size_t w = *it;
      if (!visited[w])
        stack.push_back(w);
    }
  }

  return component;
}

std::vector<size_t> KineticDelaunay::extractGraphConnectedComponent(size_t u, std::vector<bool>& visited) const
{
  std::vector<size_t> component;
  std::vector<size_t> stack;
  stack.push_back(u);

  while (!stack.empty())
  {
    const size_t v = stack.back();
    stack.pop_back();

    if (visited[v])
    {
      continue;
    }

    visited[v] = true;
    component.push_back(v);

    const auto nbrs = graph.inducedNeighborsFromLiveGraph(v);
    for (auto it = nbrs.rbegin(); it != nbrs.rend(); ++it)
    {
      const size_t w = *it;
      if (!visited[w])
      {
        stack.push_back(w);
      }
    }
  }

  return component;
}

const std::vector<glm::dvec2>& KineticDelaunay::getDummyBoundary() const { return dummy_boundary; }

std::vector<std::vector<size_t>> KineticDelaunay::checkForSplit(const std::array<int, 3>& tri_vertices, double t) const
{
  return checkForSplit(tri_vertices, face_inside, t);
}

std::vector<std::vector<size_t>> KineticDelaunay::checkForSplit(
  const std::array<int, 3>& tri_vertices, const std::vector<bool>& inside_state, double t) const
{
  if (tri_vertices[0] < 0 || tri_vertices[1] < 0 || tri_vertices[2] < 0)
  {
    return {};
  }

  const size_t seed0 = static_cast<size_t>(tri_vertices[0]);
  if (seed0 >= component_data.component_map.size())
  {
    return {};
  }
  const size_t parent_component_id = component_data.component_map[seed0];
  if (parent_component_id >= component_data.components.size())
  {
    return {};
  }

  const std::vector<size_t>& parent_strands = component_data.components[parent_component_id];
  if (parent_strands.size() < 2)
  {
    return {};
  }

  // Match isMinimalInputBranchTriangle: branch membership for the kinetic interval around t.
  const size_t section = static_cast<size_t>(std::ceil(t));
  const auto& strands_by_branch = branch_trajs.getStrandsByBranchId();
  if (section >= strands_by_branch.size())
  {
    return {};
  }

  const auto input_branch_of = [&](size_t strand_id) -> size_t
  { return branch_trajs.getBranchIndex(strand_id, section); };

  // Strands of each input branch that still belong to this kinetic component.
  std::unordered_map<size_t, std::vector<size_t>> strands_by_input_branch;
  strands_by_input_branch.reserve(parent_strands.size());
  for (size_t strand_id : parent_strands)
  {
    if (isDummyBoundary(strand_id))
    {
      continue;
    }
    strands_by_input_branch[input_branch_of(strand_id)].push_back(strand_id);
  }
  if (strands_by_input_branch.size() < 2)
  {
    // A single input branch cannot initiate a pending split.
    return {};
  }

  // Pending-split criterion (component-limited multi-source BFS / DFS):
  // For each input branch in this kinetic component, seed a search from *all* of its strands and walk only
  // parent-component induced neighbors. Trigger isolation only if that search never reaches a different input
  // branch. A lone disconnected strand of a still-mixed branch must not start a pending split / graph cut.
  std::unordered_set<size_t> parent_strand_set(parent_strands.begin(), parent_strands.end());
  std::vector<size_t> isolated_input_branches;
  isolated_input_branches.reserve(strands_by_input_branch.size());
  std::vector<bool> branch_search_visited(graph.getVertexCount(), false);
  for (const auto& [input_branch, branch_strands] : strands_by_input_branch)
  {
    if (branch_strands.empty())
    {
      continue;
    }

    std::fill(branch_search_visited.begin(), branch_search_visited.end(), false);
    std::vector<size_t> queue;
    queue.reserve(branch_strands.size());
    for (size_t strand_id : branch_strands)
    {
      if (strand_id >= branch_search_visited.size())
      {
        continue;
      }
      branch_search_visited[strand_id] = true;
      queue.push_back(strand_id);
    }

    bool found_foreign = false;
    size_t head = 0;
    while (head < queue.size() && !found_foreign)
    {
      const size_t v = queue[head++];
      for (size_t w : graph.inducedNeighbors(v, inside_state))
      {
        if (w >= graph.getVertexCount() || isDummyBoundary(w) || parent_strand_set.count(w) == 0)
        {
          continue;
        }
        if (input_branch_of(w) != input_branch)
        {
          found_foreign = true;
          break;
        }
        if (!branch_search_visited[w])
        {
          branch_search_visited[w] = true;
          queue.push_back(w);
        }
      }
    }

    if (!found_foreign)
    {
      isolated_input_branches.push_back(input_branch);
    }
  }

  if (isolated_input_branches.empty())
  {
    return {};
  }

  std::sort(isolated_input_branches.begin(), isolated_input_branches.end());
  std::unordered_set<size_t> isolated_input_branch_set(
    isolated_input_branches.begin(), isolated_input_branches.end());

  // One pending-split child per isolated input branch (all of that branch's strands in the parent, even if they
  // form multiple induced pieces). Retained parent keeps every non-isolated strand (and dummies). Runtime-branch
  // children for the same input branch are still deduplicated later in notePendingBranchSplit.
  std::vector<size_t> retained_strands;
  retained_strands.reserve(parent_strands.size());
  for (size_t strand_id : parent_strands)
  {
    if (isDummyBoundary(strand_id)
      || isolated_input_branch_set.find(input_branch_of(strand_id)) == isolated_input_branch_set.end())
    {
      retained_strands.push_back(strand_id);
    }
  }

  std::vector<std::vector<size_t>> components;
  components.reserve(isolated_input_branches.size() + 1);
  if (!retained_strands.empty())
  {
    components.push_back(std::move(retained_strands));
  }
  for (size_t input_branch : isolated_input_branches)
  {
    components.push_back(strands_by_input_branch[input_branch]);
  }

  if (components.size() < 2)
  {
    return {};
  }

  // Prefer the piece that contains the radius-triangle seed as components[0] (retained parent id).
  for (size_t i = 0; i < components.size(); ++i)
  {
    if (std::find(components[i].begin(), components[i].end(), seed0) != components[i].end())
    {
      if (i != 0)
      {
        std::swap(components[0], components[i]);
      }
      break;
    }
  }

  return components;
}

std::vector<std::vector<size_t>> KineticDelaunay::extractConnectedComponents() const
{
  std::vector<std::vector<size_t>> components;
  std::vector<bool> visited(graph.getVertexCount(), false);
  for (size_t u = 0; u < graph.getVertexCount(); u++)
  {
    if (visited[u])
    {
      continue;
    }

    auto component = extractConnectedComponent(u, visited);
    components.push_back(component);
  }

  return components;
}

std::vector<std::vector<size_t>> KineticDelaunay::extractGraphConnectedComponents() const
{
  std::vector<std::vector<size_t>> components;
  std::vector<bool> visited(graph.getVertexCount(), false);
  for (size_t u = 0; u < graph.getVertexCount(); ++u)
  {
    if (visited[u] || isDummyBoundary(u) || !isStrandLiveInGraph(u))
    {
      visited[u] = true;
      continue;
    }

    auto component = extractGraphConnectedComponent(u, visited);
    if (!component.empty())
    {
      components.push_back(std::move(component));
    }
  }

  return components;
}

bool KineticDelaunay::isStrandLiveInGraph(size_t strand_id) const
{
  if (isDummyBoundary(strand_id))
  {
    return false;
  }
  return graph.incidentEdgesBegin(strand_id) != graph.incidentEdgesEnd(strand_id);
}

std::pair<size_t, size_t> KineticDelaunay::countLiveStrandsAndBranches() const
{
  size_t strand_count = 0;
  const size_t vertex_count = graph.getVertexCount();
  for (size_t strand_id = 0; strand_id < vertex_count; ++strand_id)
  {
    if (isStrandLiveInGraph(strand_id))
    {
      ++strand_count;
    }
  }

  size_t branch_count = 0;
  for (size_t branch_id = 0; branch_id < runtime_branch_data_.alive.size(); ++branch_id)
  {
    if (runtime_branch_data_.alive[branch_id])
    {
      ++branch_count;
    }
  }
  return { strand_count, branch_count };
}

std::vector<BoundaryPoint> KineticDelaunay::traverseBoundary(
  size_t start_he_id, double t, bool apply_reference_transform, bool include_virtual_offset) const
{
  // Walk the boundary to extract the boundary half-edges
  std::vector<BoundaryPoint> boundary_points;
  size_t he_id = start_he_id;
  do
  {
    size_t origin = graph.halfEdge(he_id).origin;
    if (origin == -1)
    {
      KINDS_ERROR("Followed infinite edge.");
    }
    glm::dvec2 pos = getPointAt(origin, t, apply_reference_transform, include_virtual_offset);
    boundary_points.emplace_back(BoundaryPoint { origin, he_id, pos });
    he_id = nextOnComponentBoundaryId(he_id);
  } while (he_id != start_he_id);

  return boundary_points;
}

std::vector<std::vector<BoundaryPoint>> KineticDelaunay::extractComponentBoundaries(const std::vector<size_t>& component,
  double t, std::vector<bool>& he_visited, bool apply_reference_transform, bool include_virtual_offset) const
{
  // KINDS_DEBUG("Extracting component boundaries at t = " << t);
  if (component.size() < 3)
  {
    return { {} };
  }

  std::vector<std::vector<BoundaryPoint>> boundaries;
  double min_x = std::numeric_limits<double>::infinity();
  // TODO: this is not perfectly safe if points of the outer and an inner boundary coincide at the minimum
  size_t min_x_id = 0;
  for (size_t i = 0; i < component.size(); i++)
  {
    const size_t& v = component[i];

    for (auto it = graph.incidentEdgesBegin(v); it != graph.incidentEdgesEnd(v); it++)
    {
      auto he_id = *it;

      if (he_visited[he_id] || !isOnComponentBoundaryOutside(he_id))
      {
        continue;
      }

      if (graph.destination(he_id) == -1)
      {
        KINDS_ERROR("Destination of half-edge is invalid");
      }
      if (graph.halfEdge(he_id).origin == -1)
      {
        KINDS_ERROR("Origin of half-edge is invalid");
      }

      auto boundary_points = traverseBoundary(he_id, t, apply_reference_transform, include_virtual_offset);

      for (auto& bp : boundary_points)
      {
        he_visited[bp.he_id] = true;
        if (bp.p[0] < min_x)
        {
          min_x = bp.p[0];
          min_x_id = boundaries.size();
        }
      }

      boundaries.emplace_back(boundary_points);
    }
  }

  // swap the boundary with the minimum x to the front
  if (min_x_id != 0)
  {
    std::swap(boundaries[0], boundaries[min_x_id]);
  }

  return boundaries;
}

std::vector<BoundaryPoint> KineticDelaunay::extractComponentBoundary(
  const std::vector<size_t>& component, double t, bool apply_reference_transform, bool include_virtual_offset) const
{
  // Find an extreme point to start the boundary walk as it must be on the boundary
  // Note that merely being on the outside of the boundary is not sufficent as there can also be holes inside the
  // component

  size_t start_vertex_id = -1;
  double min_x = std::numeric_limits<double>::infinity();
  for (size_t i = 0; i < component.size(); i++)
  {
    const size_t& v = component[i];

    // Get position and check if it's the minimum x
    glm::dvec2 pos = getPointAt(v, t, apply_reference_transform, include_virtual_offset);
    if (pos[0] < min_x)
    {
      min_x = pos[0];
      start_vertex_id = v;
    }
  }

  // From the starting vertex, find a half-edge that is on the boundary
  size_t start_he_id = -1;
  for (auto it = graph.incidentEdgesBegin(start_vertex_id); it != graph.incidentEdgesEnd(start_vertex_id); it++)
  {
    if (isOnComponentBoundaryOutside(*it))
    {
      start_he_id = *it;
      break;
    }
  }

  return traverseBoundary(start_he_id, t, apply_reference_transform, include_virtual_offset);
}

const std::vector<bool>& KineticDelaunay::getFacesInside() const { return face_inside; }

bool KineticDelaunay::getFaceInside(size_t face_index) const { return face_inside[face_index]; }

size_t KineticDelaunay::getRuntimeBranchIdForStrand(size_t strand_id) const
{
  if (strand_id >= runtime_branch_data_.branch_map.size())
  {
    throw std::runtime_error(
      "getRuntimeBranchIdForStrand: strand " + std::to_string(strand_id) + " has no runtime branch entry");
  }
  return runtime_branch_data_.branch_map[strand_id];
}

size_t KineticDelaunay::getInsideFaceComponentId(size_t strand_id) const
{
  if (strand_id >= component_data.component_map.size())
  {
    throw std::runtime_error(
      "getInsideFaceComponentId: strand " + std::to_string(strand_id) + " has no component entry");
  }
  return component_data.component_map[strand_id];
}

size_t KineticDelaunay::getComponentLowestStrandId(size_t component_id) const
{
  if (component_id >= component_data.components.size() || component_data.components[component_id].empty())
  {
    throw std::runtime_error("getComponentLowestStrandId: invalid component id " + std::to_string(component_id));
  }
  return component_data.components[component_id].front();
}

size_t KineticDelaunay::getRuntimeBranchIdForFace(size_t face_id) const
{
  constexpr size_t kDefaultRuntimeBranch = 0;

  const HalfEdgeDelaunayGraph& delaunay_graph = getGraph();
  if (!delaunay_graph.isLiveFace(face_id))
  {
    return kDefaultRuntimeBranch;
  }

  const std::array<int, 3> vertices = delaunay_graph.getTriangleVertexIndices(face_id);
  std::optional<size_t> runtime_branch;
  for (int vertex : vertices)
  {
    if (vertex < 0)
    {
      continue;
    }
    if (isDummyBoundary(static_cast<size_t>(vertex)))
    {
      return kDefaultRuntimeBranch;
    }
    if (static_cast<size_t>(vertex) >= runtime_branch_data_.branch_map.size())
    {
      return kDefaultRuntimeBranch;
    }

    const size_t vertex_runtime_branch = runtime_branch_data_.branch_map[static_cast<size_t>(vertex)];
    if (!runtime_branch.has_value())
    {
      runtime_branch = vertex_runtime_branch;
    }
    else if (runtime_branch.value() != vertex_runtime_branch)
    {
      return kDefaultRuntimeBranch;
    }
  }

  return runtime_branch.value_or(kDefaultRuntimeBranch);
}

size_t KineticDelaunay::getRuntimeBranchIdForHalfEdge(size_t half_edge_id) const
{
  constexpr size_t kDefaultRuntimeBranch = 0;

  const HalfEdgeDelaunayGraph& delaunay_graph = getGraph();
  if (!delaunay_graph.isLiveHalfEdge(half_edge_id))
  {
    return kDefaultRuntimeBranch;
  }

  const size_t from_face
    = getRuntimeBranchIdForFace(static_cast<size_t>(delaunay_graph.halfEdge(half_edge_id).face));
  const size_t from_twin
    = getRuntimeBranchIdForFace(static_cast<size_t>(delaunay_graph.halfEdge(half_edge_id ^ 1).face));

  if (from_face == from_twin)
  {
    return from_face;
  }
  if (from_face != kDefaultRuntimeBranch)
  {
    return from_face;
  }
  return from_twin;
}

bool KineticDelaunay::runtimeBranchHasSingleFiniteTriangle(size_t runtime_branch_id) const
{
  const HalfEdgeDelaunayGraph& delaunay_graph = getGraph();
  size_t finite_triangle_count = 0;

  for (size_t live_face_id : delaunay_graph.liveFaces())
  {
    const std::array<int, 3> tri_vertices = delaunay_graph.getTriangleVertexIndices(live_face_id);
    if (tri_vertices[0] == -1 || tri_vertices[1] == -1 || tri_vertices[2] == -1)
    {
      continue;
    }

    if (getRuntimeBranchIdForFace(live_face_id) != runtime_branch_id)
    {
      continue;
    }

    if (++finite_triangle_count > 1)
    {
      return false;
    }
  }

  return finite_triangle_count == 1;
}

bool KineticDelaunay::isMinimalInputBranchTriangle(const std::array<int, 3>& vertices, double t) const
{
  if (vertices[0] == -1 || vertices[1] == -1 || vertices[2] == -1)
  {
    return false;
  }

  for (int v : vertices)
  {
    if (isDummyBoundary(static_cast<size_t>(v)))
    {
      return false;
    }
  }

  const size_t section = static_cast<size_t>(std::ceil(t));
  const auto& strands_by_branch = branch_trajs.getStrandsByBranchId();
  if (section >= strands_by_branch.size())
  {
    return false;
  }

  const size_t branch0 = getBranchIndex(static_cast<size_t>(vertices[0]), section);
  const size_t branch1 = getBranchIndex(static_cast<size_t>(vertices[1]), section);
  const size_t branch2 = getBranchIndex(static_cast<size_t>(vertices[2]), section);
  if (branch0 != branch1 || branch1 != branch2)
  {
    return false;
  }

  if (branch0 >= strands_by_branch[section].size())
  {
    return false;
  }

  return strands_by_branch[section][branch0].size() == 3;
}

bool KineticDelaunay::mustRemainInside(size_t face_index, double t) const
{
  const auto& tri = graph.face(face_index);
  return isMinimalInputBranchTriangle(graph.adjacentTriangleVertices(tri.half_edges[0]), t);
}

void KineticDelaunay::setFaceInside(size_t face_index, bool value, double t)
{
  if (!value && mustRemainInside(face_index, t))
  {
    return;
  }

  if (value)
  {
    auto tri_vertices = graph.adjacentTriangleVertices(graph.face(face_index).half_edges[0]);

    for (int& v : tri_vertices)
    {
      if (v == -1)
      {
        // cannot set face with infinite vertex to inside
        throw std::runtime_error("Cannot set face " + std::to_string(face_index) + " to inside!");
      }
    }
  }
  face_inside[face_index] = value;
}

bool KineticDelaunay::isOnFutureBranchSeamForComponent(
  size_t he_id, size_t component_id, const std::unordered_set<size_t>& partner_component_ids) const
{
  if (!graph.isLiveHalfEdge(he_id))
  {
    return false;
  }

  const int origin = graph.halfEdge(he_id).origin;
  const int destination = graph.destination(he_id);
  if (origin < 0 || destination < 0)
  {
    return false;
  }

  const size_t origin_id = static_cast<size_t>(origin);
  const size_t destination_id = static_cast<size_t>(destination);
  if (origin_id >= component_data.component_map.size()
    || destination_id >= component_data.component_map.size())
  {
    return false;
  }
  if (isDummyBoundary(origin_id) || isDummyBoundary(destination_id))
  {
    return false;
  }

  const size_t origin_component = component_data.component_map[origin_id];
  const size_t destination_component = component_data.component_map[destination_id];
  if (origin_component != component_id)
  {
    return false;
  }
  return partner_component_ids.count(destination_component) > 0;
}

std::vector<size_t> KineticDelaunay::collectFutureBranchSeamStartEdges(
  size_t component_id, const std::unordered_set<size_t>& partner_component_ids) const
{
  std::vector<size_t> start_edges;
  if (component_id >= component_data.components.size())
  {
    return start_edges;
  }

  std::unordered_set<size_t> seen;
  for (size_t strand_id : component_data.components[component_id])
  {
    for (auto it = graph.incidentEdgesBegin(strand_id); it != graph.incidentEdgesEnd(strand_id); ++it)
    {
      // Incident half-edges have their origin at strand_id (inside component_id). An outward seam half-edge crosses
      // into a partner component; this is exactly the tombstoning criterion (endpoints in different components).
      const size_t outward_he = *it;
      if (!isOnFutureBranchSeamForComponent(outward_he, component_id, partner_component_ids))
      {
        continue;
      }

      // The inward-pointing half-edge (destination inside component_id) is the start seed; it is not part of the
      // outline itself.
      const size_t inward_he = graph.twin(outward_he);
      if (seen.insert(inward_he).second)
      {
        start_edges.push_back(inward_he);
      }
    }
  }

  return start_edges;
}

size_t KineticDelaunay::nextOnFutureBranchSeamId(size_t he_id) const
{
  // A seam-outline half-edge has both endpoints in the same runtime branch. Rotating around the pivot vertex (the
  // destination of @p he_id, i.e. the origin of its `next`) yields the next outline edge of that branch.
  const auto within_same_runtime_branch = [&](size_t candidate_he) {
    const size_t origin = static_cast<size_t>(graph.halfEdge(candidate_he).origin);
    const int destination = graph.destination(candidate_he);
    if (destination < 0 || isDummyBoundary(origin) || isDummyBoundary(static_cast<size_t>(destination)))
    {
      return false;
    }
    const size_t origin_branch = runtime_branch_data_.branchForStrand(origin);
    const size_t destination_branch = runtime_branch_data_.branchForStrand(static_cast<size_t>(destination));
    return origin_branch != RuntimeBranchData::no_branch && origin_branch == destination_branch;
  };

  const size_t first_candidate_he_id = graph.halfEdge(he_id).next;
  size_t candidate_he_id = first_candidate_he_id;
  while (!within_same_runtime_branch(candidate_he_id))
  {
    candidate_he_id = graph.neighborEdgeId(candidate_he_id);
    if (candidate_he_id == first_candidate_he_id)
    {
      // Rotated all the way around the pivot vertex without finding a same-branch edge.
      KINDS_DEBUG("nextOnFutureBranchSeamId: no next seam edge from he_id=" << he_id);
      return he_id;
    }
  }
  return candidate_he_id;
}

std::vector<BoundaryPoint> KineticDelaunay::traverseFutureBranchSeam(size_t start_he_id, double t) const
{
  // `start_he_id` is a start seed (an inward-pointing cross-component seam half-edge whose destination is the first
  // outline vertex). The seed itself is not recorded: the walk advances with @ref nextOnFutureBranchSeamId and closes
  // once a half-edge's destination returns to the seed's destination.
  const int start_destination = graph.destination(start_he_id);
  if (start_destination < 0)
  {
    return {};
  }

  std::vector<BoundaryPoint> seam_points;
  size_t current_he_id = start_he_id;
  size_t guard = 0;
  const size_t guard_limit = graph.halfEdgeSlotCount() + 1;
  do
  {
    const size_t next_he_id = nextOnFutureBranchSeamId(current_he_id);
    if (next_he_id == current_he_id || ++guard >= guard_limit)
    {
      return {};
    }
    current_he_id = next_he_id;

    const size_t origin = static_cast<size_t>(graph.halfEdge(current_he_id).origin);
    seam_points.emplace_back(BoundaryPoint { origin, current_he_id, getPointAt(origin, t) });
  } while (graph.destination(current_he_id) != start_destination);

  return seam_points;
}

bool KineticDelaunay::isOnComponentBoundary(size_t he_id) const
{
  size_t face_id = graph.halfEdge(he_id).face;
  size_t twin_face_id = graph.halfEdge(he_id ^ 1).face;
  return (face_inside[face_id] != face_inside[twin_face_id]);
}

bool KineticDelaunay::isOnComponentBoundaryOutside(size_t he_id) const
{
  size_t face_id = graph.halfEdge(he_id).face;
  size_t twin_face_id = graph.halfEdge(he_id ^ 1).face;
  return (!face_inside[face_id] && face_inside[twin_face_id]);
}

size_t KineticDelaunay::nextOnComponentBoundaryId(size_t he_id) const
{
  size_t next_he_id = graph.halfEdge(he_id).next;

  while (!isOnComponentBoundaryOutside(next_he_id))
  {
    next_he_id = graph.twin(next_he_id);
    next_he_id = graph.halfEdge(next_he_id).next;
  }

  return next_he_id;
}

namespace
{
struct FiniteFaceInsideExpectation
{
  bool valid = false;
  double circumradius = 0.0;
  bool expected_inside = false;
  std::array<int, 3> vertices { -1, -1, -1 };
};

FiniteFaceInsideExpectation computeFiniteFaceInsideExpectation(
  const KineticDelaunay& kd, size_t face_id, double t)
{
  FiniteFaceInsideExpectation out;
  const HalfEdgeDelaunayGraph& graph = kd.getGraph();
  if (!graph.isLiveFace(face_id) || face_id >= kd.getFacesInside().size())
  {
    return out;
  }

  out.vertices = graph.getTriangleVertexIndices(face_id);
  if (out.vertices[0] < 0 || out.vertices[1] < 0 || out.vertices[2] < 0)
  {
    return out;
  }
  for (int vertex : out.vertices)
  {
    if (!kd.isDiagnosticsStrandIdValid(static_cast<size_t>(vertex)))
    {
      return out;
    }
  }

  const glm::dvec2 p0 = kd.getPointAt(static_cast<size_t>(out.vertices[0]), t);
  const glm::dvec2 p1 = kd.getPointAt(static_cast<size_t>(out.vertices[1]), t);
  const glm::dvec2 p2 = kd.getPointAt(static_cast<size_t>(out.vertices[2]), t);
  try
  {
    out.circumradius = circumradius(p0, p1, p2);
  }
  catch (const std::exception&)
  {
    return out;
  }

  out.expected_inside = (out.circumradius < kd.getCutoff()) || kd.mustRemainInside(face_id, t);
  out.valid = true;
  return out;
}

std::string formatTriangleVertices(const std::array<int, 3>& vertices)
{
  return "(" + std::to_string(vertices[0]) + "," + std::to_string(vertices[1]) + "," + std::to_string(vertices[2])
    + ")";
}

std::string formatFaceInsideStateValues(
  bool stored_inside, const FiniteFaceInsideExpectation& info, double cutoff)
{
  std::ostringstream oss;
  oss << "(stored_inside=" << stored_inside;
  if (info.valid)
  {
    oss << ", expected_inside=" << info.expected_inside << ", circumradius=" << info.circumradius << ", cutoff="
        << cutoff << ", vertices=" << formatTriangleVertices(info.vertices) << ", must_remain_inside="
        << (info.expected_inside && !(info.circumradius < cutoff)) << ")";
  }
  else
  {
    oss << ", expected_inside=unavailable, circumradius=unavailable, cutoff=" << cutoff
        << ", vertices=" << formatTriangleVertices(info.vertices) << ", must_remain_inside=unavailable)";
  }
  return oss.str();
}

[[noreturn]] void throwIncorrectFaceInsideState(const char* context, size_t face_id, bool stored_inside,
  const FiniteFaceInsideExpectation& info, double cutoff, double t, const std::string& extra_detail = {})
{
  std::ostringstream oss;
  oss << "Face inside/outside sanity check failed";
  if (context != nullptr && context[0] != '\0')
  {
    oss << " (" << context << ")";
  }
  oss << ": face " << face_id << " has incorrect inside state at t=" << t << " "
      << formatFaceInsideStateValues(stored_inside, info, cutoff);
  if (!extra_detail.empty())
  {
    oss << ". " << extra_detail;
  }
  throw std::runtime_error(oss.str());
}

void validateStoredFaceInsideAgainstExpectation(const KineticDelaunay& kd, size_t face_id, double t,
  const char* context, const std::string& extra_detail = {})
{
  if (!kd.getGraph().isLiveFace(face_id) || face_id >= kd.getFacesInside().size())
  {
    return;
  }

  const FiniteFaceInsideExpectation info = computeFiniteFaceInsideExpectation(kd, face_id, t);
  if (!info.valid)
  {
    return;
  }

  const bool stored_inside = kd.getFaceInside(face_id);
  if (stored_inside != info.expected_inside)
  {
    throwIncorrectFaceInsideState(context, face_id, stored_inside, info, kd.getCutoff(), t, extra_detail);
  }
}
} // namespace

void KineticDelaunay::validateFlipAdjacentFaceInsideConsistency(size_t half_edge_id, double t) const
{
  if (!isDiagnosticsHalfEdgeIdValid(half_edge_id) || !isDiagnosticsHalfEdgeIdValid(half_edge_id ^ 1))
  {
    return;
  }

  const size_t face_a = graph.halfEdge(half_edge_id).face;
  const size_t face_b = graph.halfEdge(half_edge_id ^ 1).face;
  if (face_a >= face_inside.size() || face_b >= face_inside.size())
  {
    return;
  }

  const bool inside_a = face_inside[face_a];
  const bool inside_b = face_inside[face_b];
  if (inside_a == inside_b)
  {
    return;
  }

  const FiniteFaceInsideExpectation info_a = computeFiniteFaceInsideExpectation(*this, face_a, t);
  const FiniteFaceInsideExpectation info_b = computeFiniteFaceInsideExpectation(*this, face_b, t);
  if (!info_a.valid || !info_b.valid)
  {
    return;
  }

  const bool incorrect_a = inside_a != info_a.expected_inside;
  const bool incorrect_b = inside_b != info_b.expected_inside;

  std::ostringstream detail;
  detail << "flip half-edge " << half_edge_id << " adjacent faces " << face_a << " (inside=" << inside_a << ", r="
         << info_a.circumradius << ", expected=" << info_a.expected_inside << ") and " << face_b << " (inside="
         << inside_b << ", r=" << info_b.circumradius << ", expected=" << info_b.expected_inside << ")";

  if (incorrect_a && incorrect_b)
  {
    detail << "; both faces disagree with circumradius/cutoff expectation";
    throwIncorrectFaceInsideState("flip_event", face_a, inside_a, info_a, getCutoff(), t, detail.str());
  }
  if (incorrect_a)
  {
    throwIncorrectFaceInsideState("flip_event", face_a, inside_a, info_a, getCutoff(), t, detail.str());
  }
  if (incorrect_b)
  {
    throwIncorrectFaceInsideState("flip_event", face_b, inside_b, info_b, getCutoff(), t, detail.str());
  }

  std::ostringstream oss;
  oss << "Face inside/outside sanity check failed (flip_event): flip half-edge " << half_edge_id << " at t=" << t
      << " has adjacent faces with mismatched inside status (face " << face_a << " inside=" << inside_a << ", r="
      << info_a.circumradius << "; face " << face_b << " inside=" << inside_b << ", r=" << info_b.circumradius
      << ", cutoff=" << getCutoff()
      << ") but neither face violates circumradius/cutoff expectation (possible mustRemainInside asymmetry)";
  throw std::runtime_error(oss.str());
}

void KineticDelaunay::validateAllFaceInsideStatesAtTime(double t, const char* context) const
{
  for (size_t face_id : graph.liveFaces())
  {
    validateStoredFaceInsideAgainstExpectation(*this, face_id, t, context);
  }
}

void KineticDelaunay::logFaceInsideStateAtTime(size_t face_id, double t, const char* context) const
{
  if (!isDiagnosticsFaceIdValid(face_id))
  {
    return;
  }

  std::ostringstream oss;
  oss << "Face inside/outside sanity check monitor";
  if (context != nullptr && context[0] != '\0')
  {
    oss << " (" << context << ")";
  }

  if (!graph.isLiveFace(face_id))
  {
    oss << ": face " << face_id << " at t=" << t << " (not live)";
    KINDS_MONITOR(oss.str());
    return;
  }

  if (face_id >= face_inside.size())
  {
    oss << ": face " << face_id << " at t=" << t << " (face_inside slot missing)";
    KINDS_MONITOR(oss.str());
    return;
  }

  const FiniteFaceInsideExpectation info = computeFiniteFaceInsideExpectation(*this, face_id, t);
  const bool stored_inside = face_inside[face_id];
  oss << ": face " << face_id << " has inside state at t=" << t << " "
      << formatFaceInsideStateValues(stored_inside, info, getCutoff());
  if (info.valid && stored_inside != info.expected_inside)
  {
    oss << " MISMATCH";
  }
  KINDS_MONITOR(oss.str());
}
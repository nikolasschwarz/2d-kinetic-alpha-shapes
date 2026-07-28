#pragma once

#include "KineticDelaunay.hpp"
#include "KineticDelaunayRadiusEvent.hpp"
#include "KineticDelaunayFlipEventTriggerDump.hpp"

#include <cmath>
#include <limits>
#include <map>
#include <algorithm>
#include <optional>
#include <string>

namespace kinDS
{
inline constexpr double flip_voronoi_vertex_distance_eps = 1e-6;
inline constexpr double flip_boundary_collinearity_eps = 1e-8;

inline double normalizedTriangleCollinearityMetric(const glm::dvec2& pa, const glm::dvec2& pb, const glm::dvec2& pc)
{
  const glm::dvec2 ab = pb - pa;
  const glm::dvec2 ac = pc - pa;
  const double area2 = std::abs(ab.x * ac.y - ab.y * ac.x);
  const double max_len2 = std::max({ glm::dot(ab, ab), glm::dot(ac, ac), glm::dot(pc - pb, pc - pb) });
  if (max_len2 <= 0.0)
  {
    return 0.0;
  }
  return area2 / max_len2;
}

inline std::vector<size_t> collectFlipQuadrilateralStrandIds(const HalfEdgeDelaunayGraph& graph, size_t he_id)
{
  std::vector<size_t> strand_ids;
  strand_ids.reserve(4);

  const auto append_strand = [&](int vertex)
  {
    if (vertex >= 0)
    {
      const size_t strand_id = static_cast<size_t>(vertex);
      if (std::find(strand_ids.begin(), strand_ids.end(), strand_id) == strand_ids.end())
      {
        strand_ids.push_back(strand_id);
      }
    }
  };

  if (graph.isOnConvexBoundary(he_id) || graph.isOutsideConvexBoundary(he_id))
  {
    size_t boundary_he_id = he_id;
    if (graph.isOutsideConvexBoundary(boundary_he_id))
    {
      boundary_he_id ^= 1;
    }

    const int indices[4] = {
      graph.halfEdge(boundary_he_id).origin,
      graph.triangleOppositeVertex(boundary_he_id ^ 1),
      graph.halfEdge(boundary_he_id ^ 1).origin,
      graph.triangleOppositeVertex(boundary_he_id),
    };
    for (int vertex : indices)
    {
      append_strand(vertex);
    }
  }
  else
  {
    append_strand(graph.halfEdge(he_id).origin);
    append_strand(graph.triangleOppositeVertex(he_id ^ 1));
    append_strand(graph.halfEdge(he_id ^ 1).origin);
    append_strand(graph.triangleOppositeVertex(he_id));
  }

  return strand_ids;
}

inline glm::dvec2 flipTriangleCircumcenterAt(const KineticDelaunay& kd, const HalfEdgeDelaunayGraph& graph, size_t he_id,
  double t, bool use_transformed_points, std::optional<size_t> shared_reference_branch = std::nullopt)
{
  const std::array<int, 3> vertices = graph.adjacentTriangleVertices(he_id);
  const StrandTree& tree = kd.getStrandTree();
  const auto point_at = [&](int vertex) -> glm::dvec2
  {
    const size_t strand_id = static_cast<size_t>(vertex);
    if (!use_transformed_points)
    {
      return tree.evaluate(strand_id, t);
    }
    if (shared_reference_branch.has_value())
    {
      return kd.getPointAtWithReferenceBranch(strand_id, t, shared_reference_branch.value());
    }
    return kd.getPointAt(strand_id, t);
  };
  return HalfEdgeDelaunayGraph::circumcenter(
    point_at(vertices[0]), point_at(vertices[1]), point_at(vertices[2]));
}

inline std::string flipUntransformedFrameMismatchNote(bool transformed_passes, bool untransformed_passes)
{
  if (!transformed_passes && untransformed_passes)
  {
    return " [untransformed geometry satisfies the flip criterion; possible frame mismatch in event polynomials]";
  }
  if (transformed_passes && !untransformed_passes)
  {
    return " [only transformed geometry satisfies the flip criterion]";
  }
  return {};
}

class KineticDelaunay::FlipEvent final : public KineticDelaunay::Event
{
 public:
  KineticDelaunay* kd_;
  size_t half_edge_id;
  glm::dvec2 position;

  FlipEvent(KineticDelaunay* kd, double t, size_t he_id, double creation_time, glm::dvec2 position)
    : KineticDelaunay::Event(t, creation_time, 20u)
    , kd_(kd)
    , half_edge_id(he_id)
    , position(position)
  {
  }

  void handleEvent() override;
};

class KineticDelaunay::FlipEventManager final : public KineticDelaunay::EventManager
{
 public:
  explicit FlipEventManager(KineticDelaunay* kd)
    : kd_(kd)
  {
  }

  void computeEvents(double t, size_t event_id) override;

 private:
  KineticDelaunay* kd_;
};

} // namespace kinDS

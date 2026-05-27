#pragma once

#include "HalfEdgeDelaunayGraph.hpp"

#include <cstddef>
#include <cstdint>
#include <optional>
#include <functional>
#include <unordered_set>

namespace kinDS
{
class KineticDelaunay;

/// Primitives to emphasize in event-triggered debug SVG exports (see @ref HalfEdgeDelaunayGraphToSVG::write).
struct VisualDebugHighlight
{
  std::unordered_set<size_t> delaunay_vertices;
  std::unordered_set<size_t> delaunay_faces;
  std::unordered_set<size_t> directed_half_edges;
  std::unordered_set<size_t> voronoi_vertices;
  /// Undirected dual / Voronoi edge id (= even directed Delaunay half-edge index / 2).
  std::unordered_set<size_t> voronoi_edges;
  /// Subset of @ref voronoi_edges drawn with a stronger focus style (e.g. the flip edge in a quad).
  std::unordered_set<size_t> primary_voronoi_edges;
  /// Label every intersection whose delaunay_edge_id is in this set (see @ref shouldLabelCrossing).
  std::unordered_set<size_t> label_crossings_on_delaunay_edges;
  /// (d,v) pairs drawn with a stronger intersection marker (@ref emphasizesCrossing).
  std::unordered_set<uint64_t> crossing_intersection_keys;

  /// Queue-driven; does not recurse on the C++ call stack.
  void addUndirectedDelaunayEdge(const HalfEdgeDelaunayGraph& graph, size_t he_id);
  void addDelaunayTriangle(const HalfEdgeDelaunayGraph& graph, size_t face_id);

  bool affectsDelaunayVertex(size_t vertex_id) const;
  bool affectsDelaunayFace(size_t face_id) const;
  bool affectsDirectedHalfEdge(size_t he_id) const;
  bool affectsVoronoiVertex(size_t voronoi_vertex_id) const;
  bool affectsVoronoiEdge(size_t voronoi_edge_id) const;
  bool affectsPrimaryVoronoiEdge(size_t voronoi_edge_id) const;
  bool shouldLabelCrossing(size_t delaunay_edge_id, size_t voronoi_edge_id) const;
  bool emphasizesCrossing(size_t delaunay_edge_id, size_t voronoi_edge_id) const;

  static VisualDebugHighlight forFlip(const HalfEdgeDelaunayGraph& graph, size_t flip_half_edge_id);
  static VisualDebugHighlight forRadius(const HalfEdgeDelaunayGraph& graph, size_t radius_half_edge_id);
  static VisualDebugHighlight forCrossing(
    const HalfEdgeDelaunayGraph& graph, size_t crossed_half_edge_id, size_t voronoi_vertex_id);
  static VisualDebugHighlight forSubdivisionStrand(const HalfEdgeDelaunayGraph& graph, size_t strand_vertex_id);
  static VisualDebugHighlight forSectionBoundary(
    const HalfEdgeDelaunayGraph& graph, const std::function<bool(size_t even_half_edge_id)>& is_boundary_edge);

  /**
   * Highlight invariant-failure context: at most one @p primary_dual_edge (magenta), label all intersections on that
   * edge, and emphasize (larger marker) the @p crossing_intersection_keys failure pair.
   */
  static VisualDebugHighlight forInvariantViolation(const HalfEdgeDelaunayGraph& graph,
    std::optional<size_t> primary_dual_edge, const std::unordered_set<size_t>& auxiliary_dual_edges,
    const std::unordered_set<uint64_t>& crossing_intersection_keys);
};

} // namespace kinDS

#pragma once

#include "HalfEdgeDelaunayGraphToSVG.hpp"

#include <cstddef>
#include <optional>
#include <string>
#include <vector>

namespace kinDS
{
/// Snapshot of a failed closing-mesh cap walk, anchored to graph elements so it overlays the kinetic debug SVG.
struct ClosingMeshWalkDebug
{
  using Anchor = HalfEdgeDelaunayGraphToSVG::OverlayAnchor;
  enum class StripState
  {
    Unused, ///< Never walked; still open at the failure.
    Walked, ///< Walked by the failing attempt (see @ref walk_strip_order).
    Closed, ///< Consumed by an earlier, successfully closed polygon.
  };
  struct Strip
  {
    size_t index = 0;
    Anchor start;
    Anchor end;
    std::string start_label;
    std::string end_label;
    std::string label;
    StripState state = StripState::Unused;
  };
  struct WalkPoint
  {
    Anchor anchor;
    std::string label;
  };
  struct BoundaryEdge
  {
    size_t half_edge = 0;
    size_t origin = 0;
    size_t destination = 0;
  };

  size_t strand_id = 0;
  double t = 0.0;
  size_t seed_index = 0;
  std::string reason;
  std::vector<Strip> strips;
  /// Strip indices in the order the failing attempt walked them.
  std::vector<size_t> walk_strip_order;
  /// Cap outline vertices in walk order.
  std::vector<WalkPoint> walk_points;
  /// Component-boundary Delaunay half-edges traversed by the failing attempt, in order.
  std::vector<BoundaryEdge> boundary_edges;
  std::optional<WalkPoint> stuck;
};

/// Overlay for @ref HalfEdgeDelaunayGraphToSVG::write: strips colored by state (unused red / walked blue / closed
/// gray) labelled s<i>, the walk as a green path with numbered steps, traversed boundary half-edges (orange, b<k>),
/// the stuck point (magenta), and notes with the reason, walk order, strip list and step legend.
HalfEdgeDelaunayGraphToSVG::DebugOverlay buildClosingMeshWalkOverlay(const ClosingMeshWalkDebug& debug);
} // namespace kinDS

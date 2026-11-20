#pragma once

#include <glm/glm.hpp>
#include <optional>
#include <vector>

namespace kinDS
{
/// Topological identity of an outline vertex or edge support: a kind tag followed by sorted site ids.
using OutlineKey = std::vector<int>;

/// Outline vertex keys of a Voronoi cell outline (see @ref outlineEdgeSupports).
OutlineKey outlineSiteKey(int site);
OutlineKey outlineVoronoiVertexKey(int a, int b, int c);
OutlineKey outlineCrossingKey(int delaunay_a, int delaunay_b, int voronoi_a, int voronoi_b);

/// Support of each edge (i, i+1) of the cell outline of @p strand with vertex @p keys: the Voronoi edge (strand, n) or
/// the boundary Delaunay edge (p, q) both endpoints lie on. Edges without a unique common support get a key that
/// matches no other edge (kind -1).
std::vector<OutlineKey> outlineEdgeSupports(const std::vector<OutlineKey>& keys, int strand);

/// Closed CCW outline: vertex keys and positions, plus the support key of each edge (i, i+1).
struct KeyedOutline
{
  std::vector<OutlineKey> keys;
  std::vector<glm::dvec2> positions;
  std::vector<OutlineKey> edge_supports;
};

/// Vertex of a difference face at lerp(vertex @c a of outline @c outline_a, vertex @c b of outline @c outline_b, @c s)
/// (outline 0 = old, 1 = new); @c s is 0 for outline vertices.
struct OutlineDifferenceVertex
{
  OutlineKey key;
  glm::dvec2 position;
  int outline_a = 0;
  size_t a = 0;
  int outline_b = 0;
  size_t b = 0;
  double s = 0.0;
};

/// CCW face of the symmetric difference; @c in_old: part of old \ new, otherwise new \ old.
struct OutlineDifferenceFace
{
  std::vector<OutlineDifferenceVertex> ring;
  bool in_old = false;
};

/// Vertex inserted into edge (@c edge, @c edge + 1) of outline @c outline (0 = old, 1 = new): a vertex of the other
/// outline on that edge, or a crossing of the outlines. Meshes built along the outlines must split that edge there.
struct OutlineEdgeSplit
{
  int outline = 0;
  size_t edge = 0;
  OutlineDifferenceVertex vertex;
};

struct OutlineDifference
{
  std::vector<OutlineDifferenceFace> faces;
  std::vector<OutlineEdgeSplit> splits;
};

/// Symmetric difference of two outlines matched by topology: shared vertex keys and equal edge supports decide where
/// the outlines coincide; only crossings between edges of different supports are computed numerically. No faces when
/// the outlines coincide; nullopt when the outlines share fewer than two vertices or visit them in different orders.
/// A point visited twice by one outline (pinched at a site) is not shared; faces may then repeat that point.
std::optional<OutlineDifference> outlineSymmetricDifference(
  const KeyedOutline& old_outline, const KeyedOutline& new_outline);
} // namespace kinDS

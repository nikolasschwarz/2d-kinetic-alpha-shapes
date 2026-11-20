#include "OutlineDifference.hpp"

#include <algorithm>
#include <cmath>
#include <map>
#include <set>

using namespace kinDS;

namespace
{
constexpr int kind_site = 0;
constexpr int kind_voronoi_vertex = 1;
constexpr int kind_crossing = 2;
constexpr int kind_outline_crossing = 3;
constexpr int support_none = -1;
constexpr int support_voronoi = 10;
constexpr int support_delaunay = 11;

OutlineKey sortedKey(int kind, std::vector<int> sites)
{
  std::sort(sites.begin(), sites.end());
  sites.insert(sites.begin(), kind);
  return sites;
}

double cross2(const glm::dvec2& a, const glm::dvec2& b) { return a.x * b.y - a.y * b.x; }

struct RingVertex
{
  OutlineDifferenceVertex vertex;
  OutlineKey support; ///< Support of the edge to the next ring vertex.
  size_t edge = 0; ///< Original outline edge the edge to the next ring vertex lies on.
};
using Ring = std::vector<RingVertex>;

struct Insertion
{
  size_t edge;
  double u;
  OutlineDifferenceVertex vertex;
};

Ring makeRing(const KeyedOutline& outline, int index)
{
  Ring ring;
  for (size_t i = 0; i < outline.keys.size(); ++i)
  {
    ring.push_back({ { outline.keys[i], outline.positions[i], index, i, index, i, 0.0 }, outline.edge_supports[i], i });
  }
  return ring;
}

bool supportsMatch(const OutlineKey& a, const OutlineKey& b) { return a[0] != support_none && a == b; }

void applyInsertions(Ring& ring, std::vector<Insertion> insertions, int outline, std::vector<OutlineEdgeSplit>& splits)
{
  std::sort(insertions.begin(), insertions.end(),
    [](const Insertion& x, const Insertion& y) { return x.edge != y.edge ? x.edge < y.edge : x.u < y.u; });
  Ring out;
  auto it = insertions.begin();
  for (size_t i = 0; i < ring.size(); ++i)
  {
    out.push_back(ring[i]);
    for (; it != insertions.end() && it->edge == i; ++it)
    {
      if (out.back().vertex.key != it->vertex.key)
      {
        out.push_back({ it->vertex, ring[i].support, ring[i].edge });
        splits.push_back({ outline, ring[i].edge, it->vertex });
      }
    }
  }
  ring = std::move(out);
}

std::set<OutlineKey> ringKeys(const Ring& ring)
{
  std::set<OutlineKey> keys;
  for (const RingVertex& v : ring)
  {
    keys.insert(v.vertex.key);
  }
  return keys;
}

/// Points of @p other on an edge of @p ring with the same support, not yet in @p ring, ordered along that edge.
std::vector<Insertion> collinearInsertions(const Ring& ring, const Ring& other)
{
  const std::set<OutlineKey> keys = ringKeys(ring);
  std::vector<Insertion> insertions;
  for (size_t i = 0; i < ring.size(); ++i)
  {
    const glm::dvec2 p0 = ring[i].vertex.position;
    const glm::dvec2 d = ring[(i + 1) % ring.size()].vertex.position - p0;
    const double len2 = glm::dot(d, d);
    for (size_t j = 0; j < other.size() && len2 > 0.0; ++j)
    {
      if (!supportsMatch(ring[i].support, other[j].support))
      {
        continue;
      }
      for (const RingVertex* q : { &other[j], &other[(j + 1) % other.size()] })
      {
        const double u = glm::dot(q->vertex.position - p0, d) / len2;
        const bool listed = std::any_of(
          insertions.begin(), insertions.end(), [&](const Insertion& ins) { return ins.vertex.key == q->vertex.key; });
        if (u > 0.0 && u < 1.0 && keys.count(q->vertex.key) == 0 && !listed)
        {
          insertions.push_back({ i, u, q->vertex });
        }
      }
    }
  }
  return insertions;
}
} // namespace

OutlineKey kinDS::outlineSiteKey(int site) { return { kind_site, site }; }

OutlineKey kinDS::outlineVoronoiVertexKey(int a, int b, int c) { return sortedKey(kind_voronoi_vertex, { a, b, c }); }

OutlineKey kinDS::outlineCrossingKey(int delaunay_a, int delaunay_b, int voronoi_a, int voronoi_b)
{
  const OutlineKey d = sortedKey(kind_crossing, { delaunay_a, delaunay_b });
  return { kind_crossing, d[1], d[2], std::min(voronoi_a, voronoi_b), std::max(voronoi_a, voronoi_b) };
}

std::vector<OutlineKey> kinDS::outlineEdgeSupports(const std::vector<OutlineKey>& keys, int strand)
{
  // Supports a vertex lies on within the strand's cell; a site lies on every boundary Delaunay edge it starts.
  const auto vertex_supports = [&](const OutlineKey& key)
  {
    std::vector<OutlineKey> supports;
    if (key[0] == kind_voronoi_vertex && std::find(key.begin() + 1, key.end(), strand) != key.end())
    {
      for (size_t k = 1; k < key.size(); ++k)
      {
        if (key[k] != strand)
        {
          supports.push_back(sortedKey(support_voronoi, { strand, key[k] }));
        }
      }
    }
    else if (key[0] == kind_crossing)
    {
      supports.push_back({ support_delaunay, key[1], key[2] });
      supports.push_back({ support_voronoi, key[3], key[4] });
    }
    return supports;
  };
  const auto on_site = [&](const OutlineKey& key, const OutlineKey& support)
  { return key[0] == kind_site && support[0] == support_delaunay && (support[1] == key[1] || support[2] == key[1]); };

  std::vector<OutlineKey> supports;
  for (size_t i = 0; i < keys.size(); ++i)
  {
    const OutlineKey& a = keys[i];
    const OutlineKey& b = keys[(i + 1) % keys.size()];
    std::vector<OutlineKey> common;
    for (const OutlineKey& support : vertex_supports(a))
    {
      const std::vector<OutlineKey> b_supports = vertex_supports(b);
      if (on_site(b, support) || std::find(b_supports.begin(), b_supports.end(), support) != b_supports.end())
      {
        common.push_back(support);
      }
    }
    for (const OutlineKey& support : vertex_supports(b))
    {
      if (on_site(a, support))
      {
        common.push_back(support);
      }
    }
    supports.push_back(common.size() == 1 ? common.front() : OutlineKey { support_none, static_cast<int>(i) });
  }
  return supports;
}

std::optional<OutlineDifference> kinDS::outlineSymmetricDifference(
  const KeyedOutline& old_outline, const KeyedOutline& new_outline)
{
  OutlineDifference result;
  Ring rings[2] = { makeRing(old_outline, 0), makeRing(new_outline, 1) };
  if (rings[0].size() < 3 || rings[1].size() < 3)
  {
    return std::nullopt;
  }
  // A point an outline visits twice (pinched at its site) is no shared vertex: each occurrence gets its own key.
  for (Ring& ring : rings)
  {
    std::map<OutlineKey, int> count;
    for (const RingVertex& v : ring)
    {
      ++count[v.vertex.key];
    }
    std::map<OutlineKey, int> seen;
    for (RingVertex& v : ring)
    {
      if (count[v.vertex.key] > 1)
      {
        const int occurrence = seen[v.vertex.key]++;
        v.vertex.key.push_back(-1 - occurrence);
      }
    }
  }

  // Topology: points of one outline on a same-support edge of the other become shared vertices.
  {
    std::vector<Insertion> into_old = collinearInsertions(rings[0], rings[1]);
    std::vector<Insertion> into_new = collinearInsertions(rings[1], rings[0]);
    applyInsertions(rings[0], std::move(into_old), 0, result.splits);
    applyInsertions(rings[1], std::move(into_new), 1, result.splits);
  }

  // Proper crossings between edges of different supports that share no endpoint.
  {
    constexpr double eps = 1e-12;
    std::vector<Insertion> insertions[2];
    for (size_t i = 0; i < rings[0].size(); ++i)
    {
      const RingVertex& a0 = rings[0][i];
      const RingVertex& a1 = rings[0][(i + 1) % rings[0].size()];
      for (size_t j = 0; j < rings[1].size(); ++j)
      {
        const RingVertex& b0 = rings[1][j];
        const RingVertex& b1 = rings[1][(j + 1) % rings[1].size()];
        if (supportsMatch(a0.support, b0.support) || a0.vertex.key == b0.vertex.key || a0.vertex.key == b1.vertex.key
          || a1.vertex.key == b0.vertex.key || a1.vertex.key == b1.vertex.key)
        {
          continue;
        }
        const glm::dvec2 r = a1.vertex.position - a0.vertex.position;
        const glm::dvec2 s = b1.vertex.position - b0.vertex.position;
        const double den = cross2(r, s);
        if (std::abs(den) <= eps * glm::length(r) * glm::length(s))
        {
          continue;
        }
        const glm::dvec2 qp = b0.vertex.position - a0.vertex.position;
        const double u = cross2(qp, s) / den;
        const double w = cross2(qp, r) / den;
        if (u <= eps || u >= 1.0 - eps || w <= eps || w >= 1.0 - eps)
        {
          continue;
        }
        const OutlineKey& lo = std::min(a0.support, b0.support);
        const OutlineKey& hi = std::max(a0.support, b0.support);
        OutlineKey key { kind_outline_crossing };
        key.insert(key.end(), lo.begin(), lo.end());
        key.insert(key.end(), hi.begin(), hi.end());
        const OutlineDifferenceVertex vertex { key, a0.vertex.position + u * r, a0.vertex.outline_a, a0.vertex.a,
          a1.vertex.outline_a, a1.vertex.a, u };
        insertions[0].push_back({ i, u, vertex });
        insertions[1].push_back({ j, w, vertex });
      }
    }
    applyInsertions(rings[0], std::move(insertions[0]), 0, result.splits);
    applyInsertions(rings[1], std::move(insertions[1]), 1, result.splits);
  }

  std::map<OutlineKey, size_t> new_index;
  for (size_t j = 0; j < rings[1].size(); ++j)
  {
    if (!new_index.emplace(rings[1][j].vertex.key, j).second)
    {
      return std::nullopt;
    }
  }
  std::vector<size_t> shared_old;
  for (size_t i = 0; i < rings[0].size(); ++i)
  {
    if (new_index.count(rings[0][i].vertex.key) != 0)
    {
      shared_old.push_back(i);
    }
  }
  if (shared_old.size() < 2 || ringKeys(rings[0]).size() != rings[0].size())
  {
    return std::nullopt;
  }

  // Between consecutive shared vertices u -> w both outlines run along chains; where they differ, the chains bound a
  // face. Its orientation (old forward, new backward) tells which side it is on.
  double scale = 0.0;
  for (const Ring& ring : rings)
  {
    for (const RingVertex& v : ring)
    {
      scale = std::max({ scale, std::abs(v.vertex.position.x), std::abs(v.vertex.position.y) });
    }
  }
  std::vector<OutlineDifferenceFace>& faces = result.faces;
  const size_t n_old = rings[0].size();
  const size_t n_new = rings[1].size();
  for (size_t k = 0; k < shared_old.size(); ++k)
  {
    const size_t u_old = shared_old[k];
    const size_t w_old = shared_old[(k + 1) % shared_old.size()];
    const size_t u_new = new_index.at(rings[0][u_old].vertex.key);
    size_t w_new = (u_new + 1) % n_new;
    while (new_index.count(rings[1][w_new].vertex.key) == 0
      || std::find_if(shared_old.begin(), shared_old.end(),
           [&](size_t i) { return rings[0][i].vertex.key == rings[1][w_new].vertex.key; })
        == shared_old.end())
    {
      w_new = (w_new + 1) % n_new;
    }
    if (rings[1][w_new].vertex.key != rings[0][w_old].vertex.key)
    {
      return std::nullopt;
    }
    const size_t old_len = (w_old + n_old - u_old) % n_old;
    const size_t new_len = (w_new + n_new - u_new) % n_new;
    if (old_len == 1 && new_len == 1)
    {
      continue;
    }
    OutlineDifferenceFace face;
    for (size_t i = u_old; i != w_old; i = (i + 1) % n_old)
    {
      face.ring.push_back(rings[0][i].vertex);
    }
    for (size_t j = w_new; j != u_new; j = (j + n_new - 1) % n_new)
    {
      face.ring.push_back(rings[1][j].vertex);
    }
    double area2 = 0.0;
    for (size_t i = 0; i < face.ring.size(); ++i)
    {
      area2 += cross2(face.ring[i].position, face.ring[(i + 1) % face.ring.size()].position);
    }
    if (std::abs(area2) <= 1e-12 * scale * scale)
    {
      continue;
    }
    face.in_old = area2 > 0.0;
    if (!face.in_old)
    {
      std::reverse(face.ring.begin(), face.ring.end());
    }
    faces.push_back(std::move(face));
  }
  return result;
}

#include "kinDS/HalfEdgeDelaunayGraph.hpp"

#include <catch2/catch_test_macros.hpp>

#include <array>
#include <cmath>
#include <map>
#include <set>
#include <vector>

using namespace kinDS;

static std::array<int, 3> face_vertex_key(
  size_t face_id, const HalfEdgeDelaunayGraph& graph)
{
  const auto& tri = graph.face(face_id);
  std::array<int, 3> key;
  for (int i = 0; i < 3; ++i)
  {
    key[i] = graph.halfEdge(tri.half_edges[i]).origin;
  }
  std::sort(key.begin(), key.end());
  return key;
}

static bool half_edge_belongs_to_face(const HalfEdgeDelaunayGraph& graph, size_t he_id, size_t face_id)
{
  if (!graph.isLiveFace(face_id))
  {
    return false;
  }
  const size_t start = graph.face(face_id).half_edges[0];
  size_t cur = start;
  for (int i = 0; i < 3; ++i)
  {
    if (cur == he_id)
    {
      return true;
    }
    cur = static_cast<size_t>(graph.halfEdge(cur).next);
  }
  return false;
}

static bool face_has_infinite_vertex(const HalfEdgeDelaunayGraph& graph, size_t face_id)
{
  const auto verts = graph.getTriangleVertexIndices(face_id);
  return verts[0] == -1 || verts[1] == -1 || verts[2] == -1;
}

static void check_live_half_edge_face_references(const HalfEdgeDelaunayGraph& graph)
{
  for (size_t he_id : graph.liveDelaunayEdges())
  {
    for (size_t directed_he : { he_id, he_id ^ 1 })
    {
      REQUIRE(graph.isLiveHalfEdge(directed_he));
      const int face_id = graph.halfEdge(directed_he).face;
      REQUIRE(face_id >= 0);
      REQUIRE(static_cast<size_t>(face_id) < graph.faceSlotCount());
      REQUIRE(graph.isLiveFace(static_cast<size_t>(face_id)));
      REQUIRE(half_edge_belongs_to_face(graph, directed_he, static_cast<size_t>(face_id)));
    }
  }

  for (size_t he_id : graph.liveDelaunayEdges())
  {
    if (!graph.isOnConvexBoundary(he_id))
    {
      continue;
    }
    const int face0 = graph.halfEdge(he_id).face;
    const int face1 = graph.halfEdge(he_id ^ 1).face;
    const bool inf0 = face_has_infinite_vertex(graph, static_cast<size_t>(face0));
    const bool inf1 = face_has_infinite_vertex(graph, static_cast<size_t>(face1));
    REQUIRE(inf0 != inf1);
  }
}

static void check_remap_invariants(const HalfEdgeDelaunayGraph& graph)
{
  REQUIRE(graph.halfEdgeSlotCount() > 0);
  REQUIRE(graph.faceSlotCount() > 0);
  REQUIRE(graph.liveFaceCount() > 0);
  REQUIRE(graph.liveDelaunayEdgeCount() > 0);

  std::set<size_t> face_ids;

  for (size_t ti : graph.liveFaces())
  {
    REQUIRE(face_ids.insert(ti).second);

    const HalfEdgeDelaunayGraph::Triangle& tri = graph.face(ti);
    for (int i = 0; i < 3; ++i)
    {
      const size_t he_id = tri.half_edges[i];
      REQUIRE(he_id < graph.halfEdgeSlotCount());
      REQUIRE(graph.halfEdge(he_id).face == static_cast<int>(ti));
      REQUIRE(graph.isLiveHalfEdge(he_id));
    }
  }

  for (size_t he_id : graph.liveDelaunayEdges())
  {
    REQUIRE(graph.isLiveHalfEdge(he_id));
    REQUIRE(graph.isLiveHalfEdge(he_id ^ 1));
    const auto& he = graph.halfEdge(he_id);
    REQUIRE(he.next >= 0);
    REQUIRE(static_cast<size_t>(he.next) < graph.halfEdgeSlotCount());
    REQUIRE(graph.isLiveHalfEdge(static_cast<size_t>(he.next)));
    if (he.face >= 0)
    {
      REQUIRE(static_cast<size_t>(he.face) < graph.faceSlotCount());
      REQUIRE(graph.isLiveFace(static_cast<size_t>(he.face)));
    }
  }

  for (size_t face_id : graph.liveFaces())
  {
    const auto& tri = graph.face(face_id);
    size_t he_id = tri.half_edges[0];
    for (int i = 0; i < 3; ++i)
    {
      REQUIRE(graph.halfEdge(he_id).face == static_cast<int>(face_id));
      he_id = static_cast<size_t>(graph.halfEdge(he_id).next);
    }
    REQUIRE(he_id == tri.half_edges[0]);
  }

  for (size_t he_id = 0; he_id < graph.halfEdgeSlotCount(); ++he_id)
  {
    if (graph.isLiveHalfEdge(he_id))
    {
      continue;
    }
    REQUIRE(HalfEdgeDelaunayGraph::isTombstone(graph.halfEdge(he_id)));
  }

  check_live_half_edge_face_references(graph);
}

TEST_CASE("HalfEdgeDelaunayGraph update preserves matched face ids", "[HalfEdgeDelaunayGraph]")
{
  std::vector<glm::dvec2> sites = {
    { 0.0, 0.0 },
    { 1.0, 0.0 },
    { 0.0, 1.0 },
    { 1.0, 1.0 },
    { 0.5, 0.5 },
  };

  HalfEdgeDelaunayGraph graph;
  graph.init(sites);

  std::vector<std::array<int, 3>> old_keys;
  old_keys.reserve(graph.liveFaceCount());
  std::vector<int> old_origins;
  old_origins.reserve(graph.halfEdgeSlotCount());
  for (size_t he_id = 0; he_id < graph.halfEdgeSlotCount(); ++he_id)
  {
    old_origins.push_back(graph.halfEdge(he_id).origin);
  }
  for (size_t ti : graph.liveFaces())
  {
    old_keys.push_back(face_vertex_key(ti, graph));
  }

  // Single-component update with identical sites should preserve every triangle id.
  graph.update(
    sites.size(),
    { { 0, 1, 2, 3, 4 } },
    [&sites](size_t v) { return sites[v]; });

  check_remap_invariants(graph);

  std::set<int> seen_face_ids;
  for (size_t ti : graph.liveFaces())
  {
    const int id = static_cast<int>(ti);
    REQUIRE(seen_face_ids.insert(id).second);

    const auto key = face_vertex_key(ti, graph);
    bool found_old = false;
    for (size_t old_ti = 0; old_ti < old_keys.size(); ++old_ti)
    {
      if (old_keys[old_ti] == key)
      {
        REQUIRE(old_ti == ti);
        found_old = true;
        break;
      }
    }
    REQUIRE(found_old);
  }

  for (size_t he_id = 0; he_id < old_origins.size(); ++he_id)
  {
    if (!graph.isLiveHalfEdge(he_id))
    {
      continue;
    }
    REQUIRE(graph.halfEdge(he_id).origin == old_origins[he_id]);
  }
}

TEST_CASE("HalfEdgeDelaunayGraph applyComponentSplit cuts cross-component topology", "[HalfEdgeDelaunayGraph]")
{
  std::vector<glm::dvec2> sites = {
    { 0.0, 0.0 },
    { 2.0, 0.0 },
    { 0.0, 2.0 },
    { 2.0, 2.0 },
    { 1.0, 1.0 },
    { 4.0, 1.0 },
    { 5.0, 0.0 },
    { 5.0, 2.0 },
  };

  HalfEdgeDelaunayGraph graph;
  graph.init(sites);

  std::vector<size_t> component_map(sites.size());
  for (size_t v : { 0, 1, 2, 4 })
  {
    component_map[v] = 0;
  }
  for (size_t v : { 3, 5, 6 })
  {
    component_map[v] = 1;
  }

  graph.applyRuntimeBranchSplit(component_map, [&sites](size_t v) { return sites[v]; });

  check_remap_invariants(graph);

  for (size_t he_id : graph.liveDelaunayEdges())
  {
    const int u = graph.halfEdge(he_id).origin;
    const int v = graph.halfEdge(he_id ^ 1).origin;
    if (u >= 0 && v >= 0)
    {
      REQUIRE(component_map[static_cast<size_t>(u)] == component_map[static_cast<size_t>(v)]);
    }
  }
}

TEST_CASE("HalfEdgeDelaunayGraph update assigns non-colliding ids for new triangles", "[HalfEdgeDelaunayGraph]")
{
  std::vector<glm::dvec2> sites = {
    { 0.0, 0.0 },
    { 2.0, 0.0 },
    { 0.0, 2.0 },
    { 2.0, 2.0 },
    { 1.0, 1.0 },
    { 4.0, 1.0 },
    { 5.0, 0.0 },
    { 5.0, 2.0 },
  };

  HalfEdgeDelaunayGraph graph;
  graph.init(sites);

  const size_t old_face_count = graph.faceSlotCount();

  graph.update(
    sites.size(),
    { { 0, 1, 2, 4 }, { 3, 5, 6 } },
    [&sites](size_t v) { return sites[v]; });

  check_remap_invariants(graph);

  REQUIRE(graph.faceSlotCount() >= old_face_count);
  REQUIRE(graph.liveFaceCount() > 0);
}

TEST_CASE("HalfEdgeDelaunayGraph update leaves tombstones for removed half-edges", "[HalfEdgeDelaunayGraph]")
{
  std::vector<glm::dvec2> sites = {
    { 0.0, 0.0 },
    { 2.0, 0.0 },
    { 0.0, 2.0 },
    { 2.0, 2.0 },
    { 1.0, 1.0 },
    { 4.0, 1.0 },
    { 5.0, 0.0 },
    { 5.0, 2.0 },
  };

  HalfEdgeDelaunayGraph graph;
  graph.init(sites);
  const size_t old_he_slots = graph.halfEdgeSlotCount();

  graph.update(
    sites.size(),
    { { 0, 1, 2, 4 }, { 3, 5, 6 } },
    [&sites](size_t v) { return sites[v]; });

  check_remap_invariants(graph);

  bool saw_tombstone = false;
  for (size_t he_id = 0; he_id < old_he_slots; ++he_id)
  {
    if (!graph.isLiveHalfEdge(he_id))
    {
      saw_tombstone = true;
      REQUIRE(HalfEdgeDelaunayGraph::isTombstone(graph.halfEdge(he_id)));
    }
  }
  REQUIRE(saw_tombstone);
}

TEST_CASE("HalfEdgeDelaunayGraph repeated update preserves half-edge topology", "[HalfEdgeDelaunayGraph]")
{
  std::vector<glm::dvec2> sites = {
    { 0.0, 0.0 },
    { 2.0, 0.0 },
    { 0.0, 2.0 },
    { 2.0, 2.0 },
    { 1.0, 1.0 },
    { 4.0, 1.0 },
    { 5.0, 0.0 },
    { 5.0, 2.0 },
  };

  HalfEdgeDelaunayGraph graph;
  graph.init(sites);

  graph.update(
    sites.size(),
    { { 0, 1, 2, 4 }, { 3, 5, 6 } },
    [&sites](size_t v) { return sites[v]; });
  check_remap_invariants(graph);

  graph.update(
    sites.size(),
    { { 0, 1, 2, 4 }, { 3, 5, 6 } },
    [&sites](size_t v) { return sites[v]; });
  check_remap_invariants(graph);
}

static std::vector<size_t> collect_infinite_faces(const HalfEdgeDelaunayGraph& graph)
{
  std::vector<size_t> infinite_faces;
  for (size_t face_id : graph.liveFaces())
  {
    const auto verts = graph.getTriangleVertexIndices(face_id);
    if (verts[0] == -1 || verts[1] == -1 || verts[2] == -1)
    {
      infinite_faces.push_back(face_id);
    }
  }
  return infinite_faces;
}

static std::pair<int, int> infinite_face_hull_edge_key(const HalfEdgeDelaunayGraph& graph, size_t face_id)
{
  const auto verts = graph.adjacentTriangleVertices(graph.face(face_id).half_edges[0]);
  int finite[2] = { -1, -1 };
  int count = 0;
  for (int v : verts)
  {
    if (v >= 0)
    {
      finite[count++] = v;
    }
  }
  if (finite[0] > finite[1])
  {
    std::swap(finite[0], finite[1]);
  }
  return { finite[0], finite[1] };
}

static std::map<std::pair<int, int>, glm::dvec2> infinite_directions_by_hull_edge(
  const HalfEdgeDelaunayGraph& graph, const std::vector<glm::dvec2>& sites)
{
  const auto circumcenters = graph.computeCircumcenters(sites);
  std::map<std::pair<int, int>, glm::dvec2> directions;
  for (size_t face_id : collect_infinite_faces(graph))
  {
    directions.emplace(infinite_face_hull_edge_key(graph, face_id), circumcenters[face_id].first);
  }
  return directions;
}

static void require_infinite_directions_point_outward(
  const HalfEdgeDelaunayGraph& graph, const std::vector<glm::dvec2>& sites)
{
  const auto circumcenters = graph.computeCircumcenters(sites);
  for (size_t face_id : collect_infinite_faces(graph))
  {
    const glm::dvec2 dir = circumcenters[face_id].first;
    REQUIRE(circumcenters[face_id].second);

    bool checked_face = false;
    for (size_t he_id : graph.face(face_id).half_edges)
    {
      if (!graph.isOnConvexBoundaryOutside(he_id))
      {
        continue;
      }

      const int origin = graph.halfEdge(he_id).origin;
      const int dest = graph.destination(he_id);
      REQUIRE(origin >= 0);
      REQUIRE(dest >= 0);

      const int interior_vertex = graph.triangleOppositeVertex(he_id ^ 1);
      REQUIRE(interior_vertex >= 0);

      const glm::dvec2 u = sites[static_cast<size_t>(origin)];
      const glm::dvec2 v = sites[static_cast<size_t>(dest)];
      const glm::dvec2 w = sites[static_cast<size_t>(interior_vertex)];
      const glm::dvec2 midpoint = 0.5 * (u + v);
      REQUIRE(glm::dot(dir, midpoint - w) > 0.0);
      checked_face = true;
      break;
    }

    REQUIRE(checked_face);
  }
}

static void require_same_infinite_directions(
  const std::map<std::pair<int, int>, glm::dvec2>& reference,
  const HalfEdgeDelaunayGraph& graph,
  const std::vector<glm::dvec2>& sites,
  bool require_all_reference_keys = true)
{
  const auto current = infinite_directions_by_hull_edge(graph, sites);
  size_t matched = 0;
  for (const auto& [key, before] : reference)
  {
    const auto after_it = current.find(key);
    if (after_it == current.end())
    {
      continue;
    }
    ++matched;
    const glm::dvec2& after = after_it->second;
    const double dot = before.x * after.x + before.y * after.y;
    const double before_len2 = before.x * before.x + before.y * before.y;
    const double after_len2 = after.x * after.x + after.y * after.y;
    REQUIRE(before_len2 > 0.0);
    REQUIRE(after_len2 > 0.0);
    REQUIRE(dot / std::sqrt(before_len2 * after_len2) > 0.999);
  }
  if (require_all_reference_keys)
  {
    REQUIRE(matched == reference.size());
  }
  else
  {
    REQUIRE(matched > 0);
  }
}

TEST_CASE("HalfEdgeDelaunayGraph update restores convex hull infinite face refs after split", "[HalfEdgeDelaunayGraph]")
{
  std::vector<glm::dvec2> sites = {
    { 0.0, 0.0 },
    { 2.0, 0.0 },
    { 0.0, 2.0 },
    { 2.0, 2.0 },
    { 1.0, 1.0 },
    { 5.0, 1.0 },
    { 6.0, 0.0 },
    { 6.0, 2.0 },
  };

  HalfEdgeDelaunayGraph graph;
  graph.init(sites);

  graph.update(
    sites.size(),
    { { 0, 1, 2, 4 }, { 3, 5, 6 } },
    [&sites](size_t v) { return sites[v]; });

  check_remap_invariants(graph);

  size_t convex_hull_edges = 0;
  for (size_t he_id : graph.liveDelaunayEdges())
  {
    if (!graph.isOnConvexBoundary(he_id))
    {
      continue;
    }
    const auto& he = graph.halfEdge(he_id);
    const auto& twin = graph.halfEdge(he_id ^ 1);
    if (he.origin == -1 || twin.origin == -1)
    {
      continue;
    }
    ++convex_hull_edges;
    const bool inf0 = face_has_infinite_vertex(graph, static_cast<size_t>(he.face));
    const bool inf1 = face_has_infinite_vertex(graph, static_cast<size_t>(twin.face));
    REQUIRE(inf0 != inf1);
  }
  REQUIRE(convex_hull_edges > 0);
}

TEST_CASE("HalfEdgeDelaunayGraph infinite Voronoi ray direction is slot-order invariant", "[HalfEdgeDelaunayGraph]")
{
  std::vector<glm::dvec2> sites = {
    { 0.0, 0.0 },
    { 2.0, 0.0 },
    { 0.0, 2.0 },
    { 2.0, 2.0 },
    { 1.0, 1.0 },
  };

  HalfEdgeDelaunayGraph graph;
  graph.init(sites);
  const auto reference_directions = infinite_directions_by_hull_edge(graph, sites);
  REQUIRE(!reference_directions.empty());
  require_infinite_directions_point_outward(graph, sites);

  // flipEdge rewrites triangle.half_edges slots; reorder_from_old can do the same.
  bool flipped_hull_edge = false;
  for (size_t he_id : graph.liveDelaunayEdges())
  {
    if (!graph.isOnConvexBoundary(he_id))
    {
      continue;
    }
    const auto& he = graph.halfEdge(he_id);
    const auto& twin = graph.halfEdge(he_id ^ 1);
    if (he.origin == -1 || twin.origin == -1)
    {
      continue;
    }
    graph.flipEdge(he_id);
    flipped_hull_edge = true;
    break;
  }
  REQUIRE(flipped_hull_edge);
  require_same_infinite_directions(reference_directions, graph, sites, false);
  require_infinite_directions_point_outward(graph, sites);

  graph.update(
    sites.size(),
    { { 0, 1, 2, 3, 4 } },
    [&sites](size_t v) { return sites[v]; });
  check_remap_invariants(graph);
  require_same_infinite_directions(reference_directions, graph, sites, true);
  require_infinite_directions_point_outward(graph, sites);
}

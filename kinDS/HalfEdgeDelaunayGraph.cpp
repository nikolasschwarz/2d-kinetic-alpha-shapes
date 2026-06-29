#include "HalfEdgeDelaunayGraph.hpp"

#include "Logger.hpp"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>
#include <vector>

using namespace kinDS;

namespace
{
void tombstoneHalfEdge(std::vector<HalfEdgeDelaunayGraph::HalfEdge>& half_edges, size_t he_id)
{
  half_edges[he_id] = HalfEdgeDelaunayGraph::HalfEdge {};
}

void killFaceSlot(HalfEdgeDelaunayGraph::Triangle& tri, size_t invalid_he)
{
  tri.half_edges[0] = invalid_he;
  tri.half_edges[1] = invalid_he;
  tri.half_edges[2] = invalid_he;
}

double ccwAngle(const glm::dvec2& from, const glm::dvec2& to)
{
  return std::atan2(to.y - from.y, to.x - from.x);
}

size_t nextBoundaryHalfEdgeLeaving(size_t vertex, size_t arrived_from, const std::vector<size_t>& boundary_edges,
  const std::vector<HalfEdgeDelaunayGraph::HalfEdge>& half_edges,
  const std::function<glm::dvec2(size_t)>& vertex_at, size_t component_id, const std::vector<size_t>& component_map)
{
  const glm::dvec2 pos_v = vertex_at(vertex);
  const double arrive_angle = ccwAngle(pos_v, vertex_at(arrived_from));
  double best_delta = std::numeric_limits<double>::infinity();
  size_t best_he = static_cast<size_t>(-1);

  for (size_t he_id : boundary_edges)
  {
    const int origin = half_edges[he_id].origin;
    const int dest = half_edges[he_id ^ 1].origin;
    if (origin < 0 || dest < 0 || static_cast<size_t>(origin) != vertex)
    {
      continue;
    }
    if (component_map[static_cast<size_t>(origin)] != component_id)
    {
      continue;
    }

    const double out_angle = ccwAngle(pos_v, vertex_at(static_cast<size_t>(dest)));
    double delta = out_angle - arrive_angle;
    while (delta <= 0.0)
    {
      delta += 2.0 * std::acos(-1.0);
    }
    if (delta < best_delta)
    {
      best_delta = delta;
      best_he = he_id;
    }
  }

  return best_he;
}

std::vector<size_t> traceBoundaryCycle(size_t start_he, const std::vector<size_t>& boundary_edges,
  const std::vector<HalfEdgeDelaunayGraph::HalfEdge>& half_edges,
  const std::function<glm::dvec2(size_t)>& vertex_at, size_t component_id, const std::vector<size_t>& component_map,
  std::unordered_set<size_t>& visited)
{
  std::vector<size_t> cycle;
  cycle.reserve(boundary_edges.size());
  size_t current = start_he;
  const size_t start_tail = static_cast<size_t>(half_edges[start_he ^ 1].origin);
  for (size_t guard = 0; guard <= boundary_edges.size(); ++guard)
  {
    if (!visited.insert(current).second)
    {
      break;
    }
    cycle.push_back(current);
    const size_t tail = static_cast<size_t>(half_edges[current ^ 1].origin);
    const size_t tip = static_cast<size_t>(half_edges[current].origin);
    current = nextBoundaryHalfEdgeLeaving(tip, tail, boundary_edges, half_edges, vertex_at, component_id, component_map);
    if (current == static_cast<size_t>(-1))
    {
      break;
    }
    if (current == start_he && static_cast<size_t>(half_edges[current ^ 1].origin) == start_tail)
    {
      break;
    }
  }

  return cycle;
}

void createInfiniteFacesFromBoundary(const std::vector<size_t>& boundary_interior,
  std::vector<HalfEdgeDelaunayGraph::HalfEdge>& half_edges, std::vector<HalfEdgeDelaunayGraph::Triangle>& triangles,
  size_t vertex_count)
{
  if (boundary_interior.empty())
  {
    return;
  }

  std::vector<int> boundary_edge_map(vertex_count, -1);
  for (size_t interior_he : boundary_interior)
  {
    const size_t exterior_he = interior_he ^ 1;
    HalfEdgeDelaunayGraph::HalfEdge& exterior = half_edges[exterior_he];
    const int u = exterior.origin;
    if (u < 0 || static_cast<size_t>(u) >= vertex_count)
    {
      throw std::runtime_error("createInfiniteFacesFromBoundary: boundary exterior has invalid finite origin");
    }
    exterior.face = -1;
    exterior.next = -1;
    if (boundary_edge_map[static_cast<size_t>(u)] != -1)
    {
      throw std::runtime_error("createInfiniteFacesFromBoundary: duplicate outgoing boundary edge at vertex "
        + std::to_string(static_cast<size_t>(u)));
    }
    boundary_edge_map[static_cast<size_t>(u)] = static_cast<int>(exterior_he);
  }

  std::vector<int> incoming_edge_map(vertex_count, -1);
  for (size_t u = 0; u < vertex_count; ++u)
  {
    const int boundary_edge_index = boundary_edge_map[u];
    if (boundary_edge_index == -1)
    {
      continue;
    }

    HalfEdgeDelaunayGraph::HalfEdge& boundary = half_edges[static_cast<size_t>(boundary_edge_index)];
    const int v = half_edges[static_cast<size_t>(boundary_edge_index ^ 1)].origin;
    if (v < 0 || static_cast<size_t>(v) >= vertex_count)
    {
      throw std::runtime_error("createInfiniteFacesFromBoundary: boundary edge has invalid finite destination");
    }
    const int next_boundary_he = boundary_edge_map[static_cast<size_t>(v)];
    if (next_boundary_he == -1)
    {
      throw std::runtime_error("createInfiniteFacesFromBoundary: missing successor boundary edge at vertex "
        + std::to_string(static_cast<size_t>(v)));
    }

    const size_t he_at_v = half_edges.size();
    half_edges.push_back(HalfEdgeDelaunayGraph::HalfEdge {});
    const size_t he_at_inf = half_edges.size();
    half_edges.push_back(HalfEdgeDelaunayGraph::HalfEdge {});

    half_edges[he_at_v].origin = v;
    half_edges[he_at_v].face = -1;
    half_edges[he_at_inf].origin = -1;
    half_edges[he_at_inf].face = -1;
    half_edges[he_at_inf].next = next_boundary_he;

    boundary.next = static_cast<int>(he_at_v);
    incoming_edge_map[static_cast<size_t>(v)] = static_cast<int>(he_at_inf);
  }

  for (size_t u = 0; u < vertex_count; ++u)
  {
    const int incoming_edge_index = incoming_edge_map[u];
    if (incoming_edge_index == -1)
    {
      continue;
    }

    HalfEdgeDelaunayGraph::HalfEdge& incoming = half_edges[static_cast<size_t>(incoming_edge_index)];
    HalfEdgeDelaunayGraph::HalfEdge& boundary = half_edges[static_cast<size_t>(incoming.next)];
    HalfEdgeDelaunayGraph::HalfEdge& outgoing = half_edges[static_cast<size_t>(boundary.next)];

    outgoing.next = incoming_edge_index;

    const size_t new_face = triangles.size();
    triangles.emplace_back();
    triangles[new_face].half_edges[0] = static_cast<size_t>(incoming_edge_index);
    triangles[new_face].half_edges[1] = static_cast<size_t>(incoming.next);
    triangles[new_face].half_edges[2] = static_cast<size_t>(boundary.next);

    incoming.face = static_cast<int>(new_face);
    boundary.face = static_cast<int>(new_face);
    outgoing.face = static_cast<int>(new_face);
  }

  for (size_t interior_he : boundary_interior)
  {
    const size_t exterior_he = interior_he ^ 1;
    if (half_edges[exterior_he].next < 0 || half_edges[exterior_he].face < 0)
    {
      throw std::runtime_error(
        "createInfiniteFacesFromBoundary: failed to assign rebuilt infinite face to boundary edge "
        + std::to_string(interior_he));
    }
  }
}
} // namespace

void HalfEdgeDelaunayGraph::build(const std::vector<size_t>& index_buffer)
{
  KINDS_DEBUG("Building half-edge mesh from triangle index buffer of size " << index_buffer.size());
  assert(index_buffer.size() % 3 == 0 && "Input must be a triangle index buffer.");
  const int num_tris = index_buffer.size() / 3;

  // TODO: update to account for components
  // note that the following reserves are off by one compared to Euler's formula because there is an implicit vertex at
  // infinity that is not counted in the vertex count
  half_edges.clear();
  half_edges.reserve(
    6 * (vertex_count - 1)); // Reserve space for half-edges, at most 6 * (V - 1) half-edges in a triangulation

  triangles.clear();
  triangles.reserve(2 * (vertex_count - 1)); // Reserve space for faces, at most 2 * (V - 1) faces in a triangulation
  triangles.resize(num_tris);

  // Store outgoing edge along boundary from vertex, -1 for non-boundary vertices
  std::vector<int> boundary_edge_map(vertex_count, -1);

  struct TmpHalfEdge
  {
    size_t v;
    size_t f;
    size_t j; // index of the half-edge in the triangle
  };

  // Sort up- and down-edges into an adjacency list by the higher vertex index
  std::vector<std::vector<TmpHalfEdge>> adjacency_list_up_in(vertex_count);
  std::vector<std::vector<TmpHalfEdge>> adjacency_list_down_out(vertex_count);

  // iterate through the triangles
  for (size_t i = 0; i < num_tris; ++i)
  {
    size_t v[3];
    v[0] = index_buffer[i * 3];
    v[1] = index_buffer[i * 3 + 1];
    v[2] = index_buffer[i * 3 + 2];

    // iterate through the triangle's edges
    for (size_t j = 0; j < 3; ++j)
    {
      size_t a = v[j];
      size_t b = v[(j + 1) % 3];

      if (a > b)
      {
        adjacency_list_down_out[a].push_back({ b, i, j });
      }
      else
      {
        adjacency_list_up_in[b].push_back({ a, i, j });
      }
    }
  }

  // now sort the opposite way to bring them into order
  std::vector<std::vector<TmpHalfEdge>> adjacency_list_up_out(vertex_count);
  std::vector<std::vector<TmpHalfEdge>> adjacency_list_down_in(vertex_count);

  for (size_t u = 0; u < vertex_count; u++)
  {
    for (TmpHalfEdge& the : adjacency_list_up_in[u])
    {
      adjacency_list_up_out[the.v].push_back({ u, the.f, the.j });
    }

    for (TmpHalfEdge& the : adjacency_list_down_out[u])
    {
      adjacency_list_down_in[the.v].push_back({ u, the.f, the.j });
    }
  }

  // finally, iterate over both adjacency lists to build the half-edges and faces
  for (size_t u = 0; u < vertex_count; u++)
  {
    // go through both adjacency lists in a merge-like way
    auto up_it = adjacency_list_up_out[u].begin();
    auto down_it = adjacency_list_down_in[u].begin();

    while (up_it != adjacency_list_up_out[u].end() || down_it != adjacency_list_down_in[u].end())
    {
      int v_up = (up_it != adjacency_list_up_out[u].end()) ? up_it->v : std::numeric_limits<int>::max();
      int v_down = (down_it != adjacency_list_down_in[u].end()) ? down_it->v : std::numeric_limits<int>::max();

      // always create both half-edges, but if one of them is not present, it will be a boundary half-edge
      HalfEdge he_up;
      he_up.origin = u;
      he_up.face = (v_up <= v_down) ? up_it->f : -1; // -1 means boundary
      he_up.next = -1; // will be set later

      if (he_up.face != -1)
      {
        triangles[he_up.face].half_edges[up_it->j] = half_edges.size(); // store the half-edge index in the face
      }
      else
      {
        boundary_edge_map[u] = half_edges.size(); // store the boundary edge for this vertex
      }
      half_edges.push_back(he_up);

      size_t v = std::min(v_up, v_down);

      HalfEdge he_down;
      he_down.origin = v;
      he_down.face = (v_down <= v_up) ? down_it->f : -1; // -1 means boundary
      he_down.next = -1; // will be set later

      if (he_down.face != -1)
      {
        triangles[he_down.face].half_edges[down_it->j] = half_edges.size(); // store the half-edge index in the face
      }
      else
      {
        boundary_edge_map[v] = half_edges.size(); // store the boundary edge for this vertex
      }
      half_edges.push_back(he_down);

      // increment the iterators
      if (v_up <= v_down)
      {
        up_it++;
      }
      if (v_down <= v_up)
      {
        down_it++;
      }
    }
  }

  // now link the half-edges together
  // triangle edges first
  for (auto& face : triangles)
  {
    for (size_t j = 0; j < 3; ++j)
    {
      size_t he_index = face.half_edges[j];
      if (he_index != -1)
      {
        HalfEdge& he = half_edges[he_index];
        he.next = face.half_edges[(j + 1) % 3]; // link to the next half-edge in the face
      }
      else
      {
        KINDS_DEBUG("Face " << (&face - &triangles[0]) << " has a missing half-edge at index " << j << ".");
      }
    }
  }

  // Store incoming edge from infinity to vertex, -1 for non-boundary vertices
  std::vector<int> incoming_edge_map(vertex_count, -1);

  // at the boundary, create additional faces and half-edges that connect to a vertex at infinity
  for (size_t u = 0; u < vertex_count; ++u)
  {
    int boundary_edge_index = boundary_edge_map[u];
    if (boundary_edge_index != -1)
    {
      HalfEdge& he = half_edges[boundary_edge_index];
      size_t v = destination(boundary_edge_index); // situated between he and next_he
      size_t next_he_id = boundary_edge_map[v]; // next boundary half-edge
      HalfEdge& next_he = half_edges[next_he_id];
      // create new pair of half-edges that connect to a vertex at infinity
      HalfEdge he_infinity;
      he_infinity.origin = v; // origin is the current vertex
      he_infinity.face = -1; // -1 means boundary, assign proper id later
      he.next = half_edges.size(); // link the current half-edge to the new one

      HalfEdge he_infinity_twin;
      he_infinity_twin.origin = -1; // vertex at infinity
      he_infinity_twin.face = -1; // -1 means boundary, assign proper id later
      he_infinity_twin.next = next_he_id;

      half_edges.push_back(he_infinity);
      // store the incoming edge from infinity so we can later find it
      incoming_edge_map[v] = half_edges.size();
      half_edges.push_back(he_infinity_twin);
    }
  }

  // need another run to create faces and link the half-edges that meet at infinity
  for (size_t u = 0; u < vertex_count; ++u)
  {
    int incoming_edge_index = incoming_edge_map[u];
    if (incoming_edge_index == -1)
    {
      continue; // not a boundary vertex
    }
    HalfEdge& incoming = half_edges[incoming_edge_index];
    HalfEdge& boundary = half_edges[incoming.next];
    assert(incoming.next == boundary_edge_map[u]);
    HalfEdge& outgoing = half_edges[boundary.next];

    outgoing.next = incoming_edge_index; // link the outgoing half-edge to the incoming one

    // create new face for the boundary edge
    Triangle new_face;
    new_face.half_edges[0] = incoming_edge_index; // first half-edge is the incoming one
    new_face.half_edges[1] = incoming.next; // second half-edge is the boundary half-edge
    new_face.half_edges[2] = boundary.next; // third half-edge is the outgoing one

    incoming.face = triangles.size(); // assign the new face index to the incoming half-edge
    boundary.face = triangles.size(); // assign the new face index to the boundary half-edge
    outgoing.face = triangles.size(); // assign the new face index to the outgoing half-edge

    triangles.push_back(new_face); // add the new face to the list
  }

  // iterate over all edges to set vertex to half-edge mapping
  for (size_t he_id = 0; he_id < half_edges.size(); ++he_id)
  {
    const HalfEdge& he = half_edges[he_id];
    if (he.origin != -1)
    {
      vertex_to_half_edge[he.origin] = he_id;
    }
  }

  KINDS_DEBUG("Half-edge mesh built with " << half_edges.size() << " half-edges and " << triangles.size() << " faces.");
  rebuildLiveIndices();
}

void HalfEdgeDelaunayGraph::rebuildLiveIndices()
{
  half_edge_live_.assign(half_edges.size(), 0);
  live_even_half_edges_.clear();
  live_even_half_edges_.reserve(half_edges.size() / 2);
  for (size_t he_id = 0; he_id < half_edges.size(); he_id += 2)
  {
    if (!isLiveHalfEdge(he_id))
    {
      continue;
    }
    half_edge_live_[he_id] = 1;
    half_edge_live_[he_id ^ 1] = 1;
    live_even_half_edges_.push_back(he_id);
  }

  live_faces_.clear();
  live_faces_.reserve(triangles.size());
  for (size_t face_id = 0; face_id < triangles.size(); ++face_id)
  {
    if (isLiveFace(face_id))
    {
      live_faces_.push_back(face_id);
    }
  }
}

bool HalfEdgeDelaunayGraph::isLiveHalfEdge(size_t he_id) const
{
  if (he_id >= half_edges.size())
  {
    return false;
  }
  const HalfEdge& he = half_edges[he_id];
  if (isTombstone(he))
  {
    return false;
  }
  // Every connected half-edge in a built mesh has its face cycle linked.
  return he.next >= 0;
}

bool HalfEdgeDelaunayGraph::isLiveFace(size_t face_id) const
{
  if (face_id >= triangles.size())
  {
    return false;
  }

  const Triangle& tri = triangles[face_id];
  if (tri.half_edges[0] >= half_edges.size())
  {
    return false;
  }

  for (int i = 0; i < 3; ++i)
  {
    const size_t he_id = tri.half_edges[i];
    if (he_id >= half_edges.size())
    {
      return false;
    }
    if (half_edges[he_id].face != static_cast<int>(face_id))
    {
      return false;
    }
  }
  return true;
}

bool HalfEdgeDelaunayGraph::halfEdgeInFaceCycle(size_t he_id, size_t face_id) const
{
  if (!isLiveFace(face_id))
  {
    return false;
  }

  const size_t start = triangles[face_id].half_edges[0];
  size_t cur = start;
  for (int i = 0; i < 3; ++i)
  {
    if (cur == he_id)
    {
      return true;
    }
    cur = static_cast<size_t>(half_edges[cur].next);
  }
  return false;
}

bool HalfEdgeDelaunayGraph::faceHasInfiniteVertex(size_t face_id) const
{
  const auto verts = getTriangleVertexIndices(face_id);
  return verts[0] == -1 || verts[1] == -1 || verts[2] == -1;
}

void HalfEdgeDelaunayGraph::validateLiveHalfEdgeFaceReferences() const
{
  for (size_t he_id = 0; he_id < half_edges.size(); ++he_id)
  {
    if (!isLiveHalfEdge(he_id))
    {
      continue;
    }

    const HalfEdge& he = half_edges[he_id];
    const int face_id = he.face;
    if (face_id < 0 || static_cast<size_t>(face_id) >= triangles.size() || !isLiveFace(static_cast<size_t>(face_id)))
    {
      std::ostringstream message;
      message << "Live half-edge " << he_id << " (origin=" << he.origin << ", destination=" << destination(he_id)
              << ") references ";
      if (face_id < 0)
      {
        message << "invalid face index " << face_id;
      }
      else
      {
        message << "dead face slot " << face_id;
      }
      throw std::runtime_error(message.str());
    }
  }
}

void HalfEdgeDelaunayGraph::applyFreshHalfEdgeFaceReferences(
  const std::vector<HalfEdge>& fresh_half_edges, const std::vector<int>& tri_remap)
{
  auto remap_face = [&](int fresh_face) -> int
  {
    if (fresh_face < 0)
    {
      return -1;
    }
    const size_t f = static_cast<size_t>(fresh_face);
    if (f >= tri_remap.size())
    {
      return -1;
    }
    return tri_remap[f];
  };

  auto directed_key = [](int origin, int dest) -> uint64_t
  {
    return (static_cast<uint64_t>(static_cast<uint32_t>(origin)) << 32)
      | static_cast<uint64_t>(static_cast<uint32_t>(dest));
  };

  // Keyed by directed site edge; authoritative face incidence from per-branch retriangulation.
  std::unordered_map<uint64_t, int> fresh_face_by_directed;
  fresh_face_by_directed.reserve(fresh_half_edges.size());
  for (size_t fresh_hi = 0; fresh_hi < fresh_half_edges.size(); ++fresh_hi)
  {
    const int origin = fresh_half_edges[fresh_hi].origin;
    const int dest = fresh_half_edges[fresh_hi ^ 1].origin;
    if (origin < 0 || dest < 0)
    {
      continue;
    }
    fresh_face_by_directed[directed_key(origin, dest)] = remap_face(fresh_half_edges[fresh_hi].face);
  }

  for (size_t he_id = 0; he_id < half_edges.size(); ++he_id)
  {
    if (!isLiveHalfEdge(he_id))
    {
      continue;
    }
    const int origin = half_edges[he_id].origin;
    const int dest = half_edges[he_id ^ 1].origin;
    if (origin < 0 || dest < 0)
    {
      continue;
    }
    const auto it = fresh_face_by_directed.find(directed_key(origin, dest));
    if (it != fresh_face_by_directed.end())
    {
      half_edges[he_id].face = it->second;
    }
  }
}

const HalfEdgeDelaunayGraph::HalfEdge& HalfEdgeDelaunayGraph::halfEdge(size_t he_id) const
{
  assert(he_id < half_edges.size());
  return half_edges[he_id];
}

const HalfEdgeDelaunayGraph::Triangle& HalfEdgeDelaunayGraph::face(size_t face_id) const
{
  assert(face_id < triangles.size());
  return triangles[face_id];
}

size_t HalfEdgeDelaunayGraph::halfEdgeSlotCount() const { return half_edges.size(); }

size_t HalfEdgeDelaunayGraph::faceSlotCount() const { return triangles.size(); }

void kinDS::HalfEdgeDelaunayGraph::flipEdge(size_t he_id)
{
  if (he_id >= half_edges.size())
  {
    KINDS_ERROR("Invalid half-edge ID for flipping: " << he_id);
    return;
  }
  HalfEdge& he = half_edges[he_id];
  HalfEdge& twin = half_edges[he_id ^ 1];

  size_t u = he.origin; // origin vertex of the half-edge
  size_t v = twin.origin; // origin vertex of the twin half-edge

  // Check if we can flip the edge
  if (he.face == -1 || twin.face == -1)
  {
    KINDS_ERROR("Cannot flip edge " << he_id << " because one of the triangles is invalid.");
    return;
  }

  // check if edge is referenced in vertex_to_half_edge and update
  if (u != size_t(-1))
  {

    if (vertex_to_half_edge[u] == he_id)
    {
      // just set to next half-edge on the vertex
      vertex_to_half_edge[u] = neighborEdgeId(he_id);
    }
  }

  if (v != size_t(-1))
  {
    if (vertex_to_half_edge[v] == HalfEdgeDelaunayGraph::twin(he_id))
    {
      vertex_to_half_edge[v] = neighborEdgeId(HalfEdgeDelaunayGraph::twin(he_id)); // update to point to the half-edge
    }
  }

  int he_next_id = he.next;
  int twin_next_id = twin.next;

  size_t he_last_id = half_edges[he_next_id].next; // The next half-edge in the triangle
  size_t twin_last_id = half_edges[twin_next_id].next; // The next half-edge in the twin triangle

  size_t he_opposite_id = triangleOppositeVertex(he_id);
  size_t twin_opposite_id = triangleOppositeVertex(he_id ^ 1);

  // Update the half-edges to flip the edge
  he.origin = twin_opposite_id;
  twin.origin = he_opposite_id;

  // Update the next pointers
  // bridge the two triangles where the edge used to be
  half_edges[he_last_id].next = twin_next_id;
  half_edges[twin_last_id].next = he_next_id;

  // intercept at the new end points of the flipped edge
  half_edges[he_next_id].next = he_id ^ 1;
  half_edges[twin_next_id].next = he_id;

  // finally, update the twin pointers
  half_edges[he_id].next = he_last_id; // the half-edge now points to the twin
  half_edges[he_id ^ 1].next = twin_last_id; // the twin of he now points to the next half-edge

  // Update the faces
  std::swap(half_edges[he_next_id].face, half_edges[twin_next_id].face);

  // Update the half-edge indices in the faces
  triangles[he.face].half_edges[0] = he_id; // Update the first half-edge of the face
  triangles[he.face].half_edges[1] = he_last_id; // Update the second half-edge of the face
  triangles[he.face].half_edges[2] = twin_next_id; // Update the third half-edge of the face

  triangles[twin.face].half_edges[0] = he_id ^ 1; // Update the first half-edge of the twin face
  triangles[twin.face].half_edges[1] = twin_last_id; // Update the second half-edge of the twin face
  triangles[twin.face].half_edges[2] = he_next_id; // Update the third half-edge of the twin face

  // KINDS_DEBUG("Flipped edge " << he_id << " between vertices " << u << " and " << v << ".");

  // printDebug();
}

void HalfEdgeDelaunayGraph::printDebug(std::ostream* out) const
{
  std::ostream& os = out != nullptr ? *out : std::cout;

  os << "Half-Edges:\n";
  for (size_t i = 0; i < half_edges.size(); ++i)
  {
    const HalfEdge& he = half_edges[i];
    os << "  [" << i << "]" << (isLiveHalfEdge(i) ? "" : " \u2020") << " origin = " << he.origin
       << ", next = " << he.next << ", face = " << he.face << ", twin = " << (i ^ 1) << "\n";
  }

  os << "\nFaces:\n";
  for (size_t i = 0; i < triangles.size(); ++i)
  {
    const Triangle& f = triangles[i];
    os << "  [" << i << "]" << (isLiveFace(i) ? "" : " \u2020") << " half_edge = " << f.half_edges[0] << "\n";

    // Walk the face's boundary
    int start = f.half_edges[0];
    if (start < 0)
      continue;
    int he = start;
    os << "    Vertices: ";
    do
    {
      os << half_edges[he].origin << " ";
      he = half_edges[he].next;
    } while (he != start && he != -1);
    os << "\n";

    he = start;
    os << "    Half-Edges: ";
    do
    {
      os << he << " ";
      he = half_edges[he].next;
    } while (he != start && he != -1);
    os << "\n";
  }

  os << "Vertex incident half-edges:\n";
  for (size_t u = 0; u < vertex_count; ++u)
  {
    os << " Outgoing edges from vertex " << u << ":\n";

    for (IncidentEdgeIterator it = incidentEdgesBegin(u); it != incidentEdgesEnd(u); ++it)
    {
      size_t he_id = *it;
      const HalfEdge& he = half_edges[he_id];
      os << "  Half-edge " << he_id << " to vertex " << destination(he_id) << " in face " << he.face << "\n";
    }
  }
}

void HalfEdgeDelaunayGraph::init(const std::vector<glm::dvec2>& site_positions)
{
  vertex_count = site_positions.size();
  vertex_to_half_edge.assign(vertex_count, -1);
  std::vector<float> coords;
  coords.reserve(site_positions.size() * 2);
  for (const auto& point : site_positions)
  {
    coords.push_back(point[0]);
    coords.push_back(point[1]);
  }

  Delaunator2D delaunator(coords);

  build(delaunator.triangles);
}

void HalfEdgeDelaunayGraph::init(const std::vector<std::vector<glm::dvec2>>& splines)
{
  std::vector<glm::dvec2> site_positions;
  site_positions.reserve(splines.size());
  for (const auto& spline : splines)
  {
    site_positions.push_back(spline.front());
  }
  init(site_positions);
}

void HalfEdgeDelaunayGraph::applyComponentSplit(const std::vector<size_t>& component_map,
  const std::function<glm::dvec2(size_t)>& vertex_at, std::optional<double> debug_time)
{
  // In-place component split (no retriangulation). Given component_map[v] = component id for each site vertex,
  // disconnect topology across components while keeping each component's interior triangulation:
  //  1) decide which original Delaunay edge pairs survive the split
  //  2) kill any finite triangle that uses a dead cross-component edge
  //  3) reuse the surviving boundary edge pairs whose exterior face died to create new infinite caps
  //  4) tombstone edge pairs with no live incident face on either side; refresh indices

  if (component_map.size() < vertex_count)
  {
    throw std::runtime_error("applyComponentSplit: component_map size mismatch");
  }

  static size_t split_debug_dump_counter = 0;
  const size_t debug_dump_id = split_debug_dump_counter++;
  std::string debug_tag = std::to_string(debug_dump_id);
  if (debug_time.has_value())
  {
    std::ostringstream oss;
    oss << "_t" << *debug_time;
    std::string time_suffix = oss.str();
    for (char& ch : time_suffix)
    {
      if (!std::isalnum(static_cast<unsigned char>(ch)))
      {
        ch = '_';
      }
    }
    debug_tag += time_suffix;
  }
  {
    std::ofstream before_dump("applyComponentSplit_before_" + debug_tag + ".txt");
    if (before_dump.is_open())
    {
      printDebug(&before_dump);
    }
  }

  // Sentinel half-edge pair appended now so killFaceSlot() can mark dead triangles without colliding with live slots.
  // invalid_he is the even index of that reserved pair; only used as a tombstone marker in Triangle::half_edges.
  const size_t invalid_he = half_edges.size();
  half_edges.emplace_back();
  half_edges.emplace_back();

  const size_t original_half_edge_count = invalid_he;
  std::vector<bool> kill_he(original_half_edge_count, false);

  // --- Pass 1: mark Delaunay edges whose endpoints lie in different components ---
  for (size_t he_id = 0; he_id < original_half_edge_count; he_id += 2)
  {
    if (!isLiveHalfEdge(he_id))
    {
      continue;
    }
    const int u = half_edges[he_id].origin;
    const int v = half_edges[he_id ^ 1].origin;
    if (u >= 0 && v >= 0 && component_map[static_cast<size_t>(u)] != component_map[static_cast<size_t>(v)])
    {
      kill_he[he_id] = true;
      kill_he[he_id ^ 1] = true;
    }
  }

  // --- Pass 2: kill any finite triangle touching a cross-component edge ---
  std::vector<bool> kill_face(triangles.size(), false);
  for (size_t face_id = 0; face_id < triangles.size(); ++face_id)
  {
    if (!isLiveFace(face_id))
    {
      continue;
    }

    if (faceHasInfiniteVertex(face_id))
    {
      continue;
    }

    const auto face_hes = getTriangleHalfEdgeIndices(triangles[face_id].half_edges[0]);
    if (kill_he[face_hes[0]] || kill_he[face_hes[1]] || kill_he[face_hes[2]])
    {
      kill_face[face_id] = true;
    }
  }

  // Kill mixed faces: triangle slots remain allocated but are marked dead (all three half-edge refs -> invalid_he).
  for (size_t face_id = 0; face_id < triangles.size(); ++face_id)
  {
    if (kill_face[face_id])
    {
      killFaceSlot(triangles[face_id], invalid_he);
    }
  }

  // The old infinity topology is no longer authoritative after the cut. Remove all original infinite faces and
  // half-edge pairs incident to infinity, then rebuild fresh caps around the surviving component boundaries.
  for (size_t face_id = 0; face_id < triangles.size(); ++face_id)
  {
    if (!isLiveFace(face_id) || !faceHasInfiniteVertex(face_id))
    {
      continue;
    }
    killFaceSlot(triangles[face_id], invalid_he);
  }
  for (size_t he_id = 0; he_id < original_half_edge_count; he_id += 2)
  {
    if (half_edges[he_id].origin < 0 || half_edges[he_id ^ 1].origin < 0)
    {
      tombstoneHalfEdge(half_edges, he_id);
      tombstoneHalfEdge(half_edges, he_id ^ 1);
    }
  }

  // --- Collect interior boundary half-edges ---
  // A directed half-edge is on the new component boundary when its own finite face survives, but the twin no longer
  // borders a live face. The twin half-edge is then reused as the exterior side of a new infinite cap.
  std::vector<size_t> boundary_interior;
  boundary_interior.reserve(half_edges.size() / 6);
  for (size_t he_id = 0; he_id < original_half_edge_count; ++he_id)
  {
    if (!isLiveHalfEdge(he_id))
    {
      continue;
    }
    if (half_edges[he_id].origin < 0 || half_edges[he_id ^ 1].origin < 0)
    {
      continue;
    }
    const int face_id = half_edges[he_id].face;
    if (face_id < 0 || !isLiveFace(static_cast<size_t>(face_id)))
    {
      continue;
    }
    if (faceHasInfiniteVertex(static_cast<size_t>(face_id)))
    {
      continue;
    }

    const size_t twin_he = he_id ^ 1;
    if (isLiveHalfEdge(twin_he) && half_edges[twin_he].face >= 0
      && isLiveFace(static_cast<size_t>(half_edges[twin_he].face)))
    {
      continue;
    }

    boundary_interior.push_back(he_id);
  }

  createInfiniteFacesFromBoundary(boundary_interior, half_edges, triangles, vertex_count);

  // --- Tombstone original edge pairs that no longer border any live face ---
  for (size_t he_id = 0; he_id < original_half_edge_count; he_id += 2)
  {
    const int face0 = half_edges[he_id].face;
    const int face1 = half_edges[he_id ^ 1].face;
    const bool face0_live = face0 >= 0 && static_cast<size_t>(face0) < triangles.size()
      && isLiveFace(static_cast<size_t>(face0));
    const bool face1_live = face1 >= 0 && static_cast<size_t>(face1) < triangles.size()
      && isLiveFace(static_cast<size_t>(face1));
    if (!kill_he[he_id] && (face0_live || face1_live))
    {
      continue;
    }
    tombstoneHalfEdge(half_edges, he_id);
    tombstoneHalfEdge(half_edges, he_id ^ 1);
  }

  // Repair pass: some half-edges may still reference a killed face slot after the above edits.
  // If a live half-edge's face is dead, search all live faces for one whose 3-cycle contains this half-edge.
  for (size_t he_id = 0; he_id < half_edges.size(); ++he_id)
  {
    if (!isLiveHalfEdge(he_id))
    {
      continue;
    }
    const int face_id = half_edges[he_id].face;
    if (face_id >= 0 && isLiveFace(static_cast<size_t>(face_id)))
    {
      continue;
    }
    for (size_t live_face = 0; live_face < triangles.size(); ++live_face)
    {
      if (!isLiveFace(live_face))
      {
        continue;
      }
      if (halfEdgeInFaceCycle(he_id, live_face))
      {
        half_edges[he_id].face = static_cast<int>(live_face);
        break;
      }
    }
  }

  // Refresh vertex -> outgoing half-edge (first live half-edge seen per origin; not necessarily CCW-minimal).
  vertex_to_half_edge.assign(vertex_count, static_cast<size_t>(-1));
  for (size_t he_id = 0; he_id < half_edges.size(); ++he_id)
  {
    const int origin = half_edges[he_id].origin;
    if (origin < 0)
    {
      continue;
    }
    const size_t v = static_cast<size_t>(origin);
    if (vertex_to_half_edge[v] == static_cast<size_t>(-1))
    {
      vertex_to_half_edge[v] = he_id;
    }
  }

  rebuildLiveIndices();
  validateLiveHalfEdgeFaceReferences();

  {
    std::ofstream after_dump("applyComponentSplit_after_" + debug_tag + ".txt");
    if (after_dump.is_open())
    {
      printDebug(&after_dump);
    }
  }
}

void HalfEdgeDelaunayGraph::killVertices(const std::vector<bool>& dead_vertex_mask)
{
  if (dead_vertex_mask.size() < vertex_count)
  {
    throw std::runtime_error("killVertices: dead_vertex_mask size mismatch");
  }

  const size_t invalid_he = half_edges.size();
  half_edges.emplace_back();
  half_edges.emplace_back();

  const size_t original_half_edge_count = invalid_he;

  for (size_t face_id = 0; face_id < triangles.size(); ++face_id)
  {
    if (!isLiveFace(face_id))
    {
      continue;
    }

    const std::array<int, 3> tri_vertices = getTriangleVertexIndices(face_id);
    bool touches_dead_vertex = false;
    for (int vertex : tri_vertices)
    {
      if (vertex >= 0 && dead_vertex_mask[static_cast<size_t>(vertex)])
      {
        touches_dead_vertex = true;
        break;
      }
    }

    if (touches_dead_vertex)
    {
      killFaceSlot(triangles[face_id], invalid_he);
    }
  }

  for (size_t he_id = 0; he_id < original_half_edge_count; he_id += 2)
  {
    const int origin0 = half_edges[he_id].origin;
    const int origin1 = half_edges[he_id ^ 1].origin;
    const bool touches_dead_vertex
      = (origin0 >= 0 && dead_vertex_mask[static_cast<size_t>(origin0)])
      || (origin1 >= 0 && dead_vertex_mask[static_cast<size_t>(origin1)]);

    const int face0 = half_edges[he_id].face;
    const int face1 = half_edges[he_id ^ 1].face;
    const bool face0_live = face0 >= 0 && static_cast<size_t>(face0) < triangles.size()
      && isLiveFace(static_cast<size_t>(face0));
    const bool face1_live = face1 >= 0 && static_cast<size_t>(face1) < triangles.size()
      && isLiveFace(static_cast<size_t>(face1));

    if (!touches_dead_vertex && (face0_live || face1_live))
    {
      continue;
    }

    tombstoneHalfEdge(half_edges, he_id);
    tombstoneHalfEdge(half_edges, he_id ^ 1);
  }

  vertex_to_half_edge.assign(vertex_count, static_cast<size_t>(-1));
  for (size_t he_id = 0; he_id < half_edges.size(); ++he_id)
  {
    const int origin = half_edges[he_id].origin;
    if (origin < 0)
    {
      continue;
    }

    const size_t vertex = static_cast<size_t>(origin);
    if (dead_vertex_mask[vertex])
    {
      continue;
    }

    if (vertex_to_half_edge[vertex] == static_cast<size_t>(-1))
    {
      vertex_to_half_edge[vertex] = he_id;
    }
  }

  rebuildLiveIndices();
  validateLiveHalfEdgeFaceReferences();
}

void kinDS::HalfEdgeDelaunayGraph::update(
  size_t vertex_count_,
  const std::vector<std::vector<size_t>>& components,
  const std::function<glm::dvec2(size_t)>& vertex_position)
{
  vertex_count = vertex_count_;
  vertex_to_half_edge.assign(vertex_count, -1);

  std::vector<size_t> index_buffer;

  for (const auto& c : components)
  {
    std::vector<float> coords;
    coords.reserve(c.size() * 2);
    for (const auto& v : c)
    {
      glm::dvec2 point = vertex_position(v);
      coords.push_back(point[0]);
      coords.push_back(point[1]);
    }

    Delaunator2D delaunator(coords);

    for (size_t i : delaunator.triangles)
    {
      index_buffer.emplace_back(c[i]);
    }
  }

  // copy old data:
  auto old_triangles = triangles;
  auto old_half_edges = half_edges;

  build(index_buffer);

  // reorder to match old triangulation
  reorder_from_old(old_triangles, old_half_edges);

  validateLiveHalfEdgeFaceReferences();
}

glm::dvec2 HalfEdgeDelaunayGraph::circumcenter(const glm::dvec2& a, const glm::dvec2& b, const glm::dvec2& c)
{
  // Calculate the circumcenter of the triangle formed by points a, b, c
  double D = 2 * (a[0] * (b[1] - c[1]) + b[0] * (c[1] - a[1]) + c[0] * (a[1] - b[1]));
  if (D == 0)
  {
    // Degenerate case, return a point at infinity
    KINDS_ERROR("Circumcenter calculation failed due to zero denominator. Points may be collinear.");
    return glm::dvec2 { std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity() };
  }
  double Ux = ((a[0] * a[0] + a[1] * a[1]) * (b[1] - c[1]) + (b[0] * b[0] + b[1] * b[1]) * (c[1] - a[1])
                + (c[0] * c[0] + c[1] * c[1]) * (a[1] - b[1]))
    / D;
  double Uy = ((a[0] * a[0] + a[1] * a[1]) * (c[0] - b[0]) + (b[0] * b[0] + b[1] * b[1]) * (a[0] - c[0])
                + (c[0] * c[0] + c[1] * c[1]) * (b[0] - a[0]))
    / D;
  return { Ux, Uy };
}

std::optional<glm::dvec2> HalfEdgeDelaunayGraph::infiniteVoronoiRayDirection(
  size_t face_id, const std::function<glm::dvec2(int vertex_index)>& vertex_at) const
{
  if (!isLiveFace(face_id))
  {
    return std::nullopt;
  }
  for (size_t he_id : face(face_id).half_edges)
  {
    if (!isOnConvexBoundaryOutside(he_id))
    {
      continue;
    }

    const int origin = half_edges[he_id].origin;
    const int dest = destination(he_id);
    if (origin < 0 || dest < 0)
    {
      continue;
    }

    const int interior_vertex = triangleOppositeVertex(he_id ^ 1);
    if (interior_vertex < 0)
    {
      continue;
    }

    const glm::dvec2 u = vertex_at(origin);
    const glm::dvec2 v = vertex_at(dest);
    const glm::dvec2 w = vertex_at(interior_vertex);
    const glm::dvec2 edge = v - u;
    const glm::dvec2 midpoint = 0.5 * (u + v);
    const glm::dvec2 outward_hint = midpoint - w;

    glm::dvec2 dir(edge.y, -edge.x);
    if (glm::dot(dir, outward_hint) < 0.0)
    {
      dir = -dir;
    }
    return dir;
  }

  return std::nullopt;
}

std::vector<std::pair<glm::dvec2, bool>> HalfEdgeDelaunayGraph::computeCircumcenters(
  const std::vector<glm::dvec2>& vertices) const
{
  // give either the position of the circumcenter or a direction vector if the triangle is infinite, the boolean
  // indicates if the circumcenter is infinite
  std::vector<std::pair<glm::dvec2, bool>> circumcenters(triangles.size());

  const auto vertex_at = [&vertices](int vertex_index) { return vertices[static_cast<size_t>(vertex_index)]; };

  for (size_t triangle_id = 0; triangle_id < triangles.size(); triangle_id++)
  {
    if (!isLiveFace(triangle_id))
    {
      continue;
    }

    if (const std::optional<glm::dvec2> infinite_dir = infiniteVoronoiRayDirection(triangle_id, vertex_at))
    {
      circumcenters[triangle_id] = { *infinite_dir, true };
      continue;
    }

    const Triangle& triangle = triangles[triangle_id];
    const std::array<int, 3> cycle_vertices = adjacentTriangleVertices(triangle.half_edges[0]);
    const glm::dvec2& v0 = vertices[static_cast<size_t>(cycle_vertices[0])];
    const glm::dvec2& v1 = vertices[static_cast<size_t>(cycle_vertices[1])];
    const glm::dvec2& v2 = vertices[static_cast<size_t>(cycle_vertices[2])];
    // Compute a finite circumcenter if the triangle is well-conditioned.
    // If the triangle is nearly colinear, the circumcenter becomes numerically unstable and can explode.
    //
    // In that case, fall back to a direction vector (marked as "infinite") derived from the longest edge.
    // The direction is chosen as the opposite of the infinite-vertex case, i.e. the opposite perpendicular.

    const glm::dvec2 e01 = v1 - v0;
    const glm::dvec2 e12 = v2 - v1;
    const glm::dvec2 e20 = v0 - v2;

    const double len01_2 = glm::dot(e01, e01);
    const double len12_2 = glm::dot(e12, e12);
    const double len20_2 = glm::dot(e20, e20);

    const double max_edge_len2 = [&]() -> double
    {
      double m = len01_2;
      if (len12_2 > m)
        m = len12_2;
      if (len20_2 > m)
        m = len20_2;
      return m;
    }();

    // area2 is proportional to the determinant in the circumcenter formula (2x triangle area).
    const double area2 = std::abs(e01.x * (v2.y - v0.y) - e01.y * (v2.x - v0.x));

    // Threshold relative to edge length scale.
    constexpr double degenerate_eps = 1e-12;
    if (area2 <= degenerate_eps * max_edge_len2)
    {
      // Pick the longest edge direction.
      glm::dvec2 longest_edge_vec = e01; // v1 - v0 by default
      if (len12_2 >= len01_2 && len12_2 >= len20_2)
      {
        longest_edge_vec = e12; // v2 - v1
      }
      else if (len20_2 >= len01_2 && len20_2 >= len12_2)
      {
        longest_edge_vec = e20; // v0 - v2
      }

      // Infinite-vertex case uses: circumcenter_dir = perp(edge_vec) = (dy, -dx)
      // Here we want the opposite direction.
      glm::dvec2 circumcenter_dir { longest_edge_vec[1], -longest_edge_vec[0] };
      circumcenter_dir = -circumcenter_dir;

      circumcenters[triangle_id] = { circumcenter_dir, true };
    }
    else
    {
      // Compute the circumcenter of the triangle (finite case).
      circumcenters[triangle_id] = { circumcenter(v0, v1, v2), false };
    }
  }

  return circumcenters;
}

// utils
bool HalfEdgeDelaunayGraph::isOutsideConvexBoundary(size_t he_id) const
{
  // walk the triangle and check if any vertex is -1
  for (size_t i = 0; i < 3; i++)
  {
    if (half_edges[he_id].origin == -1)
    {
      return true;
    }
    he_id = half_edges[he_id].next;
  }

  return false;
}

bool HalfEdgeDelaunayGraph::isOnConvexBoundary(size_t he_id) const
{
  // XOR this as one half-edge will be inside and the other one outside
  return isOutsideConvexBoundary(he_id) != isOutsideConvexBoundary(he_id ^ 1);
}

bool kinDS::HalfEdgeDelaunayGraph::isOnConvexBoundaryOutside(size_t he_id) const
{
  return isOutsideConvexBoundary(he_id) && isOnConvexBoundary(he_id);
}

bool kinDS::HalfEdgeDelaunayGraph::isInfinite(size_t he_id) const
{
  return (half_edges[he_id].origin == -1) || (destination(he_id) == -1);
}

int HalfEdgeDelaunayGraph::destination(size_t he_id) const { return half_edges[he_id ^ 1].origin; }

int HalfEdgeDelaunayGraph::triangleOppositeVertex(size_t he_id) const
{
  // For debugging, get the half-edge and triangle
  auto& half_edge = half_edges[he_id];
  auto& triangle = triangles[half_edge.face];

  // Returns the vertex opposite to the half-edge in its triangle
  size_t next_he_id = half_edges[he_id].next;
  next_he_id = half_edges[next_he_id].next;
  return half_edges[next_he_id].origin;
}

std::vector<size_t> kinDS::HalfEdgeDelaunayGraph::neighbors(size_t v)
{
  std::vector<size_t> nbrs;

  // Iterate over incident edges
  for (IncidentEdgeIterator it = incidentEdgesBegin(v); it != incidentEdgesEnd(v); ++it)
  {
    size_t he_id = *it;
    nbrs.push_back(destination(he_id));
  }

  return nbrs;
}

std::vector<size_t> kinDS::HalfEdgeDelaunayGraph::inducedNeighbors(size_t v, const std::vector<bool>& face_inside) const
{
  std::vector<size_t> nbrs;

  // Iterate over incident edges
  for (IncidentEdgeIterator it = incidentEdgesBegin(v); it != incidentEdgesEnd(v); ++it)
  {
    size_t he_id = *it;

    // Check if the face on the half-edge or its twin is inside
    if (!face_inside[half_edges[he_id].face] && !face_inside[half_edges[he_id ^ 1].face])
    {
      continue; // skip this neighbor
    }

    nbrs.push_back(destination(he_id));
  }

  return nbrs;
}

std::vector<size_t> kinDS::HalfEdgeDelaunayGraph::inducedNeighborsFromLiveGraph(size_t v) const
{
  std::vector<size_t> nbrs;

  for (IncidentEdgeIterator it = incidentEdgesBegin(v); it != incidentEdgesEnd(v); ++it)
  {
    const size_t he_id = *it;
    if (!isLiveHalfEdge(he_id))
    {
      continue;
    }

    const size_t face0 = static_cast<size_t>(half_edges[he_id].face);
    const size_t face1 = static_cast<size_t>(half_edges[he_id ^ 1].face);
    if (!isLiveFace(face0) && !isLiveFace(face1))
    {
      continue;
    }

    const int dest = destination(he_id);
    if (dest >= 0)
    {
      nbrs.push_back(static_cast<size_t>(dest));
    }
  }

  return nbrs;
}

std::array<int, 3> HalfEdgeDelaunayGraph::adjacentTriangleVertices(size_t he_id) const
{
  // Returns the vertices of the triangle that the half-edge belongs to
  std::array<int, 3> vertices;
  for (size_t i = 0; i < 3; i++)
  {
    vertices[i] = half_edges[he_id].origin;
    he_id = half_edges[he_id].next;
  }
  return vertices;
}

size_t kinDS::HalfEdgeDelaunayGraph::neighborEdgeId(size_t he_id) const { return half_edges[twin(he_id)].next; }

size_t kinDS::HalfEdgeDelaunayGraph::nextOnConvexBoundaryId(size_t he_id) const
{
  size_t next_he_id = half_edges[he_id].next;
  next_he_id = twin(next_he_id);
  next_he_id = half_edges[next_he_id].next;
  return next_he_id;
}

size_t kinDS::HalfEdgeDelaunayGraph::prevOnConvexBoundaryId(size_t he_id) const
{
  size_t prev_he_id = prev(he_id);
  prev_he_id = twin(prev_he_id);
  prev_he_id = prev(prev_he_id);
  return prev_he_id;
}

size_t HalfEdgeDelaunayGraph::twin(size_t he_id) { return he_id ^ 1; }

size_t HalfEdgeDelaunayGraph::prev(size_t he_id) const
{
  // Walk around face until we reach the edge from behind
  size_t next_he_id = half_edges[he_id].next;
  size_t prev_he_id;

  while (next_he_id != he_id)
  {
    prev_he_id = next_he_id;
    next_he_id = half_edges[next_he_id].next;
  }

  assert(half_edges[prev_he_id].next == he_id);
  return prev_he_id;
}

std::array<int, 3> kinDS::HalfEdgeDelaunayGraph::getTriangleVertexIndices(size_t face_id) const
{
  std::array<int, 3> result;

  for (size_t i = 0; i < 3; i++)
  {
    size_t he_id = triangles[face_id].half_edges[i];
    result[i] = half_edges[he_id].origin;
  }
  return result;
}

std::array<size_t, 3> kinDS::HalfEdgeDelaunayGraph::getTriangleHalfEdgeIndices(size_t start_he_id) const
{
  std::array<size_t, 3> result;
  size_t he_id = start_he_id;
  for (size_t i = 0; i < 3; i++)
  {
    result[i] = he_id;
    he_id = half_edges[he_id].next;
  }
  return result;
}

std::array<size_t, 4> kinDS::HalfEdgeDelaunayGraph::getQuadBoundaryHalfEdgeIndices(size_t quad_id) const
{
  std::array<size_t, 4> result;

  result[0] = half_edges[quad_id * 2].next; // Next half-edge in the quadrilateral
  result[1] = half_edges[result[0]].next; // Next half-edge in the quadrilateral
  result[2] = half_edges[quad_id * 2 + 1].next; // Next half-edge in the quadrilateral
  result[3] = half_edges[result[2]].next; // Next half-edge in the quadrilateral

  return result;
}

// getters
size_t HalfEdgeDelaunayGraph::getVertexCount() const { return vertex_count; }

static std::array<int, 3> triangle_vertex_key(
  size_t tri_idx, const std::vector<HalfEdgeDelaunayGraph::Triangle>& tris,
  const std::vector<HalfEdgeDelaunayGraph::HalfEdge>& hes)
{
  std::array<int, 3> key;
  for (int i = 0; i < 3; ++i)
  {
    key[i] = hes[tris[tri_idx].half_edges[i]].origin;
  }
  std::sort(key.begin(), key.end());
  return key;
}

struct TriangleVertexKeyHash
{
  size_t operator()(const std::array<int, 3>& k) const noexcept
  {
    size_t h = 1469598103934665603ull;
    for (int v : k)
    {
      h ^= static_cast<size_t>(v) + 0x9e3779b97f4a7c15ull + (h << 6) + (h >> 2);
    }
    return h;
  }
};

static bool is_live_face_slot(
  size_t face_id, const std::vector<HalfEdgeDelaunayGraph::Triangle>& tris,
  const std::vector<HalfEdgeDelaunayGraph::HalfEdge>& hes)
{
  if (face_id >= tris.size())
  {
    return false;
  }

  const HalfEdgeDelaunayGraph::Triangle& tri = tris[face_id];
  if (tri.half_edges[0] >= hes.size())
  {
    return false;
  }

  for (int i = 0; i < 3; ++i)
  {
    const size_t he_id = tri.half_edges[i];
    if (he_id >= hes.size())
    {
      return false;
    }
    if (hes[he_id].face != static_cast<int>(face_id))
    {
      return false;
    }
  }
  return true;
}

static void map_triangle_half_edges_by_directed_edge(
  size_t new_tri_idx, size_t old_tri_idx, const std::vector<HalfEdgeDelaunayGraph::Triangle>& new_tris,
  const std::vector<HalfEdgeDelaunayGraph::HalfEdge>& new_hes,
  const std::vector<HalfEdgeDelaunayGraph::Triangle>& old_tris,
  const std::vector<HalfEdgeDelaunayGraph::HalfEdge>& old_hes, std::vector<int>& he_remap)
{
  std::array<bool, 3> old_slot_used = { false, false, false };

  for (int i = 0; i < 3; ++i)
  {
    const size_t new_he = new_tris[new_tri_idx].half_edges[i];
    if (he_remap[new_he] != -1)
    {
      continue;
    }

    const int new_origin = new_hes[new_he].origin;
    const int new_destination = new_hes[new_he ^ 1].origin;

    for (int j = 0; j < 3; ++j)
    {
      if (old_slot_used[j])
      {
        continue;
      }

      const size_t old_he = old_tris[old_tri_idx].half_edges[j];
      if (old_hes[old_he].origin == new_origin && old_hes[old_he ^ 1].origin == new_destination)
      {
        he_remap[new_he] = static_cast<int>(old_he);
        old_slot_used[j] = true;
        break;
      }
    }
  }
}

void HalfEdgeDelaunayGraph::reorder_from_old(
  const std::vector<Triangle>& old_triangles, const std::vector<HalfEdge>& old_half_edges)
{
  // Fresh per-branch retriangulation is authoritative for convex-hull / infinite-face incidence.
  const std::vector<HalfEdge> fresh_half_edges = half_edges;

  // Map sorted triangle vertex triple -> old triangle index (each old triangle used at most once).
  std::unordered_map<std::array<int, 3>, size_t, TriangleVertexKeyHash> old_tri_map;

  for (size_t ti = 0; ti < old_triangles.size(); ++ti)
  {
    if (!is_live_face_slot(ti, old_triangles, old_half_edges))
    {
      continue;
    }
    old_tri_map[triangle_vertex_key(ti, old_triangles, old_half_edges)] = ti;
  }

  const size_t old_tri_count = old_triangles.size();
  const size_t old_he_count = old_half_edges.size();

  std::vector<int> tri_remap(triangles.size(), -1);
  std::vector<int> he_remap(half_edges.size(), -1);

  size_t next_new_tri = old_tri_count;
  size_t next_new_he = old_he_count;

  // Match unchanged triangles and preserve their old face / half-edge indices.
  for (size_t ti = 0; ti < triangles.size(); ++ti)
  {
    const auto key = triangle_vertex_key(ti, triangles, half_edges);
    const auto it = old_tri_map.find(key);
    if (it == old_tri_map.end())
    {
      continue;
    }

    const size_t old_ti = it->second;
    old_tri_map.erase(it);

    tri_remap[ti] = static_cast<int>(old_ti);
    map_triangle_half_edges_by_directed_edge(ti, old_ti, triangles, half_edges, old_triangles, old_half_edges, he_remap);
  }

  for (size_t ti = 0; ti < triangles.size(); ++ti)
  {
    if (tri_remap[ti] != -1)
    {
      continue;
    }

    tri_remap[ti] = static_cast<int>(next_new_tri);
    ++next_new_tri;
  }

  for (size_t ti = 0; ti < triangles.size(); ++ti)
  {
    for (int i = 0; i < 3; ++i)
    {
      const size_t new_he = triangles[ti].half_edges[i];
      if (he_remap[new_he] != -1)
      {
        continue;
      }

      const size_t twin_he = new_he ^ 1;
      if (he_remap[twin_he] != -1)
      {
        he_remap[new_he] = he_remap[twin_he] ^ 1;
        continue;
      }

      he_remap[new_he] = static_cast<int>(next_new_he);
      he_remap[twin_he] = static_cast<int>(next_new_he ^ 1);
      next_new_he += 2;
    }
  }

  std::vector<Triangle> new_triangles(next_new_tri);
  std::vector<HalfEdge> new_half_edges(next_new_he);

  for (size_t ti = 0; ti < triangles.size(); ++ti)
  {
    const int dst_tri = tri_remap[ti];
    assert(dst_tri >= 0 && static_cast<size_t>(dst_tri) < new_triangles.size());

    Triangle t;
    for (int i = 0; i < 3; ++i)
    {
      const int remapped_he = he_remap[triangles[ti].half_edges[i]];
      assert(remapped_he >= 0 && static_cast<size_t>(remapped_he) < new_half_edges.size());
      t.half_edges[i] = static_cast<size_t>(remapped_he);
    }
    new_triangles[static_cast<size_t>(dst_tri)] = t;
  }

  for (size_t hi = 0; hi < half_edges.size(); ++hi)
  {
    const int dst_he = he_remap[hi];
    if (dst_he < 0)
    {
      continue;
    }
    assert(static_cast<size_t>(dst_he) < new_half_edges.size());

    HalfEdge he = half_edges[hi];
    if (he.next >= 0)
    {
      const int remapped_next = he_remap[static_cast<size_t>(he.next)];
      assert(remapped_next >= 0);
      he.next = remapped_next;
    }
    if (he.face >= 0)
    {
      he.face = tri_remap[static_cast<size_t>(he.face)];
    }
    new_half_edges[static_cast<size_t>(dst_he)] = he;
  }

  triangles.swap(new_triangles);
  half_edges.swap(new_half_edges);

  vertex_to_half_edge.assign(vertex_count, size_t(-1));

  for (size_t he = 0; he < half_edges.size(); ++he)
  {
    const int origin = half_edges[he].origin;
    if (origin < 0)
    {
      continue;
    }
    const size_t v = static_cast<size_t>(origin);
    if (vertex_to_half_edge[v] == size_t(-1))
    {
      vertex_to_half_edge[v] = he;
    }
  }

  // Remapping preserves ids but can leave .face inconsistent with the fresh retriangulation.
  // Re-apply face indices from the per-branch build through tri_remap, matching directed site edges.
  applyFreshHalfEdgeFaceReferences(fresh_half_edges, tri_remap);

  auto remap_face = [&](int fresh_face) -> int
  {
    if (fresh_face < 0)
    {
      return -1;
    }
    const size_t f = static_cast<size_t>(fresh_face);
    return (f < tri_remap.size()) ? tri_remap[f] : -1;
  };

  // Half-edges incident to the vertex at infinity have no finite directed site edge key.
  for (size_t fresh_hi = 0; fresh_hi < fresh_half_edges.size(); ++fresh_hi)
  {
    if (fresh_half_edges[fresh_hi].origin != -1)
    {
      continue;
    }
    const int dst_he = he_remap[fresh_hi];
    if (dst_he < 0 || static_cast<size_t>(dst_he) >= half_edges.size())
    {
      continue;
    }
    half_edges[static_cast<size_t>(dst_he)].face = remap_face(fresh_half_edges[fresh_hi].face);
  }

  rebuildLiveIndices();
}
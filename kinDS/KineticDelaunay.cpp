#include "KineticDelaunay.hpp"
#include "../Logger.hpp"

using namespace kinDS;
uint64_t DelaunayGraph::edge_key(int a, int b) {
  return (static_cast<uint64_t>(a) << 32) | static_cast<uint64_t>(b);
}

void DelaunayGraph::build(const std::vector<size_t>& index_buffer) {
  logger.log(INFO, "Building half-edge mesh from triangle index buffer of size %zu", index_buffer.size());
  assert(index_buffer.size() % 3 == 0 && "Input must be a triangle index buffer.");

  half_edges.clear();
  faces.clear();

  const int num_tris = index_buffer.size() / 3;
  faces.resize(num_tris);

  std::vector<int> boundary_edge_map(vertex_count, -1);

  struct TmpHalfEdge {
    size_t v;
    size_t f;
    size_t j; // index of the half-edge in the triangle
  };

  // Sort up- and down-edges into an adjacency list by the higher vertex index
  std::vector<std::vector<TmpHalfEdge>> adjacency_list_up_in(vertex_count);
  std::vector<std::vector<TmpHalfEdge>> adjacency_list_down_out(vertex_count);

  //iterate through the triangles
  for (size_t i = 0; i < num_tris; ++i) {
    size_t v[3];
    v[0] = index_buffer[i * 3];
    v[1] = index_buffer[i * 3 + 1];
    v[2] = index_buffer[i * 3 + 2];


    // iterate through the triangle's edges
    for (size_t j = 0; j < 3; ++j) {
      size_t a = v[j];
      size_t b = v[(j + 1) % 3];

      if (a > b)
      {
        adjacency_list_down_out[a].push_back({ b,i, j });
      }
      else {
        adjacency_list_up_in[b].push_back({ a,i, j });
      }
    }
  }

  // now sort the opposite way to bring them into order
  std::vector<std::vector<TmpHalfEdge>> adjacency_list_up_out(vertex_count);
  std::vector<std::vector<TmpHalfEdge>> adjacency_list_down_in(vertex_count);

  for (size_t u = 0; u < vertex_count; u++)
  {
    for (TmpHalfEdge& the : adjacency_list_up_in[u]) {
      adjacency_list_up_out[the.v].push_back({ u, the.f, the.j });
    }

    for (TmpHalfEdge& the : adjacency_list_down_out[u]) {
      adjacency_list_down_in[the.v].push_back({ u, the.f, the.j });
    }
  }

  // finally, iterate over both adjacency lists to build the half-edges and faces
  for (size_t u = 0; u < vertex_count; u++) {
    // go through both adjacency lists in a merge-like way
    auto up_it = adjacency_list_up_out[u].begin();
    auto down_it = adjacency_list_down_in[u].begin();

    while (up_it != adjacency_list_up_out[u].end() || down_it != adjacency_list_down_in[u].end()) {
      int v_up = (up_it != adjacency_list_up_out[u].end()) ? up_it->v : std::numeric_limits<int>::max();
      int v_down = (down_it != adjacency_list_down_in[u].end()) ? down_it->v : std::numeric_limits<int>::max();

      // always create both half-edges, but if one of them is not present, it will be a boundary half-edge
      HalfEdge he_up;
      he_up.origin = u;
      he_up.face = (v_up <= v_down) ? up_it->f : -1; // -1 means boundary
      he_up.next = -1; // will be set later

      if (he_up.face != -1) {
        faces[he_up.face].half_edges[up_it->j] = half_edges.size(); // store the half-edge index in the face
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

      if (he_down.face != -1) {
        faces[he_down.face].half_edges[down_it->j] = half_edges.size(); // store the half-edge index in the face
      }
      else
      {
        boundary_edge_map[v] = half_edges.size(); // store the boundary edge for this vertex
      }
      half_edges.push_back(he_down);

      // increment the iterators
      if (v_up <= v_down) {
        up_it++;
      }
      if (v_down <= v_up) {
        down_it++;
      }
    }
  }

  // now link the half-edges together
  // triangle edges first
  for (auto& face : faces)
  {
    for (size_t j = 0; j < 3; ++j) {
      size_t he_index = face.half_edges[j];
      if (he_index != -1) {
        HalfEdge& he = half_edges[he_index];
        he.next = face.half_edges[(j + 1) % 3]; // link to the next half-edge in the face
      }
      else
      {
        logger.log(WARNING, "Face %zu has a missing half-edge at index %zu.", &face - &faces[0], j);
      }
    }
  }

  // continue with boundary edges
  for (size_t u = 0; u < vertex_count; ++u) {
    int boundary_edge_index = boundary_edge_map[u];
    if (boundary_edge_index != -1) {
      HalfEdge& he = half_edges[boundary_edge_index];
      // we need the twin to get the destination vertex
      HalfEdge& twin_he = half_edges[boundary_edge_index ^ 1]; // twin is the next half-edge in the list

      he.next = boundary_edge_map[twin_he.origin]; // link to the next boundary half-edge
      if (he.next == -1) {
        logger.log(WARNING, "Boundary edge for vertex %zu is not linked to another boundary edge.", u);
      }
    }
  }

  logger.log(INFO, "Half-edge mesh built with %zu half-edges and %zu faces.", half_edges.size(), faces.size());
}

#include <iostream>

void DelaunayGraph::print_debug() const {
  std::cout << "Half-Edges:\n";
  for (size_t i = 0; i < half_edges.size(); ++i) {
    const HalfEdge& he = half_edges[i];
    std::cout << "  [" << i << "] origin = " << he.origin
      << ", next = " << he.next
      << ", face = " << he.face
      << ", twin = " << (i ^ 1) << "\n";
  }

  std::cout << "\nFaces:\n";
  for (size_t i = 0; i < faces.size(); ++i) {
    const Face& f = faces[i];
    std::cout << "  [" << i << "] half_edge = " << f.half_edges[0] << "\n";

    // Walk the face's boundary
    int start = f.half_edges[0];
    if (start < 0) continue;
    int he = start;
    std::cout << "    Vertices: ";
    do {
      std::cout << half_edges[he].origin << " ";
      he = half_edges[he].next;
    } while (he != start && he != -1);
    std::cout << "\n";
  }
}

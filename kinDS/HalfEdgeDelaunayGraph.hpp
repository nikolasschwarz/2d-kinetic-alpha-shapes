#pragma once
#include "../delaunay/Delaunator2D.hpp"
#include "CubicHermiteSpline.hpp"
#include <array>
#include <vector>

namespace kinDS {
  // A graph that can represent the delaunay triangulation such that edges are explicitly stored and can be flipped in its quadrilateral.
  class DelaunayGraph {
  public:

    struct HalfEdge {
      int origin = -1;    // index into vertices
      int next = -1;      // index into half_edges
      int face = -1;      // index into faces; -1 means boundary
      // twin = index ^ 1
    };

    struct Face {
      std::array<size_t, 3> half_edges;
    };

  private:
    size_t vertex_count = 0; // Number of vertices in the triangulation
    std::vector<Face> faces;
    std::vector<HalfEdge> half_edges; // List of half-edges in the triangulation

    // Builds the half-edge mesh from a triangle index buffer (size must be multiple of 3)
    void build(const std::vector<size_t>& index_buffer);
    static uint64_t edge_key(int a, int b); // encodes directed edge (a->b)

    const std::vector<HalfEdge>& get_half_edges() const { return half_edges; }
    const std::vector<Face>& get_faces() const { return faces; }

    // utils
    inline size_t destination(size_t he_id) const {
      return half_edges[he_id ^ 1].origin;
    }

    size_t triangle_opposite_vertex(size_t he_id) const {
      // Returns the vertex opposite to the half-edge in its triangle
      size_t next_he_id = half_edges[he_id].next;
      next_he_id = half_edges[next_he_id].next;
      return half_edges[next_he_id].origin;
    }


  public:
    DelaunayGraph() = default;

    void init(const std::vector<CubicHermiteSpline<2>>& splines) {
      vertex_count = splines.size();
      std::vector<float> coords;
      coords.reserve(splines.size() * 2); // Reserve space for x and y coordinates
      for (const auto& spline : splines) {
        Point<2> point = spline.evaluate(0.0);
        coords.push_back(point[0]);
        coords.push_back(point[1]);
      }

      Delaunator::Delaunator2D delaunator(coords);

      build(delaunator.triangles);

    }
    // Flips an edge between two triangles
    void flipEdge(size_t he_id);
    // Other methods to manipulate and query the triangulation can be added here.
    void print_debug() const;

  };
}
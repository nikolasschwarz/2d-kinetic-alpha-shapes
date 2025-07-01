#pragma once
#include "Polynomial.hpp"
#include "CubicHermiteSpline.hpp"
#include "../delaunay/Delaunator2D.hpp"

namespace kinDS
{
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
    void flipEdge(size_t v0, size_t v1);
    // Other methods to manipulate and query the triangulation can be added here.
    void print_debug() const;

  };

class KineticDelaunay {
private:
  std::vector<CubicHermiteSpline<2>> splines;


  /* Compare to Leonidas Guibas and Jorge Stolfi. 1985. Primitives for the manipulation of general subdivisions and the computation of Voronoi. ACM Trans. Graph. 4, 2 (April 1985), 74–123. https://doi.org/10.1145/282918.282923
   */
  static Polynomial inCircle(const Polynomial& ax, const Polynomial& ay, const Polynomial& bx, const Polynomial& by, const Polynomial& cx, const Polynomial& cy,
    const Polynomial& px, const Polynomial& py) {
    const Polynomial dx = ax - px;
    const Polynomial dy = ay - py;
    const Polynomial ex = bx - px;
    const Polynomial ey = by - py;
    const Polynomial fx = cx - px;
    const Polynomial fy = cy - py;

    const Polynomial ap = dx * dx + dy * dy;
    const Polynomial bp = ex * ex + ey * ey;
    const Polynomial cp = fx * fx + fy * fy;

    return (dx * (ey * cp - bp * fy) - dy * (ex * cp - bp * fx) + ap * (ex * fy - ey * fx));
  }

  static Polynomial ccw(const Polynomial& ax, const Polynomial& ay, const Polynomial& bx, const Polynomial& by, const Polynomial& cx, const Polynomial& cy) {
    return (ax * by) + (bx * cy) + (cx * ay) - (ay * bx) - (by * cx) - (cy * ax);
  }

  void computeStep(double t) {

  }

public:
  KineticDelaunay(const std::vector<CubicHermiteSpline<2>>& splines) : splines(splines) {}
  // Computes the Delaunay triangulation of the given splines
  void compute() {
    DelaunayGraph graph;
    graph.init(splines);

    graph.print_debug();
  }
  // Other methods to manipulate and query the triangulation can be added here.
};
} // namespace kinDS
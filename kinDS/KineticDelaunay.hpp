#pragma once
#include "Polynomial.hpp"
#include "CubicHermiteSpline.hpp"
#include "HalfEdgeDelaunayGraph.hpp"

namespace kinDS
{

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

    graph.flipEdge(6); // Example of flipping the first edge

    graph.print_debug(); // Print the graph after flipping the edge
  }
  // Other methods to manipulate and query the triangulation can be added here.
};
} // namespace kinDS
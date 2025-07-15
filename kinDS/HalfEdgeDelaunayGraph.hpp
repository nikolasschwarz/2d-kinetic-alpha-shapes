#pragma once
#include "../delaunay/Delaunator2D.hpp"
#include "CubicHermiteSpline.hpp"
#include "Logger.hpp"
#include <array>
#include <vector>

namespace kinDS
{
// A graph that can represent the delaunay triangulation such that edges are explicitly stored and can be flipped in its quadrilateral.
class HalfEdgeDelaunayGraph
{
 public:
  struct HalfEdge
  {
    int origin = -1; // index into vertices
    int next = -1; // index into half_edges
    int face = -1; // index into faces; -1 means boundary
    // twin = index ^ 1
  };

  struct Triangle
  {
    std::array<size_t, 3> half_edges;
  };

 private:
  size_t vertex_count = 0; // Number of vertices in the triangulation
  std::vector<Triangle> triangles;
  std::vector<HalfEdge> half_edges; // List of half-edges in the triangulation

  /**
   * @brief Build the data structure from an index buffer of triangles.
   *
   * Constructs the half-edge data structure from a list of triangle indices. Construction is done in O(n) time.
   * Twin edges are stored implicitly by storing them next to each other. The twin index can be computed as `index ^ 1`.
   *
   * @param index_buffer index buffer of triangles, size must be multiple of 3 to be valid and each index must be in range [0, vertex_count).
   */
  void build(const std::vector<size_t>& index_buffer);

 public:
  HalfEdgeDelaunayGraph() = default;

  void init(const std::vector<CubicHermiteSpline<2>>& splines)
  {
    vertex_count = splines.size();
    std::vector<float> coords;
    coords.reserve(splines.size() * 2); // Reserve space for x and y coordinates
    for (const auto& spline : splines)
    {
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
  void printDebug() const;

  static Point<2> circumcenter(const Point<2>& a, const Point<2>& b, const Point<2>& c)
  {
    // Calculate the circumcenter of the triangle formed by points a, b, c
    double D = 2 * (a[0] * (b[1] - c[1]) + b[0] * (c[1] - a[1]) + c[0] * (a[1] - b[1]));
    if (D == 0)
    {
      // Degenerate case, return a point at infinity
      logger.log(ERROR, "Circumcenter calculation failed due to zero denominator. Points may be collinear.");
    }
    double Ux = ((a[0] * a[0] + a[1] * a[1]) * (b[1] - c[1]) + (b[0] * b[0] + b[1] * b[1]) * (c[1] - a[1]) + (c[0] * c[0] + c[1] * c[1]) * (a[1] - b[1])) / D;
    double Uy = ((a[0] * a[0] + a[1] * a[1]) * (c[0] - b[0]) + (b[0] * b[0] + b[1] * b[1]) * (a[0] - c[0]) + (c[0] * c[0] + c[1] * c[1]) * (b[0] - a[0])) / D;
    return { Ux, Uy };
  }

  std::vector<std::pair<Point<2>, bool>> computeCircumcenters(const std::vector<Point<2>>& vertices) const
  {
    // give either the position of the circumcenter or a direction vector if the triangle is infinite, the boolean indicates if the circumcenter is infinite
    std::vector<std::pair<Point<2>, bool>> circumcenters(triangles.size());

    for (size_t triangle_id = 0; triangle_id < triangles.size(); triangle_id++)
    {
      const Triangle& triangle = triangles[triangle_id];
      const HalfEdge& he0 = half_edges[triangle.half_edges[0]];
      const HalfEdge& he1 = half_edges[triangle.half_edges[1]];
      const HalfEdge& he2 = half_edges[triangle.half_edges[2]];

      // TODO:: How do we obtain the correct direction for these infinite vertices?
      // Filter for infinity vertices
      if (he0.origin == -1)
      {
        const Point<2>& v1 = vertices[he1.origin];
        const Point<2>& v2 = vertices[he2.origin];
        const Point<2> dir = v2 - v1;
        circumcenters[triangle_id] = { Point<2> { dir[1], -dir[0] }, true };
        continue;
      }

      if (he1.origin == -1)
      {
        const Point<2>& v0 = vertices[he0.origin];
        const Point<2>& v2 = vertices[he2.origin];
        const Point<2> dir = v0 - v2;
        circumcenters[triangle_id] = { Point<2> { dir[1], -dir[0] }, true };
        continue;
      }

      if (he2.origin == -1)
      {
        const Point<2>& v0 = vertices[he0.origin];
        const Point<2>& v1 = vertices[he1.origin];
        const Point<2> dir = v1 - v0;
        circumcenters[triangle_id] = { Point<2> { dir[1], -dir[0] }, true };
        continue;
      }

      // Get the vertices of the triangle
      const Point<2>& v0 = vertices[he0.origin];
      const Point<2>& v1 = vertices[he1.origin];
      const Point<2>& v2 = vertices[he2.origin];
      // Compute the circumcenter of the triangle
      circumcenters[triangle_id] = { circumcenter(v0, v1, v2), false };
    }

    return circumcenters;
  }

  // utils
  bool isOutsideBoundary(size_t he_id) const
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

  bool isOnBoundary(size_t he_id) const
  {
    // XOR this as one half-edge will be inside and the other one outside
    return isOutsideBoundary(he_id) != isOutsideBoundary(he_id ^ 1);
  }
  inline size_t destination(size_t he_id) const
  {
    return half_edges[he_id ^ 1].origin;
  }
  int triangleOppositeVertex(size_t he_id) const
  {
    // Returns the vertex opposite to the half-edge in its triangle
    size_t next_he_id = half_edges[he_id].next;
    next_he_id = half_edges[next_he_id].next;
    return half_edges[next_he_id].origin;
  }

  std::array<int, 3> adjacentTriangleVertices(size_t he_id) const
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

  // getters
  const std::vector<HalfEdge>& getHalfEdges() const { return half_edges; }
  const std::vector<Triangle>& getFaces() const { return triangles; }
  size_t getVertexCount() const { return vertex_count; }
};
}

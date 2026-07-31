#pragma once

#include <algorithm>
#include <cmath>
#include <glm/glm.hpp>

namespace kinDS
{

/// Non-instantiable helpers for common 3D geometric queries (GLM only, no CGAL).
class GeometryUtils
{
 public:
  GeometryUtils() = delete;
  ~GeometryUtils() = delete;
  GeometryUtils(const GeometryUtils&) = delete;
  GeometryUtils& operator=(const GeometryUtils&) = delete;

  /// Closest point on the closed triangle ABC to query point `p`.
  ///
  /// Algorithm:
  /// 1. Orthogonal projection of `p` onto the supporting plane. If that point
  ///    lies inside the triangle, it is the closest point.
  /// 2. Otherwise the closest point lies on the boundary: take the closest
  ///    point on each of the three closed edge segments (clamping to endpoints
  ///    covers the vertices; no separate vertex pass is needed).
  /// 3. Degenerate (near-zero area) triangles fall back to the edge-segment
  ///    stage only.
  static glm::dvec3 closestPointOnTriangle(const glm::dvec3& p, const glm::dvec3& a, const glm::dvec3& b,
                                           const glm::dvec3& c)
  {
    constexpr double kAreaEps2 = 1e-24;

    const glm::dvec3 ab = b - a;
    const glm::dvec3 ac = c - a;
    const glm::dvec3 n = glm::cross(ab, ac);
    const double n2 = glm::dot(n, n);

    if (n2 > kAreaEps2)
    {
      // Orthogonal projection onto the supporting plane of ABC.
      const glm::dvec3 q = p - n * (glm::dot(n, p - a) / n2);
      if (pointInTrianglePlane(q, a, b, c, n))
      {
        return q;
      }
    }

    // Closest point on the boundary (three closed segments).
    glm::dvec3 best = closestPointOnSegment(p, a, b);
    double best_d2 = squaredDistance(p, best);

    const glm::dvec3 on_bc = closestPointOnSegment(p, b, c);
    const double d2_bc = squaredDistance(p, on_bc);
    if (d2_bc < best_d2)
    {
      best = on_bc;
      best_d2 = d2_bc;
    }

    const glm::dvec3 on_ca = closestPointOnSegment(p, c, a);
    const double d2_ca = squaredDistance(p, on_ca);
    if (d2_ca < best_d2)
    {
      best = on_ca;
    }

    return best;
  }

  /// Euclidean distance from `p` to the closed triangle ABC.
  static double distancePointTriangle(const glm::dvec3& p, const glm::dvec3& a, const glm::dvec3& b,
                                      const glm::dvec3& c)
  {
    return glm::length(p - closestPointOnTriangle(p, a, b, c));
  }

  /// Squared Euclidean distance from `p` to the closed triangle ABC.
  static double squaredDistancePointTriangle(const glm::dvec3& p, const glm::dvec3& a, const glm::dvec3& b,
                                             const glm::dvec3& c)
  {
    return squaredDistance(p, closestPointOnTriangle(p, a, b, c));
  }

  /// Closest point on the closed segment AB to query point `p`.
  static glm::dvec3 closestPointOnSegment(const glm::dvec3& p, const glm::dvec3& a, const glm::dvec3& b)
  {
    const glm::dvec3 ab = b - a;
    const double ab2 = glm::dot(ab, ab);
    if (ab2 <= 0.0)
    {
      return a;
    }
    const double t = std::clamp(glm::dot(p - a, ab) / ab2, 0.0, 1.0);
    return a + t * ab;
  }

 private:
  static double squaredDistance(const glm::dvec3& u, const glm::dvec3& v)
  {
    const glm::dvec3 d = u - v;
    return glm::dot(d, d);
  }

  /// True if `q` (assumed coplanar with ABC) lies in the closed triangle.
  /// Uses consistent half-plane tests against the three edges, oriented by `n`.
  static bool pointInTrianglePlane(const glm::dvec3& q, const glm::dvec3& a, const glm::dvec3& b,
                                   const glm::dvec3& c, const glm::dvec3& n)
  {
    const double s_ab = glm::dot(glm::cross(b - a, q - a), n);
    const double s_bc = glm::dot(glm::cross(c - b, q - b), n);
    const double s_ca = glm::dot(glm::cross(a - c, q - c), n);
    return (s_ab >= 0.0 && s_bc >= 0.0 && s_ca >= 0.0) || (s_ab <= 0.0 && s_bc <= 0.0 && s_ca <= 0.0);
  }
};

} // namespace kinDS

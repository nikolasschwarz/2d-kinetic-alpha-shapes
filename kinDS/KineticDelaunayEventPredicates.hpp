#pragma once

#include "Polynomial.hpp"
#include "StrandTree.hpp"

#include <cmath>
#include <glm/geometric.hpp>

namespace kinDS
{

// Geometric predicate helpers used to build event trigger polynomials.
// These were historically defined inside `KineticDelaunay.cpp` but must be available
// to the per-event factory translation units.

/* Compare to Leonidas Guibas and Jorge Stolfi. 1985. Primitives for the manipulation of general subdivisions and the
 * computation of Voronoi. ACM Trans. Graph. 4, 2 (April 1985), 74-123. https://doi.org/10.1145/282918.282923
 */
static Polynomial inCircle(const Polynomial& ax, const Polynomial& ay, const Polynomial& bx, const Polynomial& by,
  const Polynomial& cx, const Polynomial& cy, const Polynomial& px, const Polynomial& py)
{
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

static Polynomial ccw(const Polynomial& ax, const Polynomial& ay, const Polynomial& bx, const Polynomial& by,
  const Polynomial& cx, const Polynomial& cy)
{
  return (ax * by) + (bx * cy) + (cx * ay) - (ay * bx) - (by * cx) - (cy * ax);
}

// Polynomial that evaluates to zero iff the distance from A to the circumcenter equals the value r
static Polynomial circumradiusEquals(const Polynomial& ax, const Polynomial& ay, const Polynomial& bx,
  const Polynomial& by, const Polynomial& cx, const Polynomial& cy, double r)
{
  // We first do the same computations as for the circumcenter
  Polynomial D = (ax * (by - cy) + bx * (cy - ay) + cx * (ay - by)) * 2.0;

  // only compute the numerators
  Polynomial Nx = ((ax * ax + ay * ay) * (by - cy) + (bx * bx + by * by) * (cy - ay) + (cx * cx + cy * cy) * (ay - by));
  Polynomial Ny = ((ax * ax + ay * ay) * (cx - bx) + (bx * bx + by * by) * (ax - cx) + (cx * cx + cy * cy) * (bx - ax));

  // By taking the distance formula between the circumcenter and the point A and setting it equal to r, we get the
  // following after rearranging:
  Polynomial circumradius_eq = (Nx - ax * D) * (Nx - ax * D) + (Ny - ay * D) * (Ny - ay * D) - (D * D * (r * r));
  return circumradius_eq;
}

// Numeric counterpart of the angular bisector predicate used for infinite / hull cases.
// Returns a value whose sign indicates on which side of the angular bisector at `a`
// (between rays a->c and a->c_prime) the query point `x` lies.
static double angularBisectorPredicate(
  const glm::dvec2& a, const glm::dvec2& c, const glm::dvec2& c_prime, const glm::dvec2& x)
{
  glm::dvec2 ac = c - a;
  glm::dvec2 ac_prime = c_prime - a;
  glm::dvec2 ax = x - a;

  double ac_sq = glm::dot(ac, ac);
  double ac_prime_sq = glm::dot(ac_prime, ac_prime);

  double term1 = glm::dot(ac, ax) * ac_prime_sq;
  double term2 = glm::dot(ac_prime, ax) * ac_sq;

  return term1 - term2;
}

static Polynomial angularBisectorHelper(
  Trajectory<3>& voronoi_homogeneous, Trajectory<2>& a, Trajectory<2>& c_i, Trajectory<2>& c_j)
{
  Trajectory<2> a_scaled = a * voronoi_homogeneous[2];
  Trajectory<2> voronoi_xy { voronoi_homogeneous[0], voronoi_homogeneous[1] };
  return (Trajectory<2>::dot(c_i - a, voronoi_xy - a_scaled) * (c_j - a).squaredNorm());
}
static Polynomial angularBisector(
  Trajectory<2>& a, Trajectory<2>& c, Trajectory<2>& c_prime, Trajectory<3> voronoi_homogeneous)
{
  return angularBisectorHelper(voronoi_homogeneous, a, c, c_prime)
    - angularBisectorHelper(voronoi_homogeneous, a, c_prime, c);
}

/// Upper fraction bound for event-root filtering; matches @ref eventIntervalUpperBound in absolute time.
constexpr double kEventIntervalFractionUpperBound = 1.0;

/// Absolute-time upper bound of the section interval events are computed for: [t, result].
inline double eventIntervalUpperBound(double t)
{
  return std::floor(t) + kEventIntervalFractionUpperBound;
}
} // namespace kinDS


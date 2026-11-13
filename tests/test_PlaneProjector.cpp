#include "kinDS/PlaneProjector.hpp"

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

using namespace kinDS;
using Catch::Approx;

namespace
{
// Profile-plane transform: local (u, 0, v) -> world, columns u / normal / v / origin.
glm::dmat4 planeTransform(const glm::dvec3& origin, const glm::dvec3& u, const glm::dvec3& v)
{
  glm::dmat4 m(1.0);
  m[0] = glm::dvec4(u, 0.0);
  m[1] = glm::dvec4(glm::normalize(glm::cross(v, u)), 0.0);
  m[2] = glm::dvec4(v, 0.0);
  m[3] = glm::dvec4(origin, 1.0);
  return m;
}

glm::dvec3 toWorld(const glm::dmat4& m, const glm::dvec2& p) { return glm::dvec3(m * glm::dvec4(p.x, 0.0, p.y, 1.0)); }

double planeDistance(const glm::dmat4& m, const glm::dvec3& x)
{
  const glm::dvec3 n = glm::normalize(glm::cross(glm::dvec3(m[0]), glm::dvec3(m[2])));
  return glm::dot(n, x - glm::dvec3(m[3]));
}

// Trunk-like horizontal plane and a tilted child plane, both far from the world origin (Elm section 11 scale).
const glm::dmat4 kA = planeTransform({ -0.003, 0.094, -0.001 }, { 0.001, 0.0, 0.0 }, { 0.0, 0.0, 0.001 });
const glm::dmat4 kB = planeTransform({ 0.034, 0.092, -0.014 }, glm::dvec3(0.000995, -0.0001, 0.0) * 1.0,
  glm::normalize(glm::dvec3(0.0, 0.04, 1.0)) * 0.001);
} // namespace

TEST_CASE("PlaneProjector hinge point lies on both planes", "[PlaneProjector]")
{
  const PlaneProjector projector(kA, kB);
  REQUIRE_FALSE(projector.isParallel());
  CHECK(planeDistance(kA, projector.intersectionPoint()) == Approx(0.0).margin(1e-12));
  CHECK(planeDistance(kB, projector.intersectionPoint()) == Approx(0.0).margin(1e-12));
}

TEST_CASE("PlaneProjector folds about the hinge onto the target plane", "[PlaneProjector]")
{
  const PlaneProjector projector(kA, kB);
  const glm::dvec3 p0 = projector.intersectionPoint();
  const glm::dvec3 axis = projector.intersectionAxis();
  const auto hinge_distance = [&](const glm::dvec3& x)
  {
    const glm::dvec3 d = x - p0;
    return glm::length(d - glm::dot(d, axis) * axis);
  };

  for (const glm::dvec2 p : { glm::dvec2(10.0, -5.0), glm::dvec2(-30.0, 20.0), glm::dvec2(45.0, 12.0) })
  {
    const glm::dvec3 on_a = toWorld(kA, p);
    const glm::dvec3 on_b = toWorld(kB, projector.project(p));
    CHECK(planeDistance(kB, on_b) == Approx(0.0).margin(1e-12));
    CHECK(hinge_distance(on_b) == Approx(hinge_distance(on_a)).epsilon(1e-9));
    CHECK(glm::length(PlaneProjector(kB, kA).project(projector.project(p)) - p) == Approx(0.0).margin(1e-9));
  }

  // Points on the hinge stay where they are.
  const auto line = projector.intersectionLineInLocalA();
  REQUIRE(line.has_value());
  const glm::dvec2 on_hinge = line->first + 7.0 * line->second;
  CHECK(glm::length(toWorld(kB, projector.project(on_hinge)) - toWorld(kA, on_hinge)) == Approx(0.0).margin(1e-12));
}

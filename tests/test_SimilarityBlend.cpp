#include "kinDS/SimilarityBlend.hpp"

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <glm/gtc/quaternion.hpp>

using namespace kinDS;
using Catch::Approx;

namespace
{
void checkNear(const glm::dvec3& a, const glm::dvec3& b) { CHECK(glm::length(a - b) == Approx(0.0).margin(1e-12)); }
} // namespace

TEST_CASE("Similarity fit recovers rotation, scale and translation of a coplanar set", "[SimilarityBlend]")
{
  const glm::dquat rotation = glm::angleAxis(0.6, glm::normalize(glm::dvec3(1.0, 2.0, 0.5)));
  const double scale = 1.3;
  const glm::dvec3 translation(0.2, -0.4, 1.5);
  const std::vector<glm::dvec3> pre { { 0.0, 0.0, 0.0 }, { 1.0, 0.0, 0.0 }, { 0.0, 2.0, 0.0 }, { 1.5, 1.0, 0.0 } };
  std::vector<glm::dvec3> post;
  for (const auto& p : pre)
  {
    post.push_back(translation + scale * (rotation * p));
  }

  const SimilarityFit fit = fitSimilarity(pre, post, { 0.0, 0.0, 1.0 }, rotation * glm::dvec3(0.0, 0.0, 1.0));
  CHECK(fit.scale == Approx(scale).epsilon(1e-12));
  checkNear(fit.rotation * glm::dvec3(0.3, -0.7, 0.9), rotation * glm::dvec3(0.3, -0.7, 0.9));

  for (size_t i = 0; i < pre.size(); ++i)
  {
    checkNear(blendSimilarity(fit, pre[i], post[i], 0.0), pre[i]);
    checkNear(blendSimilarity(fit, pre[i], post[i], 1.0), post[i]);
  }
  // Distances change linearly; a plain linear blend would shrink them under rotation.
  const double d_pre = glm::length(pre[1] - pre[2]);
  for (const double w : { 0.25, 0.5, 0.75 })
  {
    const double d = glm::length(blendSimilarity(fit, pre[1], post[1], w) - blendSimilarity(fit, pre[2], post[2], w));
    CHECK(d == Approx(d_pre * (1.0 + w * (scale - 1.0))).epsilon(1e-12));
  }
}

TEST_CASE("Similarity blend keeps exact endpoints with a non-similar residual", "[SimilarityBlend]")
{
  const std::vector<glm::dvec3> pre { { 0.0, 0.0, 0.0 }, { 1.0, 0.0, 0.0 }, { 0.0, 1.0, 0.0 } };
  const std::vector<glm::dvec3> post { { 0.0, 0.0, 0.0 }, { 2.0, 0.0, 0.0 }, { 0.0, 1.0, 0.0 } };
  const SimilarityFit fit = fitSimilarity(pre, post, { 0.0, 0.0, 1.0 }, { 0.0, 0.0, 1.0 });
  for (size_t i = 0; i < pre.size(); ++i)
  {
    checkNear(blendSimilarity(fit, pre[i], post[i], 0.0), pre[i]);
    checkNear(blendSimilarity(fit, pre[i], post[i], 1.0), post[i]);
  }
}

TEST_CASE("Similarity fit of a single point is a pure translation", "[SimilarityBlend]")
{
  const SimilarityFit fit = fitSimilarity({ { 1.0, 2.0, 3.0 } }, { { 3.0, 2.0, 1.0 } }, {}, {});
  CHECK(fit.scale == 1.0);
  checkNear(blendSimilarity(fit, { 1.0, 2.0, 3.0 }, { 3.0, 2.0, 1.0 }, 0.5), { 2.0, 2.0, 2.0 });
}

TEST_CASE("Similarity fit of two points stays finite", "[SimilarityBlend]")
{
  // Two-strand split group from a meshing run: rank-deficient covariance made the polar iteration diverge to NaN.
  const std::vector<glm::dvec3> pre { { 0.39355623741794887, 0.02729327622459192, 0.35063514112889005 },
    { 0.39377076983080961, 0.02543337250198251, 0.34990621710902364 } };
  const std::vector<glm::dvec3> post { { 0.39354495405204099, 0.027304172081389008, 0.35064938050773975 },
    { 0.39375958109237508, 0.025444367488616532, 0.34992009399066448 } };
  // Branch plane normals are cross products of scaled profile axes, so they can be tiny.
  for (const glm::dvec3 normal : { glm::dvec3(2.0163426205112579e-08, 3.9001789521209486e-07, -9.2058779977345704e-07),
         glm::dvec3(0.0, 0.0, 1.0), glm::dvec3(0.0, 0.0, 0.0) })
  {
    const SimilarityFit fit = fitSimilarity(pre, post, normal, normal);
    for (size_t i = 0; i < pre.size(); ++i)
    {
      checkNear(blendSimilarity(fit, pre[i], post[i], 0.0), pre[i]);
      checkNear(blendSimilarity(fit, pre[i], post[i], 1.0), post[i]);
      const glm::dvec3 mid = blendSimilarity(fit, pre[i], post[i], 0.5);
      CHECK((std::isfinite(mid.x) && std::isfinite(mid.y) && std::isfinite(mid.z)));
    }
  }
}

namespace
{
glm::dvec2 rotate(double angle, const glm::dvec2& v)
{
  return { std::cos(angle) * v.x - std::sin(angle) * v.y, std::sin(angle) * v.x + std::cos(angle) * v.y };
}

glm::dvec2 evaluateQuadratic(const std::array<glm::dvec2, 3>& c, double w) { return c[0] + w * (c[1] + w * c[2]); }

// Distance between two static points under the polynomial blend of a frame change by `angle` about `center`.
double blendedDistance(double angle, double w)
{
  const glm::dvec2 a(1.0, 2.0);
  const glm::dvec2 b(3.0, 1.5);
  const glm::dvec2 center(-4.0, 7.0);
  const std::complex<double> similarity = std::polar(1.0, angle);
  const auto blended = [&](const glm::dvec2& p)
  {
    const glm::dvec2 p1 = center + rotate(angle, p - center);
    return evaluateQuadratic(quadraticThroughMidpoint(p, similarityBlendMidpoint(p, p, p1, similarity), p1), w);
  };
  return glm::length(blended(a) - blended(b));
}
} // namespace

TEST_CASE("Polynomial frame blend keeps distances under a frame rotation", "[SimilarityBlend]")
{
  const double d = glm::length(glm::dvec2(1.0, 2.0) - glm::dvec2(3.0, 1.5));
  const double angle = 26.0 * 3.14159265358979323846 / 180.0;
  CHECK(blendedDistance(angle, 0.0) == Approx(d).epsilon(1e-12));
  CHECK(blendedDistance(angle, 0.5) == Approx(d).epsilon(1e-12));
  CHECK(blendedDistance(angle, 1.0) == Approx(d).epsilon(1e-12));
  CHECK(blendedDistance(angle, 0.25) == Approx(d).epsilon(1e-3));
  CHECK(quadraticSimilarityBlendDistanceError(std::polar(1.0, angle)) < 1e-3);

  // A linear blend mirrors every point through one spot halfway at 180 degrees; the polynomial blend stays close.
  for (const double w : { 0.1, 0.25, 0.5, 0.75, 0.9 })
  {
    CHECK(blendedDistance(3.14159265358979323846, w) > 0.85 * d);
  }
}

TEST_CASE("Polynomial frame blend is the linear blend without rotation or scale", "[SimilarityBlend]")
{
  const glm::dvec2 p0(1.0, 2.0);
  const glm::dvec2 p1_start(2.0, -1.0);
  const glm::dvec2 p1(5.0, 0.5); // translated end frame
  const auto c = quadraticThroughMidpoint(p0, similarityBlendMidpoint(p0, p1_start, p1, 1.0), p1);
  CHECK(glm::length(c[2]) == Approx(0.0).margin(1e-12));
  CHECK(glm::length(evaluateQuadratic(c, 0.3) - (p0 + 0.3 * (p1 - p0))) == Approx(0.0).margin(1e-12));
  CHECK(quadraticSimilarityBlendDistanceError(1.0) == Approx(0.0).margin(1e-12));
  CHECK(quadraticSimilarityBlendDistanceError(2.0) < 1e-2); // pure scale: |q(w)| vs 2^w
}

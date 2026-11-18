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

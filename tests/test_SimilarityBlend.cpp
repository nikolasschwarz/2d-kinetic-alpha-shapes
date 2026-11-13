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

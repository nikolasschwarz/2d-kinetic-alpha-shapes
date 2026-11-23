#include "SimilarityBlend.hpp"

#include <algorithm>
#include <cmath>

namespace kinDS
{
namespace
{
double squaredNorm(const glm::dmat3& m) { return glm::dot(m[0], m[0]) + glm::dot(m[1], m[1]) + glm::dot(m[2], m[2]); }

/// Rotation factor of the polar decomposition (Newton iteration); identity unless @p m is a well-conditioned proper
/// rotation-scale. Rank-deficient input (e.g. two points) has no unique rotation and makes the iteration diverge.
glm::dmat3 polarRotation(glm::dmat3 m)
{
  const double det = glm::determinant(m);
  const double norm = std::sqrt(squaredNorm(m));
  if (!(det > 1e-9 * norm * norm * norm) || !std::isfinite(det))
  {
    return glm::dmat3(1.0);
  }
  for (int i = 0; i < 64; ++i)
  {
    const glm::dmat3 next = 0.5 * (m + glm::transpose(glm::inverse(m)));
    const double change = squaredNorm(next - m);
    m = next;
    if (change < 1e-28)
    {
      break;
    }
  }
  return std::isfinite(squaredNorm(m)) ? m : glm::dmat3(1.0);
}
} // namespace

SimilarityFit fitSimilarity(const std::vector<glm::dvec3>& pre, const std::vector<glm::dvec3>& post,
  const glm::dvec3& pre_normal, const glm::dvec3& post_normal)
{
  SimilarityFit fit;
  if (pre.empty() || pre.size() != post.size())
  {
    return fit;
  }
  for (size_t i = 0; i < pre.size(); ++i)
  {
    fit.pre_centroid += pre[i];
    fit.post_centroid += post[i];
  }
  fit.pre_centroid /= static_cast<double>(pre.size());
  fit.post_centroid /= static_cast<double>(post.size());

  glm::dmat3 covariance(0.0);
  double spread = 0.0;
  for (size_t i = 0; i < pre.size(); ++i)
  {
    const glm::dvec3 p = pre[i] - fit.pre_centroid;
    covariance += glm::outerProduct(post[i] - fit.post_centroid, p);
    spread += glm::dot(p, p);
  }
  if (spread <= 0.0)
  {
    return fit;
  }
  const double pre_len = glm::length(pre_normal);
  const double post_len = glm::length(post_normal);
  if (pre_len > 0.0 && post_len > 0.0)
  {
    covariance += spread * glm::outerProduct(post_normal / post_len, pre_normal / pre_len);
  }

  const glm::dmat3 rotation = polarRotation(covariance);
  double projected = 0.0;
  for (size_t i = 0; i < pre.size(); ++i)
  {
    projected += glm::dot(post[i] - fit.post_centroid, rotation * (pre[i] - fit.pre_centroid));
  }
  fit.rotation = glm::normalize(glm::quat_cast(rotation));
  fit.scale = projected > 0.0 ? projected / spread : 1.0;
  return fit;
}

glm::dvec3 blendSimilarity(const SimilarityFit& fit, const glm::dvec3& pre, const glm::dvec3& post, const double w)
{
  const glm::dvec3 offset = pre - fit.pre_centroid;
  const glm::dvec3 residual = post - (fit.post_centroid + fit.scale * (fit.rotation * offset));
  const glm::dquat rotation = glm::slerp(glm::dquat(1.0, 0.0, 0.0, 0.0), fit.rotation, w);
  return fit.pre_centroid + w * (fit.post_centroid - fit.pre_centroid)
    + (1.0 + w * (fit.scale - 1.0)) * (rotation * offset) + w * residual;
}
namespace
{
glm::dvec2 multiply(std::complex<double> z, const glm::dvec2& v)
{
  return { z.real() * v.x - z.imag() * v.y, z.imag() * v.x + z.real() * v.y };
}
} // namespace

glm::dvec2 similarityBlendMidpoint(
  const glm::dvec2& p0, const glm::dvec2& p1_in_start_frame, const glm::dvec2& p1, std::complex<double> similarity)
{
  // Half similarity on the start-frame motion; the remaining (translation / non-similar) part is blended linearly.
  return multiply(std::sqrt(similarity), 0.5 * (p0 + p1_in_start_frame))
    + 0.5 * (p1 - multiply(similarity, p1_in_start_frame));
}

std::array<glm::dvec2, 3> quadraticThroughMidpoint(const glm::dvec2& p0, const glm::dvec2& mid, const glm::dvec2& p1)
{
  return { p0, -3.0 * p0 + 4.0 * mid - p1, 2.0 * p0 - 4.0 * mid + 2.0 * p1 };
}

double quadraticSimilarityBlendDistanceError(std::complex<double> similarity)
{
  // A start-frame difference vector d becomes q(w) * d with q through 1, sqrt(z), z; exactly it scales by |z|^w.
  const std::complex<double> half = std::sqrt(similarity);
  const std::complex<double> c1 = -3.0 + 4.0 * half - similarity;
  const std::complex<double> c2 = 2.0 - 4.0 * half + 2.0 * similarity;
  double max_error = 0.0;
  for (int i = 1; i < 100; ++i)
  {
    const double w = i / 100.0;
    const double exact = std::pow(std::abs(similarity), w);
    max_error = std::max(max_error, std::abs(std::abs(1.0 + w * c1 + w * w * c2) / exact - 1.0));
  }
  return max_error;
}
} // namespace kinDS

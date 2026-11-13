#include "SimilarityBlend.hpp"

#include <cmath>

namespace kinDS
{
namespace
{
/// Rotation factor of the polar decomposition (Newton iteration); identity unless @p m is a proper rotation-scale.
glm::dmat3 polarRotation(glm::dmat3 m)
{
  const double det = glm::determinant(m);
  if (!(det > 1e-300) || !std::isfinite(det))
  {
    return glm::dmat3(1.0);
  }
  for (int i = 0; i < 64; ++i)
  {
    const glm::dmat3 next = 0.5 * (m + glm::transpose(glm::inverse(m)));
    const glm::dmat3 diff = next - m;
    m = next;
    if (glm::dot(diff[0], diff[0]) + glm::dot(diff[1], diff[1]) + glm::dot(diff[2], diff[2]) < 1e-28)
    {
      break;
    }
  }
  return m;
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
} // namespace kinDS

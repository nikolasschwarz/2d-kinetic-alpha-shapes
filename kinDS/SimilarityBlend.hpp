#pragma once

#include <glm/glm.hpp>
#include <glm/gtc/quaternion.hpp>
#include <vector>

namespace kinDS
{
/// Similarity transform (rotation, uniform scale, translation) that best maps @c pre points onto @c post points:
/// post ≈ post_centroid + scale * rotation * (pre - pre_centroid).
struct SimilarityFit
{
  glm::dvec3 pre_centroid { 0.0 };
  glm::dvec3 post_centroid { 0.0 };
  glm::dquat rotation { 1.0, 0.0, 0.0, 0.0 };
  double scale = 1.0;
};

/// Least-squares similarity fit between corresponding point sets. The plane normals are an extra rotation
/// constraint so coplanar point sets still determine a full 3D rotation. Falls back to a pure translation when
/// the rotation is undetermined.
SimilarityFit fitSimilarity(const std::vector<glm::dvec3>& pre, const std::vector<glm::dvec3>& post,
  const glm::dvec3& pre_normal, const glm::dvec3& post_normal);

/// Moves @p pre towards @p post by @p w in [0, 1]: slerped rotation, linear scale and translation, plus the
/// linear share of the point's fit residual. Exactly @p pre at w = 0 and exactly @p post at w = 1, and
/// distances between points of one fit change linearly with @p w.
glm::dvec3 blendSimilarity(const SimilarityFit& fit, const glm::dvec3& pre, const glm::dvec3& post, double w);
} // namespace kinDS

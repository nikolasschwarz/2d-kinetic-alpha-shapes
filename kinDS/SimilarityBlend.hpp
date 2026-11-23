#pragma once

#include <array>
#include <complex>
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

/// Midpoint of the polynomial (kinetic) frame blend from @p p0 (start frame) to @p p1 (end frame). @p p1_in_start_frame
/// is @p p1 expressed in the start frame and @p similarity = scale * e^(i * angle) maps start-frame vectors to
/// end-frame vectors. The midpoint applies the half similarity, so distances there follow the similarity instead of
/// shrinking as under a linear blend (by cos(angle / 2), down to mirroring through one point at 180 degrees). Without
/// rotation or scale change it is the linear midpoint.
glm::dvec2 similarityBlendMidpoint(
  const glm::dvec2& p0, const glm::dvec2& p1_in_start_frame, const glm::dvec2& p1, std::complex<double> similarity);

/// Coefficients {c0, c1, c2} of c0 + c1 * w + c2 * w^2 through @p p0 at w = 0, @p mid at w = 1/2 and @p p1 at w = 1.
std::array<glm::dvec2, 3> quadraticThroughMidpoint(const glm::dvec2& p0, const glm::dvec2& mid, const glm::dvec2& p1);

/// Largest relative error over w in [0, 1] of distances under the quadratic blend, compared with the exact similarity
/// blend (rotation by w * angle, scale^w).
double quadraticSimilarityBlendDistanceError(std::complex<double> similarity);
} // namespace kinDS

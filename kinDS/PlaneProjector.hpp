#pragma once

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <optional>
#include <utility>

namespace kinDS
{
class PlaneProjector
{
 public:
  // Construct from plane-local -> world transforms
  // Local coordinates are (u, 0, v)
  PlaneProjector(const glm::dmat4& planeAToWorld, const glm::dmat4& planeBToWorld);

  // Project local v on plane A to local return value on plane B
  glm::dvec2 project(const glm::dvec2& v) const;

  /** True when the two profile planes are parallel (shift / normal projection; no rotation axis). */
  bool isParallel() const { return m_parallel; }

  /** World-space point on the intersection line; valid only when @ref isParallel is false.
   *  WARNING: see BUG comment in the constructor — this point may not lie on either plane. */
  const glm::dvec3& intersectionPoint() const { return m_p0; }
  /** Unit direction of the intersection line; valid only when @ref isParallel is false. */
  const glm::dvec3& intersectionAxis() const { return m_axis; }
  /** Hinge angle (radians) rotating plane A onto plane B about @ref intersectionAxis; 0 when parallel. */
  double hingeAngle() const { return m_parallel ? 0.0 : m_angle; }

  /**
   * Intersection line of the two planes expressed in plane A's local (u,v) coordinates.
   * Returns nullopt when @ref isParallel (no unique intersection line).
   * @return (point on line, direction) in A-local 2D.
   */
  std::optional<std::pair<glm::dvec2, glm::dvec2>> intersectionLineInLocalA() const;

 private:
  // Extract origin + spanning vectors from transform
  void extractPlaneFromTransform(const glm::dmat4& M, glm::dvec3& origin, glm::dvec3& u, glm::dvec3& v);

  // Plane data
  glm::dvec3 m_oA, m_uA, m_vA;
  glm::dvec3 m_oB, m_uB, m_vB;

  // Normals
  glm::dvec3 m_nA, m_nB;

  // Parallel or not
  bool m_parallel;

  // --- Non-parallel case ---
  glm::dvec3 m_axis;
  double m_angle;
  glm::dmat3 m_rot;
  glm::dvec3 m_p0; // point on intersection line

  // --- Parallel case ---
  double m_dB; // plane B offset

  // Helpers
  glm::dvec3 localAToWorld(double a, double b) const;
  glm::dvec3 applyTransform(const glm::dvec3& x) const;
  glm::dvec2 worldToLocalA(const glm::dvec3& x) const;
  glm::dvec2 worldToLocalB(const glm::dvec3& x) const;
  static glm::dvec2 worldToLocalOnPlane(
    const glm::dvec3& x, const glm::dvec3& origin, const glm::dvec3& u, const glm::dvec3& v);
};
} // namespace kinDS

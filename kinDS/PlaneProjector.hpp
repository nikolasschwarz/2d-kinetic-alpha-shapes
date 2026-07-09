#pragma once

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>

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
  glm::dvec2 worldToLocalB(const glm::dvec3& x) const;
};
} // namespace kinDS

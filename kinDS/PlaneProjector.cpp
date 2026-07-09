#include "PlaneProjector.hpp"

#include <cassert>
#include <cmath>
#include "Logger.hpp"

using namespace kinDS;

static const double EPS = 1e-8;

PlaneProjector::PlaneProjector(const glm::dmat4& planeAToWorld, const glm::dmat4& planeBToWorld)
{
  // Extract plane A
  extractPlaneFromTransform(planeAToWorld, m_oA, m_uA, m_vA);

  // Extract plane B
  extractPlaneFromTransform(planeBToWorld, m_oB, m_uB, m_vB);

  // Compute normals
  m_nA = glm::normalize(glm::cross(m_uA, m_vA));
  m_nB = glm::normalize(glm::cross(m_uB, m_vB));

  // Check parallelism
  glm::dvec3 axis = glm::cross(m_nA, m_nB);
  double axisLen = glm::length(axis);

  if (axisLen < EPS)
  {
    // Parallel planes
    m_parallel = true;
    m_dB = -glm::dot(m_nB, m_oB);
  }
  else
  {
    // Non-parallel planes
    m_parallel = false;

    m_axis = axis / axisLen;

    double cosTheta = glm::clamp(glm::dot(m_nA, m_nB), -1.0, 1.0);
    m_angle = std::acos(cosTheta);

    // Rotation matrix
    m_rot = glm::dmat3(glm::rotate(glm::dmat4(1.0), m_angle, m_axis));

    // Plane offsets
    double dA = -glm::dot(m_nA, m_oA);
    double dB = -glm::dot(m_nB, m_oB);

    // Point on intersection line
    double denom = glm::dot(m_nA, glm::cross(m_nB, m_axis));

    assert(std::abs(denom) > EPS);

    m_p0 = (dB * glm::cross(m_axis, m_nA) + dA * glm::cross(m_nB, m_axis)) / denom;
  }
}

void PlaneProjector::extractPlaneFromTransform(const glm::dmat4& M, glm::dvec3& origin, glm::dvec3& u, glm::dvec3& v)
{
  // Origin
  origin = glm::dvec3(M * glm::dvec4(0.0, 0.0, 0.0, 1.0));

  // Spanning vectors (directions)
  u = glm::dvec3(M * glm::dvec4(1.0, 0.0, 0.0, 0.0));

  v = glm::dvec3(M * glm::dvec4(0.0, 0.0, 1.0, 0.0));
}

glm::dvec3 PlaneProjector::localAToWorld(double a, double b) const { return m_oA + a * m_uA + b * m_vA; }

glm::dvec3 PlaneProjector::applyTransform(const glm::dvec3& x) const
{
  if (!m_parallel)
  {
    // Rotate around intersection line
    return m_p0 + m_rot * (x - m_p0);
  }
  else
  {
    // Project along plane normal
    double t = (glm::dot(m_nB, x) + m_dB) / glm::dot(m_nB, m_nA);
    return x - t * m_nA;
  }
}

glm::dvec2 PlaneProjector::worldToLocalB(const glm::dvec3& x) const
{
  glm::dvec3 w = x - m_oB;

  double uu = glm::dot(m_uB, m_uB);
  double uv = glm::dot(m_uB, m_vB);
  double vv = glm::dot(m_vB, m_vB);

  double wu = glm::dot(w, m_uB);
  double wv = glm::dot(w, m_vB);

  double det = uu * vv - uv * uv;
  //assert(std::abs(det) > EPS);
  if(std::abs(det) < EPS)
  {
    // Degenerate case, warn
    //KINDS_WARNING("Degenerate plane in PlaneProjector");
  }

  double c = (wu * vv - wv * uv) / det;
  double d = (wv * uu - wu * uv) / det;

  return glm::dvec2(c, d);
}

glm::dvec2 PlaneProjector::project(const glm::dvec2& v) const
{
  glm::dvec3 xA = localAToWorld(v.x, v.y);
  glm::dvec3 xW = applyTransform(xA);
  return worldToLocalB(xW);
}

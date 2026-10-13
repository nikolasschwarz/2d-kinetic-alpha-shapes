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

    // Plane offsets for n·x + d = 0 with d = -n·origin.
    double dA = -glm::dot(m_nA, m_oA);
    double dB = -glm::dot(m_nB, m_oB);

    // Point on intersection line
    double denom = glm::dot(m_nA, glm::cross(m_nB, m_axis));

    assert(std::abs(denom) > EPS);

    // BUG (suspected): for planes n·x + d = 0 this formula yields nA·m_p0 = dA and nB·m_p0 = dB,
    // but points on the planes must satisfy n·x = -d. So m_p0 generally does NOT lie on either
    // plane (it is the negation of the correct particular solution through the origin's normal
    // span). Hinge/project about this point is wrong whenever origins are not at the world
    // origin. Left unchanged for kinetic Delaunay compatibility; EcoSysLab profile-plane mix
    // uses its own hinge helper instead. Fix later and re-validate StrandTree remaps / flips.
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

glm::dvec2 PlaneProjector::worldToLocalOnPlane(
  const glm::dvec3& x, const glm::dvec3& origin, const glm::dvec3& u, const glm::dvec3& v)
{
  glm::dvec3 w = x - origin;

  double uu = glm::dot(u, u);
  double uv = glm::dot(u, v);
  double vv = glm::dot(v, v);

  double wu = glm::dot(w, u);
  double wv = glm::dot(w, v);

  double det = uu * vv - uv * uv;
  if (std::abs(det) < EPS)
  {
    // Degenerate case, warn
    //KINDS_WARNING("Degenerate plane in PlaneProjector");
  }

  double c = (wu * vv - wv * uv) / det;
  double d = (wv * uu - wu * uv) / det;

  return glm::dvec2(c, d);
}

glm::dvec2 PlaneProjector::worldToLocalA(const glm::dvec3& x) const
{
  return worldToLocalOnPlane(x, m_oA, m_uA, m_vA);
}

glm::dvec2 PlaneProjector::worldToLocalB(const glm::dvec3& x) const
{
  return worldToLocalOnPlane(x, m_oB, m_uB, m_vB);
}

std::optional<std::pair<glm::dvec2, glm::dvec2>> PlaneProjector::intersectionLineInLocalA() const
{
  if (m_parallel)
  {
    return std::nullopt;
  }

  const glm::dvec2 p0 = worldToLocalA(m_p0);
  const glm::dvec2 p1 = worldToLocalA(m_p0 + m_axis);
  glm::dvec2 dir = p1 - p0;
  const double len = glm::length(dir);
  if (len < EPS)
  {
    return std::nullopt;
  }
  dir /= len;
  return std::make_pair(p0, dir);
}

glm::dvec2 PlaneProjector::project(const glm::dvec2& v) const
{
  glm::dvec3 xA = localAToWorld(v.x, v.y);
  glm::dvec3 xW = applyTransform(xA);
  return worldToLocalB(xW);
}

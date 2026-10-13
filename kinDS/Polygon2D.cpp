#include "Polygon2D.hpp"

#include <cmath>

using namespace kinDS;

namespace
{
double cross2(const glm::dvec2& a, const glm::dvec2& b, const glm::dvec2& c)
{
  return glm::cross(b - a, c - b);
}
} // namespace

std::vector<glm::dvec2> kinDS::simplifyPolygonRing(const std::vector<glm::dvec2>& polygon, double eps)
{
  std::vector<glm::dvec2> vertices;
  vertices.reserve(polygon.size());
  for (const glm::dvec2& point : polygon)
  {
    if (vertices.empty() || glm::length(point - vertices.back()) > eps)
    {
      vertices.push_back(point);
    }
  }
  if (vertices.size() > 1 && glm::length(vertices.front() - vertices.back()) <= eps)
  {
    vertices.pop_back();
  }

  bool removed_collinear = true;
  while (removed_collinear && vertices.size() > 3)
  {
    removed_collinear = false;
    for (size_t i = 0; i < vertices.size(); ++i)
    {
      const glm::dvec2& prev = vertices[(i + vertices.size() - 1) % vertices.size()];
      const glm::dvec2& current = vertices[i];
      const glm::dvec2& next = vertices[(i + 1) % vertices.size()];
      if (std::abs(cross2(prev, current, next)) <= eps)
      {
        vertices.erase(vertices.begin() + static_cast<std::ptrdiff_t>(i));
        removed_collinear = true;
        break;
      }
    }
  }

  return vertices;
}

bool kinDS::isConvexPolygon(const std::vector<glm::dvec2>& polygon, double eps)
{
  const std::vector<glm::dvec2> vertices = simplifyPolygonRing(polygon, eps);
  if (vertices.size() < 3)
  {
    return false;
  }

  int turn_sign = 0;
  const size_t vertex_count = vertices.size();
  for (size_t i = 0; i < vertex_count; ++i)
  {
    const glm::dvec2& prev = vertices[(i + vertex_count - 1) % vertex_count];
    const glm::dvec2& current = vertices[i];
    const glm::dvec2& next = vertices[(i + 1) % vertex_count];
    const double cross = cross2(prev, current, next);
    if (std::abs(cross) <= eps)
    {
      continue;
    }

    const int sign = cross > 0.0 ? 1 : -1;
    if (turn_sign == 0)
    {
      turn_sign = sign;
    }
    else if (turn_sign != sign)
    {
      return false;
    }
  }

  if (turn_sign == 0)
  {
    return false;
  }

  double signed_area2 = 0.0;
  for (size_t i = 0; i < vertex_count; ++i)
  {
    const glm::dvec2& p0 = vertices[i];
    const glm::dvec2& p1 = vertices[(i + 1) % vertex_count];
    signed_area2 += p0.x * p1.y - p1.x * p0.y;
  }

  return std::abs(signed_area2) > eps;
}

bool kinDS::isConvexPolygonOutline(const std::vector<BoundaryPoint>& outline, double eps)
{
  std::vector<glm::dvec2> polygon;
  polygon.reserve(outline.size());
  for (const BoundaryPoint& point : outline)
  {
    polygon.push_back(point.p);
  }
  return isConvexPolygon(polygon, eps);
}

#pragma once

#include <glm/geometric.hpp>
#include <cmath>
#include <stdexcept>

namespace kinDS
{
static double circumradius(const glm::dvec2& p0, const glm::dvec2& p1, const glm::dvec2& p2)
{
  const double x0 = p0[0], y0 = p0[1];
  const double x1 = p1[0], y1 = p1[1];
  const double x2 = p2[0], y2 = p2[1];

  // Side lengths
  const double a = std::hypot(x1 - x2, y1 - y2);
  const double b = std::hypot(x0 - x2, y0 - y2);
  const double c = std::hypot(x0 - x1, y0 - y1);

  // Twice the triangle area (cross product magnitude)
  const double area2 = std::abs((x1 - x0) * (y2 - y0) - (y1 - y0) * (x2 - x0));

  if (area2 == 0.0)
  {
    throw std::runtime_error("Degenerate triangle: circumradius undefined");
  }

  // R = (a * b * c) / (4 * A), and area2 = 2 * A
  return (a * b * c) / (2.0 * area2);
}
}
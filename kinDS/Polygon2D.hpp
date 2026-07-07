#pragma once

#include "KineticDelaunay.hpp"

#include <glm/glm.hpp>
#include <vector>

namespace kinDS
{
/// Remove duplicate and collinear vertices from a closed polygon ring.
std::vector<glm::dvec2> simplifyPolygonRing(const std::vector<glm::dvec2>& polygon, double eps = 1e-12);

/// True when @p polygon is a strictly convex ring with positive area (after simplification).
bool isConvexPolygon(const std::vector<glm::dvec2>& polygon, double eps = 1e-12);

/// Convenience wrapper for @ref BoundaryPoint outlines.
bool isConvexPolygonOutline(const std::vector<BoundaryPoint>& outline, double eps = 1e-12);
} // namespace kinDS

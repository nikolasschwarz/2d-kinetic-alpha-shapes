#include "kinDS/Polygon2D.hpp"

#include <catch2/catch_test_macros.hpp>

#include <glm/glm.hpp>
#include <vector>

using namespace kinDS;

namespace
{
std::vector<glm::dvec2> square(bool reversed)
{
  std::vector<glm::dvec2> polygon {
    { 0.0, 0.0 },
    { 1.0, 0.0 },
    { 1.0, 1.0 },
    { 0.0, 1.0 },
  };
  if (reversed)
  {
    std::reverse(polygon.begin(), polygon.end());
  }
  return polygon;
}
} // namespace

TEST_CASE("isConvexPolygon accepts simple convex polygons", "[polygon][convexity]")
{
  REQUIRE(isConvexPolygon(square(false)));
  REQUIRE(isConvexPolygon(square(true)));

  const std::vector<glm::dvec2> triangle {
    { 0.0, 0.0 },
    { 2.0, 0.0 },
    { 1.0, 1.5 },
  };
  REQUIRE(isConvexPolygon(triangle));
}

TEST_CASE("isConvexPolygon rejects reflex polygons", "[polygon][convexity]")
{
  const std::vector<glm::dvec2> l_shape {
    { 0.0, 0.0 },
    { 2.0, 0.0 },
    { 2.0, 1.0 },
    { 1.0, 1.0 },
    { 1.0, 2.0 },
    { 0.0, 2.0 },
  };
  REQUIRE_FALSE(isConvexPolygon(l_shape));
}

TEST_CASE("isConvexPolygon rejects degenerate polygons", "[polygon][convexity]")
{
  REQUIRE_FALSE(isConvexPolygon({}));
  REQUIRE_FALSE(isConvexPolygon({ { 0.0, 0.0 }, { 1.0, 1.0 } }));

  const std::vector<glm::dvec2> collinear {
    { 0.0, 0.0 },
    { 1.0, 0.0 },
    { 2.0, 0.0 },
    { 3.0, 0.0 },
  };
  REQUIRE_FALSE(isConvexPolygon(collinear));
}

TEST_CASE("isConvexPolygonOutline mirrors point-ring convexity", "[polygon][convexity]")
{
  const std::vector<BoundaryPoint> convex_outline {
    BoundaryPoint { 0, 0, { 0.0, 0.0 } },
    BoundaryPoint { 1, 1, { 1.0, 0.0 } },
    BoundaryPoint { 2, 2, { 0.5, 1.0 } },
  };
  REQUIRE(isConvexPolygonOutline(convex_outline));

  const std::vector<BoundaryPoint> reflex_outline {
    BoundaryPoint { 0, 0, { 0.0, 0.0 } },
    BoundaryPoint { 1, 1, { 2.0, 0.0 } },
    BoundaryPoint { 2, 2, { 2.0, 1.0 } },
    BoundaryPoint { 3, 3, { 1.0, 1.0 } },
    BoundaryPoint { 4, 4, { 1.0, 2.0 } },
    BoundaryPoint { 5, 5, { 0.0, 2.0 } },
  };
  REQUIRE_FALSE(isConvexPolygonOutline(reflex_outline));
}

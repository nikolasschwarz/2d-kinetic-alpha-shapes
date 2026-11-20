#include "kinDS/OutlineDifference.hpp"

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

using namespace kinDS;

namespace
{
// Cell of strand 0. Voronoi vertex keys use the strand and two neighbours; edges get explicit supports.
OutlineKey vv(int a, int b) { return outlineVoronoiVertexKey(0, a, b); }
OutlineKey vor(int n) { return { 10, 0, n }; }

double area(const OutlineDifferenceFace& face)
{
  double a2 = 0.0;
  for (size_t i = 0; i < face.ring.size(); ++i)
  {
    const glm::dvec2& p = face.ring[i].position;
    const glm::dvec2& q = face.ring[(i + 1) % face.ring.size()].position;
    a2 += p.x * q.y - p.y * q.x;
  }
  return 0.5 * a2;
}
} // namespace

TEST_CASE("outline difference is empty for identical outlines", "[outline_difference]")
{
  const KeyedOutline square { { vv(1, 2), vv(2, 3), vv(3, 4), vv(4, 1) }, { { 0, 0 }, { 1, 0 }, { 1, 1 }, { 0, 1 } },
    { vor(2), vor(3), vor(4), vor(1) } };
  const auto faces = outlineSymmetricDifference(square, square);
  REQUIRE(faces.has_value());
  REQUIRE(faces->faces.empty());
  REQUIRE(faces->splits.empty());
}

TEST_CASE("outline difference finds the grown part via collinear topology", "[outline_difference]")
{
  // Old: unit square; removing neighbour 3 extends the cell to x = 2. The old vertices (1,0) / (1,1) lie on the new
  // bottom / top edges (same Voronoi supports 2 / 4), so they are inserted topologically, not by intersection.
  const KeyedOutline old_outline { { vv(1, 2), vv(2, 3), vv(3, 4), vv(4, 1) },
    { { 0, 0 }, { 1, 0 }, { 1, 1 }, { 0, 1 } }, { vor(2), vor(3), vor(4), vor(1) } };
  const KeyedOutline new_outline { { vv(1, 2), vv(2, 5), vv(5, 4), vv(4, 1) },
    { { 0, 0 }, { 2, 0 }, { 2, 1 }, { 0, 1 } }, { vor(2), vor(5), vor(4), vor(1) } };
  const auto faces = outlineSymmetricDifference(old_outline, new_outline);
  REQUIRE(faces.has_value());
  REQUIRE(faces->faces.size() == 1);
  CHECK_FALSE(faces->faces.front().in_old);
  CHECK(area(faces->faces.front()) == Catch::Approx(1.0));
  CHECK(faces->faces.front().ring.size() == 4);
  // Old vertices (1,0) and (1,1) split the new bottom and top edges.
  REQUIRE(faces->splits.size() == 2);
  CHECK(faces->splits[0].outline == 1);
  CHECK(faces->splits[1].outline == 1);
}

TEST_CASE("outline difference splits at a proper crossing", "[outline_difference]")
{
  // Right side tilts: old edge x = 1, new edge from (0.5, 0) to (1.5, 1). They cross at (1, 0.5).
  const KeyedOutline old_outline { { vv(1, 2), vv(2, 3), vv(3, 4), vv(4, 1) },
    { { 0, 0 }, { 1, 0 }, { 1, 1 }, { 0, 1 } }, { vor(2), vor(3), vor(4), vor(1) } };
  const KeyedOutline new_outline { { vv(1, 2), vv(2, 5), vv(5, 4), vv(4, 1) },
    { { 0, 0 }, { 0.5, 0 }, { 1.5, 1 }, { 0, 1 } }, { vor(2), vor(5), vor(4), vor(1) } };
  const auto faces = outlineSymmetricDifference(old_outline, new_outline);
  REQUIRE(faces.has_value());
  REQUIRE(faces->faces.size() == 2);
  double old_area = 0.0;
  double new_area = 0.0;
  for (const auto& face : faces->faces)
  {
    REQUIRE(area(face) > 0.0);
    (face.in_old ? old_area : new_area) += area(face);
  }
  CHECK(old_area == Catch::Approx(0.125));
  CHECK(new_area == Catch::Approx(0.125));
}

TEST_CASE("outline difference handles an outline pinched at its site", "[outline_difference]")
{
  // Old cell touches the boundary at its site S twice (two wedges); the new cell fills the gap between them.
  const OutlineKey s = outlineSiteKey(0);
  const KeyedOutline old_outline { { s, vv(1, 2), vv(2, 3), s, vv(3, 4), vv(4, 5) },
    { { 0, 0 }, { 2, 0 }, { 1, 1 }, { 0, 0 }, { -1, 1 }, { -2, 0 } },
    { { 11, 0, 1 }, vor(2), { 11, 0, 2 }, { 11, 0, 3 }, vor(4), { 11, 0, 5 } } };
  const KeyedOutline new_outline { { s, vv(1, 2), vv(2, 3), vv(3, 4), vv(4, 5) },
    { { 0, 0 }, { 2, 0 }, { 1, 1 }, { -1, 1 }, { -2, 0 } }, { { 11, 0, 1 }, vor(2), vor(9), vor(4), { 11, 0, 5 } } };
  const auto faces = outlineSymmetricDifference(old_outline, new_outline);
  REQUIRE(faces.has_value());
  REQUIRE(faces->faces.size() == 1);
  CHECK_FALSE(faces->faces.front().in_old);
  CHECK(area(faces->faces.front()) == Catch::Approx(1.0));
}

TEST_CASE("outline difference needs two shared vertices", "[outline_difference]")
{
  const KeyedOutline a { { vv(1, 2), vv(2, 3), vv(3, 1) }, { { 0, 0 }, { 1, 0 }, { 0, 1 } },
    { vor(2), vor(3), vor(1) } };
  const KeyedOutline b { { vv(5, 6), vv(6, 7), vv(7, 5) }, { { 0, 0 }, { 2, 0 }, { 0, 2 } },
    { vor(6), vor(7), vor(5) } };
  CHECK_FALSE(outlineSymmetricDifference(a, b).has_value());
}

TEST_CASE("outline edge supports follow shared sites", "[outline_difference]")
{
  // Strand 0: Voronoi vertex (0,1,2) -> crossing of Voronoi edge (0,2) with boundary edge (0,3) -> site 0.
  const std::vector<OutlineKey> keys { outlineVoronoiVertexKey(0, 1, 2), outlineCrossingKey(3, 0, 2, 0),
    outlineSiteKey(0) };
  const auto supports = outlineEdgeSupports(keys, 0);
  REQUIRE(supports.size() == 3);
  CHECK(supports[0] == OutlineKey { 10, 0, 2 });
  CHECK(supports[1] == OutlineKey { 11, 0, 3 });
  CHECK(supports[2][0] == -1);
}

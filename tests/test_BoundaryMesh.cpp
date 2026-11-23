#include "kinDS/BoundaryMesh.hpp"

#include <catch2/catch_test_macros.hpp>

using namespace kinDS;

namespace
{
std::vector<size_t> addVertices(BoundaryMesh& mesh, size_t count)
{
  std::vector<size_t> ids;
  for (size_t i = 0; i < count; ++i)
  {
    ids.push_back(mesh.addVertex(glm::dvec3(static_cast<double>(i), 0.0, 0.0), i, 0.0, glm::dvec2(0.0), "{}"));
  }
  return ids;
}

std::vector<size_t> triangle(const BoundaryMesh& mesh, size_t index)
{
  const auto& t = mesh.mesh().getTriangles();
  return { t[3 * index], t[3 * index + 1], t[3 * index + 2] };
}
} // namespace

TEST_CASE("BoundaryMesh stores triangles with reversed winding")
{
  BoundaryMesh mesh;
  const auto v = addVertices(mesh, 3);
  mesh.addTriangle(v[0], v[1], v[2], "{}", 0, 1.0, 1.0);
  CHECK(triangle(mesh, 0) == std::vector<size_t> { v[0], v[2], v[1] });
}

TEST_CASE("BoundaryMesh closeFront fans through a pending vertex and drops it")
{
  BoundaryMesh mesh;
  mesh.resetHalfEdgeSlots(2);
  const auto v = addVertices(mesh, 5); // left, right, pending, apex, second apex
  mesh.setLastLeftAndRight(0, v[0], v[1]);
  mesh.setHalfEdgeBoundaryVertex(0, v[2]);

  CHECK(mesh.closeFront(0, v[3], "{}", 1.0, 1.0) == std::make_pair(v[0], v[1]));
  REQUIRE(mesh.triangleCount() == 2);
  CHECK(triangle(mesh, 0) == std::vector<size_t> { v[0], v[3], v[2] });
  CHECK(triangle(mesh, 1) == std::vector<size_t> { v[2], v[3], v[1] });

  mesh.closeFront(0, v[4], "{}", 1.0, 1.0);
  CHECK(mesh.triangleCount() == 3);
}

TEST_CASE("BoundaryMesh setLastLeftAndRight drops a stale pending vertex")
{
  BoundaryMesh mesh;
  mesh.resetHalfEdgeSlots(1);
  const auto v = addVertices(mesh, 4);
  mesh.setLastLeftAndRight(0, v[0], v[1]);
  mesh.setHalfEdgeBoundaryVertex(0, v[2]);
  mesh.setLastLeftAndRight(0, v[0], v[1]);

  mesh.closeFront(0, v[3], "{}", 1.0, 1.0);
  CHECK(mesh.triangleCount() == 1);
}

TEST_CASE("BoundaryMesh closeCornerGap fans through bridged corners")
{
  BoundaryMesh mesh;
  const auto v = addVertices(mesh, 6); // a, m, b, bridge apexes, gap apex
  mesh.bridgeCorners(v[1], v[0], v[3], "{}", 1.0, 1.0);
  mesh.bridgeCorners(v[2], v[1], v[4], "{}", 1.0, 1.0);
  mesh.closeCornerGap(v[0], v[2], v[5], "{}", 1.0, 1.0);
  REQUIRE(mesh.triangleCount() == 4);
  // Stored with reversed winding: (a, m, apex) -> (a, apex, m).
  CHECK(triangle(mesh, 2) == std::vector<size_t> { v[0], v[5], v[1] });
  CHECK(triangle(mesh, 3) == std::vector<size_t> { v[1], v[5], v[2] });

  mesh.closeCornerGap(v[3], v[4], v[5], "{}", 1.0, 1.0); // not bridged: one direct triangle
  CHECK(mesh.triangleCount() == 5);
}

TEST_CASE("BoundaryMesh vertexOfStrandAt finds only the strand's vertex at that time")
{
  BoundaryMesh mesh;
  const size_t a = mesh.addVertex(glm::dvec3(0.0), 7, 3.0, glm::dvec2(0.0), "{}");
  CHECK(mesh.vertexOfStrandAt(7, 3.0) == a);
  CHECK_FALSE(mesh.vertexOfStrandAt(7, 3.01).has_value());
  CHECK_FALSE(mesh.vertexOfStrandAt(8, 3.0).has_value());
}

TEST_CASE("BoundaryMesh finalize collapses coincident vertices only when an edge joins them")
{
  BoundaryMesh mesh;
  const auto p
    = [&](double x, double y) { return mesh.addVertex(glm::dvec3(x, y, 0.0), 0, 0.0, glm::dvec2(0.0), "{}"); };
  const size_t a = p(0, 0), b = p(1, 0), c = p(0, 1), a_copy = p(0, 0), d = p(-1, 0), e = p(0, -1);
  mesh.addTriangle(a, b, c, "{}", 0, 1.0, 1.0);
  mesh.addTriangle(a_copy, d, e, "{}", 0, 1.0, 1.0); // pinch copy, no edge to a
  mesh.finalize();
  CHECK(mesh.mesh().getVertexCount() == 6);

  BoundaryMesh joined;
  const auto q
    = [&](double x, double y) { return joined.addVertex(glm::dvec3(x, y, 0.0), 0, 0.0, glm::dvec2(0.0), "{}"); };
  const size_t u = q(0, 0), w = q(1, 0), u_copy = q(0, 0), z = q(0, 1);
  joined.addTriangle(u, w, z, "{}", 0, 1.0, 1.0);
  joined.addTriangle(u, u_copy, w, "{}", 0, 1.0, 1.0); // zero-area face reconnecting the copy
  joined.addTriangle(u_copy, z, w, "{}", 0, 1.0, 1.0);
  CHECK(joined.finalize() == 1);
  CHECK(joined.mesh().getVertexCount() == 3);
  CHECK(joined.triangleCount() == 2);
}

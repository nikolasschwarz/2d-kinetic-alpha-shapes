#include "eigen/Eigen/Dense"
#include "kinDS/CubicHermiteSpline.hpp"
#include "kinDS/HalfEdgeDelaunayGraphToSVG.hpp"
#include "kinDS/KineticDelaunay.hpp"
#include "kinDS/MeshIntersection.hpp"
#include "kinDS/ObjExporter.hpp"
#include "kinDS/Polynomial.hpp"
#include "kinDS/SegmentBuilder.hpp"
#include "simple_svg.hpp"
#include <iostream>

static void eigen_example()
{
  Eigen::VectorXd a(3);
  a << 1, 2, 3; // 1 + 2x + 3x^2
  Eigen::VectorXd b(2);
  b << 4, 5; // 4 + 5x

  kinDS::Polynomial p1(a), p2(b);
  kinDS::Polynomial sum = p1 + p2;
  kinDS::Polynomial prod = p1 * p2;

  sum.print(); // e.g. 3x^2 + 7x + 5
  prod.print(); // e.g. 15x^3 + 22x^2 + 13x + 4

  // Test evaluation
  double x = 2.0;
  double result_sum = sum(x);
  double result_prod = prod(x);
  std::cout << "Sum evaluated at x = " << x << ": " << result_sum << std::endl; // e.g. 3*2^2 + 7*2 + 5
  std::cout << "Product evaluated at x = " << x << ": " << result_prod << std::endl; // e.g. 15*2^3 + 22*2^2 + 13*2 + 4

  Eigen::VectorXd c(3);
  c << 4, 5, -2; // 4 + 5x - 2x^2

  auto test = POLYNOMIAL(x ^ 2); // X^2
  kinDS::Polynomial p3
    = kinDS::Polynomial([&](kinDS::Var x) { return (kinDS::Polynomial { (10 * (x ^ 4) + 4 * (x ^ 2) + 2) }); });
  p3.print();

  Eigen::VectorXcd roots = p3.roots();

  std::cout << "Roots of the polynomial: " << roots << std::endl; // Outputs the roots of the polynomial
}

#include <queue>
#include <utility> // for std::pair
#include <vector>

std::vector<std::pair<size_t, double>> merge_sorted_vectors(const std::vector<std::vector<double>>& inputs)
{
  using Entry = std::pair<size_t, double>; // (index of input vector, value)
  std::vector<std::pair<size_t, double>> result;

  struct HeapNode
  {
    size_t vec_idx; // which input vector
    size_t elem_idx; // index inside that vector
    double value; // value itself

    bool operator>(const HeapNode& other) const
    {
      return value > other.value; // for min-heap
    }
  };

  std::priority_queue<HeapNode, std::vector<HeapNode>, std::greater<>> min_heap;

  // Initialize heap with the first element of each vector
  for (size_t i = 0; i < inputs.size(); ++i)
  {
    if (!inputs[i].empty())
    {
      min_heap.push({ i, 0, inputs[i][0] });
    }
  }

  while (!min_heap.empty())
  {
    auto node = min_heap.top();
    min_heap.pop();

    // record (vector index, value)
    result.emplace_back(node.vec_idx, node.value);

    // advance in that vector
    if (node.elem_idx + 1 < inputs[node.vec_idx].size())
    {
      min_heap.push({ node.vec_idx, node.elem_idx + 1, inputs[node.vec_idx][node.elem_idx + 1] });
    }
  }

  return result;
}

static void kinetic_delaunay_example()
{
#define TEST_TRAJECTORIES
#ifdef TEST_TRAJECTORIES
  std::vector<kinDS::Point<2>> trajectory_A = {
    //{ -0.420113, -0.558875 },
    { -0.432132, -0.426942 },
    { -0.447292, -0.580708 },
    { -0.469864, -0.531837 },
    { -0.578741, -0.494280 },
    { -0.519044, -0.496727 },
    { -0.487418, -0.587100 },
    { -0.536664, -0.465019 },
  };

  std::vector<kinDS::Point<2>> trajectory_B = {
    //{ -0.200000, -0.500000 },
    { -0.150887, -0.424968 },
    { -0.101774, -0.349936 },
    { -0.052661, -0.274904 },
    { -0.003548, -0.199872 },
    { 0.045565, -0.124840 },
    { 0.094678, -0.049808 },
    { 0.143791, 0.025224 },
  };

  std::vector<kinDS::Point<2>> trajectory_C = {
    //{ 0.100000, -0.400000 },
    { -0.048665, -0.333097 },
    { -0.197330, -0.266194 },
    { -0.345995, -0.199291 },
    { -0.494660, -0.132388 },
    { -0.643325, -0.065485 },
    { -0.791990, 0.001418 },
    { -0.940656, 0.068321 },
  };

  std::vector<kinDS::Point<2>> trajectory_D = {
    //{ 0.500000, 0.200000 },
    { 0.467745, 0.111272 },
    { 0.435490, 0.022544 },
    { 0.403235, -0.066183 },
    { 0.370980, -0.154910 },
    { 0.338725, -0.243637 },
    { 0.306470, -0.332364 },
    { 0.274215, -0.421091 },
  };

  /* std::vector<kinDS::Point<2>> trajectory_B = {
    { 0.556188, -0.476838 },
    { 0.572065, -0.428676 },
    { 0.478597, -0.594223 },
    { 0.520113, -0.439088 },
    { 0.443697, -0.404739 },
    { 0.490520, -0.404842 },
    { 0.457843, -0.405508 },
    { 0.430738, -0.451790 },
  };

  std::vector<kinDS::Point<2>> trajectory_C = {
    { 0.541622, 0.548149 },
    { 0.429803, 0.496248 },
    { 0.413522, 0.498909 },
    { 0.580101, 0.499790 },
    { 0.479949, 0.562262 },
    { 0.593821, 0.428951 },
    { 0.595407, 0.407437 },
    { 0.447700, 0.444491 },
  };

  std::vector<kinDS::Point<2>> trajectory_D = {
    { -0.559656, 0.421442 },
    { -0.464952, 0.428177 },
    { -0.578349, 0.568869 },
    { -0.453764, 0.424407 },
    { -0.420543, 0.594949 },
    { -0.489771, 0.566934 },
    { -0.431853, 0.415965 },
    { -0.428599, 0.554348 },
  };*/

  kinDS::CubicHermiteSpline<2> spline_A(trajectory_A);
  kinDS::CubicHermiteSpline<2> spline_B(trajectory_B);
  kinDS::CubicHermiteSpline<2> spline_C(trajectory_C);
  kinDS::CubicHermiteSpline<2> spline_D(trajectory_D);

  std::vector<kinDS::CubicHermiteSpline<2>> splines = { spline_A, spline_B, spline_C, spline_D };

#else
  std::vector<std::vector<kinDS::Point<2>>> strand_guide_points = {
    { kinDS::Point<2> { 4.12703, -1.57022 }, kinDS::Point<2> { 4.12703, -1.57022 },
      kinDS::Point<2> { 4.12703, -1.57022 } },
    { kinDS::Point<2> { 4.3178, -3.52816 }, kinDS::Point<2> { 4.3178, -3.52816 },
      kinDS::Point<2> { 4.3178, -3.52816 } },
    { kinDS::Point<2> { -2.97276, 4.03951 }, kinDS::Point<2> { -2.97276, 4.03951 },
      kinDS::Point<2> { -2.97276, 4.03951 } },
    { kinDS::Point<2> { 3.88387, 2.44586 }, kinDS::Point<2> { 3.88387, 2.44586 },
      kinDS::Point<2> { 3.88387, 2.44586 } },
    { kinDS::Point<2> { -2.87601, 0.0555282 }, kinDS::Point<2> { -2.87601, 0.0555282 },
      kinDS::Point<2> { -2.87601, 0.0555282 } },
    { kinDS::Point<2> { -4.58249, 0.899374 }, kinDS::Point<2> { -4.58249, 0.899374 },
      kinDS::Point<2> { -4.58249, 0.899374 } },
    { kinDS::Point<2> { 2.64728, -2.65377 }, kinDS::Point<2> { 2.64728, -2.65377 },
      kinDS::Point<2> { 2.64728, -2.65377 } },
    { kinDS::Point<2> { -0.927691, -0.936193 }, kinDS::Point<2> { -0.927691, -0.936193 },
      kinDS::Point<2> { -0.927691, -0.936193 } },
    { kinDS::Point<2> { 4.04879, 0.409984 }, kinDS::Point<2> { 4.04879, 0.409984 },
      kinDS::Point<2> { 4.04879, 0.409984 } },
    { kinDS::Point<2> { 2.49795, -0.656439 }, kinDS::Point<2> { 2.49795, -0.656439 },
      kinDS::Point<2> { 2.49795, -0.656439 } },
    { kinDS::Point<2> { 0.861075, -3.77961 }, kinDS::Point<2> { 0.861075, -3.77961 },
      kinDS::Point<2> { 0.861075, -3.77961 } },
    { kinDS::Point<2> { -2.54759, -3.88424 }, kinDS::Point<2> { -2.54759, -3.88424 },
      kinDS::Point<2> { -2.54759, -3.88423 } },
    { kinDS::Point<2> { -4.39657, -0.990622 }, kinDS::Point<2> { -4.39657, -0.990622 },
      kinDS::Point<2> { -4.39657, -0.990622 } },
    { kinDS::Point<2> { -0.764416, -4.80726 }, kinDS::Point<2> { -0.764416, -4.80726 },
      kinDS::Point<2> { -0.764416, -4.80726 } },
    { kinDS::Point<2> { -1.018, 1.12591 }, kinDS::Point<2> { -1.018, 1.12591 }, kinDS::Point<2> { -1.018, 1.12591 } },
    { kinDS::Point<2> { 0.724303, -1.71089 }, kinDS::Point<2> { 0.724303, -1.71089 },
      kinDS::Point<2> { 0.724303, -1.71089 } },
    { kinDS::Point<2> { 0.429867, 2.15062 }, kinDS::Point<2> { 0.429867, 2.15062 },
      kinDS::Point<2> { 0.429867, 2.15062 } },
    { kinDS::Point<2> { 2.19038, 3.23518 }, kinDS::Point<2> { 2.19038, 3.23517 },
      kinDS::Point<2> { 2.19038, 3.23517 } },
    { kinDS::Point<2> { -1.13884, 3.0388 }, kinDS::Point<2> { -1.13884, 3.0388 },
      kinDS::Point<2> { -1.13884, 3.0388 } },
    { kinDS::Point<2> { -4.52197, 2.9438 }, kinDS::Point<2> { -4.52197, 2.9438 },
      kinDS::Point<2> { -4.52197, 2.9438 } },
    { kinDS::Point<2> { 5.8964, -0.557665 }, kinDS::Point<2> { 5.8964, -0.557665 },
      kinDS::Point<2> { 5.8964, -0.557665 } },
    { kinDS::Point<2> { 2.44269, 1.24945 }, kinDS::Point<2> { 2.44269, 1.24945 },
      kinDS::Point<2> { 2.44269, 1.24945 } },
    { kinDS::Point<2> { 2.03784, 5.19061 }, kinDS::Point<2> { 2.03784, 5.19061 },
      kinDS::Point<2> { 2.03784, 5.19061 } },
    { kinDS::Point<2> { -0.901509, -2.7655 }, kinDS::Point<2> { -0.901509, -2.7655 },
      kinDS::Point<2> { -0.901509, -2.7655 } },
    { kinDS::Point<2> { 0.289699, 4.10103 }, kinDS::Point<2> { 0.289699, 4.10103 },
      kinDS::Point<2> { 0.289699, 4.10103 } },
    { kinDS::Point<2> { 0.59763, 0.123849 }, kinDS::Point<2> { 0.59763, 0.123849 },
      kinDS::Point<2> { 0.59763, 0.123849 } },
    { kinDS::Point<2> { -4.27136, -3.04379 }, kinDS::Point<2> { -4.27136, -3.04379 },
      kinDS::Point<2> { -4.27136, -3.04379 } },
    { kinDS::Point<2> { 5.76047, 1.54084 }, kinDS::Point<2> { 5.76047, 1.54084 },
      kinDS::Point<2> { 5.76047, 1.54084 } },
    { kinDS::Point<2> { -2.88041, 1.87929 }, kinDS::Point<2> { -2.88041, 1.87929 },
      kinDS::Point<2> { -2.88041, 1.87929 } },
    { kinDS::Point<2> { -2.59, -2.00589 }, kinDS::Point<2> { -2.59, -2.00589 }, kinDS::Point<2> { -2.59, -2.00589 } }
  };

  std::vector<kinDS::CubicHermiteSpline<2>> splines;

  for (auto& points : strand_guide_points)
  {
    splines.emplace_back(points);
  }
#endif

  // For debugging, print several samples along each spline
  std::cout << "Spline evaluations:\n";
  for (double t = 0.0; t <= 2.0; t += 0.25)
  {
    std::cout << "t = " << t << ":\n";
    for (size_t i = 0; i < splines.size(); ++i)
    {
      auto p = splines[i].evaluate(t);
      std::cout << "  Spline " << i << ": (" << p[0] << ", " << p[1] << ")\n";
    }
  }

  // Test subdivisions for 4 strands
  /*std::vector<std::vector<double>> subdivisions = { { 0.42, 1.37, 2.89, 5.46 }, { 0.15, 1.92, 3.08, 4.61, 6.73 },
    { 0.08, 2.14, 2.95, 4.22, 6.11 }, { 0.33, 1.25, 2.67, 4.19, 5.78, 6.92 } };*/

  std::vector<std::vector<double>> subdivisions
    = { { 0.2, 0.4, 0.6, 0.8 }, { 0.2, 0.4, 0.6, 0.8 }, { 0.2, 0.4, 0.6, 0.8 }, { 0.2, 0.4, 0.6, 0.8 } };

  // Sort subdivisions into pairs of strand index and parameter
  std::vector<std::pair<size_t, double>> sorted_subdivisions = merge_sorted_vectors(subdivisions);

  kinDS::KineticDelaunay kinetic_delaunay(splines);

  kinetic_delaunay.init();
  kinDS::SegmentBuilder mesh_builder(kinetic_delaunay, splines, sorted_subdivisions);
  mesh_builder.init();
  auto points = kinetic_delaunay.getPointsAt(0.0);
  kinDS::HalfEdgeDelaunayGraphToSVG::write(points, kinetic_delaunay.getGraph(), "test.svg", 0.1);
  kinDS::HalfEdgeDelaunayGraphToSVG::writeVoronoi(points, kinetic_delaunay.getGraph(), "test_voronoi.svg", 0.1);

  size_t section_count = 1;
  // kinetic_delaunay.getSectionCount();

  for (size_t i = 0; i < section_count; ++i)
  {
    if (i != 0)
      mesh_builder.betweenSections(i);
    kinetic_delaunay.advanceOneSection(mesh_builder);
    // kinetic_delaunay.getGraph().printDebug();
    points = kinetic_delaunay.getPointsAt(static_cast<double>(i + 1));
    kinDS::HalfEdgeDelaunayGraphToSVG::write(
      points, kinetic_delaunay.getGraph(), "test_" + std::to_string(i + 1) + ".svg", 0.1);
    std::cout << "Wrote " << ("test_" + std::to_string(i + 1) + ".svg") << std::endl;
    kinDS::HalfEdgeDelaunayGraphToSVG::writeVoronoi(
      points, kinetic_delaunay.getGraph(), "test_voronoi_" + std::to_string(i + 1) + ".svg", 0.1);
    std::cout << "Wrote " << ("test_voronoi_" + std::to_string(i + 1) + ".svg") << std::endl;
  }

  mesh_builder.finalize(section_count);

  // mesh_builder.printDebugInfo();

  // auto meshes = mesh_builder.extractMeshes();
  auto meshes = mesh_builder.extractSegmentMeshlets();
  //(0.1, 0.01, subdivisions);

  for (size_t i = 0; i < meshes.size(); ++i)
  {
    std::string filename = "mesh_" + std::to_string(i) + ".obj";
    kinDS::ObjExporter::writeMesh(meshes[i], filename);
    std::cout << "Mesh saved to " << filename << std::endl;
  }

  auto& boundary_mesh = mesh_builder.getBoundaryMesh();
  kinDS::ObjExporter::writeMesh(boundary_mesh, "boundary_mesh.obj");

  boundary_mesh.checkForDegenerateTriangles();

  // intersect all meshes with the boundary mesh and save the result
  /* for (size_t i = 0; i < meshes.size(); ++i)
  {
    meshes[i].checkForDegenerateTriangles();
    auto intersection = kinDS::MeshIntersection::intersect(meshes[i], boundary_mesh);
    std::string filename = "intersection_" + std::to_string(i) + ".obj";
    kinDS::ObjExporter::writeMesh(intersection, filename);
    std::cout << "Intersection mesh saved to " << filename << std::endl;
  }*/
}

int main()
{
  // eigen_example();

  // voronoi_example();

  kinetic_delaunay_example();

  // kinDS::MeshIntersection::test();
}

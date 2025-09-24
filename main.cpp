#include "eigen/Eigen/Dense"
#include "kinDS/CubicHermiteSpline.hpp"
#include "kinDS/HalfEdgeDelaunayGraphToSVG.hpp"
#include "kinDS/KineticDelaunay.hpp"
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
  kinDS::Polynomial p3 = kinDS::Polynomial([&](kinDS::Var x)
    { return (kinDS::Polynomial { (10 * (x ^ 4) + 4 * (x ^ 2) + 2) }); });
  p3.print();

  Eigen::VectorXcd roots = p3.roots();

  std::cout << "Roots of the polynomial: " << roots << std::endl; // Outputs the roots of the polynomial
}

#include <queue>
#include <utility> // for std::pair
#include <vector>

std::vector<std::pair<size_t, double>> merge_sorted_vectors(
  const std::vector<std::vector<double>>& inputs)
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
      min_heap.push({ node.vec_idx,
        node.elem_idx + 1,
        inputs[node.vec_idx][node.elem_idx + 1] });
    }
  }

  return result;
}

static void kinetic_delaunay_example()
{
  std::vector<kinDS::Point<2>> trajectory_A = {
    { -0.420113, -0.558875 },
    { -0.432132, -0.426942 },
    { -0.447292, -0.580708 },
    { -0.469864, -0.531837 },
    { -0.578741, -0.494280 },
    { -0.519044, -0.496727 },
    { -0.487418, -0.587100 },
    { -0.536664, -0.465019 },
  };

  std::vector<kinDS::Point<2>> trajectory_B = {
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
  };

  kinDS::CubicHermiteSpline<2> spline_A(trajectory_A);
  kinDS::CubicHermiteSpline<2> spline_B(trajectory_B);
  kinDS::CubicHermiteSpline<2> spline_C(trajectory_C);
  kinDS::CubicHermiteSpline<2> spline_D(trajectory_D);

  std::vector<kinDS::CubicHermiteSpline<2>> splines = {
    spline_A,
    spline_B,
    spline_C,
    spline_D
  };

  // Test subdivisions for 4 strands
  std::vector<std::vector<double>> subdivisions = {
    { 0.42, 1.37, 2.89, 5.46 },
    { 0.15, 1.92, 3.08, 4.61, 6.73 },
    { 0.08, 2.14, 2.95, 4.22, 6.11 },
    { 0.33, 1.25, 2.67, 4.19, 5.78, 6.92 }
  };

  // Sort subdivisions into pairs of strand index and parameter
  std::vector<std::pair<size_t, double>> sorted_subdivisions = merge_sorted_vectors(subdivisions);

  kinDS::KineticDelaunay kinetic_delaunay(splines);

  kinetic_delaunay.init();
  kinDS::SegmentBuilder mesh_builder(kinetic_delaunay, splines, sorted_subdivisions);
  mesh_builder.init();
  auto points = kinetic_delaunay.getPointsAt(0.0);
  kinDS::HalfEdgeDelaunayGraphToSVG::write(points, kinetic_delaunay.getGraph(), "test.svg", 0.1);
  kinDS::HalfEdgeDelaunayGraphToSVG::writeVoronoi(points, kinetic_delaunay.getGraph(), "test_voronoi.svg", 0.1);

  size_t section_count = kinetic_delaunay.getSectionCount();

  for (size_t i = 0; i < section_count; ++i)
  {
    if (i != 0)
      mesh_builder.betweenSections(i);
    kinetic_delaunay.advanceOneSection(mesh_builder);
    // kinetic_delaunay.getGraph().printDebug();
    points = kinetic_delaunay.getPointsAt(static_cast<double>(i + 1));
    kinDS::HalfEdgeDelaunayGraphToSVG::write(points, kinetic_delaunay.getGraph(), "test_" + std::to_string(i + 1) + ".svg", 0.1);
    kinDS::HalfEdgeDelaunayGraphToSVG::writeVoronoi(points, kinetic_delaunay.getGraph(), "test_voronoi_" + std::to_string(i + 1) + ".svg", 0.1);
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
}

int main()
{
  // eigen_example();

  // voronoi_example();

  kinetic_delaunay_example();
}

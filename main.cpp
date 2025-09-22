#include "eigen/Eigen/Dense"
#include "kinDS/CubicHermiteSpline.hpp"
#include "kinDS/HalfEdgeDelaunayGraphToSVG.hpp"
#include "kinDS/KineticDelaunay.hpp"
#include "kinDS/MeshBuilder.hpp"
#include "kinDS/ObjExporter.hpp"
#include "kinDS/Polynomial.hpp"
#include "kinDS/SegmentBuilder.hpp"
#include "simple_svg.hpp"
#include "voronoi/VoronoiDiagramGenerator.h"
#include <iostream>

static void voronoi_example()
{
  const size_t count = 4;
  float xValues[count] = { -22, -17, 4, 22 };
  float yValues[count] = { -9, 31, 13, -5 };

  float minX = -100, maxX = 100, minY = -100, maxY = 100;

  voronoi_diagram_generator::VoronoiDiagramGenerator vdg;
  vdg.generateVoronoi(xValues, yValues, count, -100, 100, -100, 100, 3);

  vdg.resetIterator();

  float x1, y1, x2, y2;

  printf("\n-------------------------------\n");

  {
    using namespace svg;
    float width = 3000.0f;
    float scale = width / (maxX - minX);
    Dimensions dimensions(maxX - minX, maxY - minY);
    Layout layout(dimensions, Layout::TopLeft, scale, svg::Point(-minX, -minY));

    Document doc("voronoi.svg", layout);

    // output lines of the Voronoi diagram
    while (vdg.getNext(x1, y1, x2, y2))
    {
      printf("GOT Line (%f,%f)->(%f,%f)\n", x1, y1, x2, y2);
      svg::Line line(
        svg::Point(x1, y1), svg::Point(x2, y2),
        svg::Stroke(0.5, svg::Color::Black));

      doc << line;
    }

    // output the original sites
    for (long i = 0; i < count; ++i)
    {
      printf("GOT Site (%f,%f)\n", xValues[i], yValues[i]);
      doc << svg::Circle(svg::Point(xValues[i], yValues[i]), 3, svg::Fill(svg::Color::Red), svg::Stroke(1, svg::Color::Black));
    }

    // Save document and report status
    if (doc.save())
    {
      std::cout << "File saved successfully." << std::endl;
    }
    else
    {
      std::cout << "Failed to save the file." << std::endl;
    }
  }
}

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

  kinDS::KineticDelaunay kinetic_delaunay(splines);

  kinetic_delaunay.init();
  kinDS::SegmentBuilder mesh_builder(kinetic_delaunay, splines);
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

  // Test subdivisions for 4 strands
  std::vector<std::vector<double>> subdivisions = {
    { 0.42, 1.37, 2.89, 5.46 },
    { 0.15, 1.92, 3.08, 4.61, 6.73 },
    { 0.08, 2.14, 2.95, 4.22, 6.11 },
    { 0.33, 1.25, 2.67, 4.19, 5.78, 6.92 }
  };

  auto meshes = mesh_builder.extractMeshes();
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

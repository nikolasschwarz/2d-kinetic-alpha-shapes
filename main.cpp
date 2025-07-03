#include "eigen/Eigen/Dense"
#include "kinDS/AVLTree.hpp"
#include "kinDS/CubicHermiteSpline.hpp"
#include "kinDS/HalfEdgeDelaunayGraphToSVG.hpp"
#include "kinDS/KineticDelaunay.hpp"
#include "kinDS/Polynomial.hpp"
#include "simple_svg.hpp"
#include "voronoi/VoronoiDiagramGenerator.h"
#include <iostream>

void voronoi_example()
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

        /*Polygon test_cell(Stroke(1, Color::Red));

        test_cell << svg::Point(0, 0) << svg::Point(layout.dimensions.width, 0)
                << svg::Point(layout.dimensions.width, layout.dimensions.height)
                << svg::Point(0, layout.dimensions.height);
        doc << test_cell;

        doc << Circle(svg::Point(0, 0), 20, Fill(Color::Red), Stroke(1, Color::Black));
        doc << Circle(
                svg::Point(layout.dimensions.width - 10, layout.dimensions.height - 10), 20,
                Fill(Color::Red), Stroke(1, Color::Black));
                */

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

void eigen_example()
{

    /*using Eigen::MatrixXd;
    MatrixXd m(2, 2);
    m(0, 0) = 3;
    m(1, 0) = 2.5;
    m(0, 1) = -1;
    m(1, 1) = m(1, 0) + m(0, 1);
    std::cout << m << std::endl;*/
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
    double result_sum = sum.eval(x);
    double result_prod = prod.eval(x);
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

void kinetic_delaunay_example()
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

    kinDS::KineticDelaunay kinetic_delaunay({ spline_A,
        spline_B,
        spline_C,
        spline_D });

    kinetic_delaunay.init();
    auto points = kinetic_delaunay.getPointsAt(0.0);
    HalfEdgeDelaunayGraphToSVG::write(points, kinetic_delaunay.getGraph(), "test.svg");

    size_t section_count = kinetic_delaunay.getSectionCount();

    for (size_t i = 0; i < section_count; ++i)
    {
        kinetic_delaunay.advanceOneSection();
        points = kinetic_delaunay.getPointsAt(static_cast<double>(i + 1));
        HalfEdgeDelaunayGraphToSVG::write(points, kinetic_delaunay.getGraph(), "test_" + std::to_string(i + 1) + ".svg");
    }
}

void avl_tree_example()
{
    kinDS::AVLTree<int, int> test;

    test.insert(1, 1);
    test.insert(3, 3);
    test.insert(4, 4);
    test.insert(5, 5);
    test.insert(2, 2);

    test.printTreeStructure();
}

int main()
{
    // eigen_example();

    // voronoi_example();

    kinetic_delaunay_example();
}

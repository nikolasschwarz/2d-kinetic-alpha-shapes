#include <iostream>
#include "kinDS/AVLTree.hpp"
#include "voronoi/VoronoiDiagramGenerator.h"
#include "simple_svg.hpp"
#include "eigen/Eigen/Dense"
#include "Polynomial.hpp"

void voronoi_example()
{
	const size_t count = 4;
	float xValues[count] = { -22, -17, 4, 22};
	float yValues[count] = { -9, 31,13,-5};

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
	Eigen::VectorXd a(3); a << 1, 2, 3;  // 1 + 2x + 3x^2
	Eigen::VectorXd b(2); b << 4, 5;     // 4 + 5x

	Polynomial p1(a), p2(b);
	Polynomial sum = p1 + p2;
	Polynomial prod = p1 * p2;

	sum.print();   // e.g. 3x^2 + 7x + 5
	prod.print();  // e.g. 15x^3 + 22x^2 + 13x + 4

	// Test evaluation
  double x = 2.0;
  double result_sum = sum.eval(x);
  double result_prod = prod.eval(x);
  std::cout << "Sum evaluated at x = " << x << ": " << result_sum << std::endl;   // e.g. 3*2^2 + 7*2 + 5
  std::cout << "Product evaluated at x = " << x << ": " << result_prod << std::endl; // e.g. 15*2^3 + 22*2^2 + 13*2 + 4
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
	eigen_example();

	//voronoi_example();
}

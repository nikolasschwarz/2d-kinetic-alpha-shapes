#include <iostream>
#include "kinDS/AVLTree.hpp"
#include "voronoi/VoronoiDiagramGenerator.h"

void voronoi_example()
{
	float xValues[4] = { -22, -17, 4,22 };
	float yValues[4] = { -9, 31,13,-5 };

	long count = 4;

	VoronoiDiagramGenerator vdg;
	vdg.generateVoronoi(xValues, yValues, count, -100, 100, -100, 100, 3);

	vdg.resetIterator();

	float x1, y1, x2, y2;

	printf("\n-------------------------------\n");
	while (vdg.getNext(x1, y1, x2, y2))
	{
		printf("GOT Line (%f,%f)->(%f,%f)\n", x1, y1, x2, y2);

	}
}

int main()
{
    kinDS::AVLTree<int, int> test;

    test.insert(1, 1);
    test.insert(3, 3);
    test.insert(4, 4);
    test.insert(5, 5);
    test.insert(2, 2);
    
    test.printTreeStructure();

		voronoi_example();
}

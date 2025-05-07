#include <iostream>
#include "kinDS/AVLTree.hpp"
#include "voronoi/VoronoiDiagram2D.hpp"

void voronoi_example()
{
  std::vector<Point2D> points = {
    {0.0, 0.0},
    {1.1, 0.0},
    {0.5, 1.0},
    {2.0, 2.0},
    {3.0, 1.5}
  };

  VoronoiDiagram2D voronoi(points);
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

#include "VoronoiDiagram2D.hpp"

#include "libqhullcpp/Qhull.h"

using namespace orgQhull;


VoronoiDiagram2D::VoronoiDiagram2D(const std::vector<Point2D>& inputPoints) {
  computeVoronoi(inputPoints);
}

VoronoiDiagram2D::~VoronoiDiagram2D() {
  cleanup();
}

void VoronoiDiagram2D::cleanup() {
 
}

void VoronoiDiagram2D::computeVoronoi(const std::vector<Point2D>& inputPoints) {
  cleanup();  // in case re-used

  Qhull qh;
}

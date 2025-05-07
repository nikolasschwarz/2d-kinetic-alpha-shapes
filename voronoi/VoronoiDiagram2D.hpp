#pragma once

#include <array>
#include <vector>
#include "qhull_ra.h"

using Point2D = std::array<double, 2>;

class VoronoiDiagram2D {
public:
  

  struct Cell {
    std::vector<Point2D> vertices;  // Convex polygon
    std::vector<int> neighbors;     // Adjacent input sites
  };

  explicit VoronoiDiagram2D(const std::vector<Point2D>& inputPoints);
  ~VoronoiDiagram2D();

  const std::vector<Cell>& cells() const {
    return _cells;
  }

private:
  std::vector<Cell> _cells;
  qhT _qh;  // Reentrant Qhull context
  bool _qhInitialized = false;

  void computeVoronoi(const std::vector<Point2D>& inputPoints);
  void cleanup();
};

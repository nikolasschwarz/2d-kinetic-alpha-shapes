#include "VoronoiDiagram2D.hpp"
#include <cstring>  // for memset
#include <stdexcept>
#include <iostream>

VoronoiDiagram2D::VoronoiDiagram2D(const std::vector<Point2D>& inputPoints) {
  computeVoronoi(inputPoints);
}

VoronoiDiagram2D::~VoronoiDiagram2D() {
  cleanup();
}

void VoronoiDiagram2D::cleanup() {
  if (_qhInitialized) {
    qh_freeqhull(&_qh, !qh_ALL);  // Free memory
    int curlong, totlong;
    qh_memfreeshort(&_qh, &curlong, &totlong);
    _qhInitialized = false;
    if (curlong || totlong) {
      std::cerr << "Qhull: memory not fully freed: "
        << curlong << " blocks remaining, "
        << totlong << " bytes in total\n";
    }
  }
}

void VoronoiDiagram2D::computeVoronoi(const std::vector<Point2D>& inputPoints) {
  cleanup();  // in case re-used

  if (inputPoints.size() < 3)
    throw std::runtime_error("At least 3 points are needed for a Voronoi diagram");

  int dim = 2;
  int numPoints = static_cast<int>(inputPoints.size());
  std::vector<coordT> points(dim * numPoints);
  for (int i = 0; i < numPoints; ++i) {
    points[2 * i] = inputPoints[i][0];
    points[2 * i + 1] = inputPoints[i][1];
  }

  // Initialize Qhull
  memset(&_qh, 0, sizeof(qhT));
  qh_init_A(&_qh, nullptr, stdout, stderr, 0, nullptr);  // output to stdout/stderr
  _qhInitialized = true;

  char* options = (char*)"qhull v Qbb Qc";  // v=Voronoi, Qbb=scale, Qc=center
  int exitcode = qh_new_qhull(&_qh, dim, numPoints, points.data(), false, options, nullptr, nullptr);

  if (exitcode != 0)
    throw std::runtime_error("Qhull failed to compute the Voronoi diagram");

  // Extract Voronoi vertices and cells
  _cells.clear();
  _cells.resize(numPoints);

  // TODO: Handle unbounded cells
  for (facetT* facet = _qh.facet_list; facet && facet->next; facet = facet->next) {
    if (!facet->upperdelaunay && facet->center) {
      for (int i = 0; i < facet->vertices->maxsize; ++i) {

        for (vertexT* vertex = _qh.vertex_list; vertex && vertex->next; vertex = vertex->next) {
          setT* neighbors = vertex->neighbors;
          int neighborCount = neighbors->maxsize;

          for (int i = 0; i < neighborCount; ++i) {
            vertexT* neighbor = (vertexT*)SETelem_(neighbors, i);
            if (!neighbor)
              continue;  // Skip NULL entries
            // Get associated Delaunay input site for this facet
            int siteIndex = qh_pointid(&_qh, neighbor->point);
            if (siteIndex >= 0 && siteIndex < numPoints) {
              auto& cell = _cells[siteIndex];
              coordT* center = facet->center;
              cell.vertices.push_back({ center[0], center[1] });
              // optionally collect neighbors here
              cell.neighbors.push_back(siteIndex);
            }
          }
          break;  // Only need to process each facet once
        }
      }
    }
  }
}

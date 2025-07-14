#pragma once

#include "../simple_svg.hpp"
#include "HalfEdgeDelaunayGraph.hpp"
#include "Logger.hpp"

namespace kinDS
{

class HalfEdgeDelaunayGraphToSVG
{
 public:
  HalfEdgeDelaunayGraphToSVG() = delete;
  ~HalfEdgeDelaunayGraphToSVG() = delete;

  struct BoundingBox
  {
    double min_x;
    double min_y;
    double max_x;
    double max_y;
    BoundingBox(double min_x, double min_y, double max_x, double max_y)
      : min_x(min_x)
      , min_y(min_y)
      , max_x(max_x)
      , max_y(max_y)
    {
    }
  };

  static BoundingBox computeBoundingBox(const std::vector<Point<2>>& points, double margin = 0.0)
  {
    double min_x = std::numeric_limits<double>::max();
    double min_y = std::numeric_limits<double>::max();
    double max_x = std::numeric_limits<double>::lowest();
    double max_y = std::numeric_limits<double>::lowest();
    for (const auto& point : points)
    {
      if (point[0] < min_x)
        min_x = point[0];
      if (point[1] < min_y)
        min_y = point[1];
      if (point[0] > max_x)
        max_x = point[0];
      if (point[1] > max_y)
        max_y = point[1];
    }
    // Apply margin to the bounding box
    min_x -= margin;
    min_y -= margin;
    max_x += margin;
    max_y += margin;
    return BoundingBox(min_x, min_y, max_x, max_y);
  }

  static svg::Document setupDocument(const std::vector<Point<2>> points, const std::string& filename, const BoundingBox& bb)
  {
    // Create an SVG document with the bounding box
    double width = bb.max_x - bb.min_x;
    double height = bb.max_y - bb.min_y;

    svg::Dimensions dimensions(width, height);

    return svg::Document(filename, svg::Layout(dimensions, svg::Layout::TopLeft, 1.0, svg::Point(-bb.min_x, -bb.min_y)));
  }

  /**
   * @brief Converts the half-edge Delaunay graph to an SVG representation.
   *
   * @param graph The half-edge Delaunay graph to convert.
   * @param filename The name of the output SVG file.
   */
  static void write(const std::vector<Point<2>> points, const HalfEdgeDelaunayGraph& graph, const std::string& filename, double margin = 0.0)
  {
    BoundingBox bb = computeBoundingBox(points, margin);
    svg::Document doc = setupDocument(points, filename, bb);
    // Draw edges
    for (size_t he_id = 0; he_id < graph.getHalfEdges().size(); he_id += 2)
    {
      const HalfEdgeDelaunayGraph::HalfEdge& he = graph.getHalfEdges()[he_id];
      if (he.origin != -1 && graph.getHalfEdges()[he_id ^ 1].origin != -1) // Only draw edges that are not boundary edges
      {
        Point<2> start = points[graph.getHalfEdges()[he_id].origin];
        Point<2> end = points[graph.getHalfEdges()[he_id ^ 1].origin];
        doc << svg::Line(svg::Point(start[0], start[1]), svg::Point(end[0], end[1]),
          svg::Stroke(0.01, svg::Color::Black));
      }
    }

    // Draw vertices
    for (auto& point : points)
    {
      doc << svg::Circle(svg::Point(point[0], point[1]), 0.02, svg::Fill(svg::Color::Blue), svg::Stroke(0.0, svg::Color::Black));
    }
    doc.save();
  }

  static void writeVoronoi(const std::vector<Point<2>> points, const HalfEdgeDelaunayGraph& graph, const std::string& filename, double margin = 0.0)
  {
    auto circumcenters = graph.computeCircumcenters(points);

    BoundingBox bb = computeBoundingBox(points, margin);
    svg::Document doc = setupDocument(points, filename, bb);

    // Draw vertices
    for (auto& point : points)
    {
      doc << svg::Circle(svg::Point(point[0], point[1]), 0.02, svg::Fill(svg::Color::Blue), svg::Stroke(0.0, svg::Color::Black));
    }

    // Draw Voronoi edges
    for (size_t he_id = 0; he_id < graph.getHalfEdges().size(); he_id += 2)
    {
      const HalfEdgeDelaunayGraph::HalfEdge& he = graph.getHalfEdges()[he_id];
      if (he.origin != -1 && graph.getHalfEdges()[he_id ^ 1].origin != -1) // Only draw edges that are not boundary edges
      {
        auto start = circumcenters[graph.getHalfEdges()[he_id].face];
        auto end = circumcenters[graph.getHalfEdges()[he_id ^ 1].face];

        if (!start.second && !end.second) // Only draw finite lines
        {

          doc << svg::Line(svg::Point(start.first[0], start.first[1]), svg::Point(end.first[0], end.first[1]),
            svg::Stroke(0.01, svg::Color::Black));
        }
        else if (start.second && !end.second)
        {
          // Draw a line to infinity in the direction of the circumcenter
          svg::Point end_point(end.first[0], end.first[1]);
          svg::Point start_point(end.first[0] + 1000 * start.first[0], end.first[1] + 1000 * start.first[1]);
          doc << svg::Line(start_point, end_point, svg::Stroke(0.01, svg::Color::Black));
        }
        else if (!start.second && end.second)
        {
          // Draw a line to infinity in the direction of the circumcenter
          svg::Point start_point(start.first[0], start.first[1]);
          svg::Point end_point(start.first[0] + 1000 * end.first[0], start.first[1] + 1000 * end.first[1]);
          doc << svg::Line(start_point, end_point, svg::Stroke(0.01, svg::Color::Black));
        }
        else
        {
          // This should only happen if all points lie on one line, maybe we should implement it later
          logger.log(WARNING, "Both circumcenters are infinite, skipping edge.");
        }

        // TODO:: Handle infinite edges
      }
    }

    // Draw Voronoi vertices
    for (const auto& circumcenter : circumcenters)
    {
      if (!circumcenter.second) // Only draw finite vertices
      {
        doc << svg::Circle(svg::Point(circumcenter.first[0], circumcenter.first[1]), 0.02, svg::Fill(svg::Color::Red), svg::Stroke(0.0, svg::Color::Black));
      }
    }

    doc.save();
  }
};
} // namespace kinDS

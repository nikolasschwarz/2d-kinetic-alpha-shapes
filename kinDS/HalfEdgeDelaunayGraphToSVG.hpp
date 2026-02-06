#pragma once

#include "HalfEdgeDelaunayGraph.hpp"
#include "Logger.hpp"
#include "simple_svg.hpp"
#include <glm/glm.hpp>

namespace kinDS
{

// use this for font placement
static glm::dvec2 triangleIncenter(const glm::dvec2& A, const glm::dvec2& B, const glm::dvec2& C)
{
  const double a = glm::length(B - C); // opposite A
  const double b = glm::length(C - A); // opposite B
  const double c = glm::length(A - B); // opposite C

  const double sum = a + b + c;

  // Degenerate triangle guard (optional but recommended)
  if (sum == 0.0)
  {
    return A; // or any reasonable fallback
  }

  return (a * A + b * B + c * C) / sum;
}

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

  static BoundingBox computeBoundingBox(const std::vector<glm::dvec2>& points, double margin = 0.0)
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

  static svg::Document setupDocument(
    const std::vector<glm::dvec2> points, const std::string& filename, const BoundingBox& bb)
  {
    // Create an SVG document with the bounding box
    double width = bb.max_x - bb.min_x;
    double height = bb.max_y - bb.min_y;

    svg::Dimensions dimensions(width, height);

    return svg::Document(
      filename, svg::Layout(dimensions, svg::Layout::TopLeft, 1.0, svg::Point(-bb.min_x, -bb.min_y)));
  }

  /**
   * @brief Converts the half-edge Delaunay graph to an SVG representation.
   *
   * @param graph The half-edge Delaunay graph to convert.
   * @param filename The name of the output SVG file.
   */
  static void write(const std::vector<glm::dvec2> points, const HalfEdgeDelaunayGraph& graph,
    const std::string& filename, double margin = 0.0, std::vector<bool>* face_inside = nullptr)
  {
    BoundingBox bb = computeBoundingBox(points, margin);
    svg::Document doc = setupDocument(points, filename, bb);

    // Draw faces
    for (size_t face_id = 0; face_id < graph.getFaces().size(); face_id++)
    {
      auto face_vertex_indices = graph.getTriangleVertexIndices(face_id);

      // If any vertex is infinite, skip the face
      if (face_vertex_indices[0] == -1 || face_vertex_indices[1] == -1 || face_vertex_indices[2] == -1)
      {
        // assert that face is outside
        if (face_inside)
        {
          assert(!(*face_inside)[face_id] && "Face is inside despite being infinite!");
        }
        continue;
      }

      svg::Color face_color { svg::Color::Green };

      if (face_inside && !(*face_inside)[face_id])
      {
        face_color = svg::Color { svg::Color::Red };
      }

      std::array<glm::dvec2, 3> face_vertices
        = { points[face_vertex_indices[0]], points[face_vertex_indices[1]], points[face_vertex_indices[2]] };

      svg::Polygon face { svg::Fill(face_color) };
      face << svg::Point(face_vertices[0][0], face_vertices[0][1])
           << svg::Point(face_vertices[1][0], face_vertices[1][1])
           << svg::Point(face_vertices[2][0], face_vertices[2][1]);
      doc << face;

      // Draw face id at incenter
      glm::dvec2 incenter = triangleIncenter(face_vertices[0], face_vertices[1], face_vertices[2]);
      doc << svg::Text(
        svg::Point(incenter[0], incenter[1]), std::to_string(face_id), svg::Fill(svg::Color::White), svg::Font(0.01));
    }

    // Draw edges
    for (size_t he_id = 0; he_id < graph.getHalfEdges().size(); he_id += 2)
    {
      const HalfEdgeDelaunayGraph::HalfEdge& he = graph.getHalfEdges()[he_id];
      if (he.origin != -1
        && graph.getHalfEdges()[he_id ^ 1].origin != -1) // Only draw edges that are not boundary edges
      {
        glm::dvec2 start = points[graph.getHalfEdges()[he_id].origin];
        glm::dvec2 end = points[graph.getHalfEdges()[he_id ^ 1].origin];
        doc << svg::Line(
          svg::Point(start[0], start[1]), svg::Point(end[0], end[1]), svg::Stroke(0.01, svg::Color::Black));
      }
    }

    // Draw vertices
    for (size_t v = 0; v < points.size(); v++)
    {
      auto& point = points[v];
      doc << svg::Circle(
        svg::Point(point[0], point[1]), 0.02, svg::Fill(svg::Color::Blue), svg::Stroke(0.0, svg::Color::Black));

      // Draw vertex id
      doc << svg::Text(svg::Point(point[0] - 0.005, point[1] - 0.005), std::to_string(v), svg::Fill(svg::Color::White),
        svg::Font(0.01));
    }

    for (size_t he_id = 0; he_id < graph.getHalfEdges().size(); he_id += 2)
    {
      const HalfEdgeDelaunayGraph::HalfEdge& he = graph.getHalfEdges()[he_id];
      // Draw half-edge id at midpoint but slightly offset to the left in the direction of the edge normal

      if (he.origin != -1 && graph.getHalfEdges()[he_id ^ 1].origin != -1)
      {
        // Do this for both half-edges
        for (size_t i = 0; i < 2; i++)
        {

          glm::dvec2 start = points[graph.getHalfEdges()[he_id].origin];
          glm::dvec2 end = points[graph.getHalfEdges()[he_id ^ 1].origin];
          glm::dvec2 midpoint = (start + end) / 2.0;
          glm::dvec2 edge_dir = glm::normalize(end - start);
          glm::dvec2 edge_normal(-edge_dir[1], edge_dir[0]); // Rotate 90 degrees to get normal
          glm::dvec2 label_pos
            = midpoint + std::pow(-1, i) * 0.01 * edge_normal - glm::dvec2(0.005, 0.005); // Offset by 0.02 units
          doc << svg::Text(svg::Point(label_pos[0], label_pos[1]), std::to_string(he_id + i),
            svg::Fill(svg::Color::Yellow), svg::Font(0.01));
        }
      }
    }
    doc.save();
  }

  static void writeVoronoi(const std::vector<glm::dvec2> points, const HalfEdgeDelaunayGraph& graph,
    const std::string& filename, bool also_draw_delaunay, double margin = 0.0)
  {
    auto circumcenters = graph.computeCircumcenters(points);

    const std::vector<glm::dvec2> allFinitePoints = [&]()
    {
      std::vector<glm::dvec2> finitePoints;

      std::copy(points.begin(), points.end(), std::back_inserter(finitePoints));

      for (const auto& cc : circumcenters)
      {
        if (!cc.second) // only finite points
        {
          finitePoints.push_back(cc.first);
        }
      }
      return finitePoints;
    }();

    BoundingBox bb = computeBoundingBox(allFinitePoints, margin);
    svg::Document doc = setupDocument(points, filename, bb);

    // Draw vertices
    for (auto& point : points)
    {
      doc << svg::Circle(
        svg::Point(point[0], point[1]), 0.02, svg::Fill(svg::Color::Blue), svg::Stroke(0.0, svg::Color::Black));
    }

    // Draw Voronoi edges
    for (size_t he_id = 0; he_id < graph.getHalfEdges().size(); he_id += 2)
    {
      const HalfEdgeDelaunayGraph::HalfEdge& he = graph.getHalfEdges()[he_id];
      if (he.origin != -1
        && graph.getHalfEdges()[he_id ^ 1].origin != -1) // Only draw edges that are not boundary edges
      {
        auto start = circumcenters[graph.getHalfEdges()[he_id].face];
        auto end = circumcenters[graph.getHalfEdges()[he_id ^ 1].face];

        if (!start.second && !end.second) // Only draw finite lines
        {

          doc << svg::Line(svg::Point(start.first[0], start.first[1]), svg::Point(end.first[0], end.first[1]),
            svg::Stroke(0.01, svg::Color::Red));
        }
        else if (start.second && !end.second)
        {
          // Draw a line to infinity in the direction of the circumcenter
          svg::Point end_point(end.first[0], end.first[1]);
          svg::Point start_point(end.first[0] + 1000 * start.first[0], end.first[1] + 1000 * start.first[1]);
          doc << svg::Line(start_point, end_point, svg::Stroke(0.01, svg::Color::Red));
        }
        else if (!start.second && end.second)
        {
          // Draw a line to infinity in the direction of the circumcenter
          svg::Point start_point(start.first[0], start.first[1]);
          svg::Point end_point(start.first[0] + 1000 * end.first[0], start.first[1] + 1000 * end.first[1]);
          doc << svg::Line(start_point, end_point, svg::Stroke(0.01, svg::Color::Red));
        }
        else
        {
          // This should only happen if all points lie on one line, maybe we should implement it later
          KINDS_WARNING("Both circumcenters are infinite, skipping edge.");
        }

        // TODO:: Handle infinite edges
      }
    }

    // Draw Voronoi vertices
    for (const auto& circumcenter : circumcenters)
    {
      if (!circumcenter.second) // Only draw finite vertices
      {
        doc << svg::Circle(svg::Point(circumcenter.first[0], circumcenter.first[1]), 0.02, svg::Fill(svg::Color::Red),
          svg::Stroke(0.0, svg::Color::Black));
      }
    }

    if (also_draw_delaunay)
    {
      // Draw edges
      for (size_t he_id = 0; he_id < graph.getHalfEdges().size(); he_id += 2)
      {
        const HalfEdgeDelaunayGraph::HalfEdge& he = graph.getHalfEdges()[he_id];
        if (he.origin != -1
          && graph.getHalfEdges()[he_id ^ 1].origin != -1) // Only draw edges that are not boundary edges
        {
          glm::dvec2 start = points[graph.getHalfEdges()[he_id].origin];
          glm::dvec2 end = points[graph.getHalfEdges()[he_id ^ 1].origin];
          doc << svg::Line(
            svg::Point(start[0], start[1]), svg::Point(end[0], end[1]), svg::Stroke(0.01, svg::Color::Black));
        }
      }

      // Draw vertices
      for (auto& point : points)
      {
        doc << svg::Circle(
          svg::Point(point[0], point[1]), 0.02, svg::Fill(svg::Color::Blue), svg::Stroke(0.0, svg::Color::Black));
      }
    }

    doc.save();
  }
};
} // namespace kinDS

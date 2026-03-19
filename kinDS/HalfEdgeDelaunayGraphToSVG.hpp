#pragma once

#include "HalfEdgeDelaunayGraph.hpp"
#include "Logger.hpp"
#include "simple_svg.hpp"
#include <glm/glm.hpp>
#include <array>

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

  // Check if a point lies within the (optionally slightly expanded) bounding box.
  static bool isWithinBoundingBox(const glm::dvec2& p, const BoundingBox& bb, double tolerance = 0.0)
  {
    const double min_x = bb.min_x - tolerance;
    const double max_x = bb.max_x + tolerance;
    const double min_y = bb.min_y - tolerance;
    const double max_y = bb.max_y + tolerance;
    return (p.x >= min_x && p.x <= max_x && p.y >= min_y && p.y <= max_y);
  }

  // Clip a line segment to the bounding box using Liang–Barsky. Returns false if it lies completely outside.
  static bool clipSegmentToBoundingBox(glm::dvec2& p0, glm::dvec2& p1, const BoundingBox& bb, double tolerance = 0.0)
  {
    double x_min = bb.min_x - tolerance;
    double x_max = bb.max_x + tolerance;
    double y_min = bb.min_y - tolerance;
    double y_max = bb.max_y + tolerance;

    double dx = p1.x - p0.x;
    double dy = p1.y - p0.y;

    double u1 = 0.0;
    double u2 = 1.0;

    auto clip = [&](double p, double q) -> bool
    {
      if (p == 0.0)
      {
        if (q < 0.0)
        {
          return false;
        }
        return true;
      }
      double r = q / p;
      if (p < 0.0)
      {
        if (r > u2)
        {
          return false;
        }
        if (r > u1)
        {
          u1 = r;
        }
      }
      else if (p > 0.0)
      {
        if (r < u1)
        {
          return false;
        }
        if (r < u2)
        {
          u2 = r;
        }
      }
      return true;
    };

    if (!clip(-dx, p0.x - x_min))
      return false;
    if (!clip(dx, x_max - p0.x))
      return false;
    if (!clip(-dy, p0.y - y_min))
      return false;
    if (!clip(dy, y_max - p0.y))
      return false;
    if (u2 < u1)
      return false;

    glm::dvec2 new_p0 = p0 + u1 * glm::dvec2(dx, dy);
    glm::dvec2 new_p1 = p0 + u2 * glm::dvec2(dx, dy);
    p0 = new_p0;
    p1 = new_p1;
    return true;
  }

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

  // [delaunay_edge_id, voronoi_edge_id, delaunay_list_index, voronoi_list_index]
  typedef std::array<size_t, 4> IntersectionDebugInfo;

  static bool lineIntersection(
    const glm::dvec2& a0, const glm::dvec2& a1, const glm::dvec2& b0, const glm::dvec2& b1, glm::dvec2& out)
  {
    auto cross2d = [](const glm::dvec2& u, const glm::dvec2& v) { return u.x * v.y - u.y * v.x; };
    glm::dvec2 r = a1 - a0;
    glm::dvec2 s = b1 - b0;
    double denom = cross2d(r, s);
    if (std::abs(denom) < 1e-12)
    {
      return false;
    }
    double t = cross2d(b0 - a0, s) / denom;
    out = a0 + t * r;
    return true;
  }

  static bool getDelaunayEdgeEndpoints(
    const std::vector<glm::dvec2>& points, const HalfEdgeDelaunayGraph& graph, size_t delaunay_edge_id, glm::dvec2& p0, glm::dvec2& p1)
  {
    size_t he0 = 2 * delaunay_edge_id;
    if (he0 + 1 >= graph.getHalfEdges().size())
    {
      return false;
    }
    const auto& he_a = graph.getHalfEdges()[he0];
    const auto& he_b = graph.getHalfEdges()[he0 ^ 1];
    if (he_a.origin == static_cast<size_t>(-1) || he_b.origin == static_cast<size_t>(-1))
    {
      return false;
    }
    p0 = points[he_a.origin];
    p1 = points[he_b.origin];
    return true;
  }

  static bool getVoronoiEdgeEndpoints(const HalfEdgeDelaunayGraph& graph,
    const std::vector<std::pair<glm::dvec2, bool>>& circumcenters, size_t voronoi_edge_id, glm::dvec2& q0, glm::dvec2& q1)
  {
    size_t he0 = 2 * voronoi_edge_id;
    if (he0 + 1 >= graph.getHalfEdges().size())
    {
      return false;
    }
    size_t start_face = graph.getHalfEdges()[he0].face;
    size_t end_face = graph.getHalfEdges()[he0 ^ 1].face;
    const auto& start = circumcenters[start_face];
    const auto& end = circumcenters[end_face];

    if (!start.second && !end.second)
    {
      q0 = start.first;
      q1 = end.first;
      return true;
    }
    if (start.second && !end.second)
    {
      q1 = end.first;
      q0 = end.first + 1000.0 * start.first;
      return true;
    }
    if (!start.second && end.second)
    {
      q0 = start.first;
      q1 = start.first + 1000.0 * end.first;
      return true;
    }
    return false;
  }

  static std::vector<std::pair<glm::dvec2, std::pair<size_t, size_t>>> computeIntersectionMarkerData(
    const std::vector<glm::dvec2>& points, const HalfEdgeDelaunayGraph& graph,
    const std::vector<IntersectionDebugInfo>& intersection_debug_info)
  {
    std::vector<std::pair<glm::dvec2, std::pair<size_t, size_t>>> markers;
    auto circumcenters = graph.computeCircumcenters(points);
    markers.reserve(intersection_debug_info.size());

    for (const auto& item : intersection_debug_info)
    {
      size_t d_edge_id = item[0];
      size_t v_edge_id = item[1];
      size_t d_index = item[2];
      size_t v_index = item[3];

      glm::dvec2 p0, p1, q0, q1;
      if (!getDelaunayEdgeEndpoints(points, graph, d_edge_id, p0, p1))
      {
        continue;
      }
      if (!getVoronoiEdgeEndpoints(graph, circumcenters, v_edge_id, q0, q1))
      {
        continue;
      }

      glm::dvec2 intersection;
      if (!lineIntersection(p0, p1, q0, q1, intersection))
      {
        continue;
      }
      markers.push_back({ intersection, { d_index, v_index } });
    }
    return markers;
  }

  /**
   * @brief Converts the half-edge Delaunay graph to an SVG representation.
   *
   * @param points Vertex positions.
   * @param graph The half-edge Delaunay graph to convert.
   * @param filename The name of the output SVG file.
   * @param margin Margin around the bounding box (from points only).
   * @param face_inside Optional: per-face inside flag for coloring.
   * @param draw_voronoi_edges If true, draw Voronoi edges in red (bbox unchanged; edges may extend beyond it).
   * @param voronoi_vertex_to_tri Optional: mapping from Voronoi vertex id (index into circumcenters)
   *        to containing triangle id. If provided, Voronoi vertices are labeled as "id/tri".
   *        If not provided, only "id" is used (no slash).
   */
  static void write(const std::vector<glm::dvec2> points, const HalfEdgeDelaunayGraph& graph,
    const std::string& filename, double margin = 0.0, const std::vector<bool>* face_inside = nullptr,
    bool draw_voronoi_edges = false, const std::vector<size_t>* voronoi_vertex_to_tri = nullptr,
    const std::vector<IntersectionDebugInfo>* intersection_debug_info = nullptr)
  {
    BoundingBox bb = computeBoundingBox(points, margin);
    svg::Document doc = setupDocument(points, filename, bb);

    struct Label
    {
      double x;
      double y;
      std::string text;
      svg::Color color;
      double font_size;
      Label(double x_, double y_, std::string text_, const svg::Color& color_, double font_size_)
        : x(x_), y(y_), text(std::move(text_)), color(color_), font_size(font_size_)
      {
      }
    };
    std::vector<Label> labels;

    // Draw faces
    for (size_t face_id = 0; face_id < graph.getFaces().size(); face_id++)
    {
      auto face_vertex_indices = graph.getTriangleVertexIndices(face_id);

      // If any vertex is infinite, we still want to label the (infinite) face,
      // but we do not draw a filled triangle.
      if (face_vertex_indices[0] == -1 || face_vertex_indices[1] == -1 || face_vertex_indices[2] == -1)
      {
        // assert that face is outside
        if (face_inside)
        {
          assert(!(*face_inside)[face_id] && "Face is inside despite being infinite!");
        }

        // Find the two finite vertices (the edge opposite the infinite vertex).
        int finite_indices[2];
        int count = 0;
        for (int k = 0; k < 3; ++k)
        {
          if (face_vertex_indices[k] != -1)
          {
            finite_indices[count++] = face_vertex_indices[k];
          }
        }
        if (count == 2)
        {

          if(face_vertex_indices[1] == -1)
          {
            std::swap(finite_indices[0], finite_indices[1]);
          }
          glm::dvec2 a = points[finite_indices[0]];
          glm::dvec2 b = points[finite_indices[1]];
          glm::dvec2 midpoint = 0.5 * (a + b);

          // Offset the label along a vector orthogonal to the edge (consistent side).
          glm::dvec2 edge_dir = b - a;
          double len2 = glm::dot(edge_dir, edge_dir);
          if (len2 == 0.0)
          {
            edge_dir = glm::dvec2(1.0, 0.0);
          }
          else
          {
            edge_dir /= std::sqrt(len2);
          }
          glm::dvec2 edge_normal(edge_dir.y, -edge_dir.x);
          glm::dvec2 label_pos = midpoint + 0.03 * edge_normal;

          labels.push_back(Label(
            label_pos.x,
            label_pos.y,
            std::to_string(face_id),
            svg::Color(svg::Color::Black),
            0.01));
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
      labels.push_back(Label(incenter[0], incenter[1], std::to_string(face_id), svg::Color(svg::Color::White), 0.01));
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
        if (clipSegmentToBoundingBox(start, end, bb))
        {
          doc << svg::Line(
            svg::Point(start[0], start[1]), svg::Point(end[0], end[1]), svg::Stroke(0.01, svg::Color::Black));
        }
      }
    }

    // Draw vertices (filter out extreme outliers)
    for (size_t v = 0; v < points.size(); v++)
    {
      glm::dvec2 point = points[v];
      if (!isWithinBoundingBox(point, bb))
      {
        continue;
      }
      doc << svg::Circle(
        svg::Point(point[0], point[1]), 0.02, svg::Fill(svg::Color::Blue), svg::Stroke(0.0, svg::Color::Black));

      // Draw vertex id
      labels.push_back(Label(point[0] - 0.005, point[1] - 0.005, std::to_string(v), svg::Color(svg::Color::White), 0.01));
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
          glm::dvec2 edge_normal(edge_dir[1], -edge_dir[0]); // Rotate 90 degrees to get normal
          glm::dvec2 label_pos
            = midpoint + std::pow(-1, i) * 0.01 * edge_normal - glm::dvec2(0.005, 0.005); // Offset by 0.02 units
          labels.push_back(
            Label(label_pos[0], label_pos[1], std::to_string(he_id + i), svg::Color(svg::Color::Yellow), 0.01));
        }
      }
    }

    if (draw_voronoi_edges)
    {
      auto circumcenters = graph.computeCircumcenters(points);
      svg::Color purple(128, 0, 128);
      for (size_t he_id = 0; he_id < graph.getHalfEdges().size(); he_id += 2)
      {
        const HalfEdgeDelaunayGraph::HalfEdge& he = graph.getHalfEdges()[he_id];
        if (he.origin != -1 && graph.getHalfEdges()[he_id ^ 1].origin != -1)
        {
          auto start = circumcenters[graph.getHalfEdges()[he_id].face];
          auto end = circumcenters[graph.getHalfEdges()[he_id ^ 1].face];
          glm::dvec2 p0, p1;
          if (!start.second && !end.second)
          {
            p0 = start.first;
            p1 = end.first;
          }
          else if (start.second && !end.second)
          {
            // start.first is a direction, end.first is a finite point
            p1 = end.first;
            p0 = end.first + 1000.0 * start.first;
          }
          else if (!start.second && end.second)
          {
            // start.first is a finite point, end.first is a direction
            p0 = start.first;
            p1 = start.first + 1000.0 * end.first;
          }
          else
          {
            KINDS_WARNING("Both circumcenters are infinite, skipping edge.");
            continue;
          }

          if (clipSegmentToBoundingBox(p0, p1, bb))
          {
            doc << svg::Line(
              svg::Point(p0.x, p0.y), svg::Point(p1.x, p1.y), svg::Stroke(0.005, svg::Color::Red));
          }
        }
      }

      constexpr size_t invalid_id = static_cast<size_t>(-1);
      for (size_t i = 0; i < circumcenters.size(); ++i)
      {
        const auto& cc = circumcenters[i];
        if (!cc.second && isWithinBoundingBox(cc.first, bb))
        {
          doc << svg::Circle(
            svg::Point(cc.first[0], cc.first[1]), 0.02, svg::Fill(purple), svg::Stroke(0.0, svg::Color::Black));

          // Label Voronoi vertex. If crossing data is provided, use "id/containingTriId".
          // Otherwise, just use "id" (no slash).
          std::string label_text;
          if (voronoi_vertex_to_tri && i < voronoi_vertex_to_tri->size())
          {
            size_t tri_id = (*voronoi_vertex_to_tri)[i];
            if (tri_id != invalid_id)
            {
              label_text = std::to_string(i) + "/" + std::to_string(tri_id);
            }
            else
            {
              label_text = std::to_string(i);
            }
          }
          else
          {
            label_text = std::to_string(i);
          }

          labels.push_back(Label(
            cc.first[0] - 0.005,
            cc.first[1] - 0.005,
            label_text,
            svg::Color(svg::Color::White),
            0.01));
        }
      }
    }

    // Helper to draw intersection markers and labels, if the caller provides them.
    auto drawIntersections = [&](const std::vector<std::pair<glm::dvec2, std::pair<size_t, size_t>>>& intersections)
    {
      svg::Color light_blue(173, 216, 230);
      svg::Color pale_pink(255, 182, 193);
      for (const auto& inter : intersections)
      {
        glm::dvec2 p = inter.first;
        if (!isWithinBoundingBox(p, bb))
        {
          continue;
        }

        // Draw a small dot at the intersection.
        doc << svg::Circle(
          svg::Point(p.x, p.y), 0.01, svg::Fill(light_blue), svg::Stroke(0.0, light_blue));

        // Label as "(d,v)" with d in light blue and v in pale pink.
        const double x0 = p.x + 0.005;
        const double y0 = p.y + 0.005;
        const double font_size = 0.01;
        // crude glyph advance for this SVG font size/layout
        const double dx = 0.0055;

        std::string d_text = std::to_string(inter.second.first);
        std::string v_text = std::to_string(inter.second.second);

        // "(d," in light blue
        std::string prefix_text = "(" + d_text + ",";
        doc << svg::Text(svg::Point(x0, y0), prefix_text, svg::Fill(light_blue), svg::Font(font_size));

        // "v" in pale pink
        double vx = x0 + dx * static_cast<double>(prefix_text.size());
        doc << svg::Text(svg::Point(vx, y0), v_text, svg::Fill(pale_pink), svg::Font(font_size));

        // ")" in light blue
        double suffix_x = vx + dx * static_cast<double>(v_text.size());
        doc << svg::Text(svg::Point(suffix_x, y0), ")", svg::Fill(light_blue), svg::Font(font_size));
      }
    };

    if (intersection_debug_info)
    {
      auto marker_data = computeIntersectionMarkerData(points, graph, *intersection_debug_info);
      drawIntersections(marker_data);
    }

    for (const auto& label : labels)
    {
      if (!isWithinBoundingBox(glm::dvec2(label.x, label.y), bb, 0.01))
      {
        continue;
      }
      doc << svg::Text(
        svg::Point(label.x, label.y), label.text, svg::Fill(label.color), svg::Font(label.font_size));
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
        && graph.getHalfEdges()[he_id ^ 1].origin != -1)
      {
        auto start = circumcenters[graph.getHalfEdges()[he_id].face];
        auto end = circumcenters[graph.getHalfEdges()[he_id ^ 1].face];

        if (!start.second && !end.second)
        {
          doc << svg::Line(svg::Point(start.first[0], start.first[1]), svg::Point(end.first[0], end.first[1]),
            svg::Stroke(0.005, svg::Color::Red));
        }
        else if (start.second && !end.second)
        {
          svg::Point end_point(end.first[0], end.first[1]);
          svg::Point start_point(end.first[0] + 1000 * start.first[0], end.first[1] + 1000 * start.first[1]);
          doc << svg::Line(start_point, end_point, svg::Stroke(0.005, svg::Color::Red));
        }
        else if (!start.second && end.second)
        {
          svg::Point start_point(start.first[0], start.first[1]);
          svg::Point end_point(start.first[0] + 1000 * end.first[0], start.first[1] + 1000 * end.first[1]);
          doc << svg::Line(start_point, end_point, svg::Stroke(0.005, svg::Color::Red));
        }
        else
        {
          KINDS_WARNING("Both circumcenters are infinite, skipping edge.");
        }
      }
    }

    // Draw Voronoi vertices
    svg::Color purple(128, 0, 128);
    for (size_t i = 0; i < circumcenters.size(); ++i)
    {
      const auto& circumcenter = circumcenters[i];
      if (!circumcenter.second)
      {
        doc << svg::Circle(svg::Point(circumcenter.first[0], circumcenter.first[1]), 0.02,
          svg::Fill(purple), svg::Stroke(0.0, svg::Color::Black));

        // Label Voronoi vertex with "id/triId". For the generic Voronoi writer, we label with "i/i".
        std::string label_text = std::to_string(i) + "/" + std::to_string(i);
        doc << svg::Text(svg::Point(circumcenter.first[0] - 0.005, circumcenter.first[1] - 0.005),
          label_text, svg::Fill(svg::Color::White), svg::Font(0.01));
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

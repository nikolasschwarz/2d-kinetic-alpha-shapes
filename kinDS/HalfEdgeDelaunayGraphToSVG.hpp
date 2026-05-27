#pragma once

#include "HalfEdgeDelaunayGraph.hpp"
#include "Logger.hpp"
#include "VisualDebugHighlight.hpp"
#include "simple_svg.hpp"
#include <array>
#include <glm/glm.hpp>
#include <sstream>
#include <string>
#include <unordered_set>
#include <vector>

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

  /// Multiplier for text size and label placement offsets in @ref write and @ref writeVoronoi.
  static constexpr double label_scale = 15.0;

  /// Scales all SVG user units (geometry, strokes, fonts) so viewers can zoom in further.
  static constexpr double coordinate_scale = 10.0;

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
    const double width = (bb.max_x - bb.min_x) * coordinate_scale;
    const double height = (bb.max_y - bb.min_y) * coordinate_scale;

    return svg::Document(filename,
      svg::Layout(svg::Dimensions(width, height), svg::Layout::TopLeft, coordinate_scale,
        svg::Point(-bb.min_x, -bb.min_y)));
  }

  struct IntersectionDebugInfo
  {
    size_t delaunay_edge_id = 0;
    size_t voronoi_edge_id = 0;
    size_t delaunay_list_index = 0;
    size_t voronoi_list_index = 0;
    size_t prev_segment_mesh_pair_index = static_cast<size_t>(-1);
    size_t next_segment_mesh_pair_index = static_cast<size_t>(-1);
    double delaunay_edge_param = 0.0;
  };

  using IntersectionMarker = std::pair<glm::dvec2, IntersectionDebugInfo>;

  static std::string formatIntersectionParam(double delaunay_edge_param)
  {
    std::ostringstream oss;
    oss.setf(std::ios::fixed);
    oss.precision(6);
    oss << delaunay_edge_param;
    return oss.str();
  }

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

  static bool getDelaunayEdgeEndpoints(const std::vector<glm::dvec2>& points, const HalfEdgeDelaunayGraph& graph,
    size_t delaunay_edge_id, glm::dvec2& p0, glm::dvec2& p1)
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
    const std::vector<std::pair<glm::dvec2, bool>>& circumcenters, size_t voronoi_edge_id, glm::dvec2& q0,
    glm::dvec2& q1)
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

  static std::vector<IntersectionMarker> computeIntersectionMarkerData(const std::vector<glm::dvec2>& points,
    const HalfEdgeDelaunayGraph& graph, const std::vector<IntersectionDebugInfo>& intersection_debug_info)
  {
    std::vector<IntersectionMarker> markers;
    auto circumcenters = graph.computeCircumcenters(points);
    markers.reserve(intersection_debug_info.size());

    for (const IntersectionDebugInfo& item : intersection_debug_info)
    {
      glm::dvec2 p0, p1, q0, q1;
      if (!getDelaunayEdgeEndpoints(points, graph, item.delaunay_edge_id, p0, p1))
      {
        KINDS_DEBUG("Intersection marker skip: invalid Delaunay endpoints for (d,v)=("
          << item.delaunay_edge_id << "," << item.voronoi_edge_id << ") listIdx(d,v)=(" << item.delaunay_list_index
          << "," << item.voronoi_list_index << ")");
        continue;
      }
      if (!getVoronoiEdgeEndpoints(graph, circumcenters, item.voronoi_edge_id, q0, q1))
      {
        KINDS_DEBUG("Intersection marker skip: invalid Voronoi endpoints for (d,v)=("
          << item.delaunay_edge_id << "," << item.voronoi_edge_id << ") listIdx(d,v)=(" << item.delaunay_list_index
          << "," << item.voronoi_list_index << ")");
        continue;
      }

      glm::dvec2 intersection;
      if (!lineIntersection(p0, p1, q0, q1, intersection))
      {
        KINDS_DEBUG("Intersection marker skip: line intersection ill-defined for (d,v)=("
          << item.delaunay_edge_id << "," << item.voronoi_edge_id << ") listIdx(d,v)=(" << item.delaunay_list_index
          << "," << item.voronoi_list_index << ")");
        continue;
      }
      markers.push_back({ intersection, item });
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
    const std::vector<IntersectionDebugInfo>* intersection_debug_info = nullptr,
    const VisualDebugHighlight* highlight = nullptr)
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
        : x(x_)
        , y(y_)
        , text(std::move(text_))
        , color(color_)
        , font_size(font_size_)
      {
      }
    };
    std::vector<Label> labels;

    constexpr double label_unit = 0.01 * label_scale;
    const double label_font_size = label_unit;
    const double label_pad_sm = 0.5 * label_unit;
    const double label_pad_inf_face = 3.0 * label_unit;
    const double label_pad_he_along = label_unit;
    const double label_glyph_advance = 0.55 * label_unit;
    const double intersection_label_glyph_advance = 0.42 * label_unit;
    const double label_secondary_line_dy = 1.2 * label_unit;
    const double label_bbox_margin = label_unit;

    const bool selective = highlight != nullptr;

    std::unordered_set<size_t> delaunay_face_label_ids;
    if (selective)
    {
      for (size_t face_id : highlight->delaunay_faces)
      {
        delaunay_face_label_ids.insert(face_id);
        for (size_t he : graph.getFaces()[face_id].half_edges)
        {
          const int neighbor_face = graph.getHalfEdges()[he ^ 1].face;
          if (neighbor_face >= 0)
          {
            delaunay_face_label_ids.insert(static_cast<size_t>(neighbor_face));
          }
        }
      }
    }

    const svg::Color dim_inside_face(210, 235, 210);
    const svg::Color dim_outside_face(240, 215, 215);
    const svg::Color dim_edge(170, 170, 170);
    const svg::Color site_vertex_circle(173, 216, 230);
    const svg::Color site_vertex_circle_hi(135, 206, 250);
    const svg::Color dim_voronoi_edge(200, 160, 200);
    const svg::Color dim_voronoi_vertex(160, 120, 160);
    const svg::Color hi_inside_face(255, 210, 60);
    const svg::Color hi_outside_face(255, 90, 90);
    const svg::Color hi_edge(255, 120, 0);
    const svg::Color hi_primary_edge(220, 0, 120);
    const svg::Color hi_voronoi_edge(200, 0, 200);
    const svg::Color hi_voronoi_vertex(120, 0, 180);

    auto triangle_inside = [&](size_t face_id) -> bool
    {
      return face_inside && face_id < face_inside->size() && (*face_inside)[face_id];
    };

    auto face_fill_color = [&](size_t face_id) -> svg::Color
    {
      const bool inside = triangle_inside(face_id);
      if (!selective)
      {
        return inside ? svg::Color(svg::Color::Green) : svg::Color(svg::Color::Red);
      }
      const bool affected = highlight->affectsDelaunayFace(face_id);
      if (inside)
      {
        return affected ? hi_inside_face : dim_inside_face;
      }
      return affected ? hi_outside_face : dim_outside_face;
    };

    auto label_delaunay_face = [&](size_t face_id, double x, double y, const svg::Color& color)
    {
      if (!selective || delaunay_face_label_ids.find(face_id) != delaunay_face_label_ids.end())
      {
        labels.push_back(Label(x, y, std::to_string(face_id), color, label_font_size));
      }
    };

    constexpr double delaunay_edge_stroke_w = 0.006;
    constexpr double voronoi_edge_stroke_w = 0.004;

    auto delaunay_edge_color = [&](size_t voronoi_edge_id) -> svg::Color
    {
      if (!selective)
      {
        return svg::Color(svg::Color::Black);
      }
      if (highlight->affectsPrimaryVoronoiEdge(voronoi_edge_id))
      {
        return hi_primary_edge;
      }
      if (highlight->affectsVoronoiEdge(voronoi_edge_id))
      {
        return hi_edge;
      }
      return dim_edge;
    };

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

          if (face_vertex_indices[1] == -1)
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
          glm::dvec2 label_pos = midpoint + label_pad_inf_face * edge_normal;

          label_delaunay_face(face_id, label_pos.x, label_pos.y, svg::Color(svg::Color::Black));
        }

        continue;
      }

      const svg::Color face_color = face_fill_color(face_id);

      std::array<glm::dvec2, 3> face_vertices
        = { points[face_vertex_indices[0]], points[face_vertex_indices[1]], points[face_vertex_indices[2]] };

      svg::Polygon face { svg::Fill(face_color) };
      face << svg::Point(face_vertices[0][0], face_vertices[0][1])
           << svg::Point(face_vertices[1][0], face_vertices[1][1])
           << svg::Point(face_vertices[2][0], face_vertices[2][1]);
      doc << face;

      glm::dvec2 incenter = triangleIncenter(face_vertices[0], face_vertices[1], face_vertices[2]);
      label_delaunay_face(face_id, incenter[0], incenter[1], svg::Color(svg::Color::White));
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
          const size_t voronoi_edge_id = he_id / 2;
          const double stroke_w = selective ? delaunay_edge_stroke_w : 0.01;
          doc << svg::Line(svg::Point(start[0], start[1]), svg::Point(end[0], end[1]),
            svg::Stroke(stroke_w, delaunay_edge_color(voronoi_edge_id)));
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
      const bool affected = !selective || highlight->affectsDelaunayVertex(v);
      const double radius = affected ? 0.028 : 0.012;
      const svg::Color fill_c = affected ? site_vertex_circle_hi : site_vertex_circle;
      doc << svg::Circle(svg::Point(point[0], point[1]), radius, svg::Fill(fill_c), svg::Stroke(0.0, svg::Color::Black));

      if (!selective || highlight->affectsDelaunayVertex(v))
      {
        labels.push_back(Label(point[0] - label_pad_sm, point[1] - label_pad_sm, std::to_string(v),
          svg::Color(svg::Color::Black), label_font_size));
      }
    }

    for (size_t he_id = 0; he_id < graph.getHalfEdges().size(); he_id += 2)
    {
      const HalfEdgeDelaunayGraph::HalfEdge& he = graph.getHalfEdges()[he_id];
      // Draw half-edge id at midpoint but slightly offset to the left in the direction of the edge normal

      if (he.origin != -1 && graph.getHalfEdges()[he_id ^ 1].origin != -1)
      {
        glm::dvec2 start = points[graph.getHalfEdges()[he_id].origin];
        glm::dvec2 end = points[graph.getHalfEdges()[he_id ^ 1].origin];
        glm::dvec2 midpoint = (start + end) / 2.0;
        glm::dvec2 edge_dir = glm::normalize(end - start);
        glm::dvec2 edge_normal(edge_dir[1], -edge_dir[0]);

        // Do this for both half-edges
        for (size_t i = 0; i < 2; i++)
        {
          const size_t current_he_id = he_id + i;
          const int source = graph.getHalfEdges()[current_he_id].origin;
          const int destination = graph.getHalfEdges()[current_he_id ^ 1].origin;
          const std::string label_text = std::to_string(current_he_id) + " (" + std::to_string(source) + " --> "
            + std::to_string(destination) + ")";
          // Rotate 90 degrees to get normal
          glm::dvec2 label_pos = midpoint + std::pow(-1, i) * label_pad_he_along * edge_normal
            - glm::dvec2(label_pad_sm, label_pad_sm);
          const size_t undirected_edge_id = he_id / 2;
          if (!selective || highlight->affectsDirectedHalfEdge(current_he_id)
            || highlight->affectsPrimaryVoronoiEdge(undirected_edge_id))
          {
            const svg::Color label_color = selective && highlight->affectsPrimaryVoronoiEdge(undirected_edge_id)
              ? hi_primary_edge
              : hi_edge;
            labels.push_back(Label(label_pos[0], label_pos[1], label_text, label_color, label_font_size));
          }
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
            // KINDS_WARNING("Both circumcenters are infinite, skipping edge.");
            continue;
          }

          if (clipSegmentToBoundingBox(p0, p1, bb))
          {
            const size_t voronoi_edge_id = he_id / 2;
            svg::Color stroke_c = dim_voronoi_edge;
            if (selective)
            {
              if (highlight->affectsPrimaryVoronoiEdge(voronoi_edge_id))
              {
                stroke_c = hi_primary_edge;
              }
              else if (highlight->affectsVoronoiEdge(voronoi_edge_id))
              {
                stroke_c = hi_voronoi_edge;
              }
            }
            doc << svg::Line(svg::Point(p0.x, p0.y), svg::Point(p1.x, p1.y),
              svg::Stroke(voronoi_edge_stroke_w, stroke_c));
          }
        }
      }

      constexpr size_t invalid_id = static_cast<size_t>(-1);
      for (size_t i = 0; i < circumcenters.size(); ++i)
      {
        const auto& cc = circumcenters[i];
        if (!cc.second && isWithinBoundingBox(cc.first, bb))
        {
          const bool affected = !selective || highlight->affectsVoronoiVertex(i);
          const double radius = affected ? 0.026 : 0.012;
          const svg::Color fill_c = affected ? hi_voronoi_vertex : dim_voronoi_vertex;
          doc << svg::Circle(
            svg::Point(cc.first[0], cc.first[1]), radius, svg::Fill(fill_c), svg::Stroke(0.0, svg::Color::Black));

          if (!selective || highlight->affectsVoronoiVertex(i))
          {
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

            labels.push_back(Label(cc.first[0] - label_pad_sm, cc.first[1] - label_pad_sm, label_text,
              svg::Color(svg::Color::White), label_font_size));
          }
        }
      }
    }

    // Helper to draw intersection markers and labels, if the caller provides them.
    auto drawIntersections = [&](const std::vector<IntersectionMarker>& intersections)
    {
      svg::Color light_blue(173, 216, 230);
      svg::Color pale_pink(255, 182, 193);
      svg::Color mint_green(170, 255, 170);
      for (const auto& inter : intersections)
      {
        const IntersectionDebugInfo& info = inter.second;
        const size_t d_list_index = info.delaunay_list_index;
        const size_t v_list_index = info.voronoi_list_index;
        const size_t voronoi_edge_id = info.voronoi_edge_id;
        const size_t delaunay_edge_id = info.delaunay_edge_id;

        glm::dvec2 p = inter.first;
        if (!std::isfinite(p.x) || !std::isfinite(p.y))
        {
          KINDS_DEBUG("Intersection marker skip: non-finite point for (d,v)=("
            << delaunay_edge_id << "," << voronoi_edge_id << ") listIdx(" << d_list_index << "," << v_list_index
            << ")");
          continue;
        }
        if (!isWithinBoundingBox(p, bb))
        {
          KINDS_DEBUG("Intersection marker skip: outside bbox for (d,v)=("
            << delaunay_edge_id << "," << voronoi_edge_id << ") at (" << p.x << "," << p.y << ")");
          continue;
        }

        const bool emphasized = selective && highlight->emphasizesCrossing(delaunay_edge_id, voronoi_edge_id);
        const double dot_radius = emphasized ? 0.022 : 0.01;
        const svg::Color dot_fill = emphasized ? hi_voronoi_vertex : light_blue;

        const bool show_intersection_labels
          = !selective || highlight->shouldLabelCrossing(delaunay_edge_id, voronoi_edge_id);

        const std::string group_id = "intersection_dedge" + std::to_string(delaunay_edge_id) + "_vedge"
          + std::to_string(voronoi_edge_id) + "_di" + std::to_string(d_list_index) + "_vi"
          + std::to_string(v_list_index);
        svg::Group intersection_group(group_id);
        intersection_group << svg::Circle(
          svg::Point(p.x, p.y), dot_radius, svg::Fill(dot_fill), svg::Stroke(0.0, dot_fill));

        if (show_intersection_labels)
        {
          // Keep edge ids and per-edge list positions explicit; these labels are read in SVG editors.
          const double x0 = p.x + label_pad_sm;
          const double y0 = p.y + label_pad_sm;
          const double dx = intersection_label_glyph_advance;

          const std::string d_edge_text = "d=" + std::to_string(delaunay_edge_id) + ",";
          const std::string v_edge_text = "v=" + std::to_string(voronoi_edge_id) + ",";
          const std::string param_text = "t=" + formatIntersectionParam(info.delaunay_edge_param);
          const std::string d_index_text = "dIdx=" + std::to_string(d_list_index) + ",";
          const std::string v_index_text = "vIdx=" + std::to_string(v_list_index);

          intersection_group << svg::Text(
            svg::Point(x0, y0), d_edge_text, svg::Fill(light_blue), svg::Font(label_font_size));

          double label_x = x0 + dx * static_cast<double>(d_edge_text.size());
          intersection_group << svg::Text(
            svg::Point(label_x, y0), v_edge_text, svg::Fill(pale_pink), svg::Font(label_font_size));

          label_x += dx * static_cast<double>(v_edge_text.size());
          intersection_group << svg::Text(
            svg::Point(label_x, y0), param_text, svg::Fill(mint_green), svg::Font(label_font_size));

          const double y1 = y0 + label_secondary_line_dy;
          intersection_group << svg::Text(
            svg::Point(x0, y1), d_index_text, svg::Fill(light_blue), svg::Font(label_font_size));

          label_x = x0 + dx * static_cast<double>(d_index_text.size());
          intersection_group << svg::Text(
            svg::Point(label_x, y1), v_index_text, svg::Fill(pale_pink), svg::Font(label_font_size));

          if (!(info.prev_segment_mesh_pair_index == static_cast<size_t>(-1)
                && info.next_segment_mesh_pair_index == static_cast<size_t>(-1)))
          {
            const double y2 = y0 + 2.0 * label_secondary_line_dy;
            const std::string prev_text = (info.prev_segment_mesh_pair_index == static_cast<size_t>(-1))
              ? "X"
              : std::to_string(info.prev_segment_mesh_pair_index);
            const std::string next_text = (info.next_segment_mesh_pair_index == static_cast<size_t>(-1))
              ? "X"
              : std::to_string(info.next_segment_mesh_pair_index);
            const std::string mesh_pair_text = "m(" + prev_text + "," + next_text + ")";
            intersection_group << svg::Text(svg::Point(x0, y2), mesh_pair_text,
              svg::Fill(mint_green), svg::Font(label_font_size));
          }
        }

        doc << intersection_group;
      }
    };

    if (intersection_debug_info)
    {
      auto marker_data = computeIntersectionMarkerData(points, graph, *intersection_debug_info);
      drawIntersections(marker_data);
    }

    for (const auto& label : labels)
    {
      if (!isWithinBoundingBox(glm::dvec2(label.x, label.y), bb, label_bbox_margin))
      {
        continue;
      }
      doc << svg::Text(svg::Point(label.x, label.y), label.text, svg::Fill(label.color), svg::Font(label.font_size));
    }

    KINDS_DEBUG("Wrote SVG: " << filename);
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

    constexpr double label_unit = 0.01 * label_scale;
    const double label_font_size = label_unit;
    const double label_pad_sm = 0.5 * label_unit;
    const svg::Color site_vertex_circle(173, 216, 230);

    // Draw vertices
    for (auto& point : points)
    {
      doc << svg::Circle(svg::Point(point[0], point[1]), 0.02, svg::Fill(site_vertex_circle),
        svg::Stroke(0.0, svg::Color::Black));
    }

    // Draw Voronoi edges
    for (size_t he_id = 0; he_id < graph.getHalfEdges().size(); he_id += 2)
    {
      const HalfEdgeDelaunayGraph::HalfEdge& he = graph.getHalfEdges()[he_id];
      if (he.origin != -1 && graph.getHalfEdges()[he_id ^ 1].origin != -1)
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
        doc << svg::Circle(svg::Point(circumcenter.first[0], circumcenter.first[1]), 0.02, svg::Fill(purple),
          svg::Stroke(0.0, svg::Color::Black));

        // Label Voronoi vertex with "id/triId". For the generic Voronoi writer, we label with "i/i".
        std::string label_text = std::to_string(i) + "/" + std::to_string(i);
        doc << svg::Text(svg::Point(circumcenter.first[0] - label_pad_sm, circumcenter.first[1] - label_pad_sm),
          label_text, svg::Fill(svg::Color::White), svg::Font(label_font_size));
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
        doc << svg::Circle(svg::Point(point[0], point[1]), 0.02, svg::Fill(site_vertex_circle),
          svg::Stroke(0.0, svg::Color::Black));
      }
    }

    KINDS_DEBUG("Wrote SVG: " << filename);
    doc.save();
  }
};
} // namespace kinDS

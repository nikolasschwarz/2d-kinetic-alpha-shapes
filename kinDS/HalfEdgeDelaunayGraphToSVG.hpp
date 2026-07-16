#pragma once

#include "HalfEdgeDelaunayGraph.hpp"
#include "Logger.hpp"
#include "VisualDebugHighlight.hpp"
#include "simple_svg.hpp"
#include <algorithm>
#include <array>
#include <glm/glm.hpp>
#include <limits>
#include <sstream>
#include <string>
#include <unordered_map>
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

  /// Positions for the two directed labels on one Delaunay edge; enforces @p min_twin_label_dy on y separation.
  static std::array<glm::dvec2, 2> twinHalfEdgeLabelPositions(const glm::dvec2& midpoint,
    const glm::dvec2& edge_dir_normalized, double label_pad_he_along, double label_pad_sm, double min_twin_label_dy)
  {
    const glm::dvec2 edge_normal(edge_dir_normalized.y, -edge_dir_normalized.x);
    const glm::dvec2 corner_offset(-label_pad_sm, -label_pad_sm);
    std::array<glm::dvec2, 2> positions {
      midpoint + label_pad_he_along * edge_normal + corner_offset,
      midpoint - label_pad_he_along * edge_normal + corner_offset,
    };

    const double dy = positions[1].y - positions[0].y;
    if (std::abs(dy) < min_twin_label_dy)
    {
      const double center_y = 0.5 * (positions[0].y + positions[1].y);
      positions[0].y = center_y - 0.5 * min_twin_label_dy;
      positions[1].y = center_y + 0.5 * min_twin_label_dy;
    }

    return positions;
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

  /// Base position and offset end position for a separated strand site.
  struct SeparationOffsetSegment
  {
    glm::dvec2 base { 0.0, 0.0 };
    glm::dvec2 end { 0.0, 0.0 };
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

  static std::string formatShortMeshPairIndex(size_t mesh_pair_index)
  {
    if (mesh_pair_index == static_cast<size_t>(-1))
    {
      return "-";
    }
    return std::to_string(mesh_pair_index);
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
    size_t delaunay_edge_id, glm::dvec2& p0, glm::dvec2& p1,
    const std::unordered_set<size_t>* positioned_strands = nullptr)
  {
    size_t he0 = 2 * delaunay_edge_id;
    if (he0 + 1 >= graph.halfEdgeSlotCount())
    {
      return false;
    }
    const auto& he_a = graph.halfEdge(he0);
    const auto& he_b = graph.halfEdge(he0 ^ 1);
    if (he_a.origin == static_cast<size_t>(-1) || he_b.origin == static_cast<size_t>(-1))
    {
      return false;
    }
    const size_t origin_a = static_cast<size_t>(he_a.origin);
    const size_t origin_b = static_cast<size_t>(he_b.origin);
    if (positioned_strands != nullptr
      && (positioned_strands->count(origin_a) == 0 || positioned_strands->count(origin_b) == 0))
    {
      return false;
    }
    if (origin_a >= points.size() || origin_b >= points.size())
    {
      return false;
    }
    p0 = points[origin_a];
    p1 = points[origin_b];
    return true;
  }

  static bool getVoronoiEdgeEndpoints(const HalfEdgeDelaunayGraph& graph,
    const std::vector<std::pair<glm::dvec2, bool>>& circumcenters, size_t voronoi_edge_id, glm::dvec2& q0,
    glm::dvec2& q1)
  {
    size_t he0 = 2 * voronoi_edge_id;
    if (he0 + 1 >= graph.halfEdgeSlotCount())
    {
      return false;
    }
    const size_t start_face = static_cast<size_t>(graph.halfEdge(he0).face);
    const size_t end_face = static_cast<size_t>(graph.halfEdge(he0 ^ 1).face);
    if (start_face >= circumcenters.size() || end_face >= circumcenters.size())
    {
      return false;
    }
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
    const HalfEdgeDelaunayGraph& graph, const std::vector<IntersectionDebugInfo>& intersection_debug_info,
    const std::unordered_set<size_t>* positioned_strands = nullptr)
  {
    std::vector<IntersectionMarker> markers;
    auto circumcenters = graph.computeCircumcenters(points);
    markers.reserve(intersection_debug_info.size());

    for (const IntersectionDebugInfo& item : intersection_debug_info)
    {
      glm::dvec2 p0, p1, q0, q1;
      if (!getDelaunayEdgeEndpoints(points, graph, item.delaunay_edge_id, p0, p1, positioned_strands))
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
   * @param component_strands Optional: when set, only geometry belonging to this connected component is drawn
   *        and the view is cropped to its strand sites.
   * @param site_input_branch_labels Optional: strand id -> input branch id (from @ref StrandTree::getBranchIndex at
   *        the debug time). Appended to the site label as `branch=X`. Whole-number @p occurrence_time uses that
   *        section index; fractional times use the next section (ceil).
   * @param positioned_strands Optional: strand ids for which @p points contains valid coordinates. When set,
   *        geometry that needs an unpositioned strand is skipped instead of indexing outside the live set.
   * @param seam_outlines Optional: closed loops (world/profile positions) drawn as bright magenta polygons, used to
   *        visualize the extracted future-branch seam outlines that feed the pending-split convexity check.
   */
  static void write(const std::vector<glm::dvec2> points, const HalfEdgeDelaunayGraph& graph,
    const std::string& filename, double margin = 0.0, const std::vector<bool>* face_inside = nullptr,
    bool draw_voronoi_edges = false, const std::vector<size_t>* voronoi_vertex_to_tri = nullptr,
    const std::vector<IntersectionDebugInfo>* intersection_debug_info = nullptr,
    const VisualDebugHighlight* highlight = nullptr, const std::unordered_set<size_t>* component_strands = nullptr,
    const std::unordered_map<size_t, size_t>* site_input_branch_labels = nullptr,
    const std::unordered_set<size_t>* positioned_strands = nullptr,
    const std::unordered_map<size_t, glm::dvec3>* site_world_positions = nullptr,
    const std::unordered_map<size_t, glm::dvec3>* voronoi_vertex_world_positions = nullptr,
    const std::vector<SeparationOffsetSegment>* separation_offset_segments = nullptr,
    const std::vector<std::vector<glm::dvec2>>* seam_outlines = nullptr)
  {
    auto in_component = [&](size_t strand_id) -> bool
    {
      return component_strands == nullptr || component_strands->count(strand_id) != 0;
    };

    auto try_point = [&](size_t strand_id) -> const glm::dvec2*
    {
      if (strand_id >= points.size())
      {
        return nullptr;
      }
      if (positioned_strands != nullptr && positioned_strands->count(strand_id) == 0)
      {
        return nullptr;
      }
      return &points[strand_id];
    };

    auto triangle_in_component = [&](size_t face_id) -> bool
    {
      if (component_strands == nullptr)
      {
        return true;
      }
      const auto face_vertex_indices = graph.getTriangleVertexIndices(face_id);
      bool has_member = false;
      for (int k = 0; k < 3; ++k)
      {
        if (face_vertex_indices[k] == -1)
        {
          continue;
        }
        if (!in_component(static_cast<size_t>(face_vertex_indices[k])))
        {
          return false;
        }
        has_member = true;
      }
      return has_member;
    };

    auto delaunay_edge_in_component = [&](size_t he_id) -> bool
    {
      if (component_strands == nullptr)
      {
        return true;
      }
      const HalfEdgeDelaunayGraph::HalfEdge& he = graph.halfEdge(he_id);
      const HalfEdgeDelaunayGraph::HalfEdge& twin = graph.halfEdge(he_id ^ 1);
      if (he.origin == -1 || twin.origin == -1)
      {
        const size_t finite_origin
          = he.origin != -1 ? static_cast<size_t>(he.origin) : static_cast<size_t>(twin.origin);
        return in_component(finite_origin) && triangle_in_component(he.face);
      }
      return in_component(static_cast<size_t>(he.origin)) && in_component(static_cast<size_t>(twin.origin));
    };

    BoundingBox bb = [&]() -> BoundingBox
    {
      if (component_strands != nullptr && !component_strands->empty())
      {
        double min_x = std::numeric_limits<double>::max();
        double min_y = std::numeric_limits<double>::max();
        double max_x = std::numeric_limits<double>::lowest();
        double max_y = std::numeric_limits<double>::lowest();
        for (size_t strand_id : *component_strands)
        {
          const glm::dvec2* point = try_point(strand_id);
          if (point == nullptr)
          {
            continue;
          }
          min_x = std::min(min_x, (*point)[0]);
          min_y = std::min(min_y, (*point)[1]);
          max_x = std::max(max_x, (*point)[0]);
          max_y = std::max(max_y, (*point)[1]);
        }
        if (min_x > max_x || min_y > max_y)
        {
          return computeBoundingBox(points, margin);
        }
        return BoundingBox(min_x - margin, min_y - margin, max_x + margin, max_y + margin);
      }
      if (positioned_strands != nullptr && !positioned_strands->empty())
      {
        double min_x = std::numeric_limits<double>::max();
        double min_y = std::numeric_limits<double>::max();
        double max_x = std::numeric_limits<double>::lowest();
        double max_y = std::numeric_limits<double>::lowest();
        for (size_t strand_id : *positioned_strands)
        {
          const glm::dvec2* point = try_point(strand_id);
          if (point == nullptr)
          {
            continue;
          }
          min_x = std::min(min_x, (*point)[0]);
          min_y = std::min(min_y, (*point)[1]);
          max_x = std::max(max_x, (*point)[0]);
          max_y = std::max(max_y, (*point)[1]);
        }
        if (min_x > max_x || min_y > max_y)
        {
          return computeBoundingBox(points, margin);
        }
        return BoundingBox(min_x - margin, min_y - margin, max_x + margin, max_y + margin);
      }
      return computeBoundingBox(points, margin);
    }();
    if (separation_offset_segments != nullptr)
    {
      for (const SeparationOffsetSegment& segment : *separation_offset_segments)
      {
        bb.min_x = std::min({ bb.min_x, segment.base.x, segment.end.x });
        bb.min_y = std::min({ bb.min_y, segment.base.y, segment.end.y });
        bb.max_x = std::max({ bb.max_x, segment.base.x, segment.end.x });
        bb.max_y = std::max({ bb.max_y, segment.base.y, segment.end.y });
      }
    }
    if (seam_outlines != nullptr)
    {
      for (const std::vector<glm::dvec2>& outline : *seam_outlines)
      {
        for (const glm::dvec2& p : outline)
        {
          bb.min_x = std::min(bb.min_x, p.x);
          bb.min_y = std::min(bb.min_y, p.y);
          bb.max_x = std::max(bb.max_x, p.x);
          bb.max_y = std::max(bb.max_y, p.y);
        }
      }
    }
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
    std::unordered_set<size_t> labeled_delaunay_vertices;

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
        for (size_t he : graph.face(face_id).half_edges)
        {
          const int neighbor_face = graph.halfEdge(he ^ 1).face;
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
    const svg::Color dim_voronoi_edge(220, 0, 0);
    const svg::Color dim_voronoi_vertex(160, 120, 160);
    const svg::Color hi_inside_face(255, 210, 60);
    const svg::Color hi_outside_face(255, 90, 90);
    const svg::Color hi_edge(255, 120, 0);
    const svg::Color hi_primary_edge(220, 0, 120);
    const svg::Color hi_voronoi_edge(255, 0, 0);
    const svg::Color hi_voronoi_vertex(120, 0, 180);

    auto triangle_inside = [&](size_t face_id) -> bool
    {
      return face_inside && face_id < face_inside->size() && (*face_inside)[face_id];
    };

    // Delaunay edge on the convex hull: one incident face is infinite (mesh boundary).
    auto is_convex_hull_edge = [&](size_t he_id) -> bool
    {
      return graph.isOnConvexBoundary(he_id);
    };

    const svg::Color convex_hull_edge_color(0, 0, 255);
    constexpr double convex_hull_edge_stroke_w = 0.012;

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

    auto label_delaunay_vertex = [&](size_t vertex_id, const svg::Color& color)
    {
      if (!labeled_delaunay_vertices.insert(vertex_id).second)
      {
        return;
      }
      const glm::dvec2* point = try_point(vertex_id);
      if (point == nullptr)
      {
        return;
      }
      std::ostringstream label_text;
      label_text.setf(std::ios::fixed);
      label_text.precision(6);
      label_text << vertex_id << " (" << (*point)[0] << "," << (*point)[1] << ")";
      if (site_world_positions != nullptr)
      {
        const auto world_it = site_world_positions->find(vertex_id);
        if (world_it != site_world_positions->end())
        {
          label_text << " world=(" << world_it->second.x << "," << world_it->second.y << "," << world_it->second.z
                     << ")";
        }
      }
      if (site_input_branch_labels != nullptr)
      {
        const auto branch_it = site_input_branch_labels->find(vertex_id);
        if (branch_it != site_input_branch_labels->end())
        {
          label_text << " branch=" << branch_it->second;
        }
      }
      labels.push_back(Label((*point)[0] - label_pad_sm, (*point)[1] - label_pad_sm, label_text.str(), color, label_font_size));
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
      if (highlight->affectsDirectedHalfEdge(voronoi_edge_id * 2)
        || highlight->affectsDirectedHalfEdge(voronoi_edge_id * 2 + 1))
      {
        return hi_edge;
      }
      return dim_edge;
    };

    // Draw faces
    for (size_t face_id : graph.liveFaces())
    {
      if (!triangle_in_component(face_id))
      {
        continue;
      }

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
          const glm::dvec2* a = try_point(static_cast<size_t>(finite_indices[0]));
          const glm::dvec2* b = try_point(static_cast<size_t>(finite_indices[1]));
          if (a == nullptr || b == nullptr)
          {
            continue;
          }
          glm::dvec2 midpoint = 0.5 * ((*a) + (*b));

          // Offset the label along a vector orthogonal to the edge (consistent side).
          glm::dvec2 edge_dir = (*b) - (*a);
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

      const glm::dvec2* v0 = try_point(static_cast<size_t>(face_vertex_indices[0]));
      const glm::dvec2* v1 = try_point(static_cast<size_t>(face_vertex_indices[1]));
      const glm::dvec2* v2 = try_point(static_cast<size_t>(face_vertex_indices[2]));
      if (v0 == nullptr || v1 == nullptr || v2 == nullptr)
      {
        continue;
      }

      std::array<glm::dvec2, 3> face_vertices = { *v0, *v1, *v2 };

      svg::Polygon face { svg::Fill(face_color) };
      face << svg::Point(face_vertices[0][0], face_vertices[0][1])
           << svg::Point(face_vertices[1][0], face_vertices[1][1])
           << svg::Point(face_vertices[2][0], face_vertices[2][1]);
      doc << face;

      glm::dvec2 incenter = triangleIncenter(face_vertices[0], face_vertices[1], face_vertices[2]);
      label_delaunay_face(face_id, incenter[0], incenter[1], svg::Color(svg::Color::White));
    }

    // Draw edges
    for (size_t he_id : graph.liveDelaunayEdges())
    {
      if (!delaunay_edge_in_component(he_id))
      {
        continue;
      }

      const HalfEdgeDelaunayGraph::HalfEdge& he = graph.halfEdge(he_id);
      if (he.origin != -1
        && graph.halfEdge(he_id ^ 1).origin != -1) // Only draw edges with two finite endpoints
      {
        const glm::dvec2* start = try_point(static_cast<size_t>(graph.halfEdge(he_id).origin));
        const glm::dvec2* end = try_point(static_cast<size_t>(graph.halfEdge(he_id ^ 1).origin));
        if (start == nullptr || end == nullptr)
        {
          continue;
        }
        glm::dvec2 clipped_start = *start;
        glm::dvec2 clipped_end = *end;
        if (clipSegmentToBoundingBox(clipped_start, clipped_end, bb))
        {
          const size_t voronoi_edge_id = he_id / 2;
          const bool convex_hull = is_convex_hull_edge(he_id);
          const double stroke_w = convex_hull ? convex_hull_edge_stroke_w
                                              : (selective ? delaunay_edge_stroke_w : 0.01);
          const svg::Color stroke_color
            = convex_hull ? convex_hull_edge_color : delaunay_edge_color(voronoi_edge_id);
          doc << svg::Line(svg::Point(clipped_start[0], clipped_start[1]), svg::Point(clipped_end[0], clipped_end[1]),
            svg::Stroke(stroke_w, stroke_color));
        }
      }
    }

    // Draw vertices (filter out extreme outliers)
    const auto draw_site_vertex = [&](size_t v)
    {
      if (!in_component(v))
      {
        return;
      }

      const glm::dvec2* point = try_point(v);
      if (point == nullptr || !isWithinBoundingBox(*point, bb))
      {
        return;
      }
      const bool affected = !selective || highlight->affectsDelaunayVertex(v);
      const double radius = affected ? 0.028 : 0.012;
      const svg::Color fill_c = affected ? site_vertex_circle_hi : site_vertex_circle;
      doc << svg::Circle(svg::Point((*point)[0], (*point)[1]), radius, svg::Fill(fill_c), svg::Stroke(0.0, svg::Color::Black));

      if (!selective || highlight->affectsDelaunayVertex(v))
      {
        label_delaunay_vertex(v, svg::Color(svg::Color::Black));
      }
    };

    if (positioned_strands != nullptr)
    {
      for (size_t v : *positioned_strands)
      {
        draw_site_vertex(v);
      }
    }
    else
    {
      for (size_t v = 0; v < points.size(); v++)
      {
        draw_site_vertex(v);
      }
    }

    for (size_t he_id : graph.liveDelaunayEdges())
    {
      if (!delaunay_edge_in_component(he_id))
      {
        continue;
      }

      const HalfEdgeDelaunayGraph::HalfEdge& he = graph.halfEdge(he_id);
      // Draw half-edge id at midpoint but slightly offset to the left in the direction of the edge normal

      if (he.origin != -1 && graph.halfEdge(he_id ^ 1).origin != -1)
      {
        const glm::dvec2* start = try_point(static_cast<size_t>(graph.halfEdge(he_id).origin));
        const glm::dvec2* end = try_point(static_cast<size_t>(graph.halfEdge(he_id ^ 1).origin));
        if (start == nullptr || end == nullptr)
        {
          continue;
        }
        glm::dvec2 midpoint = (*start + *end) / 2.0;
        glm::dvec2 edge_dir = (*end) - (*start);
        const double edge_len2 = glm::dot(edge_dir, edge_dir);
        if (edge_len2 > 0.0)
        {
          edge_dir /= std::sqrt(edge_len2);
        }
        else
        {
          edge_dir = glm::dvec2(1.0, 0.0);
        }

        const std::array<glm::dvec2, 2> label_positions
          = twinHalfEdgeLabelPositions(midpoint, edge_dir, label_pad_he_along, label_pad_sm, label_secondary_line_dy);

        for (size_t i = 0; i < 2; i++)
        {
          const size_t current_he_id = he_id + i;
          const int source = graph.halfEdge(current_he_id).origin;
          const int destination = graph.halfEdge(current_he_id ^ 1).origin;
          const int next = graph.halfEdge(current_he_id).next;
          const int face = graph.halfEdge(current_he_id).face;
          const std::string label_text = std::to_string(current_he_id) + " ("
            + (source < 0 ? "inf" : std::to_string(source)) + " --> "
            + (destination < 0 ? "inf" : std::to_string(destination)) + ") next=" + std::to_string(next)
            + " face=" + std::to_string(face);
          const glm::dvec2& label_pos = label_positions[i];
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
      for (size_t he_id : graph.liveDelaunayEdges())
      {
        if (!delaunay_edge_in_component(he_id))
        {
          continue;
        }

        const HalfEdgeDelaunayGraph::HalfEdge& he = graph.halfEdge(he_id);
        if (he.origin != -1 && graph.halfEdge(he_id ^ 1).origin != -1)
        {
          const size_t start_face = static_cast<size_t>(graph.halfEdge(he_id).face);
          const size_t end_face = static_cast<size_t>(graph.halfEdge(he_id ^ 1).face);
          if (start_face >= circumcenters.size() || end_face >= circumcenters.size())
          {
            continue;
          }

          auto start = circumcenters[start_face];
          auto end = circumcenters[end_face];
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
        if (!triangle_in_component(i))
        {
          continue;
        }

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
            std::ostringstream label_text;
            label_text.setf(std::ios::fixed);
            label_text.precision(6);
            if (voronoi_vertex_to_tri && i < voronoi_vertex_to_tri->size())
            {
              size_t tri_id = (*voronoi_vertex_to_tri)[i];
              if (tri_id != invalid_id)
              {
                label_text << i << "/" << tri_id;
              }
              else
              {
                label_text << i;
              }
            }
            else
            {
              label_text << i;
            }
            label_text << " (" << cc.first[0] << "," << cc.first[1] << ")";

            if (voronoi_vertex_world_positions != nullptr)
            {
              const auto world_it = voronoi_vertex_world_positions->find(i);
              if (world_it != voronoi_vertex_world_positions->end())
              {
                label_text << " world=(" << world_it->second.x << "," << world_it->second.y << ","
                           << world_it->second.z << ")";
              }
            }

            labels.push_back(Label(cc.first[0] - label_pad_sm, cc.first[1] - label_pad_sm, label_text.str(),
              svg::Color(svg::Color::White), label_font_size));
          }
        }
      }
    }

    // Helper to draw intersection markers and labels, if the caller provides them.
    auto drawIntersections = [&](const std::vector<IntersectionMarker>& intersections)
    {
      svg::Color marker_light_blue(173, 216, 230);
      svg::Color label_blue(30, 95, 120);
      svg::Color label_pink(165, 45, 75);
      svg::Color label_green(45, 130, 45);
      for (const auto& inter : intersections)
      {
        const IntersectionDebugInfo& info = inter.second;
        const size_t d_list_index = info.delaunay_list_index;
        const size_t v_list_index = info.voronoi_list_index;
        const size_t voronoi_edge_id = info.voronoi_edge_id;
        const size_t delaunay_edge_id = info.delaunay_edge_id;

        if (!delaunay_edge_in_component(2 * delaunay_edge_id))
        {
          continue;
        }

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
        const svg::Color dot_fill = emphasized ? hi_voronoi_vertex : marker_light_blue;

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
          if (voronoi_edge_id < graph.halfEdgeSlotCount() / 2)
          {
            for (size_t quad_he_id : graph.getQuadBoundaryHalfEdgeIndices(voronoi_edge_id))
            {
              const int vertex_id = graph.halfEdge(quad_he_id).origin;
              if (vertex_id >= 0)
              {
                label_delaunay_vertex(static_cast<size_t>(vertex_id), svg::Color(svg::Color::Black));
              }
            }
          }

          // Keep edge ids and per-edge list positions explicit; these labels are read in SVG editors.
          const double x0 = p.x + label_pad_sm;
          const double y0 = p.y + label_pad_sm;
          const double dx = intersection_label_glyph_advance;

          const std::string d_edge_text = "d=" + std::to_string(delaunay_edge_id) + ",";
          const std::string v_edge_text = "v=" + std::to_string(voronoi_edge_id) + ",";
          const std::string param_text = "dParam=" + formatIntersectionParam(info.delaunay_edge_param);
          const std::string d_index_text = "dIdx=" + std::to_string(d_list_index) + ",";
          const std::string v_index_text = "vIdx=" + std::to_string(v_list_index);
          int delaunay_face0 = -1;
          int delaunay_face1 = -1;
          const size_t d_he0 = 2 * delaunay_edge_id;
          if (d_he0 + 1 < graph.halfEdgeSlotCount())
          {
            delaunay_face0 = graph.halfEdge(d_he0).face;
            delaunay_face1 = graph.halfEdge(d_he0 ^ 1).face;
          }
          const std::string faces_text
            = "faces=(" + std::to_string(delaunay_face0) + "," + std::to_string(delaunay_face1) + ")";

          intersection_group << svg::Text(
            svg::Point(x0, y0), d_edge_text, svg::Fill(label_blue), svg::Font(label_font_size));

          double label_x = x0 + dx * static_cast<double>(d_edge_text.size());
          intersection_group << svg::Text(
            svg::Point(label_x, y0), v_edge_text, svg::Fill(label_pink), svg::Font(label_font_size));

          label_x += dx * static_cast<double>(v_edge_text.size());
          intersection_group << svg::Text(
            svg::Point(label_x, y0), param_text, svg::Fill(label_green), svg::Font(label_font_size));

          const double y1 = y0 + label_secondary_line_dy;
          intersection_group << svg::Text(
            svg::Point(x0, y1), d_index_text, svg::Fill(label_blue), svg::Font(label_font_size));

          label_x = x0 + dx * static_cast<double>(d_index_text.size());
          intersection_group << svg::Text(
            svg::Point(label_x, y1), v_index_text, svg::Fill(label_pink), svg::Font(label_font_size));

          const double y2 = y0 + 2.0 * label_secondary_line_dy;
          intersection_group << svg::Text(svg::Point(x0, y2), faces_text,
            svg::Fill(svg::Color(svg::Color::Black)), svg::Font(label_font_size));

          const double y3 = y0 + 3.0 * label_secondary_line_dy;
          const std::string mesh_pair_text = "p=" + formatShortMeshPairIndex(info.prev_segment_mesh_pair_index) + ",n="
            + formatShortMeshPairIndex(info.next_segment_mesh_pair_index);
          intersection_group << svg::Text(svg::Point(x0, y3), mesh_pair_text, svg::Fill(label_green),
            svg::Font(label_font_size));
        }

        doc << intersection_group;
      }
    };

    if (intersection_debug_info)
    {
      auto marker_data = computeIntersectionMarkerData(points, graph, *intersection_debug_info, positioned_strands);
      drawIntersections(marker_data);
    }

    if (separation_offset_segments != nullptr)
    {
      const svg::Color offset_line_color(0, 128, 255);
      const svg::Color offset_base_color(255, 140, 0);
      constexpr double offset_line_stroke_w = 0.01;
      constexpr double offset_base_radius = 0.018;
      for (const SeparationOffsetSegment& segment : *separation_offset_segments)
      {
        glm::dvec2 clipped_start = segment.base;
        glm::dvec2 clipped_end = segment.end;
        if (!clipSegmentToBoundingBox(clipped_start, clipped_end, bb))
        {
          continue;
        }
        doc << svg::Line(svg::Point(clipped_start[0], clipped_start[1]), svg::Point(clipped_end[0], clipped_end[1]),
          svg::Stroke(offset_line_stroke_w, offset_line_color));
        if (isWithinBoundingBox(segment.base, bb))
        {
          doc << svg::Circle(svg::Point(segment.base[0], segment.base[1]), offset_base_radius,
            svg::Fill(svg::Color::White), svg::Stroke(offset_line_stroke_w, offset_base_color));
        }
      }
    }

    if (seam_outlines != nullptr)
    {
      const svg::Color seam_outline_color(255, 0, 255);
      constexpr double seam_outline_stroke_w = 0.022;
      constexpr double seam_outline_vertex_radius = 0.03;
      for (const std::vector<glm::dvec2>& outline : *seam_outlines)
      {
        const size_t n = outline.size();
        if (n < 2)
        {
          continue;
        }
        for (size_t i = 0; i < n; ++i)
        {
          glm::dvec2 clipped_start = outline[i];
          glm::dvec2 clipped_end = outline[(i + 1) % n];
          if (clipSegmentToBoundingBox(clipped_start, clipped_end, bb))
          {
            doc << svg::Line(svg::Point(clipped_start[0], clipped_start[1]),
              svg::Point(clipped_end[0], clipped_end[1]), svg::Stroke(seam_outline_stroke_w, seam_outline_color));
          }
        }
        for (const glm::dvec2& p : outline)
        {
          if (isWithinBoundingBox(p, bb))
          {
            doc << svg::Circle(svg::Point(p[0], p[1]), seam_outline_vertex_radius, svg::Fill(seam_outline_color),
              svg::Stroke(0.0, seam_outline_color));
          }
        }
      }
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
    for (size_t he_id : graph.liveDelaunayEdges())
    {
      const HalfEdgeDelaunayGraph::HalfEdge& he = graph.halfEdge(he_id);
      if (he.origin != -1 && graph.halfEdge(he_id ^ 1).origin != -1)
      {
        const size_t start_face = static_cast<size_t>(graph.halfEdge(he_id).face);
        const size_t end_face = static_cast<size_t>(graph.halfEdge(he_id ^ 1).face);
        if (start_face >= circumcenters.size() || end_face >= circumcenters.size())
        {
          continue;
        }

        auto start = circumcenters[start_face];
        auto end = circumcenters[end_face];

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

        // Label Voronoi vertex with "id/triId" plus local coordinates, matching site labels.
        std::ostringstream label_text;
        label_text.setf(std::ios::fixed);
        label_text.precision(6);
        label_text << i << "/" << i << " (" << circumcenter.first[0] << "," << circumcenter.first[1] << ")";
        doc << svg::Text(svg::Point(circumcenter.first[0] - label_pad_sm, circumcenter.first[1] - label_pad_sm),
          label_text.str(), svg::Fill(svg::Color::White), svg::Font(label_font_size));
      }
    }

    if (also_draw_delaunay)
    {
      // Draw edges
      for (size_t he_id : graph.liveDelaunayEdges())
      {
        const HalfEdgeDelaunayGraph::HalfEdge& he = graph.halfEdge(he_id);
        if (he.origin != -1
          && graph.halfEdge(he_id ^ 1).origin != -1) // Only draw edges that are not boundary edges
        {
          glm::dvec2 start = points[graph.halfEdge(he_id).origin];
          glm::dvec2 end = points[graph.halfEdge(he_id ^ 1).origin];
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

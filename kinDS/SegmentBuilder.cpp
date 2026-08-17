#include "SegmentBuilder.hpp"

#include "DebugExportFormatting.hpp"
#include "GeometryUtils.hpp"
#include "KineticDelaunayEventPredicates.hpp"
#include "KineticDelaunayFlipEvent.hpp"
#include "KineticDelaunayRadiusEvent.hpp"
#include "Logger.hpp"
#include "SegmentBuilderCrossingCallback.hpp"
#include "SegmentBuilderFlipCallback.hpp"
#include "SegmentBuilderRadiusCallback.hpp"
#include "SegmentBuilderSectionCallback.hpp"
#include "SegmentBuilderSeparationCallback.hpp"
#include "SegmentBuilderSubdivisionCallback.hpp"
#include "SegmentBuilderVisualDebug.hpp"
#include "Validator.hpp"
#include "VisualDebugHighlight.hpp"

#include <algorithm>
#include <array>
#include <cassert>
#include <cctype>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <glm/gtx/exterior_product.hpp>
#include <glm/gtx/norm.hpp>
#include <iomanip>
#include <limits>
#include <map>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>

using namespace kinDS;

namespace
{
std::string jsonStringLiteral(const std::string& value)
{
  std::string out = "\"";
  for (char ch : value)
  {
    switch (ch)
    {
    case '\\':
      out += "\\\\";
      break;
    case '"':
      out += "\\\"";
      break;
    case '\n':
      out += "\\n";
      break;
    case '\r':
      out += "\\r";
      break;
    case '\t':
      out += "\\t";
      break;
    default:
      out += ch;
      break;
    }
  }
  out += '"';
  return out;
}

std::string numberLiteral(double value)
{
  std::ostringstream o;
  o << std::setprecision(kDebugExportTimePrecision) << std::showpoint << value;
  return o.str();
}

std::string svgEscape(const std::string& text)
{
  std::string out;
  out.reserve(text.size());
  for (char ch : text)
  {
    switch (ch)
    {
    case '&':
      out += "&amp;";
      break;
    case '<':
      out += "&lt;";
      break;
    case '>':
      out += "&gt;";
      break;
    case '"':
      out += "&quot;";
      break;
    default:
      out += ch;
      break;
    }
  }
  return out;
}

std::optional<size_t> metadataSizeField(const std::string& metadata, const char* key);

std::optional<std::string> metadataStringField(const std::string& metadata, const char* key)
{
  const std::string needle = jsonStringLiteral(key) + ":";
  const size_t key_pos = metadata.find(needle);
  if (key_pos == std::string::npos)
  {
    return std::nullopt;
  }
  size_t value_pos = key_pos + needle.size();
  while (value_pos < metadata.size() && std::isspace(static_cast<unsigned char>(metadata[value_pos])))
  {
    ++value_pos;
  }
  if (value_pos >= metadata.size() || metadata[value_pos] != '"')
  {
    return std::nullopt;
  }
  ++value_pos;
  std::string value;
  for (; value_pos < metadata.size(); ++value_pos)
  {
    const char ch = metadata[value_pos];
    if (ch == '"')
    {
      return value;
    }
    if (ch == '\\' && value_pos + 1 < metadata.size())
    {
      ++value_pos;
      value += metadata[value_pos];
      continue;
    }
    value += ch;
  }
  return std::nullopt;
}

template<typename ListIt> bool listConstIteratorPrecedes(ListIt begin, ListIt end, ListIt earlier, ListIt later)
{
  for (ListIt it = begin; it != end; ++it)
  {
    if (it == earlier)
    {
      return true;
    }
    if (it == later)
    {
      return false;
    }
  }
  return false;
}

bool triangulationPlaneXYEqual(const VoronoiMesh& mesh, size_t vertex_a, size_t vertex_b)
{
  const glm::dvec2 a = mesh.triangulationPlaneXY(vertex_a);
  const glm::dvec2 b = mesh.triangulationPlaneXY(vertex_b);
  return a.x == b.x && a.y == b.y;
}

double resolveTriangulateSimplePolygonDebugTime(const VoronoiMesh& mesh, std::optional<double> occurrence_time)
{
  if (occurrence_time.has_value() && std::isfinite(occurrence_time.value()))
  {
    return occurrence_time.value();
  }
  return mesh.getCreationKineticTime();
}

std::filesystem::path makeTriangulateSimplePolygonDebugPath(const KineticDelaunay& kin_del, const VoronoiMesh& mesh,
  const char* tag, const char* extension, std::optional<double> occurrence_time,
  std::optional<size_t> runtime_branch_id)
{
  static size_t debug_counter = 0;
  ++debug_counter;

  const double kinetic_time = resolveTriangulateSimplePolygonDebugTime(mesh, occurrence_time);
  const std::string filename = formatDebugExportTimeToken(kinetic_time) + "_triangulateSimplePolygon_" + tag + "_"
    + std::to_string(debug_counter) + extension;

  const std::string branch_folder = runtime_branch_id.has_value()
    ? ("branch" + std::to_string(runtime_branch_id.value()))
    : kVisualDebugUnresolvedBranchFolder;
  std::filesystem::path filepath = std::filesystem::path(branch_folder) / filename;
  if (const std::optional<std::filesystem::path>& output_root = kin_del.getVisualDebugOutputRoot();
    output_root.has_value())
  {
    filepath = *output_root / filepath;
  }
  if (filepath.has_parent_path())
  {
    std::filesystem::create_directories(filepath.parent_path());
  }
  return filepath;
}

void appendTriangulateSimplePolygonVertexSourceFields(std::ostream& out, const std::string& metadata)
{
  const std::optional<std::string> source = metadataStringField(metadata, "source");
  out << " source=" << (source.has_value() ? source.value() : "unknown");

  if (const std::optional<std::string> event_type = metadataStringField(metadata, "event_type"); event_type.has_value())
  {
    out << " event_type=" << event_type.value();
  }

  // Associated IDs depend on source: site → strand_id; Voronoi vertex → voronoi_vertex_id;
  // intersection → delaunay_edge_id / voronoi_edge_id (and conceptual_* when present).
  if (const auto strand_id = metadataSizeField(metadata, "strand_id"); strand_id.has_value())
  {
    out << " strand_id=" << strand_id.value();
  }
  if (const auto voronoi_vertex_id = metadataSizeField(metadata, "voronoi_vertex_id"); voronoi_vertex_id.has_value())
  {
    out << " voronoi_vertex_id=" << voronoi_vertex_id.value();
  }
  if (const auto delaunay_edge_id = metadataSizeField(metadata, "delaunay_edge_id"); delaunay_edge_id.has_value())
  {
    out << " delaunay_edge_id=" << delaunay_edge_id.value();
  }
  if (const auto voronoi_edge_id = metadataSizeField(metadata, "voronoi_edge_id"); voronoi_edge_id.has_value())
  {
    out << " voronoi_edge_id=" << voronoi_edge_id.value();
  }
  if (const auto conceptual_delaunay_edge_id = metadataSizeField(metadata, "conceptual_delaunay_edge_id");
    conceptual_delaunay_edge_id.has_value())
  {
    out << " conceptual_delaunay_edge_id=" << conceptual_delaunay_edge_id.value();
  }
  if (const auto conceptual_voronoi_edge_id = metadataSizeField(metadata, "conceptual_voronoi_edge_id");
    conceptual_voronoi_edge_id.has_value())
  {
    out << " conceptual_voronoi_edge_id=" << conceptual_voronoi_edge_id.value();
  }
  if (const auto strand_cell_id = metadataSizeField(metadata, "strand_cell_id"); strand_cell_id.has_value())
  {
    out << " strand_cell_id=" << strand_cell_id.value();
  }
  if (const auto delaunay_face_id = metadataSizeField(metadata, "delaunay_face_id"); delaunay_face_id.has_value())
  {
    out << " delaunay_face_id=" << delaunay_face_id.value();
  }
}

void writeTriangulateSimplePolygonDebugTxt(const KineticDelaunay& kin_del, const VoronoiMesh& mesh,
  const std::filesystem::path& filepath, const char* tag,
  const std::vector<std::pair<std::string, std::vector<size_t>>>& rings, std::optional<double> occurrence_time,
  std::optional<size_t> runtime_branch_id)
{
  std::ofstream out(filepath);
  if (!out)
  {
    KINDS_WARNING("triangulateSimplePolygon: failed to open debug TXT " << filepath.generic_string());
    return;
  }

  const double kinetic_time = resolveTriangulateSimplePolygonDebugTime(mesh, occurrence_time);
  const auto& stored_vertices = mesh.getVertices();
  const auto& vertex_metadata = mesh.getVertexMetadata();

  out << "# tag=" << tag << '\n';
  out << "# occurrence_time=" << formatDebugExportTime(kinetic_time) << '\n';
  out << "# runtime_branch_id="
      << (runtime_branch_id.has_value() ? std::to_string(runtime_branch_id.value()) : "unresolved") << '\n';
  out << "# equality=exact profile_xy (triangulationPlaneXY), not mesh vertex id\n";
  out << "# ring_count=" << rings.size() << "\n\n";

  for (const auto& [ring_label, polygon_vertices] : rings)
  {
    out << "## ring=" << ring_label << " vertex_count=" << polygon_vertices.size() << '\n';

    for (size_t i = 0; i < polygon_vertices.size(); ++i)
    {
      for (size_t j = 0; j < i; ++j)
      {
        if (triangulationPlaneXYEqual(mesh, polygon_vertices[i], polygon_vertices[j]))
        {
          out << "# duplicate_xy: index=" << i << " matches index=" << j << " (ids " << polygon_vertices[i] << " vs "
              << polygon_vertices[j] << ")\n";
        }
      }
    }

    out << "# columns: index vertex_id profile_x profile_y object_x object_y object_z <source meta...>\n";
    for (size_t i = 0; i < polygon_vertices.size(); ++i)
    {
      const size_t vertex_id = polygon_vertices[i];
      const glm::dvec2 profile = mesh.triangulationPlaneXY(vertex_id);
      const glm::dvec3 object = vertex_id < stored_vertices.size() ? stored_vertices[vertex_id] : glm::dvec3(0.0);
      out << i << ' ' << vertex_id << ' ' << numberLiteral(profile.x) << ' ' << numberLiteral(profile.y) << ' '
          << numberLiteral(object.x) << ' ' << numberLiteral(object.y) << ' ' << numberLiteral(object.z);
      if (vertex_id < vertex_metadata.size())
      {
        appendTriangulateSimplePolygonVertexSourceFields(out, vertex_metadata[vertex_id]);
      }
      else
      {
        out << " source=unknown";
      }
      out << '\n';
    }
    out << '\n';
  }

  KINDS_WARNING("triangulateSimplePolygon: wrote debug TXT to " << filepath.generic_string());
}

std::vector<std::vector<size_t>> splitPolygonAtRepeatedVertices(
  const VoronoiMesh& mesh, const std::vector<size_t>& polygon)
{
  std::vector<std::vector<size_t>> worklist;
  worklist.push_back(polygon);
  std::vector<std::vector<size_t>> simple_rings;

  while (!worklist.empty())
  {
    std::vector<size_t> ring = std::move(worklist.back());
    worklist.pop_back();
    if (ring.size() < 3)
    {
      continue;
    }

    bool split_any = false;
    for (size_t i = 0; i < ring.size() && !split_any; ++i)
    {
      for (size_t j = 0; j < i; ++j)
      {
        if (!triangulationPlaneXYEqual(mesh, ring[i], ring[j]))
        {
          continue;
        }

        const size_t start = j;
        const size_t end = i;

        // Inner loop: [start, end) — shared pinch point only at the start index.
        std::vector<size_t> sub1(
          ring.begin() + static_cast<std::ptrdiff_t>(start), ring.begin() + static_cast<std::ptrdiff_t>(end));

        // Outer loop: [end, n) + [0, start) — shared pinch point only at sub2[0].
        std::vector<size_t> sub2;
        sub2.reserve(ring.size() - (end - start));
        sub2.insert(sub2.end(), ring.begin() + static_cast<std::ptrdiff_t>(end), ring.end());
        sub2.insert(sub2.end(), ring.begin(), ring.begin() + static_cast<std::ptrdiff_t>(start));

        if (sub1.size() >= 3)
        {
          worklist.push_back(std::move(sub1));
        }
        if (sub2.size() >= 3)
        {
          worklist.push_back(std::move(sub2));
        }
        split_any = true;
        break;
      }
    }

    if (!split_any)
    {
      simple_rings.push_back(std::move(ring));
    }
  }

  return simple_rings;
}

void writeTriangulateSimplePolygonFailSvg(KineticDelaunay& kin_del, VoronoiMesh& mesh,
  const std::vector<size_t>& polygon_vertices, std::optional<double> occurrence_time,
  std::optional<size_t> runtime_branch_id)
{
  if (polygon_vertices.size() < 3)
  {
    return;
  }

  const double kinetic_time = resolveTriangulateSimplePolygonDebugTime(mesh, occurrence_time);

  // Common event-style kinetic SVG with sites / intersections / VVs from the failed polygon highlighted.
  {
    const auto& graph = kin_del.getGraph();
    VisualDebugHighlight highlight;
    const auto& vertex_metadata = mesh.getVertexMetadata();
    std::optional<size_t> delaunay_face_id;
    for (size_t vertex_id : polygon_vertices)
    {
      if (vertex_id >= vertex_metadata.size())
      {
        continue;
      }
      const std::string& metadata = vertex_metadata[vertex_id];
      if (const auto strand_id = metadataSizeField(metadata, "strand_id"); strand_id.has_value())
      {
        highlight.delaunay_vertices.insert(strand_id.value());
      }
      if (const auto voronoi_vertex_id = metadataSizeField(metadata, "voronoi_vertex_id");
        voronoi_vertex_id.has_value())
      {
        highlight.voronoi_vertices.insert(voronoi_vertex_id.value());
      }
      if (const auto delaunay_edge_id = metadataSizeField(metadata, "delaunay_edge_id"); delaunay_edge_id.has_value())
      {
        const size_t de = delaunay_edge_id.value();
        highlight.label_crossings_on_delaunay_edges.insert(de);
        highlight.addUndirectedDelaunayEdge(graph, 2 * de);
      }
      if (const auto voronoi_edge_id = metadataSizeField(metadata, "voronoi_edge_id"); voronoi_edge_id.has_value())
      {
        const size_t ve = voronoi_edge_id.value();
        highlight.voronoi_edges.insert(ve);
        highlight.label_crossings_on_voronoi_edges.insert(ve);
      }
      if (const auto face_id = metadataSizeField(metadata, "delaunay_face_id"); face_id.has_value())
      {
        delaunay_face_id = face_id;
      }
      const auto de = metadataSizeField(metadata, "delaunay_edge_id");
      const auto ve = metadataSizeField(metadata, "voronoi_edge_id");
      if (de.has_value() && ve.has_value())
      {
        highlight.crossing_intersection_keys.insert((static_cast<uint64_t>(de.value()) << 32) | ve.value());
      }
    }
    if (delaunay_face_id.has_value() && graph.isLiveFace(delaunay_face_id.value()))
    {
      highlight.addDelaunayTriangle(graph, delaunay_face_id.value());
    }
    writeSegmentBuilderVisualDebugSvg(
      true, kin_del, graph, kinetic_time, "error", "triangulateSimplePolygon_FAIL", highlight, runtime_branch_id);
  }

  const std::filesystem::path filepath
    = makeTriangulateSimplePolygonDebugPath(kin_del, mesh, "FAIL", ".svg", occurrence_time, runtime_branch_id);

  std::vector<glm::dvec2> ring_xy;
  ring_xy.reserve(polygon_vertices.size());
  double min_x = std::numeric_limits<double>::infinity();
  double min_y = std::numeric_limits<double>::infinity();
  double max_x = -std::numeric_limits<double>::infinity();
  double max_y = -std::numeric_limits<double>::infinity();
  for (size_t vertex_id : polygon_vertices)
  {
    const glm::dvec2 xy = mesh.triangulationPlaneXY(vertex_id);
    ring_xy.push_back(xy);
    min_x = std::min(min_x, xy.x);
    min_y = std::min(min_y, xy.y);
    max_x = std::max(max_x, xy.x);
    max_y = std::max(max_y, xy.y);
  }

  constexpr double pad = 0.05;
  const double span_x = std::max(max_x - min_x, 1e-6);
  const double span_y = std::max(max_y - min_y, 1e-6);
  min_x -= pad * span_x;
  min_y -= pad * span_y;
  max_x += pad * span_x;
  max_y += pad * span_y;
  const double width = max_x - min_x;
  const double height = max_y - min_y;

  auto svg_x = [&](double x) { return x - min_x; };
  auto svg_y = [&](double y) { return max_y - y; };

  std::ostringstream poly_points;
  for (size_t i = 0; i < ring_xy.size(); ++i)
  {
    if (i > 0)
    {
      poly_points << ' ';
    }
    poly_points << svg_x(ring_xy[i].x) << ',' << svg_y(ring_xy[i].y);
  }

  std::ofstream out(filepath);
  if (!out)
  {
    KINDS_WARNING("triangulateSimplePolygon: failed to open debug SVG " << filepath.generic_string());
    return;
  }

  out << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n";
  out << "<svg xmlns=\"http://www.w3.org/2000/svg\" viewBox=\"0 0 " << width << ' ' << height << "\" width=\""
      << width * 100.0 << "\" height=\"" << height * 100.0 << "\">\n";
  out << "<rect x=\"0\" y=\"0\" width=\"" << width << "\" height=\"" << height << "\" fill=\"#f8f8f8\"/>\n";
  out << "<polygon points=\"" << poly_points.str()
      << "\" fill=\"rgba(255,80,80,0.15)\" stroke=\"#cc0000\" stroke-width=\"" << span_x * 0.002
      << "\" fill-rule=\"evenodd\"/>\n";

  const auto& stored_vertices = mesh.getVertices();
  const auto& vertex_metadata = mesh.getVertexMetadata();
  constexpr double eps = 1e-12;
  for (size_t i = 0; i < polygon_vertices.size(); ++i)
  {
    const size_t vertex_id = polygon_vertices[i];
    const size_t prev_id = polygon_vertices[(i + polygon_vertices.size() - 1) % polygon_vertices.size()];
    const size_t next_id = polygon_vertices[(i + 1) % polygon_vertices.size()];
    const glm::dvec2 a = mesh.triangulationPlaneXY(prev_id);
    const glm::dvec2 b = ring_xy[i];
    const glm::dvec2 c = mesh.triangulationPlaneXY(next_id);
    const double turn = glm::cross(b - a, c - b);
    const bool reflex = turn <= eps;
    const bool finite_xy = std::isfinite(b.x) && std::isfinite(b.y);
    const glm::dvec3 stored = vertex_id < stored_vertices.size() ? stored_vertices[vertex_id] : glm::dvec3(0.0);
    const bool finite_stored = std::isfinite(stored.x) && std::isfinite(stored.y) && std::isfinite(stored.z);

    const double cx = svg_x(b.x);
    const double cy = svg_y(b.y);
    const double r = span_x * 0.008;
    const char* fill = !finite_xy || !finite_stored ? "#ff00ff" : (reflex ? "#ffaa00" : "#00aa00");
    out << "<circle cx=\"" << cx << "\" cy=\"" << cy << "\" r=\"" << r << "\" fill=\"" << fill
        << "\" stroke=\"#222\" stroke-width=\"" << span_x * 0.0008 << "\"/>\n";

    std::ostringstream label;
    label << "i=" << i << " id=" << vertex_id << " (" << numberLiteral(b.x) << "," << numberLiteral(b.y) << ")";
    if (!finite_xy || !finite_stored)
    {
      label << " NON_FINITE";
    }
    if (reflex)
    {
      label << " reflex";
    }
    if (vertex_id < vertex_metadata.size())
    {
      const std::string& metadata = vertex_metadata[vertex_id];
      if (const std::optional<std::string> source = metadataStringField(metadata, "source"); source.has_value())
      {
        label << " " << source.value();
      }
      if (const auto voronoi_vertex_id = metadataSizeField(metadata, "voronoi_vertex_id");
        voronoi_vertex_id.has_value())
      {
        label << " vv=" << voronoi_vertex_id.value();
      }
      if (const auto strand_id = metadataSizeField(metadata, "strand_id"); strand_id.has_value())
      {
        label << " strand=" << strand_id.value();
      }
      if (const auto delaunay_edge_id = metadataSizeField(metadata, "delaunay_edge_id"); delaunay_edge_id.has_value())
      {
        label << " de=" << delaunay_edge_id.value();
      }
      if (const auto voronoi_edge_id = metadataSizeField(metadata, "voronoi_edge_id"); voronoi_edge_id.has_value())
      {
        label << " ve=" << voronoi_edge_id.value();
      }
    }
    out << "<text x=\"" << (cx + r * 1.5) << "\" y=\"" << (cy - r * 1.5) << "\" font-size=\"" << span_x * 0.02
        << "\" fill=\"#111\">" << svgEscape(label.str()) << "</text>\n";
  }

  out << "<text x=\"" << span_x * 0.02 << "\" y=\"" << span_y * 0.05 << "\" font-size=\"" << span_x * 0.025
      << "\" fill=\"#111\">triangulateSimplePolygon ear-clip failure (" << polygon_vertices.size()
      << " vertices)</text>\n";
  out << "</svg>\n";

  KINDS_WARNING("triangulateSimplePolygon: wrote debug SVG to " << filepath.generic_string());
}

std::optional<std::array<double, 3>> barycentricCoordinates2D(
  const glm::dvec2& a, const glm::dvec2& b, const glm::dvec2& c, const glm::dvec2& p)
{
  const glm::dvec2 v0 = b - a;
  const glm::dvec2 v1 = c - a;
  const glm::dvec2 v2 = p - a;

  const double d00 = glm::dot(v0, v0);
  const double d01 = glm::dot(v0, v1);
  const double d11 = glm::dot(v1, v1);
  const double d20 = glm::dot(v2, v0);
  const double d21 = glm::dot(v2, v1);
  const double denom = d00 * d11 - d01 * d01;
  const double scale = std::max(d00, d11);
  if (std::abs(denom) <= 1e-24 * scale * scale)
  {
    return std::nullopt;
  }

  const double v = (d11 * d20 - d01 * d21) / denom;
  const double w = (d00 * d21 - d01 * d20) / denom;
  const double u = 1.0 - v - w;
  return std::array<double, 3> { u, v, w };
}

std::optional<size_t> metadataSizeField(const std::string& metadata, const char* key)
{
  const std::string needle = jsonStringLiteral(key) + ":";
  const size_t key_pos = metadata.find(needle);
  if (key_pos == std::string::npos)
  {
    return std::nullopt;
  }
  size_t value_pos = key_pos + needle.size();
  while (value_pos < metadata.size() && std::isspace(static_cast<unsigned char>(metadata[value_pos])))
  {
    ++value_pos;
  }
  if (value_pos >= metadata.size() || !std::isdigit(static_cast<unsigned char>(metadata[value_pos])))
  {
    return std::nullopt;
  }
  size_t value = 0;
  for (; value_pos < metadata.size() && std::isdigit(static_cast<unsigned char>(metadata[value_pos])); ++value_pos)
  {
    value = value * 10 + static_cast<size_t>(metadata[value_pos] - '0');
  }
  return value;
}

std::optional<double> metadataDoubleField(const std::string& metadata, const char* key)
{
  const std::string needle = jsonStringLiteral(key) + ":";
  const size_t key_pos = metadata.find(needle);
  if (key_pos == std::string::npos)
  {
    return std::nullopt;
  }
  size_t value_pos = key_pos + needle.size();
  while (value_pos < metadata.size() && std::isspace(static_cast<unsigned char>(metadata[value_pos])))
  {
    ++value_pos;
  }
  if (value_pos >= metadata.size())
  {
    return std::nullopt;
  }

  const char* start = metadata.c_str() + value_pos;
  char* end = nullptr;
  const double value = std::strtod(start, &end);
  if (end == start)
  {
    return std::nullopt;
  }
  return value;
}

bool isKnownVertexMetadataSource(const std::string& source)
{
  return source == "Voronoi vertex" || source == "intersection" || source == "site";
}
} // namespace

SegmentBuilder::MetadataBuilder SegmentBuilder::MetadataBuilder::fromObject(const std::string& metadata)
{
  MetadataBuilder builder;
  const size_t begin = metadata.find('{');
  const size_t end = metadata.rfind('}');
  if (begin != std::string::npos && end != std::string::npos && end > begin + 1)
  {
    builder.raw_body_ = metadata.substr(begin + 1, end - begin - 1);
  }
  return builder;
}

SegmentBuilder::MetadataBuilder& SegmentBuilder::MetadataBuilder::addString(const char* key, const char* value)
{
  return addRaw(key, jsonStringLiteral(value != nullptr ? std::string(value) : std::string()));
}

SegmentBuilder::MetadataBuilder& SegmentBuilder::MetadataBuilder::addString(const char* key, const std::string& value)
{
  return addRaw(key, jsonStringLiteral(value));
}

SegmentBuilder::MetadataBuilder& SegmentBuilder::MetadataBuilder::addSize(const char* key, size_t value)
{
  return addRaw(key, std::to_string(value));
}

SegmentBuilder::MetadataBuilder& SegmentBuilder::MetadataBuilder::addInt(const char* key, int value)
{
  return addRaw(key, std::to_string(value));
}

SegmentBuilder::MetadataBuilder& SegmentBuilder::MetadataBuilder::addDouble(const char* key, double value)
{
  return addRaw(key, numberLiteral(value));
}

SegmentBuilder::MetadataBuilder& SegmentBuilder::MetadataBuilder::addBool(const char* key, bool value)
{
  return addRaw(key, value ? "true" : "false");
}

SegmentBuilder::MetadataBuilder& SegmentBuilder::MetadataBuilder::addRaw(
  const char* key, const std::string& raw_json_value)
{
  fields_.emplace_back(key != nullptr ? std::string(key) : std::string(), raw_json_value);
  return *this;
}

SegmentBuilder::MetadataBuilder& SegmentBuilder::MetadataBuilder::ensureString(const char* key, const char* value)
{
  if (!hasKey(key))
  {
    addString(key, value);
  }
  return *this;
}

SegmentBuilder::MetadataBuilder& SegmentBuilder::MetadataBuilder::ensureBool(const char* key, bool value)
{
  if (!hasKey(key))
  {
    addBool(key, value);
  }
  return *this;
}

SegmentBuilder::MetadataBuilder& SegmentBuilder::MetadataBuilder::ensureSize(const char* key, size_t value)
{
  if (!hasKey(key))
  {
    addSize(key, value);
  }
  return *this;
}

SegmentBuilder::MetadataBuilder& SegmentBuilder::MetadataBuilder::ensureDouble(const char* key, double value)
{
  if (!hasKey(key))
  {
    addDouble(key, value);
  }
  return *this;
}

std::string SegmentBuilder::MetadataBuilder::build() const
{
  std::string out = "{";
  bool needs_comma = false;
  if (!raw_body_.empty())
  {
    out += raw_body_;
    needs_comma = true;
  }
  for (const auto& [key, value] : fields_)
  {
    if (needs_comma)
    {
      out += ",";
    }
    out += jsonStringLiteral(key);
    out += ":";
    out += value;
    needs_comma = true;
  }
  out += "}";
  return out;
}

bool SegmentBuilder::MetadataBuilder::hasKey(const char* key) const
{
  const std::string needle = jsonStringLiteral(key != nullptr ? std::string(key) : std::string()) + ":";
  if (raw_body_.find(needle) != std::string::npos)
  {
    return true;
  }
  return std::any_of(fields_.begin(), fields_.end(), [&](const auto& field) { return field.first == key; });
}

SegmentBuilder::ScopedMetadataCallbackPhase::ScopedMetadataCallbackPhase(
  SegmentBuilder& segment_builder, const char* callback_phase)
  : segment_builder_(segment_builder)
  , previous_phase_(segment_builder.metadata_callback_phase_)
{
  segment_builder_.metadata_callback_phase_ = callback_phase != nullptr ? callback_phase : std::string {};
}

SegmentBuilder::ScopedMetadataCallbackPhase::~ScopedMetadataCallbackPhase()
{
  segment_builder_.metadata_callback_phase_ = std::move(previous_phase_);
}

SegmentBuilder::ScopedEndActiveCrossingEvent::ScopedEndActiveCrossingEvent(SegmentBuilder& segment_builder)
  : segment_builder_(segment_builder)
{
}

SegmentBuilder::ScopedEndActiveCrossingEvent::~ScopedEndActiveCrossingEvent()
{
  segment_builder_.endActiveCrossingEvent();
}

void SegmentBuilder::endActiveCrossingEvent()
{
  active_crossing_voronoi_vertex_id_.reset();
  active_crossing_delaunay_edge_id_.reset();
  active_crossing_voronoi_edge_ids_ = { static_cast<size_t>(-1), static_cast<size_t>(-1), static_cast<size_t>(-1) };
  active_crossing_mesh_position_.reset();
  active_crossing_delaunay_xy_.reset();
}

bool SegmentBuilder::activeCrossingEventBufferApplies(std::optional<size_t> voronoi_vertex_id,
  const std::optional<KineticDelaunay::CrossingData::EdgeIntersectionRef>& intersection) const
{
  if (!active_crossing_mesh_position_.has_value())
  {
    return false;
  }
  if (voronoi_vertex_id.has_value() && active_crossing_voronoi_vertex_id_.has_value()
    && voronoi_vertex_id.value() == active_crossing_voronoi_vertex_id_.value())
  {
    return true;
  }
  if (!intersection.has_value() || !active_crossing_delaunay_edge_id_.has_value())
  {
    return false;
  }
  if (intersection.value()->delaunay_edge_id != active_crossing_delaunay_edge_id_.value())
  {
    return false;
  }
  const size_t voronoi_edge_id = intersection.value()->voronoi_edge_id;
  for (size_t incident_edge_id : active_crossing_voronoi_edge_ids_)
  {
    if (incident_edge_id == voronoi_edge_id)
    {
      return true;
    }
  }
  return false;
}

namespace
{
bool vertexPositionFinite(const glm::dvec3& p);
}

void SegmentBuilder::beginActiveCrossingEvent(
  size_t voronoi_vertex_id, size_t delaunay_edge_id, double t, size_t strand_id)
{
  endActiveCrossingEvent();
  active_crossing_voronoi_vertex_id_ = voronoi_vertex_id;
  active_crossing_delaunay_edge_id_ = delaunay_edge_id;

  const auto& graph = kin_del.getGraph();
  const glm::dvec3 fallback_profile(0.0, 0.0, t);
  if (voronoi_vertex_id < graph.faceSlotCount() && graph.isLiveFace(voronoi_vertex_id))
  {
    const std::array<size_t, 3> half_edges = graph.face(voronoi_vertex_id).half_edges;
    for (size_t i = 0; i < 3; ++i)
    {
      active_crossing_voronoi_edge_ids_[i] = half_edges[i] / 2;
    }
    active_crossing_delaunay_xy_ = glm::dvec2(computeVoronoiVertex(half_edges[0], t));
  }

  // Prefer one pre-event incident intersection so the shared point uses Delaunay-edge interpolation.
  const auto coincidence = findIncidentIntersectionOnDelaunayEdge(voronoi_vertex_id, delaunay_edge_id);
  if (coincidence.has_value())
  {
    MeshletVertexRuntimeInfo runtime_info;
    runtime_info.position_intersection = coincidence;
    runtime_info.conceptual_intersection = coincidence;
    const MeshIntersectionObjectSpaceResult inter
      = computeMeshIntersectionObjectSpace(runtime_info, fallback_profile, strand_id, t);
    if (vertexPositionFinite(inter.position))
    {
      active_crossing_mesh_position_ = inter.position;
      active_crossing_delaunay_xy_
        = glm::dvec2(getCrossingCoordsInDelaunaySpace(kin_del, coincidence.value(), t));
      KINDS_DEBUG("beginActiveCrossingEvent: buffered vv=" << voronoi_vertex_id << " de=" << delaunay_edge_id
        << " from intersection ve=" << coincidence.value()->voronoi_edge_id
        << " d_param=" << coincidence.value()->delaunay_edge_param << " t=" << t);
      return;
    }
    KINDS_WARNING("beginActiveCrossingEvent: incident intersection placement was non-finite for Voronoi vertex "
      << voronoi_vertex_id << " at t=" << t << "; falling back to barycentric.");
  }

  const MeshVoronoiVertexObjectSpaceResult voronoi
    = computeMeshVoronoiVertexObjectSpace(voronoi_vertex_id, fallback_profile, strand_id, t);
  if (vertexPositionFinite(voronoi.position))
  {
    active_crossing_mesh_position_ = voronoi.position;
    KINDS_DEBUG("beginActiveCrossingEvent: buffered vv=" << voronoi_vertex_id << " de=" << delaunay_edge_id
      << " from barycentric Voronoi vertex t=" << t);
    return;
  }
  KINDS_WARNING("beginActiveCrossingEvent: failed to buffer a finite crossing vertex for Voronoi vertex "
    << voronoi_vertex_id << " at t=" << t);
}

namespace
{
std::optional<glm::dvec3> meshVertexUv(const VoronoiMesh& mesh, size_t vertex_index)
{
  if (const std::optional<glm::dvec3> semantic_uv = mesh.vertexSemanticUv(vertex_index); semantic_uv.has_value())
  {
    return semantic_uv;
  }
  const std::vector<glm::dvec3>& uvs = mesh.getUVs();
  if (vertex_index < uvs.size())
  {
    const glm::dvec3& pooled_uv = uvs[vertex_index];
    if (std::isfinite(pooled_uv.x) && std::isfinite(pooled_uv.y) && std::isfinite(pooled_uv.z))
    {
      return pooled_uv;
    }
  }
  for (size_t triangle_corner : mesh.findTriangleCorners(vertex_index, true))
  {
    if (mesh.hasValidUVIndex(triangle_corner))
    {
      return mesh.getUV(triangle_corner);
    }
  }
  return std::nullopt;
}

void warnIfTriangleKineticTimesNotInUnitSection(
  size_t v1, size_t v2, size_t v3, const std::vector<glm::dvec3>& vertices, const char* mesh_kind, int material_id)
{
  const double t0 = vertices[v1][2];
  const double t1 = vertices[v2][2];
  const double t2 = vertices[v3][2];
  const double min_t = std::min({ t0, t1, t2 });
  const double max_t = std::max({ t0, t1, t2 });
  const double n = std::floor(min_t);
  constexpr double eps = 1e-9;
  if (max_t > n + 1.0 + eps)
  {
    KINDS_WARNING("[" << mesh_kind << "] Triangle kinetic times not contained in any unit section [n,n+1]: vertices=("
                      << v1 << "," << v2 << "," << v3 << ") times=(" << t0 << "," << t1 << "," << t2 << ") min="
                      << min_t << " max=" << max_t << " floor(min)=" << n << " material_id=" << material_id);
  }
}

void setMeshVertexUv(VoronoiMesh& mesh, size_t vertex_index, const glm::dvec3& uv)
{
  mesh.setVertexSemanticUv(vertex_index, uv);

  for (size_t triangle_corner : mesh.findTriangleCorners(vertex_index))
  {
    if (mesh.hasValidUVIndex(triangle_corner))
    {
      mesh.setUV(uv, triangle_corner);
    }
    else if (triangle_corner < mesh.getUVIndices().size())
    {
      mesh.getUVIndices()[triangle_corner] = mesh.addUV(uv);
    }
  }
}

struct AdjustedBoundaryTriangleUvs
{
  glm::dvec3 u {};
  glm::dvec3 v {};
  glm::dvec3 w {};
};

AdjustedBoundaryTriangleUvs adjustedBoundaryTriangleUvs(const glm::dvec2& raw_u, const glm::dvec2& raw_v,
  const glm::dvec2& raw_w, double uv_circum_factor, double uv_height_factor)
{
  AdjustedBoundaryTriangleUvs adjusted;
  adjusted.u = glm::dvec3(raw_u.x, raw_u.y, 0.0);
  adjusted.v = glm::dvec3(raw_v.x, raw_v.y, 0.0);
  adjusted.w = glm::dvec3(raw_w.x, raw_w.y, 0.0);

  const double base_angle = adjusted.u.x;
  double& angle_v = adjusted.v.x;
  angle_v -= std::round(angle_v - base_angle);
  double& angle_w = adjusted.w.x;
  angle_w -= std::round(angle_w - base_angle);

  adjusted.u.x *= uv_circum_factor;
  adjusted.v.x *= uv_circum_factor;
  adjusted.w.x *= uv_circum_factor;
  adjusted.u.y *= uv_height_factor;
  adjusted.v.y *= uv_height_factor;
  adjusted.w.y *= uv_height_factor;
  return adjusted;
}

bool vertexPositionFinite(const glm::dvec3& p)
{
  return std::isfinite(p.x) && std::isfinite(p.y) && std::isfinite(p.z);
}

void syncResolvedFlexibleVertexMetadata(VoronoiMesh& mesh, size_t vertex_index)
{
  if (!mesh.storeMetadata() || vertex_index >= mesh.getVertices().size())
  {
    return;
  }
  const glm::dvec2 profile = mesh.triangulationPlaneXY(vertex_index);
  const double kinetic_time = mesh.vertexKineticTime(vertex_index);
  std::string existing = "{}";
  if (vertex_index < mesh.getVertexMetadata().size())
  {
    existing = mesh.getVertexMetadata()[vertex_index];
  }
  const std::string parent_event_type = metadataStringField(existing, "parent_event_type").value_or("unknown_event");
  const std::string pos = metadataStringField(existing, "pos").value_or("unknown");
  SegmentBuilder::MetadataBuilder builder;
  builder.addString("event_type", parent_event_type)
    .addString("source", "site")
    .addString("pos", pos)
    .addBool("intersection_flexible_placeholder", false)
    .addBool("shift", false)
    .addDouble("x", profile.x)
    .addDouble("y", profile.y)
    .addDouble("t", kinetic_time);
  if (const auto strand_id = metadataSizeField(existing, "strand_id"); strand_id.has_value())
  {
    builder.addSize("strand_id", strand_id.value());
  }
  mesh.setVertexMetadata(vertex_index, builder.build());
}

const char* boundaryEventTypeToString(SegmentBuilder::BoundaryEventType event_type)
{
  switch (event_type)
  {
  case SegmentBuilder::BoundaryEventType::Init:
    return "init";
  case SegmentBuilder::BoundaryEventType::Section:
    return "section_event";
  case SegmentBuilder::BoundaryEventType::Radius:
    return "radius_event";
  case SegmentBuilder::BoundaryEventType::Crossing:
    return "crossing_event";
  case SegmentBuilder::BoundaryEventType::Subdivision:
    return "subdivision_event";
  case SegmentBuilder::BoundaryEventType::Separation:
    return "separation_event";
  default:
    return "unknown_event";
  }
}

const char* boundarySegmentActionToString(SegmentBuilder::BoundarySegmentAction segment_action)
{
  switch (segment_action)
  {
  case SegmentBuilder::BoundarySegmentAction::NewSegment:
    return "new_segment";
  case SegmentBuilder::BoundarySegmentAction::SegmentCompleted:
    return "segment_completed";
  case SegmentBuilder::BoundarySegmentAction::SegmentRemoved:
    return "segment_removed";
  case SegmentBuilder::BoundarySegmentAction::SegmentRemapped:
    return "segment_remapped";
  default:
    return "unknown_action";
  }
}

std::string makeBoundaryMeshMetadata(
  SegmentBuilder::BoundaryEventType event_type, SegmentBuilder::BoundarySegmentAction segment_action)
{
  return SegmentBuilder::MetadataBuilder()
    .addString("event_type", boundaryEventTypeToString(event_type))
    .addString("segment_action", boundarySegmentActionToString(segment_action))
    .build();
}

std::string makeClosingMeshVertexMetadata(const char* source,
  std::optional<KineticDelaunay::CrossingData::EdgeIntersectionRef> crossing = std::nullopt,
  std::optional<size_t> voronoi_vertex_id = std::nullopt, std::optional<size_t> strand_id = std::nullopt)
{
  SegmentBuilder::MetadataBuilder builder;
  builder.addString("event_type", "closing_mesh").addString("source", source);
  const std::string source_string(source);
  if (source_string == "intersection" && crossing.has_value())
  {
    builder.addSize("delaunay_edge_id", crossing.value()->delaunay_edge_id)
      .addSize("voronoi_edge_id", crossing.value()->voronoi_edge_id)
      .addDouble("stored_d_param", crossing.value()->delaunay_edge_param);
  }
  else if (source_string == "Voronoi vertex" && voronoi_vertex_id.has_value())
  {
    builder.addSize("voronoi_vertex_id", voronoi_vertex_id.value());
  }
  else if (source_string == "site" && strand_id.has_value())
  {
    builder.addSize("strand_id", strand_id.value());
  }
  return builder.build();
}

std::optional<std::pair<size_t, size_t>> closingMeshIntersectionListIndices(
  const KineticDelaunay& kd, KineticDelaunay::CrossingData::EdgeIntersectionRef ref)
{
  const auto& crossing_data = kd.getCrossingData();
  const auto& ir = *ref;
  const size_t d_edge = ir.delaunay_edge_id;
  const size_t v_edge = ir.voronoi_edge_id;
  if (d_edge >= crossing_data.delaunay_edge_intersections.size()
    || v_edge >= crossing_data.voronoi_edge_intersections.size())
  {
    return std::nullopt;
  }
  size_t d_list_idx = 0;
  for (auto d_it = crossing_data.delaunay_edge_intersections[d_edge].begin();
    d_it != crossing_data.delaunay_edge_intersections[d_edge].end(); ++d_it, ++d_list_idx)
  {
    const KineticDelaunay::CrossingData::EdgeIntersectionRef r = *d_it;
    if (&(*r) != &ir)
    {
      continue;
    }
    size_t v_list_idx = 0;
    for (auto v_it = crossing_data.voronoi_edge_intersections[v_edge].begin();
      v_it != crossing_data.voronoi_edge_intersections[v_edge].end(); ++v_it, ++v_list_idx)
    {
      if (&(**v_it) == &ir)
      {
        return std::make_pair(d_list_idx, v_list_idx);
      }
    }
    break;
  }
  return std::nullopt;
}

void logClosingMeshIntersectionVertex(const KineticDelaunay& kd, KineticDelaunay::CrossingData::EdgeIntersectionRef ref)
{
  const auto& ir = *ref;
  const auto idx = closingMeshIntersectionListIndices(kd, ref);
  if (idx.has_value())
  {
    KINDS_DEBUG("Adding vertex at intersection of Delaunay edge "
      << ir.delaunay_edge_id << " and Voronoi edge " << ir.voronoi_edge_id << " with intersection indices ("
      << idx->first << "," << idx->second << ")");
  }
  else
  {
    KINDS_DEBUG("Adding vertex at intersection of Delaunay edge "
      << ir.delaunay_edge_id << " and Voronoi edge " << ir.voronoi_edge_id << " with intersection indices (?,?)");
  }
}

void logClosingMeshStripEndpointVertex(const KineticDelaunay& kd, size_t segment_index,
  const std::optional<KineticDelaunay::CrossingData::EdgeIntersectionRef>& endpoint_crossing,
  size_t closing_voronoi_edge_id)
{
  if (endpoint_crossing.has_value())
  {
    logClosingMeshIntersectionVertex(kd, endpoint_crossing.value());
  }
  else if (closing_voronoi_edge_id == static_cast<size_t>(-1))
  {
    KINDS_DEBUG(
      "Polygon walk: using mesh vertex on segment " << segment_index << " of Voronoi edge ? (strip circumcenter).");
  }
  else
  {
    KINDS_DEBUG("Polygon walk: using mesh vertex on segment " << segment_index << " of Voronoi edge "
                                                              << closing_voronoi_edge_id << " (strip circumcenter).");
  }
}

void logExtractionClosingMeshStripEndpoint(const KineticDelaunay& kd, double t, size_t extraction_index,
  const std::optional<KineticDelaunay::CrossingData::EdgeIntersectionRef>& endpoint_crossing,
  size_t closing_voronoi_edge_id)
{
  if (endpoint_crossing.has_value())
  {
    const auto ref = endpoint_crossing.value();
    const auto& ir = *ref;
    const auto idx = closingMeshIntersectionListIndices(kd, ref);
    glm::dvec2 xy(0.0, 0.0);
    const bool have_xy = tryComputeCrossingIntersectionPosition2D(kd, endpoint_crossing, t, xy, false, false);
    if (idx.has_value())
    {
      if (have_xy)
      {
        KINDS_DEBUG("Closing mesh extraction: strip endpoint at Delaunay edge "
          << ir.delaunay_edge_id << " / Voronoi edge " << ir.voronoi_edge_id << ", list indices (" << idx->first << ","
          << idx->second << "), 2D=(" << xy.x << "," << xy.y << ")");
      }
      else
      {
        KINDS_DEBUG("Closing mesh extraction: strip endpoint at Delaunay edge "
          << ir.delaunay_edge_id << " / Voronoi edge " << ir.voronoi_edge_id << ", list indices (" << idx->first << ","
          << idx->second << ")");
      }
    }
    else if (have_xy)
    {
      KINDS_DEBUG("Closing mesh extraction: strip endpoint at Delaunay edge "
        << ir.delaunay_edge_id << " / Voronoi edge " << ir.voronoi_edge_id << ", list indices (?,?), 2D=(" << xy.x
        << "," << xy.y << ")");
    }
    else
    {
      KINDS_DEBUG("Closing mesh extraction: strip endpoint at Delaunay edge "
        << ir.delaunay_edge_id << " / Voronoi edge " << ir.voronoi_edge_id << ", list indices (?,?)");
    }
  }
  else if (closing_voronoi_edge_id == static_cast<size_t>(-1))
  {
    KINDS_DEBUG("Closing mesh extraction: strip endpoint on extraction segment " << extraction_index
                                                                                 << " of Voronoi edge ? (circumcenter, "
                                                                                    "no crossing ref).");
  }
  else
  {
    KINDS_DEBUG("Closing mesh extraction: strip endpoint on extraction segment "
      << extraction_index << " of Voronoi edge " << closing_voronoi_edge_id << " (circumcenter, no crossing ref).");
  }
}

std::string makeRegularStripVertexMetadataJson(double kinetic_time, size_t voronoi_edge_id, size_t even_half_edge_id,
  int strand_even_origin, int strand_odd_origin, SegmentBuilder::BoundaryEventType event_type,
  SegmentBuilder::BoundarySegmentAction segment_action,
  const std::optional<KineticDelaunay::CrossingData::EdgeIntersectionRef>& crossing, const char* pos, const char* op,
  const char* source, const std::optional<KineticDelaunay::CrossingData::EdgeIntersectionRef>& position_crossing,
  bool radius_shift_explicit_profile_position)
{
  (void)kinetic_time;
  (void)voronoi_edge_id;
  (void)even_half_edge_id;
  (void)strand_even_origin;
  (void)strand_odd_origin;
  (void)segment_action;
  (void)pos;
  (void)op;
  SegmentBuilder::MetadataBuilder builder;
  builder.addString("event_type", boundaryEventTypeToString(event_type))
    .addString("source", source != nullptr ? source : (crossing.has_value() ? "intersection" : "site"));
  if (crossing.has_value() || position_crossing.has_value())
  {
    const KineticDelaunay::CrossingData::EdgeIntersectionRef geom_ref
      = position_crossing.has_value() ? position_crossing.value() : crossing.value();
    const KineticDelaunay::CrossingData::EdgeIntersectionRef conceptual_ref
      = crossing.has_value() ? crossing.value() : geom_ref;
    builder.addSize("delaunay_edge_id", geom_ref->delaunay_edge_id)
      .addSize("voronoi_edge_id", geom_ref->voronoi_edge_id)
      .addDouble("stored_d_param", geom_ref->delaunay_edge_param);
    if (conceptual_ref != geom_ref)
    {
      builder.addSize("conceptual_delaunay_edge_id", conceptual_ref->delaunay_edge_id)
        .addSize("conceptual_voronoi_edge_id", conceptual_ref->voronoi_edge_id);
    }
    if (radius_shift_explicit_profile_position)
    {
      builder.addBool("radius_shift_explicit_profile_position", true);
    }
  }
  return builder.build();
}

std::string makeRegularStripFaceMetadataJson(double kinetic_time, size_t voronoi_edge_id, size_t even_half_edge_id,
  int strand_even_origin, int strand_odd_origin, SegmentBuilder::BoundaryEventType event_type,
  SegmentBuilder::BoundarySegmentAction segment_action, const char* op)
{
  return SegmentBuilder::MetadataBuilder()
    .addString("event_type", boundaryEventTypeToString(event_type))
    .addString("mesh_type", "regular")
    .addString("segment_action", boundarySegmentActionToString(segment_action))
    .addDouble("time", kinetic_time)
    .addDouble("t", kinetic_time)
    .addSize("voronoi_edge_id", voronoi_edge_id)
    .addSize("even_half_edge_id", even_half_edge_id)
    .addInt("strand_even_origin", strand_even_origin)
    .addInt("strand_odd_origin", strand_odd_origin)
    .addString("op", op)
    .build();
}

std::string makeBoundaryMeshFaceMetadataJson(
  double kinetic_time, const char* mesh_type, size_t half_edge_id, size_t delaunay_face_id, size_t input_branch_id)
{
  SegmentBuilder::MetadataBuilder builder;
  builder.addString("event_type", "boundary_mesh")
    .addString("mesh_type", mesh_type)
    .addDouble("time", kinetic_time)
    .addDouble("t", kinetic_time);
  if (half_edge_id != static_cast<size_t>(-1))
  {
    builder.addSize("half_edge_id", half_edge_id);
  }
  if (delaunay_face_id != static_cast<size_t>(-1))
  {
    builder.addSize("delaunay_face_id", delaunay_face_id);
  }
  if (input_branch_id != static_cast<size_t>(-1))
  {
    builder.addSize("input_branch_id", input_branch_id);
  }
  return builder.build();
}

std::string makeClosingMeshFaceMetadataJson(double kinetic_time, size_t strand_id)
{
  return SegmentBuilder::MetadataBuilder()
    .addString("event_type", "closing_mesh")
    .addString("mesh_type", "closing_cap")
    .addDouble("time", kinetic_time)
    .addDouble("t", kinetic_time)
    .addSize("strand_id", strand_id)
    .build();
}
} // namespace

bool SegmentBuilder::interpolateFlexibleVerticesAlongEdge(
  VoronoiMesh& mesh, std::vector<int>& flex, size_t anchor_old_vertex, size_t anchor_new_vertex)
{
  if (flex.empty())
  {
    return true;
  }
  auto& verts = mesh.getVertices();
  if (anchor_old_vertex >= verts.size() || anchor_new_vertex >= verts.size())
  {
    KINDS_WARNING("interpolateFlexibleVerticesAlongEdge: anchor index out of range (old="
      << anchor_old_vertex << " new=" << anchor_new_vertex << " verts=" << verts.size() << ").");
    return false;
  }
  const glm::dvec3 p0 = verts[anchor_old_vertex];
  const glm::dvec3 p1 = verts[anchor_new_vertex];
  if (!vertexPositionFinite(p0) || !vertexPositionFinite(p1))
  {
    KINDS_WARNING("interpolateFlexibleVerticesAlongEdge: non-finite anchor position (old="
      << anchor_old_vertex << " new=" << anchor_new_vertex << ").");
    return false;
  }
  const double z0 = mesh.vertexKineticTime(anchor_old_vertex);
  const double z1 = mesh.vertexKineticTime(anchor_new_vertex);
  const double denom = z1 - z0;
  const bool is_boundary_interval_meshlet = intersectionMeshletIndexForMesh(mesh).has_value();
  const std::optional<glm::dvec2> raw_uv0
    = is_boundary_interval_meshlet ? boundaryIntervalRawUvAtVertex(mesh, anchor_old_vertex) : std::nullopt;
  const std::optional<glm::dvec2> raw_uv1
    = is_boundary_interval_meshlet ? boundaryIntervalRawUvAtVertex(mesh, anchor_new_vertex) : std::nullopt;
  const std::optional<glm::dvec3> uv0
    = is_boundary_interval_meshlet ? std::nullopt : meshVertexUv(mesh, anchor_old_vertex);
  const std::optional<glm::dvec3> uv1
    = is_boundary_interval_meshlet ? std::nullopt : meshVertexUv(mesh, anchor_new_vertex);
  const bool interpolate_boundary_raw_uv = raw_uv0.has_value() && raw_uv1.has_value();
  const bool interpolate_interior_uv = uv0.has_value() && uv1.has_value();
  const glm::dvec2 profile0 = mesh.triangulationPlaneXY(anchor_old_vertex);
  const glm::dvec2 profile1 = mesh.triangulationPlaneXY(anchor_new_vertex);
  const size_t k = flex.size();
  for (size_t j = 0; j < k; ++j)
  {
    const int fj = flex[j];
    if (fj < 0)
    {
      continue;
    }
    const size_t fju = static_cast<size_t>(fj);
    if (fju >= verts.size())
    {
      continue;
    }
    const double fz = mesh.vertexKineticTime(fju);
    double s = 0.0;
    if (std::abs(denom) > 1e-18)
    {
      s = (fz - z0) / denom;
      if (s < 0.0)
      {
        s = 0.0;
      }
      else if (s > 1.0)
      {
        s = 1.0;
      }
    }
    else
    {
      s = static_cast<double>(j + 1) / static_cast<double>(k + 1);
    }
    glm::dvec3 resolved = p0 + s * (p1 - p0);
    if (!vertexPositionFinite(resolved))
    {
      resolved = (s < 0.5) ? p0 : p1;
      KINDS_WARNING("interpolateFlexibleVerticesAlongEdge: non-finite interpolation result for flex vertex "
        << fj << "; snapping to nearest anchor.");
    }
    mesh.replaceVertex(fju, resolved);
    mesh.setVertexKineticTime(fju, fz);
    mesh.setProfilePlanePosition(fju, profile0 + s * (profile1 - profile0));
    syncResolvedFlexibleVertexMetadata(mesh, fju);
    if (interpolate_boundary_raw_uv)
    {
      const glm::dvec2 raw_interp = *raw_uv0 + s * (*raw_uv1 - *raw_uv0);
      setBoundaryIntervalRawUv(mesh, fju, raw_interp);
      refreshBoundaryIntervalTrianglesIncidentToVertex(mesh, fju);
    }
    else if (interpolate_interior_uv)
    {
      const glm::dvec3 uv_interp = *uv0 + s * (*uv1 - *uv0);
      setMeshVertexUv(mesh, fju, uv_interp);
    }
  }
  return true;
}

void SegmentBuilder::snapFlexibleVerticesToAnchor(VoronoiMesh& mesh, const std::vector<int>& flex, size_t anchor_vertex)
{
  if (flex.empty() || anchor_vertex >= mesh.getVertices().size())
  {
    return;
  }
  const glm::dvec3 anchor_position = mesh.getVertices()[anchor_vertex];
  if (!vertexPositionFinite(anchor_position))
  {
    KINDS_WARNING("snapFlexibleVerticesToAnchor: anchor " << anchor_vertex << " is non-finite.");
    return;
  }
  const glm::dvec2 anchor_profile = mesh.triangulationPlaneXY(anchor_vertex);
  const double anchor_time = mesh.vertexKineticTime(anchor_vertex);
  const bool is_boundary_interval_meshlet = intersectionMeshletIndexForMesh(mesh).has_value();
  const std::optional<glm::dvec2> anchor_raw_uv
    = is_boundary_interval_meshlet ? boundaryIntervalRawUvAtVertex(mesh, anchor_vertex) : std::nullopt;
  const std::optional<glm::dvec3> anchor_uv
    = is_boundary_interval_meshlet ? std::nullopt : meshVertexUv(mesh, anchor_vertex);
  for (int flex_vertex_id : flex)
  {
    if (flex_vertex_id < 0)
    {
      continue;
    }
    const size_t flex_index = static_cast<size_t>(flex_vertex_id);
    if (flex_index >= mesh.getVertices().size())
    {
      continue;
    }
    mesh.replaceVertex(flex_index, anchor_position);
    mesh.setVertexKineticTime(flex_index, anchor_time);
    mesh.setProfilePlanePosition(flex_index, anchor_profile);
    syncResolvedFlexibleVertexMetadata(mesh, flex_index);
    if (anchor_raw_uv.has_value())
    {
      setBoundaryIntervalRawUv(mesh, flex_index, anchor_raw_uv.value());
      refreshBoundaryIntervalTrianglesIncidentToVertex(mesh, flex_index);
    }
    else if (anchor_uv.has_value())
    {
      setMeshVertexUv(mesh, flex_index, anchor_uv.value());
    }
  }
}

static bool raySegmentIntersection(
  const glm::dvec2& C, const glm::dvec2& D, const glm::dvec2& A, const glm::dvec2& B, double& t_out)
{
  glm::dvec2 E = B - A;
  glm::dvec2 AC = A - C;

  double det = glm::cross(D, E);
  if (std::abs(det) < 1e-12)
    return false;

  double t = glm::cross(AC, E) / det;
  double u = glm::cross(AC, D) / det;

  if (t >= 0.0 && u >= 0.0 && u <= 1.0)
  {
    t_out = t;
    return true;
  }
  return false;
}

static std::vector<double> rayCast(
  const std::vector<BoundaryPoint>& polygon, const glm::dvec2& origin, const glm::dvec2& dir)
{
  double lenCP = glm::length(dir);

  if (lenCP < 1e-12)
    return {};

  std::vector<double> hits;
  const size_t n = polygon.size();
  for (size_t i = 0; i < n; ++i)
  {
    const glm::dvec2& A = polygon[i].p;
    const glm::dvec2& B = polygon[(i + 1) % n].p;

    double t;
    if (raySegmentIntersection(origin, dir, A, B, t))
    {
      hits.emplace_back(t);
    }
  }

  return hits;
}

static std::vector<double> rayCast(
  const std::vector<glm::dvec2>& polygon, const glm::dvec2& origin, const glm::dvec2& dir)
{
  double lenCP = glm::length(dir);

  if (lenCP < 1e-12)
    return {};

  std::vector<double> hits;
  const size_t n = polygon.size();
  for (size_t i = 0; i < n; ++i)
  {
    const glm::dvec2& A = polygon[i];
    const glm::dvec2& B = polygon[(i + 1) % n];

    double t;
    if (raySegmentIntersection(origin, dir, A, B, t))
    {
      hits.emplace_back(t);
    }
  }

  return hits;
}

static double relativeDistanceFromCenter(
  const std::vector<BoundaryPoint>& polygon, const glm::dvec2& center, const glm::dvec2& point)
{
  glm::dvec2 dir = point - center;
  auto hits = rayCast(polygon, center, dir);

  if (hits.empty())
    return std::numeric_limits<double>::quiet_NaN();

  double t_max = -std::numeric_limits<double>::infinity();

  for (double t : hits)
  {
    t_max = std::max(t_max, t);
  }

  if (t_max < 0)
  {
    return std::numeric_limits<double>::quiet_NaN();
  }

  // |B - C| = t_max * |D|
  return 1.0 / t_max;
}

static bool isInside(const std::vector<BoundaryPoint>& polygon, const glm::dvec2& center, const glm::dvec2& point)
{
  double rel_dist = relativeDistanceFromCenter(polygon, center, point);
  if (std::isnan(rel_dist))
  {
    throw std::runtime_error("Center lies outside of polygon");
  }
  return rel_dist <= 1.0;
}

[[nodiscard]] glm::dvec3 kinDS::SegmentBuilder::computeVoronoiVertex(size_t half_edge_id, double t) const
{
  return kin_del.computeVoronoiVertexClampedInfinity(half_edge_id, t, false, false);
}

bool clampVoronoiVertices(glm::dvec3& left_vertex, glm::dvec3& right_vertex,
  const std::vector<BoundaryPoint>& boundary_points, const glm::dvec2& centroid)
{
  // don't do this for now
  return true;

  bool left_inside = isInside(boundary_points, centroid, glm::dvec2 { left_vertex[0], left_vertex[1] });
  bool right_inside = isInside(boundary_points, centroid, glm::dvec2 { right_vertex[0], right_vertex[1] });

  if (left_inside && !right_inside)
  {
    KINDS_DEBUG(
      "Clamping right vertex: (" << right_vertex[0] << ", " << right_vertex[1] << ", " << right_vertex[2] << ")");
    // clamp right to the boundary
    glm::dvec2 origin { left_vertex[0], left_vertex[1] };
    glm::dvec2 dir { right_vertex[0] - left_vertex[0], right_vertex[1] - left_vertex[1] };
    auto hits = rayCast(boundary_points, origin, dir);

    double s_min = std::numeric_limits<double>::infinity();

    for (double s : hits)
    {
      if (s >= 0)
      {
        s_min = std::min(s_min, s);
      }
    }

    auto right_vertex_2d = origin + s_min * dir;
    right_vertex = glm::dvec3 { right_vertex_2d[0], right_vertex_2d[1], right_vertex[2] };
    KINDS_DEBUG(
      "Clamped right vertex to: (" << right_vertex[0] << ", " << right_vertex[1] << ", " << right_vertex[2] << ")");
  }
  else if (!left_inside && right_inside)
  {
    KINDS_DEBUG("Clamping left vertex: (" << left_vertex[0] << ", " << left_vertex[1] << ", " << left_vertex[2] << ")");
    // clamp left to the boundary
    glm::dvec2 origin { right_vertex[0], right_vertex[1] };
    glm::dvec2 dir { left_vertex[0] - right_vertex[0], left_vertex[1] - right_vertex[1] };
    auto hits = rayCast(boundary_points, origin, dir);

    double s_min = std::numeric_limits<double>::infinity();

    for (double s : hits)
    {
      if (s >= 0)
      {
        s_min = std::min(s_min, s);
      }
    }

    auto left_vertex_2d = origin + s_min * dir;
    left_vertex = glm::dvec3 { left_vertex_2d[0], left_vertex_2d[1], left_vertex[2] };
    KINDS_DEBUG(
      "Clamped left vertex to: (" << left_vertex[0] << ", " << left_vertex[1] << ", " << left_vertex[2] << ")");
  }
  else if (!left_inside && !right_inside)
  {
    // I'm not sure yet if this will work, especially for connections, but for now we will just discard it
    return false;
  }
  return true;
}

void kinDS::SegmentBuilder::logDiagnosticsMonitoredFaceInsideState(double t, const char* event_context) const
{
  if (!diagnostics || !kin_del.isDiagnosticsMonitoredFaceValid())
  {
    return;
  }
  kin_del.logFaceInsideStateAtTime(KineticDelaunay::kDiagnosticsMonitoredFaceId, t, event_context);
}

namespace
{
std::string formatMeshPairIndex(size_t pair_idx)
{
  if (pair_idx == static_cast<size_t>(-1))
  {
    return "-1";
  }
  return std::to_string(pair_idx);
}
} // namespace

void kinDS::SegmentBuilder::maybeLogDiagnosticsMonitoredDelaunayEdgeTrigger(double t, const char* event_context,
  std::optional<size_t> delaunay_edge_id, std::optional<size_t> mesh_pair_index) const
{
  if (!diagnostics)
  {
    return;
  }
  const bool edge_hit = delaunay_edge_id.has_value()
    && KineticDelaunay::matchesDiagnosticsMonitorId(*delaunay_edge_id, kDiagnosticsMonitoredDelaunayEdgeId);
  const bool pair_hit = mesh_pair_index.has_value()
    && KineticDelaunay::matchesDiagnosticsMonitorId(*mesh_pair_index, kDiagnosticsMonitoredMeshPairId);
  if (!edge_hit && !pair_hit)
  {
    return;
  }
  logDiagnosticsMonitoredDelaunayEdgeState(t, event_context);
}

void kinDS::SegmentBuilder::logDiagnosticsMonitoredDelaunayEdgeState(
  double t, const char* event_context, size_t delaunay_edge_id) const
{
  if (!diagnostics || !KineticDelaunay::isDiagnosticsMonitorIdEnabled(delaunay_edge_id))
  {
    return;
  }

  const auto format_intersection_mesh_pair_metadata = [this](size_t pair_idx) -> std::string
  {
    if (pair_idx == static_cast<size_t>(-1))
    {
      return "(unset)";
    }
    if (pair_idx >= intersection_mesh_pair_metadata.size())
    {
      return "(no metadata slot)";
    }
    const auto& md = intersection_mesh_pair_metadata[pair_idx];
    std::ostringstream oss;
    oss << "{cell=" << md.voronoi_cell_id << " owner_seg=" << md.owner_segment_id
        << " start_d=" << md.start_delaunay_edge_id << " end_d=" << md.end_delaunay_edge_id << "}";
    return oss.str();
  };

  const auto format_intersection_mesh_pair_strip_state = [this](size_t pair_idx) -> std::string
  {
    if (pair_idx == static_cast<size_t>(-1) || pair_idx >= intersection_mesh_pair_last_left_and_right_vertex.size())
    {
      return "(no strip state)";
    }
    const auto& segs = intersection_mesh_pair_last_left_and_right_vertex[pair_idx];
    if (segs.empty())
    {
      return "(empty strip list)";
    }

    // MeshingData::{start,end}_crossing are list iterators into CrossingData::edge_intersections. After erase/insert,
    // std::optional can still be engaged while the iterator is singular ("value-initialized") or orphaned — never
    // dereference the stored iterator for diagnostics. Resolve from live CrossingData via the endpoint half-edge.
    const auto& crossing_data = kin_del.getCrossingData();
    const auto format_crossing_endpoint
      = [&](const char* label, const std::optional<KineticDelaunay::CrossingData::EdgeIntersectionRef>& crossing_opt,
          int half_edge_id)
    {
      if (!crossing_opt.has_value())
      {
        return std::string {};
      }

      std::ostringstream part;
      part << " " << label << "{";
      if (half_edge_id < 0)
      {
        part << "stale_or_unset he=-1}";
        return part.str();
      }

      const size_t d_edge = static_cast<size_t>(half_edge_id) / 2;
      if (d_edge >= crossing_data.delaunay_edge_intersections.size())
      {
        part << "stale he=" << half_edge_id << " (no d-slot)}";
        return part.str();
      }

      const auto& live_refs = crossing_data.delaunay_edge_intersections[d_edge];
      if (live_refs.empty())
      {
        part << "stale he=" << half_edge_id << " d=" << d_edge << " (no live crossing)}";
        return part.str();
      }

      if (live_refs.size() == 1)
      {
        const auto& ref = live_refs.front();
        part << "de=" << ref->delaunay_edge_id << ",ve=" << ref->voronoi_edge_id << "}";
        return part.str();
      }

      part << "he=" << half_edge_id << " d=" << d_edge << " live=[";
      bool first = true;
      for (const auto& ref : live_refs)
      {
        if (!first)
        {
          part << ";";
        }
        first = false;
        part << "de=" << ref->delaunay_edge_id << ",ve=" << ref->voronoi_edge_id;
      }
      part << "]}";
      return part.str();
    };

    std::ostringstream oss;
    size_t seg_i = 0;
    for (const auto& seg : segs)
    {
      if (seg_i > 0)
      {
        oss << " | ";
      }
      oss << "seg" << seg_i << "[v" << seg.mesh_start_vertex_id << "->" << seg.mesh_end_vertex_id << " he("
          << seg.start_half_edge_id << "," << seg.end_half_edge_id << ")";
      oss << format_crossing_endpoint("start_x", seg.start_crossing, seg.start_half_edge_id);
      oss << format_crossing_endpoint("end_x", seg.end_crossing, seg.end_half_edge_id);
      ++seg_i;
    }
    return oss.str();
  };

  const size_t d_edge = delaunay_edge_id;
  const size_t he_even = 2 * d_edge;
  const size_t he_odd = he_even + 1;
  const auto& graph = kin_del.getGraph();

  std::ostringstream header;
  header << "Monitored Delaunay edge d=" << d_edge << " (he " << he_even << "/" << he_odd << ")";
  if (event_context != nullptr && event_context[0] != '\0')
  {
    header << " [" << event_context << "]";
  }
  header << " at t=" << t;
  KINDS_MONITOR(header.str());

  if (he_odd >= graph.halfEdgeSlotCount())
  {
    KINDS_MONITOR("  edge slot out of bounds (he_slots=" << graph.halfEdgeSlotCount() << ")");
    return;
  }

  const bool on_boundary = kin_del.isOnComponentBoundary(he_even);
  const bool even_outside = on_boundary && kin_del.isOnComponentBoundaryOutside(he_even);
  KINDS_MONITOR("  alpha_boundary=" << (on_boundary ? "yes" : "no")
                                    << " even_is_outside=" << (even_outside ? "yes" : "no") << " even_origin="
                                    << graph.halfEdge(he_even).origin << " odd_origin=" << graph.halfEdge(he_odd).origin
                                    << " even_face=" << graph.halfEdge(he_even).face
                                    << " odd_face=" << graph.halfEdge(he_odd).face);

  const auto& crossing_data = kin_del.getCrossingData();
  if (d_edge >= crossing_data.delaunay_edge_intersections.size())
  {
    KINDS_MONITOR("  no delaunay_edge_intersections slot");
    return;
  }

  const auto& refs = crossing_data.delaunay_edge_intersections[d_edge];
  KINDS_MONITOR("  crossing_count=" << refs.size());
  size_t list_idx = 0;
  for (const auto& ref : refs)
  {
    std::ostringstream line;
    line << "  [" << list_idx << "] ve=" << ref->voronoi_edge_id << " param=" << ref->delaunay_edge_param
         << " prev_pair=" << formatMeshPairIndex(ref->prev_segment_mesh_pair_index)
         << " next_pair=" << formatMeshPairIndex(ref->next_segment_mesh_pair_index);
    KINDS_MONITOR(line.str());
    if (ref->prev_segment_mesh_pair_index != static_cast<size_t>(-1))
    {
      KINDS_MONITOR("    prev_pair_meta "
        << format_intersection_mesh_pair_metadata(ref->prev_segment_mesh_pair_index)
        << " strip=" << format_intersection_mesh_pair_strip_state(ref->prev_segment_mesh_pair_index));
    }
    if (ref->next_segment_mesh_pair_index != static_cast<size_t>(-1))
    {
      KINDS_MONITOR("    next_pair_meta "
        << format_intersection_mesh_pair_metadata(ref->next_segment_mesh_pair_index)
        << " strip=" << format_intersection_mesh_pair_strip_state(ref->next_segment_mesh_pair_index));
    }
    if (KineticDelaunay::matchesDiagnosticsMonitorId(ref->prev_segment_mesh_pair_index, kDiagnosticsMonitoredMeshPairId)
      || KineticDelaunay::matchesDiagnosticsMonitorId(
        ref->next_segment_mesh_pair_index, kDiagnosticsMonitoredMeshPairId))
    {
      KINDS_MONITOR("    ** references monitored mesh pair " << kDiagnosticsMonitoredMeshPairId << " **");
    }
    ++list_idx;
  }

  if (KineticDelaunay::isDiagnosticsMonitorIdEnabled(kDiagnosticsMonitoredMeshPairId)
    && kDiagnosticsMonitoredMeshPairId < intersection_mesh_pair_metadata.size())
  {
    KINDS_MONITOR("  monitored_pair_" << kDiagnosticsMonitoredMeshPairId << "_meta="
                                      << format_intersection_mesh_pair_metadata(kDiagnosticsMonitoredMeshPairId));
    KINDS_MONITOR("  monitored_pair_" << kDiagnosticsMonitoredMeshPairId << "_strip="
                                      << format_intersection_mesh_pair_strip_state(kDiagnosticsMonitoredMeshPairId));
    if (kDiagnosticsMonitoredMeshPairId < intersection_meshes.size())
    {
      KINDS_MONITOR(
        "  monitored_pair_" << kDiagnosticsMonitoredMeshPairId
                            << "_verts=" << intersection_meshes[kDiagnosticsMonitoredMeshPairId].getVertexCount()
                            << " tris=" << intersection_meshes[kDiagnosticsMonitoredMeshPairId].getTriangleCount()
                            << " completed="
                            << (kDiagnosticsMonitoredMeshPairId < boundary_meshlet_completed_.size()
                                     && boundary_meshlet_completed_[kDiagnosticsMonitoredMeshPairId]
                                   ? "yes"
                                   : "no"));
    }
  }
}

void kinDS::SegmentBuilder::meshletDiagnosticLogLine(
  const char* tag, size_t half_edge_id, double t, const char* extra_note) const
{
  if (!diagnostics)
  {
    return;
  }

  const size_t even = half_edge_id & ~1;
  const size_t dual_edge = even / 2;
  size_t pi = static_cast<size_t>(-1);
  if (even < half_edge_index_to_segment_mesh_pair_index.size())
  {
    pi = half_edge_index_to_segment_mesh_pair_index[even];
  }
  size_t nv = 0;
  size_t nt = 0;
  size_t nstrip = 0;
  size_t nclosed = 0;
  const bool pi_ok = (pi < meshes.size()) && (pi < segment_mesh_pair_last_left_and_right_vertex.size());
  if (pi_ok)
  {
    nv = meshes[pi].getVertices().size();
    nt = meshes[pi].getTriangleCount();
    nstrip = segment_mesh_pair_last_left_and_right_vertex[pi].size();
    for (const auto& s : segment_mesh_pair_last_left_and_right_vertex[pi])
    {
      if (s.mesh_start_vertex_id >= 0 && s.mesh_end_vertex_id >= 0)
      {
        ++nclosed;
      }
    }
  }
  std::ostringstream oss;
  oss << "meshlet_diag " << tag << " he=" << half_edge_id << " even=" << even << " dual_edge=" << dual_edge
      << " t=" << t;
  if (even < half_edge_index_to_segment_mesh_pair_index.size() && !pi_ok)
  {
    oss << " pair=INVALID(raw_pi=" << pi << ")";
  }
  else
  {
    oss << " pair=" << (pi_ok ? static_cast<long long>(pi) : -1);
  }
  oss << " verts=" << nv << " tris=" << nt << " strips=" << nstrip << " closed_strips=" << nclosed;
  if (extra_note != nullptr && extra_note[0] != '\0')
  {
    oss << ' ' << extra_note;
  }
  KINDS_DEBUG(oss.str());
}

void kinDS::SegmentBuilder::strandInitDiagnosticLogLine(
  const char* tag, size_t strand_id, double t, const char* extra_note) const
{
  if (!diagnostics)
  {
    return;
  }

  std::ostringstream oss;
  oss << "strand_init_diag " << tag << " strand=" << strand_id << " t=" << t;
  if (extra_note != nullptr && extra_note[0] != '\0')
  {
    oss << ' ' << extra_note;
  }
  KINDS_DEBUG(oss.str());
}

void kinDS::SegmentBuilder::logStrandInitDiagnosticsSummary(double t) const
{
  if (!diagnostics)
  {
    return;
  }

  const size_t focus_strand = 0;
  strandInitDiagnosticLogLine("init_summary_begin", focus_strand, t, "post_init_snapshot");

  if (focus_strand >= strand_to_segment_indices.size())
  {
    strandInitDiagnosticLogLine("init_summary_missing_strand", focus_strand, t, "strand_to_segment_indices_too_small");
    return;
  }

  const auto& seg_ids = strand_to_segment_indices[focus_strand];
  std::ostringstream seg_oss;
  seg_oss << "segment_ids=[";
  for (size_t i = 0; i < seg_ids.size(); ++i)
  {
    if (i > 0)
    {
      seg_oss << ',';
    }
    seg_oss << seg_ids[i];
  }
  seg_oss << ']';
  strandInitDiagnosticLogLine("init_summary_segments", focus_strand, t, seg_oss.str().c_str());

  if (!seg_ids.empty())
  {
    const size_t segment_id = seg_ids.front();
    if (segment_id < segment_properties.size())
    {
      std::ostringstream prop_oss;
      prop_oss << "segment_id=" << segment_id << " neighbor_count=" << segment_properties[segment_id].neighbor_count
               << " (0 until finalize/accumulateSegmentProperties)";
      strandInitDiagnosticLogLine("init_summary_segment_props", focus_strand, t, prop_oss.str().c_str());
    }
  }

  if (focus_strand < meshes.size())
  {
    std::ostringstream mesh_oss;
    mesh_oss << "initial_cap_meshlet_index=" << focus_strand << " verts=" << meshes[focus_strand].getVertexCount()
             << " tris=" << meshes[focus_strand].getTriangleCount();
    strandInitDiagnosticLogLine("init_summary_initial_cap_mesh", focus_strand, t, mesh_oss.str().c_str());
    if (meshes[focus_strand].getVertexCount() == 0 && meshes[focus_strand].getTriangleCount() == 0)
    {
      strandInitDiagnosticLogLine("init_summary_empty_initial_cap", focus_strand, t,
        "initial closing cap meshlet is empty — check createClosingMesh trace/triangulation");
    }
  }
  else
  {
    strandInitDiagnosticLogLine("init_summary_no_initial_cap_mesh", focus_strand, t,
      "meshes[] has no entry at strand index (createClosingMesh may not have registered)");
  }

  if (focus_strand < segment_mesh_pairs.size())
  {
    const auto& pair = segment_mesh_pairs[focus_strand];
    std::ostringstream pair_oss;
    pair_oss << "mesh_pair_index=" << focus_strand << " segment_index0=" << pair.segment_index0
             << " segment_index1=" << pair.segment_index1;
    strandInitDiagnosticLogLine("init_summary_cap_pair", focus_strand, t, pair_oss.str().c_str());
  }

  const bool dummy = kin_del.isDummyBoundary(focus_strand);
  const bool live = kin_del.isStrandLiveInGraph(focus_strand);
  std::ostringstream flags_oss;
  flags_oss << "is_dummy_boundary=" << (dummy ? "true" : "false") << " is_live_in_graph=" << (live ? "true" : "false");
  if (focus_strand < kin_del.component_data.component_map.size())
  {
    flags_oss << " component_id=" << kin_del.component_data.component_map[focus_strand];
  }
  strandInitDiagnosticLogLine("init_summary_strand_flags", focus_strand, t, flags_oss.str().c_str());

  size_t incident_strip_pairs = 0;
  auto& graph = kin_del.getGraph();
  for (HalfEdgeDelaunayGraph::IncidentEdgeIterator it = graph.incidentEdgesBegin(focus_strand),
                                                   end = graph.incidentEdgesEnd(focus_strand);
    it != end; ++it)
  {
    const size_t he = *it;
    const size_t even = he & ~1;
    if (even < half_edge_index_to_segment_mesh_pair_index.size())
    {
      const size_t pair_idx = half_edge_index_to_segment_mesh_pair_index[even];
      if (pair_idx < meshes.size())
      {
        ++incident_strip_pairs;
        std::ostringstream strip_oss;
        strip_oss << "incident_he=" << he << " pair_idx=" << pair_idx << " verts=" << meshes[pair_idx].getVertexCount()
                  << " tris=" << meshes[pair_idx].getTriangleCount();
        strandInitDiagnosticLogLine("init_summary_incident_strip", focus_strand, t, strip_oss.str().c_str());
      }
    }
  }
  if (incident_strip_pairs == 0)
  {
    strandInitDiagnosticLogLine(
      "init_summary_no_incident_strips", focus_strand, t, "no Voronoi strip meshlets indexed for incident half-edges");
  }

  const size_t section = static_cast<size_t>(t);
  const auto& branches_at_section = kin_del.getStrandTree().getStrandBranchesByHeight(section);
  for (size_t input_branch_id = 0; input_branch_id < branches_at_section.size(); ++input_branch_id)
  {
    const auto& branch_strands = branches_at_section[input_branch_id];
    const bool in_branch
      = std::find(branch_strands.begin(), branch_strands.end(), focus_strand) != branch_strands.end();
    if (!in_branch)
    {
      continue;
    }
    std::ostringstream branch_oss;
    branch_oss << "input_branch_id=" << input_branch_id << " branch_size=" << branch_strands.size();
    strandInitDiagnosticLogLine("init_summary_input_branch", focus_strand, t, branch_oss.str().c_str());
  }

  strandInitDiagnosticLogLine("init_summary_end", focus_strand, t, nullptr);
}

void kinDS::SegmentBuilder::meshletDiagnosticWarnIfUnexpectedEmptyAfterStartNewMesh(size_t half_edge_even, double t,
  bool initial_left_inside, const VoronoiMesh& mesh, const std::list<MeshingData>& strips) const
{
  if (!diagnostics)
  {
    return;
  }

  const size_t dual_edge = half_edge_even / 2;
  size_t with_any_mesh_id = 0;
  for (const auto& s : strips)
  {
    if (s.mesh_start_vertex_id >= 0 || s.mesh_end_vertex_id >= 0)
    {
      ++with_any_mesh_id;
    }
  }
  const size_t nv = mesh.getVertices().size();
  if (nv == 0 && initial_left_inside)
  {
    KINDS_WARNING("meshlet_diag unexpected_empty after startNewMesh: dual_edge="
      << dual_edge << " t=" << t
      << " — left circumcenter face was inside alpha-shape but mesh has no vertices (strips=" << strips.size() << ").");
  }
  if (nv == 0 && with_any_mesh_id > 0)
  {
    KINDS_WARNING("meshlet_diag unexpected_empty after startNewMesh: dual_edge="
      << dual_edge << " t=" << t
      << " — strip entries reference mesh vertex ids but mesh has no vertices (strips=" << strips.size() << ").");
  }
  // disable for now, probably related to the vertex shift and harmless, but could be a real bug if it happens in other
  // contexts
  /*if (nv > 0 && strips.empty())
  {
    KINDS_WARNING("meshlet_diag inconsistent after startNewMesh: dual_edge="
      << dual_edge << " t=" << t << " - mesh has " << nv << " vertices but strip list is empty.");
  }*/
}

SegmentBuilder::RegularMeshStripIntervalEndpoints SegmentBuilder::regularMeshStripIntervalFromMeshingData(
  const MeshingData& segment, size_t even_half_edge_id, size_t odd_half_edge_id)
{
  RegularMeshStripIntervalEndpoints interval;
  if (segment.start_half_edge_id < 0)
  {
    interval.start_open_voronoi_half_edge_id = even_half_edge_id;
  }
  else
  {
    interval.start_crossing = segment.start_crossing;
    interval.start_crossed_inside_half_edge_id = segment.start_half_edge_id;
  }
  if (segment.end_half_edge_id < 0)
  {
    interval.end_open_voronoi_half_edge_id = odd_half_edge_id;
  }
  else
  {
    interval.end_crossing = segment.end_crossing;
    interval.end_crossed_inside_half_edge_id = segment.end_half_edge_id;
  }
  return interval;
}

std::tuple<size_t, size_t> SegmentBuilder::finishRegularMeshStripInterval(VoronoiMesh& mesh,
  const std::vector<BoundaryPoint>& boundary_polygon, const glm::dvec2& centroid, size_t even_half_edge_id,
  size_t voronoi_edge_id, double t, size_t strand_vertex_id, int strand_even_origin_i, int strand_odd_origin_i,
  BoundaryEventType event_type, BoundarySegmentAction segment_action, const RegularMeshStripIntervalEndpoints& interval,
  size_t last_start_vertex_index, size_t last_end_vertex_index, const std::string& finish_face_metadata,
  const RadiusBoundaryTransitionShiftContext* boundary_transition_shift)
{
  const size_t odd_half_edge_id = even_half_edge_id + 1;
  const std::optional<RadiusBoundaryTransitionCrossingPlacement> start_placement = interval.start_crossing.has_value()
    ? std::make_optional(
        resolveRadiusBoundaryTransitionCrossingPlacement(t, interval.start_crossing.value(), boundary_transition_shift))
    : std::nullopt;
  const std::optional<RadiusBoundaryTransitionCrossingPlacement> end_placement = interval.end_crossing.has_value()
    ? std::make_optional(
        resolveRadiusBoundaryTransitionCrossingPlacement(t, interval.end_crossing.value(), boundary_transition_shift))
    : std::nullopt;

  const auto endpoint_position
    = [&](bool at_start, const std::optional<RadiusBoundaryTransitionCrossingPlacement>& placement)
  {
    const auto crossing = at_start ? interval.start_crossing : interval.end_crossing;
    if (crossing.has_value())
    {
      // In particular for section events, use the crossing stored in MeshingData as the coordinate source. Radius
      // events may replace position_intersection with the one-event destination crossing.
      const auto position_ref = placement.has_value() ? placement->position_intersection : crossing.value();
      const glm::dvec3 p = getCrossingCoordsInMeshSpace(position_ref, t);
      if (placement.has_value() && placement->positionDiffersFromConceptual())
      {
        logRadiusBoundaryTransitionVertexShift("finishRegularMeshStripInterval_endpoint", t,
          placement->conceptual_intersection, placement->position_intersection,
          getCrossingCoordsInMeshSpace(placement->conceptual_intersection, t), p);
      }
      return p;
    }
    const auto open_half_edge
      = at_start ? interval.start_open_voronoi_half_edge_id : interval.end_open_voronoi_half_edge_id;
    return open_half_edge.has_value() ? computeVoronoiVertex(open_half_edge.value(), t) : glm::dvec3(0.0);
  };
  const glm::dvec3 new_start_pos = endpoint_position(true, start_placement);
  const glm::dvec3 new_end_pos = endpoint_position(false, end_placement);

  const auto& graph = kin_del.getGraph();
  const std::optional<size_t> start_vv = interval.start_open_voronoi_half_edge_id.has_value()
    ? std::optional<size_t>(graph.halfEdge(even_half_edge_id).face)
    : std::nullopt;
  const std::optional<size_t> end_vv = interval.end_open_voronoi_half_edge_id.has_value()
    ? std::optional<size_t>(graph.halfEdge(odd_half_edge_id).face)
    : std::nullopt;

  const char* start_source
    = interval.start_crossing.has_value() ? "intersection" : (start_vv.has_value() ? "Voronoi vertex" : "site");
  const char* end_source
    = interval.end_crossing.has_value() ? "intersection" : (end_vv.has_value() ? "Voronoi vertex" : "site");

  auto interval_endpoint_metadata = [&](bool at_start) -> std::string
  {
    const char* pos = at_start ? "left" : "right";
    const char* source = at_start ? start_source : end_source;
    const std::optional<KineticDelaunay::CrossingData::EdgeIntersectionRef> crossing
      = at_start ? interval.start_crossing : interval.end_crossing;
    if (!store_mesh_metadata)
    {
      return {};
    }
    if (source == std::string("intersection") && crossing.has_value())
    {
      const auto& placement = (at_start ? start_placement : end_placement).value();
      MetadataBuilder builder;
      builder.addString("event_type", boundaryEventTypeToString(event_type))
        .addString("segment_action", boundarySegmentActionToString(segment_action));
      return intersectionCrossingVertexMetadata(builder.build(), placement.conceptual_intersection,
        placement.position_intersection, pos, placement.explicit_profile_position.has_value());
    }
    return composeRegularStripVertexMetadata(t, voronoi_edge_id, even_half_edge_id, strand_even_origin_i,
      strand_odd_origin_i, event_type, segment_action, crossing, pos, nullptr, source);
  };

  const std::string meta_start = interval_endpoint_metadata(true);
  const std::string meta_end = interval_endpoint_metadata(false);
  const std::optional<size_t> start_placement_vv
    = start_placement.has_value() && start_placement->snap_voronoi_vertex_id.has_value()
    ? start_placement->snap_voronoi_vertex_id
    : start_vv;
  const std::optional<size_t> end_placement_vv
    = end_placement.has_value() && end_placement->snap_voronoi_vertex_id.has_value()
    ? end_placement->snap_voronoi_vertex_id
    : end_vv;

  const size_t new_start_vertex_index = addMeshletVertex(mesh, boundary_polygon, centroid, new_start_pos,
    strand_vertex_id, t, false, start_placement_vv, meta_start, std::nullopt,
    start_placement.has_value()
      ? MeshletVertexRuntimeInfo { false, start_placement->explicit_profile_position.has_value(),
          start_placement->position_intersection, start_placement->conceptual_intersection,
          start_placement->projection }
      : MeshletVertexRuntimeInfo {});
  const size_t new_end_vertex_index = addMeshletVertex(mesh, boundary_polygon, centroid, new_end_pos, strand_vertex_id,
    t, false, end_placement_vv, meta_end, std::nullopt,
    end_placement.has_value()
      ? MeshletVertexRuntimeInfo { false, end_placement->explicit_profile_position.has_value(),
          end_placement->position_intersection, end_placement->conceptual_intersection, end_placement->projection }
      : MeshletVertexRuntimeInfo {});

  if (last_start_vertex_index == last_end_vertex_index)
  {
    addMeshletTriangle(mesh, new_start_vertex_index, last_end_vertex_index, new_end_vertex_index, finish_face_metadata);
  }
  else if (mesh.vertexKineticTime(last_start_vertex_index) < mesh.vertexKineticTime(last_end_vertex_index))
  {
    addMeshletTriangle(
      mesh, last_start_vertex_index, last_end_vertex_index, new_start_vertex_index, finish_face_metadata);
    addMeshletTriangle(mesh, new_start_vertex_index, last_end_vertex_index, new_end_vertex_index, finish_face_metadata);
  }
  else
  {
    addMeshletTriangle(
      mesh, last_start_vertex_index, last_end_vertex_index, new_end_vertex_index, finish_face_metadata);
    addMeshletTriangle(
      mesh, last_start_vertex_index, new_end_vertex_index, new_start_vertex_index, finish_face_metadata);
  }

  return { new_start_vertex_index, new_end_vertex_index };
}

bool SegmentBuilder::regularMeshStripCrossingTouchesVoronoiEdge(
  const MeshingData& segment, size_t voronoi_edge_id) const
{
  /*if (!kin_del.computeBoundaryOnTheFly())
  {
    const bool start_on = segment.start_crossing.has_value()
      && segment.start_crossing.value()->voronoi_edge_id == voronoi_edge_id;
    const bool end_on = segment.end_crossing.has_value()
      && segment.end_crossing.value()->voronoi_edge_id == voronoi_edge_id;
    return start_on || end_on;
  }*/

  const KineticDelaunay::CrossingData& crossing_data = kin_del.getCrossingData();
  if (voronoi_edge_id >= crossing_data.voronoi_edge_intersections.size())
  {
    return false;
  }

  const auto& v_list = crossing_data.voronoi_edge_intersections[voronoi_edge_id];
  if (v_list.empty())
  {
    return false;
  }

  size_t idx_lo = 0;
  size_t idx_hi = v_list.size() - 1;

  if (segment.start_crossing.has_value())
  {
    const KineticDelaunay::CrossingData::EdgeIntersectionRef start_ref = segment.start_crossing.value();
    if (start_ref->voronoi_edge_id == voronoi_edge_id)
    {
      if (const auto idx = closingMeshIntersectionListIndices(kin_del, start_ref))
      {
        idx_lo = idx->second;
      }
    }
  }

  if (segment.end_crossing.has_value())
  {
    const KineticDelaunay::CrossingData::EdgeIntersectionRef end_ref = segment.end_crossing.value();
    if (end_ref->voronoi_edge_id == voronoi_edge_id)
    {
      if (const auto idx = closingMeshIntersectionListIndices(kin_del, end_ref))
      {
        idx_hi = idx->second;
      }
    }
  }

  if (idx_lo > idx_hi)
  {
    std::swap(idx_lo, idx_hi);
  }

  size_t list_idx = 0;
  for (const KineticDelaunay::CrossingData::EdgeIntersectionRef& ref : v_list)
  {
    if (list_idx >= idx_lo && list_idx <= idx_hi && ref->voronoi_edge_id == voronoi_edge_id)
    {
      return true;
    }
    ++list_idx;
  }

  return false;
}

/**
 * @brief Extends every active strip on one Voronoi edge to time @p t (emit quads as two triangles per strip).
 *
 * @details Expects @ref startNewMesh to have run earlier on this edge so @c last_segments holds seeded corners.
 * Per strip, the “before” state is (@c last_left, @c last_right) = (@c mesh_start_vertex_id, @c mesh_end_vertex_id).
 * @ref finishRegularMeshStripInterval places new corners at the interval endpoints at @p t, connects them to the old
 * corners, then overwrites @c mesh_start_vertex_id / @c mesh_end_vertex_id so a later finish continues from the new
 * front. Strips with either corner unset (@c -1) are skipped (e.g. open strip never fully seeded).
 */
std::vector<SegmentBuilder::MeshingData> kinDS::SegmentBuilder::finishMesh(size_t he_id, double t,
  const std::vector<BoundaryPoint>& boundary_points, BoundaryEventType event_type, BoundarySegmentAction segment_action,
  const RadiusBoundaryTransitionShiftContext* boundary_transition_shift, bool only_adjacent_segment)
{
  std::vector<MeshingData> operated_segments;
  if (kin_del.getGraph().isInfinite(he_id) && kin_del.computeBoundaryOnTheFly())
  {
    meshletDiagnosticLogLine("finish_mesh_skip", he_id, t, "reason=infinite_boundary_on_the_fly");
    return operated_segments;
  }

  size_t segment_mesh_pair_index = half_edge_index_to_segment_mesh_pair_index[he_id];
  if (segment_mesh_pair_index >= meshes.size()
    || segment_mesh_pair_index >= segment_mesh_pair_last_left_and_right_vertex.size())
  {
    if (diagnostics)
    {
      KINDS_WARNING("meshlet_diag finish_mesh: invalid pair index he=" << he_id << " t=" << t
                                                                       << " pair=" << segment_mesh_pair_index);
    }
    return operated_segments;
  }

  const bool retire_meshlet_at_subdivision
    = event_type == BoundaryEventType::Subdivision && segment_action == BoundarySegmentAction::SegmentCompleted;
  const auto mark_retired_if_subdivision = [&]()
  {
    if (retire_meshlet_at_subdivision)
    {
      markInteriorMeshletCompleted(segment_mesh_pair_index);
    }
  };

  meshletDiagnosticLogLine("finish_mesh_enter", he_id, t, "");

  // Live strip state for this dual edge (same list @ref startNewMesh populated).
  auto& last_segments = segment_mesh_pair_last_left_and_right_vertex[segment_mesh_pair_index];
  if (last_segments.empty())
  {
    meshletDiagnosticLogLine("finish_mesh_noop", he_id, t, "last_segments empty (no extension)");
    mark_retired_if_subdivision();
    return operated_segments;
  }

  VoronoiMesh& mesh = meshes[segment_mesh_pair_index];
  const size_t even_id = he_id & ~1;
  const size_t odd_id = even_id + 1;
  const size_t voronoi_edge_id = even_id / 2;
  auto& he = kin_del.getGraph().halfEdge(even_id);
  const auto& twin_he_finish = kin_del.getGraph().halfEdge(odd_id);
  const int strand_even_origin_i = static_cast<int>(he.origin);
  const int strand_odd_origin_i = static_cast<int>(twin_he_finish.origin);
  glm::dvec2 centroid = polygonCentroid(boundary_points);

  const std::string finish_face_meta = composeRegularStripFaceMetadata(
    t, voronoi_edge_id, even_id, strand_even_origin_i, strand_odd_origin_i, event_type, segment_action, "finish_strip");

  if (mesh.getVertices().empty())
  {
    size_t closed_strip_refs = 0;
    for (const auto& s : last_segments)
    {
      if (s.mesh_start_vertex_id >= 0 && s.mesh_end_vertex_id >= 0)
      {
        ++closed_strip_refs;
      }
    }
    if (closed_strip_refs > 0 && diagnostics)
    {
      KINDS_WARNING("meshlet_diag finish_mesh_enter: mesh has no vertices but "
        << closed_strip_refs << " strip(s) have mesh_start/end set (expected non-empty mesh) dual_edge="
        << voronoi_edge_id << " he=" << he_id << " t=" << t << ".");
    }
  }

  // Strand site used as @c origin for addMeshletVertex on this edge (finite endpoint of the even directed edge).
  size_t v = he.origin;
  if (v == -1)
  {
    v = kin_del.getGraph().destination(even_id);
  }

  // Snapshot interval topology from current MeshingData (crossing refs refreshed after any CrossingData rebuild).
  std::vector<RegularMeshStripIntervalEndpoints> finish_intervals;
  finish_intervals.reserve(last_segments.size());
  for (auto& segment : last_segments)
  {
    if (segment.mesh_start_vertex_id < 0 || segment.mesh_end_vertex_id < 0)
    {
      continue;
    }
    refreshMeshingDataCrossingRefs(segment, voronoi_edge_id);
    finish_intervals.push_back(regularMeshStripIntervalFromMeshingData(segment, even_id, odd_id));
  }

  const size_t processable_strips = finish_intervals.size();
  size_t loops_ran = 0;
  const size_t tris_before = diagnostics ? mesh.getTriangleCount() : 0;

  // Parallel walk: one finish_interval per processable strip, in list order.
  auto interval_it = finish_intervals.begin();
  for (auto& segment : last_segments)
  {
    if (segment.mesh_start_vertex_id < 0 || segment.mesh_end_vertex_id < 0)
    {
      continue;
    }

    if (only_adjacent_segment)
    {
      // TODO: describe purpose of this assertion in comment
      assert(boundary_transition_shift != nullptr);
      if (!regularMeshStripCrossingTouchesVoronoiEdge(segment, voronoi_edge_id))
      {
        if (diagnostics)
        {
          const size_t crossing_voronoi_edge_id_start = segment.start_crossing.has_value()
            ? segment.start_crossing.value()->voronoi_edge_id
            : static_cast<size_t>(-1);
          const size_t crossing_voronoi_edge_id_end = segment.end_crossing.has_value()
            ? segment.end_crossing.value()->voronoi_edge_id
            : static_cast<size_t>(-1);
          meshletDiagnosticLogLine("finish_mesh_skip_interval", he_id, t,
            ("reason=only_adjacent_segment crossing_voronoi_edge_ids=" + std::to_string(crossing_voronoi_edge_id_start)
              + "," + std::to_string(crossing_voronoi_edge_id_end))
              .c_str());
        }
        continue;
      }
    }

    ++loops_ran;
    const RegularMeshStripIntervalEndpoints& interval = *interval_it++;
    const size_t last_left = static_cast<size_t>(segment.mesh_start_vertex_id);
    const size_t last_right = static_cast<size_t>(segment.mesh_end_vertex_id);
    const auto [new_start_vertex_index, new_end_vertex_index] = finishRegularMeshStripInterval(mesh, boundary_points,
      centroid, even_id, voronoi_edge_id, t, v, strand_even_origin_i, strand_odd_origin_i, event_type, segment_action,
      interval, last_left, last_right, finish_face_meta, boundary_transition_shift);
    // After finish: strip front is at the new corners; next event’s finish will use these as last_left / last_right.
    segment.mesh_start_vertex_id = static_cast<int>(new_start_vertex_index);
    segment.mesh_end_vertex_id = static_cast<int>(new_end_vertex_index);
    operated_segments.push_back(segment);
  }

  if (diagnostics)
  {
    const size_t tris_after = mesh.getTriangleCount();
    const long long d_tris = static_cast<long long>(tris_after) - static_cast<long long>(tris_before);
    {
      std::ostringstream oss;
      oss << "tris_before=" << tris_before << " tris_after=" << tris_after << " d_tris=" << d_tris
          << " processable_strips=" << processable_strips << " loops_ran=" << loops_ran;
      meshletDiagnosticLogLine("finish_mesh_exit", he_id, t, oss.str().c_str());
    }
    if (processable_strips > 0 && tris_after == tris_before)
    {
      KINDS_WARNING("meshlet_diag finish_mesh: expected triangle extension (processable_strips="
        << processable_strips << ") but triangle count unchanged for dual_edge=" << voronoi_edge_id << " he=" << he_id
        << " t=" << t << ".");
    }
    // Strips may exist in the list but be “open” on one side (only one seeded corner) — nothing to extrude this step.
    if (!last_segments.empty() && processable_strips == 0)
    {
      KINDS_DEBUG("meshlet_diag finish_mesh: strips present but none had both mesh_start/end set — no extension for he="
        << he_id << " dual_edge=" << voronoi_edge_id << " t=" << t << ".");
    }
  }
  mark_retired_if_subdivision();
  return operated_segments;
}

SegmentBuilder::SegmentBuilder(KineticDelaunay& kin_del, std::vector<std::pair<size_t, double>> subdivisions,
  bool create_transformed_mesh, bool visual_debug,
  std::function<void(size_t, std::function<void(size_t)>)> parallel_for)
  : kin_del(kin_del)
  , parallel_for(std::move(parallel_for))
  , section_callback_(std::make_unique<SegmentBuilderSectionCallback>(*this))
  , flip_callback_(std::make_unique<SegmentBuilderFlipCallback>(*this))
  , radius_callback_(std::make_unique<SegmentBuilderRadiusCallback>(*this))
  , crossing_callback_(std::make_unique<SegmentBuilderCrossingCallback>(*this))
  , subdivision_callback_(std::make_unique<SegmentBuilderSubdivisionCallback>(*this))
  , separation_callback_(std::make_unique<SegmentBuilderSeparationCallback>(*this))
  , visual_debug(visual_debug)
  , create_transformed_mesh(create_transformed_mesh)
  , boundary_mesh({ "bark", "interior" })
{
  if (!this->parallel_for)
  {
    this->parallel_for = [](size_t count, const std::function<void(size_t)>& func)
    {
      for (size_t i = 0; i < count; ++i)
      {
        func(i);
      }
    };
  }
  assert(std::is_sorted(
    subdivisions.begin(), subdivisions.end(), [](const auto& a, const auto& b) { return a.second < b.second; }));
  kin_del.setSubdivisionSchedule(std::move(subdivisions));
}

SegmentBuilder::SegmentBuilder(KineticDelaunay& kin_del, bool create_transformed_mesh, bool visual_debug,
  std::function<void(size_t, std::function<void(size_t)>)> parallel_for)
  : kin_del(kin_del)
  , parallel_for(std::move(parallel_for))
  , section_callback_(std::make_unique<SegmentBuilderSectionCallback>(*this))
  , flip_callback_(std::make_unique<SegmentBuilderFlipCallback>(*this))
  , radius_callback_(std::make_unique<SegmentBuilderRadiusCallback>(*this))
  , crossing_callback_(std::make_unique<SegmentBuilderCrossingCallback>(*this))
  , subdivision_callback_(std::make_unique<SegmentBuilderSubdivisionCallback>(*this))
  , separation_callback_(std::make_unique<SegmentBuilderSeparationCallback>(*this))
  , visual_debug(visual_debug)
  , create_transformed_mesh(create_transformed_mesh)
  , boundary_mesh({ "bark", "interior" })
{
  if (!this->parallel_for)
  {
    this->parallel_for = [](size_t count, const std::function<void(size_t)>& func)
    {
      for (size_t i = 0; i < count; ++i)
      {
        func(i);
      }
    };
  }
  kin_del.setSubdivisionSchedule({});
}

SegmentBuilder::~SegmentBuilder() = default;

/// Infinite-line intersection along the first segment direction (matches @c intersect_segments_2d_closing).
static glm::dvec2 intersectAlongFirstSegmentWithLine2D(
  const glm::dvec2& p, const glm::dvec2& p2, const glm::dvec2& q, const glm::dvec2& q2)
{
  const glm::dvec2 r = p2 - p;
  const glm::dvec2 s = q2 - q;
  const double rxs = r.x * s.y - r.y * s.x;
  if (std::abs(rxs) < 1e-30)
  {
    return glm::dvec2(std::numeric_limits<double>::quiet_NaN());
  }
  const double tt = ((q - p).x * s.y - (q - p).y * s.x) / rxs;
  return p + tt * r;
}

static glm::dvec2 intersectSegments(
  const glm::dvec2& p, const glm::dvec2& p2, const glm::dvec2& q, const glm::dvec2& q2)
{
  glm::dvec2 r = p2 - p;
  glm::dvec2 s = q2 - q;

  double rxs = r.x * s.y - r.y * s.x;
  double qpxr = (q - p).x * r.y - (q - p).y * r.x;

  if (rxs == 0.0 && qpxr == 0.0)
    return glm::dvec2(std::numeric_limits<double>::quiet_NaN()); // collinear

  if (rxs == 0.0 && qpxr != 0.0)
    return glm::dvec2(std::numeric_limits<double>::infinity()); // parallel

  double t = ((q - p).x * s.y - (q - p).y * s.x) / rxs;
  double u = ((q - p).x * r.y - (q - p).y * r.x) / rxs;

  if (t < 0 || t > 1 || u < 0 || u > 1)
  {
    KINDS_WARNING("Segments do not intersect within the segment bounds");
  }

  return p + t * r;
}

std::vector<SegmentBuilder::DirectedVoronoiEdgeCrossing> SegmentBuilder::orientCrossingsAlongVoronoiEdge(
  size_t voronoi_edge_id, size_t start_containing_tri_id, bool reverse_traversal,
  size_t even_half_edge_id_for_diagnostics) const
{
  const auto& graph = kin_del.getGraph();
  const KineticDelaunay::CrossingData& crossing_data = kin_del.getCrossingData();
  if (voronoi_edge_id >= crossing_data.voronoi_edge_intersections.size())
  {
    std::ostringstream oss;
    oss << "orientCrossingsAlongVoronoiEdge: voronoi_edge_intersections has no slot for voronoi_edge_id="
        << voronoi_edge_id;
    if (even_half_edge_id_for_diagnostics != static_cast<size_t>(-1))
    {
      oss << " (even_half_edge_id=" << even_half_edge_id_for_diagnostics << ")";
    }
    oss << ", start_containing_tri_id=" << start_containing_tri_id
        << ", voronoi_edge_intersections.size=" << crossing_data.voronoi_edge_intersections.size();
    throw std::runtime_error(oss.str());
  }

  const auto& v_intersections = crossing_data.voronoi_edge_intersections[voronoi_edge_id];
  std::vector<DirectedVoronoiEdgeCrossing> oriented;
  oriented.reserve(v_intersections.size());

  size_t current_face_id = start_containing_tri_id;
  size_t step_index = 0;
  const auto orient_one = [&](const KineticDelaunay::CrossingData::EdgeIntersectionRef& iref)
  {
    const size_t d = iref->delaunay_edge_id;
    const size_t he0 = d * 2;
    const size_t he1 = he0 + 1;
    size_t crossed_he_id = static_cast<size_t>(-1);
    if (graph.halfEdge(he0).face == current_face_id)
    {
      crossed_he_id = he0;
    }
    else if (graph.halfEdge(he1).face == current_face_id)
    {
      crossed_he_id = he1;
    }
    else
    {
      const size_t voronoi_he_even = 2 * voronoi_edge_id;
      const size_t left_voronoi_vertex_id = graph.halfEdge(voronoi_he_even).face;
      const size_t right_voronoi_vertex_id = graph.halfEdge(voronoi_he_even + 1).face;
      KINDS_WARNING("orientCrossingsAlongVoronoiEdge: skipping disconnected crossing at step "
        << step_index << " (voronoi_edge_id=" << voronoi_edge_id << ", delaunay_edge_id=" << d << ", current_face_id="
        << current_face_id << ", face(he0)=" << graph.halfEdge(he0).face << ", face(he1)=" << graph.halfEdge(he1).face
        << ", voronoi_vertices=[" << left_voronoi_vertex_id << "," << right_voronoi_vertex_id << "])");
      return;
    }
    oriented.push_back(DirectedVoronoiEdgeCrossing { crossed_he_id, iref });
    current_face_id = graph.halfEdge(crossed_he_id ^ 1).face;
    ++step_index;
  };

  if (reverse_traversal)
  {
    for (auto it = v_intersections.rbegin(); it != v_intersections.rend(); ++it)
    {
      orient_one(*it);
    }
  }
  else
  {
    for (const KineticDelaunay::CrossingData::EdgeIntersectionRef& iref : v_intersections)
    {
      orient_one(iref);
    }
  }

  return oriented;
}

std::vector<SegmentBuilder::RegularMeshStripIntervalEndpoints>
SegmentBuilder::collectRegularMeshStripIntervalsOnVoronoiEdge(
  size_t even_half_edge_id, size_t voronoi_edge_id, size_t left_containing_tri_id) const
{
  std::vector<RegularMeshStripIntervalEndpoints> intervals;
  const size_t odd_half_edge_id = even_half_edge_id + 1;
  const auto& graph = kin_del.getGraph();
  const size_t left_voronoi_vertex_id = graph.halfEdge(even_half_edge_id).face;
  const size_t right_voronoi_vertex_id = graph.halfEdge(odd_half_edge_id).face;

  bool inside = kin_del.getFaceInside(left_containing_tri_id);
  std::optional<RegularMeshStripIntervalEndpoints> open_interval;

  if (inside || !kin_del.computeBoundaryOnTheFly())
  {
    RegularMeshStripIntervalEndpoints interval;
    interval.start_open_voronoi_half_edge_id = even_half_edge_id;
    open_interval = interval;
  }

  if (kin_del.computeBoundaryOnTheFly())
  {
    const std::vector<DirectedVoronoiEdgeCrossing> directed_crossings
      = orientCrossingsAlongVoronoiEdge(voronoi_edge_id, left_containing_tri_id, false, even_half_edge_id);

    for (const DirectedVoronoiEdgeCrossing& directed : directed_crossings)
    {
      const size_t crossed_he_id = directed.crossed_half_edge_id;
      const auto& crossing_ref = directed.ref;
      const size_t next_face_id = graph.halfEdge(crossed_he_id ^ 1).face;
      const bool next_inside = kin_del.getFaceInside(next_face_id);

      if (inside != next_inside)
      {
        if (inside)
        {
          if (open_interval.has_value())
          {
            open_interval->end_crossing = crossing_ref;
            open_interval->end_crossed_inside_half_edge_id = static_cast<int>(crossed_he_id);
            intervals.push_back(open_interval.value());
            open_interval.reset();
          }
        }
        else
        {
          RegularMeshStripIntervalEndpoints interval;
          interval.start_crossing = crossing_ref;
          interval.start_crossed_inside_half_edge_id = static_cast<int>(crossed_he_id ^ 1);
          open_interval = interval;
        }
      }

      inside = next_inside;
    }
  }

  if (inside || !kin_del.computeBoundaryOnTheFly())
  {
    if (open_interval.has_value())
    {
      open_interval->end_open_voronoi_half_edge_id = odd_half_edge_id;
      intervals.push_back(open_interval.value());
    }
  }

  return intervals;
}

SegmentBuilder::MeshingData SegmentBuilder::meshRegularStripInterval(VoronoiMesh& mesh,
  const std::vector<BoundaryPoint>& boundary_polygon, const glm::dvec2& centroid, size_t even_half_edge_id,
  size_t voronoi_edge_id, double t, int strand_even_origin_i, int strand_odd_origin_i, BoundaryEventType event_type,
  BoundarySegmentAction segment_action, const RegularMeshStripIntervalEndpoints& interval,
  const RadiusBoundaryTransitionShiftContext* boundary_transition_shift)
{
  const auto& graph = kin_del.getGraph();

  auto crossing_endpoint = [&](KineticDelaunay::CrossingData::EdgeIntersectionRef crossing_ref)
    -> std::pair<RadiusBoundaryTransitionCrossingPlacement, glm::dvec3>
  {
    const RadiusBoundaryTransitionCrossingPlacement placement
      = resolveRadiusBoundaryTransitionCrossingPlacement(t, crossing_ref, boundary_transition_shift);
    const glm::dvec3 old_chord_pos = crossingProfilePosition(t, placement.conceptual_intersection);
    const glm::dvec3 intersection_point = crossingProfilePositionFromPlacement(t, placement);
    if (placement.positionDiffersFromConceptual())
    {
      logRadiusBoundaryTransitionVertexShift("meshRegularStripInterval_crossing", t, placement.conceptual_intersection,
        placement.position_intersection, old_chord_pos, intersection_point);
    }
    return { placement, intersection_point };
  };

  const auto strand_id_for_inside_half_edge = [&](int inside_half_edge_id) -> size_t
  {
    if (inside_half_edge_id >= 0)
    {
      const int origin = graph.halfEdge(static_cast<size_t>(inside_half_edge_id)).origin;
      if (origin >= 0)
      {
        return static_cast<size_t>(origin);
      }
    }
    if (strand_even_origin_i >= 0)
    {
      return static_cast<size_t>(strand_even_origin_i);
    }
    if (strand_odd_origin_i >= 0)
    {
      return static_cast<size_t>(strand_odd_origin_i);
    }
    throw std::runtime_error("meshRegularStripInterval: no strand id for crossing vertex transform.");
  };

  int start_half_edge_id = -1;
  size_t start_vertex_index = 0;
  if (interval.start_crossing.has_value())
  {
    start_half_edge_id = interval.start_crossed_inside_half_edge_id;
    const auto [placement, pos] = crossing_endpoint(interval.start_crossing.value());
    const std::string meta_start = store_mesh_metadata
      ? intersectionCrossingVertexMetadata(MetadataBuilder()
                                             .addString("event_type", boundaryEventTypeToString(event_type))
                                             .addString("segment_action", boundarySegmentActionToString(segment_action))
                                             .build(),
          placement.conceptual_intersection, placement.position_intersection, "left",
          placement.explicit_profile_position.has_value())
      : std::string {};
    start_vertex_index
      = addMeshletVertex(mesh, boundary_polygon, centroid, pos, strand_id_for_inside_half_edge(start_half_edge_id), t,
        false, placement.snap_voronoi_vertex_id, meta_start, std::nullopt,
        MeshletVertexRuntimeInfo { false, placement.explicit_profile_position.has_value(),
          placement.position_intersection, placement.conceptual_intersection, placement.projection });
  }
  else if (interval.start_open_voronoi_half_edge_id.has_value())
  {
    const size_t open_he = interval.start_open_voronoi_half_edge_id.value();
    const glm::dvec3 pos = computeVoronoiVertex(open_he, t);
    const size_t voronoi_vertex_id = graph.halfEdge(open_he).face;
    const int origin = graph.halfEdge(open_he).origin;
    const std::string meta_start
      = composeRegularStripVertexMetadata(t, voronoi_edge_id, even_half_edge_id, strand_even_origin_i,
        strand_odd_origin_i, event_type, segment_action, std::nullopt, "left", nullptr, "Voronoi vertex");
    start_vertex_index = addMeshletVertex(
      mesh, boundary_polygon, centroid, pos, origin, t, false, std::optional<size_t>(voronoi_vertex_id), meta_start);
  }

  int end_half_edge_id = -1;
  size_t end_vertex_index = 0;
  if (interval.end_crossing.has_value())
  {
    end_half_edge_id = interval.end_crossed_inside_half_edge_id;
    const auto [placement, pos] = crossing_endpoint(interval.end_crossing.value());
    const std::string meta_end = store_mesh_metadata
      ? intersectionCrossingVertexMetadata(MetadataBuilder()
                                             .addString("event_type", boundaryEventTypeToString(event_type))
                                             .addString("segment_action", boundarySegmentActionToString(segment_action))
                                             .build(),
          placement.conceptual_intersection, placement.position_intersection, "right",
          placement.explicit_profile_position.has_value())
      : std::string {};
    end_vertex_index
      = addMeshletVertex(mesh, boundary_polygon, centroid, pos, strand_id_for_inside_half_edge(end_half_edge_id), t,
        false, placement.snap_voronoi_vertex_id, meta_end, std::nullopt,
        MeshletVertexRuntimeInfo { false, placement.explicit_profile_position.has_value(),
          placement.position_intersection, placement.conceptual_intersection, placement.projection });
  }
  else if (interval.end_open_voronoi_half_edge_id.has_value())
  {
    const size_t open_he = interval.end_open_voronoi_half_edge_id.value();
    const glm::dvec3 pos = computeVoronoiVertex(open_he, t);
    const size_t voronoi_vertex_id = graph.halfEdge(open_he).face;
    const int origin = graph.halfEdge(open_he).origin;
    const std::string meta_end
      = composeRegularStripVertexMetadata(t, voronoi_edge_id, even_half_edge_id, strand_even_origin_i,
        strand_odd_origin_i, event_type, segment_action, std::nullopt, "right", nullptr, "Voronoi vertex");
    end_vertex_index = addMeshletVertex(
      mesh, boundary_polygon, centroid, pos, origin, t, false, std::optional<size_t>(voronoi_vertex_id), meta_end);
  }

  MeshingData segment { static_cast<int>(start_vertex_index), static_cast<int>(end_vertex_index), start_half_edge_id,
    end_half_edge_id };
  segment.start_crossing = interval.start_crossing;
  segment.end_crossing = interval.end_crossing;
  return segment;
}

/**
 * @brief (Re)builds the regular strip meshlet(s) on one undirected Voronoi / Delaunay dual edge at time @p t.
 *
 * @details One call corresponds to “open” or “re-open” meshing on this edge. Each inside-alpha interval along the edge
 * becomes one @ref MeshingData strip stored in @c segment_mesh_pair_last_left_and_right_vertex[pair]. A strip’s state
 * after this function (before any @ref finishMesh):
 *   - @c mesh_start_vertex_id / @c mesh_end_vertex_id: mesh indices of the **left** and **right** strip corners seeded
 *     at interval endpoints (circumcenter or boundary crossing).
 *   - @c start_half_edge_id / @c end_half_edge_id: inside-oriented Delaunay half-edge at each boundary endpoint, or @c
 * -1 for an open circumcenter end.
 *   - @c start_crossing / @c end_crossing: iterators into @c CrossingData when the endpoint is a stored crossing
 * (refreshed after collection). No quads are emitted here; @ref finishMesh advances the corners and adds triangles
 * between old and new positions.
 *
 * @param half_edge_id Directed Delaunay half-edge on the strand side used for lookup (even or odd of the dual edge).
 * @param reuse_existing_pair_and_mesh If true, append vertices into the mesh already tied to this edge and replace the
 *   strip list; if false, allocate a new @c SegmentMeshPair and mesh when none exists yet.
 */
std::vector<SegmentBuilder::MeshingData> SegmentBuilder::startNewMesh(size_t half_edge_id, double t,
  bool reuse_existing_pair_and_mesh, BoundaryEventType event_type, BoundarySegmentAction segment_action,
  const RadiusBoundaryTransitionShiftContext* boundary_transition_shift, bool only_adjacent_segment)
{
  std::vector<MeshingData> operated_segments;
  if (kin_del.getGraph().isInfinite(half_edge_id) && kin_del.computeBoundaryOnTheFly())
  {
    meshletDiagnosticLogLine("start_new_mesh_skip", half_edge_id, t, "reason=infinite_boundary_on_the_fly");
    return operated_segments;
  }

  // Canonical even/odd pair for the undirected dual edge: even = left Voronoi vertex, odd = right Voronoi vertex.
  size_t even_id = half_edge_id & ~1;
  size_t odd_id = even_id + 1;
  if (diagnostics)
  {
    std::ostringstream oss;
    oss << "reuse_existing_pair=" << (reuse_existing_pair_and_mesh ? "true" : "false");
    meshletDiagnosticLogLine("start_new_mesh_begin", half_edge_id, t, oss.str().c_str());
  }

  const auto& graph = kin_del.getGraph();
  const auto& he = graph.halfEdge(even_id);
  const auto& twin_he = graph.halfEdge(odd_id);

  // --- SegmentMeshPair: book-keeping linking this Voronoi-edge mesh to the two incident strand segments (sites). ---
  size_t segment_mesh_pair_index = static_cast<size_t>(-1);
  bool created_new_pair = false;
  if (reuse_existing_pair_and_mesh && even_id < half_edge_index_to_segment_mesh_pair_index.size())
  {
    segment_mesh_pair_index = half_edge_index_to_segment_mesh_pair_index[even_id];
  }

  MeshStructure::SegmentMeshPair segment_mesh_pair;
  segment_mesh_pair.segment_index0 = he.origin == -1 ? -1 : strand_to_segment_indices[he.origin].back();
  segment_mesh_pair.segment_index1 = twin_he.origin == -1 ? -1 : strand_to_segment_indices[twin_he.origin].back();
  if (segment_mesh_pair_index == static_cast<size_t>(-1))
  {
    created_new_pair = true;
    segment_mesh_pair_index = segment_mesh_pairs.size();
    half_edge_index_to_segment_mesh_pair_index[even_id] = segment_mesh_pair_index;
    half_edge_index_to_segment_mesh_pair_index[odd_id] = segment_mesh_pair_index;
    segment_mesh_pairs.push_back(segment_mesh_pair);
  }
  else
  {
    half_edge_index_to_segment_mesh_pair_index[even_id] = segment_mesh_pair_index;
    half_edge_index_to_segment_mesh_pair_index[odd_id] = segment_mesh_pair_index;
    segment_mesh_pairs[segment_mesh_pair_index] = segment_mesh_pair;
  }

  KINDS_DEBUG("Using mesh pair index " << segment_mesh_pair_index << " for half-edge " << even_id
                                       << " (segment indices " << segment_mesh_pair.segment_index0 << ", "
                                       << segment_mesh_pair.segment_index1 << ")");

  // Component boundary + centroid for alpha-shape checks and texture coordinates on new vertices.
  size_t vertex = std::max(he.origin, twin_he.origin);
  size_t component_id = kin_del.component_data.component_map[vertex];

  std::vector<bool> he_visited(graph.halfEdgeSlotCount(), false);
  updateBoundary(t, he_visited, component_id);

  auto& boundary_polygon = kin_del.component_data.component_boundaries[component_id][0];
  auto& centroid = kin_del.component_data.component_centroids[component_id];

  // --- Target mesh: new pair → local mesh registered at end; reuse → append to existing VoronoiMesh for this pair. ---
  VoronoiMesh mesh_local(MeshletExportMaterialNames);
  configureMeshletStorage(mesh_local);
  const bool reuse_in_place = !created_new_pair && segment_mesh_pair_index < meshes.size();
  VoronoiMesh& mesh = reuse_in_place ? meshes[segment_mesh_pair_index] : mesh_local;
  const size_t vertex_count_before = mesh.getVertexCount();

  const size_t voronoi_edge_id = even_id / 2;
  const int strand_even_origin_i = static_cast<int>(he.origin);
  const int strand_odd_origin_i = static_cast<int>(twin_he.origin);

  // Walk the dual edge from the left circumcenter: inside intervals are strips bounded by crossings or open ends.
  const size_t left_voronoi_vertex_id = kin_del.getGraph().halfEdge(even_id).face;
  const size_t left_containing_tri_id = kin_del.getCrossingDataContainingTriId(left_voronoi_vertex_id);
  const size_t right_voronoi_vertex_id = kin_del.getGraph().halfEdge(odd_id).face;
  const bool initial_left_inside = kin_del.getFaceInside(left_containing_tri_id);

  const std::vector<RegularMeshStripIntervalEndpoints> strip_intervals
    = collectRegularMeshStripIntervalsOnVoronoiEdge(even_id, voronoi_edge_id, left_containing_tri_id);
  const size_t voronoi_edge_crossing_count
    = kin_del.computeBoundaryOnTheFly() && voronoi_edge_id < kin_del.getCrossingData().voronoi_edge_intersections.size()
    ? kin_del.getCrossingData().voronoi_edge_intersections[voronoi_edge_id].size()
    : 0;
  const bool ending_inside
    = !strip_intervals.empty() && strip_intervals.back().end_open_voronoi_half_edge_id.has_value();

  // Replace strip list for this pair: previous MeshingData (if any) is discarded; new strips hold seeded corners only.
  if (segment_mesh_pair_index >= segment_mesh_pair_last_left_and_right_vertex.size())
  {
    segment_mesh_pair_last_left_and_right_vertex.resize(segment_mesh_pair_index + 1);
  }
  segment_mesh_pair_last_left_and_right_vertex[segment_mesh_pair_index].clear();
  auto& segments_for_pair = segment_mesh_pair_last_left_and_right_vertex[segment_mesh_pair_index];
  const size_t operated_segment_begin = segments_for_pair.size();

  for (const RegularMeshStripIntervalEndpoints& interval : strip_intervals)
  {
    const bool interval_has_start
      = interval.start_crossing.has_value() || interval.start_open_voronoi_half_edge_id.has_value();
    const bool interval_has_end
      = interval.end_crossing.has_value() || interval.end_open_voronoi_half_edge_id.has_value();
    if (!interval_has_start || !interval_has_end)
    {
      continue;
    }

    if (only_adjacent_segment && reuse_existing_pair_and_mesh)
    {
      assert(boundary_transition_shift != nullptr);
      MeshingData crossing_probe { -1, -1, interval.start_crossed_inside_half_edge_id,
        interval.end_crossed_inside_half_edge_id };
      crossing_probe.start_crossing = interval.start_crossing;
      crossing_probe.end_crossing = interval.end_crossing;
      refreshMeshingDataCrossingRefs(crossing_probe, voronoi_edge_id);
      if (!regularMeshStripCrossingTouchesVoronoiEdge(crossing_probe, voronoi_edge_id))
      {
        if (diagnostics)
        {
          const size_t crossing_voronoi_edge_id_start = crossing_probe.start_crossing.has_value()
            ? crossing_probe.start_crossing.value()->voronoi_edge_id
            : static_cast<size_t>(-1);
          const size_t crossing_voronoi_edge_id_end = crossing_probe.end_crossing.has_value()
            ? crossing_probe.end_crossing.value()->voronoi_edge_id
            : static_cast<size_t>(-1);
          meshletDiagnosticLogLine("start_new_mesh_skip_interval", half_edge_id, t,
            ("reason=only_adjacent_segment crossing_voronoi_edge_ids=" + std::to_string(crossing_voronoi_edge_id_start)
              + "," + std::to_string(crossing_voronoi_edge_id_end))
              .c_str());
        }
        continue;
      }
    }

    segments_for_pair.emplace_back(meshRegularStripInterval(mesh, boundary_polygon, centroid, even_id, voronoi_edge_id,
      t, strand_even_origin_i, strand_odd_origin_i, event_type, segment_action, interval, boundary_transition_shift));
  }

  for (auto& seg : segments_for_pair)
  {
    refreshMeshingDataCrossingRefs(seg, voronoi_edge_id);
  }

  size_t segment_index = 0;
  for (const MeshingData& seg : segments_for_pair)
  {
    if (segment_index++ >= operated_segment_begin)
    {
      operated_segments.push_back(seg);
    }
  }

  meshletDiagnosticWarnIfUnexpectedEmptyAfterStartNewMesh(even_id, t, initial_left_inside, mesh, segments_for_pair);

  const size_t vertex_count_after = mesh.getVertexCount();
  const size_t vertices_added = vertex_count_after - vertex_count_before;

  if (reuse_in_place)
  {
    if (!std::isfinite(mesh.getCreationKineticTime()))
    {
      mesh.setCreationKineticTime(t);
    }
  }
  else
  {
    registerMeshletWithSuffix(std::move(mesh_local), std::string("_voronoi") + std::to_string(voronoi_edge_id), t);
  }

  if (diagnostics)
  {
    std::ostringstream oss;
    oss << "created_new_pair=" << (created_new_pair ? "true" : "false") << " pair_idx=" << segment_mesh_pair_index;
    meshletDiagnosticLogLine("start_new_mesh_end", even_id, t, oss.str().c_str());
  }

  KINDS_DEBUG("startNewMesh summary: vertices_added="
    << vertices_added << " mesh_vertex_count_after=" << vertex_count_after
    << " mesh_vertex_count_before=" << vertex_count_before << " half_edge_id=" << half_edge_id << " even_id=" << even_id
    << " voronoi_edge_id=" << voronoi_edge_id << " t=" << t << " pair_idx=" << segment_mesh_pair_index
    << " created_new_pair=" << (created_new_pair ? "true" : "false")
    << " reuse_in_place=" << (reuse_in_place ? "true" : "false") << " boundary_on_the_fly="
    << (kin_del.computeBoundaryOnTheFly() ? "true" : "false") << " strip_segments=" << segments_for_pair.size()
    << " voronoi_edge_crossing_list_size=" << voronoi_edge_crossing_count
    << " initial_left_inside=" << (initial_left_inside ? "true" : "false")
    << " ending_inside=" << (ending_inside ? "true" : "false") << " component_id=" << component_id
    << " voronoi_vertices=[" << left_voronoi_vertex_id << "," << right_voronoi_vertex_id << "]"
    << " segment_indices=[" << segment_mesh_pair.segment_index0 << "," << segment_mesh_pair.segment_index1 << "]");

  if (segment_mesh_pair_index >= segment_mesh_pair_last_left_and_right_vertex.size())
  {
    KINDS_WARNING("startNewMesh: strip storage out of bounds for pair index " << segment_mesh_pair_index);
  }

  assert(segment_mesh_pairs.size() == meshes.size());
  return operated_segments;
}

std::string SegmentBuilder::composeBoundaryMetadata(
  BoundaryEventType event_type, BoundarySegmentAction segment_action) const
{
  if (!store_mesh_metadata)
  {
    return {};
  }
  return makeBoundaryMeshMetadata(event_type, segment_action);
}

std::string SegmentBuilder::composeRegularStripVertexMetadata(double kinetic_time, size_t voronoi_edge_id,
  size_t even_half_edge_id, int strand_even_origin, int strand_odd_origin, BoundaryEventType event_type,
  BoundarySegmentAction segment_action,
  const std::optional<KineticDelaunay::CrossingData::EdgeIntersectionRef>& crossing, const char* pos, const char* op,
  const char* source, const std::optional<KineticDelaunay::CrossingData::EdgeIntersectionRef>& position_crossing,
  bool radius_shift_explicit_profile_position) const
{
  if (!store_mesh_metadata)
  {
    return {};
  }
  return makeRegularStripVertexMetadataJson(kinetic_time, voronoi_edge_id, even_half_edge_id, strand_even_origin,
    strand_odd_origin, event_type, segment_action, crossing, pos, op, source, position_crossing,
    radius_shift_explicit_profile_position);
}

std::string SegmentBuilder::intersectionCrossingVertexMetadata(const std::string& base_metadata,
  KineticDelaunay::CrossingData::EdgeIntersectionRef conceptual_ref,
  KineticDelaunay::CrossingData::EdgeIntersectionRef position_ref, const char* pos_label,
  bool radius_shift_explicit_profile_position)
{
  MetadataBuilder builder = MetadataBuilder::fromObject(base_metadata).ensureString("source", "intersection");
  builder.addSize("delaunay_edge_id", position_ref->delaunay_edge_id)
    .addSize("voronoi_edge_id", position_ref->voronoi_edge_id)
    .addDouble("stored_d_param", position_ref->delaunay_edge_param);
  if (conceptual_ref != position_ref)
  {
    builder.addSize("conceptual_delaunay_edge_id", conceptual_ref->delaunay_edge_id)
      .addSize("conceptual_voronoi_edge_id", conceptual_ref->voronoi_edge_id);
  }
  if (pos_label != nullptr)
  {
    builder.addString("pos", pos_label);
  }
  if (radius_shift_explicit_profile_position)
  {
    builder.addBool("radius_shift_explicit_profile_position", true);
  }
  return builder.build();
}

void SegmentBuilder::appendIntersectionInterpolationDebugToMetadata(
  MetadataBuilder& builder, const IntersectionInterpolationDebug& debug)
{
  if (!debug.valid())
  {
    return;
  }
  builder.addDouble("interp_d_param", debug.param)
    .addDouble("interp_p0_x", debug.endpoint0.x)
    .addDouble("interp_p0_y", debug.endpoint0.y)
    .addDouble("interp_p0_z", debug.endpoint0.z)
    .addDouble("interp_p1_x", debug.endpoint1.x)
    .addDouble("interp_p1_y", debug.endpoint1.y)
    .addDouble("interp_p1_z", debug.endpoint1.z);
}

std::string SegmentBuilder::composeRegularStripFaceMetadata(double kinetic_time, size_t voronoi_edge_id,
  size_t even_half_edge_id, int strand_even_origin, int strand_odd_origin, BoundaryEventType event_type,
  BoundarySegmentAction segment_action, const char* op) const
{
  if (!store_mesh_metadata)
  {
    return {};
  }
  return makeRegularStripFaceMetadataJson(kinetic_time, voronoi_edge_id, even_half_edge_id, strand_even_origin,
    strand_odd_origin, event_type, segment_action, op);
}

std::string SegmentBuilder::composeBoundaryMeshFaceMetadata(double kinetic_time, const char* mesh_type,
  size_t half_edge_id, size_t delaunay_face_id, size_t input_branch_id) const
{
  if (!store_mesh_metadata)
  {
    return {};
  }
  return makeBoundaryMeshFaceMetadataJson(kinetic_time, mesh_type, half_edge_id, delaunay_face_id, input_branch_id);
}

std::string SegmentBuilder::composeClosingMeshFaceMetadata(double kinetic_time, size_t strand_id) const
{
  if (!store_mesh_metadata)
  {
    return {};
  }
  return makeClosingMeshFaceMetadataJson(kinetic_time, strand_id);
}

void SegmentBuilder::configureMeshletStorage(VoronoiMesh& mesh) const { mesh.setStoreMetadata(store_mesh_metadata); }

std::optional<size_t> SegmentBuilder::regularMeshletIndexForMesh(const VoronoiMesh& mesh) const
{
  for (size_t i = 0; i < meshes.size(); ++i)
  {
    if (&meshes[i] == &mesh)
    {
      return i;
    }
  }
  return std::nullopt;
}

std::optional<size_t> SegmentBuilder::intersectionMeshletIndexForMesh(const VoronoiMesh& mesh) const
{
  for (size_t i = 0; i < intersection_meshes.size(); ++i)
  {
    if (&intersection_meshes[i] == &mesh)
    {
      return i;
    }
  }
  return std::nullopt;
}

void SegmentBuilder::markInteriorMeshletCompleted(size_t meshlet_index)
{
  if (meshlet_index >= meshes.size())
  {
    return;
  }
  if (interior_meshlet_completed_.size() < meshes.size())
  {
    interior_meshlet_completed_.resize(meshes.size(), false);
  }
  interior_meshlet_completed_[meshlet_index] = true;
}

void SegmentBuilder::markInteriorMeshletsCompletedForVoronoiCell(size_t voronoi_cell_id)
{
  auto& graph = kin_del.getGraph();
  if (voronoi_cell_id >= graph.getVertexCount())
  {
    return;
  }
  for (HalfEdgeDelaunayGraph::IncidentEdgeIterator it = graph.incidentEdgesBegin(voronoi_cell_id),
                                                   end = graph.incidentEdgesEnd(voronoi_cell_id);
    it != end; ++it)
  {
    const size_t he_id = *it;
    if (he_id >= half_edge_index_to_segment_mesh_pair_index.size())
    {
      continue;
    }
    const size_t pair_index = half_edge_index_to_segment_mesh_pair_index[he_id];
    if (pair_index != static_cast<size_t>(-1) && pair_index < meshes.size())
    {
      markInteriorMeshletCompleted(pair_index);
    }
  }
}

void SegmentBuilder::markBoundaryMeshletCompleted(size_t meshlet_index)
{
  if (meshlet_index >= intersection_meshes.size())
  {
    return;
  }
  if (boundary_meshlet_completed_.size() < intersection_meshes.size())
  {
    boundary_meshlet_completed_.resize(intersection_meshes.size(), false);
  }
  boundary_meshlet_completed_[meshlet_index] = true;
}

bool SegmentBuilder::isBoundaryMeshletCompleted(size_t meshlet_index) const
{
  return meshlet_index < boundary_meshlet_completed_.size() && boundary_meshlet_completed_[meshlet_index];
}

size_t SegmentBuilder::preferLiveBoundaryMeshPair(size_t pair_idx) const
{
  if (pair_idx == static_cast<size_t>(-1) || pair_idx >= intersection_meshes.size())
  {
    return pair_idx;
  }
  if (!isBoundaryMeshletCompleted(pair_idx))
  {
    return pair_idx;
  }

  size_t cell = static_cast<size_t>(-1);
  size_t start_d = static_cast<size_t>(-1);
  size_t end_d = static_cast<size_t>(-1);
  if (pair_idx < intersection_mesh_pair_metadata.size())
  {
    cell = intersection_mesh_pair_metadata[pair_idx].voronoi_cell_id;
    start_d = intersection_mesh_pair_metadata[pair_idx].start_delaunay_edge_id;
    end_d = intersection_mesh_pair_metadata[pair_idx].end_delaunay_edge_id;
  }

  std::optional<KineticDelaunay::CrossingData::EdgeIntersectionRef> start_crossing;
  std::optional<KineticDelaunay::CrossingData::EdgeIntersectionRef> end_crossing;
  if (pair_idx < intersection_mesh_pair_last_left_and_right_vertex.size()
    && !intersection_mesh_pair_last_left_and_right_vertex[pair_idx].empty())
  {
    const MeshingData& stale = intersection_mesh_pair_last_left_and_right_vertex[pair_idx].front();
    start_crossing = stale.start_crossing;
    end_crossing = stale.end_crossing;
  }

  auto crossings_match = [](const std::optional<KineticDelaunay::CrossingData::EdgeIntersectionRef>& a,
                           const std::optional<KineticDelaunay::CrossingData::EdgeIntersectionRef>& b) -> bool
  {
    if (!a.has_value() && !b.has_value())
    {
      return true;
    }
    return a.has_value() && b.has_value() && a.value() == b.value();
  };
  auto interval_matches = [&](const MeshingData& seg) -> bool
  {
    return (crossings_match(seg.start_crossing, start_crossing) && crossings_match(seg.end_crossing, end_crossing))
      || (crossings_match(seg.start_crossing, end_crossing) && crossings_match(seg.end_crossing, start_crossing));
  };

  for (size_t i = intersection_meshes.size(); i-- > 0;)
  {
    if (i == pair_idx || isBoundaryMeshletCompleted(i))
    {
      continue;
    }
    bool meta_ok = false;
    if (i < intersection_mesh_pair_metadata.size())
    {
      const auto& meta = intersection_mesh_pair_metadata[i];
      meta_ok = meta.voronoi_cell_id == cell && meta.start_delaunay_edge_id == start_d
        && meta.end_delaunay_edge_id == end_d;
    }
    bool interval_ok = false;
    if (i < intersection_mesh_pair_last_left_and_right_vertex.size()
      && !intersection_mesh_pair_last_left_and_right_vertex[i].empty())
    {
      interval_ok = interval_matches(intersection_mesh_pair_last_left_and_right_vertex[i].front());
    }
    if (meta_ok || interval_ok)
    {
      KINDS_DEBUG("preferLiveBoundaryMeshPair: pair " << pair_idx << " retired after subdivision; using live pair " << i);
      return i;
    }
  }
  return pair_idx;
}

bool SegmentBuilder::warnAndSkipIfMeshletCompleted(const VoronoiMesh& mesh, const char* operation, double t) const
{
  if (const std::optional<size_t> regular_idx = regularMeshletIndexForMesh(mesh); regular_idx.has_value())
  {
    if (regular_idx.value() < interior_meshlet_completed_.size() && interior_meshlet_completed_[regular_idx.value()])
    {
      std::ostringstream oss;
      oss << "SegmentBuilder: skipping " << (operation != nullptr ? operation : "meshlet operation")
          << " on completed interior meshlet " << regular_idx.value();
      if (std::isfinite(t))
      {
        oss << " at t=" << t;
      }
      oss << ".";
      KINDS_WARNING(oss.str());
      return true;
    }
    return false;
  }

  if (const std::optional<size_t> boundary_idx = intersectionMeshletIndexForMesh(mesh); boundary_idx.has_value())
  {
    if (boundary_idx.value() < boundary_meshlet_completed_.size() && boundary_meshlet_completed_[boundary_idx.value()])
    {
      std::ostringstream oss;
      oss << "SegmentBuilder: skipping " << (operation != nullptr ? operation : "meshlet operation")
          << " on completed boundary meshlet " << boundary_idx.value();
      if (std::isfinite(t))
      {
        oss << " at t=" << t;
      }
      oss << ".";
      KINDS_WARNING(oss.str());
      return true;
    }
  }
  return false;
}

void SegmentBuilder::logRadiusBoundaryTransitionVertexShift(const char* context, double t,
  KineticDelaunay::CrossingData::EdgeIntersectionRef from_ref,
  KineticDelaunay::CrossingData::EdgeIntersectionRef to_ref, const glm::dvec3& old_pos, const glm::dvec3& new_pos) const
{
  if (from_ref == to_ref)
  {
    return;
  }
  KINDS_DEBUG("Radius boundary transition ["
    << context << "]: vertex shifted from Delaunay edge " << from_ref->delaunay_edge_id << " (Voronoi edge "
    << from_ref->voronoi_edge_id << ") to Delaunay edge " << to_ref->delaunay_edge_id << " (Voronoi edge "
    << to_ref->voronoi_edge_id << ") at t=" << t << "; old=(" << old_pos.x << "," << old_pos.y << "," << old_pos.z
    << ") new=(" << new_pos.x << "," << new_pos.y << "," << new_pos.z << ")");
}

std::optional<KineticDelaunay::CrossingData::EdgeIntersectionRef>
SegmentBuilder::neighborIntersectionOnTargetAlongVoronoiEdge(
  KineticDelaunay::CrossingData::EdgeIntersectionRef crossing_ref, size_t target_delaunay_edge) const
{
  if (!crossing_ref->voronoi_ref.has_value())
  {
    return std::nullopt;
  }

  const size_t voronoi_edge_id = crossing_ref->voronoi_edge_id;
  const auto& crossing_data = kin_del.getCrossingData();
  if (voronoi_edge_id >= crossing_data.voronoi_edge_intersections.size())
  {
    return std::nullopt;
  }

  const auto& v_list = crossing_data.voronoi_edge_intersections[voronoi_edge_id];
  const auto list_it = crossing_ref->voronoi_ref.value();
  if (list_it == v_list.end() || *list_it != crossing_ref)
  {
    KINDS_WARNING(
      "neighborIntersectionOnTargetAlongVoronoiEdge: voronoi_ref stale for voronoi_edge=" << voronoi_edge_id << ".");
    return std::nullopt;
  }

  if (list_it != v_list.begin())
  {
    const auto pit = std::prev(list_it);
    if ((*pit)->delaunay_edge_id == target_delaunay_edge)
    {
      return *pit;
    }
  }

  const auto nit = std::next(list_it);
  if (nit != v_list.end() && (*nit)->delaunay_edge_id == target_delaunay_edge)
  {
    return *nit;
  }

  return std::nullopt;
}

glm::dvec3 SegmentBuilder::crossingProfilePosition(
  double t, KineticDelaunay::CrossingData::EdgeIntersectionRef intersection_ref) const
{
  return getCrossingCoordsInMeshSpace(intersection_ref, t);
}

glm::dvec3 SegmentBuilder::crossingProfilePositionFromPlacement(
  double t, const RadiusBoundaryTransitionCrossingPlacement& placement) const
{
  // A shifted crossing is still exactly another crossing. Use the destination reference as the sole coordinate source
  // so profile/object-space selection always follows getCrossingCoordsInMeshSpace().
  return crossingProfilePosition(t, placement.position_intersection);
}

RadiusBoundaryTransitionCrossingPlacement SegmentBuilder::resolveRadiusBoundaryTransitionCrossingPlacement(double t,
  KineticDelaunay::CrossingData::EdgeIntersectionRef conceptual_ref,
  const RadiusBoundaryTransitionShiftContext* boundary_transition_shift) const
{
  (void)t;
  RadiusBoundaryTransitionCrossingPlacement placement { conceptual_ref, conceptual_ref, std::nullopt };

  if (radius_boundary_transition_shift_enabled && boundary_transition_shift != nullptr
    && boundary_transition_shift->roles_valid)
  {
    const size_t d_edge = conceptual_ref->delaunay_edge_id;
    if (d_edge == boundary_transition_shift->source_delaunay_edges[0]
      || d_edge == boundary_transition_shift->source_delaunay_edges[1])
    {
      if (auto neighbor_opt
        = neighborIntersectionOnTargetAlongVoronoiEdge(conceptual_ref, boundary_transition_shift->target_delaunay_edge);
        neighbor_opt.has_value())
      {
        // Keep the source crossing as the conceptual/topological identity. The adjacent crossing on the target edge is
        // only the coordinate source; addMeshletVertex reconstructs it through getCrossingCoordsInMeshSpace().
        placement.position_intersection = neighbor_opt.value();
      }
    }
  }

  return placement;
}

bool SegmentBuilder::delaunayUndirectedEdgeHasVertex(
  const HalfEdgeDelaunayGraph& graph, size_t delaunay_edge_id, size_t vertex_id)
{
  const size_t he_even = 2 * delaunay_edge_id;
  if (he_even + 1 >= graph.halfEdgeSlotCount())
  {
    return false;
  }
  const int a = graph.halfEdge(he_even).origin;
  const int b = graph.halfEdge(he_even + 1).origin;
  return (a >= 0 && static_cast<size_t>(a) == vertex_id) || (b >= 0 && static_cast<size_t>(b) == vertex_id);
}

std::optional<size_t> SegmentBuilder::oppositeFiniteDelaunayVertexOnUndirectedEdge(
  const HalfEdgeDelaunayGraph& graph, size_t delaunay_edge_id, size_t site_vertex_id)
{
  const size_t he_even = 2 * delaunay_edge_id;
  if (he_even + 1 >= graph.halfEdgeSlotCount())
  {
    return std::nullopt;
  }
  const int a = graph.halfEdge(he_even).origin;
  const int b = graph.halfEdge(he_even + 1).origin;
  if (a >= 0 && static_cast<size_t>(a) == site_vertex_id)
  {
    if (b >= 0)
    {
      return static_cast<size_t>(b);
    }
    return std::nullopt;
  }
  if (b >= 0 && static_cast<size_t>(b) == site_vertex_id)
  {
    if (a >= 0)
    {
      return static_cast<size_t>(a);
    }
    return std::nullopt;
  }
  return std::nullopt;
}

glm::dvec3 SegmentBuilder::radiusTransitionProjectionPosition(
  const RadiusTransitionProjection& projection, double t, bool mesh_space) const
{
  const auto endpoint_position = [&](const RadiusTransitionProjectionEndpoint& endpoint) -> glm::dvec3
  {
    if (endpoint.intersection.has_value())
    {
      const auto ref = endpoint.intersection.value();
      if (mesh_space)
      {
        return getCrossingCoordsInMeshSpace(ref, t);
      }

      kinDS::ensureCrossingIntersectionParamUpToDate(const_cast<KineticDelaunay&>(kin_del), ref, t);
      const auto& graph = kin_del.getGraph();
      const size_t he_even = 2 * ref->delaunay_edge_id;
      if (he_even + 1 < graph.halfEdgeSlotCount())
      {
        const int a = graph.halfEdge(he_even).origin;
        const int b = graph.halfEdge(he_even + 1).origin;
        if (a >= 0 && b >= 0 && std::isfinite(ref->delaunay_edge_param))
        {
          const glm::dvec2 pa = kin_del.getPointAt(t, static_cast<size_t>(a), false, false);
          const glm::dvec2 pb = kin_del.getPointAt(t, static_cast<size_t>(b), false, false);
          const glm::dvec2 p = pa * (1.0 - ref->delaunay_edge_param) + pb * ref->delaunay_edge_param;
          return glm::dvec3(p, t);
        }
      }
    }
    else if (endpoint.site_vertex_id.has_value())
    {
      if (mesh_space)
      {
        return getPointInMeshSpace(endpoint.site_vertex_id.value(), t);
      }
      const glm::dvec2 p = kin_del.getPointAt(t, endpoint.site_vertex_id.value(), false, false);
      return glm::dvec3(p, t);
    }
    else if (endpoint.voronoi_vertex_id.has_value())
    {
      const size_t vv = endpoint.voronoi_vertex_id.value();
      if (mesh_space)
      {
        return computeMeshVoronoiVertexObjectSpace(vv, glm::dvec3(0.0, 0.0, t), 0, t).position;
      }
      const auto& graph = kin_del.getGraph();
      if (vv < graph.faceSlotCount() && graph.isLiveFace(vv))
      {
        return computeVoronoiVertex(graph.face(vv).half_edges[0], t);
      }
    }
    return glm::dvec3(std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN(), t);
  };

  const glm::dvec3 p0 = endpoint_position(projection.endpoints[0]);
  const glm::dvec3 p1 = endpoint_position(projection.endpoints[1]);
  return glm::dvec3(p0 * (1.0 - projection.param) + p1 * projection.param);
}

std::optional<SegmentBuilder::RadiusTransitionSitePlacement> SegmentBuilder::radiusTransitionInterpolatedSitePosition(
  double t, size_t site_vertex_id, size_t strip_delaunay_edge_id,
  const RadiusBoundaryTransitionShiftContext* boundary_transition_shift) const
{
  if (!radius_boundary_transition_shift_enabled || boundary_transition_shift == nullptr
    || !boundary_transition_shift->roles_valid)
  {
    return std::nullopt;
  }
  const auto& graph = kin_del.getGraph();
  const size_t s0 = boundary_transition_shift->source_delaunay_edges[0];
  const size_t s1 = boundary_transition_shift->source_delaunay_edges[1];

  if (boundary_transition_shift->interior_voronoi_vertex_id.has_value())
  {
    // Interior Voronoi vertices reject shifting entirely (see radiusBoundaryTransitionShiftApplicable).
    // Keep this guard so a stale shift context cannot snap corner sites to the VV.
    return std::nullopt;
  }

  if (!delaunayUndirectedEdgeHasVertex(graph, s0, site_vertex_id)
    || !delaunayUndirectedEdgeHasVertex(graph, s1, site_vertex_id))
  {
    return std::nullopt;
  }
  const size_t d_strip = strip_delaunay_edge_id;
  if (d_strip != s0 && d_strip != s1)
  {
    return std::nullopt;
  }
  const size_t d_other = (d_strip == s0) ? s1 : s0;

  struct ProjectionAnchor
  {
    glm::dvec3 old_profile_position {};
    RadiusTransitionProjectionEndpoint shifted_endpoint {};
  };
  const auto anchor = [&](size_t delaunay_edge_id) -> std::optional<ProjectionAnchor>
  {
    const size_t he_even = 2 * delaunay_edge_id;
    if (he_even + 1 >= graph.halfEdgeSlotCount())
    {
      return std::nullopt;
    }
    const std::vector<KineticDelaunay::CrossingData::EdgeIntersectionRef> refs
      = getBoundaryIntersectionsInBoundaryOrder(delaunay_edge_id);
    if (!refs.empty())
    {
      const int o0 = graph.halfEdge(he_even).origin;
      const int o1 = graph.halfEdge(he_even + 1).origin;
      std::optional<KineticDelaunay::CrossingData::EdgeIntersectionRef> ref;
      if (o0 >= 0 && static_cast<size_t>(o0) == site_vertex_id)
      {
        ref = refs.front();
      }
      else if (o1 >= 0 && static_cast<size_t>(o1) == site_vertex_id)
      {
        ref = refs.back();
      }
      if (!ref.has_value())
      {
        return std::nullopt;
      }
      const RadiusBoundaryTransitionCrossingPlacement placement
        = resolveRadiusBoundaryTransitionCrossingPlacement(t, ref.value(), boundary_transition_shift);
      if (placement.explicit_profile_position.has_value() || placement.snap_voronoi_vertex_id.has_value())
      {
        return std::nullopt;
      }
      RadiusTransitionProjectionEndpoint conceptual_endpoint;
      conceptual_endpoint.intersection = ref.value();
      const RadiusTransitionProjection conceptual_projection { { conceptual_endpoint, conceptual_endpoint }, 0.0 };
      ProjectionAnchor out;
      out.old_profile_position = radiusTransitionProjectionPosition(conceptual_projection, t, false);
      out.shifted_endpoint.intersection = placement.position_intersection;
      return out;
    }
    if (!delaunayUndirectedEdgeHasVertex(graph, delaunay_edge_id, site_vertex_id))
    {
      return std::nullopt;
    }
    if (auto v_opp = oppositeFiniteDelaunayVertexOnUndirectedEdge(graph, delaunay_edge_id, site_vertex_id);
      v_opp.has_value())
    {
      const glm::dvec2 xy = kin_del.getPointAt(t, v_opp.value(), false, false);
      ProjectionAnchor out;
      out.old_profile_position = glm::dvec3(xy, t);
      out.shifted_endpoint.site_vertex_id = v_opp.value();
      return out;
    }
    return std::nullopt;
  };

  const auto a0 = anchor(d_strip);
  const auto a1 = anchor(d_other);
  if (!a0.has_value() || !a1.has_value())
  {
    return std::nullopt;
  }

  const glm::dvec3 p0_old = a0->old_profile_position;
  const glm::dvec3 p1_old = a1->old_profile_position;
  const glm::dvec2 site_xy = kin_del.getPointAt(t, site_vertex_id, false, false);
  const glm::dvec2 projection_axis = glm::dvec2(p1_old) - glm::dvec2(p0_old);
  const double projection_axis_len2 = glm::dot(projection_axis, projection_axis);
  if (!(projection_axis_len2 > 1e-30))
  {
    return std::nullopt;
  }

  RadiusTransitionProjection projection;
  projection.endpoints = { a0->shifted_endpoint, a1->shifted_endpoint };
  projection.param = glm::dot(site_xy - glm::dvec2(p0_old), projection_axis) / projection_axis_len2;
  const glm::dvec3 out = radiusTransitionProjectionPosition(projection, t, false);
  KINDS_DEBUG("Radius boundary transition corner site: cell="
    << site_vertex_id << " strip_d=" << d_strip << " other_d=" << d_other << " projection_param=" << projection.param
    << " t=" << t << " out=(" << out.x << "," << out.y << ")");
  return RadiusTransitionSitePlacement { out, std::nullopt, projection };
}

void SegmentBuilder::clearRadiusShiftedSiteCache()
{
  radius_shifted_site_cache_.clear();
}

SegmentBuilder::RadiusShiftedSiteCacheKey SegmentBuilder::makeRadiusShiftedSiteCacheKey(
  size_t site_vertex_id, const RadiusBoundaryTransitionShiftContext& ctx)
{
  size_t lo = ctx.source_delaunay_edges[0];
  size_t hi = ctx.source_delaunay_edges[1];
  if (hi < lo)
  {
    std::swap(lo, hi);
  }
  return RadiusShiftedSiteCacheKey { site_vertex_id, lo, hi };
}

SegmentBuilder::RadiusShiftedSiteCacheEntry* SegmentBuilder::getOrFillRadiusShiftedSiteCache(double t,
  size_t site_vertex_id, const RadiusBoundaryTransitionShiftContext* boundary_transition_shift)
{
  if (boundary_transition_shift == nullptr || !boundary_transition_shift->roles_valid
    || !radius_boundary_transition_shift_enabled)
  {
    return nullptr;
  }

  const RadiusShiftedSiteCacheKey key = makeRadiusShiftedSiteCacheKey(site_vertex_id, *boundary_transition_shift);
  if (auto it = radius_shifted_site_cache_.find(key); it != radius_shifted_site_cache_.end())
  {
    KINDS_DEBUG("RadiusShiftedSiteCache reuse: site=" << site_vertex_id << " sources=(" << key.source_edge_lo << ","
                                                      << key.source_edge_hi << ") has_mesh_pos="
                                                      << (it->second.explicit_mesh_position.has_value() ? "true" : "false")
                                                      << " t=" << t);
    return &it->second;
  }

  const size_t canonical_strip = key.source_edge_lo;
  auto placement
    = radiusTransitionInterpolatedSitePosition(t, site_vertex_id, canonical_strip, boundary_transition_shift);
  if (!placement.has_value())
  {
    return nullptr;
  }

  RadiusShiftedSiteCacheEntry entry;
  entry.placement = std::move(placement.value());
  KINDS_DEBUG("RadiusShiftedSiteCache fill: site=" << site_vertex_id << " sources=(" << key.source_edge_lo << ","
                                                   << key.source_edge_hi << ") canonical_strip=" << canonical_strip
                                                   << " param="
                                                   << (entry.placement.projection.has_value()
                                                         ? entry.placement.projection->param
                                                         : std::numeric_limits<double>::quiet_NaN())
                                                   << " profile=(" << entry.placement.position.x << ","
                                                   << entry.placement.position.y << "," << entry.placement.position.z
                                                   << ") t=" << t);
  auto [ins_it, inserted] = radius_shifted_site_cache_.emplace(key, std::move(entry));
  (void)inserted;
  return &ins_it->second;
}

void SegmentBuilder::noteRadiusShiftedSiteMeshPosition(
  RadiusShiftedSiteCacheEntry& entry, const glm::dvec3& mesh_position, const char* consumer_context)
{
  if (!entry.explicit_mesh_position.has_value())
  {
    entry.explicit_mesh_position = mesh_position;
    KINDS_DEBUG("RadiusShiftedSiteCache mesh_pos store (" << (consumer_context != nullptr ? consumer_context : "?")
                                                          << "): (" << mesh_position.x << "," << mesh_position.y << ","
                                                          << mesh_position.z << ")");
    return;
  }
  const glm::dvec3& cached = entry.explicit_mesh_position.value();
  const bool identical
    = cached.x == mesh_position.x && cached.y == mesh_position.y && cached.z == mesh_position.z;
  KINDS_DEBUG("RadiusShiftedSiteCache mesh_pos reuse (" << (consumer_context != nullptr ? consumer_context : "?")
                                                        << "): identical=" << (identical ? "true" : "false")
                                                        << " cached=(" << cached.x << "," << cached.y << "," << cached.z
                                                        << ") got=(" << mesh_position.x << "," << mesh_position.y << ","
                                                        << mesh_position.z << ")");
}

bool SegmentBuilder::midIntervalMatchesRadiusShiftAnchors(
  KineticDelaunay::CrossingData::EdgeIntersectionRef start, KineticDelaunay::CrossingData::EdgeIntersectionRef end,
  size_t site_vertex_id, const RadiusBoundaryTransitionShiftContext* boundary_transition_shift, double t)
{
  RadiusShiftedSiteCacheEntry* cache_entry
    = getOrFillRadiusShiftedSiteCache(t, site_vertex_id, boundary_transition_shift);
  if (cache_entry == nullptr || !cache_entry->placement.projection.has_value())
  {
    return false;
  }
  const RadiusTransitionProjection& projection = cache_entry->placement.projection.value();
  if (!projection.endpoints[0].intersection.has_value() || !projection.endpoints[1].intersection.has_value())
  {
    return false;
  }
  const auto remapped0 = projection.endpoints[0].intersection.value();
  const auto remapped1 = projection.endpoints[1].intersection.value();
  if (remapped0->delaunay_edge_id != boundary_transition_shift->target_delaunay_edge
    || remapped1->delaunay_edge_id != boundary_transition_shift->target_delaunay_edge)
  {
    return false;
  }
  if (start->delaunay_edge_id != boundary_transition_shift->target_delaunay_edge
    || end->delaunay_edge_id != boundary_transition_shift->target_delaunay_edge)
  {
    return false;
  }
  return (start == remapped0 && end == remapped1) || (start == remapped1 && end == remapped0);
}

void SegmentBuilder::queuePendingRadiusComplementarySplit(size_t intersection_pair_index, size_t edge_a, size_t edge_b,
  const glm::dvec3& mesh_position, double t, size_t site_vertex_id, size_t target_delaunay_edge, bool from_shift,
  std::optional<size_t> snap_voronoi_vertex_id, bool split_last_triangle)
{
  intersection_pair_index = preferLiveBoundaryMeshPair(intersection_pair_index);
  if (intersection_pair_index == static_cast<size_t>(-1))
  {
    return;
  }
  for (const PendingRadiusComplementarySplit& pending : pending_radius_complementary_splits_)
  {
    if (pending.intersection_pair_index == intersection_pair_index)
    {
      KINDS_DEBUG("Radius complementary queue skip: pair=" << intersection_pair_index
                                                             << " already pending site=" << site_vertex_id << " t=" << t);
      return;
    }
  }
  PendingRadiusComplementarySplit pending;
  pending.intersection_pair_index = intersection_pair_index;
  pending.edge_vertex_a = edge_a;
  pending.edge_vertex_b = edge_b;
  pending.mesh_position = mesh_position;
  pending.t = t;
  pending.site_vertex_id = site_vertex_id;
  pending.target_delaunay_edge = target_delaunay_edge;
  pending.from_shift = from_shift;
  pending.split_last_triangle = split_last_triangle;
  pending.snap_voronoi_vertex_id = snap_voronoi_vertex_id;
  pending_radius_complementary_splits_.push_back(pending);
  KINDS_DEBUG("Radius complementary queue: pair=" << intersection_pair_index << " edge=(" << edge_a << "," << edge_b
                                                   << ") site=" << site_vertex_id << " target_d=" << target_delaunay_edge
                                                   << " from_shift=" << (from_shift ? "true" : "false")
                                                   << " split_last=" << (split_last_triangle ? "true" : "false")
                                                   << " mesh_pos=(" << mesh_position.x << "," << mesh_position.y << ","
                                                   << mesh_position.z << ") t=" << t);
}

void SegmentBuilder::maybeQueueRadiusComplementarySplitForStartedMid(size_t intersection_pair_index,
  KineticDelaunay::CrossingData::EdgeIntersectionRef start, KineticDelaunay::CrossingData::EdgeIntersectionRef end,
  double t, size_t site_vertex_id, const RadiusBoundaryTransitionShiftContext* boundary_transition_shift)
{
  intersection_pair_index = preferLiveBoundaryMeshPair(intersection_pair_index);
  if (intersection_pair_index == static_cast<size_t>(-1) || boundary_transition_shift == nullptr)
  {
    return;
  }
  if (!midIntervalMatchesRadiusShiftAnchors(start, end, site_vertex_id, boundary_transition_shift, t))
  {
    return;
  }
  RadiusShiftedSiteCacheEntry* cache_entry
    = getOrFillRadiusShiftedSiteCache(t, site_vertex_id, boundary_transition_shift);
  if (cache_entry == nullptr || !cache_entry->placement.projection.has_value())
  {
    return;
  }
  glm::dvec3 insert_pos;
  if (cache_entry->explicit_mesh_position.has_value())
  {
    insert_pos = cache_entry->explicit_mesh_position.value();
  }
  else
  {
    insert_pos = radiusTransitionProjectionPosition(cache_entry->placement.projection.value(), t, true);
    noteRadiusShiftedSiteMeshPosition(*cache_entry, insert_pos, "complementary_queue_start");
  }
  // Fresh mid-interval from startNewMeshFromIntersections: seed verts are always 0 and 1.
  queuePendingRadiusComplementarySplit(intersection_pair_index, 0, 1, insert_pos, t, site_vertex_id,
    boundary_transition_shift->target_delaunay_edge);
}

void SegmentBuilder::maybeQueueRadiusComplementarySplitForExistingMid(
  double t, size_t site_vertex_id, const RadiusBoundaryTransitionShiftContext* boundary_transition_shift)
{
  if (boundary_transition_shift == nullptr)
  {
    return;
  }
  RadiusShiftedSiteCacheEntry* cache_entry
    = getOrFillRadiusShiftedSiteCache(t, site_vertex_id, boundary_transition_shift);
  if (cache_entry == nullptr || !cache_entry->placement.projection.has_value())
  {
    return;
  }
  const RadiusTransitionProjection& projection = cache_entry->placement.projection.value();
  if (!projection.endpoints[0].intersection.has_value() || !projection.endpoints[1].intersection.has_value())
  {
    return;
  }
  const auto remapped0 = projection.endpoints[0].intersection.value();
  const auto remapped1 = projection.endpoints[1].intersection.value();
  const size_t target_d = boundary_transition_shift->target_delaunay_edge;
  if (remapped0->delaunay_edge_id != target_d || remapped1->delaunay_edge_id != target_d)
  {
    return;
  }

  auto crossings_match = [](const std::optional<KineticDelaunay::CrossingData::EdgeIntersectionRef>& a,
                           KineticDelaunay::CrossingData::EdgeIntersectionRef b) -> bool
  {
    return a.has_value() && a.value() == b;
  };

  size_t pair_idx = static_cast<size_t>(-1);
  if (remapped0->next_segment_mesh_pair_index != static_cast<size_t>(-1)
    && remapped0->next_segment_mesh_pair_index == remapped1->prev_segment_mesh_pair_index)
  {
    pair_idx = remapped0->next_segment_mesh_pair_index;
  }
  else if (remapped1->next_segment_mesh_pair_index != static_cast<size_t>(-1)
    && remapped1->next_segment_mesh_pair_index == remapped0->prev_segment_mesh_pair_index)
  {
    pair_idx = remapped1->next_segment_mesh_pair_index;
  }

  auto segment_matches = [&](const MeshingData& seg) -> bool
  {
    if (!seg.start_crossing.has_value() || !seg.end_crossing.has_value())
    {
      return false;
    }
    if (seg.start_crossing.value()->delaunay_edge_id != target_d
      || seg.end_crossing.value()->delaunay_edge_id != target_d)
    {
      return false;
    }
    return (crossings_match(seg.start_crossing, remapped0) && crossings_match(seg.end_crossing, remapped1))
      || (crossings_match(seg.start_crossing, remapped1) && crossings_match(seg.end_crossing, remapped0));
  };

  if (pair_idx != static_cast<size_t>(-1))
  {
    pair_idx = preferLiveBoundaryMeshPair(pair_idx);
    if (pair_idx >= intersection_mesh_pair_last_left_and_right_vertex.size()
      || intersection_mesh_pair_last_left_and_right_vertex[pair_idx].empty()
      || !segment_matches(intersection_mesh_pair_last_left_and_right_vertex[pair_idx].front()))
    {
      pair_idx = static_cast<size_t>(-1);
    }
  }
  if (pair_idx == static_cast<size_t>(-1))
  {
    for (size_t i = intersection_mesh_pair_last_left_and_right_vertex.size(); i-- > 0;)
    {
      if (isBoundaryMeshletCompleted(i) || intersection_mesh_pair_last_left_and_right_vertex[i].empty())
      {
        continue;
      }
      if (segment_matches(intersection_mesh_pair_last_left_and_right_vertex[i].front()))
      {
        pair_idx = i;
        break;
      }
    }
  }
  if (pair_idx == static_cast<size_t>(-1) || pair_idx >= intersection_meshes.size())
  {
    return;
  }

  // Already queued from a fresh start of this mid-interval.
  for (const PendingRadiusComplementarySplit& pending : pending_radius_complementary_splits_)
  {
    if (pending.intersection_pair_index == pair_idx)
    {
      return;
    }
  }

  const MeshingData& seg = intersection_mesh_pair_last_left_and_right_vertex[pair_idx].front();
  if (seg.mesh_start_vertex_id < 0 || seg.mesh_end_vertex_id < 0)
  {
    return;
  }
  const size_t v0 = static_cast<size_t>(seg.mesh_start_vertex_id);
  const size_t v1 = static_cast<size_t>(seg.mesh_end_vertex_id);
  if (v0 == v1)
  {
    return;
  }

  glm::dvec3 insert_pos;
  if (cache_entry->explicit_mesh_position.has_value())
  {
    insert_pos = cache_entry->explicit_mesh_position.value();
  }
  else
  {
    insert_pos = radiusTransitionProjectionPosition(projection, t, true);
    noteRadiusShiftedSiteMeshPosition(*cache_entry, insert_pos, "complementary_queue_existing");
  }

  // Finished complementary (e.g. 1->2): split the closing edge of the finished strip (last triangle).
  queuePendingRadiusComplementarySplit(pair_idx, v0, v1, insert_pos, t, site_vertex_id, target_d, true, std::nullopt,
    true);
}

std::pair<bool, bool> SegmentBuilder::delaunayEdgeAdjacentInsideFlags(
  size_t delaunay_edge_id, std::optional<size_t> override_face_id, bool override_inside) const
{
  const auto& graph = kin_del.getGraph();
  const size_t he0 = 2 * delaunay_edge_id;
  const size_t he1 = he0 + 1;
  const size_t face_a = graph.halfEdge(he0).face;
  const size_t face_b = graph.halfEdge(he1).face;
  auto inside_of = [&](size_t face_id) -> bool
  {
    if (override_face_id.has_value() && face_id == override_face_id.value())
    {
      return override_inside;
    }
    return kin_del.getFaceInside(face_id);
  };
  return { inside_of(face_a), inside_of(face_b) };
}

SegmentBuilder::RadiusEdgeInsideTransition SegmentBuilder::classifyRadiusEdgeInsideTransition(
  std::pair<bool, bool> pre_inside, std::pair<bool, bool> post_inside)
{
  const bool pre_io = pre_inside.first != pre_inside.second;
  const bool pre_ii = pre_inside.first && pre_inside.second;
  const bool post_io = post_inside.first != post_inside.second;
  const bool post_ii = post_inside.first && post_inside.second;
  if (pre_io && post_ii)
  {
    return RadiusEdgeInsideTransition::IoToIi;
  }
  if (pre_ii && post_io)
  {
    return RadiusEdgeInsideTransition::IiToIo;
  }
  return RadiusEdgeInsideTransition::None;
}

std::optional<glm::dvec3> SegmentBuilder::resolveRadiusNonShiftComplementaryInsertPosition(double t,
  size_t site_vertex_id, const RadiusBoundaryTransitionShiftContext* ctx) const
{
  if (ctx != nullptr && ctx->interior_voronoi_vertex_id.has_value())
  {
    const size_t vv = ctx->interior_voronoi_vertex_id.value();
    return computeMeshVoronoiVertexObjectSpace(vv, glm::dvec3(0.0, 0.0, t), site_vertex_id, t).position;
  }
  return computeMeshSiteVertexPosition(glm::dvec3(0.0, 0.0, t), site_vertex_id, t);
}

void SegmentBuilder::maybeQueueRadiusNonShiftSplitForStartedMid(size_t intersection_pair_index, size_t delaunay_edge_id,
  size_t mid_cell, double t, size_t affected_face_id, bool pre_affected_inside, bool post_affected_inside,
  const RadiusBoundaryTransitionShiftContext* roles_ctx)
{
  intersection_pair_index = preferLiveBoundaryMeshPair(intersection_pair_index);
  if (intersection_pair_index == static_cast<size_t>(-1) || mid_cell == static_cast<size_t>(-1))
  {
    return;
  }
  const auto pre = delaunayEdgeAdjacentInsideFlags(delaunay_edge_id, affected_face_id, pre_affected_inside);
  const auto post = delaunayEdgeAdjacentInsideFlags(delaunay_edge_id, affected_face_id, post_affected_inside);
  if (classifyRadiusEdgeInsideTransition(pre, post) != RadiusEdgeInsideTransition::IoToIi)
  {
    return;
  }

  const size_t site_id = mid_cell;
  std::optional<size_t> snap_vv;
  if (roles_ctx != nullptr && roles_ctx->interior_voronoi_vertex_id.has_value())
  {
    snap_vv = roles_ctx->interior_voronoi_vertex_id;
  }
  const auto insert_pos = resolveRadiusNonShiftComplementaryInsertPosition(t, site_id, roles_ctx);
  if (!insert_pos.has_value())
  {
    return;
  }
  KINDS_DEBUG("Radius non-shift IO→II: queue started mid pair=" << intersection_pair_index << " d_edge="
                                                                << delaunay_edge_id << " site=" << site_id << " t=" << t);
  queuePendingRadiusComplementarySplit(intersection_pair_index, 0, 1, insert_pos.value(), t, site_id, delaunay_edge_id,
    false, snap_vv);
}

void SegmentBuilder::maybeQueueRadiusNonShiftSplitForFinishedMid(size_t intersection_pair_index, size_t delaunay_edge_id,
  size_t mid_cell, double t, size_t affected_face_id, bool pre_affected_inside, bool post_affected_inside,
  const RadiusBoundaryTransitionShiftContext* roles_ctx)
{
  intersection_pair_index = preferLiveBoundaryMeshPair(intersection_pair_index);
  if (intersection_pair_index == static_cast<size_t>(-1) || mid_cell == static_cast<size_t>(-1))
  {
    return;
  }
  if (intersection_pair_index >= intersection_mesh_pair_last_left_and_right_vertex.size()
    || intersection_mesh_pair_last_left_and_right_vertex[intersection_pair_index].empty())
  {
    return;
  }
  const auto pre = delaunayEdgeAdjacentInsideFlags(delaunay_edge_id, affected_face_id, pre_affected_inside);
  const auto post = delaunayEdgeAdjacentInsideFlags(delaunay_edge_id, affected_face_id, post_affected_inside);
  if (classifyRadiusEdgeInsideTransition(pre, post) != RadiusEdgeInsideTransition::IiToIo)
  {
    return;
  }

  const MeshingData& seg = intersection_mesh_pair_last_left_and_right_vertex[intersection_pair_index].front();
  if (seg.mesh_start_vertex_id < 0 || seg.mesh_end_vertex_id < 0)
  {
    return;
  }
  const size_t v0 = static_cast<size_t>(seg.mesh_start_vertex_id);
  const size_t v1 = static_cast<size_t>(seg.mesh_end_vertex_id);
  if (v0 == v1)
  {
    return;
  }

  const size_t site_id = mid_cell;
  std::optional<size_t> snap_vv;
  if (roles_ctx != nullptr && roles_ctx->interior_voronoi_vertex_id.has_value())
  {
    snap_vv = roles_ctx->interior_voronoi_vertex_id;
  }
  const auto insert_pos = resolveRadiusNonShiftComplementaryInsertPosition(t, site_id, roles_ctx);
  if (!insert_pos.has_value())
  {
    return;
  }
  KINDS_DEBUG("Radius non-shift II→IO: queue finished mid pair=" << intersection_pair_index << " d_edge="
                                                                  << delaunay_edge_id << " site=" << site_id
                                                                  << " edge=(" << v0 << "," << v1 << ") t=" << t);
  queuePendingRadiusComplementarySplit(intersection_pair_index, v0, v1, insert_pos.value(), t, site_id, delaunay_edge_id,
    false, snap_vv, true);
}

void SegmentBuilder::maybeQueueRadiusNonShiftSplitForExistingMid(
  double t, size_t site_vertex_id, const RadiusBoundaryTransitionShiftContext* roles_ctx)
{
  if (roles_ctx == nullptr || !roles_ctx->roles_valid)
  {
    return;
  }
  // Reuse remapped-anchor mid lookup from the shift path (neighbors still exist without shifting XY).
  const size_t s0 = roles_ctx->source_delaunay_edges[0];
  const size_t s1 = roles_ctx->source_delaunay_edges[1];
  const size_t target_d = roles_ctx->target_delaunay_edge;

  auto site_side_crossing = [&](size_t delaunay_edge_id) -> std::optional<KineticDelaunay::CrossingData::EdgeIntersectionRef>
  {
    const auto& graph = kin_del.getGraph();
    const size_t he_even = 2 * delaunay_edge_id;
    const std::vector<KineticDelaunay::CrossingData::EdgeIntersectionRef> refs
      = getBoundaryIntersectionsInBoundaryOrder(delaunay_edge_id);
    if (refs.empty())
    {
      return std::nullopt;
    }
    const int o0 = graph.halfEdge(he_even).origin;
    const int o1 = graph.halfEdge(he_even + 1).origin;
    if (o0 >= 0 && static_cast<size_t>(o0) == site_vertex_id)
    {
      return refs.front();
    }
    if (o1 >= 0 && static_cast<size_t>(o1) == site_vertex_id)
    {
      return refs.back();
    }
    return std::nullopt;
  };

  const auto c0 = site_side_crossing(s0);
  const auto c1 = site_side_crossing(s1);
  if (!c0.has_value() || !c1.has_value())
  {
    return;
  }
  const auto n0 = neighborIntersectionOnTargetAlongVoronoiEdge(c0.value(), target_d);
  const auto n1 = neighborIntersectionOnTargetAlongVoronoiEdge(c1.value(), target_d);
  if (!n0.has_value() || !n1.has_value())
  {
    return;
  }

  size_t pair_idx = static_cast<size_t>(-1);
  if (n0.value()->next_segment_mesh_pair_index != static_cast<size_t>(-1)
    && n0.value()->next_segment_mesh_pair_index == n1.value()->prev_segment_mesh_pair_index)
  {
    pair_idx = n0.value()->next_segment_mesh_pair_index;
  }
  else if (n1.value()->next_segment_mesh_pair_index != static_cast<size_t>(-1)
    && n1.value()->next_segment_mesh_pair_index == n0.value()->prev_segment_mesh_pair_index)
  {
    pair_idx = n1.value()->next_segment_mesh_pair_index;
  }
  if (pair_idx != static_cast<size_t>(-1))
  {
    pair_idx = preferLiveBoundaryMeshPair(pair_idx);
  }
  if (pair_idx == static_cast<size_t>(-1) || pair_idx >= intersection_mesh_pair_last_left_and_right_vertex.size()
    || intersection_mesh_pair_last_left_and_right_vertex[pair_idx].empty())
  {
    return;
  }
  const MeshingData& seg = intersection_mesh_pair_last_left_and_right_vertex[pair_idx].front();
  if (seg.mesh_start_vertex_id < 0 || seg.mesh_end_vertex_id < 0)
  {
    return;
  }
  const size_t v0 = static_cast<size_t>(seg.mesh_start_vertex_id);
  const size_t v1 = static_cast<size_t>(seg.mesh_end_vertex_id);
  if (v0 == v1)
  {
    return;
  }
  const auto insert_pos = resolveRadiusNonShiftComplementaryInsertPosition(t, site_vertex_id, roles_ctx);
  if (!insert_pos.has_value())
  {
    return;
  }
  KINDS_DEBUG("Radius non-shift existing mid queue: pair=" << pair_idx << " site=" << site_vertex_id << " target_d="
                                                             << target_d << " t=" << t);
  queuePendingRadiusComplementarySplit(pair_idx, v0, v1, insert_pos.value(), t, site_vertex_id, target_d, false,
    roles_ctx->interior_voronoi_vertex_id, true);
}

bool SegmentBuilder::applyPendingRadiusComplementarySplit(
  const PendingRadiusComplementarySplit& pending, bool allow_completed)
{
  if (pending.intersection_pair_index >= intersection_meshes.size())
  {
    KINDS_WARNING("flushPendingRadiusComplementarySplits: invalid pair=" << pending.intersection_pair_index);
    return false;
  }
  VoronoiMesh& mesh = intersection_meshes[pending.intersection_pair_index];
  if (!allow_completed
    && warnAndSkipIfMeshletCompleted(mesh, "flushPendingRadiusComplementarySplits/splitTriangle", pending.t))
  {
    return false;
  }

  size_t v0 = pending.edge_vertex_a;
  size_t v1 = pending.edge_vertex_b;
  bool prefer_last_triangle = pending.split_last_triangle;
  if (prefer_last_triangle && pending.intersection_pair_index < intersection_mesh_pair_last_left_and_right_vertex.size()
    && !intersection_mesh_pair_last_left_and_right_vertex[pending.intersection_pair_index].empty())
  {
    const MeshingData& seg
      = intersection_mesh_pair_last_left_and_right_vertex[pending.intersection_pair_index].front();
    if (seg.mesh_start_vertex_id >= 0 && seg.mesh_end_vertex_id >= 0
      && seg.mesh_start_vertex_id != seg.mesh_end_vertex_id)
    {
      v0 = static_cast<size_t>(seg.mesh_start_vertex_id);
      v1 = static_cast<size_t>(seg.mesh_end_vertex_id);
    }
  }

  if (v0 >= mesh.getVertices().size() || v1 >= mesh.getVertices().size() || v0 == v1)
  {
    KINDS_WARNING("flushPendingRadiusComplementarySplits: bad edge verts pair=" << pending.intersection_pair_index
                                                                                << " (" << v0 << "," << v1 << ")");
    return false;
  }

  const std::string split_meta = [&]() -> std::string
  {
    if (!store_mesh_metadata)
    {
      return {};
    }
    MetadataBuilder builder;
    builder.addString("event_type", boundaryEventTypeToString(BoundaryEventType::Radius))
      .addString("segment_action", boundarySegmentActionToString(BoundarySegmentAction::SegmentRemapped))
      .addString("source", pending.snap_voronoi_vertex_id.has_value() ? "Voronoi vertex" : "site")
      .addString("pos", pending.from_shift ? "complementary_shift" : "complementary_noshift")
      .addBool("split_triangle", true)
      .addBool("radius_shifted_corner", pending.from_shift)
      .addBool("radius_noshift_corner", !pending.from_shift)
      .addSize("strand_id", pending.site_vertex_id)
      .addSize("target_delaunay_edge", pending.target_delaunay_edge)
      .addSize("intersection_pair_index", pending.intersection_pair_index)
      .addDouble("time", pending.t)
      .addDouble("t", pending.t);
    if (pending.snap_voronoi_vertex_id.has_value())
    {
      builder.addSize("voronoi_vertex_id", pending.snap_voronoi_vertex_id.value());
    }
    return builder.build();
  }();

  auto find_tri = [&](size_t a, size_t b) -> std::optional<size_t>
  {
    const auto& tris = mesh.getTriangles();
    const size_t n_tri = tris.size() / 3;
    if (prefer_last_triangle)
    {
      for (size_t tri = n_tri; tri-- > 0;)
      {
        for (size_t e = 0; e < 3; ++e)
        {
          if (tris[3 * tri + e] == a && tris[3 * tri + ((e + 1) % 3)] == b)
          {
            return tri;
          }
        }
      }
      return std::nullopt;
    }
    // Prefer early strip faces (seed edge 0-1 lives on the first triangles).
    const size_t scan_n = n_tri < 2 ? n_tri : 2;
    for (size_t pass = 0; pass < 2; ++pass)
    {
      const size_t begin = (pass == 0) ? 0 : scan_n;
      const size_t end = (pass == 0) ? scan_n : n_tri;
      for (size_t tri = begin; tri < end; ++tri)
      {
        for (size_t e = 0; e < 3; ++e)
        {
          if (tris[3 * tri + e] == a && tris[3 * tri + ((e + 1) % 3)] == b)
          {
            return tri;
          }
        }
      }
    }
    return std::nullopt;
  };

  auto try_edge = [&](size_t a, size_t b) -> std::optional<std::pair<size_t, size_t>>
  {
    std::optional<size_t> split_tri = find_tri(a, b);
    size_t edge_a = a;
    size_t edge_b = b;
    if (!split_tri.has_value())
    {
      split_tri = find_tri(b, a);
      edge_a = b;
      edge_b = a;
    }
    if (!split_tri.has_value())
    {
      return std::nullopt;
    }
    return std::make_pair(edge_a, edge_b);
  };

  std::optional<std::pair<size_t, size_t>> edge = try_edge(v0, v1);
  if (!edge.has_value() && prefer_last_triangle)
  {
    prefer_last_triangle = false;
    v0 = pending.edge_vertex_a;
    v1 = pending.edge_vertex_b;
    if (v0 < mesh.getVertices().size() && v1 < mesh.getVertices().size() && v0 != v1)
    {
      edge = try_edge(v0, v1);
    }
  }
  if (!edge.has_value())
  {
    KINDS_WARNING("flushPendingRadiusComplementarySplits: no triangle on edge pair="
      << pending.intersection_pair_index << " verts=(" << v0 << "," << v1 << ") tri_count="
      << mesh.getTriangleCount() << " t=" << pending.t);
    return false;
  }

  std::optional<size_t> split_tri = find_tri(edge->first, edge->second);
  if (!split_tri.has_value())
  {
    return false;
  }
  const size_t tri_corner0 = mesh.triangleCornerIndex(split_tri.value(), edge->first);
  const size_t tri_corner1 = mesh.triangleCornerIndex(split_tri.value(), edge->second);
  if (tri_corner0 == static_cast<size_t>(-1) || tri_corner1 == static_cast<size_t>(-1))
  {
    KINDS_WARNING("flushPendingRadiusComplementarySplits: triangleCornerIndex failed pair="
      << pending.intersection_pair_index << " tri=" << split_tri.value());
    return false;
  }

  const auto [new_vid, new_tri]
    = mesh.splitTriangle(tri_corner0, tri_corner1, pending.mesh_position, split_meta, pending.t);
  if (new_vid == static_cast<size_t>(-1) || new_tri == static_cast<size_t>(-1))
  {
    KINDS_WARNING("flushPendingRadiusComplementarySplits: splitTriangle failed pair="
      << pending.intersection_pair_index);
    return false;
  }
  KINDS_DEBUG("flushPendingRadiusComplementarySplits: pair=" << pending.intersection_pair_index << " tri="
                                                               << split_tri.value() << " -> new_vid=" << new_vid
                                                               << " new_tri=" << new_tri << " t=" << pending.t);
  return true;
}

void SegmentBuilder::flushPendingRadiusComplementarySplitsForPair(size_t intersection_pair_index)
{
  if (intersection_pair_index == static_cast<size_t>(-1) || pending_radius_complementary_splits_.empty())
  {
    return;
  }
  std::vector<PendingRadiusComplementarySplit> remaining;
  remaining.reserve(pending_radius_complementary_splits_.size());
  for (const PendingRadiusComplementarySplit& pending : pending_radius_complementary_splits_)
  {
    if (pending.intersection_pair_index != intersection_pair_index)
    {
      remaining.push_back(pending);
      continue;
    }
    applyPendingRadiusComplementarySplit(pending, true);
  }
  pending_radius_complementary_splits_.swap(remaining);
}

void SegmentBuilder::flushPendingRadiusComplementarySplits()
{
  if (pending_radius_complementary_splits_.empty())
  {
    return;
  }
  KINDS_INFO("flushPendingRadiusComplementarySplits: applying " << pending_radius_complementary_splits_.size()
                                                                << " deferred complementary split(s)");

  for (const PendingRadiusComplementarySplit& pending : pending_radius_complementary_splits_)
  {
    applyPendingRadiusComplementarySplit(pending, true);
  }

  pending_radius_complementary_splits_.clear();
}

size_t SegmentBuilder::resolveIntersectionMeshPairIndex(size_t voronoi_cell_id,
  std::optional<KineticDelaunay::CrossingData::EdgeIntersectionRef> start_intersection,
  std::optional<KineticDelaunay::CrossingData::EdgeIntersectionRef> end_intersection, double event_time) const
{
  const auto live_pair = [this](size_t pair_idx) -> size_t { return preferLiveBoundaryMeshPair(pair_idx); };

  if (!start_intersection.has_value() && !end_intersection.has_value())
  {
    const std::string msg = "resolveIntersectionMeshPairIndex: both start_intersection and end_intersection are null.";
    KINDS_ERROR(msg);
    throw std::runtime_error(msg);
  }

  if (start_intersection.has_value() && end_intersection.has_value())
  {
    if (start_intersection.value()->delaunay_edge_id != end_intersection.value()->delaunay_edge_id)
    {
      std::ostringstream oss;
      oss << "resolveIntersectionMeshPairIndex: start/end intersection Delaunay edge mismatch (start="
          << start_intersection.value()->delaunay_edge_id << ", end=" << end_intersection.value()->delaunay_edge_id
          << ", voronoi_cell_id=" << voronoi_cell_id << ", event_time=" << event_time << ").";
      KINDS_ERROR(oss.str());
      throw std::runtime_error(oss.str());
    }

    // We assume correct order here, check that
    auto& start_value = start_intersection.value();
    auto& end_value = end_intersection.value();

    const size_t d_edge_id = start_value->delaunay_edge_id;
    const auto& d_list = kin_del.getCrossingData().delaunay_edge_intersections[d_edge_id];
    const auto it_start = std::find(d_list.begin(), d_list.end(), start_value);
    const auto it_end = std::find(d_list.begin(), d_list.end(), end_value);
    if (it_start != d_list.end() && it_end != d_list.end() && it_start != it_end)
    {
      KineticDelaunay::CrossingData::EdgeIntersectionRef list_first;
      KineticDelaunay::CrossingData::EdgeIntersectionRef list_second;
      if (listConstIteratorPrecedes(d_list.begin(), d_list.end(), it_start, it_end))
      {
        list_first = start_value;
        list_second = end_value;
      }
      else
      {
        list_first = end_value;
        list_second = start_value;
      }
      if (list_first->next_segment_mesh_pair_index == list_second->prev_segment_mesh_pair_index)
      {
        return live_pair(list_first->next_segment_mesh_pair_index);
      }
    }

    if (start_value->next_segment_mesh_pair_index == end_value->prev_segment_mesh_pair_index)
    {
      return live_pair(start_value->next_segment_mesh_pair_index);
    }

    // Try to recover, perhaps they are swapped, but issue a warning since this should not happen
    if (start_value->prev_segment_mesh_pair_index == end_value->next_segment_mesh_pair_index)
    {
      KINDS_WARNING("resolveIntersectionMeshPairIndex: start/end intersection mesh pair index mismatch but recoverable "
                    "by swapping start/end (start_next="
        << start_value->next_segment_mesh_pair_index << ", start_prev=" << start_value->prev_segment_mesh_pair_index
        << ", end_next=" << end_value->next_segment_mesh_pair_index
        << ", end_prev=" << end_value->prev_segment_mesh_pair_index << ", voronoi_cell_id=" << voronoi_cell_id
        << ", event_time=" << event_time << ").");
      return live_pair(start_value->prev_segment_mesh_pair_index);
    }

    // One-sided link after radius boundary reset or stale crossing cleanup: trust the populated side.
    if (start_value->next_segment_mesh_pair_index != static_cast<size_t>(-1)
      && end_value->prev_segment_mesh_pair_index == static_cast<size_t>(-1))
    {
      KINDS_WARNING("resolveIntersectionMeshPairIndex: end prev unset; recovering from start next (start_next="
        << start_value->next_segment_mesh_pair_index << ", voronoi_cell_id=" << voronoi_cell_id
        << ", event_time=" << event_time << ").");
      return live_pair(start_value->next_segment_mesh_pair_index);
    }
    if (end_value->prev_segment_mesh_pair_index != static_cast<size_t>(-1)
      && start_value->next_segment_mesh_pair_index == static_cast<size_t>(-1))
    {
      KINDS_WARNING("resolveIntersectionMeshPairIndex: start next unset; recovering from end prev (end_prev="
        << end_value->prev_segment_mesh_pair_index << ", voronoi_cell_id=" << voronoi_cell_id
        << ", event_time=" << event_time << ").");
      return live_pair(end_value->prev_segment_mesh_pair_index);
    }

    std::ostringstream oss;
    oss << "resolveIntersectionMeshPairIndex: start/end intersection mesh pair index mismatch (start_next="
        << start_value->next_segment_mesh_pair_index << ", end_prev=" << end_value->prev_segment_mesh_pair_index
        << ", voronoi_cell_id=" << voronoi_cell_id << ", event_time=" << event_time << ", de=" << d_edge_id
        << ", s_ve=" << start_value->voronoi_edge_id << ", e_ve=" << end_value->voronoi_edge_id
        << ", edge_links=" << formatDelaunayEdgeCrossingMeshPairLinkSequence(d_edge_id) << ").";
    KINDS_ERROR(oss.str());
    throw std::runtime_error(oss.str());
  }

  // Now handle cases where one is null
  if (start_intersection.has_value())
  {
    // Verify that this is indeed the last intersection in the delaunay edge
    size_t delaunay_edge_id = start_intersection.value()->delaunay_edge_id;
    auto d_ref = start_intersection.value()->delaunay_ref;
    if (!d_ref.has_value())
    {
      throw std::runtime_error("resolveIntersectionMeshPairIndex: start_intersection has unset delaunay_ref");
    }
    auto next_ref = std::next(*d_ref);
    if (next_ref != kin_del.getCrossingData().delaunay_edge_intersections[delaunay_edge_id].end())
    {
      std::ostringstream oss;
      oss << "resolveIntersectionMeshPairIndex: start_intersection is not the last one on its Delaunay edge "
             "(total_on_edge="
          << kin_del.getCrossingData().delaunay_edge_intersections[delaunay_edge_id].size()
          << ", voronoi_cell_id=" << voronoi_cell_id << ", event_time=" << event_time << ").";
      KINDS_ERROR(oss.str());
      throw std::runtime_error(oss.str());
    }

    // check if value is valid
    if (start_intersection.value()->next_segment_mesh_pair_index == static_cast<size_t>(-1))
    {
      size_t prev_index = start_intersection.value()->prev_segment_mesh_pair_index;
      KINDS_WARNING("resolveIntersectionMeshPairIndex: start_intersection has no linked mesh pair "
                    "(voronoi_cell_id="
        << voronoi_cell_id << ", event_time=" << event_time << ", prev_segment_mesh_pair_index=" << prev_index
        << "); skipping interval.");
      return static_cast<size_t>(-1);
    }

    return live_pair(start_intersection.value()->next_segment_mesh_pair_index);
  }
  else // if(end_intersection.has_value())
  {
    // Verify that this is indeed the first intersection in the delaunay edge
    size_t delaunay_edge_id = end_intersection.value()->delaunay_edge_id;
    auto d_ref = end_intersection.value()->delaunay_ref;
    if (!d_ref.has_value())
    {
      throw std::runtime_error("resolveIntersectionMeshPairIndex: end_intersection has unset delaunay_ref");
    }
    if (*d_ref != kin_del.getCrossingData().delaunay_edge_intersections[delaunay_edge_id].begin())
    {
      std::ostringstream oss;
      oss << "resolveIntersectionMeshPairIndex: end_intersection is not the first one on its Delaunay edge "
             "(voronoi_cell_id="
          << voronoi_cell_id << ", event_time=" << event_time << ").";
      KINDS_ERROR(oss.str());
      throw std::runtime_error(oss.str());
    }

    // check if value is valid
    if (end_intersection.value()->prev_segment_mesh_pair_index == static_cast<size_t>(-1))
    {
      size_t next_index = end_intersection.value()->next_segment_mesh_pair_index;
      KINDS_WARNING("resolveIntersectionMeshPairIndex: end_intersection has no linked mesh pair "
                    "(voronoi_cell_id="
        << voronoi_cell_id << ", event_time=" << event_time << ", next_segment_mesh_pair_index=" << next_index
        << "); skipping interval.");
      return static_cast<size_t>(-1);
    }

    return live_pair(end_intersection.value()->prev_segment_mesh_pair_index);
  }
}

size_t SegmentBuilder::startNewMeshFromIntersections(size_t voronoi_cell_id, double t,
  std::optional<KineticDelaunay::CrossingData::EdgeIntersectionRef> start_intersection,
  std::optional<KineticDelaunay::CrossingData::EdgeIntersectionRef> end_intersection, bool reuse_existing_pair_and_mesh,
  BoundaryEventType event_type, BoundarySegmentAction segment_action, bool force_single_seed_vertex,
  const RadiusBoundaryTransitionShiftContext* boundary_transition_shift)
{
  if (!start_intersection.has_value() && !end_intersection.has_value())
  {
    throw std::runtime_error("startNewMeshFromIntersections requires at least one intersection reference.");
  }

  const auto& graph = kin_del.getGraph();

  // Interval semantics (start/end of the strip) and mesh topology are passed explicitly; only `SegmentBuilder`
  // state is captured via `[this]`.
  const auto finish_for_delaunay_edge
    = [this, boundary_transition_shift](const HalfEdgeDelaunayGraph& graph, size_t delaunay_edge_id,
        size_t voronoi_cell_id, double t,
        const std::optional<KineticDelaunay::CrossingData::EdgeIntersectionRef>& interval_start_crossing,
        const std::optional<KineticDelaunay::CrossingData::EdgeIntersectionRef>& interval_end_crossing,
        bool reuse_existing_pair_and_mesh, BoundaryEventType event_type, BoundarySegmentAction segment_action,
        bool force_single_seed_vertex) -> size_t
  {
    // Directed half-edges for this undirected Delaunay edge: even index is one orientation, odd is the twin.
    const size_t he_even = 2 * delaunay_edge_id;
    const size_t he_odd = he_even + 1;
    if (he_odd >= graph.halfEdgeSlotCount())
    {
      return static_cast<size_t>(-1);
    }

    // `MeshingData` stores the *inside* boundary half-edge for each endpoint that is an actual crossing (not an open
    // site). On a boundary Delaunay edge, one directed half-edge points outward from the meshed component; we store the
    // opposite.
    int inside_boundary_he_id = -1;
    if (kin_del.isOnComponentBoundary(he_even))
    {
      const bool even_is_outside = kin_del.isOnComponentBoundaryOutside(he_even);
      inside_boundary_he_id = static_cast<int>(even_is_outside ? he_odd : he_even);
    }

    size_t intersection_pair_index = static_cast<size_t>(-1);
    bool created_new_pair = false;
    if (reuse_existing_pair_and_mesh)
    {
      // Pair lookup keys off the same interval tuple the caller passed (cell, start/end refs, time).
      intersection_pair_index
        = resolveIntersectionMeshPairIndex(voronoi_cell_id, interval_start_crossing, interval_end_crossing, t);
      if (intersection_pair_index != static_cast<size_t>(-1)
        && intersection_pair_index >= intersection_segment_mesh_pairs.size())
      {
        intersection_pair_index = static_cast<size_t>(-1);
      }
      if (intersection_pair_index == static_cast<size_t>(-1))
      {
        return static_cast<size_t>(-1);
      }
    }

    const auto& he = graph.halfEdge(he_even);
    const auto& twin_he = graph.halfEdge(he_odd);
    const size_t owner_segment_id
      = (voronoi_cell_id < strand_to_segment_indices.size() && !strand_to_segment_indices[voronoi_cell_id].empty())
      ? strand_to_segment_indices[voronoi_cell_id].back()
      : static_cast<size_t>(-1);
    MeshStructure::SegmentMeshPair segment_mesh_pair;
    segment_mesh_pair.segment_index0 = he.origin == -1 ? -1 : strand_to_segment_indices[he.origin].back();
    segment_mesh_pair.segment_index1 = twin_he.origin == -1 ? -1 : strand_to_segment_indices[twin_he.origin].back();
    if (intersection_pair_index == static_cast<size_t>(-1))
    {
      created_new_pair = true;
      intersection_pair_index = intersection_segment_mesh_pairs.size();
      intersection_segment_mesh_pairs.push_back(segment_mesh_pair);
    }
    else
    {
      intersection_segment_mesh_pairs[intersection_pair_index] = segment_mesh_pair;
    }

    if (intersection_pair_index < intersection_mesh_pair_metadata.size())
    {
      intersection_mesh_pair_metadata[intersection_pair_index]
        = MeshStructure::IntersectionMeshPairMetadata { voronoi_cell_id, owner_segment_id,
            interval_start_crossing.has_value() ? interval_start_crossing.value()->delaunay_edge_id
                                                : static_cast<size_t>(-1),
            interval_end_crossing.has_value() ? interval_end_crossing.value()->delaunay_edge_id
                                              : static_cast<size_t>(-1) };
    }
    else
    {
      intersection_mesh_pair_metadata.resize(intersection_pair_index + 1);
      intersection_mesh_pair_metadata[intersection_pair_index]
        = MeshStructure::IntersectionMeshPairMetadata { voronoi_cell_id, owner_segment_id,
            interval_start_crossing.has_value() ? interval_start_crossing.value()->delaunay_edge_id
                                                : static_cast<size_t>(-1),
            interval_end_crossing.has_value() ? interval_end_crossing.value()->delaunay_edge_id
                                              : static_cast<size_t>(-1) };
    }

    const size_t component_vertex = std::max(he.origin, twin_he.origin);
    const size_t component_id = kin_del.component_data.component_map[component_vertex];
    std::vector<bool> he_visited(graph.halfEdgeSlotCount(), false);
    updateBoundary(t, he_visited, component_id);
    auto& boundary_polygon = kin_del.component_data.component_boundaries[component_id][0];
    auto& centroid = kin_del.component_data.component_centroids[component_id];

    const bool reuse_in_place = !created_new_pair && intersection_pair_index < intersection_meshes.size();
    if (!reuse_in_place)
    {
      VoronoiMesh mesh_local(MeshletExportMaterialNames);
      configureMeshletStorage(mesh_local);
      intersection_meshes.push_back(std::move(mesh_local));
      intersection_mesh_raw_uvs.emplace_back();
      boundary_meshlet_completed_.push_back(false);
      intersection_meshlet_export_suffixes.push_back(std::string("_intersection_d") + std::to_string(delaunay_edge_id));
    }
    VoronoiMesh& mesh = intersection_meshes[intersection_pair_index];

    const std::string base_boundary_meta = composeBoundaryMetadata(event_type, segment_action);

    struct CrossingOrSiteEndpoint
    {
      glm::dvec3 position {};
      std::optional<RadiusBoundaryTransitionCrossingPlacement> placement {};
      std::optional<size_t> snap_voronoi_vertex_id {};
      bool explicit_profile_position = false;
      std::optional<RadiusTransitionProjection> projection {};
      std::optional<glm::dvec3> explicit_mesh_position {};
      RadiusShiftedSiteCacheEntry* shifted_site_cache = nullptr;
    };

    auto crossing_or_site_endpoint
      = [&](const std::optional<KineticDelaunay::CrossingData::EdgeIntersectionRef>& ref) -> CrossingOrSiteEndpoint
    {
      if (ref.has_value())
      {
        const RadiusBoundaryTransitionCrossingPlacement placement
          = resolveRadiusBoundaryTransitionCrossingPlacement(t, ref.value(), boundary_transition_shift);
        const glm::dvec3 old_pos = crossingProfilePosition(t, placement.conceptual_intersection);
        const glm::dvec3 pos = crossingProfilePositionFromPlacement(t, placement);
        if (placement.positionDiffersFromConceptual())
        {
          logRadiusBoundaryTransitionVertexShift("startNewMeshFromIntersections_interval", t,
            placement.conceptual_intersection, placement.position_intersection, old_pos, pos);
        }
        return { pos, placement, placement.snap_voronoi_vertex_id, placement.explicit_profile_position.has_value(),
          placement.projection, std::nullopt, nullptr };
      }
      if (RadiusShiftedSiteCacheEntry* cache_entry
        = getOrFillRadiusShiftedSiteCache(t, voronoi_cell_id, boundary_transition_shift);
        cache_entry != nullptr)
      {
        const glm::dvec2 p_site = kin_del.getPointAt(t, voronoi_cell_id, false, false);
        KINDS_DEBUG("Radius boundary transition [startNewMeshFromIntersections_site]: cell="
          << voronoi_cell_id << " strip_d=" << delaunay_edge_id << " t=" << t << " old=(" << p_site.x << "," << p_site.y
          << ") new=(" << cache_entry->placement.position.x << "," << cache_entry->placement.position.y
          << ") cache_mesh=" << (cache_entry->explicit_mesh_position.has_value() ? "yes" : "no"));
        return { cache_entry->placement.position, std::nullopt, cache_entry->placement.snap_voronoi_vertex_id, false,
          cache_entry->placement.projection, cache_entry->explicit_mesh_position, cache_entry };
      }
      const glm::dvec2 p = kin_del.getPointAt(t, voronoi_cell_id, false, false);
      return { glm::dvec3(p.x, p.y, t), std::nullopt, std::nullopt, false, std::nullopt, std::nullopt, nullptr };
    };

    const CrossingOrSiteEndpoint start_endpoint = crossing_or_site_endpoint(interval_start_crossing);
    const CrossingOrSiteEndpoint end_endpoint = crossing_or_site_endpoint(interval_end_crossing);
    const glm::dvec3& start_pos = start_endpoint.position;
    const glm::dvec3& end_pos = end_endpoint.position;
    const auto& start_placement_opt = start_endpoint.placement;
    const auto& end_placement_opt = end_endpoint.placement;
    const auto endpoint_runtime_info = [](const CrossingOrSiteEndpoint& endpoint)
    {
      MeshletVertexRuntimeInfo info;
      info.radius_shift_explicit_profile_position = endpoint.explicit_profile_position;
      info.radius_transition_projection = endpoint.projection;
      info.explicit_mesh_position = endpoint.explicit_mesh_position;
      if (endpoint.placement.has_value())
      {
        info.position_intersection = endpoint.placement->position_intersection;
        info.conceptual_intersection = endpoint.placement->conceptual_intersection;
      }
      return info;
    };

    auto remember_shifted_site_mesh = [&](const CrossingOrSiteEndpoint& endpoint, size_t vertex_index)
    {
      if (endpoint.shifted_site_cache == nullptr)
      {
        return;
      }
      noteRadiusShiftedSiteMeshPosition(
        *endpoint.shifted_site_cache, mesh.getVertices()[vertex_index], "startNewMeshFromIntersections");
    };

    std::string boundary_start_meta_left;
    std::string boundary_start_meta_right;
    std::string boundary_start_meta_uniform;
    if (store_mesh_metadata)
    {
      const char* start_source = interval_start_crossing.has_value() ? "intersection" : "site";
      const char* end_source = interval_end_crossing.has_value() ? "intersection" : "site";
      auto make_endpoint_meta = [&](const char* source, const char* pos,
                                  const std::optional<RadiusBoundaryTransitionCrossingPlacement>& placement_opt)
      {
        MetadataBuilder builder;
        builder.addString("event_type", boundaryEventTypeToString(event_type)).addString("source", source);
        if (source == std::string("intersection") && placement_opt.has_value())
        {
          const auto& placement = placement_opt.value();
          return intersectionCrossingVertexMetadata(builder.build(), placement.conceptual_intersection,
            placement.position_intersection, pos, placement.explicit_profile_position.has_value());
        }
        return builder.build();
      };
      boundary_start_meta_left = make_endpoint_meta(start_source, "left", start_placement_opt);
      boundary_start_meta_right = make_endpoint_meta(end_source, "right", end_placement_opt);
      boundary_start_meta_uniform = make_endpoint_meta(start_source, "uniform", start_placement_opt);
    }
    const double dx = start_pos.x - end_pos.x;
    const double dy = start_pos.y - end_pos.y;
    const double dz = start_pos.z - end_pos.z;
    const bool same_endpoint = force_single_seed_vertex || (dx * dx + dy * dy + dz * dz) <= 1e-20;

    // First vertex = interval start (left tag); second = interval end (right tag). Uniform when the strip collapses to
    // one vertex.
    const std::string& start_vertex_meta = same_endpoint ? boundary_start_meta_uniform : boundary_start_meta_left;
    const size_t start_vertex_index = addMeshletVertex(mesh, boundary_polygon, centroid, start_pos, voronoi_cell_id, t,
      false, start_endpoint.snap_voronoi_vertex_id, start_vertex_meta,
      same_endpoint ? std::optional<glm::dvec3>(glm::dvec3(1.0, 0.0, 1.0))
                    : std::optional<glm::dvec3>(glm::dvec3(1.0, 0.0, 0.0)),
      endpoint_runtime_info(start_endpoint));
    remember_shifted_site_mesh(start_endpoint, start_vertex_index);
    const size_t end_vertex_index = same_endpoint
      ? start_vertex_index
      : addMeshletVertex(mesh, boundary_polygon, centroid, end_pos, voronoi_cell_id, t, false,
          end_endpoint.snap_voronoi_vertex_id, boundary_start_meta_right, glm::dvec3(0.0, 0.0, 1.0),
          endpoint_runtime_info(end_endpoint));
    if (!same_endpoint)
    {
      remember_shifted_site_mesh(end_endpoint, end_vertex_index);
    }

    std::list<MeshingData> local_segments;
    MeshingData seg { static_cast<int>(start_vertex_index), static_cast<int>(end_vertex_index),
      interval_start_crossing.has_value() ? inside_boundary_he_id : -1,
      interval_end_crossing.has_value() ? inside_boundary_he_id : -1 };
    seg.start_crossing = interval_start_crossing;
    seg.end_crossing = interval_end_crossing;
    if (reuse_in_place)
    {
      if (intersection_pair_index >= intersection_mesh_pair_last_left_and_right_vertex.size())
      {
        intersection_mesh_pair_last_left_and_right_vertex.resize(intersection_pair_index + 1);
      }
      auto& segments_for_pair = intersection_mesh_pair_last_left_and_right_vertex[intersection_pair_index];
      segments_for_pair.clear();
      segments_for_pair.emplace_back(std::move(seg));
    }
    else
    {
      local_segments.emplace_back(std::move(seg));
    }

    // Prev/next on `CrossingData` follow list order on the Delaunay edge for `[ref,ref]`; one-null cases use
    // `voronoi_cell_id` vs even-half-edge origin inside `writeIntersectionPairLinks` (see there). Arguments mirror
    // interval start/end crossings.
    writeIntersectionPairLinks(
      intersection_pair_index, voronoi_cell_id, interval_start_crossing, interval_end_crossing, t);

    if (reuse_in_place)
    {
      mesh.setCreationKineticTime(t);
    }
    else
    {
      intersection_mesh_pair_last_left_and_right_vertex.emplace_back(std::move(local_segments));
      mesh.setCreationKineticTime(t);
    }

    return intersection_pair_index;
  };

  // Both crossings `[ref, ref]` — must agree on the same Delaunay edge.
  if (start_intersection.has_value() && end_intersection.has_value())
  {
    if (start_intersection.value()->delaunay_edge_id != end_intersection.value()->delaunay_edge_id)
    {
      return static_cast<size_t>(-1);
    }
    return finish_for_delaunay_edge(graph, start_intersection.value()->delaunay_edge_id, voronoi_cell_id, t,
      start_intersection, end_intersection, reuse_existing_pair_and_mesh, event_type, segment_action,
      force_single_seed_vertex);
  }

  // Crossing at start, open site at end `[ref, null]`.
  if (start_intersection.has_value())
  {
    return finish_for_delaunay_edge(graph, start_intersection.value()->delaunay_edge_id, voronoi_cell_id, t,
      start_intersection, end_intersection, reuse_existing_pair_and_mesh, event_type, segment_action,
      force_single_seed_vertex);
  }

  // Open site at start, crossing at end `[null, ref]`.
  return finish_for_delaunay_edge(graph, end_intersection.value()->delaunay_edge_id, voronoi_cell_id, t,
    start_intersection, end_intersection, reuse_existing_pair_and_mesh, event_type, segment_action,
    force_single_seed_vertex);
}

void SegmentBuilder::finishMeshFromIntersections(size_t voronoi_cell_id, double t,
  std::optional<KineticDelaunay::CrossingData::EdgeIntersectionRef> start_intersection,
  std::optional<KineticDelaunay::CrossingData::EdgeIntersectionRef> end_intersection, BoundaryEventType event_type,
  BoundarySegmentAction segment_action, const RadiusBoundaryTransitionShiftContext* boundary_transition_shift)
{
  // voronoi_cell_id is the strand site (Delaunay vertex), not a Voronoi vertex (Delaunay face).
  // Validate circumcenters only via the intersection Voronoi-edge endpoints below.
  if (start_intersection.has_value())
  {
    requireLiveRegisteredVoronoiEdgeEndpoints(
      start_intersection.value()->voronoi_edge_id, "finishMeshFromIntersections(start_intersection)");
  }
  if (end_intersection.has_value())
  {
    requireLiveRegisteredVoronoiEdgeEndpoints(
      end_intersection.value()->voronoi_edge_id, "finishMeshFromIntersections(end_intersection)");
  }

  if (!start_intersection.has_value() && !end_intersection.has_value())
  {
    throw std::runtime_error("finishMeshFromIntersections requires at least one intersection reference.");
  }

  const size_t intersection_pair_index
    = resolveIntersectionMeshPairIndex(voronoi_cell_id, start_intersection, end_intersection, t);

  if (intersection_pair_index == static_cast<size_t>(-1) || intersection_pair_index >= intersection_meshes.size()
    || intersection_pair_index >= intersection_mesh_pair_last_left_and_right_vertex.size())
  {
    std::ostringstream oss;
    oss << "finishMeshFromIntersections: failed to resolve mesh pair for extension (cell=" << voronoi_cell_id
        << ", t=" << t << ", start_ref=";
    if (start_intersection.has_value())
    {
      oss << "{d_edge=" << start_intersection.value()->delaunay_edge_id
          << ", v_edge=" << start_intersection.value()->voronoi_edge_id
          << ", param=" << start_intersection.value()->delaunay_edge_param << "}";
    }
    else
    {
      oss << "null";
    }
    oss << ", end_ref=";
    if (end_intersection.has_value())
    {
      oss << "{d_edge=" << end_intersection.value()->delaunay_edge_id
          << ", v_edge=" << end_intersection.value()->voronoi_edge_id
          << ", param=" << end_intersection.value()->delaunay_edge_param << "}";
    }
    else
    {
      oss << "null";
    }
    oss << ").";
    KINDS_WARNING(oss.str());
    return;
  }

  auto& segs = intersection_mesh_pair_last_left_and_right_vertex[intersection_pair_index];
  if (segs.empty())
  {
    std::ostringstream oss;
    oss << "finishMeshFromIntersections: missing segment state for pair " << intersection_pair_index
        << " (cell=" << voronoi_cell_id << ", t=" << t << ").";
    throw std::runtime_error(oss.str());
  }
  auto& seg = segs.front();
  if (seg.mesh_start_vertex_id < 0 || seg.mesh_end_vertex_id < 0)
  {
    std::ostringstream oss;
    oss << "finishMeshFromIntersections: missing last-left/last-right indices for pair " << intersection_pair_index
        << " (mesh_start_vertex_id=" << seg.mesh_start_vertex_id << ", mesh_end_vertex_id=" << seg.mesh_end_vertex_id
        << ", cell=" << voronoi_cell_id << ", t=" << t << ").";
    throw std::runtime_error(oss.str());
  }

  auto& mesh = intersection_meshes[intersection_pair_index];

  size_t component_id = kin_del.component_data.component_map[voronoi_cell_id];
  std::vector<bool> he_visited(kin_del.getGraph().halfEdgeSlotCount(), false);
  updateBoundary(t, he_visited, component_id);
  auto& boundary_polygon = kin_del.component_data.component_boundaries[component_id][0];
  auto centroid = kin_del.component_data.component_centroids[component_id];

  const size_t strip_delaunay_edge_id_for_site = start_intersection.has_value()
    ? start_intersection.value()->delaunay_edge_id
    : end_intersection.value()->delaunay_edge_id;

  std::optional<RadiusBoundaryTransitionCrossingPlacement> start_placement_opt;
  std::optional<RadiusBoundaryTransitionCrossingPlacement> end_placement_opt;
  std::optional<size_t> start_snap_voronoi_vertex_id;
  std::optional<size_t> end_snap_voronoi_vertex_id;

  struct CrossingOrSiteEndpoint
  {
    glm::dvec3 position {};
    std::optional<RadiusBoundaryTransitionCrossingPlacement> placement {};
    std::optional<size_t> snap_voronoi_vertex_id {};
    bool explicit_profile_position = false;
    std::optional<RadiusTransitionProjection> projection {};
    std::optional<glm::dvec3> explicit_mesh_position {};
    RadiusShiftedSiteCacheEntry* shifted_site_cache = nullptr;
  };

  auto crossing_or_site_endpoint
    = [&](const std::optional<KineticDelaunay::CrossingData::EdgeIntersectionRef>& input_ref) -> CrossingOrSiteEndpoint
  {
    if (!input_ref.has_value())
    {
      if (RadiusShiftedSiteCacheEntry* cache_entry
        = getOrFillRadiusShiftedSiteCache(t, voronoi_cell_id, boundary_transition_shift);
        cache_entry != nullptr)
      {
        const glm::dvec2 p_site = kin_del.getPointAt(t, voronoi_cell_id, false, false);
        KINDS_DEBUG("Radius boundary transition [finishMeshFromIntersections_site]: cell="
          << voronoi_cell_id << " strip_d=" << strip_delaunay_edge_id_for_site << " t=" << t << " old=(" << p_site.x
          << "," << p_site.y << ") new=(" << cache_entry->placement.position.x << ","
          << cache_entry->placement.position.y
          << ") cache_mesh=" << (cache_entry->explicit_mesh_position.has_value() ? "yes" : "no"));
        return { cache_entry->placement.position, std::nullopt, cache_entry->placement.snap_voronoi_vertex_id, false,
          cache_entry->placement.projection, cache_entry->explicit_mesh_position, cache_entry };
      }
      const glm::dvec2 p = kin_del.getPointAt(t, voronoi_cell_id, false, false);
      return { glm::dvec3(p.x, p.y, t), std::nullopt, std::nullopt, false, std::nullopt, std::nullopt, nullptr };
    }

    const RadiusBoundaryTransitionCrossingPlacement placement
      = resolveRadiusBoundaryTransitionCrossingPlacement(t, input_ref.value(), boundary_transition_shift);
    const glm::dvec3 old_pos = crossingProfilePosition(t, placement.conceptual_intersection);
    const glm::dvec3 pos = crossingProfilePositionFromPlacement(t, placement);
    if (placement.positionDiffersFromConceptual())
    {
      logRadiusBoundaryTransitionVertexShift("finishMeshFromIntersections_interval", t,
        placement.conceptual_intersection, placement.position_intersection, old_pos, pos);
    }
    return { pos, placement, placement.snap_voronoi_vertex_id, placement.explicit_profile_position.has_value(),
      placement.projection, std::nullopt, nullptr };
  };

  const CrossingOrSiteEndpoint start_endpoint = crossing_or_site_endpoint(start_intersection);
  const CrossingOrSiteEndpoint end_endpoint = crossing_or_site_endpoint(end_intersection);
  const glm::dvec3 new_start_pos = start_endpoint.position;
  const glm::dvec3 new_end_pos = end_endpoint.position;
  start_placement_opt = start_endpoint.placement;
  end_placement_opt = end_endpoint.placement;
  start_snap_voronoi_vertex_id = start_endpoint.snap_voronoi_vertex_id;
  end_snap_voronoi_vertex_id = end_endpoint.snap_voronoi_vertex_id;
  const auto endpoint_runtime_info = [](const CrossingOrSiteEndpoint& endpoint)
  {
    MeshletVertexRuntimeInfo info;
    info.radius_shift_explicit_profile_position = endpoint.explicit_profile_position;
    info.radius_transition_projection = endpoint.projection;
    info.explicit_mesh_position = endpoint.explicit_mesh_position;
    if (endpoint.placement.has_value())
    {
      info.position_intersection = endpoint.placement->position_intersection;
      info.conceptual_intersection = endpoint.placement->conceptual_intersection;
    }
    return info;
  };

  auto remember_shifted_site_mesh = [&](const CrossingOrSiteEndpoint& endpoint, size_t vertex_index)
  {
    if (endpoint.shifted_site_cache == nullptr)
    {
      return;
    }
    noteRadiusShiftedSiteMeshPosition(
      *endpoint.shifted_site_cache, mesh.getVertices()[vertex_index], "finishMeshFromIntersections");
  };

  std::string boundary_finish_meta;
  std::string boundary_finish_meta_left;
  std::string boundary_finish_meta_right;
  std::string boundary_finish_meta_uniform;
  if (store_mesh_metadata)
  {
    const char* start_source = start_intersection.has_value() ? "intersection" : "site";
    const char* end_source = end_intersection.has_value() ? "intersection" : "site";
    auto make_intersection_meta = [&](const char* source, const char* pos,
                                    const std::optional<RadiusBoundaryTransitionCrossingPlacement>& placement_opt)
    {
      MetadataBuilder builder;
      builder.addString("event_type", boundaryEventTypeToString(event_type)).addString("source", source);
      if (source == std::string("intersection") && placement_opt.has_value())
      {
        const auto& placement = placement_opt.value();
        return intersectionCrossingVertexMetadata(builder.build(), placement.conceptual_intersection,
          placement.position_intersection, pos, placement.explicit_profile_position.has_value());
      }
      return builder.build();
    };
    boundary_finish_meta = make_intersection_meta(start_source, "uniform", start_placement_opt);
    boundary_finish_meta_left = make_intersection_meta(start_source, "left", start_placement_opt);
    boundary_finish_meta_right = make_intersection_meta(end_source, "right", end_placement_opt);
    boundary_finish_meta_uniform = boundary_finish_meta;
  }
  const double fdx = new_start_pos.x - new_end_pos.x;
  const double fdy = new_start_pos.y - new_end_pos.y;
  const double fdz = new_start_pos.z - new_end_pos.z;
  const bool collapsed_finish_endpoints = (fdx * fdx + fdy * fdy + fdz * fdz) <= 1e-20;
  const size_t new_start_vertex_index = addMeshletVertex(mesh, boundary_polygon, centroid, new_start_pos,
    voronoi_cell_id, t, false, start_snap_voronoi_vertex_id,
    collapsed_finish_endpoints ? boundary_finish_meta_uniform : boundary_finish_meta_left,
    collapsed_finish_endpoints ? std::optional<glm::dvec3>(glm::dvec3(1.0, 0.0, 1.0))
                               : std::optional<glm::dvec3>(glm::dvec3(1.0, 0.0, 0.0)),
    endpoint_runtime_info(start_endpoint));
  remember_shifted_site_mesh(start_endpoint, new_start_vertex_index);
  const size_t new_end_vertex_index = collapsed_finish_endpoints
    ? new_start_vertex_index
    : addMeshletVertex(mesh, boundary_polygon, centroid, new_end_pos, voronoi_cell_id, t, false,
        end_snap_voronoi_vertex_id, boundary_finish_meta_right, glm::dvec3(0.0, 0.0, 1.0),
        endpoint_runtime_info(end_endpoint));
  if (!collapsed_finish_endpoints)
  {
    remember_shifted_site_mesh(end_endpoint, new_end_vertex_index);
  }

  size_t ordered_new_start_vertex_index = new_start_vertex_index;
  size_t ordered_new_end_vertex_index = new_end_vertex_index;
  const int old_fixed_start_id = seg.mesh_start_vertex_id;
  const int old_fixed_end_id = seg.mesh_end_vertex_id;
  const size_t eff_left = intersectionStripEffectiveVertexIndex(seg, true);
  const size_t eff_right = intersectionStripEffectiveVertexIndex(seg, false);

  // Single new anchor (collapsed interval / one mesh vertex for both ends): both flex chains meet that point.
  const bool uniform_finish_targets = (ordered_new_start_vertex_index == ordered_new_end_vertex_index);
  const size_t flex_interp_target = ordered_new_start_vertex_index;
  if (!interpolateFlexibleVerticesAlongEdge(mesh, seg.flexible_left_vertex_ids, static_cast<size_t>(old_fixed_start_id),
        uniform_finish_targets ? flex_interp_target : ordered_new_start_vertex_index))
  {
    KINDS_WARNING(
      "finishMeshFromIntersections: left flexible interpolation failed for pair " << intersection_pair_index << ".");
  }
  else
  {
    seg.flexible_left_vertex_ids.clear();
  }
  if (!interpolateFlexibleVerticesAlongEdge(mesh, seg.flexible_right_vertex_ids, static_cast<size_t>(old_fixed_end_id),
        uniform_finish_targets ? flex_interp_target : ordered_new_end_vertex_index))
  {
    KINDS_WARNING(
      "finishMeshFromIntersections: right flexible interpolation failed for pair " << intersection_pair_index << ".");
  }
  else
  {
    seg.flexible_right_vertex_ids.clear();
  }
  resolveRemainingFlexibleVertices(mesh, seg, "finishMeshFromIntersections");

  const int inside_boundary_he_id = (seg.start_half_edge_id >= 0) ? seg.start_half_edge_id : seg.end_half_edge_id;
  const auto& verts_after = mesh.getVertices();
  bool start_apex_first = false;
  if (eff_left != eff_right)
  {
    const double zl = verts_after[eff_left].z;
    const double zr = verts_after[eff_right].z;
    constexpr double z_eps = 1e-12;
    if (std::abs(zl - zr) > z_eps)
    {
      start_apex_first = zl < zr;
    }
    else
    {
      // Flexible placeholders often share z with strip corners; pick the diagonal from projected order along
      // eff_l–eff_r.
      const glm::dvec2 pl(verts_after[eff_left].x, verts_after[eff_left].y);
      const glm::dvec2 pr(verts_after[eff_right].x, verts_after[eff_right].y);
      const glm::dvec2 ps(verts_after[ordered_new_start_vertex_index].x, verts_after[ordered_new_start_vertex_index].y);
      const glm::dvec2 pe(verts_after[ordered_new_end_vertex_index].x, verts_after[ordered_new_end_vertex_index].y);
      const glm::dvec2 d = pr - pl;
      const double strip_len2 = glm::dot(d, d);
      if (strip_len2 > 1e-24)
      {
        start_apex_first = glm::dot(ps - pl, d) < glm::dot(pe - pl, d);
      }
      else
      {
        start_apex_first = zl < zr;
      }
    }
  }
  if (eff_left == eff_right)
  {
    addBoundaryIntervalTriangleOriented(mesh, ordered_new_start_vertex_index, eff_right, ordered_new_end_vertex_index,
      inside_boundary_he_id, t, boundary_finish_meta);
  }
  else if (start_apex_first)
  {
    addBoundaryIntervalTriangleOriented(
      mesh, eff_left, eff_right, ordered_new_start_vertex_index, inside_boundary_he_id, t, boundary_finish_meta);
    addBoundaryIntervalTriangleOriented(mesh, ordered_new_start_vertex_index, eff_right, ordered_new_end_vertex_index,
      inside_boundary_he_id, t, boundary_finish_meta);
  }
  else
  {
    addBoundaryIntervalTriangleOriented(
      mesh, eff_left, eff_right, ordered_new_end_vertex_index, inside_boundary_he_id, t, boundary_finish_meta);
    addBoundaryIntervalTriangleOriented(mesh, eff_left, ordered_new_end_vertex_index, ordered_new_start_vertex_index,
      inside_boundary_he_id, t, boundary_finish_meta);
  }

  seg.mesh_start_vertex_id = static_cast<int>(ordered_new_start_vertex_index);
  seg.mesh_end_vertex_id = static_cast<int>(ordered_new_end_vertex_index);
}

std::vector<KineticDelaunay::CrossingData::EdgeIntersectionRef> SegmentBuilder::getBoundaryIntersectionsInBoundaryOrder(
  size_t delaunay_edge_id) const
{
  const auto& crossing_data = kin_del.getCrossingData();
  if (delaunay_edge_id >= crossing_data.delaunay_edge_intersections.size())
  {
    std::ostringstream oss;
    oss << "getBoundaryIntersectionsInBoundaryOrder: Delaunay edge out of bounds (edge=" << delaunay_edge_id
        << ", crossing_data_size=" << crossing_data.delaunay_edge_intersections.size() << ").";
    throw std::runtime_error(oss.str());
  }

  std::vector<KineticDelaunay::CrossingData::EdgeIntersectionRef> refs;
  const auto& d_intersections = crossing_data.delaunay_edge_intersections[delaunay_edge_id];
  for (const auto& ref : d_intersections)
  {
    refs.push_back(ref);
  }
  return refs;
}

std::string SegmentBuilder::formatCrossingMeshPairLinkSequence(
  const std::vector<std::pair<size_t, size_t>>& prev_next_pairs)
{
  std::ostringstream oss;
  oss << '(';
  bool first = true;
  for (const auto& [prev_pair, next_pair] : prev_next_pairs)
  {
    if (!first)
    {
      oss << ", ";
    }
    first = false;
    oss << formatMeshPairIndex(prev_pair) << ", " << formatMeshPairIndex(next_pair);
  }
  oss << ')';
  return oss.str();
}

std::string SegmentBuilder::formatDelaunayEdgeCrossingMeshPairLinkSequence(size_t delaunay_edge_id) const
{
  const auto& crossing_data = kin_del.getCrossingData();
  if (delaunay_edge_id >= crossing_data.delaunay_edge_intersections.size())
  {
    return "(out_of_range)";
  }

  const auto& d_list = crossing_data.delaunay_edge_intersections[delaunay_edge_id];
  std::vector<std::pair<size_t, size_t>> prev_next_pairs;
  prev_next_pairs.reserve(d_list.size());
  for (const KineticDelaunay::CrossingData::EdgeIntersectionRef& ref : d_list)
  {
    prev_next_pairs.emplace_back(ref->prev_segment_mesh_pair_index, ref->next_segment_mesh_pair_index);
  }
  return formatCrossingMeshPairLinkSequence(prev_next_pairs);
}

void SegmentBuilder::clearIntersectionMeshPairLinksOnDelaunayEdge(size_t delaunay_edge_id)
{
  auto& crossing_data = kin_del.getCrossingDataMutable();
  if (delaunay_edge_id >= crossing_data.delaunay_edge_intersections.size())
  {
    return;
  }

  for (auto& inter : crossing_data.delaunay_edge_intersections[delaunay_edge_id])
  {
    assignIntersectionMeshPairLink(
      inter, true, static_cast<size_t>(-1), "clearIntersectionMeshPairLinksOnDelaunayEdge:prev");
    assignIntersectionMeshPairLink(
      inter, false, static_cast<size_t>(-1), "clearIntersectionMeshPairLinksOnDelaunayEdge:next");
  }
  maybeLogDiagnosticsMonitoredDelaunayEdgeTrigger(std::numeric_limits<double>::quiet_NaN(),
    "clearIntersectionMeshPairLinksOnDelaunayEdge", delaunay_edge_id, std::nullopt);
}

void SegmentBuilder::validateDelaunayEdgeIntersectionMeshPairLinks(
  size_t delaunay_edge_id, double event_time, const char* context) const
{
  const auto& crossing_data = kin_del.getCrossingData();
  if (delaunay_edge_id >= crossing_data.delaunay_edge_intersections.size())
  {
    return;
  }

  const auto& d_list = crossing_data.delaunay_edge_intersections[delaunay_edge_id];
  if (d_list.size() < 2)
  {
    return;
  }

  const auto format_mesh_pair_index = [](size_t mesh_pair_index) -> std::string
  {
    if (mesh_pair_index == static_cast<size_t>(-1))
    {
      return "unset";
    }
    return std::to_string(mesh_pair_index);
  };

  const auto format_crossing
    = [&](size_t list_index, KineticDelaunay::CrossingData::EdgeIntersectionRef ref) -> std::string
  {
    std::ostringstream oss;
    oss << "[" << list_index << "] v_edge=" << ref->voronoi_edge_id << " d_edge=" << ref->delaunay_edge_id
        << " param=" << ref->delaunay_edge_param
        << " prev=" << format_mesh_pair_index(ref->prev_segment_mesh_pair_index)
        << " next=" << format_mesh_pair_index(ref->next_segment_mesh_pair_index);
    return oss.str();
  };

  size_t list_index = 0;
  for (auto it = d_list.begin(); it != d_list.end(); ++it, ++list_index)
  {
    const auto next_it = std::next(it);
    if (next_it == d_list.end())
    {
      break;
    }

    const size_t next_index = list_index + 1;
    const size_t cur_next = (*it)->next_segment_mesh_pair_index;
    const size_t next_prev = (*next_it)->prev_segment_mesh_pair_index;
    if (cur_next == next_prev)
    {
      continue;
    }

    std::ostringstream oss;
    oss << context << ": Delaunay edge " << delaunay_edge_id
        << " intersection mesh-pair link mismatch at t=" << event_time << " (crossing[" << list_index
        << "].next=" << format_mesh_pair_index(cur_next) << " != crossing[" << next_index
        << "].prev=" << format_mesh_pair_index(next_prev)
        << ") edge_links=" << formatDelaunayEdgeCrossingMeshPairLinkSequence(delaunay_edge_id) << ":\n";

    size_t dump_index = 0;
    for (KineticDelaunay::CrossingData::EdgeIntersectionRef ref : d_list)
    {
      oss << "  " << format_crossing(dump_index, ref) << '\n';
      ++dump_index;
    }

    const std::string msg = oss.str();
    KINDS_ERROR(msg);
    throw std::runtime_error(msg);
  }
}

void SegmentBuilder::assignIntersectionMeshPairLink(KineticDelaunay::CrossingData::EdgeIntersectionRef ref,
  bool is_prev_link, size_t new_pair_index, const char* context, double t)
{
  size_t& slot = is_prev_link ? ref->prev_segment_mesh_pair_index : ref->next_segment_mesh_pair_index;
  const size_t old_pair_index = slot;
  if (old_pair_index == new_pair_index)
  {
    return;
  }
  slot = new_pair_index;

  if (!diagnostics)
  {
    return;
  }

  // Disabled monitor id must never match cleared/unset pair links (also -1).
  if (!KineticDelaunay::isDiagnosticsMonitorIdEnabled(kDiagnosticsMonitoredMeshPairId))
  {
    return;
  }
  const bool assigned_monitored
    = KineticDelaunay::matchesDiagnosticsMonitorId(new_pair_index, kDiagnosticsMonitoredMeshPairId);
  const bool cleared_monitored
    = KineticDelaunay::matchesDiagnosticsMonitorId(old_pair_index, kDiagnosticsMonitoredMeshPairId)
    && !KineticDelaunay::matchesDiagnosticsMonitorId(new_pair_index, kDiagnosticsMonitoredMeshPairId);
  if (!assigned_monitored && !cleared_monitored)
  {
    return;
  }

  std::ostringstream oss;
  oss << "Monitored mesh-pair link " << (assigned_monitored ? "ASSIGNED" : "CLEARED") << " "
      << (is_prev_link ? "prev" : "next") << "_segment_mesh_pair_index";
  if (assigned_monitored)
  {
    oss << "=" << kDiagnosticsMonitoredMeshPairId;
    if (old_pair_index != static_cast<size_t>(-1))
    {
      oss << " (was " << old_pair_index << ")";
    }
  }
  else
  {
    oss << " (was " << kDiagnosticsMonitoredMeshPairId << " -> "
        << (new_pair_index == static_cast<size_t>(-1) ? -1 : static_cast<long long>(new_pair_index)) << ")";
  }
  oss << " de=" << ref->delaunay_edge_id << " ve=" << ref->voronoi_edge_id << " param=" << ref->delaunay_edge_param;
  if (context != nullptr && context[0] != '\0')
  {
    oss << " [" << context << "]";
  }
  if (std::isfinite(t))
  {
    oss << " t=" << t;
  }
  if (!metadata_callback_phase_.empty())
  {
    oss << " callback=" << metadata_callback_phase_;
  }
  KINDS_DEBUG(oss.str());
  maybeLogDiagnosticsMonitoredDelaunayEdgeTrigger(t, context, ref->delaunay_edge_id, new_pair_index);
}

void SegmentBuilder::invalidateStaleAlphaBoundaryMeshLinksOnEvenHalfEdge(size_t even_half_edge_id)
{
  if (kin_del.isOnComponentBoundary(even_half_edge_id))
  {
    return;
  }

  const size_t odd_half_edge_id = even_half_edge_id ^ 1;
  clearIntersectionMeshPairLinksOnDelaunayEdge(even_half_edge_id / 2);

  if (even_half_edge_id < boundary_mesh_last_left_and_right_vertex.size())
  {
    boundary_mesh_last_left_and_right_vertex[even_half_edge_id] = std::make_pair(-1, -1);
  }
  if (odd_half_edge_id < boundary_mesh_last_left_and_right_vertex.size())
  {
    boundary_mesh_last_left_and_right_vertex[odd_half_edge_id] = std::make_pair(-1, -1);
  }
  if (even_half_edge_id < half_edge_to_boundary_vertex_index.size())
  {
    half_edge_to_boundary_vertex_index[even_half_edge_id] = -1;
  }
  if (odd_half_edge_id < half_edge_to_boundary_vertex_index.size())
  {
    half_edge_to_boundary_vertex_index[odd_half_edge_id] = -1;
  }
}

void SegmentBuilder::invalidateStaleAlphaBoundaryMeshLinksOnTriangleEdgesLeftBoundary(
  const std::array<size_t, 3>& triangle_half_edges)
{
  std::unordered_set<size_t> seen_even_half_edges;
  for (size_t he_id : triangle_half_edges)
  {
    const size_t even_half_edge_id = he_id & ~size_t(1);
    if (!seen_even_half_edges.insert(even_half_edge_id).second)
    {
      continue;
    }
    if (kin_del.isOnComponentBoundary(even_half_edge_id))
    {
      continue;
    }
    const size_t delaunay_edge_id = even_half_edge_id / 2;
    clearIntersectionMeshPairLinksOnDelaunayEdge(delaunay_edge_id);
    KINDS_DEBUG("Radius: cleared stale intersection mesh-pair links on Delaunay edge " << delaunay_edge_id
                                                                                       << " (left alpha boundary).");
  }
}

bool SegmentBuilder::writeOneNullIntersectionPairLinkByNullVertex(size_t intersection_pair_index, size_t null_vertex_id,
  KineticDelaunay::CrossingData::EdgeIntersectionRef ref, bool interval_is_ref_to_null, double t)
{
  const size_t d_edge_id = ref->delaunay_edge_id;
  const size_t he_even = 2 * d_edge_id;
  const size_t he_odd = he_even + 1;
  const auto& graph = kin_del.getGraph();
  if (he_odd >= graph.halfEdgeSlotCount())
  {
    throw std::runtime_error("writeOneNullIntersectionPairLinkByNullVertex: Delaunay half-edge out of bounds.");
  }

  const int even_origin = graph.halfEdge(he_even).origin;
  const bool write_prev = even_origin >= 0 && null_vertex_id == static_cast<size_t>(even_origin);
  const char* link_context = write_prev ? "writeOneNullIntersectionPairLinkByNullVertex:prev"
                                        : "writeOneNullIntersectionPairLinkByNullVertex:next";
  assignIntersectionMeshPairLink(ref, write_prev, intersection_pair_index, link_context, t);
  KINDS_DEBUG("one-null link write " << (interval_is_ref_to_null ? "[ref,null]" : "[null,ref]") << ": de=" << d_edge_id
                                     << " pair=" << intersection_pair_index << " null_vertex=" << null_vertex_id
                                     << " even_origin=" << even_origin << " wrote=" << (write_prev ? "prev" : "next"));
  return write_prev;
}

void SegmentBuilder::writeIntersectionPairLinks(size_t intersection_pair_index, size_t voronoi_cell_id,
  std::optional<KineticDelaunay::CrossingData::EdgeIntersectionRef> start_intersection,
  std::optional<KineticDelaunay::CrossingData::EdgeIntersectionRef> end_intersection, double t)
{
  const auto edge_intersection_ref_desc
    = [](const std::optional<KineticDelaunay::CrossingData::EdgeIntersectionRef>& ref) -> std::string
  {
    if (!ref.has_value())
    {
      return "null";
    }
    return "{d_edge=" + std::to_string(ref.value()->delaunay_edge_id)
      + ", v_edge=" + std::to_string(ref.value()->voronoi_edge_id)
      + ", param=" + std::to_string(ref.value()->delaunay_edge_param) + "}";
  };

  std::string wrote_pair_index_to;

  // Keep prev/next links consistent with the actual order in delaunay_edge_intersections:
  // if y follows x in list order -> x.next = pair, y.prev = pair.
  if (start_intersection.has_value() && end_intersection.has_value())
  {
    const size_t d_edge_id = start_intersection.value()->delaunay_edge_id;
    if (d_edge_id < kin_del.getCrossingData().delaunay_edge_intersections.size())
    {
      const auto& d_list = kin_del.getCrossingData().delaunay_edge_intersections[d_edge_id];
      const auto it_start = std::find(d_list.begin(), d_list.end(), start_intersection.value());
      const auto it_end = std::find(d_list.begin(), d_list.end(), end_intersection.value());
      if (it_start != d_list.end() && it_end != d_list.end() && it_start != it_end)
      {
        KineticDelaunay::CrossingData::EdgeIntersectionRef list_first;
        KineticDelaunay::CrossingData::EdgeIntersectionRef list_second;
        if (listConstIteratorPrecedes(d_list.begin(), d_list.end(), it_start, it_end))
        {
          list_first = start_intersection.value();
          list_second = end_intersection.value();
        }
        else
        {
          list_first = end_intersection.value();
          list_second = start_intersection.value();
        }
        assignIntersectionMeshPairLink(
          list_first, false, intersection_pair_index, "writeIntersectionPairLinks:list_first_next", t);
        assignIntersectionMeshPairLink(
          list_second, true, intersection_pair_index, "writeIntersectionPairLinks:list_second_prev", t);
        wrote_pair_index_to = "list_first->next_segment_mesh_pair_index, list_second->prev_segment_mesh_pair_index";
      }
      else
      {
        assignIntersectionMeshPairLink(start_intersection.value(), false, intersection_pair_index,
          "writeIntersectionPairLinks:start_next_fallback", t);
        assignIntersectionMeshPairLink(
          end_intersection.value(), true, intersection_pair_index, "writeIntersectionPairLinks:end_prev_fallback", t);
        wrote_pair_index_to
          = "start_ref->next_segment_mesh_pair_index, end_ref->prev_segment_mesh_pair_index (list lookup fallback)";
      }
    }
    else
    {
      assignIntersectionMeshPairLink(start_intersection.value(), false, intersection_pair_index,
        "writeIntersectionPairLinks:start_next_missing_list", t);
      assignIntersectionMeshPairLink(
        end_intersection.value(), true, intersection_pair_index, "writeIntersectionPairLinks:end_prev_missing_list", t);
      wrote_pair_index_to = "start_ref->next_segment_mesh_pair_index, end_ref->prev_segment_mesh_pair_index (missing "
                            "crossing list fallback)";
    }
  }
  else
  {
    std::ostringstream w;
    if (start_intersection.has_value())
    {
      const bool prev = writeOneNullIntersectionPairLinkByNullVertex(
        intersection_pair_index, voronoi_cell_id, start_intersection.value(), true, t);
      w << "start_ref->" << (prev ? "prev_segment_mesh_pair_index" : "next_segment_mesh_pair_index");
    }
    if (end_intersection.has_value())
    {
      const bool prev = writeOneNullIntersectionPairLinkByNullVertex(
        intersection_pair_index, voronoi_cell_id, end_intersection.value(), false, t);
      if (start_intersection.has_value())
      {
        w << ", ";
      }
      w << "end_ref->" << (prev ? "prev_segment_mesh_pair_index" : "next_segment_mesh_pair_index");
    }
    wrote_pair_index_to = w.str();
  }

  std::ostringstream log;
  log << "writeIntersectionPairLinks: pair=" << intersection_pair_index << " cell=" << voronoi_cell_id
      << " start_ref=" << edge_intersection_ref_desc(start_intersection)
      << " end_ref=" << edge_intersection_ref_desc(end_intersection) << " wrote_pair_index_to=" << wrote_pair_index_to;
  KINDS_DEBUG(log.str());

  std::optional<size_t> trigger_d_edge;
  if (start_intersection.has_value())
  {
    trigger_d_edge = start_intersection.value()->delaunay_edge_id;
  }
  else if (end_intersection.has_value())
  {
    trigger_d_edge = end_intersection.value()->delaunay_edge_id;
  }
  maybeLogDiagnosticsMonitoredDelaunayEdgeTrigger(
    t, "writeIntersectionPairLinks", trigger_d_edge, intersection_pair_index);
}

size_t SegmentBuilder::determineVoronoiCellForBoundaryIntersectionInterval(size_t delaunay_edge_id,
  std::optional<KineticDelaunay::CrossingData::EdgeIntersectionRef> start_intersection,
  std::optional<KineticDelaunay::CrossingData::EdgeIntersectionRef> end_intersection) const
{
  if (!start_intersection.has_value() && !end_intersection.has_value())
  {
    throw std::runtime_error(
      "determineVoronoiCellForBoundaryIntersectionInterval requires at least one intersection reference.");
  }

  const auto& graph = kin_del.getGraph();
  const size_t he_even = 2 * delaunay_edge_id;
  const size_t he_odd = he_even + 1;
  if (he_odd >= graph.halfEdgeSlotCount())
  {
    std::ostringstream oss;
    oss << "determineVoronoiCellForBoundaryIntersectionInterval: Delaunay edge out of bounds (edge=" << delaunay_edge_id
        << ", half_edge_count=" << graph.halfEdgeSlotCount() << ").";
    throw std::runtime_error(oss.str());
  }

  if (!start_intersection.has_value() && end_intersection.has_value())
  {
    return static_cast<size_t>(graph.halfEdge(he_even).origin);
  }
  if (start_intersection.has_value() && !end_intersection.has_value())
  {
    return static_cast<size_t>(graph.halfEdge(he_odd).origin);
  }

  const size_t ve0 = start_intersection.value()->voronoi_edge_id;
  const size_t ve1 = end_intersection.value()->voronoi_edge_id;
  const size_t ve0_he_even = 2 * ve0;
  const size_t ve0_he_odd = ve0_he_even + 1;
  const size_t ve1_he_even = 2 * ve1;
  const size_t ve1_he_odd = ve1_he_even + 1;
  if (ve0_he_odd >= graph.halfEdgeSlotCount() || ve1_he_odd >= graph.halfEdgeSlotCount())
  {
    std::ostringstream oss;
    oss << "determineVoronoiCellForBoundaryIntersectionInterval: Voronoi edge out of bounds (ve0=" << ve0
        << ", ve1=" << ve1 << ", half_edge_count=" << graph.halfEdgeSlotCount() << ").";
    throw std::runtime_error(oss.str());
  }

  const size_t a0 = static_cast<size_t>(graph.halfEdge(ve0_he_even).origin);
  const size_t b0 = static_cast<size_t>(graph.halfEdge(ve0_he_odd).origin);
  const size_t a1 = static_cast<size_t>(graph.halfEdge(ve1_he_even).origin);
  const size_t b1 = static_cast<size_t>(graph.halfEdge(ve1_he_odd).origin);
  if (a0 == a1 || a0 == b1)
  {
    return a0;
  }
  if (b0 == a1 || b0 == b1)
  {
    return b0;
  }

  std::ostringstream oss;
  oss
    << "determineVoronoiCellForBoundaryIntersectionInterval: no shared Voronoi cell between intersection Voronoi edges "
    << ve0 << " and " << ve1 << " on Delaunay edge " << delaunay_edge_id << ".";
  KINDS_WARNING(oss.str());
  return static_cast<size_t>(-1);
}

size_t kinDS::SegmentBuilder::registerMeshletWithSuffix(
  VoronoiMesh&& mesh, std::string suffix, double creation_kinetic_time)
{
  configureMeshletStorage(mesh);
  if (std::isfinite(creation_kinetic_time))
  {
    mesh.setCreationKineticTime(creation_kinetic_time);
  }
  const size_t index = meshes.size();
  meshes.push_back(std::move(mesh));
  interior_meshlet_completed_.push_back(false);
  meshlet_export_suffixes.push_back(std::move(suffix));
  return index;
}

void kinDS::SegmentBuilder::completeBoundaryMeshSection(size_t he_id, size_t new_left, size_t new_right, double t)
{
  const std::string face_metadata = composeBoundaryMeshFaceMetadata(t, "boundary_section", he_id);
  const auto& last_left_and_right = boundary_mesh_last_left_and_right_vertex[he_id];
  if (last_left_and_right.first != -1)
  {
    // distinguish the case that we have previously flipped an infinite edge that became a boundary edge
    if (half_edge_to_boundary_vertex_index[he_id] == -1)
    {
      // no edge flip
      addBoundaryTriangle(last_left_and_right.first, new_right, new_left, face_metadata);
      if (last_left_and_right.second != -1)
      {
        addBoundaryTriangle(last_left_and_right.second, new_right, last_left_and_right.first, face_metadata);
      }
    }
    else
    {
      assert(last_left_and_right.second != -1);
      // he_id was previously flipped, it's corresponding vertex is no longer part of the boundary
      addBoundaryTriangle(
        last_left_and_right.second, new_right, half_edge_to_boundary_vertex_index[he_id], face_metadata);
      addBoundaryTriangle(
        new_left, last_left_and_right.first, half_edge_to_boundary_vertex_index[he_id], face_metadata);
      addBoundaryTriangle(new_left, half_edge_to_boundary_vertex_index[he_id], new_right, face_metadata);

      // reset the half-edge to boundary vertex index
      half_edge_to_boundary_vertex_index[he_id] = -1;
    }
  }
  else
  {
    // assert(last_left_and_right.second == -1);
  }
}

size_t kinDS::SegmentBuilder::addBoundaryTriangle(
  size_t u, size_t v, size_t w, const std::string& metadata, int material_id)
{
  // check bounds
  if (u >= boundary_mesh.getVertexCount() || v >= boundary_mesh.getVertexCount() || w >= boundary_mesh.getVertexCount())
  {
    KINDS_ERROR("Vertex index " + std::to_string(u) + ", " + std::to_string(v) + ", " + std::to_string(w)
      + " out of boundary mesh range.");
    return -1;
  }

  if (u >= boundary_mesh_raw_uvs.size() || v >= boundary_mesh_raw_uvs.size() || w >= boundary_mesh_raw_uvs.size())
  {
    KINDS_ERROR("Vertex index " + std::to_string(u) + ", " + std::to_string(v) + ", " + std::to_string(w)
      + " out of raw uv range.");
    return -1;
  }

  const AdjustedBoundaryTriangleUvs adjusted = adjustedBoundaryTriangleUvs(
    boundary_mesh_raw_uvs[u], boundary_mesh_raw_uvs[v], boundary_mesh_raw_uvs[w], uv_circum_factor, uv_height_factor);

  size_t uv_index_u = boundary_mesh.addUV(adjusted.u);
  size_t uv_index_v = boundary_mesh.addUV(adjusted.v);
  size_t uv_index_w = boundary_mesh.addUV(adjusted.w);

  /*KINDS_DEBUG("UVs after adjustment: u(" + std::to_string(uv_u[0]) + ", " + std::to_string(uv_u[1]) + "), v(" +
                std::to_string(uv_v[0]) + ", " + std::to_string(uv_v[1]) + "), w(" + std::to_string(uv_w[0]) + ", " +
                std::to_string(uv_w[1]) + ")");*/
  // warnIfTriangleKineticTimesNotInUnitSection(u, v, w, boundary_mesh.getVertices(), "boundary_mesh", 0);
  return boundary_mesh.addTriangle(u, v, w, uv_index_u, uv_index_v, uv_index_w, material_id, metadata);
}

size_t kinDS::SegmentBuilder::addBoundaryIntervalTriangle(
  VoronoiMesh& mesh, size_t u, size_t v, size_t w, const std::string& metadata, int material_id)
{
  if (warnAndSkipIfMeshletCompleted(mesh, "addBoundaryIntervalTriangle", std::numeric_limits<double>::quiet_NaN()))
  {
    return static_cast<size_t>(-1);
  }
  if (u >= mesh.getVertexCount() || v >= mesh.getVertexCount() || w >= mesh.getVertexCount())
  {
    KINDS_ERROR("Vertex index out of boundary-interval mesh range.");
    return static_cast<size_t>(-1);
  }

  const std::vector<glm::dvec2>& raw_uvs = boundaryIntervalRawUvs(mesh);
  if (u >= raw_uvs.size() || v >= raw_uvs.size() || w >= raw_uvs.size())
  {
    KINDS_ERROR("Vertex index out of boundary-interval raw uv range.");
    return static_cast<size_t>(-1);
  }

  const AdjustedBoundaryTriangleUvs adjusted
    = adjustedBoundaryTriangleUvs(raw_uvs[u], raw_uvs[v], raw_uvs[w], uv_circum_factor, uv_height_factor);
  const size_t uv_index_u = mesh.addUV(adjusted.u);
  const size_t uv_index_v = mesh.addUV(adjusted.v);
  const size_t uv_index_w = mesh.addUV(adjusted.w);
  return mesh.addTriangle(u, v, w, uv_index_u, uv_index_v, uv_index_w, material_id, metadata);
}

std::vector<glm::dvec2>& SegmentBuilder::boundaryIntervalRawUvs(VoronoiMesh& mesh)
{
  const std::optional<size_t> mesh_index = intersectionMeshletIndexForMesh(mesh);
  if (!mesh_index.has_value())
  {
    throw std::runtime_error("boundaryIntervalRawUvs: mesh is not a boundary-interval meshlet.");
  }
  if (intersection_mesh_raw_uvs.size() <= mesh_index.value())
  {
    intersection_mesh_raw_uvs.resize(mesh_index.value() + 1);
  }
  return intersection_mesh_raw_uvs[mesh_index.value()];
}

const std::vector<glm::dvec2>& SegmentBuilder::boundaryIntervalRawUvs(const VoronoiMesh& mesh) const
{
  const std::optional<size_t> mesh_index = intersectionMeshletIndexForMesh(mesh);
  if (!mesh_index.has_value() || mesh_index.value() >= intersection_mesh_raw_uvs.size())
  {
    throw std::runtime_error("boundaryIntervalRawUvs: mesh is not a boundary-interval meshlet.");
  }
  return intersection_mesh_raw_uvs[mesh_index.value()];
}

std::optional<glm::dvec2> SegmentBuilder::boundaryIntervalRawUvAtVertex(
  const VoronoiMesh& mesh, size_t vertex_index) const
{
  const std::optional<size_t> mesh_index = intersectionMeshletIndexForMesh(mesh);
  if (!mesh_index.has_value() || mesh_index.value() >= intersection_mesh_raw_uvs.size())
  {
    return std::nullopt;
  }
  const std::vector<glm::dvec2>& raw_uvs = intersection_mesh_raw_uvs[mesh_index.value()];
  if (vertex_index >= raw_uvs.size())
  {
    return std::nullopt;
  }
  const glm::dvec2& raw_uv = raw_uvs[vertex_index];
  if (!std::isfinite(raw_uv.x) || !std::isfinite(raw_uv.y))
  {
    return std::nullopt;
  }
  return raw_uv;
}

void SegmentBuilder::setBoundaryIntervalRawUv(VoronoiMesh& mesh, size_t vertex_index, const glm::dvec2& raw_uv)
{
  std::vector<glm::dvec2>& raw_uvs = boundaryIntervalRawUvs(mesh);
  if (vertex_index >= raw_uvs.size())
  {
    raw_uvs.resize(vertex_index + 1, glm::dvec2(0.0));
  }
  raw_uvs[vertex_index] = raw_uv;
}

void SegmentBuilder::refreshBoundaryIntervalTrianglesIncidentToVertex(VoronoiMesh& mesh, size_t vertex_index)
{
  const std::vector<glm::dvec2>& raw_uvs = boundaryIntervalRawUvs(mesh);
  const auto raw_uv_at_vertex
    = [&](size_t vi) -> glm::dvec2 { return vi < raw_uvs.size() ? raw_uvs[vi] : glm::dvec2(0.0); };

  const std::vector<size_t>& flat_triangles = mesh.getTriangles();
  for (const size_t corner : mesh.findTriangleCorners(vertex_index, true))
  {
    const size_t tri_base = corner - (corner % 3);
    if (tri_base + 2 >= flat_triangles.size())
    {
      continue;
    }
    const size_t u = flat_triangles[tri_base];
    const size_t v = flat_triangles[tri_base + 1];
    const size_t w = flat_triangles[tri_base + 2];
    const AdjustedBoundaryTriangleUvs adjusted = adjustedBoundaryTriangleUvs(
      raw_uv_at_vertex(u), raw_uv_at_vertex(v), raw_uv_at_vertex(w), uv_circum_factor, uv_height_factor);
    mesh.setUV(adjusted.u, tri_base);
    mesh.setUV(adjusted.v, tri_base + 1);
    mesh.setUV(adjusted.w, tri_base + 2);
  }
}

glm::dvec2 SegmentBuilder::boundaryRawUv(const glm::dvec2& delaunay_xy, const glm::dvec2& centroid, double t)
{
  const double angle = std::atan2(centroid.y - delaunay_xy.y, centroid.x - delaunay_xy.x);
  return glm::dvec2(angle / (2.0 * glm::pi<double>()), t);
}

glm::dvec3 SegmentBuilder::interiorMeshUv(const std::vector<BoundaryPoint>& boundary_polygon,
  const glm::dvec2& centroid, const glm::dvec2& delaunay_xy, double t) const
{
  const double relative_distance = relativeDistanceFromCenter(boundary_polygon, centroid, delaunay_xy);
  const double angle = std::atan2(centroid.y - delaunay_xy.y, centroid.x - delaunay_xy.x);
  const double radial_scale = texture_diameter * relative_distance * 0.5;
  return glm::dvec3(0.5 + radial_scale * std::cos(angle), 0.5 + radial_scale * std::sin(angle), t * uv_height_factor);
}

glm::dvec2 SegmentBuilder::meshVirtualShiftForStrand(size_t strand_id, double t) const
{
  // Mesh vertices stay unshifted. Kinetic separation offsets remain in getPointAt / SVG only.
  (void)strand_id;
  (void)t;
  return glm::dvec2(0.0, 0.0);
}

void SegmentBuilder::applyMeshVirtualShiftToProfileVertex(
  glm::dvec3& vertex, glm::dvec2& profile_xy, size_t strand_id, double t, bool& includes_virtual_shift) const
{
  const glm::dvec2 shift = meshVirtualShiftForStrand(strand_id, t);
  if (shift.x == 0.0 && shift.y == 0.0)
  {
    return;
  }

  vertex.x += shift.x;
  vertex.y += shift.y;
  profile_xy += shift;
  includes_virtual_shift = true;
}

void SegmentBuilder::applyUntransformedMeshViewTransform()
{
  const glm::dmat4 transform = VoronoiMesh::profileSpaceSwapYAndZTransform();
  boundary_mesh.applyTransform(transform);
  for (VoronoiMesh& mesh : meshes)
  {
    mesh.applyTransform(transform);
  }
  for (VoronoiMesh& mesh : intersection_meshes)
  {
    mesh.applyTransform(transform);
  }
}

size_t kinDS::SegmentBuilder::addBoundaryVertex(
  glm::dvec3 vertex, glm::dvec2 centroid, size_t strand_id, double t, bool includes_virtual_shift)
{
  const glm::dvec2 delaunay_xy(vertex.x, vertex.y);
  glm::dvec2 profile_xy(vertex.x, vertex.y);
  if (!create_transformed_mesh)
  {
    vertex = computeMeshSiteVertexPosition(vertex, strand_id, t);
    profile_xy = glm::dvec2(vertex.x, vertex.y);
    includes_virtual_shift = false;
  }
  else
  {
    applyMeshVirtualShiftToProfileVertex(vertex, profile_xy, strand_id, t, includes_virtual_shift);
    vertex = transformFromInputBranchToObjectSpace(vertex, strand_id, t);
  }
  const glm::dvec2 raw_uv = boundaryRawUv(delaunay_xy, centroid, t);

  const std::string metadata = store_mesh_metadata ? [&]()
  {
    MetadataBuilder builder;
    builder.addString("event_type", "boundary_vertex").addString("source", "site").addSize("strand_id", strand_id);
    if (!metadata_callback_phase_.empty())
    {
      builder.addString("callback", metadata_callback_phase_);
    }
    return builder.addBool("shift", includes_virtual_shift)
      .addDouble("x", profile_xy.x)
      .addDouble("y", profile_xy.y)
      .addDouble("t", t)
      .build();
  }()
                                                   : std::string {};
  size_t index = boundary_mesh.addVertex(vertex, metadata);
  boundary_vertex_to_strand_id.push_back(strand_id);
  boundary_mesh_raw_uvs.resize(index + 1, glm::dvec2 {});
  boundary_mesh_raw_uvs[index] = raw_uv;

  return index;
}

size_t kinDS::SegmentBuilder::addMeshletTriangle(
  VoronoiMesh& mesh, size_t u, size_t v, size_t w, const std::string& metadata, int material_id)
{
  if (warnAndSkipIfMeshletCompleted(mesh, "addMeshletTriangle", std::numeric_limits<double>::quiet_NaN()))
  {
    return static_cast<size_t>(-1);
  }

  if (mesh.getMaterialNames().empty())
  {
    mesh.setMaterialNames(MeshletExportMaterialNames);
  }
  // warnIfTriangleKineticTimesNotInUnitSection(u, v, w, mesh.getVertices(), "meshlet", material_id);
  const auto vertex_uv = [&](size_t vertex_index) -> glm::dvec3
  {
    if (const std::optional<glm::dvec3> semantic_uv = mesh.vertexSemanticUv(vertex_index); semantic_uv.has_value())
    {
      return semantic_uv.value();
    }
    return glm::dvec3(0.0);
  };
  const size_t uv_index_u = mesh.addUV(vertex_uv(u));
  const size_t uv_index_v = mesh.addUV(vertex_uv(v));
  const size_t uv_index_w = mesh.addUV(vertex_uv(w));
  return mesh.addTriangle(u, v, w, uv_index_u, uv_index_v, uv_index_w, material_id, metadata);
}

size_t kinDS::SegmentBuilder::addBoundaryIntervalTriangleOriented(VoronoiMesh& mesh, size_t u, size_t v, size_t w,
  int inside_boundary_he_id, double t, const std::string& metadata, std::optional<size_t> boundary_meshlet_id)
{
  std::string triangle_metadata = metadata;
  if (store_mesh_metadata)
  {
    MetadataBuilder builder = MetadataBuilder::fromObject(metadata).ensureDouble("t", t);
    std::optional<size_t> meshlet_id = boundary_meshlet_id;
    if (!meshlet_id.has_value())
    {
      meshlet_id = intersectionMeshletIndexForMesh(mesh);
    }
    if (meshlet_id.has_value())
    {
      builder.ensureSize("boundary_meshlet_id", meshlet_id.value());
    }
    triangle_metadata = builder.build();
  }

  if (mesh.getMaterialNames().empty())
  {
    mesh.setMaterialNames(MeshletExportMaterialNames);
  }
  const int boundary_material_id = BoundaryIntervalMeshletMaterialId;
  if (inside_boundary_he_id < 0 || static_cast<size_t>(inside_boundary_he_id) >= kin_del.getGraph().halfEdgeSlotCount())
  {
    return addBoundaryIntervalTriangle(mesh, u, v, w, triangle_metadata, boundary_material_id);
  }

  // `inside_boundary_he_id` is the inside-directed boundary half-edge; its twin is the outside one on the same Delaunay
  // edge.
  const size_t outside_he = static_cast<size_t>(inside_boundary_he_id) ^ 1u;
  if (outside_he >= kin_del.getGraph().halfEdgeSlotCount())
  {
    return addBoundaryIntervalTriangle(mesh, u, v, w, triangle_metadata, boundary_material_id);
  }

  if ((outside_he & 1u) != 0u)
  {
    std::swap(v, w);
  }
  return addBoundaryIntervalTriangle(mesh, u, v, w, triangle_metadata, boundary_material_id);
}

void kinDS::SegmentBuilder::requireLiveRegisteredVoronoiEdgeEndpoints(size_t voronoi_edge_id, const char* context) const
{
  const size_t voronoi_he0 = 2 * voronoi_edge_id;
  const size_t voronoi_he1 = voronoi_he0 + 1;
  const int left_face = kin_del.getGraph().halfEdge(voronoi_he0).face;
  const int right_face = kin_del.getGraph().halfEdge(voronoi_he1).face;
  if (left_face >= 0)
  {
    kin_del.requireLiveRegisteredVoronoiVertex(static_cast<size_t>(left_face), context);
  }
  if (right_face >= 0)
  {
    kin_del.requireLiveRegisteredVoronoiVertex(static_cast<size_t>(right_face), context);
  }
}

void kinDS::SegmentBuilder::warnIfVoronoiVertexOutsideAlphaShape(
  const char* context, size_t voronoi_vertex_id, const glm::dvec3& position, size_t strand_id, double t) const
{
  const size_t containing_tri_id = kin_del.getCrossingDataContainingTriId(voronoi_vertex_id);
  if (kin_del.getFaceInside(containing_tri_id))
  {
    return;
  }
  // This spams too much, so we comment it out for now. It is useful for debugging, but not for normal operation.
  /*KINDS_WARNING("SegmentBuilder: " << context << " - Voronoi vertex " << voronoi_vertex_id
                                   << " (containing Delaunay triangle " << containing_tri_id
                                   << ") is outside the alpha-shape; position (" << position.x << ", " << position.y
                                   << ", " << position.z << "); " << formatStrandBranchLogInfo(strand_id, t));*/
}

std::string kinDS::SegmentBuilder::formatStrandBranchLogInfo(size_t strand_id, double t) const
{
  std::ostringstream oss;
  oss << "strand=" << strand_id;
  if (strand_id == static_cast<size_t>(-1) || kin_del.isDummyBoundary(strand_id))
  {
    return oss.str();
  }

  const size_t runtime_branch = kin_del.getRuntimeBranchIdForStrand(strand_id);
  oss << " runtime_branch=";
  if (runtime_branch == KineticDelaunay::RuntimeBranchData::no_branch)
  {
    oss << "none";
  }
  else
  {
    oss << runtime_branch;
  }

  const size_t height = kin_del.getStrandTree().getHeight();
  size_t section = 0;
  if (height > 0)
  {
    section = static_cast<size_t>(std::ceil(t));
    if (section >= height)
    {
      section = height - 1;
    }
  }
  oss << " input_branch=" << kin_del.getStrandTree().getBranchIndex(strand_id, section);

  if (strand_id < kin_del.component_data.component_map.size())
  {
    oss << " component=" << kin_del.component_data.component_map[strand_id];
  }
  return oss.str();
}

glm::dvec3 SegmentBuilder::transformFromInputBranchToObjectSpace(glm::dvec3 vertex, size_t strand_id, double t) const
{
  if (strand_id == static_cast<size_t>(-1))
  {
    throw std::runtime_error("transformFromInputBranchToObjectSpace: invalid strand id.");
  }

  return kin_del.getStrandTree().transformToObjectSpace(vertex, strand_id, t);
}

glm::dvec3 SegmentBuilder::getPointInMeshSpace(size_t strand_id, double t) const
{
  if (strand_id == static_cast<size_t>(-1))
  {
    throw std::runtime_error("getPointInMeshSpace: invalid strand id.");
  }

  // Never bake kinetic separation into mesh coordinates.
  if (!create_transformed_mesh)
  {
    // CLI --untransformed: local profile coords only (no reference-branch remap, no separation).
    const glm::dvec2 xy
      = kin_del.getPointAt(strand_id, t, /*apply_reference_transform=*/false, /*include_virtual_offset=*/false);
    return glm::dvec3(xy.x, xy.y, t);
  }

  // Object-space mesh: local profile → object space (also without separation).
  const glm::dvec2 local_xy = kin_del.getStrandTree().evaluate(strand_id, t);
  return transformFromInputBranchToObjectSpace(glm::dvec3(local_xy.x, local_xy.y, t), strand_id, t);
}

glm::dvec3 SegmentBuilder::getCrossingCoordsInMeshSpace(
  KineticDelaunay::CrossingData::EdgeIntersectionRef intersection, double t) const
{
  kinDS::ensureCrossingIntersectionParamUpToDate(const_cast<KineticDelaunay&>(kin_del), intersection, t);
  const auto& ir = *intersection;
  const size_t d_he0 = 2 * ir.delaunay_edge_id;
  const size_t d_he1 = d_he0 + 1;
  const auto& graph = kin_del.getGraph();
  if (d_he1 >= graph.halfEdgeSlotCount())
  {
    return glm::dvec3(std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN(), t);
  }

  const int a = graph.halfEdge(d_he0).origin;
  const int b = graph.halfEdge(d_he1).origin;
  if (a < 0 || b < 0)
  {
    return glm::dvec3(std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN(), t);
  }

  const double param = ir.delaunay_edge_param;
  if (!std::isfinite(param))
  {
    return glm::dvec3(std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN(), t);
  }

  const glm::dvec3 a_mesh = getPointInMeshSpace(static_cast<size_t>(a), t);
  const glm::dvec3 b_mesh = getPointInMeshSpace(static_cast<size_t>(b), t);
  const glm::dvec3 pos = a_mesh * (1.0 - param) + b_mesh * param;
  if (!std::isfinite(pos.x) || !std::isfinite(pos.y))
  {
    return glm::dvec3(std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN(), t);
  }
  return pos;
}

glm::dvec3 SegmentBuilder::computeMeshSiteVertexPosition(glm::dvec3 profile_vertex, size_t strand_id, double t) const
{
  (void)profile_vertex;
  return getPointInMeshSpace(strand_id, t);
}

SegmentBuilder::MeshIntersectionObjectSpaceResult SegmentBuilder::computeMeshIntersectionObjectSpace(
  const MeshletVertexRuntimeInfo& runtime_info, glm::dvec3 fallback_profile_vertex, size_t fallback_strand_id,
  double t) const
{
  (void)fallback_profile_vertex;
  (void)fallback_strand_id;
  MeshIntersectionObjectSpaceResult result;
  if (!runtime_info.position_intersection.has_value())
  {
    throw std::runtime_error(
      "computeMeshIntersectionObjectSpace: intersection placement has no position_intersection.");
  }
  const auto ref = runtime_info.position_intersection.value();
  result.position = getCrossingCoordsInMeshSpace(ref, t);
  if (!vertexPositionFinite(result.position))
  {
    return result;
  }

  const size_t delaunay_edge_id = ref->delaunay_edge_id;
  const auto& graph = kin_del.getGraph();
  const size_t d_he0 = 2 * delaunay_edge_id;
  const size_t d_he1 = d_he0 + 1;
  if (d_he1 < graph.halfEdgeSlotCount())
  {
    const int a = graph.halfEdge(d_he0).origin;
    const int b = graph.halfEdge(d_he1).origin;
    if (a >= 0 && b >= 0)
    {
      const double param = ref->delaunay_edge_param;
      const glm::dvec3 a_mesh = getPointInMeshSpace(static_cast<size_t>(a), t);
      const glm::dvec3 b_mesh = getPointInMeshSpace(static_cast<size_t>(b), t);
      result.mesh_interpolation = IntersectionInterpolationDebug { a_mesh, b_mesh, param };
      // Keep the actual placement tied to the same captured endpoints/parameter reported in metadata. This also avoids
      // a second mutable CrossingData parameter read disagreeing with the diagnostic values during parallel callbacks.
      result.position = a_mesh * (1.0 - param) + b_mesh * param;
    }
  }

  return result;
}

std::optional<KineticDelaunay::CrossingData::EdgeIntersectionRef>
SegmentBuilder::findIncidentIntersectionOnDelaunayEdge(size_t voronoi_vertex_id, size_t delaunay_edge_id) const
{
  const auto& graph = kin_del.getGraph();
  if (voronoi_vertex_id >= graph.faceSlotCount() || !graph.isLiveFace(voronoi_vertex_id))
  {
    return std::nullopt;
  }
  const auto& crossing_data = kin_del.getCrossingData();
  for (size_t he : graph.face(voronoi_vertex_id).half_edges)
  {
    const size_t voronoi_edge_id = he / 2;
    if (voronoi_edge_id >= crossing_data.voronoi_edge_intersections.size())
    {
      continue;
    }
    for (const auto& ref : crossing_data.voronoi_edge_intersections[voronoi_edge_id])
    {
      if (ref->delaunay_edge_id == delaunay_edge_id)
      {
        return ref;
      }
    }
  }
  return std::nullopt;
}

SegmentBuilder::MeshVoronoiVertexObjectSpaceResult SegmentBuilder::computeMeshVoronoiVertexObjectSpace(
  size_t voronoi_vertex_id, glm::dvec3 fallback_profile_vertex, size_t fallback_strand_id, double t) const
{
  MeshVoronoiVertexObjectSpaceResult result;
  result.position = computeMeshSiteVertexPosition(fallback_profile_vertex, fallback_strand_id, t);

  const auto& graph = kin_del.getGraph();
  if (voronoi_vertex_id >= graph.faceSlotCount() || !graph.isLiveFace(voronoi_vertex_id))
  {
    return result;
  }

  // Crossing event: reuse the single position buffered in beforeEvent for this Voronoi vertex.
  if (active_crossing_mesh_position_.has_value() && active_crossing_voronoi_vertex_id_.has_value()
    && active_crossing_voronoi_vertex_id_.value() == voronoi_vertex_id)
  {
    result.position = active_crossing_mesh_position_.value();
    result.from_crossing_event_buffer = true;
    return result;
  }

  // Containing Delaunay triangle (may differ from the dual triangle of this Voronoi vertex).
  const std::optional<size_t> containing_tri_opt = kin_del.getCrossingData().peekContainingTriId(voronoi_vertex_id);
  if (!containing_tri_opt.has_value())
  {
    KINDS_WARNING("computeMeshVoronoiVertexObjectSpace: Voronoi vertex "
      << voronoi_vertex_id << " has no registered containing triangle at t=" << t << "; falling back.");
    return result;
  }
  const size_t containing_tri_id = containing_tri_opt.value();
  result.containing_tri_id = containing_tri_id;
  if (containing_tri_id >= graph.faceSlotCount() || !graph.isLiveFace(containing_tri_id))
  {
    KINDS_WARNING("computeMeshVoronoiVertexObjectSpace: containing triangle "
      << containing_tri_id << " for Voronoi vertex " << voronoi_vertex_id << " is not live at t=" << t
      << "; falling back.");
    return result;
  }

  const std::array<int, 3> containing_vertices = graph.getTriangleVertexIndices(containing_tri_id);
  bool containing_tri_infinite = false;
  for (int vertex_index : containing_vertices)
  {
    if (vertex_index < 0)
    {
      containing_tri_infinite = true;
      break;
    }
  }
  if (!containing_tri_infinite)
  {
    for (size_t he_id : graph.face(containing_tri_id).half_edges)
    {
      if (graph.isInfinite(he_id))
      {
        containing_tri_infinite = true;
        break;
      }
    }
  }
  result.containing_tri_infinite = containing_tri_infinite;

  // Collect finite sites (and their slots in the containing-triangle vertex order).
  std::array<size_t, 3> finite_site_ids {};
  std::array<size_t, 3> finite_slots {};
  size_t finite_count = 0;
  for (size_t i = 0; i < containing_vertices.size(); ++i)
  {
    if (containing_vertices[i] < 0)
    {
      continue;
    }
    if (finite_count >= finite_site_ids.size())
    {
      break;
    }
    finite_slots[finite_count] = i;
    finite_site_ids[finite_count] = static_cast<size_t>(containing_vertices[i]);
    ++finite_count;
  }

  if (containing_tri_infinite)
  {
    // By construction the Voronoi vertex lies on the boundary edge of the two finite sites.
    if (finite_count != 2)
    {
      KINDS_WARNING("computeMeshVoronoiVertexObjectSpace: infinite containing triangle "
        << containing_tri_id << " for Voronoi vertex " << voronoi_vertex_id << " has " << finite_count
        << " finite sites at t=" << t << "; falling back.");
      return result;
    }

    const std::vector<size_t> finite_strands = { finite_site_ids[0], finite_site_ids[1] };
    const size_t shared_ref = kin_del.sharedReferenceBranchForEventTrigger(finite_strands, eventIntervalUpperBound(t));

    const size_t dual_he = graph.face(voronoi_vertex_id).half_edges[0];
    const glm::dvec3 vv_shifted = kin_del.computeVoronoiVertexClampedInfinityWithReferenceBranch(
      dual_he, t, shared_ref, /*apply_reference_transform=*/true, /*include_virtual_offset=*/true);
    if (!vertexPositionFinite(vv_shifted))
    {
      KINDS_WARNING("computeMeshVoronoiVertexObjectSpace: non-finite shifted Delaunay Voronoi vertex "
        << voronoi_vertex_id << " at t=" << t << " (infinite containing tri); falling back.");
      return result;
    }
    const glm::dvec2 vv_xy(vv_shifted.x, vv_shifted.y);

    const glm::dvec2 p0_shifted = kin_del.getPointInDelaunaySpace(finite_site_ids[0], t, shared_ref);
    const glm::dvec2 p1_shifted = kin_del.getPointInDelaunaySpace(finite_site_ids[1], t, shared_ref);
    const glm::dvec3 p0_mesh = getPointInMeshSpace(finite_site_ids[0], t);
    const glm::dvec3 p1_mesh = getPointInMeshSpace(finite_site_ids[1], t);
    result.site_mesh_positions[0] = { finite_site_ids[0], p0_mesh };
    result.site_mesh_positions[1] = { finite_site_ids[1], p1_mesh };
    result.site_count = 2;

    const glm::dvec2 edge = p1_shifted - p0_shifted;
    const double len2 = glm::dot(edge, edge);
    if (len2 < 1e-24)
    {
      KINDS_WARNING("computeMeshVoronoiVertexObjectSpace: degenerate finite edge on infinite containing triangle "
        << containing_tri_id << " for Voronoi vertex " << voronoi_vertex_id << " at t=" << t << "; falling back.");
      return result;
    }
    // Same 1D parameter as intersection mesh placement: project onto the finite edge in Delaunay space.
    const double param = glm::dot(vv_xy - p0_shifted, edge) / len2;

    std::array<double, 3> bary { 0.0, 0.0, 0.0 };
    bary[finite_slots[0]] = 1.0 - param;
    bary[finite_slots[1]] = param;
    result.barycentric = bary;

    const glm::dvec3 mesh_pos = p0_mesh * (1.0 - param) + p1_mesh * param;
    if (!vertexPositionFinite(mesh_pos))
    {
      KINDS_WARNING(
        "computeMeshVoronoiVertexObjectSpace: non-finite edge-interpolated mesh position for Voronoi vertex "
        << voronoi_vertex_id << " at t=" << t << "; falling back.");
      return result;
    }

    result.position = mesh_pos;
    return result;
  }

  if (finite_count != 3)
  {
    KINDS_WARNING("computeMeshVoronoiVertexObjectSpace: containing triangle "
      << containing_tri_id << " for Voronoi vertex " << voronoi_vertex_id << " has " << finite_count
      << " finite sites at t=" << t << "; falling back.");
    return result;
  }

  const std::array<size_t, 3> site_ids = { finite_site_ids[0], finite_site_ids[1], finite_site_ids[2] };

  // One shared Delaunay frame for barycentric (containing-triangle strands).
  const std::vector<size_t> containing_strands = { site_ids[0], site_ids[1], site_ids[2] };
  const size_t shared_ref
    = kin_del.sharedReferenceBranchForEventTrigger(containing_strands, eventIntervalUpperBound(t));

  // 1) Voronoi vertex in that frame WITH kinetic separation shifts.
  const size_t dual_he = graph.face(voronoi_vertex_id).half_edges[0];
  const glm::dvec3 vv_shifted = kin_del.computeVoronoiVertexClampedInfinityWithReferenceBranch(
    dual_he, t, shared_ref, /*apply_reference_transform=*/true, /*include_virtual_offset=*/true);
  if (!vertexPositionFinite(vv_shifted))
  {
    KINDS_WARNING("computeMeshVoronoiVertexObjectSpace: non-finite shifted Delaunay Voronoi vertex "
      << voronoi_vertex_id << " at t=" << t << "; falling back.");
    return result;
  }
  const glm::dvec2 vv_xy(vv_shifted.x, vv_shifted.y);

  // 2) Containing sites: shifted (barycentric) vs unshifted (mesh interpolation).
  std::array<glm::dvec2, 3> shifted_sites {};
  std::array<glm::dvec3, 3> mesh_sites {};
  for (size_t i = 0; i < site_ids.size(); ++i)
  {
    const size_t strand_id = site_ids[i];
    // WITH shift — same shared frame as the Voronoi vertex above.
    shifted_sites[i] = kin_del.getPointInDelaunaySpace(strand_id, t, shared_ref);

    // WITHOUT shift — same mesh-space site placement as ordinary sites.
    mesh_sites[i] = getPointInMeshSpace(strand_id, t);
    result.site_mesh_positions[i] = { strand_id, mesh_sites[i] };
    result.site_count = i + 1;
  }

  // 3) Barycentrics from shifted Delaunay coordinates only.
  const std::optional<std::array<double, 3>> bary
    = barycentricCoordinates2D(shifted_sites[0], shifted_sites[1], shifted_sites[2], vv_xy);
  if (!bary.has_value())
  {
    KINDS_WARNING("computeMeshVoronoiVertexObjectSpace: degenerate containing triangle "
      << containing_tri_id << " for Voronoi vertex " << voronoi_vertex_id << " at t=" << t << "; falling back.");
    return result;
  }
  result.barycentric = bary;

  // 4) Interpolate unshifted mesh site positions with those weights.
  const glm::dvec3 mesh_pos
    = bary.value()[0] * mesh_sites[0] + bary.value()[1] * mesh_sites[1] + bary.value()[2] * mesh_sites[2];
  if (!vertexPositionFinite(mesh_pos))
  {
    KINDS_WARNING("computeMeshVoronoiVertexObjectSpace: non-finite barycentric mesh position for Voronoi vertex "
      << voronoi_vertex_id << " at t=" << t << "; falling back.");
    return result;
  }

  result.position = mesh_pos;
  return result;
}

size_t kinDS::SegmentBuilder::addMeshletVertex(VoronoiMesh& mesh, const std::vector<BoundaryPoint>& boundary_polygon,
  const glm::dvec2& centroid, glm::dvec3 vertex, size_t strand_id, double t, bool includes_virtual_shift,
  std::optional<size_t> meshlet_voronoi_vertex_for_alpha_check, const std::string& metadata,
  const std::optional<glm::dvec3>& debug_color, const MeshletVertexRuntimeInfo& runtime_info)
{
  if (warnAndSkipIfMeshletCompleted(mesh, "addMeshletVertex", t))
  {
    return static_cast<size_t>(-1);
  }

  const auto warn_degenerate_or_non_finite = [&](const glm::dvec3& p, const char* stage)
  {
    if (!std::isfinite(p.x) || !std::isfinite(p.y) || !std::isfinite(p.z))
    {
      KINDS_WARNING("addMeshletVertex(" << stage << "): non-finite vertex at t=" << t << " -> (" << p.x << ", " << p.y
                                        << ", " << p.z << "); " << formatStrandBranchLogInfo(strand_id, t));
      return;
    }
    if (p.x == 0.0 && p.y == 0.0)
    {
      KINDS_WARNING("addMeshletVertex(" << stage << "): degenerate XY vertex (0,0,z) at t=" << t << " -> (" << p.x
                                        << ", " << p.y << ", " << p.z << "); "
                                        << formatStrandBranchLogInfo(strand_id, t));
    }
  };

  const bool is_flexible_placeholder = runtime_info.is_flexible_placeholder;
  const bool radius_shift_explicit_profile_position = runtime_info.radius_shift_explicit_profile_position;
  const bool radius_transition_projection = runtime_info.radius_transition_projection.has_value();
  if (!is_flexible_placeholder)
  {
    warn_degenerate_or_non_finite(vertex, "input");
  }
  glm::dvec2 profile_xy(vertex.x, vertex.y);
  // Never bake kinetic/mesh virtual separation offsets into stored mesh vertices.
  includes_virtual_shift = false;

  const bool is_intersection_vertex = runtime_info.isIntersectionVertex();
  const bool crossing_event_buffer = !is_flexible_placeholder
    && activeCrossingEventBufferApplies(meshlet_voronoi_vertex_for_alpha_check, runtime_info.position_intersection);
  glm::dvec2 delaunay_xy(std::numeric_limits<double>::quiet_NaN());
  if (!is_flexible_placeholder)
  {
    if (runtime_info.explicit_delaunay_xy.has_value())
    {
      delaunay_xy = runtime_info.explicit_delaunay_xy.value();
    }
    else if (crossing_event_buffer && active_crossing_delaunay_xy_.has_value())
    {
      delaunay_xy = active_crossing_delaunay_xy_.value();
    }
    else if (meshlet_voronoi_vertex_for_alpha_check.has_value())
    {
      const size_t voronoi_vertex_id = meshlet_voronoi_vertex_for_alpha_check.value();
      const auto& graph = kin_del.getGraph();
      if (voronoi_vertex_id < graph.faceSlotCount() && graph.isLiveFace(voronoi_vertex_id))
      {
        delaunay_xy = glm::dvec2(computeVoronoiVertex(graph.face(voronoi_vertex_id).half_edges[0], t));
      }
    }
    else if (radius_transition_projection)
    {
      delaunay_xy
        = glm::dvec2(radiusTransitionProjectionPosition(runtime_info.radius_transition_projection.value(), t, false));
    }
    else if (is_intersection_vertex)
    {
      delaunay_xy
        = glm::dvec2(getCrossingCoordsInDelaunaySpace(kin_del, runtime_info.position_intersection.value(), t));
    }
    else
    {
      delaunay_xy = kin_del.getPointInDelaunaySpace(strand_id, t);
    }
  }

  std::optional<MeshIntersectionObjectSpaceResult> intersection_object_space;
  if (!is_flexible_placeholder && is_intersection_vertex && !radius_transition_projection
    && !meshlet_voronoi_vertex_for_alpha_check.has_value() && !crossing_event_buffer)
  {
    intersection_object_space
      = computeMeshIntersectionObjectSpace(runtime_info, glm::dvec3(profile_xy.x, profile_xy.y, t), strand_id, t);
  }

  std::optional<MeshVoronoiVertexObjectSpaceResult> voronoi_object_space;
  if (!is_flexible_placeholder && meshlet_voronoi_vertex_for_alpha_check.has_value())
  {
    if (runtime_info.explicit_mesh_position.has_value())
    {
      // Flip (and similar): reuse mesh-space coords buffered before topology/strand reassignment.
      vertex = runtime_info.explicit_mesh_position.value();
    }
    else
    {
      // Prefer barycentric Voronoi placement; the event VV reuses the crossing buffer.
      voronoi_object_space = computeMeshVoronoiVertexObjectSpace(
        meshlet_voronoi_vertex_for_alpha_check.value(), glm::dvec3(profile_xy.x, profile_xy.y, t), strand_id, t);
      vertex = voronoi_object_space.value().position;
      if (voronoi_object_space.value().from_crossing_event_buffer && active_crossing_delaunay_xy_.has_value())
      {
        delaunay_xy = active_crossing_delaunay_xy_.value();
      }
    }
    if (!vertexPositionFinite(vertex))
    {
      KINDS_WARNING("addMeshletVertex: non-finite mesh-space Voronoi vertex at t="
        << t << "; falling back to site-equivalent profile placement; " << formatStrandBranchLogInfo(strand_id, t));
      vertex = computeMeshSiteVertexPosition(glm::dvec3(profile_xy.x, profile_xy.y, t), strand_id, t);
    }
  }
  else if (!is_flexible_placeholder && radius_transition_projection)
  {
    if (runtime_info.explicit_mesh_position.has_value())
    {
      vertex = runtime_info.explicit_mesh_position.value();
    }
    else
    {
      vertex = radiusTransitionProjectionPosition(runtime_info.radius_transition_projection.value(), t, true);
    }
    if (!vertexPositionFinite(vertex))
    {
      KINDS_WARNING("addMeshletVertex: non-finite radius-transition projection at t="
        << t << "; falling back to site placement; " << formatStrandBranchLogInfo(strand_id, t));
      vertex = computeMeshSiteVertexPosition(glm::dvec3(profile_xy.x, profile_xy.y, t), strand_id, t);
    }
  }
  else if (!is_flexible_placeholder && is_intersection_vertex && crossing_event_buffer
    && active_crossing_mesh_position_.has_value())
  {
    vertex = active_crossing_mesh_position_.value();
  }
  else if (!is_flexible_placeholder && is_intersection_vertex && intersection_object_space.has_value())
  {
    // Every ordinary/remapped intersection is reconstructed from its Delaunay edge parameter by
    // getCrossingCoordsInMeshSpace(), which obtains both endpoints through getPointInMeshSpace().
    vertex = intersection_object_space.value().position;
    if (!vertexPositionFinite(vertex))
    {
      const auto position_ref = runtime_info.position_intersection.value();
      std::ostringstream oss;
      oss << "addMeshletVertex: getCrossingCoordsInMeshSpace returned a non-finite intersection position; refusing "
             "site fallback (strand="
          << strand_id << ", t=" << t << ", input=(" << profile_xy.x << "," << profile_xy.y
          << "), position_de=" << position_ref->delaunay_edge_id << ", position_ve=" << position_ref->voronoi_edge_id
          << ", position_d_param=" << position_ref->delaunay_edge_param;
      if (runtime_info.conceptual_intersection.has_value())
      {
        const auto conceptual_ref = runtime_info.conceptual_intersection.value();
        oss << ", conceptual_de=" << conceptual_ref->delaunay_edge_id
            << ", conceptual_ve=" << conceptual_ref->voronoi_edge_id
            << ", conceptual_d_param=" << conceptual_ref->delaunay_edge_param;
      }
      oss << ", computed=(" << vertex.x << "," << vertex.y << "," << vertex.z << ")).";
      throw std::runtime_error(oss.str());
    }
  }
  else if (!is_flexible_placeholder && !is_intersection_vertex)
  {
    // Sites use their unshifted mesh placement. Parameterized radius-transition projections were handled above.
    // Keep profile_xy / delaunay_xy as Delaunay-plane coordinates for ear-clip triangulation; only @c vertex
    // becomes mesh space.
    vertex = computeMeshSiteVertexPosition(glm::dvec3(profile_xy.x, profile_xy.y, t), strand_id, t);
    includes_virtual_shift = false;
  }
  else if (create_transformed_mesh && !is_flexible_placeholder)
  {
    vertex = transformFromInputBranchToObjectSpace(vertex, strand_id, t);
    if (!vertexPositionFinite(vertex))
    {
      KINDS_WARNING("addMeshletVertex: non-finite transformed vertex at t="
        << t << "; falling back to input-branch transform; " << formatStrandBranchLogInfo(strand_id, t));
      vertex = transformFromInputBranchToObjectSpace(glm::dvec3(profile_xy.x, profile_xy.y, t), strand_id, t);
    }
  }
  if (!is_flexible_placeholder)
  {
    warn_degenerate_or_non_finite(vertex, "stored");
  }
  if (meshlet_voronoi_vertex_for_alpha_check.has_value())
  {
    warnIfVoronoiVertexOutsideAlphaShape("addMeshletVertex", meshlet_voronoi_vertex_for_alpha_check.value(),
      glm::dvec3(profile_xy.x, profile_xy.y, t), strand_id, t);
  }
  std::string vertex_metadata = metadata;
  if (store_mesh_metadata)
  {
    const std::optional<std::string> event_type_field = metadataStringField(metadata, "event_type");
    const std::string event_type = event_type_field.value_or("unknown_event");
    if (!event_type_field.has_value())
    {
      KINDS_WARNING("addMeshletVertex: missing metadata event_type at t=" << t << "; using 'unknown_event'; "
                                                                          << formatStrandBranchLogInfo(strand_id, t));
    }
    else if (event_type == "unknown_event")
    {
      KINDS_WARNING(
        "addMeshletVertex: unknown metadata event_type at t=" << t << "; " << formatStrandBranchLogInfo(strand_id, t));
    }
    // Source follows how the vertex was placed. VV snap wins over intersection/site labels.
    const std::string source = meshlet_voronoi_vertex_for_alpha_check.has_value()
      ? "Voronoi vertex"
      : (is_intersection_vertex ? "intersection" : "site");
    const std::optional<std::string> source_field = metadataStringField(metadata, "source");
    if (source_field.has_value() && source_field.value() != source)
    {
      KINDS_WARNING("addMeshletVertex: caller metadata source '"
        << source_field.value() << "' disagrees with placement source '" << source << "' at t=" << t
        << "; using placement source; " << formatStrandBranchLogInfo(strand_id, t));
    }
    else if (!source_field.has_value())
    {
      KINDS_WARNING("addMeshletVertex: missing metadata source at t="
        << t << "; using placement source '" << source << "'; " << formatStrandBranchLogInfo(strand_id, t));
    }
    MetadataBuilder builder;
    builder.addString("event_type", event_type).addString("source", source);

    if (source == "Voronoi vertex")
    {
      if (meshlet_voronoi_vertex_for_alpha_check.has_value())
      {
        builder.addSize("voronoi_vertex_id", meshlet_voronoi_vertex_for_alpha_check.value());
      }
      else if (const auto voronoi_vertex_id = metadataSizeField(metadata, "voronoi_vertex_id");
        voronoi_vertex_id.has_value())
      {
        builder.addSize("voronoi_vertex_id", voronoi_vertex_id.value());
      }
      if (const auto pos = metadataStringField(metadata, "pos"); pos.has_value())
      {
        builder.addString("pos", pos.value());
      }
      if (voronoi_object_space.has_value())
      {
        if (voronoi_object_space.value().containing_tri_id.has_value())
        {
          builder.addSize("containing_tri_id", voronoi_object_space.value().containing_tri_id.value());
        }
        if (voronoi_object_space.value().containing_tri_infinite.has_value())
        {
          builder.addBool("containing_tri_infinite", voronoi_object_space.value().containing_tri_infinite.value());
        }
        if (voronoi_object_space.value().barycentric.has_value())
        {
          const auto& b = voronoi_object_space.value().barycentric.value();
          builder.addDouble("bary_0", b[0]).addDouble("bary_1", b[1]).addDouble("bary_2", b[2]);
        }
        for (size_t i = 0; i < voronoi_object_space.value().site_count; ++i)
        {
          const auto& [site_id, site_pos] = voronoi_object_space.value().site_mesh_positions[i];
          const std::string key = "p_" + std::to_string(site_id);
          std::ostringstream value;
          value.setf(std::ios::fixed);
          value.precision(17);
          value << "(" << site_pos.x << "," << site_pos.y << "," << site_pos.z << ")";
          builder.addString(key.c_str(), value.str());
        }
      }
    }
    else if (is_intersection_vertex)
    {
      const auto position_ref = runtime_info.position_intersection.value();
      builder.addSize("delaunay_edge_id", position_ref->delaunay_edge_id)
        .addSize("voronoi_edge_id", position_ref->voronoi_edge_id)
        .addDouble("stored_d_param", position_ref->delaunay_edge_param);
      if (runtime_info.conceptual_intersection.has_value()
        && runtime_info.conceptual_intersection.value() != position_ref)
      {
        builder.addSize("conceptual_delaunay_edge_id", runtime_info.conceptual_intersection.value()->delaunay_edge_id)
          .addSize("conceptual_voronoi_edge_id", runtime_info.conceptual_intersection.value()->voronoi_edge_id);
      }
      if (const auto pos = metadataStringField(metadata, "pos"); pos.has_value())
      {
        builder.addString("pos", pos.value());
      }
      if (radius_shift_explicit_profile_position)
      {
        builder.addBool("radius_shift_explicit_profile_position", true);
      }
      if (intersection_object_space.has_value() && intersection_object_space.value().mesh_interpolation.has_value())
      {
        appendIntersectionInterpolationDebugToMetadata(
          builder, intersection_object_space.value().mesh_interpolation.value());
        builder.addDouble("interp_result_x", vertex.x)
          .addDouble("interp_result_y", vertex.y)
          .addDouble("interp_result_z", vertex.z);
      }
    }
    else
    {
      builder.addSize("strand_id", strand_id);
    }

    if (crossing_event_buffer)
    {
      builder.addBool("crossing_event_buffered_vertex", true);
    }

    if (!metadata_callback_phase_.empty())
    {
      builder.addString("callback", metadata_callback_phase_);
    }

    const glm::dvec2 metadata_delaunay_xy
      = (std::isfinite(delaunay_xy.x) && std::isfinite(delaunay_xy.y)) ? delaunay_xy : profile_xy;
    vertex_metadata = builder.addBool("shift", includes_virtual_shift)
                        .addDouble("x", metadata_delaunay_xy.x)
                        .addDouble("y", metadata_delaunay_xy.y)
                        .addDouble("mesh_x", vertex.x)
                        .addDouble("mesh_y", vertex.y)
                        .addDouble("mesh_z", vertex.z)
                        .addDouble("t", t)
                        .build();
  }
  size_t index
    = mesh.addVertex(vertex, vertex_metadata, debug_color.has_value() ? debug_color.value() : glm::dvec3(1.0));
  if (is_flexible_placeholder)
  {
    mesh.setVertexFlexible(index, true);
  }
  if (is_intersection_vertex && mesh.getVertices()[index] != vertex)
  {
    const glm::dvec3& stored = mesh.getVertices()[index];
    std::ostringstream oss;
    oss << "addMeshletVertex: intersection coordinates changed during mesh insertion (index=" << index << ", expected=("
        << vertex.x << "," << vertex.y << "," << vertex.z << "), stored=(" << stored.x << "," << stored.y << ","
        << stored.z << ")).";
    throw std::runtime_error(oss.str());
  }
  mesh.setVertexKineticTime(index, t);
  if (!is_flexible_placeholder)
  {
    // Triangulation plane must match the polygon ring's construction space (usually Delaunay). Prefer an
    // explicit caller-provided XY; otherwise fall back to recomputed Delaunay coords. Never use mesh/object XY.
    glm::dvec2 plane_xy = runtime_info.triangulation_plane_xy.value_or(delaunay_xy);
    if (!std::isfinite(plane_xy.x) || !std::isfinite(plane_xy.y))
    {
      plane_xy = delaunay_xy;
    }
    if (!std::isfinite(plane_xy.x) || !std::isfinite(plane_xy.y))
    {
      plane_xy = kin_del.getPointInDelaunaySpace(strand_id, t);
    }
    mesh.setProfilePlanePosition(index, plane_xy);
  }
  const bool is_boundary_interval_meshlet = intersectionMeshletIndexForMesh(mesh).has_value();
  if (is_boundary_interval_meshlet)
  {
    std::vector<glm::dvec2>& raw_uvs = boundaryIntervalRawUvs(mesh);
    if (raw_uvs.size() <= index)
    {
      raw_uvs.resize(index + 1, glm::dvec2(0.0));
    }
    if (is_flexible_placeholder)
    {
      // Placeholder raw UV is filled when the flex vertex is interpolated along its anchor edge.
      raw_uvs[index] = glm::dvec2(0.0, t);
    }
    else
    {
      if (!std::isfinite(delaunay_xy.x) || !std::isfinite(delaunay_xy.y))
      {
        delaunay_xy = kin_del.getPointInDelaunaySpace(strand_id, t);
      }
      raw_uvs[index] = boundaryRawUv(delaunay_xy, centroid, t);
    }
  }
  else if (is_flexible_placeholder)
  {
    // Placeholder UV is filled when the flex vertex is interpolated along its anchor edge.
    mesh.setVertexSemanticUv(index, glm::dvec3(0.0, 0.0, t * uv_height_factor));
  }
  else
  {
    if (!std::isfinite(delaunay_xy.x) || !std::isfinite(delaunay_xy.y))
    {
      delaunay_xy = kin_del.getPointInDelaunaySpace(strand_id, t);
    }
    mesh.setVertexSemanticUv(index, interiorMeshUv(boundary_polygon, centroid, delaunay_xy, t));
  }
  return index;
}

size_t SegmentBuilder::intersectionStripEffectiveVertexIndex(const MeshingData& seg, bool left_side) const
{
  if (left_side)
  {
    if (!seg.flexible_left_vertex_ids.empty())
    {
      const int v = seg.flexible_left_vertex_ids.back();
      if (v >= 0)
      {
        return static_cast<size_t>(v);
      }
    }
    if (seg.mesh_start_vertex_id < 0)
    {
      throw std::runtime_error("intersectionStripEffectiveVertexIndex: invalid mesh_start_vertex_id.");
    }
    return static_cast<size_t>(seg.mesh_start_vertex_id);
  }
  if (!seg.flexible_right_vertex_ids.empty())
  {
    const int v = seg.flexible_right_vertex_ids.back();
    if (v >= 0)
    {
      return static_cast<size_t>(v);
    }
  }
  if (seg.mesh_end_vertex_id < 0)
  {
    throw std::runtime_error("intersectionStripEffectiveVertexIndex: invalid mesh_end_vertex_id.");
  }
  return static_cast<size_t>(seg.mesh_end_vertex_id);
}

void SegmentBuilder::addFlexibleVertexToIntersectionMesh(VoronoiMesh& mesh, MeshingData& seg,
  bool flexible_on_left_side, const std::vector<BoundaryPoint>& boundary_polygon, const glm::dvec2& centroid,
  size_t strand_id, double t, const std::string& metadata)
{
  if (!intersection_strip_flexible_vertices_enabled)
  {
    return;
  }
  if (seg.mesh_start_vertex_id < 0 || seg.mesh_end_vertex_id < 0)
  {
    return;
  }
  // Wedge base must match crossing/finish: latest flex on a side wins over fixed mesh_start/mesh_end. Snapshot before
  // appending the new placeholder so the base is (eff_l, eff_r), not always (mesh_start, mesh_end).
  const size_t eff_l = intersectionStripEffectiveVertexIndex(seg, true);
  const size_t eff_r = intersectionStripEffectiveVertexIndex(seg, false);
  if (eff_l == eff_r)
  {
    return;
  }

  const glm::dvec3 placeholder(std::nan(""), std::nan(""), t);
  const std::string vertex_meta = [&]()
  {
    const std::string parent_event_type = metadataStringField(metadata, "event_type").value_or("unknown_event");
    MetadataBuilder builder;
    builder.addString("event_type", "intersection_flexible_placeholder")
      .addString("source", "site")
      .addString("parent_event_type", parent_event_type)
      .addString("pos", flexible_on_left_side ? "left" : "right")
      .addBool("intersection_flexible_placeholder", true)
      .addSize("strand_id", strand_id);
    return builder.build();
  }();
  const size_t idx = addMeshletVertex(mesh, boundary_polygon, centroid, placeholder, strand_id, t, false, std::nullopt,
    vertex_meta, std::nullopt, MeshletVertexRuntimeInfo { true, false, std::nullopt, std::nullopt });
  if (flexible_on_left_side)
  {
    seg.flexible_left_vertex_ids.push_back(static_cast<int>(idx));
  }
  else
  {
    seg.flexible_right_vertex_ids.push_back(static_cast<int>(idx));
  }

  const int inside_he = seg.start_half_edge_id >= 0 ? seg.start_half_edge_id : seg.end_half_edge_id;
  if (inside_he < 0)
  {
    return;
  }
  addBoundaryIntervalTriangleOriented(mesh, eff_l, eff_r, idx, inside_he, t, vertex_meta);
}

void SegmentBuilder::extendIntersectionMeshAtSharedCrossing(size_t neighbor_pair_idx,
  KineticDelaunay::CrossingData::EdgeIntersectionRef shared_ref, bool update_start_on_neighbor, double t,
  BoundaryEventType event_type, BoundarySegmentAction segment_action,
  const RadiusBoundaryTransitionShiftContext* boundary_transition_shift, bool append_flexible_placeholder,
  std::optional<size_t> skip_pair_idx)
{
  if (neighbor_pair_idx == static_cast<size_t>(-1)
    || (skip_pair_idx.has_value() && neighbor_pair_idx == skip_pair_idx.value()))
  {
    return;
  }
  neighbor_pair_idx = preferLiveBoundaryMeshPair(neighbor_pair_idx);
  if (neighbor_pair_idx == static_cast<size_t>(-1)
    || (skip_pair_idx.has_value() && neighbor_pair_idx == skip_pair_idx.value()))
  {
    return;
  }
  if (neighbor_pair_idx >= intersection_meshes.size()
    || neighbor_pair_idx >= intersection_mesh_pair_last_left_and_right_vertex.size())
  {
    return;
  }

  const size_t shared_d_edge = shared_ref->delaunay_edge_id;
  const size_t shared_he_even = 2 * shared_d_edge;
  if (shared_he_even + 1 >= kin_del.getGraph().halfEdgeSlotCount() || !kin_del.isOnComponentBoundary(shared_he_even))
  {
    return;
  }

  auto& segs = intersection_mesh_pair_last_left_and_right_vertex[neighbor_pair_idx];
  if (segs.empty())
  {
    return;
  }
  auto& seg = segs.front();
  if (seg.mesh_start_vertex_id < 0 || seg.mesh_end_vertex_id < 0)
  {
    return;
  }

  size_t neighbor_cell = static_cast<size_t>(-1);
  if (neighbor_pair_idx < intersection_mesh_pair_metadata.size())
  {
    const size_t cid = intersection_mesh_pair_metadata[neighbor_pair_idx].voronoi_cell_id;
    if (cid != static_cast<size_t>(-1))
    {
      neighbor_cell = cid;
    }
  }
  if (neighbor_cell == static_cast<size_t>(-1))
  {
    neighbor_cell = determineVoronoiCellForBoundaryIntersectionInterval(shared_ref->delaunay_edge_id,
      update_start_on_neighbor ? std::make_optional(shared_ref) : std::nullopt,
      update_start_on_neighbor ? std::nullopt : std::make_optional(shared_ref));
  }

  const size_t d_edge_id = shared_ref->delaunay_edge_id;
  const size_t he_even = 2 * d_edge_id;
  auto& graph = kin_del.getGraph();
  int inside_boundary_he_id = -1;
  if (he_even + 1 < graph.halfEdgeSlotCount() && kin_del.isOnComponentBoundary(he_even))
  {
    const bool boundary_even_out = kin_del.isOnComponentBoundaryOutside(he_even);
    inside_boundary_he_id = static_cast<int>(boundary_even_out ? he_even + 1 : he_even);
  }
  if (inside_boundary_he_id < 0)
  {
    inside_boundary_he_id = update_start_on_neighbor
      ? (seg.start_half_edge_id >= 0 ? seg.start_half_edge_id : seg.end_half_edge_id)
      : (seg.end_half_edge_id >= 0 ? seg.end_half_edge_id : seg.start_half_edge_id);
  }
  if (inside_boundary_he_id < 0)
  {
    return;
  }

  const size_t neighbor_component = kin_del.component_data.component_map[neighbor_cell];
  std::vector<bool> he_vis(graph.halfEdgeSlotCount(), false);
  updateBoundary(t, he_vis, neighbor_component);
  auto& neighbor_boundary = kin_del.component_data.component_boundaries[neighbor_component][0];
  const auto neighbor_centroid = polygonCentroid(neighbor_boundary);

  auto& mesh = intersection_meshes[neighbor_pair_idx];
  const RadiusBoundaryTransitionCrossingPlacement placement
    = resolveRadiusBoundaryTransitionCrossingPlacement(t, shared_ref, boundary_transition_shift);
  const glm::dvec3 crossing_pos = crossingProfilePositionFromPlacement(t, placement);
  if (placement.positionDiffersFromConceptual())
  {
    const glm::dvec3 old_pos = crossingProfilePosition(t, placement.conceptual_intersection);
    logRadiusBoundaryTransitionVertexShift("extendIntersectionMeshAtSharedCrossing", t,
      placement.conceptual_intersection, placement.position_intersection, old_pos, crossing_pos);
  }

  const std::string base_meta = composeBoundaryMetadata(event_type, segment_action);
  const std::string vertex_meta = store_mesh_metadata
    ? intersectionCrossingVertexMetadata(base_meta, placement.conceptual_intersection, placement.position_intersection,
        update_start_on_neighbor ? "left" : "right", placement.explicit_profile_position.has_value())
    : std::string {};
  const glm::dvec3 vertex_color = update_start_on_neighbor ? glm::dvec3(1.0, 0.0, 0.0) : glm::dvec3(0.0, 0.0, 1.0);

  const size_t eff_l = intersectionStripEffectiveVertexIndex(seg, true);
  const size_t eff_r = intersectionStripEffectiveVertexIndex(seg, false);
  const size_t new_vid = addMeshletVertex(mesh, neighbor_boundary, neighbor_centroid, crossing_pos, neighbor_cell, t,
    false, placement.snap_voronoi_vertex_id, vertex_meta, vertex_color,
    MeshletVertexRuntimeInfo { false, placement.explicit_profile_position.has_value(), placement.position_intersection,
      placement.conceptual_intersection, placement.projection });
  addBoundaryIntervalTriangleOriented(
    mesh, eff_l, eff_r, new_vid, inside_boundary_he_id, t, base_meta, neighbor_pair_idx);
  applyIntersectionStripOneSidedFixedVertex(mesh, seg, update_start_on_neighbor, new_vid, inside_boundary_he_id,
    std::make_optional(shared_ref), neighbor_boundary, neighbor_centroid, neighbor_cell, t, true,
    append_flexible_placeholder, base_meta);
  maybeLogDiagnosticsMonitoredDelaunayEdgeTrigger(
    t, "extendIntersectionMeshAtSharedCrossing", shared_ref->delaunay_edge_id, neighbor_pair_idx);
}

void SegmentBuilder::applyIntersectionStripOneSidedFixedVertex(VoronoiMesh& mesh, MeshingData& seg,
  bool fixed_start_side, size_t new_fixed_vertex_index, int inside_half_edge_id,
  const std::optional<KineticDelaunay::CrossingData::EdgeIntersectionRef>& new_crossing_for_updated_side,
  const std::vector<BoundaryPoint>& boundary_polygon, const glm::dvec2& centroid, size_t strand_id, double t,
  bool keep_strip_alive, bool append_flexible_placeholder, const std::string& flexible_base_metadata)
{
  const int old_fixed = fixed_start_side ? seg.mesh_start_vertex_id : seg.mesh_end_vertex_id;
  if (old_fixed < 0)
  {
    return;
  }
  std::vector<int>& flex = fixed_start_side ? seg.flexible_left_vertex_ids : seg.flexible_right_vertex_ids;
  if (interpolateFlexibleVerticesAlongEdge(mesh, flex, static_cast<size_t>(old_fixed), new_fixed_vertex_index))
  {
    flex.clear();
  }
  if (!keep_strip_alive)
  {
    return;
  }
  if (fixed_start_side)
  {
    seg.mesh_start_vertex_id = static_cast<int>(new_fixed_vertex_index);
    seg.start_half_edge_id = inside_half_edge_id;
    seg.start_crossing = new_crossing_for_updated_side;
  }
  else
  {
    seg.mesh_end_vertex_id = static_cast<int>(new_fixed_vertex_index);
    seg.end_half_edge_id = inside_half_edge_id;
    seg.end_crossing = new_crossing_for_updated_side;
  }
  if (append_flexible_placeholder)
  {
    addFlexibleVertexToIntersectionMesh(
      mesh, seg, !fixed_start_side, boundary_polygon, centroid, strand_id, t, flexible_base_metadata);
  }
}

void SegmentBuilder::applyIntersectionStripUniformClosureVertex(
  VoronoiMesh& mesh, MeshingData& seg, size_t closure_vertex_index)
{
  if (seg.mesh_start_vertex_id < 0 || seg.mesh_end_vertex_id < 0)
  {
    return;
  }
  if (interpolateFlexibleVerticesAlongEdge(
        mesh, seg.flexible_left_vertex_ids, static_cast<size_t>(seg.mesh_start_vertex_id), closure_vertex_index))
  {
    seg.flexible_left_vertex_ids.clear();
  }
  if (interpolateFlexibleVerticesAlongEdge(
        mesh, seg.flexible_right_vertex_ids, static_cast<size_t>(seg.mesh_end_vertex_id), closure_vertex_index))
  {
    seg.flexible_right_vertex_ids.clear();
  }
}

void SegmentBuilder::resolveRemainingFlexibleVertices(VoronoiMesh& mesh, MeshingData& seg, const char* context)
{
  if (seg.mesh_start_vertex_id < 0 || seg.mesh_end_vertex_id < 0)
  {
    if (!seg.flexible_left_vertex_ids.empty() || !seg.flexible_right_vertex_ids.empty())
    {
      KINDS_WARNING("Unresolved flexible vertices in "
        << (context != nullptr ? context : "unknown") << " without both fixed anchors: left="
        << seg.flexible_left_vertex_ids.size() << " right=" << seg.flexible_right_vertex_ids.size() << ".");
    }
    return;
  }

  const size_t start_anchor = static_cast<size_t>(seg.mesh_start_vertex_id);
  const size_t end_anchor = static_cast<size_t>(seg.mesh_end_vertex_id);

  if (!seg.flexible_left_vertex_ids.empty())
  {
    KINDS_WARNING("Resolving leftover left flexible vertices in "
      << (context != nullptr ? context : "unknown") << ": count=" << seg.flexible_left_vertex_ids.size() << ".");
    if (!interpolateFlexibleVerticesAlongEdge(mesh, seg.flexible_left_vertex_ids, start_anchor, end_anchor))
    {
      snapFlexibleVerticesToAnchor(mesh, seg.flexible_left_vertex_ids, start_anchor);
    }
    seg.flexible_left_vertex_ids.clear();
  }
  if (!seg.flexible_right_vertex_ids.empty())
  {
    KINDS_WARNING("Resolving leftover right flexible vertices in "
      << (context != nullptr ? context : "unknown") << ": count=" << seg.flexible_right_vertex_ids.size() << ".");
    if (!interpolateFlexibleVerticesAlongEdge(mesh, seg.flexible_right_vertex_ids, end_anchor, start_anchor))
    {
      snapFlexibleVerticesToAnchor(mesh, seg.flexible_right_vertex_ids, end_anchor);
    }
    seg.flexible_right_vertex_ids.clear();
  }

  for (int flex_vertex_id : { seg.mesh_start_vertex_id, seg.mesh_end_vertex_id })
  {
    if (flex_vertex_id < 0)
    {
      continue;
    }
    const size_t flex_index = static_cast<size_t>(flex_vertex_id);
    if (flex_index < mesh.getVertices().size() && !vertexPositionFinite(mesh.getVertices()[flex_index]))
    {
      KINDS_WARNING("resolveRemainingFlexibleVertices(" << (context != nullptr ? context : "unknown")
                                                        << "): anchor vertex " << flex_index << " remains non-finite.");
    }
  }
}

void SegmentBuilder::resolveAllIntersectionFlexibleVertices(const char* context)
{
  for (size_t mesh_id = 0;
    mesh_id < intersection_meshes.size() && mesh_id < intersection_mesh_pair_last_left_and_right_vertex.size();
    ++mesh_id)
  {
    for (MeshingData& seg : intersection_mesh_pair_last_left_and_right_vertex[mesh_id])
    {
      resolveRemainingFlexibleVertices(intersection_meshes[mesh_id], seg, context);
    }
  }

  for (size_t mesh_id = 0; mesh_id < intersection_meshes.size(); ++mesh_id)
  {
    const VoronoiMesh& mesh = intersection_meshes[mesh_id];
    const auto& vertices = mesh.getVertices();
    for (size_t vertex_index = 0; vertex_index < vertices.size(); ++vertex_index)
    {
      if (vertexPositionFinite(vertices[vertex_index]))
      {
        continue;
      }
      const std::string event_type = (vertex_index < mesh.getVertexMetadata().size())
        ? metadataStringField(mesh.getVertexMetadata()[vertex_index], "event_type").value_or("unknown")
        : std::string("unknown");
      KINDS_WARNING("Non-finite intersection mesh vertex " << vertex_index << " on mesh " << mesh_id << " in "
                                                           << (context != nullptr ? context : "unknown")
                                                           << " (event_type=" << event_type << ").");
    }
  }
}

void SegmentBuilder::collapseDegreeTwoFlexibleVerticesOnIntersectionMeshes()
{
  if (!collapse_degree_two_flexible_vertices_postprocess_enabled)
  {
    return;
  }
  for (VoronoiMesh& mesh : intersection_meshes)
  {
    mesh.collapseDegreeTwoFlexibleVertices();
  }
}

void kinDS::SegmentBuilder::addDelaunayTriangulationToBoundaryMesh(
  double t, size_t input_branch_id, bool invert_orientation, double offset)
{
  auto& graph = kin_del.getGraph();
  const size_t section = static_cast<size_t>(t);
  const auto& branch_strands = kin_del.getStrandTree().getStrandsByBranch(section, input_branch_id);
  const std::unordered_set<size_t> branch_strand_set(branch_strands.begin(), branch_strands.end());

  size_t skip_dummy = 0;
  size_t strands_vertex_added = 0;
  size_t skip_infinite_face = 0;
  size_t skip_not_in_branch = 0;
  size_t skip_no_mesh_vertex = 0;
  size_t skip_outside_face = 0;
  size_t triangles_added = 0;
  const size_t boundary_tris_before = diagnostics ? boundary_mesh.getTriangleCount() : 0;

  std::vector<int> strand_to_mesh_vertex(graph.getVertexCount(), -1);

  for (size_t strand_id : branch_strands)
  {
    if (kin_del.isDummyBoundary(strand_id))
    {
      ++skip_dummy;
      if (diagnostics && strand_id == 0)
      {
        strandInitDiagnosticLogLine("boundary_mesh_strand_skip", strand_id, t,
          "reason=is_dummy_boundary during addDelaunayTriangulationToBoundaryMesh");
      }
      continue;
    }

    glm::dvec2 vertex = kin_del.getPointAt(t, strand_id, false, false);

    auto component_index = kin_del.component_data.component_map[strand_id];
    auto& centroid = kin_del.component_data.component_centroids[component_index];

    const size_t mesh_vertex_index
      = addBoundaryVertex(glm::dvec3 { vertex[0], vertex[1], t + offset }, centroid, strand_id, t, false);

    strand_to_mesh_vertex[strand_id] = static_cast<int>(mesh_vertex_index);
    ++strands_vertex_added;

    if (diagnostics && strand_id == 0)
    {
      std::ostringstream oss;
      oss << "input_branch_id=" << input_branch_id << " boundary_mesh_vertex=" << mesh_vertex_index
          << " offset=" << offset << " invert_orientation=" << (invert_orientation ? "true" : "false");
      strandInitDiagnosticLogLine("boundary_mesh_strand_vertex", strand_id, t, oss.str().c_str());
    }
  }

  for (size_t face_index : graph.liveFaces())
  {
    const auto& triangle = graph.face(face_index);
    const auto& he_ids = triangle.half_edges;
    auto vertices = graph.adjacentTriangleVertices(triangle.half_edges[0]);

    if (vertices[0] == -1 || vertices[1] == -1 || vertices[2] == -1)
    {
      ++skip_infinite_face;
      continue;
    }

    if (!branch_strand_set.contains(static_cast<size_t>(vertices[0]))
      || !branch_strand_set.contains(static_cast<size_t>(vertices[1]))
      || !branch_strand_set.contains(static_cast<size_t>(vertices[2])))
    {
      ++skip_not_in_branch;
      continue;
    }

    if (strand_to_mesh_vertex[vertices[0]] < 0 || strand_to_mesh_vertex[vertices[1]] < 0
      || strand_to_mesh_vertex[vertices[2]] < 0)
    {
      ++skip_no_mesh_vertex;
      if (diagnostics
        && (static_cast<size_t>(vertices[0]) == 0 || static_cast<size_t>(vertices[1]) == 0
          || static_cast<size_t>(vertices[2]) == 0))
      {
        std::ostringstream oss;
        oss << "input_branch_id=" << input_branch_id << " face=" << face_index << " verts=[" << vertices[0] << ','
            << vertices[1] << ',' << vertices[2] << "] mesh_v=[" << strand_to_mesh_vertex[vertices[0]] << ','
            << strand_to_mesh_vertex[vertices[1]] << ',' << strand_to_mesh_vertex[vertices[2]]
            << "] reason=missing_boundary_mesh_vertex";
        strandInitDiagnosticLogLine("boundary_mesh_triangle_skip", 0, t, oss.str().c_str());
      }
      continue;
    }

    for (size_t i = 0; i < 3; i++)
    {
      const int left_mesh_vertex = strand_to_mesh_vertex[vertices[i]];
      const int right_mesh_vertex = strand_to_mesh_vertex[vertices[(i + 1) % 3]];
      if (left_mesh_vertex < 0 || right_mesh_vertex < 0)
      {
        continue;
      }

      if (kin_del.isOnComponentBoundaryOutside(he_ids[i]))
      {
        completeBoundaryMeshSection(
          he_ids[i], static_cast<size_t>(left_mesh_vertex), static_cast<size_t>(right_mesh_vertex), t);
        boundary_mesh_last_left_and_right_vertex[he_ids[i]]
          = std::make_pair(static_cast<size_t>(left_mesh_vertex), static_cast<size_t>(right_mesh_vertex));
      }
    }

    if (!kin_del.getFaceInside(face_index))
    {
      ++skip_outside_face;
      if (diagnostics
        && (static_cast<size_t>(vertices[0]) == 0 || static_cast<size_t>(vertices[1]) == 0
          || static_cast<size_t>(vertices[2]) == 0))
      {
        std::ostringstream oss;
        oss << "input_branch_id=" << input_branch_id << " face=" << face_index << " verts=[" << vertices[0] << ','
            << vertices[1] << ',' << vertices[2] << "] reason=face_not_inside";
        strandInitDiagnosticLogLine("boundary_mesh_triangle_skip", 0, t, oss.str().c_str());
      }
      continue;
    }

    if (invert_orientation)
    {
      std::swap(vertices[1], vertices[2]);
    }

    const size_t tri_v0 = static_cast<size_t>(strand_to_mesh_vertex[vertices[0]]);
    const size_t tri_v1 = static_cast<size_t>(strand_to_mesh_vertex[vertices[1]]);
    const size_t tri_v2 = static_cast<size_t>(strand_to_mesh_vertex[vertices[2]]);
    const std::string face_metadata
      = composeBoundaryMeshFaceMetadata(t, "boundary_delaunay", static_cast<size_t>(-1), face_index, input_branch_id);
    addBoundaryTriangle(tri_v0, tri_v1, tri_v2, face_metadata, 1);
    ++triangles_added;
  }

  if (diagnostics)
  {
    const bool branch_contains_strand_0 = branch_strand_set.contains(0);
    std::ostringstream oss;
    oss << "input_branch_id=" << input_branch_id << " invert_orientation=" << (invert_orientation ? "true" : "false")
        << " offset=" << offset << " branch_contains_strand_0=" << (branch_contains_strand_0 ? "true" : "false")
        << " branch_strand_count=" << branch_strands.size() << " skip_dummy=" << skip_dummy
        << " strands_vertex_added=" << strands_vertex_added << " skip_infinite_face=" << skip_infinite_face
        << " skip_not_in_branch=" << skip_not_in_branch << " skip_no_mesh_vertex=" << skip_no_mesh_vertex
        << " skip_outside_face=" << skip_outside_face << " triangles_added=" << triangles_added
        << " boundary_tris_before=" << boundary_tris_before
        << " boundary_tris_after=" << boundary_mesh.getTriangleCount();
    strandInitDiagnosticLogLine("boundary_mesh_summary",
      branch_contains_strand_0   ? 0
        : branch_strands.empty() ? 0
                                 : branch_strands.front(),
      t, oss.str().c_str());
    if (branch_contains_strand_0 && strands_vertex_added == 0)
    {
      strandInitDiagnosticLogLine(
        "boundary_mesh_strand0_missing", 0, t, "strand 0 is in branch but no boundary mesh vertex was added");
    }
  }
}

std::vector<BoundaryPoint> kinDS::SegmentBuilder::traceConvexHull(double t) const
{
  const auto& graph = kin_del.getGraph();
  std::vector<BoundaryPoint> convex_hull_points;
  for (HalfEdgeDelaunayGraph::ConvexHullEdgeIterator it = graph.boundaryEdgesBegin(), end = graph.boundaryEdgesEnd();
    it != end; ++it)
  {
    size_t he_id = *it;
    size_t strand_index = graph.halfEdge(he_id).origin;

    glm::dvec2 convex_hull_point = kin_del.getPointAt(t, strand_index, false, false);
    convex_hull_points.push_back({ strand_index, he_id, convex_hull_point });
  }

  return convex_hull_points;
}

void kinDS::SegmentBuilder::advanceBoundaryMesh(
  double t, const std::vector<BoundaryPoint>& boundary_points, const glm::dvec2& centroid)
{
  std::vector<size_t> new_vertex_indices;

  for (size_t i = 0; i < boundary_points.size(); i++)
  {
    glm::dvec2 boundary_point = boundary_points[i].p;

    new_vertex_indices.push_back(addBoundaryVertex(
      glm::dvec3 { boundary_point[0], boundary_point[1], t }, centroid, boundary_points[i].vertex_id, t, false));
  }

  for (size_t i = 0; i < boundary_points.size(); i++)
  {
    size_t he_id = boundary_points[i].he_id;
    size_t left_vertex_index = new_vertex_indices[i];
    size_t right_vertex_index = new_vertex_indices[(i + 1) % boundary_points.size()];
    auto& left_and_right = boundary_mesh_last_left_and_right_vertex[he_id];
    completeBoundaryMeshSection(he_id, left_vertex_index, right_vertex_index, t);

    // KINDS_DEBUG("Assigning boundary last vertices for he_id " << he_id);
    left_and_right.first = left_vertex_index;
    left_and_right.second = right_vertex_index;
  }
}

bool kinDS::SegmentBuilder::isComponentLive(size_t component_index) const
{
  if (component_index >= kin_del.component_data.components.size())
  {
    return false;
  }

  const auto& component = kin_del.component_data.components[component_index];
  if (component.empty())
  {
    return false;
  }

  const HalfEdgeDelaunayGraph& graph = kin_del.getGraph();
  for (size_t vertex : component)
  {
    if (graph.incidentEdgesBegin(vertex) != graph.incidentEdgesEnd(vertex))
    {
      return true;
    }
  }

  return false;
}

std::vector<size_t> kinDS::SegmentBuilder::collectLiveComponentIndices() const
{
  std::vector<size_t> live_component_indices;
  live_component_indices.reserve(kin_del.component_data.components.size());

  for (size_t component_index = 0; component_index < kin_del.component_data.components.size(); ++component_index)
  {
    if (isComponentLive(component_index))
    {
      live_component_indices.push_back(component_index);
    }
  }

  return live_component_indices;
}

void kinDS::SegmentBuilder::updateBoundary(double t, std::vector<bool>& visited, size_t component_index)
{
  if (component_index >= kin_del.component_data.components.size()
    || kin_del.component_data.components[component_index].empty() || !isComponentLive(component_index))
  {
    return;
  }

  if (kin_del.component_data.component_last_updated[component_index] != t)
  {
    kin_del.component_data.component_boundaries[component_index] = kin_del.extractComponentBoundaries(
      kin_del.component_data.components[component_index], t, visited, false, false);
    if (!kin_del.component_data.component_boundaries[component_index].empty()
      && !kin_del.component_data.component_boundaries[component_index][0].empty())
    {
      kin_del.component_data.component_centroids[component_index]
        = polygonCentroid(kin_del.component_data.component_boundaries[component_index][0]);
    }
    kin_del.component_data.component_last_updated[component_index] = t;
  }
}

void kinDS::SegmentBuilder::updateBoundaries(double t, const std::vector<size_t>& component_indices)
{
  std::vector<bool> visited(kin_del.getGraph().halfEdgeSlotCount(), false);

  for (size_t component_index : component_indices)
  {
    updateBoundary(t, visited, component_index);
  }
}

void kinDS::SegmentBuilder::advanceBoundaryMeshes(double t)
{
  const std::vector<size_t> live_component_indices = collectLiveComponentIndices();
  updateBoundaries(t, live_component_indices);

  for (size_t component_index : live_component_indices)
  {
    auto& boundaries = kin_del.component_data.component_boundaries[component_index];
    auto& centroid = kin_del.component_data.component_centroids[component_index];

    for (auto& boundary_points : boundaries)
    {
      advanceBoundaryMesh(t, boundary_points, centroid);
    }
  }
}

std::vector<glm::dvec3> SegmentBuilder::computeClampedVoronoiVertices(
  size_t strand_id, double t, const std::vector<BoundaryPoint>& boundary_polygon, const glm::dvec2& centroid)
{
  auto& graph = kin_del.getGraph();
  std::vector<glm::dvec3> voronoi_vertices;

  // we just create a triangle fan because the Voronoi cell is convex
  // iterate over all segment indices of the strand
  for (HalfEdgeDelaunayGraph::IncidentEdgeIterator it = graph.incidentEdgesBegin(strand_id),
                                                   end = graph.incidentEdgesEnd(strand_id);
    it != end; ++it)
  {
    glm::dvec3 voronoi_vertex = computeVoronoiVertex(*it, t);
    voronoi_vertices.emplace_back(voronoi_vertex);
    // addMeshletVertex(mesh, boundary_polygon, centroid, voronoi_vertex);
  }

  std::vector<glm::dvec3> clamped_voronoi_vertices;
  // go through the vertices and clamp them
  for (size_t i = 0; i < voronoi_vertices.size(); i++)
  {
    auto left = voronoi_vertices[i];
    auto right = voronoi_vertices[(i + 1) % voronoi_vertices.size()];

    bool inside = clampVoronoiVertices(left, right, boundary_polygon, centroid);

    if (inside)
    {
      clamped_voronoi_vertices.push_back(left);
      clamped_voronoi_vertices.push_back(right);
    }
  }

  voronoi_vertices.clear();
  // Put back into original vector and remove duplicates
  voronoi_vertices.emplace_back(clamped_voronoi_vertices.front());
  for (size_t i = 1; i < clamped_voronoi_vertices.size() - 1; i++)
  {
    if (clamped_voronoi_vertices[i] != voronoi_vertices.back())
    {
      voronoi_vertices.push_back(clamped_voronoi_vertices[i]);
    }
  }
  if (clamped_voronoi_vertices.front() != clamped_voronoi_vertices.back())
  {
    voronoi_vertices.push_back(clamped_voronoi_vertices.back());
  }

  return voronoi_vertices;
}

size_t kinDS::SegmentBuilder::closingMeshCountStrandIncidentEdges(size_t strand_id) const
{
  auto& graph = kin_del.getGraph();
  size_t num = 0;
  for (HalfEdgeDelaunayGraph::IncidentEdgeIterator it = graph.incidentEdgesBegin(strand_id),
                                                   end = graph.incidentEdgesEnd(strand_id);
    it != end; ++it)
  {
    ++num;
  }
  return num;
}

int kinDS::SegmentBuilder::closingMeshAppendVertex(VoronoiMesh& mesh,
  const std::vector<BoundaryPoint>& boundary_polygon, const glm::dvec2& centroid, size_t strand_id, double t,
  const glm::dvec3& position, bool includes_virtual_shift, std::optional<size_t> voronoi_vertex_for_alpha_check,
  const std::string& metadata, const MeshletVertexRuntimeInfo& runtime_info)
{
  const size_t vertex_id = addMeshletVertex(mesh, boundary_polygon, centroid, position, strand_id, t,
    includes_virtual_shift, voronoi_vertex_for_alpha_check, metadata, std::nullopt, runtime_info);
  return static_cast<int>(vertex_id);
}

std::optional<KineticDelaunay::CrossingData::EdgeIntersectionRef>
kinDS::SegmentBuilder::closingMeshFindVoronoiEdgeIntersection(
  size_t voronoi_edge_id, size_t crossed_delaunay_he_id) const
{
  const size_t d = crossed_delaunay_he_id / 2;
  const auto& crossing_data = kin_del.getCrossingData();
  if (voronoi_edge_id >= crossing_data.voronoi_edge_intersections.size())
  {
    return std::nullopt;
  }
  for (const auto& ref : crossing_data.voronoi_edge_intersections[voronoi_edge_id])
  {
    if (ref->delaunay_edge_id == d)
    {
      return ref;
    }
  }
  return std::nullopt;
}

void kinDS::SegmentBuilder::refreshMeshingDataCrossingRefs(MeshingData& seg, size_t voronoi_edge_id) const
{
  if (seg.start_half_edge_id >= 0)
  {
    seg.start_crossing
      = closingMeshFindVoronoiEdgeIntersection(voronoi_edge_id, static_cast<size_t>(seg.start_half_edge_id));
  }
  else
  {
    seg.start_crossing.reset();
  }
  if (seg.end_half_edge_id >= 0)
  {
    seg.end_crossing
      = closingMeshFindVoronoiEdgeIntersection(voronoi_edge_id, static_cast<size_t>(seg.end_half_edge_id));
  }
  else
  {
    seg.end_crossing.reset();
  }
}

void kinDS::SegmentBuilder::refreshStripCrossingRefsIncidentToVoronoiVertex(size_t voronoi_vertex_id)
{
  auto& graph = kin_del.getGraph();
  const auto& face = graph.face(voronoi_vertex_id);
  std::unordered_set<size_t> refreshed_pairs;
  for (size_t vhe : face.half_edges)
  {
    const size_t even_he = vhe & ~1;
    if (even_he >= half_edge_index_to_segment_mesh_pair_index.size())
    {
      continue;
    }
    const size_t pair_idx = half_edge_index_to_segment_mesh_pair_index[even_he];
    if (pair_idx >= segment_mesh_pair_last_left_and_right_vertex.size())
    {
      continue;
    }
    if (!refreshed_pairs.insert(pair_idx).second)
    {
      continue;
    }
    const size_t strip_ve = even_he / 2;
    for (auto& seg : segment_mesh_pair_last_left_and_right_vertex[pair_idx])
    {
      refreshMeshingDataCrossingRefs(seg, strip_ve);
    }
  }
}

void kinDS::SegmentBuilder::refreshCrossingRefsForAllStrips()
{
  for (size_t pair_idx = 0; pair_idx < segment_mesh_pair_last_left_and_right_vertex.size(); ++pair_idx)
  {
    size_t voronoi_edge_id = 0;
    bool found = false;
    for (size_t he = 0; he < half_edge_index_to_segment_mesh_pair_index.size(); ++he)
    {
      if (half_edge_index_to_segment_mesh_pair_index[he] == pair_idx)
      {
        voronoi_edge_id = (he & ~1) / 2;
        found = true;
        break;
      }
    }
    if (!found)
    {
      continue;
    }
    for (auto& seg : segment_mesh_pair_last_left_and_right_vertex[pair_idx])
    {
      refreshMeshingDataCrossingRefs(seg, voronoi_edge_id);
    }
  }
}

void kinDS::SegmentBuilder::refreshCrossingRefsForAllIntersectionStrips()
{
  for (size_t pair_idx = 0; pair_idx < intersection_mesh_pair_last_left_and_right_vertex.size(); ++pair_idx)
  {
    if (pair_idx >= intersection_mesh_pair_metadata.size())
    {
      continue;
    }

    const auto& md = intersection_mesh_pair_metadata[pair_idx];
    size_t voronoi_edge_id = md.start_delaunay_edge_id;
    if (voronoi_edge_id == static_cast<size_t>(-1))
    {
      voronoi_edge_id = md.end_delaunay_edge_id;
    }
    if (voronoi_edge_id == static_cast<size_t>(-1))
    {
      continue;
    }

    for (auto& seg : intersection_mesh_pair_last_left_and_right_vertex[pair_idx])
    {
      refreshMeshingDataCrossingRefs(seg, voronoi_edge_id);
    }
  }
}

void kinDS::SegmentBuilder::growGraphSlotArrays()
{
  const size_t he_slots = kin_del.getGraph().halfEdgeSlotCount();
  if (half_edge_index_to_segment_mesh_pair_index.size() < he_slots)
  {
    half_edge_index_to_segment_mesh_pair_index.resize(he_slots, static_cast<size_t>(-1));
  }
  if (corner_to_cutoff_mesh_indices.size() < he_slots)
  {
    corner_to_cutoff_mesh_indices.resize(he_slots, -1);
  }
  if (boundary_mesh_last_left_and_right_vertex.size() < he_slots)
  {
    boundary_mesh_last_left_and_right_vertex.resize(
      he_slots, std::make_pair(static_cast<size_t>(-1), static_cast<size_t>(-1)));
  }
  if (half_edge_to_boundary_vertex_index.size() < he_slots)
  {
    half_edge_to_boundary_vertex_index.resize(he_slots, -1);
  }
}

void kinDS::SegmentBuilder::clearDeadHalfEdgeState()
{
  const auto& graph = kin_del.getGraph();
  const size_t he_slots = graph.halfEdgeSlotCount();
  for (size_t he_id = 0; he_id < he_slots; ++he_id)
  {
    if (graph.isLiveHalfEdge(he_id))
    {
      continue;
    }
    if (he_id < half_edge_index_to_segment_mesh_pair_index.size())
    {
      half_edge_index_to_segment_mesh_pair_index[he_id] = static_cast<size_t>(-1);
    }
    if (he_id < corner_to_cutoff_mesh_indices.size())
    {
      corner_to_cutoff_mesh_indices[he_id] = -1;
    }
    if (he_id < boundary_mesh_last_left_and_right_vertex.size())
    {
      boundary_mesh_last_left_and_right_vertex[he_id]
        = std::make_pair(static_cast<size_t>(-1), static_cast<size_t>(-1));
    }
    if (he_id < half_edge_to_boundary_vertex_index.size())
    {
      half_edge_to_boundary_vertex_index[he_id] = -1;
    }
  }
}

void kinDS::SegmentBuilder::initializeNewHalfEdgesAfterGraphUpdate(double t, size_t first_new_he_slot)
{
  auto& graph = kin_del.getGraph();
  const auto& d_intersections_all = kin_del.getCrossingData().delaunay_edge_intersections;

  for (size_t he_id : graph.liveDelaunayEdges())
  {
    if (he_id < first_new_he_slot)
    {
      continue;
    }

    const bool has_pair = he_id < half_edge_index_to_segment_mesh_pair_index.size()
      && half_edge_index_to_segment_mesh_pair_index[he_id] != static_cast<size_t>(-1);
    if (!has_pair)
    {
      startNewMesh(he_id, t);
    }
  }

  for (size_t he_id : graph.liveDelaunayEdges())
  {
    if (he_id < first_new_he_slot)
    {
      continue;
    }
    if (!kin_del.isOnComponentBoundary(he_id))
    {
      continue;
    }

    const size_t d_edge_id = he_id / 2;
    if (d_edge_id >= d_intersections_all.size() || d_intersections_all[d_edge_id].empty())
    {
      continue;
    }

    const std::vector<KineticDelaunay::CrossingData::EdgeIntersectionRef> refs
      = getBoundaryIntersectionsInBoundaryOrder(d_edge_id);

    {
      const size_t first_cell
        = determineVoronoiCellForBoundaryIntersectionInterval(d_edge_id, std::nullopt, refs.front());
      startNewMeshFromIntersections(first_cell, t, std::nullopt, refs.front(), false, BoundaryEventType::Section,
        BoundarySegmentAction::NewSegment);
    }
    for (size_t k = 0; k + 1 < refs.size(); ++k)
    {
      const size_t mid_cell = determineVoronoiCellForBoundaryIntersectionInterval(d_edge_id, refs[k], refs[k + 1]);
      if (mid_cell == static_cast<size_t>(-1))
      {
        continue;
      }
      startNewMeshFromIntersections(
        mid_cell, t, refs[k], refs[k + 1], false, BoundaryEventType::Section, BoundarySegmentAction::NewSegment);
    }
    {
      const size_t last_cell
        = determineVoronoiCellForBoundaryIntersectionInterval(d_edge_id, refs.back(), std::nullopt);
      startNewMeshFromIntersections(
        last_cell, t, refs.back(), std::nullopt, false, BoundaryEventType::Section, BoundarySegmentAction::NewSegment);
    }
  }
}

void kinDS::SegmentBuilder::onGraphRetriangulated(double t, size_t prev_face_slots, size_t prev_he_slots)
{
  (void)prev_face_slots;
  growGraphSlotArrays();
  clearDeadHalfEdgeState();
  initializeNewHalfEdgesAfterGraphUpdate(t, prev_he_slots);
  refreshCrossingRefsForAllStrips();
  refreshCrossingRefsForAllIntersectionStrips();
}

void kinDS::SegmentBuilder::onGraphCutApplied(double t, size_t prev_face_slots, size_t prev_he_slots)
{
  (void)t;
  (void)prev_face_slots;
  (void)prev_he_slots;
  growGraphSlotArrays();
  clearDeadHalfEdgeState();
  refreshCrossingRefsForAllStrips();
  refreshCrossingRefsForAllIntersectionStrips();
}

void kinDS::SegmentBuilder::onBeforeComponentGraphSplit(double /*t*/) { }

glm::dvec3 kinDS::SegmentBuilder::closingMeshVoronoiDelaunayCrossingPosition(
  double t, size_t voronoi_edge_id, size_t delaunay_edge_id) const
{
  const auto& crossing_data = kin_del.getCrossingData();
  if (voronoi_edge_id < crossing_data.voronoi_edge_intersections.size())
  {
    for (const KineticDelaunay::CrossingData::EdgeIntersectionRef ref :
      crossing_data.voronoi_edge_intersections[voronoi_edge_id])
    {
      if (ref->delaunay_edge_id == delaunay_edge_id)
      {
        return getCrossingCoordsInMeshSpace(ref, t);
      }
    }
  }

  KINDS_WARNING("closingMeshVoronoiDelaunayCrossingPosition: no crossing ref for voronoi_edge="
    << voronoi_edge_id << " delaunay_edge=" << delaunay_edge_id << " t=" << t);
  return glm::dvec3(std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN(), t);
}

auto kinDS::SegmentBuilder::extractSegmentsForVoronoiEdge(double t, int incident_edge_index, size_t voronoi_edge_id,
  const std::function<int(
    const glm::dvec3&, std::optional<size_t>, const std::string&, const MeshletVertexRuntimeInfo&)>& track_vertex,
  bool reverse) -> std::vector<MeshingData>
{
  std::vector<MeshingData> out_segments;
  auto& graph = kin_del.getGraph();

  const size_t voronoi_he_even = 2 * voronoi_edge_id;
  const size_t voronoi_he_odd = voronoi_he_even + 1;

  if (graph.isInfinite(voronoi_he_even))
  {
    return out_segments;
  }

  // `reverse`: when false, walk even circumcenter → odd; when true, odd → even (strand-side-first in closingMesh).
  const size_t voronoi_he_strand_side = reverse ? voronoi_he_odd : voronoi_he_even;
  const size_t voronoi_he_other_side = voronoi_he_strand_side ^ 1;
  const bool strand_at_voronoi_even_he = !reverse;

  const glm::dvec3 strand_cm_vertex = computeVoronoiVertex(voronoi_he_strand_side, t);
  const glm::dvec3 other_cm_vertex = computeVoronoiVertex(voronoi_he_other_side, t);
  const glm::dvec2 strand2 = glm::dvec2(strand_cm_vertex);
  const glm::dvec2 other2 = glm::dvec2(other_cm_vertex);
  const glm::dvec2 vor_dir = other2 - strand2;
  const double vor_len2 = glm::dot(vor_dir, vor_dir);
  if (vor_len2 < 1e-14)
  {
    return out_segments;
  }

  const size_t strand_face_id = graph.halfEdge(voronoi_he_strand_side).face;
  const size_t strand_containing_tri_id = kin_del.getCrossingDataContainingTriId(strand_face_id);
  bool inside = kin_del.getFaceInside(strand_containing_tri_id);
  if (inside)
  {
    const size_t voronoi_vertex_id = graph.halfEdge(voronoi_he_strand_side).face;
    out_segments.push_back(MeshingData {
      track_vertex(strand_cm_vertex, voronoi_vertex_id,
        makeClosingMeshVertexMetadata("Voronoi vertex", std::nullopt, voronoi_vertex_id), MeshletVertexRuntimeInfo {}),
      -1, -1, -1 });
    out_segments.back().closing_incident_edge_index = incident_edge_index;
    out_segments.back().closing_voronoi_edge_id = voronoi_edge_id;
    out_segments.back().closing_strand_at_voronoi_even_he = strand_at_voronoi_even_he;
  }

  if (!kin_del.computeBoundaryOnTheFly())
  {
    if (inside && !out_segments.empty())
    {
      const size_t voronoi_vertex_id = graph.halfEdge(voronoi_he_other_side).face;
      out_segments.back().mesh_end_vertex_id = track_vertex(other_cm_vertex, voronoi_vertex_id,
        makeClosingMeshVertexMetadata("Voronoi vertex", std::nullopt, voronoi_vertex_id), MeshletVertexRuntimeInfo {});
    }
    for (auto& seg : out_segments)
    {
      refreshMeshingDataCrossingRefs(seg, voronoi_edge_id);
    }
    return out_segments;
  }

  const std::vector<DirectedVoronoiEdgeCrossing> directed_crossings
    = orientCrossingsAlongVoronoiEdge(voronoi_edge_id, strand_containing_tri_id, reverse);

  for (const DirectedVoronoiEdgeCrossing& directed : directed_crossings)
  {
    const size_t crossed_he_id = directed.crossed_half_edge_id;
    const auto& crossing_ref = directed.ref;
    const size_t next_face_id = graph.halfEdge(crossed_he_id ^ 1).face;
    const bool next_inside = kin_del.getFaceInside(next_face_id);

    if (inside != next_inside)
    {
      const glm::dvec3 crossing_pos = getCrossingCoordsInMeshSpace(crossing_ref, t);
      if (!std::isfinite(crossing_pos.x) || !std::isfinite(crossing_pos.y))
      {
        continue;
      }

      const int mesh_vertex_id
        = track_vertex(crossing_pos, std::nullopt, makeClosingMeshVertexMetadata("intersection", crossing_ref),
          MeshletVertexRuntimeInfo { false, false, crossing_ref, crossing_ref });

      if (inside)
      {
        if (!out_segments.empty())
        {
          auto& segment = out_segments.back();
          segment.mesh_end_vertex_id = mesh_vertex_id;
          segment.end_half_edge_id = static_cast<int>(crossed_he_id);
          segment.end_crossing = crossing_ref;
        }
      }
      else
      {
        MeshingData segment { mesh_vertex_id, -1, static_cast<int>(crossed_he_id ^ 1), -1 };
        segment.start_crossing = crossing_ref;
        segment.closing_incident_edge_index = incident_edge_index;
        segment.closing_voronoi_edge_id = voronoi_edge_id;
        segment.closing_strand_at_voronoi_even_he = strand_at_voronoi_even_he;
        out_segments.push_back(segment);
      }
    }

    inside = next_inside;
  }

  if (inside && !out_segments.empty())
  {
    const size_t voronoi_vertex_id = graph.halfEdge(voronoi_he_other_side).face;
    out_segments.back().mesh_end_vertex_id = track_vertex(other_cm_vertex, voronoi_vertex_id,
      makeClosingMeshVertexMetadata("Voronoi vertex", std::nullopt, voronoi_vertex_id), MeshletVertexRuntimeInfo {});
  }

  return out_segments;
}

auto kinDS::SegmentBuilder::closingMeshExtractRawSegmentsForVoronoiEdge(size_t strand_id, double t, VoronoiMesh& mesh,
  const std::vector<BoundaryPoint>& boundary_polygon, const glm::dvec2& centroid, int incident_edge_index,
  size_t incident_he,
  const std::function<int(const glm::dvec3&, std::optional<size_t>, const std::string&,
    const MeshletVertexRuntimeInfo&)>& track_vertex) -> std::vector<MeshingData>
{
  auto& graph = kin_del.getGraph();

  // Undirected Delaunay edge id (matches `CrossingData` / `computeEdgeIntersections`: voronoi_edge_id == he/2).
  const size_t voronoi_edge_id = incident_he / 2;
  const size_t voronoi_he_even = 2 * voronoi_edge_id;
  const size_t voronoi_he_odd = voronoi_he_even + 1;

  if (graph.isInfinite(voronoi_he_even))
  {
    return {};
  }

  const int origin_even = graph.halfEdge(voronoi_he_even).origin;
  const int origin_odd = graph.halfEdge(voronoi_he_odd).origin;
  const auto finite_origin_is_strand = [&](int o) { return o >= 0 && static_cast<size_t>(o) == strand_id; };
  if (!finite_origin_is_strand(origin_even) && !finite_origin_is_strand(origin_odd))
  {
    KINDS_ERROR("Closing mesh: Voronoi edge for dual Delaunay edge_id "
      << voronoi_edge_id << " (half-edges " << voronoi_he_even << ", " << voronoi_he_odd << ", incident_he "
      << incident_he << ") has no finite half-edge with origin strand_id " << strand_id << " (origins " << origin_even
      << ", " << origin_odd << ")");
  }

  // Orient the inside trace from the Voronoi half-edge whose origin is the strand. `voronoi_edge_intersections` is
  // stored in even→odd order; when the strand is on the odd half-edge we traverse that list in reverse.
  const bool strand_at_voronoi_even_he = finite_origin_is_strand(origin_even);
  const bool reverse = !strand_at_voronoi_even_he;

  return extractSegmentsForVoronoiEdge(t, incident_edge_index, voronoi_edge_id, track_vertex, reverse);
}

SegmentBuilder::ClosingMeshOrderedIndex kinDS::SegmentBuilder::closingMeshBuildOrderedIndex(
  std::list<MeshingData>& closing_segments)
{
  ClosingMeshOrderedIndex out;
  out.ordered_segments.reserve(closing_segments.size());
  for (auto& seg : closing_segments)
  {
    if (seg.mesh_start_vertex_id == -1 || seg.mesh_end_vertex_id == -1)
    {
      continue;
    }
    out.ordered_segments.push_back(&seg);
  }

  out.start_crossing_to_segment.reserve(out.ordered_segments.size());
  for (size_t i = 0; i < out.ordered_segments.size(); ++i)
  {
    MeshingData* s = out.ordered_segments[i];
    if (s->start_crossing.has_value())
    {
      const auto ref = s->start_crossing.value();
      const auto ins = out.start_crossing_to_segment.insert({ ref, i });
      if (!ins.second)
      {
        KINDS_ERROR("Closing mesh: duplicate start_crossing: ordered_seg[" << ins.first->second << "] and ordered_seg["
                                                                           << i << "] share the same start crossing.");
      }
    }
  }
  return out;
}

void kinDS::SegmentBuilder::closingMeshValidateOrderedSegmentGeometry(
  double t, const VoronoiMesh& mesh, const std::vector<MeshingData*>& ordered_segments) const
{
  return; // disable for now
  constexpr double k_closing_cap_geom_eps = 1e-5;
  auto mesh_vertex_xy = [&](int mesh_vid) -> glm::dvec2
  {
    const glm::dvec3& v = mesh.getVertices()[static_cast<size_t>(mesh_vid)];
    return glm::dvec2(v.x, v.y);
  };

  for (size_t si = 0; si < ordered_segments.size(); ++si)
  {
    const MeshingData* s = ordered_segments[si];
    if (s->closing_voronoi_edge_id == static_cast<size_t>(-1))
    {
      continue;
    }
    const size_t ve = s->closing_voronoi_edge_id;
    {
      const std::string ctxs = std::string("ordered_seg[") + std::to_string(si)
        + "] incident=" + std::to_string(s->closing_incident_edge_index);
      validateClosingCapCrossingRef(
        kin_del, (ctxs + " start_ref").c_str(), s->start_crossing, ve, s->start_half_edge_id);
      validateClosingCapCrossingRef(kin_del, (ctxs + " end_ref").c_str(), s->end_crossing, ve, s->end_half_edge_id);
    }

    const glm::dvec3 L3 = kin_del.computeVoronoiVertexClampedInfinity(2 * ve, t, false, false);
    const glm::dvec3 R3 = kin_del.computeVoronoiVertexClampedInfinity(2 * ve + 1, t, false, false);
    const glm::dvec2 L(L3.x, L3.y);
    const glm::dvec2 R(R3.x, R3.y);
    const glm::dvec2 axis = R - L;
    const double ax2 = glm::dot(axis, axis);
    if (ax2 < 1e-24)
    {
      KINDS_ERROR("Closing mesh ordered_seg[" << si << "]: degenerate Voronoi edge " << ve << " (zero length).");
      continue;
    }
    const glm::dvec2 ps = mesh_vertex_xy(s->mesh_start_vertex_id);
    const glm::dvec2 pe = mesh_vertex_xy(s->mesh_end_vertex_id);
    const double axis_len = glm::length(axis);
    const auto perp_dist
      = [&](const glm::dvec2& p) -> double { return std::abs((p.x - L.x) * axis.y - (p.y - L.y) * axis.x) / axis_len; };
    if (perp_dist(ps) > k_closing_cap_geom_eps || perp_dist(pe) > k_closing_cap_geom_eps)
    {
      KINDS_ERROR("Closing mesh ordered_seg[" << si << "]: mesh_start/end not collinear with canonical Voronoi "
                                              << "edge " << ve
                                              << " (even->odd circumcenters); check segment direction vs inside walk.");
    }
    const double ts = glm::dot(ps - L, axis) / ax2;
    const double te = glm::dot(pe - L, axis) / ax2;
    // mesh_start is always the strand-side circumcenter; along the canonical even->odd axis that is the smaller t iff
    // the strand Voronoi half-edge is even.
    if (s->closing_strand_at_voronoi_even_he)
    {
      if (ts > te + 1e-8)
      {
        KINDS_ERROR("Closing mesh ordered_seg[" << si << "]: mesh vertices reversed along Voronoi edge " << ve
                                                << " (t_start=" << ts << " > t_end=" << te
                                                << "); expected mesh_start nearer voronoi_he_even.");
      }
    }
    else if (te > ts + 1e-8)
    {
      KINDS_ERROR("Closing mesh ordered_seg[" << si << "]: mesh vertices reversed along Voronoi edge " << ve
                                              << " (t_start=" << ts << " > t_end=" << te
                                              << " along even->odd); expected mesh_start nearer "
                                              << "voronoi_he_odd (strand-side).");
    }

    glm::dvec2 ip;
    if (s->start_crossing.has_value())
    {
      const glm::dvec3 crossing_pos = getCrossingCoordsInMeshSpace(s->start_crossing.value(), t);
      if (std::isfinite(crossing_pos.x) && std::isfinite(crossing_pos.y))
      {
        ip = glm::dvec2(crossing_pos);
        if (glm::distance(ip, ps) > k_closing_cap_geom_eps)
        {
          KINDS_ERROR("Closing mesh ordered_seg[" << si << "]: start_ref 2D position does not match mesh_start.");
        }
      }
    }
    if (s->end_crossing.has_value())
    {
      const glm::dvec3 crossing_pos = getCrossingCoordsInMeshSpace(s->end_crossing.value(), t);
      if (std::isfinite(crossing_pos.x) && std::isfinite(crossing_pos.y))
      {
        ip = glm::dvec2(crossing_pos);
        if (glm::distance(ip, pe) > k_closing_cap_geom_eps)
        {
          KINDS_ERROR("Closing mesh ordered_seg[" << si << "]: end_ref 2D position does not match mesh_end.");
        }
      }
    }
  }
}

SegmentBuilder::ClosingMeshPolygonsTraceResult kinDS::SegmentBuilder::closingMeshTraceCapPolygons(size_t strand_id,
  double t, size_t num_incident_edges, VoronoiMesh& mesh, const std::vector<BoundaryPoint>& boundary_polygon,
  const glm::dvec2& centroid, std::vector<size_t> mesh_vertex_ids, const std::vector<MeshingData*>& ordered_segments,
  const std::unordered_map<KineticDelaunay::CrossingData::EdgeIntersectionRef, size_t, ClosingMeshCrossingIteratorHash,
    ClosingMeshCrossingIteratorEq>& start_crossing_to_segment)
{
  ClosingMeshPolygonsTraceResult result;
  result.mesh_vertex_ids = std::move(mesh_vertex_ids);
  result.segment_used.assign(ordered_segments.size(), false);

  auto& graph = kin_del.getGraph();
  const auto& crossing_data = kin_del.getCrossingData();

  struct CapOutlineVertexInfo
  {
    const char* source = "unknown";
    size_t primary_id = static_cast<size_t>(-1);
    size_t secondary_id = static_cast<size_t>(-1);
    glm::dvec2 delaunay_xy { std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN() };
  };
  std::unordered_map<size_t, CapOutlineVertexInfo> outline_vertex_info;
  const auto record_intersection = [&](size_t mesh_vertex_id, KineticDelaunay::CrossingData::EdgeIntersectionRef ref)
  {
    const glm::dvec3 p = getCrossingCoordsInDelaunaySpace(kin_del, ref, t);
    outline_vertex_info[mesh_vertex_id]
      = { "intersection", ref->delaunay_edge_id, ref->voronoi_edge_id, glm::dvec2(p.x, p.y) };
  };
  const auto record_voronoi_vertex = [&](size_t mesh_vertex_id, size_t voronoi_vertex_id)
  {
    glm::dvec2 p(std::numeric_limits<double>::quiet_NaN());
    if (voronoi_vertex_id < graph.faceSlotCount() && graph.isLiveFace(voronoi_vertex_id))
    {
      const glm::dvec3 p3 = kin_del.computeVoronoiVertexClampedInfinity(graph.face(voronoi_vertex_id).half_edges[0], t);
      p = glm::dvec2(p3.x, p3.y);
    }
    outline_vertex_info[mesh_vertex_id] = { "Voronoi vertex", voronoi_vertex_id, static_cast<size_t>(-1), p };
  };
  const auto record_site = [&](size_t mesh_vertex_id, size_t site_id)
  {
    const glm::dvec2 p = kin_del.getPointInDelaunaySpace(site_id, t);
    outline_vertex_info[mesh_vertex_id] = { "site", site_id, static_cast<size_t>(-1), p };
  };

  auto& cap_vertex_ids = result.mesh_vertex_ids;
  auto add_mesh_vertex
    = [&](const glm::dvec3& v, const std::string& metadata, const MeshletVertexRuntimeInfo& runtime_info) -> int
  {
    const int id = closingMeshAppendVertex(
      mesh, boundary_polygon, centroid, strand_id, t, v, false, std::nullopt, metadata, runtime_info);
    cap_vertex_ids.push_back(static_cast<size_t>(id));
    return id;
  };

  auto polygon_push = [&](std::vector<size_t>& poly, size_t vertex_id)
  {
    if (poly.empty() || poly.back() != vertex_id)
    {
      poly.push_back(vertex_id);
    }
  };

  // Only connect along a Delaunay boundary edge to intersections whose Voronoi edge is dual to an edge incident to
  // this strand; otherwise the first CrossingData hit on that Delaunay edge can belong to another cell and we
  // incorrectly skip the strand's Voronoi continuation.
  auto voronoi_edge_incident_to_strand = [&](size_t voronoi_edge_id) -> bool
  {
    const size_t he0 = 2 * voronoi_edge_id;
    const size_t he1 = he0 + 1;
    if (he1 >= graph.halfEdgeSlotCount())
    {
      return false;
    }
    const int o0 = graph.halfEdge(he0).origin;
    const int o1 = graph.halfEdge(he1).origin;
    auto finite_origin_is_strand = [&](int o) { return o >= 0 && static_cast<size_t>(o) == strand_id; };
    return finite_origin_is_strand(o0) || finite_origin_is_strand(o1);
  };

  // Dual graph: Voronoi half-edge `.face` is the incident Voronoi vertex (circumcenter) id — stable across strips;
  // `mesh_*_vertex_id` are mesh buffer indices and must not be used to match circumcenter joins.
  auto closing_strip_strand_voronoi_half_edge = [](const MeshingData& s) -> size_t
  {
    if (s.closing_voronoi_edge_id == static_cast<size_t>(-1))
    {
      return static_cast<size_t>(-1);
    }
    const size_t ve = s.closing_voronoi_edge_id;
    return s.closing_strand_at_voronoi_even_he ? (2 * ve) : (2 * ve + 1);
  };
  auto closing_strip_voronoi_vertex_at_strand_endpoint = [&](const MeshingData& s) -> size_t
  {
    const size_t he_strand = closing_strip_strand_voronoi_half_edge(s);
    if (he_strand == static_cast<size_t>(-1) || he_strand >= graph.halfEdgeSlotCount())
    {
      return static_cast<size_t>(-1);
    }
    return graph.halfEdge(he_strand).face;
  };
  auto closing_strip_voronoi_vertex_at_other_endpoint = [&](const MeshingData& s) -> size_t
  {
    const size_t he_strand = closing_strip_strand_voronoi_half_edge(s);
    if (he_strand == static_cast<size_t>(-1) || he_strand >= graph.halfEdgeSlotCount())
    {
      return static_cast<size_t>(-1);
    }
    return graph.halfEdge(he_strand ^ 1).face;
  };

  for (;;)
  {
    size_t seed_index = static_cast<size_t>(-1);
    for (size_t i = 0; i < ordered_segments.size(); ++i)
    {
      if (!result.segment_used[i])
      {
        seed_index = i;
        break;
      }
    }
    if (seed_index == static_cast<size_t>(-1))
    {
      break;
    }

    std::vector<size_t> polygon;
    MeshingData* seed_seg = ordered_segments[seed_index];
    {
      const size_t ve_seed = seed_seg->closing_voronoi_edge_id;
      if (ve_seed == static_cast<size_t>(-1))
      {
        KINDS_DEBUG("Starting polygon walk from ordered segment " << seed_index << " of Voronoi edge ?");
      }
      else
      {
        KINDS_DEBUG("Starting polygon walk from ordered segment " << seed_index << " of Voronoi edge " << ve_seed);
      }
    }

    size_t current_segment_index = seed_index;
    size_t last_walk_segment_index = seed_index;
    size_t guard = 0;

    std::vector<int> cap_start_hits(ordered_segments.size(), 0);
    std::vector<int> cap_end_handed_off(ordered_segments.size(), 0);
    std::vector<size_t> walk_segment_order;
    walk_segment_order.reserve(ordered_segments.size());
    bool exited_fail = false;
    bool exited_closed_loop = false;
    std::optional<KineticDelaunay::CrossingData::EdgeIntersectionRef> walk_closure_crossing_ref;

    auto record_cap_handoff = [&](size_t from_idx, size_t to_idx)
    {
      cap_end_handed_off[from_idx]++;
      cap_start_hits[to_idx]++;
    };

    while (guard++ < ordered_segments.size() * 4 && !exited_fail)
    {
      bool closed_walk_at_seed = false;
      MeshingData* current_segment = ordered_segments[current_segment_index];
      last_walk_segment_index = current_segment_index;
      result.segment_used[current_segment_index] = true;
      walk_segment_order.push_back(current_segment_index);

      const size_t start_mesh_vertex_id = static_cast<size_t>(current_segment->mesh_start_vertex_id);
      polygon_push(polygon, start_mesh_vertex_id);
      if (current_segment->start_crossing.has_value())
      {
        record_intersection(start_mesh_vertex_id, current_segment->start_crossing.value());
      }
      else
      {
        record_voronoi_vertex(start_mesh_vertex_id, closing_strip_voronoi_vertex_at_strand_endpoint(*current_segment));
      }
      logClosingMeshStripEndpointVertex(
        kin_del, current_segment_index, current_segment->start_crossing, current_segment->closing_voronoi_edge_id);
      const size_t end_mesh_vertex_id = static_cast<size_t>(current_segment->mesh_end_vertex_id);
      polygon_push(polygon, end_mesh_vertex_id);
      if (current_segment->end_crossing.has_value())
      {
        record_intersection(end_mesh_vertex_id, current_segment->end_crossing.value());
      }
      else
      {
        record_voronoi_vertex(end_mesh_vertex_id, closing_strip_voronoi_vertex_at_other_endpoint(*current_segment));
      }
      logClosingMeshStripEndpointVertex(
        kin_del, current_segment_index, current_segment->end_crossing, current_segment->closing_voronoi_edge_id);

      const bool ends_at_boundary_crossing
        = current_segment->end_crossing.has_value() && current_segment->end_half_edge_id >= 0;

      if (!ends_at_boundary_crossing)
      {
        const size_t joint_vv = closing_strip_voronoi_vertex_at_other_endpoint(*current_segment);
        if (joint_vv == static_cast<size_t>(-1))
        {
          KINDS_ERROR("Closing mesh: ordered_seg[" << current_segment_index
                                                   << "] has no boundary end_crossing but closing_voronoi_edge_id is "
                                                      "unset; cannot resolve circumcenter Voronoi "
                                                      "vertex for join.");
          exited_fail = true;
          break;
        }
        size_t next_j = static_cast<size_t>(-1);
        for (size_t j = 0; j < ordered_segments.size(); ++j)
        {
          if (j == current_segment_index || result.segment_used[j])
          {
            continue;
          }
          if (closing_strip_voronoi_vertex_at_strand_endpoint(*ordered_segments[j]) == joint_vv)
          {
            next_j = j;
            break;
          }
        }
        if (next_j != static_cast<size_t>(-1))
        {
          record_cap_handoff(current_segment_index, next_j);
          KINDS_DEBUG("Polygon walk: Voronoi-vertex join from ordered_seg["
            << current_segment_index << "] to ordered_seg[" << next_j << "] at Voronoi vertex (circumcenter) id "
            << joint_vv << " (strip end at circumcenter, no boundary end_crossing).");
          current_segment_index = next_j;
          continue;
        }
        const size_t seed_start_vv = closing_strip_voronoi_vertex_at_strand_endpoint(*seed_seg);
        if (joint_vv == seed_start_vv && current_segment_index != seed_index
          && seed_start_vv != static_cast<size_t>(-1))
        {
          walk_closure_crossing_ref.reset();
          record_cap_handoff(current_segment_index, seed_index);
          KINDS_DEBUG("Polygon walk: Voronoi-vertex join closes loop to seed ordered_seg["
            << seed_index << "] at Voronoi vertex (circumcenter) id " << joint_vv << ".");
          exited_closed_loop = true;
          break;
        }
        KINDS_ERROR("Closing mesh: ordered_seg["
          << current_segment_index << "] ends at circumcenter Voronoi vertex id " << joint_vv
          << " with no boundary end_crossing, but no unused strip has that id as its strand-side Voronoi vertex (and "
             "it "
             "does not match seed's strand-side Voronoi vertex); cannot continue cap walk.");
        exited_fail = true;
        break;
      }

      const int exit_he_id = current_segment->end_half_edge_id;
      std::optional<KineticDelaunay::CrossingData::EdgeIntersectionRef> current_ref_opt = current_segment->end_crossing;

      const size_t start_boundary_he = kin_del.isOnComponentBoundaryOutside(static_cast<size_t>(exit_he_id))
        ? static_cast<size_t>(exit_he_id)
        : static_cast<size_t>(exit_he_id) ^ 1;

      size_t boundary_he = start_boundary_he;
      bool found_next_segment = false;
      size_t boundary_guard = 0;

      constexpr double k_strand_corner_dist_eps = 1e-9;
      auto append_strand_corner_if_needed = [&]()
      {
        const glm::dvec2 xy = kin_del.getPointAt(t, strand_id, false, false);
        const glm::dvec3 strand_pos(xy, t);
        if (!polygon.empty())
        {
          const glm::dvec3& last = mesh.getVertices()[polygon.back()];
          if (glm::distance(glm::dvec2(last.x, last.y), xy) <= k_strand_corner_dist_eps)
          {
            return;
          }
        }
        const int vid = add_mesh_vertex(strand_pos,
          makeClosingMeshVertexMetadata("site", std::nullopt, std::nullopt, strand_id), MeshletVertexRuntimeInfo {});
        polygon_push(polygon, static_cast<size_t>(vid));
        record_site(static_cast<size_t>(vid), strand_id);
        KINDS_DEBUG("Adding vertex at strand " << strand_id);
      };

      while (boundary_guard++ < graph.halfEdgeSlotCount() * 2)
      {
        const int he_origin = graph.halfEdge(boundary_he).origin;
        if (he_origin >= 0 && static_cast<size_t>(he_origin) == strand_id)
        {
          append_strand_corner_if_needed();
        }

        const auto& d_intersections = crossing_data.delaunay_edge_intersections[boundary_he / 2];
        // `delaunay_edge_intersections[e]` is ordered along the even Delaunay half-edge (2*e).
        // For a boundary edge, select the inside-oriented half-edge; traverse forward iff that inside half-edge is
        // even, otherwise traverse backward.
        const size_t inside_boundary_he
          = kin_del.isOnComponentBoundaryOutside(boundary_he) ? (boundary_he ^ 1) : boundary_he;
        const bool effective_list_forward = ((inside_boundary_he % 2) != 0);

        // `traced_boundary_intervals` record endpoints in canonical interval order (even-half-edge / increasing list
        // index). The walk may list crossings backward when `effective_list_forward` is false — then swap walk order
        // (a_walk -> b_walk) to (start_intersection, end_intersection).
        const auto append_traced_boundary_interval
          = [&](const std::optional<KineticDelaunay::CrossingData::EdgeIntersectionRef>& a_walk,
              const std::optional<KineticDelaunay::CrossingData::EdgeIntersectionRef>& b_walk)
        {
          if (effective_list_forward)
          {
            result.traced_boundary_intervals.push_back(BoundaryIntersectionInterval { strand_id, a_walk, b_walk });
          }
          else
          {
            result.traced_boundary_intervals.push_back(BoundaryIntersectionInterval { strand_id, b_walk, a_walk });
          }
        };

        // Walk every crossing along this directed Delaunay half-edge after `current_ref_opt`: each crossing gets a
        // mesh vertex. Strand-incident crossings that match an ordered segment start_crossing hand off to that segment;
        // others only advance along the same Delaunay edge (no skipping to the corner).
        bool exit_boundary_chain_to_voronoi = false;
        if (!d_intersections.empty())
        {
          std::optional<KineticDelaunay::CrossingData::EdgeIntersectionRef> last_ref_on_edge = std::nullopt;
          if (current_ref_opt.has_value() && current_ref_opt.value()->delaunay_edge_id == (boundary_he / 2))
          {
            last_ref_on_edge = current_ref_opt;
          }

          auto find_next_intersection = [&](std::optional<KineticDelaunay::CrossingData::EdgeIntersectionRef> ref_opt)
            -> std::list<KineticDelaunay::CrossingData::EdgeIntersectionRef>::const_iterator
          {
            if (d_intersections.empty())
            {
              return d_intersections.end();
            }
            if (!ref_opt.has_value())
            {
              if (effective_list_forward)
              {
                return d_intersections.begin();
              }
              return std::prev(d_intersections.end());
            }
            auto it = std::find(d_intersections.begin(), d_intersections.end(), ref_opt.value());
            if (it == d_intersections.end())
            {
              return d_intersections.end();
            }
            if (effective_list_forward)
            {
              return std::next(it);
            }
            if (it == d_intersections.begin())
            {
              return d_intersections.end();
            }
            return std::prev(it);
          };

          for (size_t hop = 0; hop <= d_intersections.size(); ++hop)
          {
            const auto next_it = find_next_intersection(current_ref_opt);
            if (next_it == d_intersections.end())
            {
              break;
            }

            const KineticDelaunay::CrossingData::EdgeIntersectionRef candidate_ref = *next_it;
            const std::string candidate_log
              = formatCrossingIntersectionForLog(kin_del, std::make_optional(candidate_ref));

            append_traced_boundary_interval(last_ref_on_edge, candidate_ref);
            last_ref_on_edge = candidate_ref;

            const bool strand_inc = voronoi_edge_incident_to_strand(candidate_ref->voronoi_edge_id);
            const auto seg_it_start
              = strand_inc ? start_crossing_to_segment.find(candidate_ref) : start_crossing_to_segment.end();
            size_t next_segment_index = static_cast<size_t>(-1);
            if (seg_it_start != start_crossing_to_segment.end())
            {
              next_segment_index = seg_it_start->second;
            }
            if (next_segment_index == current_segment_index && next_segment_index != seed_index)
            {
              KINDS_DEBUG("Skipping intersection "
                << candidate_log << " during walk: handoff would target current segment's own start.");
              KINDS_ERROR("Closing mesh: Delaunay handoff from ordered_seg["
                << current_segment_index
                << "] targets its own start_crossing (only the seed loop may close onto seed).");
              exited_fail = true;
              exit_boundary_chain_to_voronoi = true;
              break;
            }
            // Close the loop only when Delaunay walk reaches the seed strip's *start* crossing (end of prior strip
            // hands off to start of next; the seed's start receives the last handoff).
            const bool closes_polygon_at_seed_crossing = (next_segment_index == seed_index);

            // Walking a new boundary half-edge resets current_ref_opt to nullopt; the first listed crossing on that
            // edge can be the seed segment's start crossing — already emitted as mesh_start. Closing the loop must not
            // add another mesh vertex there.
            if (!closes_polygon_at_seed_crossing)
            {
              const glm::dvec3 inter_pos = getCrossingCoordsInMeshSpace(candidate_ref, t);
              if (std::isfinite(inter_pos.x) && std::isfinite(inter_pos.y))
              {
                const int nv = add_mesh_vertex(inter_pos, makeClosingMeshVertexMetadata("intersection", candidate_ref),
                  MeshletVertexRuntimeInfo { false, false, candidate_ref, candidate_ref });
                polygon_push(polygon, static_cast<size_t>(nv));
                record_intersection(static_cast<size_t>(nv), candidate_ref);
                logClosingMeshIntersectionVertex(kin_del, candidate_ref);
              }
              else
              {
                KINDS_DEBUG("Skipping intersection " << candidate_log
                                                     << " during walk: computed intersection position is non-finite.");
                KINDS_WARNING("Closing mesh: crossing on Delaunay edge has non-finite 2D position (boundary_he="
                  << boundary_he << ", crossing=" << candidate_log << ").");
              }
            }
            else
            {
              KINDS_DEBUG("Skipping intersection " << candidate_log
                                                   << " during walk: closes loop at seed start crossing; "
                                                      "seed start vertex already on polygon.");
            }

            if (strand_inc)
            {
              if (next_segment_index != static_cast<size_t>(-1))
              {
                if (closes_polygon_at_seed_crossing)
                {
                  walk_closure_crossing_ref = candidate_ref;
                  record_cap_handoff(current_segment_index, seed_index);
                  found_next_segment = true;
                  closed_walk_at_seed = true;
                  exit_boundary_chain_to_voronoi = true;
                  break;
                }
                record_cap_handoff(current_segment_index, next_segment_index);
                found_next_segment = true;
                if (next_segment_index == seed_index)
                {
                  closed_walk_at_seed = true;
                }
                else
                {
                  current_segment_index = next_segment_index;
                }
                exit_boundary_chain_to_voronoi = true;
                break;
              }
              KINDS_ERROR(
                "Closing mesh: strand-incident Voronoi crossing on Delaunay boundary matches no ordered segment "
                "start_crossing at boundary_he="
                << boundary_he << ".");
              KINDS_DEBUG("Skipping intersection " << candidate_log
                                                   << " during walk: strand-incident crossing has no matching segment "
                                                      "start_crossing.");
              exited_fail = true;
              exit_boundary_chain_to_voronoi = true;
              break;
            }

            current_ref_opt = candidate_ref;
          }
        }

        if (exit_boundary_chain_to_voronoi)
        {
          break;
        }

        const int corner_vertex_id = graph.destination(boundary_he);
        if (!d_intersections.empty())
        {
          std::optional<KineticDelaunay::CrossingData::EdgeIntersectionRef> tail_ref = std::nullopt;
          if (current_ref_opt.has_value() && current_ref_opt.value()->delaunay_edge_id == (boundary_he / 2))
          {
            tail_ref = current_ref_opt;
          }
          if (tail_ref.has_value())
          {
            append_traced_boundary_interval(tail_ref, std::nullopt);
          }
        }
        if (corner_vertex_id >= 0)
        {
          if (static_cast<size_t>(corner_vertex_id) != strand_id)
          {
            KINDS_ERROR("Closing mesh: Delaunay boundary walk (half-edge "
              << boundary_he << ") reached corner vertex " << corner_vertex_id << " (expected strand " << strand_id
              << ").");
            exited_fail = true;
            const glm::dvec2 corner = kin_del.getPointAt(t, static_cast<size_t>(corner_vertex_id), false, false);
            const int cv = add_mesh_vertex(glm::dvec3(corner, t),
              makeClosingMeshVertexMetadata("site", std::nullopt, std::nullopt, static_cast<size_t>(corner_vertex_id)),
              MeshletVertexRuntimeInfo {});
            polygon_push(polygon, static_cast<size_t>(cv));
            record_site(static_cast<size_t>(cv), static_cast<size_t>(corner_vertex_id));
            break;
          }
          else
          {
            append_strand_corner_if_needed();
          }
        }

        boundary_he = kin_del.nextOnComponentBoundaryId(boundary_he);
        current_ref_opt.reset();
      }

      if (exited_fail)
      {
        break;
      }
      if (closed_walk_at_seed)
      {
        exited_closed_loop = true;
        break;
      }
      if (!found_next_segment)
      {
        KINDS_ERROR("Closing mesh: Delaunay boundary walk from end of ordered_seg["
          << current_segment_index << "] did not reach any strip's start_crossing (nor close the loop to ordered_seg["
          << seed_index << "]).");
        exited_fail = true;
        break;
      }
    }

    if (exited_closed_loop)
    {
      for (size_t idx : walk_segment_order)
      {
        if (cap_start_hits[idx] != 1)
        {
          KINDS_ERROR("Closing mesh: ordered_seg[" << idx
                                                   << "] start: expected exactly one incoming handoff (Delaunay or "
                                                      "Voronoi-vertex join), got "
                                                   << cap_start_hits[idx] << ".");
          exited_fail = true;
        }
        if (cap_end_handed_off[idx] != 1)
        {
          KINDS_ERROR("Closing mesh: ordered_seg[" << idx
                                                   << "] end: expected exactly one outgoing handoff (Delaunay or "
                                                      "Voronoi-vertex join), got "
                                                   << cap_end_handed_off[idx] << ".");
          exited_fail = true;
        }
      }
    }
    else if (!exited_fail && !walk_segment_order.empty())
    {
      KINDS_ERROR("Closing mesh: polygon walk for seed ordered_seg["
        << seed_index
        << "] ended without closing the loop (each strip end must hand off via a Delaunay boundary walk or a Voronoi "
           "circumcenter join matching dual Voronoi vertex ids to another strip's strand end or seed's strand-side "
           "vertex).");
      exited_fail = true;
    }

    {
      const MeshingData* end_seg = ordered_segments[last_walk_segment_index];
      const size_t ve_end = end_seg->closing_voronoi_edge_id;
      std::ostringstream head;
      head << "Ended polygon walk at ordered segment " << last_walk_segment_index << " of Voronoi edge ";
      if (ve_end == static_cast<size_t>(-1))
      {
        head << "?";
      }
      else
      {
        head << ve_end;
      }
      head << ".";
      if (exited_closed_loop && walk_closure_crossing_ref.has_value())
      {
        const glm::dvec3 mesh_pos = getCrossingCoordsInMeshSpace(walk_closure_crossing_ref.value(), t);
        const glm::dvec3 delaunay_pos = getCrossingCoordsInDelaunaySpace(kin_del, walk_closure_crossing_ref.value(), t);
        head << " Loop closed at crossing " << formatCrossingIntersectionForLog(kin_del, walk_closure_crossing_ref)
             << ".";
        if (std::isfinite(delaunay_pos.x) && std::isfinite(delaunay_pos.y))
        {
          head << " Delaunay 2D=(" << delaunay_pos.x << "," << delaunay_pos.y << ")";
        }
        else
        {
          head << " Delaunay 2D=(unavailable)";
        }
        if (std::isfinite(mesh_pos.x) && std::isfinite(mesh_pos.y))
        {
          head << "; mesh 2D=(" << mesh_pos.x << "," << mesh_pos.y << ").";
        }
        else
        {
          head << "; mesh 2D=(non-finite).";
        }
      }
      else if (exited_closed_loop)
      {
        head << " Loop closed at seed mesh_start (Voronoi circumcenter join; no Delaunay closure crossing).";
      }
      else if (!exited_fail)
      {
        head << " (walk did not close to seed start crossing).";
      }
      KINDS_DEBUG(head.str());
    }

    if (polygon.size() >= 3 && exited_closed_loop && !exited_fail)
    {
      const auto consecutive_delaunay_duplicate = [&](size_t a_id, size_t b_id)
      {
        const auto a_it = outline_vertex_info.find(a_id);
        const auto b_it = outline_vertex_info.find(b_id);
        if (a_it == outline_vertex_info.end() || b_it == outline_vertex_info.end())
        {
          return false;
        }
        const glm::dvec2& a = a_it->second.delaunay_xy;
        const glm::dvec2& b = b_it->second.delaunay_xy;
        if (!std::isfinite(a.x) || !std::isfinite(a.y) || !std::isfinite(b.x) || !std::isfinite(b.y))
        {
          return false;
        }
        const double scale = std::max({ 1.0, glm::length(a), glm::length(b) });
        return glm::distance(a, b) <= 1e-12 * scale;
      };

      // Filter only consecutive duplicates in the cyclic Delaunay-space outline. Deliberately do not use a set:
      // the same point appearing at non-adjacent positions can encode meaningful polygon topology.
      std::vector<size_t> filtered_polygon;
      filtered_polygon.reserve(polygon.size());
      for (size_t vertex_id : polygon)
      {
        if (filtered_polygon.empty() || !consecutive_delaunay_duplicate(filtered_polygon.back(), vertex_id))
        {
          filtered_polygon.push_back(vertex_id);
        }
      }
      // The outline is cyclic, so its last and first entries are consecutive as well.
      if (filtered_polygon.size() > 1
        && consecutive_delaunay_duplicate(filtered_polygon.back(), filtered_polygon.front()))
      {
        filtered_polygon.pop_back();
      }
      if (filtered_polygon.size() != polygon.size())
      {
        KINDS_DEBUG("Closing cap outline removed " << (polygon.size() - filtered_polygon.size())
                                                   << " consecutive Delaunay-space duplicate(s) for strand "
                                                   << strand_id << " at t=" << t);
      }
      polygon = std::move(filtered_polygon);
      if (polygon.size() < 3)
      {
        KINDS_WARNING("Closing cap outline has fewer than three vertices after consecutive Delaunay-space duplicate "
                      "filtering for strand "
          << strand_id << " at t=" << t << "; skipping triangulation.");
        continue;
      }

      // Closing-cap topology is defined by the Delaunay-space outline. Store that outline as the triangulation plane
      // before orientation so splitPolygonAtRepeatedVertices(), collinearity removal, and ear clipping all use the same
      // coordinates even when the emitted mesh vertices have been transformed to object space.
      for (size_t vertex_id : polygon)
      {
        const auto info_it = outline_vertex_info.find(vertex_id);
        if (info_it == outline_vertex_info.end() || !std::isfinite(info_it->second.delaunay_xy.x)
          || !std::isfinite(info_it->second.delaunay_xy.y))
        {
          std::ostringstream oss;
          oss << "Closing cap outline has no finite Delaunay-space coordinates for mesh vertex " << vertex_id
              << " (strand=" << strand_id << ", t=" << t << ").";
          throw std::runtime_error(oss.str());
        }
        mesh.setProfilePlanePosition(vertex_id, info_it->second.delaunay_xy);
      }

      double area2 = 0.0;
      for (size_t i = 0; i < polygon.size(); ++i)
      {
        const glm::dvec2 p0 = mesh.triangulationPlaneXY(polygon[i]);
        const glm::dvec2 p1 = mesh.triangulationPlaneXY(polygon[(i + 1) % polygon.size()]);
        area2 += p0.x * p1.y - p1.x * p0.y;
      }
      if (area2 < 0.0)
      {
        std::reverse(polygon.begin() + 1, polygon.end());
      }
      result.polygons.push_back(std::move(polygon));
    }
  }

  const auto append_outline_vertex
    = [](std::ostringstream& out, const CapOutlineVertexInfo& info, const glm::dvec3& position, bool include_z)
  {
    out << "(\"" << info.source << "\", ";
    if (std::string(info.source) == "intersection")
    {
      out << info.primary_id << ", " << info.secondary_id << ", ";
    }
    else
    {
      out << info.primary_id << ", ";
    }
    out << position.x << ", " << position.y;
    if (include_z)
    {
      out << ", " << position.z;
    }
    out << ")";
  };

  for (size_t outline_index = 0; outline_index < result.polygons.size(); ++outline_index)
  {
    const std::vector<size_t>& outline = result.polygons[outline_index];
    std::ostringstream delaunay_log;
    std::ostringstream mesh_log;
    delaunay_log << "Closing cap outline " << outline_index << " Delaunay space: [";
    mesh_log << "Closing cap outline " << outline_index << " mesh space: [";
    for (size_t i = 0; i < outline.size(); ++i)
    {
      if (i != 0)
      {
        delaunay_log << ", ";
        mesh_log << ", ";
      }
      const size_t vertex_id = outline[i];
      const auto info_it = outline_vertex_info.find(vertex_id);
      const CapOutlineVertexInfo unknown {};
      const CapOutlineVertexInfo& info = info_it != outline_vertex_info.end() ? info_it->second : unknown;
      append_outline_vertex(
        delaunay_log, info, glm::dvec3(info.delaunay_xy, std::numeric_limits<double>::quiet_NaN()), false);

      glm::dvec3 mesh_position(std::numeric_limits<double>::quiet_NaN());
      if (vertex_id < mesh.getVertices().size())
      {
        mesh_position = mesh.getVertices()[vertex_id];
      }
      append_outline_vertex(mesh_log, info, mesh_position, true);
    }
    delaunay_log << "]";
    mesh_log << "]";
    KINDS_DEBUG(delaunay_log.str());
    KINDS_DEBUG(mesh_log.str());
  }

  return result;
}

void kinDS::SegmentBuilder::closingMeshLogUnmatchedOrderedSegments(
  size_t strand_id, double t, const std::vector<MeshingData*>& ordered_segments, const std::vector<bool>& segment_used)
{
  size_t unmatched_count = 0;
  for (size_t i = 0; i < segment_used.size(); ++i)
  {
    if (!segment_used[i])
    {
      ++unmatched_count;
    }
  }
  if (unmatched_count > 0)
  {
    KINDS_WARNING("Closing mesh: " << unmatched_count << " of " << segment_used.size()
                                   << " ordered Voronoi segment(s) not reached by polygon tracing (strand " << strand_id
                                   << " t=" << t << ").");
    for (size_t i = 0; i < segment_used.size(); ++i)
    {
      if (!segment_used[i])
      {
        const MeshingData* s = ordered_segments[i];
        KINDS_WARNING("Closing mesh: unmatched ordered segment ["
          << i << "] incident_edge=" << s->closing_incident_edge_index << " voronoi_edge="
          << (s->closing_voronoi_edge_id == static_cast<size_t>(-1) ? -1 : static_cast<int>(s->closing_voronoi_edge_id))
          << ".");
      }
    }
  }
}

void kinDS::SegmentBuilder::triangulateSimplePolygon(VoronoiMesh& mesh, const std::vector<size_t>& polygon,
  const std::string& metadata, int material_id, bool orient_upwards, std::optional<double> occurrence_time,
  std::optional<size_t> runtime_branch_id)
{
  constexpr double eps = 1e-12;
  const auto cross2 = [](const glm::dvec2& u, const glm::dvec2& v) { return u.x * v.y - u.y * v.x; };

  std::vector<size_t> vertices;
  vertices.reserve(polygon.size());
  for (size_t vertex_id : polygon)
  {
    if (vertex_id >= mesh.getVertices().size())
    {
      throw std::runtime_error("triangulateSimplePolygon: polygon vertex index out of range.");
    }
    if (vertices.empty() || !triangulationPlaneXYEqual(mesh, vertices.back(), vertex_id))
    {
      vertices.push_back(vertex_id);
    }
  }
  if (vertices.size() > 1 && triangulationPlaneXYEqual(mesh, vertices.front(), vertices.back()))
  {
    vertices.pop_back();
  }
  if (vertices.size() < 3)
  {
    return;
  }

  const std::vector<std::vector<size_t>> split_polygons = splitPolygonAtRepeatedVertices(mesh, vertices);
  if (split_polygons.empty())
  {
    return;
  }
  if (split_polygons.size() > 1)
  {
    std::vector<std::pair<std::string, std::vector<size_t>>> debug_rings;
    debug_rings.reserve(split_polygons.size() + 1);
    debug_rings.emplace_back("original", vertices);
    for (size_t i = 0; i < split_polygons.size(); ++i)
    {
      debug_rings.emplace_back("sub_" + std::to_string(i), split_polygons[i]);
    }
    if (shouldDumpErrorFiles())
    {
      writeTriangulateSimplePolygonDebugTxt(kin_del, mesh,
        makeTriangulateSimplePolygonDebugPath(kin_del, mesh, "SPLIT", ".txt", occurrence_time, runtime_branch_id),
        "SPLIT", debug_rings, occurrence_time, runtime_branch_id);
    }

    for (const std::vector<size_t>& sub_polygon : split_polygons)
    {
      triangulateSimplePolygon(
        mesh, sub_polygon, metadata, material_id, orient_upwards, occurrence_time, runtime_branch_id);
    }
    return;
  }
  vertices = split_polygons.front();

  // Snapshot plane coordinates once. Ear-clip must not re-query mesh positions (object-space verts) mid-pass.
  std::vector<glm::dvec2> plane(vertices.size());
  for (size_t i = 0; i < vertices.size(); ++i)
  {
    plane[i] = mesh.triangulationPlaneXY(vertices[i]);
    if (!std::isfinite(plane[i].x) || !std::isfinite(plane[i].y))
    {
      throw std::runtime_error("triangulateSimplePolygon: non-finite triangulation-plane coordinate.");
    }
  }

  auto cross_at = [&](size_t prev_i, size_t current_i, size_t next_i)
  { return cross2(plane[current_i] - plane[prev_i], plane[next_i] - plane[current_i]); };

  bool removed_collinear = true;
  while (removed_collinear && vertices.size() > 3)
  {
    removed_collinear = false;
    for (size_t i = 0; i < vertices.size(); ++i)
    {
      const size_t prev_i = (i + vertices.size() - 1) % vertices.size();
      const size_t next_i = (i + 1) % vertices.size();
      if (std::abs(cross_at(prev_i, i, next_i)) <= eps)
      {
        vertices.erase(vertices.begin() + static_cast<std::ptrdiff_t>(i));
        plane.erase(plane.begin() + static_cast<std::ptrdiff_t>(i));
        removed_collinear = true;
        break;
      }
    }
  }
  if (vertices.size() < 3)
  {
    return;
  }

  auto signed_area2 = [&]()
  {
    double area2 = 0.0;
    for (size_t i = 0; i < plane.size(); ++i)
    {
      const glm::dvec2& p0 = plane[i];
      const glm::dvec2& p1 = plane[(i + 1) % plane.size()];
      area2 += p0.x * p1.y - p1.x * p0.y;
    }
    return area2;
  };

  const double area2 = signed_area2();
  if (std::abs(area2) <= eps)
  {
    throw std::runtime_error("triangulateSimplePolygon: polygon area is degenerate.");
  }
  if (area2 < 0.0)
  {
    std::reverse(vertices.begin(), vertices.end());
    std::reverse(plane.begin(), plane.end());
  }

  auto point_in_triangle = [&](size_t p_i, size_t a_i, size_t b_i, size_t c_i)
  {
    const glm::dvec2& p = plane[p_i];
    const glm::dvec2& a = plane[a_i];
    const glm::dvec2& b = plane[b_i];
    const glm::dvec2& c = plane[c_i];
    const double ab = cross2(b - a, p - a);
    const double bc = cross2(c - b, p - b);
    const double ca = cross2(a - c, p - c);
    return ab >= -eps && bc >= -eps && ca >= -eps;
  };

  auto emit_triangle = [&](size_t i0, size_t i1, size_t i2)
  {
    if (orient_upwards)
    {
      addMeshletTriangle(mesh, vertices[i0], vertices[i1], vertices[i2], metadata, material_id);
    }
    else
    {
      addMeshletTriangle(mesh, vertices[i0], vertices[i2], vertices[i1], metadata, material_id);
    }
  };

  auto polygon_is_convex = [&]()
  {
    bool any_pos = false;
    bool any_neg = false;
    for (size_t i = 0; i < plane.size(); ++i)
    {
      const double cr = cross_at((i + plane.size() - 1) % plane.size(), i, (i + 1) % plane.size());
      if (cr > eps)
      {
        any_pos = true;
      }
      if (cr < -eps)
      {
        any_neg = true;
      }
    }
    return !(any_pos && any_neg);
  };

  while (vertices.size() > 3)
  {
    bool clipped_ear = false;
    for (size_t i = 0; i < vertices.size(); ++i)
    {
      const size_t prev_i = (i + vertices.size() - 1) % vertices.size();
      const size_t next_i = (i + 1) % vertices.size();
      if (cross_at(prev_i, i, next_i) <= eps)
      {
        continue;
      }

      bool contains_other_vertex = false;
      for (size_t candidate_i = 0; candidate_i < vertices.size(); ++candidate_i)
      {
        if (candidate_i == prev_i || candidate_i == i || candidate_i == next_i)
        {
          continue;
        }
        if (point_in_triangle(candidate_i, prev_i, i, next_i))
        {
          contains_other_vertex = true;
          break;
        }
      }
      if (contains_other_vertex)
      {
        continue;
      }

      emit_triangle(prev_i, i, next_i);
      vertices.erase(vertices.begin() + static_cast<std::ptrdiff_t>(i));
      plane.erase(plane.begin() + static_cast<std::ptrdiff_t>(i));
      clipped_ear = true;
      break;
    }

    if (!clipped_ear)
    {
      // Convex rings (including the common radius cell quads) can always fan-triangulate.
      if (polygon_is_convex())
      {
        for (size_t i = 1; i + 1 < vertices.size(); ++i)
        {
          emit_triangle(0, i, i + 1);
        }
        return;
      }

      if (shouldDumpErrorFiles())
      {
        writeTriangulateSimplePolygonDebugTxt(kin_del, mesh,
          makeTriangulateSimplePolygonDebugPath(kin_del, mesh, "FAIL", ".txt", occurrence_time, runtime_branch_id),
          "FAIL", { { "fail", vertices } }, occurrence_time, runtime_branch_id);
        writeTriangulateSimplePolygonFailSvg(kin_del, mesh, vertices, occurrence_time, runtime_branch_id);
      }
      KINDS_ERROR(
        "triangulateSimplePolygon: failed to find an ear; polygon may be non-simple. Omitting polygon from mesh.");
      return;
    }
  }

  emit_triangle(0, 1, 2);
}

void kinDS::SegmentBuilder::fanTriangulateConvexPolygon(VoronoiMesh& mesh, const std::vector<size_t>& polygon,
  const std::string& metadata, int material_id, bool orient_upwards)
{
  std::vector<size_t> vertices;
  vertices.reserve(polygon.size());
  for (size_t vertex_id : polygon)
  {
    if (vertex_id >= mesh.getVertices().size())
    {
      throw std::runtime_error("fanTriangulateConvexPolygon: polygon vertex index out of range.");
    }
    if (vertices.empty() || vertices.back() != vertex_id)
    {
      vertices.push_back(vertex_id);
    }
  }
  if (vertices.size() > 1 && vertices.front() == vertices.back())
  {
    vertices.pop_back();
  }
  if (vertices.size() < 3)
  {
    return;
  }

  for (size_t i = 1; i + 1 < vertices.size(); ++i)
  {
    if (orient_upwards)
    {
      addMeshletTriangle(mesh, vertices[0], vertices[i], vertices[i + 1], metadata, material_id);
    }
    else
    {
      addMeshletTriangle(mesh, vertices[0], vertices[i + 1], vertices[i], metadata, material_id);
    }
  }
}

void kinDS::SegmentBuilder::closingMeshTriangulatePolygons(
  VoronoiMesh& mesh, const std::vector<std::vector<size_t>>& polygons, double t, size_t strand_id)
{
  const std::string face_metadata = composeClosingMeshFaceMetadata(t, strand_id);
  const size_t runtime_branch_id = kin_del.getRuntimeBranchIdForStrand(strand_id);
  for (const auto& polygon : polygons)
  {
    triangulateSimplePolygon(mesh, polygon, face_metadata, RegularMeshletMaterialId, true, t, runtime_branch_id);
  }
}

void kinDS::SegmentBuilder::createClosingCapForStrand(size_t strand_id, double t)
{
  if (kin_del.isDummyBoundary(strand_id) || !kin_del.isStrandLiveInGraph(strand_id))
  {
    return;
  }

  if (strand_id >= kin_del.component_data.component_map.size() || strand_id >= strand_to_segment_indices.size())
  {
    return;
  }

  if (strand_to_segment_indices[strand_id].empty())
  {
    return;
  }

  const size_t component_index = kin_del.component_data.component_map[strand_id];
  if (component_index >= kin_del.component_data.component_boundaries.size())
  {
    return;
  }

  const auto& boundary_polygon = kin_del.component_data.component_boundaries[component_index][0];
  const auto& centroid = kin_del.component_data.component_centroids[component_index];

  const size_t closing_mesh_index = createClosingMesh(strand_id, t, boundary_polygon, centroid);
  if (closing_mesh_index >= segment_mesh_pairs.size())
  {
    return;
  }

  MeshStructure::SegmentMeshPair& segment_mesh_pair = segment_mesh_pairs[closing_mesh_index];
  segment_mesh_pair.segment_index0 = static_cast<int>(strand_to_segment_indices[strand_id].back());
  segment_mesh_pair.segment_index1 = -1;
}

void kinDS::SegmentBuilder::createClosingCapsForInputBranchFinishingAtSection(double t, size_t input_branch_id)
{
  const size_t section = static_cast<size_t>(t);
  const auto& branch_strands = kin_del.getStrandTree().getStrandsByBranch(section, input_branch_id);
  for (size_t strand_id : branch_strands)
  {
    createClosingCapForStrand(strand_id, t);
  }
}

void kinDS::SegmentBuilder::createClosingCapsForInputBranchesFinishingAtSection(double t)
{
  for (size_t input_branch_id : kin_del.inputBranchesFinishingAtSection(t))
  {
    createClosingCapsForInputBranchFinishingAtSection(t, input_branch_id);
  }
}

void kinDS::SegmentBuilder::finishIncidentStripMeshesForStrandAtSection(size_t strand_id, double t)
{
  if (kin_del.isDummyBoundary(strand_id) || !kin_del.isStrandLiveInGraph(strand_id))
  {
    return;
  }

  if (strand_id >= kin_del.component_data.component_map.size())
  {
    return;
  }

  auto& graph = kin_del.getGraph();
  const size_t component_index = kin_del.component_data.component_map[strand_id];
  auto& boundary_polygon = kin_del.component_data.component_boundaries[component_index][0];

  for (HalfEdgeDelaunayGraph::IncidentEdgeIterator it = graph.incidentEdgesBegin(strand_id),
                                                   end = graph.incidentEdgesEnd(strand_id);
    it != end; ++it)
  {
    finishMesh(*it, t, boundary_polygon, BoundaryEventType::Section, BoundarySegmentAction::SegmentCompleted);
  }
}

size_t kinDS::SegmentBuilder::createClosingMesh(size_t strand_id, double t,
  const std::vector<BoundaryPoint>& boundary_polygon, const glm::dvec2& centroid,
  std::vector<BoundaryIntersectionInterval>* traced_boundary_intervals)
{
  KINDS_DEBUG("createClosingMesh strand " << strand_id << " t=" << t);
  const size_t num_incident_edges = closingMeshCountStrandIncidentEdges(strand_id);
  const size_t mesh_pair_index_before = segment_mesh_pairs.size();

  if (diagnostics)
  {
    std::ostringstream oss;
    oss << "num_incident_edges=" << num_incident_edges << " boundary_polygon_size=" << boundary_polygon.size()
        << " is_dummy_boundary=" << (kin_del.isDummyBoundary(strand_id) ? "true" : "false")
        << " is_live_in_graph=" << (kin_del.isStrandLiveInGraph(strand_id) ? "true" : "false");
    strandInitDiagnosticLogLine("create_closing_mesh_begin", strand_id, t, oss.str().c_str());
  }

  if (traced_boundary_intervals != nullptr)
  {
    traced_boundary_intervals->clear();
  }

  segment_mesh_pairs.push_back(MeshStructure::SegmentMeshPair {});

  VoronoiMesh mesh;
  configureMeshletStorage(mesh);
  std::vector<size_t> mesh_vertex_ids;
  mesh_vertex_ids.reserve(32);
  auto track_vertex = [&](const glm::dvec3& pos, std::optional<size_t> vv, const std::string& metadata,
                        const MeshletVertexRuntimeInfo& runtime_info) -> int
  {
    const int id
      = closingMeshAppendVertex(mesh, boundary_polygon, centroid, strand_id, t, pos, false, vv, metadata, runtime_info);
    mesh_vertex_ids.push_back(static_cast<size_t>(id));
    return id;
  };

  auto& graph = kin_del.getGraph();
  // Raw segments keyed by undirected Voronoi (dual) edge id; one bucket per incident finite edge.
  std::map<size_t, std::vector<MeshingData>> segments_by_voronoi_edge;
  std::list<MeshingData> closing_segments;

  // Trace every Voronoi edge incident to this strand and keep only portions that lie inside.
  // Walk stored `voronoi_edge_intersections` (re-oriented and with recomputed crossing positions) so inside/outside
  // toggles stay aligned with `CrossingData` without re-running `computeCrossedHalfEdges`.
  int incident_edge_index = -1;
  for (HalfEdgeDelaunayGraph::IncidentEdgeIterator it = graph.incidentEdgesBegin(strand_id),
                                                   end = graph.incidentEdgesEnd(strand_id);
    it != end; ++it)
  {
    ++incident_edge_index;
    const size_t incident_he = *it;
    const size_t voronoi_edge_id = incident_he / 2;
    std::vector<MeshingData> bucket = closingMeshExtractRawSegmentsForVoronoiEdge(
      strand_id, t, mesh, boundary_polygon, centroid, incident_edge_index, incident_he, track_vertex);
    segments_by_voronoi_edge[voronoi_edge_id] = bucket;
    closing_segments.insert(
      closing_segments.end(), std::make_move_iterator(bucket.begin()), std::make_move_iterator(bucket.end()));
  }

  {
    size_t extraction_index = 0;
    for (const auto& seg : closing_segments)
    {
      const int ve_log
        = seg.closing_voronoi_edge_id == static_cast<size_t>(-1) ? -1 : static_cast<int>(seg.closing_voronoi_edge_id);
      KINDS_DEBUG("Closing mesh extraction, segment "
        << extraction_index << " (incident_edge=" << seg.closing_incident_edge_index << " voronoi_edge=" << ve_log
        << " mesh_v " << seg.mesh_start_vertex_id << "->" << seg.mesh_end_vertex_id << "):");
      logExtractionClosingMeshStripEndpoint(
        kin_del, t, extraction_index, seg.start_crossing, seg.closing_voronoi_edge_id);
      logExtractionClosingMeshStripEndpoint(
        kin_del, t, extraction_index, seg.end_crossing, seg.closing_voronoi_edge_id);
      ++extraction_index;
    }
  }

  ClosingMeshOrderedIndex index_data = closingMeshBuildOrderedIndex(closing_segments);

  closingMeshValidateOrderedSegmentGeometry(t, mesh, index_data.ordered_segments);

  ClosingMeshPolygonsTraceResult trace
    = closingMeshTraceCapPolygons(strand_id, t, num_incident_edges, mesh, boundary_polygon, centroid,
      std::move(mesh_vertex_ids), index_data.ordered_segments, index_data.start_crossing_to_segment);
  if (traced_boundary_intervals != nullptr)
  {
    *traced_boundary_intervals = trace.traced_boundary_intervals;
  }

  closingMeshLogUnmatchedOrderedSegments(strand_id, t, index_data.ordered_segments, trace.segment_used);

  const size_t cap_vertices_before_tri = mesh.getVertexCount();
  closingMeshTriangulatePolygons(mesh, trace.polygons, t, strand_id);
  const size_t cap_vertices_after_tri = mesh.getVertexCount();
  const size_t cap_triangles_after_tri = mesh.getTriangleCount();

  if (diagnostics)
  {
    size_t polygon_vertex_sum = 0;
    for (const auto& polygon : trace.polygons)
    {
      polygon_vertex_sum += polygon.size();
    }
    std::ostringstream oss;
    oss << "mesh_pair_index=" << mesh_pair_index_before << " raw_segments=" << closing_segments.size()
        << " ordered_segments=" << index_data.ordered_segments.size() << " polygons=" << trace.polygons.size()
        << " polygon_vertex_sum=" << polygon_vertex_sum << " cap_verts_before_tri=" << cap_vertices_before_tri
        << " cap_verts_after_tri=" << cap_vertices_after_tri << " cap_tris=" << cap_triangles_after_tri;
    strandInitDiagnosticLogLine("create_closing_mesh_after_trace", strand_id, t, oss.str().c_str());
    if (cap_vertices_after_tri == 0 && cap_triangles_after_tri == 0)
    {
      strandInitDiagnosticLogLine(
        "create_closing_mesh_empty", strand_id, t, "closing cap has no vertices/triangles after trace+triangulation");
    }
  }

  const size_t index
    = registerMeshletWithSuffix(std::move(mesh), std::string("_strand") + std::to_string(strand_id), t);
  segment_mesh_pair_last_left_and_right_vertex.emplace_back();
  segment_mesh_pair_last_left_and_right_vertex.back() = std::move(closing_segments);

  if (diagnostics)
  {
    std::ostringstream oss;
    oss << "meshlet_index=" << index << " mesh_pair_index=" << mesh_pair_index_before
        << " pair_index_matches_meshlet=" << (index == mesh_pair_index_before ? "true" : "false");
    strandInitDiagnosticLogLine("create_closing_mesh_end", strand_id, t, oss.str().c_str());
  }

  return index;
}

void kinDS::SegmentBuilder::accumulateSegmentProperties()
{
  // Iterate through all pairs and accumulate properties
  for (size_t pair_id = 0; pair_id < segment_mesh_pairs.size(); ++pair_id)
  {
    auto& pair = segment_mesh_pairs[pair_id];
    if (pair.segment_index0 != -1)
    {
      // make sure there is space left
      if (segment_properties[pair.segment_index0].neighbor_count >= MeshStructure::SegmentProperties::MAX_NEIGHBORS)
      {
        KINDS_ERROR("Exceeded maximum number of neighbors for segment: "
          << segment_properties[pair.segment_index0].neighbor_count
          << " >= " << MeshStructure::SegmentProperties::MAX_NEIGHBORS);
        // throw std::runtime_error("Exceeded maximum number of neighbors for segment.");
      }
      else
      {
        segment_properties[pair.segment_index0]
          .mesh_pair_indices[segment_properties[pair.segment_index0].neighbor_count] = pair_id; // add mesh pair index
        segment_properties[pair.segment_index0].neighbor_indices[segment_properties[pair.segment_index0].neighbor_count]
          = pair.segment_index1; // add neighbor
        segment_properties[pair.segment_index0].neighbor_count++;
      }
    }

    if (pair.segment_index1 != -1)
    {
      // make sure there is space left
      if (segment_properties[pair.segment_index1].neighbor_count >= MeshStructure::SegmentProperties::MAX_NEIGHBORS)
      {
        KINDS_ERROR("Exceeded maximum number of neighbors for segment: "
          << segment_properties[pair.segment_index1].neighbor_count
          << " >= " << MeshStructure::SegmentProperties::MAX_NEIGHBORS);
        // throw std::runtime_error("Exceeded maximum number of neighbors for segment.");
      }
      else
      {
        segment_properties[pair.segment_index1]
          .mesh_pair_indices[segment_properties[pair.segment_index1].neighbor_count] = pair_id; // add mesh pair index
        segment_properties[pair.segment_index1].neighbor_indices[segment_properties[pair.segment_index1].neighbor_count]
          = pair.segment_index0; // add neighbor
        segment_properties[pair.segment_index1].neighbor_count++;
      }
    }
  }
}

void kinDS::SegmentBuilder::splitComponent(
  size_t component_id, const std::vector<std::vector<size_t>>& new_components, double t)
{
  if (new_components.empty())
  {
    return;
  }

  const std::vector<size_t> pre_split_parent_strands = kin_del.component_data.components[component_id];

  // Update component data
  std::vector<size_t> component_ids(new_components.size(), -1);

  component_ids[0] = component_id;
  size_t new_component_id = kin_del.component_data.components.size();
  for (size_t i = 1; i < new_components.size(); i++)
  {
    component_ids[i] = new_component_id;
    new_component_id++;
  }

  size_t new_size = kin_del.component_data.components.size() + new_components.size() - 1;
  kin_del.component_data.components.resize(new_size);
  kin_del.component_data.component_boundaries.resize(new_size);
  kin_del.component_data.component_centroids.resize(new_size);
  kin_del.component_data.component_last_updated.resize(new_size);

  // Assign strand membership first and note the pending split (writes branches.txt) before boundary
  // extraction, which can throw on malformed split-off pieces.
  for (size_t i = 0; i < new_components.size(); i++)
  {
    size_t cid = component_ids[i];

    for (size_t v : new_components[i])
    {
      kin_del.component_data.component_map[v] = cid;
    }

    kin_del.component_data.components[cid] = new_components[i];

    glm::dvec2 provisional_centroid { 0.0, 0.0 };
    double provisional_weight = 0.0;
    for (size_t strand_id : new_components[i])
    {
      if (kin_del.isDummyBoundary(strand_id))
      {
        continue;
      }
      provisional_centroid += kin_del.getPointAt(strand_id, t, false, false);
      provisional_weight += 1.0;
    }
    if (provisional_weight > 0.0)
    {
      provisional_centroid /= provisional_weight;
    }
    kin_del.component_data.component_centroids[cid] = provisional_centroid;
    kin_del.component_data.component_last_updated[cid] = t;
  }

  kin_del.notePendingBranchSplit(component_id, t, pre_split_parent_strands, new_components, component_ids);

  std::vector<bool> he_visited(kin_del.getGraph().halfEdgeSlotCount(), false);
  for (size_t i = 0; i < new_components.size(); i++)
  {
    size_t cid = component_ids[i];
    kin_del.component_data.component_boundaries[cid]
      = kin_del.extractComponentBoundaries(new_components[i], t, he_visited, false, false);
    if (!kin_del.component_data.component_boundaries[cid].empty()
      && !kin_del.component_data.component_boundaries[cid][0].empty())
    {
      kin_del.component_data.component_centroids[cid]
        = polygonCentroid(kin_del.component_data.component_boundaries[cid][0]);
    }
  }

  // Refresh separation centroids now that polygon centroids are available, then enqueue separation.
  kin_del.refreshPendingSplitSeparationCentroids(component_id);
  kin_del.maybeScheduleSeparationOrApplyPendingSplit(component_id, t);
}

void SegmentBuilder::init()
{
  configureMeshletStorage(boundary_mesh);

  kin_del.registerEventCallbacks(
    section_callback_.get(), flip_callback_.get(), radius_callback_.get(), crossing_callback_.get());
  kin_del.registerSubdivisionEventCallback(subdivision_callback_.get());
  kin_del.registerSeparationEventCallback(separation_callback_.get());

  auto& graph = kin_del.getGraph();

  size_t strand_count = graph.getVertexCount();
  strand_to_segment_indices.resize(strand_count);
  half_edge_index_to_segment_mesh_pair_index.resize(graph.halfEdgeSlotCount(), -1);
  corner_to_cutoff_mesh_indices.resize(graph.halfEdgeSlotCount(), -1);

  // Initialize the strand geometries at the kinetic bootstrap section.
  double t = static_cast<double>(kin_del.getStartSection());

  // We need a ruled surface for each half-edge in the graph except those having the infinite vertex as
  // origin
  size_t half_edge_count = graph.halfEdgeSlotCount();

  // initialize segment mesh properties for each strand
  for (size_t strand_id = 0; strand_id < strand_count; ++strand_id)
  {
    size_t new_segment_id = segment_properties.size();
    MeshStructure::SegmentProperties properties;
    segment_properties.push_back(properties);
    strand_to_segment_indices[strand_id].push_back(new_segment_id);

    size_t component_index = kin_del.component_data.component_map[strand_id];
    const auto& component = kin_del.component_data.components[component_index];
    const auto& boundary_polygon = kin_del.component_data.component_boundaries[component_index][0];
    const auto& centroid = kin_del.component_data.component_centroids[component_index];

    if (diagnostics)
    {
      std::ostringstream oss;
      oss << "segment_id=" << new_segment_id << " component_id=" << component_index
          << " component_size=" << component.size() << " boundary_polygon_size=" << boundary_polygon.size()
          << " mesh_cap_at_start=" << (mesh_cap_at_start ? "true" : "false");
      strandInitDiagnosticLogLine("init_cap_begin", strand_id, t, oss.str().c_str());
    }

    if (mesh_cap_at_start)
    {
      // create a closing mesh
      size_t closing_mesh_index = createClosingMesh(strand_id, t, boundary_polygon, centroid);
      MeshStructure::SegmentMeshPair& segment_mesh_pair = segment_mesh_pairs[new_segment_id];
      segment_mesh_pair.segment_index0 = -1;
      segment_mesh_pair.segment_index1 = strand_to_segment_indices[strand_id].back();

      if (diagnostics)
      {
        std::ostringstream oss;
        oss << "closing_meshlet_index=" << closing_mesh_index << " mesh_pair_slot=" << new_segment_id
            << " segment_index1=" << segment_mesh_pair.segment_index1;
        if (closing_mesh_index < meshes.size())
        {
          oss << " cap_verts=" << meshes[closing_mesh_index].getVertexCount()
              << " cap_tris=" << meshes[closing_mesh_index].getTriangleCount();
        }
        strandInitDiagnosticLogLine("init_cap_end", strand_id, t, oss.str().c_str());
      }
    }
    else if (diagnostics)
    {
      strandInitDiagnosticLogLine("init_cap_skipped", strand_id, t, "mesh_cap_at_start=false");
    }
  }

  // now go through all half-edges and create a segment mesh pair
  for (size_t i = 0; i < half_edge_count; i += 2)
  {
    startNewMesh(i, t);
  }

  // Create boundary interval meshes at bootstrap time from precomputed crossing data.
  for (size_t i = 0; i < half_edge_count; i += 2)
  {
    if (!kin_del.isOnComponentBoundary(i))
    {
      continue;
    }

    const size_t d_edge_id = i / 2;
    const auto& d_intersections = kin_del.getCrossingData().delaunay_edge_intersections[d_edge_id];
    if (d_intersections.empty())
    {
      continue;
    }

    const std::vector<KineticDelaunay::CrossingData::EdgeIntersectionRef> refs
      = getBoundaryIntersectionsInBoundaryOrder(d_edge_id);

    {
      const size_t first_cell
        = determineVoronoiCellForBoundaryIntersectionInterval(d_edge_id, std::nullopt, refs.front());
      startNewMeshFromIntersections(
        first_cell, t, std::nullopt, refs.front(), false, BoundaryEventType::Init, BoundarySegmentAction::NewSegment);
    }

    for (size_t k = 0; k + 1 < refs.size(); ++k)
    {
      const size_t mid_cell = determineVoronoiCellForBoundaryIntersectionInterval(d_edge_id, refs[k], refs[k + 1]);
      if (mid_cell == static_cast<size_t>(-1))
      {
        continue;
      }
      startNewMeshFromIntersections(
        mid_cell, t, refs[k], refs[k + 1], false, BoundaryEventType::Init, BoundarySegmentAction::NewSegment);
    }

    {
      const size_t last_cell
        = determineVoronoiCellForBoundaryIntersectionInterval(d_edge_id, refs.back(), std::nullopt);
      startNewMeshFromIntersections(
        last_cell, t, refs.back(), std::nullopt, false, BoundaryEventType::Init, BoundarySegmentAction::NewSegment);
    }
  }

  // initialize boundary mesh
  boundary_mesh_last_left_and_right_vertex.resize(half_edge_count, std::make_pair(-1, -1));
  half_edge_to_boundary_vertex_index.resize(half_edge_count, -1);
  const size_t init_section = static_cast<size_t>(t);
  const auto& branches_at_section = kin_del.getStrandTree().getStrandBranchesByHeight(init_section);
  for (size_t input_branch_id = 0; input_branch_id < branches_at_section.size(); ++input_branch_id)
  {
    if (branches_at_section[input_branch_id].empty())
    {
      if (diagnostics)
      {
        std::ostringstream oss;
        oss << "input_branch_id=" << input_branch_id << " reason=empty_branch_at_section";
        strandInitDiagnosticLogLine("init_boundary_branch_skip", 0, t, oss.str().c_str());
      }
      continue;
    }

    bool any_real_strand = false;
    bool contains_strand_0 = false;
    for (size_t strand_id : branches_at_section[input_branch_id])
    {
      if (strand_id == 0)
      {
        contains_strand_0 = true;
      }
      if (!kin_del.isDummyBoundary(strand_id))
      {
        any_real_strand = true;
      }
    }
    if (!any_real_strand)
    {
      if (diagnostics)
      {
        std::ostringstream oss;
        oss << "input_branch_id=" << input_branch_id << " reason=only_dummy_strands"
            << " contains_strand_0=" << (contains_strand_0 ? "true" : "false");
        strandInitDiagnosticLogLine("init_boundary_branch_skip", 0, t, oss.str().c_str());
      }
      continue;
    }

    if (diagnostics)
    {
      std::ostringstream oss;
      oss << "input_branch_id=" << input_branch_id << " invert_orientation=false offset=-0.01"
          << " contains_strand_0=" << (contains_strand_0 ? "true" : "false");
      strandInitDiagnosticLogLine("init_boundary_branch_triangulate", 0, t, oss.str().c_str());
    }
    addDelaunayTriangulationToBoundaryMesh(t, input_branch_id, false, -0.01);
  }

  for (size_t input_branch_id : kin_del.inputBranchesFinishingAtSection(t))
  {
    const auto& branch_strands = kin_del.getStrandTree().getStrandsByBranch(static_cast<size_t>(t), input_branch_id);
    for (size_t strand_id : branch_strands)
    {
      finishIncidentStripMeshesForStrandAtSection(strand_id, t);
    }
    createClosingCapsForInputBranchFinishingAtSection(t, input_branch_id);

    if (diagnostics)
    {
      std::ostringstream oss;
      oss << "input_branch_id=" << input_branch_id << " invert_orientation=true offset=0.01 reason=branch_finishing";
      strandInitDiagnosticLogLine("init_boundary_branch_triangulate", 0, t, oss.str().c_str());
    }
    addDelaunayTriangulationToBoundaryMesh(t, input_branch_id, true, 0.01);
  }

  logStrandInitDiagnosticsSummary(t);
}

void SegmentBuilder::finalize(double t)
{
  updateBoundaries(t, collectLiveComponentIndices());

  // Finalize the segments by finishing all meshes
  auto& graph = kin_del.getGraph();

  for (size_t he_id : graph.liveDelaunayEdges())
  {
    auto vertex = graph.halfEdge(he_id).origin;

    // fall back for infinite vertices
    if (vertex == -1)
    {
      vertex = graph.destination(he_id);
    }

    size_t component_index = kin_del.component_data.component_map[vertex];
    auto& boundary_points = kin_del.component_data.component_boundaries[component_index][0];

    finishMesh(he_id, t, boundary_points);
  }

  // Finalize boundary-interval meshes once more at the final time for all boundary Delaunay-edge sections.
  for (size_t he_id : graph.liveDelaunayEdges())
  {
    if (!kin_del.isOnComponentBoundary(he_id))
    {
      continue;
    }
    const size_t d_edge_id = he_id / 2;
    const auto& d_intersections = kin_del.getCrossingData().delaunay_edge_intersections[d_edge_id];
    if (d_intersections.empty())
    {
      continue;
    }

    const std::vector<KineticDelaunay::CrossingData::EdgeIntersectionRef> refs
      = getBoundaryIntersectionsInBoundaryOrder(d_edge_id);

    {
      const size_t first_cell
        = determineVoronoiCellForBoundaryIntersectionInterval(d_edge_id, std::nullopt, refs.front());
      finishMeshFromIntersections(
        first_cell, t, std::nullopt, refs.front(), BoundaryEventType::Section, BoundarySegmentAction::SegmentCompleted);
    }
    for (size_t k = 0; k + 1 < refs.size(); ++k)
    {
      const size_t mid_cell = determineVoronoiCellForBoundaryIntersectionInterval(d_edge_id, refs[k], refs[k + 1]);
      if (mid_cell == static_cast<size_t>(-1))
      {
        continue;
      }
      finishMeshFromIntersections(
        mid_cell, t, refs[k], refs[k + 1], BoundaryEventType::Section, BoundarySegmentAction::SegmentCompleted);
    }
    {
      const size_t last_cell
        = determineVoronoiCellForBoundaryIntersectionInterval(d_edge_id, refs.back(), std::nullopt);
      finishMeshFromIntersections(
        last_cell, t, refs.back(), std::nullopt, BoundaryEventType::Section, BoundarySegmentAction::SegmentCompleted);
    }
  }

  // Radius complementary splits were queued when the target mid-interval meshlet was identified
  // (often with only seed verts 0/1 and no triangles yet). Apply them now that strips have faces,
  // before flexible resolution / dedupe / glue postprocess.
  flushPendingRadiusComplementarySplits();

  // Natural branch endings at finalize time (including the tree top). Section events run on
  // [start, end), so finishers at the exclusive end_section are never capped in the section
  // callback — always seal them here. mesh_cap_at_end only seals strands truncated by a premature --end.
  createClosingCapsForInputBranchesFinishingAtSection(t);
  for (size_t input_branch_id : kin_del.inputBranchesFinishingAtSection(t))
  {
    addDelaunayTriangulationToBoundaryMesh(t, input_branch_id, true, 0.01);
  }

  if (mesh_cap_at_end)
  {
    const size_t section = static_cast<size_t>(t);
    for (size_t strand_id = 0; strand_id < graph.getVertexCount(); ++strand_id)
    {
      if (kin_del.isDummyBoundary(strand_id) || !kin_del.isStrandLiveInGraph(strand_id))
      {
        continue;
      }
      // Skip natural finishers (already capped above); only seal strands that continue past --end.
      if (kin_del.getStrandTree().getSupportPoints(strand_id).size() <= section + 1)
      {
        continue;
      }
      createClosingCapForStrand(strand_id, t);
    }
  }

  accumulateSegmentProperties();

  resolveAllIntersectionFlexibleVertices("finalize intersection mesh");
  collapseDegreeTwoFlexibleVerticesOnIntersectionMeshes();

  // Exact duplicate verts break boundary walks / glue matching; dedupe before normals and before any
  // later triangle-merge / segment-assembly passes (0 tolerance). Then drop index-degenerate faces.
  size_t degenerate_triangles_removed = 0;
  for (auto& meshlet : meshes)
  {
    meshlet.mergeDuplicateVertices(0.0);
    degenerate_triangles_removed += meshlet.removeDegenerateTriangles();
  }
  for (auto& meshlet : intersection_meshes)
  {
    meshlet.mergeDuplicateVertices(0.0);
    degenerate_triangles_removed += meshlet.removeDegenerateTriangles();
  }

  // compute normals
  for (auto& meshlet : meshes)
  {
    meshlet.ensureFaceMetadataSize();
    meshlet.computeNormals(NormalMode::PerTriangleCorner);
    meshlet.validateUVLayout("finalized interior meshlet");
  }
  for (auto& meshlet : intersection_meshes)
  {
    meshlet.ensureFaceMetadataSize();
    meshlet.computeNormals(NormalMode::PerTriangleCorner);
    meshlet.validateUVLayout("finalized intersection meshlet");
  }

  auto remap1 = boundary_mesh.mergeDuplicateVertices();
  degenerate_triangles_removed += boundary_mesh.removeDegenerateTriangles();
  auto remap2 = boundary_mesh.removeIsolatedVertices();
  boundary_mesh.ensureFaceMetadataSize();
  boundary_mesh.computeNormals(NormalMode::PerTriangleCorner);
  boundary_mesh.validateUVLayout("finalized boundary mesh");

  KINDS_INFO("finalize: removed " << degenerate_triangles_removed
                                  << " degenerate triangle(s) across interior, intersection, and boundary meshes");

  // Update boundary vertex to strand id mapping
  std::vector<size_t> new_boundary_vertex_to_strand_id;
  new_boundary_vertex_to_strand_id.resize(boundary_mesh.getVertices().size());
  for (size_t old_index = 0; old_index < boundary_vertex_to_strand_id.size(); ++old_index)
  {
    size_t new_index = remap2[remap1[old_index]];
    if (new_index != size_t(-1))
    {
      new_boundary_vertex_to_strand_id[new_index] = boundary_vertex_to_strand_id[old_index];
    }
  }

  boundary_vertex_to_strand_id.swap(new_boundary_vertex_to_strand_id);

  if (store_mesh_metadata && validate_mesh_vertex_sources)
  {
    std::vector<VoronoiMesh*> all_meshlet_ptrs;
    all_meshlet_ptrs.reserve(meshes.size() + intersection_meshes.size());
    std::vector<std::string> labels;
    labels.reserve(all_meshlet_ptrs.capacity());

    for (size_t mesh_index = 0; mesh_index < meshes.size(); ++mesh_index)
    {
      all_meshlet_ptrs.push_back(&meshes[mesh_index]);
      labels.push_back("mesh_" + std::to_string(mesh_index));
    }
    Validator::validateAndReportInteriorMeshletUvHeights(
      all_meshlet_ptrs, "segment builder finalize interior", uv_height_factor, &labels);

    std::vector<VoronoiMesh*> bark_meshlet_ptrs;
    std::vector<std::string> bark_labels;
    bark_meshlet_ptrs.reserve(intersection_meshes.size() + 1);
    bark_labels.reserve(intersection_meshes.size() + 1);
    for (size_t mesh_index = 0; mesh_index < intersection_meshes.size(); ++mesh_index)
    {
      bark_meshlet_ptrs.push_back(&intersection_meshes[mesh_index]);
      bark_labels.push_back("intersection_mesh_" + std::to_string(mesh_index));
    }
    bark_meshlet_ptrs.push_back(&boundary_mesh);
    bark_labels.push_back("boundary_mesh");
    Validator::validateAndReportBarkMeshletUvHeights(
      bark_meshlet_ptrs, "segment builder finalize bark", uv_height_factor, &bark_labels);

    for (size_t mesh_index = 0; mesh_index < intersection_meshes.size(); ++mesh_index)
    {
      all_meshlet_ptrs.push_back(&intersection_meshes[mesh_index]);
      labels.push_back("intersection_mesh_" + std::to_string(mesh_index));
    }

    Validator::validateAndReportMeshVertexSources(all_meshlet_ptrs, "segment builder finalize", &labels);
  }

  if (!create_transformed_mesh)
  {
    applyUntransformedMeshViewTransform();
  }

  finalized = true; // Set the finalized flag to true
}

std::vector<VoronoiMesh> kinDS::SegmentBuilder::extractMeshes() const { return meshes; }

namespace
{
struct SegmentAssemblyContrib
{
  VoronoiMesh mesh;
  int neighbor = -1;
  /// Index into @c meshes (interior) or @c intersection_meshes (boundary); spaces overlap without prefix.
  size_t source_index = 0;
  bool is_boundary = false;

  std::string contributorLabel() const { return (is_boundary ? "b" : "i") + std::to_string(source_index); }
};

struct Vec3ExactHash
{
  std::size_t operator()(const glm::dvec3& v) const noexcept
  {
    std::size_t h = std::hash<double> {}(v.x);
    h ^= std::hash<double> {}(v.y) + 0x9e3779b97f4a7c15ull + (h << 6) + (h >> 2);
    h ^= std::hash<double> {}(v.z) + 0x9e3779b97f4a7c15ull + (h << 6) + (h >> 2);
    return h;
  }
};

struct Vec3ExactEq
{
  bool operator()(const glm::dvec3& a, const glm::dvec3& b) const noexcept
  {
    return a.x == b.x && a.y == b.y && a.z == b.z;
  }
};

std::string markPostprocessSplitMetadata(const std::string& metadata)
{
  if (metadata.empty())
  {
    return SegmentBuilder::MetadataBuilder().addBool("postprocess_split_flexible_glue_edge", true).build();
  }
  return SegmentBuilder::MetadataBuilder::fromObject(metadata)
    .addBool("postprocess_split_flexible_glue_edge", true)
    .build();
}

struct BoundaryMatch
{
  size_t mesh_index = 0;
  size_t vertex_id = 0;
};

struct MeshBoundaryInfo
{
  std::unordered_map<size_t, size_t> successor;
  std::vector<std::vector<size_t>> cycles;
  std::unordered_set<size_t> boundary_verts;
  /// Per boundary vertex: all position matches on other contrib meshes.
  std::unordered_map<size_t, std::vector<BoundaryMatch>> matches;
};

/// Directed boundary edges as successor map: each boundary vertex -> next along the boundary.
std::unordered_map<size_t, size_t> buildBoundarySuccessorMap(const VoronoiMesh& mesh)
{
  const auto& tris = mesh.getTriangles();
  std::map<std::pair<size_t, size_t>, int> directed;
  for (size_t t = 0; t + 2 < tris.size(); t += 3)
  {
    const size_t a = tris[t];
    const size_t b = tris[t + 1];
    const size_t c = tris[t + 2];
    ++directed[{ a, b }];
    ++directed[{ b, c }];
    ++directed[{ c, a }];
  }

  std::unordered_map<size_t, size_t> successor;
  for (const auto& [edge, count] : directed)
  {
    if (count <= 0)
    {
      continue;
    }
    const auto rev = directed.find({ edge.second, edge.first });
    if (rev != directed.end() && rev->second > 0)
    {
      continue; // interior (or opposite-oriented pair)
    }
    // Prefer a single outgoing edge per vertex; last write wins if non-manifold.
    successor[edge.first] = edge.second;
  }
  return successor;
}

std::vector<std::vector<size_t>> traceBoundaryCycles(const VoronoiMesh& mesh)
{
  const auto successor = buildBoundarySuccessorMap(mesh);
  std::unordered_set<size_t> visited;
  std::vector<std::vector<size_t>> cycles;
  for (const auto& [start, first_next] : successor)
  {
    (void)first_next;
    if (visited.count(start))
    {
      continue;
    }
    std::vector<size_t> cycle;
    size_t cur = start;
    while (!visited.count(cur))
    {
      visited.insert(cur);
      cycle.push_back(cur);
      const auto it = successor.find(cur);
      if (it == successor.end())
      {
        cycle.clear();
        break;
      }
      cur = it->second;
      if (cur == start)
      {
        break;
      }
    }
    if (cycle.size() >= 2 && cur == start)
    {
      cycles.push_back(std::move(cycle));
    }
  }
  return cycles;
}

std::optional<size_t> findTriangleWithDirectedEdge(const VoronoiMesh& mesh, size_t v0, size_t v1, size_t& local_edge)
{
  const auto& tris = mesh.getTriangles();
  const size_t n_tri = tris.size() / 3;
  for (size_t t = 0; t < n_tri; ++t)
  {
    for (size_t e = 0; e < 3; ++e)
    {
      if (tris[3 * t + e] == v0 && tris[3 * t + ((e + 1) % 3)] == v1)
      {
        local_edge = e;
        return t;
      }
    }
  }
  return std::nullopt;
}

/// Remove a top-level JSON object field @p key (bool / string / array / number) if present.
std::string removeJsonObjectField(std::string metadata, const char* key)
{
  if (key == nullptr || metadata.empty())
  {
    return metadata;
  }
  const std::string needle = std::string("\"") + key + "\":";
  const size_t pos = metadata.find(needle);
  if (pos == std::string::npos)
  {
    return metadata;
  }

  size_t erase_begin = pos;
  if (erase_begin > 0)
  {
    size_t p = erase_begin;
    do
    {
      --p;
    } while (p > 0 && (metadata[p] == ' ' || metadata[p] == '\t' || metadata[p] == '\n' || metadata[p] == '\r'));
    if (metadata[p] == ',')
    {
      erase_begin = p;
    }
  }

  size_t value_start = pos + needle.size();
  while (value_start < metadata.size()
    && (metadata[value_start] == ' ' || metadata[value_start] == '\t' || metadata[value_start] == '\n'
      || metadata[value_start] == '\r'))
  {
    ++value_start;
  }
  if (value_start >= metadata.size())
  {
    return metadata;
  }

  size_t value_end = value_start;
  if (metadata[value_start] == '"')
  {
    value_end = metadata.find('"', value_start + 1);
    if (value_end == std::string::npos)
    {
      return metadata;
    }
    ++value_end;
  }
  else if (metadata[value_start] == '[')
  {
    int depth = 0;
    for (size_t i = value_start; i < metadata.size(); ++i)
    {
      if (metadata[i] == '[')
      {
        ++depth;
      }
      else if (metadata[i] == ']')
      {
        --depth;
        if (depth == 0)
        {
          value_end = i + 1;
          break;
        }
      }
    }
    if (depth != 0)
    {
      return metadata;
    }
  }
  else if (metadata.compare(value_start, 4, "true") == 0)
  {
    value_end = value_start + 4;
  }
  else if (metadata.compare(value_start, 5, "false") == 0)
  {
    value_end = value_start + 5;
  }
  else
  {
    while (value_end < metadata.size() && metadata[value_end] != ',' && metadata[value_end] != '}'
      && metadata[value_end] != ' ' && metadata[value_end] != '\t' && metadata[value_end] != '\n')
    {
      ++value_end;
    }
  }

  size_t erase_end = value_end;
  if (erase_begin == pos)
  {
    size_t t = value_end;
    while (
      t < metadata.size() && (metadata[t] == ' ' || metadata[t] == '\t' || metadata[t] == '\n' || metadata[t] == '\r'))
    {
      ++t;
    }
    if (t < metadata.size() && metadata[t] == ',')
    {
      erase_end = t + 1;
    }
  }

  metadata.erase(erase_begin, erase_end - erase_begin);
  return metadata;
}

std::string formatContributorsJson(std::vector<std::string> labels)
{
  std::sort(labels.begin(), labels.end());
  labels.erase(std::unique(labels.begin(), labels.end()), labels.end());
  std::ostringstream oss;
  oss << "[";
  for (size_t i = 0; i < labels.size(); ++i)
  {
    if (i > 0)
    {
      oss << ",";
    }
    oss << jsonStringLiteral(labels[i]);
  }
  oss << "]";
  return oss.str();
}

/// Set @c matched and @c contributors on a boundary vertex (no-op when metadata storage is disabled).
void setVertexMatchedAndContributors(
  VoronoiMesh& mesh, size_t vertex_index, const std::string& matched_raw_json, const std::string& contributors_raw_json)
{
  if (!mesh.storeMetadata() || vertex_index >= mesh.getVertices().size())
  {
    return;
  }
  std::string existing = "{}";
  if (vertex_index < mesh.getVertexMetadata().size() && !mesh.getVertexMetadata()[vertex_index].empty())
  {
    existing = mesh.getVertexMetadata()[vertex_index];
  }
  existing = removeJsonObjectField(std::move(existing), "matched");
  existing = removeJsonObjectField(std::move(existing), "contributors");
  mesh.setVertexMetadata(vertex_index,
    SegmentBuilder::MetadataBuilder::fromObject(existing)
      .addRaw("matched", matched_raw_json)
      .addRaw("contributors", contributors_raw_json)
      .build());
}

void setVertexMatchedMetadata(
  VoronoiMesh& mesh, size_t vertex_index, bool matched, const std::string& contributors_raw_json)
{
  setVertexMatchedAndContributors(mesh, vertex_index, matched ? "true" : "false", contributors_raw_json);
}

// Temporary debug filter: only emit glue-align warnings/info for this meshing segment id.
constexpr size_t kGlueAlignDebugMeshletIndex = 1962;

bool glueAlignLogEnabled(size_t meshlet_index) { return meshlet_index == kGlueAlignDebugMeshletIndex; }

/// Insert @p position on boundary edge (v0—v1), splitting the unique adjacent triangle. Marks flexible when requested.
/// @p meshlet_index matches PerSegment OBJ / failed_meshlet index (meshing segment id).
size_t insertVertexOnBoundaryEdge(VoronoiMesh& mesh, size_t v0, size_t v1, const glm::dvec3& position,
  double kinetic_time, bool mark_flexible, size_t meshlet_index, int physics_segment_id,
  const std::string& contributors_raw_json)
{
  size_t local_edge = 0;
  std::optional<size_t> tri_opt = findTriangleWithDirectedEdge(mesh, v0, v1, local_edge);
  if (!tri_opt.has_value())
  {
    tri_opt = findTriangleWithDirectedEdge(mesh, v1, v0, local_edge);
  }
  if (!tri_opt.has_value())
  {
    if (glueAlignLogEnabled(meshlet_index))
    {
      const auto& verts = mesh.getVertices();
      std::ostringstream oss;
      oss << "insertVertexOnBoundaryEdge: no triangle for boundary edge on meshlet_index=" << meshlet_index
          << " physics_segment_id=" << physics_segment_id << " vids=" << v0 << "<->" << v1 << " insert_pos=("
          << position.x << ", " << position.y << ", " << position.z << ")";
      if (v0 < verts.size())
      {
        oss << " edge_v0_pos=(" << verts[v0].x << ", " << verts[v0].y << ", " << verts[v0].z << ")";
      }
      else
      {
        oss << " edge_v0_pos=<oob>";
      }
      if (v1 < verts.size())
      {
        oss << " edge_v1_pos=(" << verts[v1].x << ", " << verts[v1].y << ", " << verts[v1].z << ")";
      }
      else
      {
        oss << " edge_v1_pos=<oob>";
      }
      KINDS_WARNING(oss.str());
    }
    return static_cast<size_t>(-1);
  }

  const size_t tri = tri_opt.value();
  auto& tris = mesh.getTriangles();
  auto& uv_indices = mesh.getUVIndices();
  const size_t a = tris[3 * tri + local_edge];
  const size_t b = tris[3 * tri + ((local_edge + 1) % 3)];
  const size_t c = tris[3 * tri + ((local_edge + 2) % 3)];

  const size_t new_vid = mesh.addVertex(position);
  mesh.setVertexKineticTime(new_vid, kinetic_time);
  if (mark_flexible)
  {
    mesh.setVertexFlexible(new_vid, true);
  }
  if (mesh.storeMetadata())
  {
    setVertexMatchedAndContributors(mesh, new_vid, "\"by_propagation\"", contributors_raw_json);
  }
  if (const auto uv0 = mesh.vertexSemanticUv(a); uv0.has_value())
  {
    if (const auto uv1 = mesh.vertexSemanticUv(b); uv1.has_value())
    {
      const double len = glm::length(glm::dvec2(mesh.getVertices()[b]) - glm::dvec2(mesh.getVertices()[a]));
      double s = 0.5;
      if (len > 1e-18)
      {
        s = glm::length(glm::dvec2(position) - glm::dvec2(mesh.getVertices()[a])) / len;
        s = std::clamp(s, 0.0, 1.0);
      }
      mesh.setVertexSemanticUv(new_vid, uv0.value() * (1.0 - s) + uv1.value() * s);
    }
  }

  size_t uv_a = std::numeric_limits<size_t>::max();
  size_t uv_b = std::numeric_limits<size_t>::max();
  size_t uv_c = std::numeric_limits<size_t>::max();
  size_t uv_f = std::numeric_limits<size_t>::max();
  const bool has_uvs = uv_indices.size() == tris.size();
  if (has_uvs)
  {
    uv_a = uv_indices[3 * tri + local_edge];
    uv_b = uv_indices[3 * tri + ((local_edge + 1) % 3)];
    uv_c = uv_indices[3 * tri + ((local_edge + 2) % 3)];
    if (uv_a < mesh.getUVs().size() && uv_b < mesh.getUVs().size())
    {
      const glm::dvec3& ua = mesh.getUVs()[uv_a];
      const glm::dvec3& ub = mesh.getUVs()[uv_b];
      uv_f = mesh.addUV(0.5 * (ua + ub));
    }
  }

  const int material_id = tri < mesh.getMaterialIDs().size() ? mesh.getMaterialIDs()[tri] : -1;
  std::string face_meta = "{}";
  if (mesh.storeMetadata() && tri < mesh.getFaceMetadata().size())
  {
    face_meta = mesh.getFaceMetadata()[tri];
  }
  if (mesh.storeMetadata())
  {
    face_meta = markPostprocessSplitMetadata(face_meta);
  }

  // Replace original triangle (a,b,c) with (a,F,c); append (F,b,c).
  tris[3 * tri + local_edge] = a;
  tris[3 * tri + ((local_edge + 1) % 3)] = new_vid;
  tris[3 * tri + ((local_edge + 2) % 3)] = c;
  if (has_uvs)
  {
    uv_indices[3 * tri + local_edge] = uv_a;
    uv_indices[3 * tri + ((local_edge + 1) % 3)] = uv_f;
    uv_indices[3 * tri + ((local_edge + 2) % 3)] = uv_c;
  }
  mesh.addTriangle(new_vid, b, c, uv_f, uv_b, uv_c, material_id, face_meta);
  return new_vid;
}

MeshBoundaryInfo buildMeshBoundaryInfo(const VoronoiMesh& mesh)
{
  MeshBoundaryInfo info;
  info.successor = buildBoundarySuccessorMap(mesh);
  info.cycles = traceBoundaryCycles(mesh);
  for (const auto& cycle : info.cycles)
  {
    for (size_t vid : cycle)
    {
      info.boundary_verts.insert(vid);
    }
  }
  for (const auto& [v, next] : info.successor)
  {
    (void)next;
    info.boundary_verts.insert(v);
  }
  return info;
}

void initBoundaryMatchedMetadataFalse(
  VoronoiMesh& mesh, const MeshBoundaryInfo& info, const SegmentAssemblyContrib& contrib)
{
  if (!mesh.storeMetadata())
  {
    return;
  }
  const std::string contributors = formatContributorsJson({ contrib.contributorLabel() });
  for (size_t vid : info.boundary_verts)
  {
    setVertexMatchedMetadata(mesh, vid, false, contributors);
  }
}

bool sameContributor(const SegmentAssemblyContrib& a, const SegmentAssemblyContrib& b)
{
  return a.is_boundary == b.is_boundary && a.source_index == b.source_index;
}

void buildBoundaryMultiMatches(std::vector<MeshBoundaryInfo>& boundaries, std::vector<SegmentAssemblyContrib>& contribs)
{
  const size_t n = contribs.size();
  std::unordered_map<glm::dvec3, std::vector<BoundaryMatch>, Vec3ExactHash, Vec3ExactEq> by_position;
  for (size_t i = 0; i < n; ++i)
  {
    const auto& verts = contribs[i].mesh.getVertices();
    for (size_t vid : boundaries[i].boundary_verts)
    {
      if (vid >= verts.size())
      {
        continue;
      }
      by_position[verts[vid]].push_back(BoundaryMatch { i, vid });
    }
  }

  for (size_t i = 0; i < n; ++i)
  {
    const auto& verts = contribs[i].mesh.getVertices();
    for (size_t vid : boundaries[i].boundary_verts)
    {
      if (vid >= verts.size())
      {
        continue;
      }
      const auto it = by_position.find(verts[vid]);
      if (it == by_position.end())
      {
        continue;
      }
      std::vector<BoundaryMatch>& out = boundaries[i].matches[vid];
      const SegmentAssemblyContrib& self = contribs[i];
      for (const BoundaryMatch& m : it->second)
      {
        if (m.mesh_index >= contribs.size())
        {
          continue;
        }
        const SegmentAssemblyContrib& other = contribs[m.mesh_index];
        // Same contributor only (i/b + source index). i9 and b9 are distinct meshes.
        if (sameContributor(other, self))
        {
          continue;
        }
        out.push_back(m);
      }
      if (!out.empty())
      {
        std::vector<std::string> labels;
        labels.reserve(out.size() + 1);
        labels.push_back(self.contributorLabel());
        for (const BoundaryMatch& m : out)
        {
          if (m.mesh_index < contribs.size())
          {
            labels.push_back(contribs[m.mesh_index].contributorLabel());
          }
        }
        setVertexMatchedMetadata(contribs[i].mesh, vid, true, formatContributorsJson(std::move(labels)));
      }
    }
  }
}

std::optional<size_t> partnerVertexOnContributor(const MeshBoundaryInfo& info, size_t vid,
  const SegmentAssemblyContrib& partner, const std::vector<SegmentAssemblyContrib>& contribs)
{
  const auto it = info.matches.find(vid);
  if (it == info.matches.end())
  {
    return std::nullopt;
  }
  for (const BoundaryMatch& m : it->second)
  {
    if (m.mesh_index >= contribs.size())
    {
      continue;
    }
    if (sameContributor(contribs[m.mesh_index], partner))
    {
      return m.vertex_id;
    }
  }
  return std::nullopt;
}

std::optional<size_t> partnerVertexOnMesh(
  const MeshBoundaryInfo& info, size_t vid, size_t partner_mesh, const std::vector<SegmentAssemblyContrib>& contribs)
{
  if (partner_mesh >= contribs.size())
  {
    return std::nullopt;
  }
  return partnerVertexOnContributor(info, vid, contribs[partner_mesh], contribs);
}

std::vector<size_t> cycleArcPathInclusive(const std::vector<size_t>& cycle, size_t c0, size_t c1)
{
  const size_t n = cycle.size();
  std::vector<size_t> path;
  if (n == 0)
  {
    return path;
  }
  size_t cur = c0;
  path.push_back(cycle[cur]);
  while (cur != c1)
  {
    cur = (cur + 1) % n;
    path.push_back(cycle[cur]);
    if (path.size() > n + 1)
    {
      return {};
    }
  }
  return path;
}

/// For each consecutive edge on @p path, record the unique adjacent triangle (either winding).
bool pathEdgeTriangles(const VoronoiMesh& mesh, const std::vector<size_t>& path, std::vector<size_t>& out_triangles)
{
  out_triangles.clear();
  if (path.size() < 2)
  {
    return false;
  }
  out_triangles.reserve(path.size() - 1);
  for (size_t i = 0; i + 1 < path.size(); ++i)
  {
    size_t local_edge = 0;
    std::optional<size_t> tri = findTriangleWithDirectedEdge(mesh, path[i], path[i + 1], local_edge);
    if (!tri.has_value())
    {
      tri = findTriangleWithDirectedEdge(mesh, path[i + 1], path[i], local_edge);
    }
    if (!tri.has_value())
    {
      out_triangles.clear();
      return false;
    }
    out_triangles.push_back(tri.value());
  }
  return true;
}

bool cycleArcPathWithTriangles(const VoronoiMesh& mesh, const std::vector<size_t>& cycle, size_t c0, size_t c1,
  std::vector<size_t>& out_path, std::vector<size_t>& out_triangles)
{
  out_path = cycleArcPathInclusive(cycle, c0, c1);
  if (out_path.size() < 2)
  {
    out_triangles.clear();
    return false;
  }
  return pathEdgeTriangles(mesh, out_path, out_triangles);
}

/// True if @p vid has at least one position match on any other contrib mesh.
bool boundaryVertexHasAnyMatch(const MeshBoundaryInfo& info, size_t vid)
{
  const auto it = info.matches.find(vid);
  return it != info.matches.end() && !it->second.empty();
}

/// True if every strict interior of @p path is unmatched to every contrib (endpoints unchecked).
bool pathInteriorsUnmatchedOnly(const MeshBoundaryInfo& info, const std::vector<size_t>& path)
{
  for (size_t t = 1; t + 1 < path.size(); ++t)
  {
    if (boundaryVertexHasAnyMatch(info, path[t]))
    {
      return false;
    }
  }
  return true;
}

/// Locate @p vid on @p cycle; returns false if absent.
bool findCycleIndex(const std::vector<size_t>& cycle, size_t vid, size_t& out_index)
{
  for (size_t i = 0; i < cycle.size(); ++i)
  {
    if (cycle[i] == vid)
    {
      out_index = i;
      return true;
    }
  }
  return false;
}

/// Partner-cycle arc between two vids.
/// Priority: (1) arcs whose interiors are unmatched to every contrib — prefer opposite order
/// (@p right → @p left), else same order; among ties pick length closest to @p preferred_path_len.
/// (2) Only if no unmatched-only arc exists (and @p require_unmatched_interiors is false): the
/// shorter of the two cycle arcs, with opposite winning length ties.
/// When @p require_unmatched_interiors is true, step (2) is skipped (failure if no unmatched-only arc).
bool findPartnerBoundaryArcPreferOpposite(const MeshBoundaryInfo& partner_info, size_t left_partner_vid,
  size_t right_partner_vid, size_t preferred_path_len, std::vector<size_t>& out_path,
  bool require_unmatched_interiors)
{
  out_path.clear();
  if (left_partner_vid == right_partner_vid)
  {
    return false;
  }

  auto collect_oriented_paths = [&](bool opposite_order, bool unmatched_interiors_only)
  {
    std::vector<std::vector<size_t>> paths;
    for (const auto& cycle : partner_info.cycles)
    {
      if (cycle.size() < 2)
      {
        continue;
      }
      size_t left_ci = 0;
      size_t right_ci = 0;
      if (!findCycleIndex(cycle, left_partner_vid, left_ci) || !findCycleIndex(cycle, right_partner_vid, right_ci))
      {
        continue;
      }
      const size_t c0 = opposite_order ? right_ci : left_ci;
      const size_t c1 = opposite_order ? left_ci : right_ci;
      std::vector<size_t> path = cycleArcPathInclusive(cycle, c0, c1);
      if (path.size() < 2)
      {
        continue;
      }
      if (unmatched_interiors_only && !pathInteriorsUnmatchedOnly(partner_info, path))
      {
        continue;
      }
      paths.push_back(std::move(path));
    }
    return paths;
  };

  auto pick_closest_len = [&](std::vector<std::vector<size_t>>& paths) -> bool
  {
    if (paths.empty())
    {
      return false;
    }
    size_t best = 0;
    size_t best_dist = static_cast<size_t>(-1);
    for (size_t i = 0; i < paths.size(); ++i)
    {
      const size_t len = paths[i].size();
      const size_t d = len > preferred_path_len ? len - preferred_path_len : preferred_path_len - len;
      if (d < best_dist)
      {
        best_dist = d;
        best = i;
      }
    }
    out_path = std::move(paths[best]);
    return true;
  };

  auto pick_shortest_prefer_opposite = [&](std::vector<std::vector<size_t>>& opposite_paths,
                                         std::vector<std::vector<size_t>>& same_paths) -> bool
  {
    auto best_of = [](std::vector<std::vector<size_t>>& paths) -> size_t
    {
      size_t best = 0;
      for (size_t i = 1; i < paths.size(); ++i)
      {
        if (paths[i].size() < paths[best].size())
        {
          best = i;
        }
      }
      return best;
    };

    const bool has_opp = !opposite_paths.empty();
    const bool has_same = !same_paths.empty();
    if (!has_opp && !has_same)
    {
      return false;
    }
    if (has_opp && !has_same)
    {
      out_path = std::move(opposite_paths[best_of(opposite_paths)]);
      return true;
    }
    if (!has_opp && has_same)
    {
      out_path = std::move(same_paths[best_of(same_paths)]);
      return true;
    }
    const size_t oi = best_of(opposite_paths);
    const size_t si = best_of(same_paths);
    if (same_paths[si].size() < opposite_paths[oi].size())
    {
      out_path = std::move(same_paths[si]);
    }
    else
    {
      // Prefer opposite on a tie or when opposite is strictly shorter.
      out_path = std::move(opposite_paths[oi]);
    }
    return true;
  };

  // Unmatched-only arcs always take precedence over a shorter arc that has matched interiors.
  std::vector<std::vector<size_t>> opposite_unmatched
    = collect_oriented_paths(/*opposite_order=*/true, /*unmatched_interiors_only=*/true);
  if (pick_closest_len(opposite_unmatched))
  {
    return true;
  }
  std::vector<std::vector<size_t>> same_unmatched
    = collect_oriented_paths(/*opposite_order=*/false, /*unmatched_interiors_only=*/true);
  if (pick_closest_len(same_unmatched))
  {
    return true;
  }

  if (require_unmatched_interiors)
  {
    return false;
  }

  // Debug dump fallback: no unmatched-only arc — pick the shorter cycle way.
  std::vector<std::vector<size_t>> opposite_any
    = collect_oriented_paths(/*opposite_order=*/true, /*unmatched_interiors_only=*/false);
  std::vector<std::vector<size_t>> same_any
    = collect_oriented_paths(/*opposite_order=*/false, /*unmatched_interiors_only=*/false);
  return pick_shortest_prefer_opposite(opposite_any, same_any);
}

/// True if the forward inclusive cycle arc is strictly longer than the complementary arc.
/// Used to skip the long way around when identifying consecutive match pairs on a cycle.
bool isStrictlyLongerCycleArc(size_t inclusive_path_size, size_t cycle_len)
{
  if (inclusive_path_size < 2 || cycle_len < 2)
  {
    return false;
  }
  const size_t edges = inclusive_path_size - 1;
  return edges > cycle_len - edges;
}

/// Partner seam arc between the two partner vids that match a source consecutive pair.
/// Interiors on that arc must be unmatched to every mesh (not merely unmatched to @p source_mesh).
/// Uses the frozen @p partner_info cycles (not updated after splits).
///
/// Glue boundaries typically run in opposite directions, so the preferred arc is partner cycle
/// order @p right_partner_vid → @p left_partner_vid. Same order
/// (@p left_partner_vid → @p right_partner_vid) is only accepted if no opposite-order arc with
/// unmatched-only interiors exists. When several candidates exist in the chosen orientation, the
/// arc whose length is closest to @p preferred_path_len wins.
bool findReciprocalConsecutivePartnerPath(const VoronoiMesh& partner_mesh, const MeshBoundaryInfo& partner_info,
  size_t source_mesh, size_t left_partner_vid, size_t right_partner_vid, size_t preferred_path_len,
  const std::vector<SegmentAssemblyContrib>& contribs, std::vector<size_t>& out_path,
  std::vector<size_t>& out_triangles)
{
  (void)source_mesh;
  (void)contribs;
  out_path.clear();
  out_triangles.clear();
  if (!findPartnerBoundaryArcPreferOpposite(
        partner_info, left_partner_vid, right_partner_vid, preferred_path_len, out_path,
        /*require_unmatched_interiors=*/true))
  {
    return false;
  }
  return pathEdgeTriangles(partner_mesh, out_path, out_triangles);
}

std::string formatVec3(const glm::dvec3& v)
{
  std::ostringstream oss;
  oss << std::setprecision(17) << "(" << v.x << ", " << v.y << ", " << v.z << ")";
  return oss.str();
}

std::string formatContribLabel(const std::vector<SegmentAssemblyContrib>& contribs, size_t contrib_index)
{
  if (contrib_index >= contribs.size())
  {
    return "?(contrib=" + std::to_string(contrib_index) + ")";
  }
  return contribs[contrib_index].contributorLabel();
}

/// TXT-style seam dump: matched endpoints as (src, partner); then source interiors as (src, _),
/// then partner interiors as (_, partner). When @p mark_non_pair_matches is set, interior indices that
/// match some contrib (i.e. not unmatched) are written with a trailing '?'.
std::string formatGlueSeamPathScheme(const std::vector<size_t>& source_path, const std::vector<size_t>& partner_path,
  const MeshBoundaryInfo* source_info = nullptr, const MeshBoundaryInfo* partner_info = nullptr,
  bool mark_non_pair_matches = false)
{
  if (source_path.size() < 2 || partner_path.size() < 2)
  {
    return {};
  }
  auto format_index = [&](size_t vid, const MeshBoundaryInfo* info)
  {
    std::string s = std::to_string(vid);
    if (mark_non_pair_matches && info != nullptr && boundaryVertexHasAnyMatch(*info, vid))
    {
      s.push_back('?');
    }
    return s;
  };

  std::ostringstream oss;
  oss << "(" << source_path.front() << ", " << partner_path.front() << ")";
  for (size_t t = 1; t + 1 < source_path.size(); ++t)
  {
    oss << " - (" << format_index(source_path[t], source_info) << ", _)";
  }
  for (size_t t = 1; t + 1 < partner_path.size(); ++t)
  {
    oss << " - (_, " << format_index(partner_path[t], partner_info) << ")";
  }
  oss << " - (" << source_path.back() << ", " << partner_path.back() << ")";
  return oss.str();
}

std::filesystem::path glueAlignBoundaryVertsTxtPath(
  const KineticDelaunay& kin_del, size_t meshlet_index, int physics_segment_id)
{
  std::filesystem::path filepath = std::filesystem::path("glue_align")
    / ("boundary_verts_meshlet" + std::to_string(meshlet_index) + "_phys" + std::to_string(physics_segment_id)
      + ".txt");
  if (const std::optional<std::filesystem::path>& output_root = kin_del.getVisualDebugOutputRoot();
    output_root.has_value())
  {
    filepath = *output_root / filepath;
  }
  return filepath;
}

bool pathHasPosition(const VoronoiMesh& mesh, const std::vector<size_t>& path, const glm::dvec3& pos)
{
  const auto& verts = mesh.getVertices();
  for (const size_t vid : path)
  {
    if (vid < verts.size() && verts[vid] == pos)
    {
      return true;
    }
  }
  return false;
}

bool meshHasExactVertexPosition(const VoronoiMesh& mesh, const glm::dvec3& pos)
{
  const auto& verts = mesh.getVertices();
  for (const glm::dvec3& v : verts)
  {
    if (v == pos)
    {
      return true;
    }
  }
  return false;
}

void warnGlueMatchedPair(size_t meshlet_index, int physics_segment_id, const char* reason,
  const std::vector<SegmentAssemblyContrib>& contribs, size_t source_mesh, size_t partner_mesh, size_t left_vid,
  size_t right_vid, const glm::dvec3& left_pos, const glm::dvec3& right_pos, size_t left_partner_vid,
  size_t right_partner_vid, const glm::dvec3& left_partner_pos, const glm::dvec3& right_partner_pos,
  const std::string& path_scheme = {})
{
  if (!glueAlignLogEnabled(meshlet_index))
  {
    return;
  }
  std::ostringstream oss;
  oss << "glue align: " << reason << "\n"
      << "\n"
      << "  meshlet_index=" << meshlet_index << " physics_segment_id=" << physics_segment_id
      << " source=" << formatContribLabel(contribs, source_mesh)
      << " partner=" << formatContribLabel(contribs, partner_mesh) << "\n";
  if (!path_scheme.empty())
  {
    oss << "\n"
        << "  path: " << path_scheme << "\n";
  }
  oss << "\n"
      << "  matched0: source=" << formatContribLabel(contribs, source_mesh) << " vid=" << left_vid
      << " pos=" << formatVec3(left_pos) << " | partner=" << formatContribLabel(contribs, partner_mesh)
      << " vid=" << left_partner_vid << " pos=" << formatVec3(left_partner_pos) << "\n"
      << "\n"
      << "  matched1: source=" << formatContribLabel(contribs, source_mesh) << " vid=" << right_vid
      << " pos=" << formatVec3(right_pos) << " | partner=" << formatContribLabel(contribs, partner_mesh)
      << " vid=" << right_partner_vid << " pos=" << formatVec3(right_partner_pos) << "\n";
  KINDS_WARNING(oss.str());
}

/**
 * @brief Propagate unmatched boundary interiors both ways along a reciprocal seam.
 *
 * First inserts source-path interiors onto the partner path, then inserts partner-path interiors
 * that were not created by that first pass onto the source path.
 *
 * @param partner_path
 *   Inclusive vertex walk on the partner mesh between the two partner vids that match the
 *   source matched pair. Built as a forward arc on the partner's own boundary cycle between those
 *   two vids whose strict interiors are unmatched to every contrib. First and last entries are
 *   therefore the matched partner endpoints; interiors (if any) lie strictly between them.
 *   Endpoint order prefers the opposite of @p source_path along the partner cycle
 *   (@c right_partner → … → @c left_partner); same order is only used if that arc with
 *   unmatched-only interiors does not exist. This function may still reverse @p partner_path by
 *   endpoint coordinates before splitting so both walks share the same geometric start/end.
 *
 * @param partner_triangles
 *   One adjacent triangle id per consecutive edge of @p partner_path
 *   (@c size == partner_path.size() - 1).
 *
 * @param source_path
 *   Inclusive vertex walk on the source mesh along its boundary cycle from one consecutive
 *   partner-matched vertex to the next (@c left … right in source cycle order). First and last
 *   are those two matched vertices; entries @c [1, size-2] are the boundary verts strictly between
 *   them on that arc. Those interiors are not matched to this partner (by construction of
 *   consecutive same-partner matches); they may still match other contributors. A pure matched
 *   edge is @c size == 2 (no interiors) and is skipped by the caller.
 *
 * @param source_triangles
 *   One adjacent triangle id per consecutive edge of @p source_path
 *   (@c size == source_path.size() - 1).
 */
bool propagateUnmatchedInteriorsOntoPartner(VoronoiMesh& partner_mesh, std::vector<size_t>& partner_path,
  std::vector<size_t> partner_triangles, VoronoiMesh& source_mesh, std::vector<size_t>& source_path,
  std::vector<size_t>& source_triangles, size_t meshlet_index, int physics_segment_id,
  const std::string& contributors_raw_json, const std::vector<SegmentAssemblyContrib>& contribs,
  size_t source_mesh_index, size_t partner_mesh_index, size_t left_vid, size_t right_vid,
  std::vector<std::string>* propagation_logs)
{
  (void)contributors_raw_json;
  (void)left_vid;
  (void)right_vid;

  auto log_propagate_call = [&](size_t count, const std::string& path_scheme, const char* note)
  {
    if (!glueAlignLogEnabled(meshlet_index))
    {
      return;
    }
    std::ostringstream line;
    line << "(" << formatContribLabel(contribs, source_mesh_index) << ", "
         << formatContribLabel(contribs, partner_mesh_index) << ") propagated x" << count << ": ";
    if (!path_scheme.empty())
    {
      line << path_scheme;
    }
    else
    {
      line << "<no scheme>";
    }
    if (note != nullptr && note[0] != '\0')
    {
      line << " [" << note << "]";
    }
    KINDS_INFO("glue align: " << line.str());
    if (propagation_logs != nullptr)
    {
      propagation_logs->push_back(line.str());
    }
  };

  if (glueAlignLogEnabled(meshlet_index))
  {
    const std::string enter_scheme = formatGlueSeamPathScheme(source_path, partner_path);
    std::ostringstream enter;
    enter << "glue align: enter propagateUnmatchedInteriorsOntoPartner"
          << " source=" << formatContribLabel(contribs, source_mesh_index)
          << " partner=" << formatContribLabel(contribs, partner_mesh_index)
          << " source_path_len=" << source_path.size() << " partner_path_len=" << partner_path.size();
    if (!enter_scheme.empty())
    {
      enter << " path=" << enter_scheme;
    }
    KINDS_INFO(enter.str());
  }

  if (source_path.size() < 2 || partner_path.size() < 2)
  {
    log_propagate_call(0, {}, "paths too short");
    KINDS_ERROR("Source and partner paths are too short");
    return false;
  }

  const auto& source_verts = source_mesh.getVertices();
  const auto& partner_verts = partner_mesh.getVertices();
  auto path_endpoint_ok = [](const std::vector<size_t>& path, const std::vector<glm::dvec3>& verts)
  { return !path.empty() && path.front() < verts.size() && path.back() < verts.size(); };
  if (!path_endpoint_ok(source_path, source_verts) || !path_endpoint_ok(partner_path, partner_verts))
  {
    log_propagate_call(0, {}, "endpoint vertex index out of range");
    KINDS_ERROR("Source or partner path endpoint vertex index out of range");
    return false;
  }

  const glm::dvec3& source_front = source_verts[source_path.front()];
  const glm::dvec3& source_back = source_verts[source_path.back()];
  const glm::dvec3& partner_front = partner_verts[partner_path.front()];
  const glm::dvec3& partner_back = partner_verts[partner_path.back()];

  // Indices are local to each mesh; align by endpoint coordinates (exact match, as in glue matching).
  if (source_front == partner_back && source_back == partner_front)
  {
    partner_path = std::vector<size_t>(partner_path.rbegin(), partner_path.rend());
    partner_triangles = std::vector<size_t>(partner_triangles.rbegin(), partner_triangles.rend());
  }

  const glm::dvec3& partner_front_aligned = partner_verts[partner_path.front()];
  const glm::dvec3& partner_back_aligned = partner_verts[partner_path.back()];
  if (source_front != partner_front_aligned || source_back != partner_back_aligned)
  {
    auto format_path = [](const std::vector<size_t>& path)
    {
      std::ostringstream oss;
      oss << "[";
      for (size_t i = 0; i < path.size(); ++i)
      {
        if (i > 0)
        {
          oss << ", ";
        }
        oss << path[i];
      }
      oss << "]";
      return oss.str();
    };
    auto format_pos = [](const glm::dvec3& p)
    {
      std::ostringstream oss;
      oss << "(" << p.x << ", " << p.y << ", " << p.z << ")";
      return oss.str();
    };
    const std::string scheme = formatGlueSeamPathScheme(
      source_path, std::vector<size_t> { partner_path.front(), partner_path.back() });
    log_propagate_call(0, scheme, "paths not aligned");
    KINDS_ERROR("Source and partner paths are not aligned"
      << " source_path=" << format_path(source_path) << " partner_path=" << format_path(partner_path)
      << " source_front=" << format_pos(source_front) << " source_back=" << format_pos(source_back) << " partner_front="
      << format_pos(partner_front_aligned) << " partner_back=" << format_pos(partner_back_aligned));
    return false;
  }

  const std::string seam_scheme = formatGlueSeamPathScheme(source_path, partner_path);

  std::vector<bool> partner_path_propagated(partner_path.size(), false);

  auto find_closest_triangle_index
    = [](const VoronoiMesh& mesh, const std::vector<size_t>& target_triangles, const glm::dvec3& vertex)
  {
    // TODO: could be sped up with squared distance
    double min_dist = std::numeric_limits<double>::infinity();
    size_t closest_tri_index = -1;
    for (size_t i = 0; i < target_triangles.size(); ++i)
    {
      size_t tri_id = target_triangles[i];
      const auto& tris = mesh.getTriangles();
      if (tri_id * 3 + 2 >= tris.size())
      {
        continue;
      }

      const glm::dvec3& a = mesh.getVertices()[tris[tri_id * 3]];
      const glm::dvec3& b = mesh.getVertices()[tris[tri_id * 3 + 1]];
      const glm::dvec3& c = mesh.getVertices()[tris[tri_id * 3 + 2]];

      double distance = GeometryUtils::distancePointTriangle(vertex, a, b, c);
      if (distance < min_dist)
      {
        min_dist = distance;
        closest_tri_index = i;
      }
    }
    return closest_tri_index;
  };

  /// Insert interiors from @p from_path onto @p target_path via closest-edge triangle splits.
  /// Skips endpoints. When @p from_skip is set, skips interior @c i with @c (*from_skip)[i] == true
  /// (already present via the opposite-direction propagation). On success, inserts into
  /// @p target_path / @p target_triangles and, if provided, @c true into @p target_propagated.
  auto propagate_path_interiors
    = [&](VoronoiMesh& target_mesh, std::vector<size_t>& target_path, std::vector<size_t>& target_triangles,
        std::vector<bool>* target_propagated, const VoronoiMesh& from_mesh, const std::vector<size_t>& from_path,
        const std::vector<bool>* from_skip, const char* direction_label) -> size_t
  {
    size_t count = 0;
    for (size_t i = 1; i + 1 < from_path.size(); ++i)
    {
      if (from_skip != nullptr && i < from_skip->size() && (*from_skip)[i])
      {
        continue;
      }

      const glm::dvec3& vertex = from_mesh.getVertices()[from_path[i]];
      const double split_t = from_mesh.vertexKineticTime(from_path[i]);
      const std::string split_meta = target_mesh.storeMetadata()
        ? SegmentBuilder::MetadataBuilder()
            .addString("event_type", "glue_align")
            .addString("source", "propagated_interior")
            .addString("pos", direction_label)
            .addBool("split_triangle", true)
            .addSize("meshlet_index", meshlet_index)
            .addInt("physics_segment_id", physics_segment_id)
            .addDouble("t", std::isfinite(split_t) ? split_t : vertex.z)
            .addDouble("time", std::isfinite(split_t) ? split_t : vertex.z)
            .build()
        : std::string {};
      size_t closest_tri_index = find_closest_triangle_index(target_mesh, target_triangles, vertex);
      if (closest_tri_index == size_t(-1) || closest_tri_index + 1 >= target_path.size()
        || closest_tri_index >= target_triangles.size())
      {
        KINDS_ERROR("No closest target path triangle for interior"
          << " direction=" << direction_label << " from_path_i=" << i << " from_vid=" << from_path[i]
          << " target_path_len=" << target_path.size() << " target_triangles_len=" << target_triangles.size());
        continue;
      }

      size_t split_tri_id = target_triangles[closest_tri_index];
      size_t vertex_id0 = target_path[closest_tri_index];
      size_t vertex_id1 = target_path[closest_tri_index + 1];

      size_t tri_vertex_id0 = target_mesh.triangleCornerIndex(split_tri_id, vertex_id0);
      size_t tri_vertex_id1 = target_mesh.triangleCornerIndex(split_tri_id, vertex_id1);
      if (tri_vertex_id0 == size_t(-1) || tri_vertex_id1 == size_t(-1))
      {
        KINDS_ERROR("triangleCornerIndex failed for split"
          << " direction=" << direction_label << " from_path_i=" << i << " from_vid=" << from_path[i]
          << " split_tri_id=" << split_tri_id << " vertex_id0=" << vertex_id0 << " vertex_id1=" << vertex_id1
          << " tri_vertex_id0=" << tri_vertex_id0 << " tri_vertex_id1=" << tri_vertex_id1);
        continue;
      }

      const auto [new_vid, new_tri] = target_mesh.splitTriangle(tri_vertex_id0, tri_vertex_id1, vertex, split_meta,
        std::isfinite(split_t) ? std::optional<double>(split_t) : std::nullopt);
      if (new_vid == size_t(-1) || new_tri == size_t(-1))
      {
        KINDS_ERROR("splitTriangle failed"
          << " direction=" << direction_label << " from_path_i=" << i << " from_vid=" << from_path[i]
          << " split_tri_id=" << split_tri_id << " vertex_id0=" << vertex_id0 << " vertex_id1=" << vertex_id1);
        continue;
      }

      // splitTriangle orients to face winding: original keeps (va→F), new gets (F→vb).
      // Path order vertex_id0→vertex_id1 may be opposite, so assign edge tris from which face owns v0.
      const auto insert_at = static_cast<std::ptrdiff_t>(closest_tri_index + 1);
      target_path.insert(target_path.begin() + insert_at, new_vid);
      if (target_mesh.triangleCornerIndex(split_tri_id, vertex_id0) != size_t(-1))
      {
        // Original face has vertex_id0→F; new face has F→vertex_id1.
        target_triangles.insert(target_triangles.begin() + insert_at, new_tri);
      }
      else
      {
        // Path opposed winding: new face has vertex_id0→F; original has F→vertex_id1.
        target_triangles[closest_tri_index] = new_tri;
        target_triangles.insert(target_triangles.begin() + insert_at, split_tri_id);
      }
      if (target_propagated != nullptr)
      {
        target_propagated->insert(target_propagated->begin() + insert_at, true);
      }
      ++count;
    }
    return count;
  };

  size_t propagation_count = 0;

  // Source interiors → partner seam.
  propagation_count += propagate_path_interiors(partner_mesh, partner_path, partner_triangles, &partner_path_propagated,
    source_mesh, source_path, /*from_skip=*/nullptr, "source-->partner");

  // Partner interiors that were not inserted from the source → source seam.
  propagation_count += propagate_path_interiors(source_mesh, source_path, source_triangles, /*target_propagated=*/nullptr,
    partner_mesh, partner_path, &partner_path_propagated, "partner-->source");

  log_propagate_call(propagation_count, seam_scheme, propagation_count == 0 ? "no inserts" : nullptr);
  return propagation_count > 0;
}

/// For every consecutive pair of verts on @p source_mesh matched to @p partner_mesh (along each
/// boundary cycle), require a reciprocal consecutive matched pair on the partner toward the source.
/// If that reciprocal pair is missing, warn and skip; otherwise propagate unmatched interiors onto
/// the partner arc. Warns again if propagation leaves an interior missing on the partner.
bool processConsecutiveMatchedPairsTowardPartner(size_t source_mesh, size_t partner_mesh,
  const MeshBoundaryInfo& source_info, const MeshBoundaryInfo& partner_info,
  std::vector<SegmentAssemblyContrib>& contribs, size_t meshlet_index, int physics_segment_id,
  std::vector<std::string>* propagation_logs)
{
  if (source_mesh >= contribs.size() || partner_mesh >= contribs.size())
  {
    return false;
  }

  bool edited = false;
  const auto& sverts = contribs[source_mesh].mesh.getVertices();
  const std::string contributors_ij
    = formatContributorsJson({ contribs[source_mesh].contributorLabel(), contribs[partner_mesh].contributorLabel() });

  for (const auto& cycle : source_info.cycles)
  {
    if (cycle.size() < 2)
    {
      continue;
    }
    const size_t cn = cycle.size();
    std::vector<size_t> matched_cycle_idx;
    std::vector<size_t> matched_partner_vid;
    for (size_t ci = 0; ci < cn; ++ci)
    {
      if (const auto partner_vid = partnerVertexOnMesh(source_info, cycle[ci], partner_mesh, contribs);
        partner_vid.has_value())
      {
        matched_cycle_idx.push_back(ci);
        matched_partner_vid.push_back(partner_vid.value());
      }
    }
    const size_t m = matched_cycle_idx.size();
    if (m < 2)
    {
      continue;
    }

    for (size_t k = 0; k < m; ++k)
    {
      const size_t k1 = (k + 1) % m;
      const size_t c0 = matched_cycle_idx[k];
      const size_t c1 = matched_cycle_idx[k1];
      const size_t left_vid = cycle[c0];
      const size_t right_vid = cycle[c1];
      const size_t left_partner_vid = matched_partner_vid[k];
      const size_t right_partner_vid = matched_partner_vid[k1];
      if (left_vid >= sverts.size() || right_vid >= sverts.size())
      {
        continue;
      }
      const glm::dvec3 left_pos = sverts[left_vid];
      const glm::dvec3 right_pos = sverts[right_vid];
      const auto& pverts = contribs[partner_mesh].mesh.getVertices();
      const glm::dvec3 left_partner_pos = left_partner_vid < pverts.size() ? pverts[left_partner_vid] : glm::dvec3(0.0);
      const glm::dvec3 right_partner_pos
        = right_partner_vid < pverts.size() ? pverts[right_partner_vid] : glm::dvec3(0.0);

      std::vector<size_t> source_path;
      std::vector<size_t> source_triangles;
      if (!cycleArcPathWithTriangles(contribs[source_mesh].mesh, cycle, c0, c1, source_path, source_triangles))
      {
        continue;
      }

      // Skip the long way around the cycle; the complementary consecutive match pair covers the short seam.
      if (isStrictlyLongerCycleArc(source_path.size(), cn))
      {
        continue;
      }

      std::vector<size_t> partner_path;
      std::vector<size_t> partner_triangles;
      if (!findReciprocalConsecutivePartnerPath(contribs[partner_mesh].mesh, partner_info, source_mesh,
            left_partner_vid, right_partner_vid, source_path.size(), contribs, partner_path, partner_triangles))
      {
        // No reciprocal partner arc: still dump both boundary arcs (source first, partner second).
        // Partner arc is collected without the unmatched-interior requirement so interrupters appear.
        std::vector<size_t> partner_display_path;
        if (!findPartnerBoundaryArcPreferOpposite(partner_info, left_partner_vid, right_partner_vid, source_path.size(),
              partner_display_path, /*require_unmatched_interiors=*/false))
        {
          partner_display_path = { left_partner_vid, right_partner_vid };
        }
        const std::string path_scheme = formatGlueSeamPathScheme(source_path, partner_display_path, &source_info,
          &partner_info, /*mark_non_pair_matches=*/true);
        warnGlueMatchedPair(meshlet_index, physics_segment_id,
          "cannot process matched pair (no reciprocal consecutive matched pair on partner seam)", contribs, source_mesh,
          partner_mesh, left_vid, right_vid, left_pos, right_pos, left_partner_vid, right_partner_vid, left_partner_pos,
          right_partner_pos, path_scheme);
        continue;
      }

      // Nothing to propagate when both sides are a pure matched edge (no unmatched interiors either way).
      if (source_path.size() == 2 && partner_path.size() == 2)
      {
        continue;
      }

      if (propagateUnmatchedInteriorsOntoPartner(contribs[partner_mesh].mesh, partner_path, partner_triangles,
            contribs[source_mesh].mesh, source_path, source_triangles, meshlet_index, physics_segment_id,
            contributors_ij, contribs, source_mesh, partner_mesh, left_vid, right_vid, propagation_logs))
      {
        edited = true;
      }
    }
  }
  return edited;
}

void logGlueBoundaryWalks(const KineticDelaunay& kin_del, size_t meshlet_index, int physics_segment_id,
  const std::vector<SegmentAssemblyContrib>& contribs, const std::vector<MeshBoundaryInfo>& boundaries)
{
  if (!glueAlignLogEnabled(meshlet_index))
  {
    return;
  }

  std::ostringstream oss;
  oss << "glue align: boundary walks meshlet_index=" << meshlet_index << " physics_segment_id=" << physics_segment_id
      << "\n";

  std::filesystem::path filepath = glueAlignBoundaryVertsTxtPath(kin_del, meshlet_index, physics_segment_id);
  if (filepath.has_parent_path())
  {
    std::filesystem::create_directories(filepath.parent_path());
  }
  std::ofstream out(filepath);
  if (!out)
  {
    KINDS_WARNING("glue align: failed to open boundary vertex TXT " << filepath.generic_string());
  }
  else
  {
    out << std::setprecision(17);
    out << "# meshlet_index=" << meshlet_index << " physics_segment_id=" << physics_segment_id
        << " contrib_count=" << contribs.size() << "\n";
  }

  for (size_t i = 0; i < contribs.size() && i < boundaries.size(); ++i)
  {
    const std::string label = contribs[i].contributorLabel();
    const auto& cycles = boundaries[i].cycles;
    const auto& verts = contribs[i].mesh.getVertices();
    if (cycles.empty())
    {
      oss << "  " << label << ": <no boundary cycle>\n";
      if (out)
      {
        out << "# " << label << ": <no boundary cycle> mesh_vertex_count=" << verts.size() << "\n";
      }
      continue;
    }
    for (size_t c = 0; c < cycles.size(); ++c)
    {
      const auto& cycle = cycles[c];
      std::ostringstream cycle_line;
      cycle_line << label;
      if (cycles.size() > 1)
      {
        cycle_line << "[cycle" << c << "]";
      }
      cycle_line << ": ";
      for (size_t k = 0; k < cycle.size(); ++k)
      {
        if (k > 0)
        {
          cycle_line << ", ";
        }
        cycle_line << cycle[k];
      }
      if (!cycle.empty())
      {
        cycle_line << ", " << cycle.front();
      }
      oss << "  " << cycle_line.str() << "\n";

      if (out)
      {
        out << "# " << cycle_line.str() << " (cycle_len=" << cycle.size() << " mesh_vertex_count=" << verts.size()
            << ")\n";
        for (const size_t vid : cycle)
        {
          out << vid;
          if (vid < verts.size())
          {
            out << " " << verts[vid].x << " " << verts[vid].y << " " << verts[vid].z;
          }
          else
          {
            out << " <oob>";
          }
          out << "\n";
        }
      }
    }
  }

  if (out)
  {
    out << "\n# matching segments (reciprocal consecutive matched pairs)\n";
    out << "# format: (source_vid, partner_vid); then source interiors (src, _), then partner interiors (_, partner)\n";
    out << "# reciprocal = partner arc between endpoints has unmatched-only interiors (no match to any mesh);\n";
    out << "# partner cycle order prefers opposite of source (right→left), else same order\n";
    out << "# closed cycle (first endpoint repeated) only when every vertex on the source cycle is matched\n";
    out << "# wrap-around consecutive match pairs are always considered; long-way cycle arcs are skipped\n";
    out << "# each unordered contrib pair is logged once (source index < partner index)\n";
    out << "# unmatched consecutive pairs: consecutive on source matches, but partner arc is interrupted\n";
    out << "#   by at least one vertex matched to some contrib; those interrupting indices are marked with '?'\n";
  }

  const size_t n = contribs.size();
  for (size_t i = 0; i < n; ++i)
  {
    for (size_t j = i + 1; j < n; ++j)
    {
      if (i >= boundaries.size() || j >= boundaries.size())
      {
        continue;
      }
      const std::string label_i = contribs[i].contributorLabel();
      const std::string label_j = contribs[j].contributorLabel();
      const MeshBoundaryInfo& source_info = boundaries[i];
      const MeshBoundaryInfo& partner_info = boundaries[j];

      std::vector<std::string> match_chains;
      std::vector<std::string> unmatched_pairs;

      for (const auto& cycle : source_info.cycles)
      {
        if (cycle.size() < 2)
        {
          continue;
        }
        const size_t cn = cycle.size();
        std::vector<size_t> matched_cycle_idx;
        std::vector<size_t> matched_src_vid;
        std::vector<size_t> matched_partner_vid;
        for (size_t ci = 0; ci < cn; ++ci)
        {
          if (const auto partner_vid = partnerVertexOnMesh(source_info, cycle[ci], j, contribs);
            partner_vid.has_value())
          {
            matched_cycle_idx.push_back(ci);
            matched_src_vid.push_back(cycle[ci]);
            matched_partner_vid.push_back(partner_vid.value());
          }
        }
        const size_t m = matched_src_vid.size();
        if (m < 2)
        {
          continue;
        }

        // Always consider wrap-around consecutive match pairs. Only close a match chain by repeating
        // the first endpoint when every cycle vertex is matched (full-cycle seam).
        const bool close_full_cycle = (m == cn);

        std::ostringstream chain;
        size_t chain_len = 0;
        auto flush_chain = [&]()
        {
          if (chain_len >= 2)
          {
            match_chains.push_back(chain.str());
          }
          chain.str("");
          chain.clear();
          chain_len = 0;
        };
        auto append_token = [&](const std::string& token)
        {
          if (chain_len > 0)
          {
            chain << " - ";
          }
          chain << token;
          ++chain_len;
        };
        auto append_matched = [&](size_t src_vid, size_t partner_vid)
        {
          append_token("(" + std::to_string(src_vid) + ", " + std::to_string(partner_vid) + ")");
        };
        auto append_source_unmatched = [&](size_t src_vid)
        { append_token("(" + std::to_string(src_vid) + ", _)"); };
        auto append_partner_unmatched = [&](size_t partner_vid)
        { append_token("(_, " + std::to_string(partner_vid) + ")"); };

        auto append_arc_interiors = [&](const std::vector<size_t>& source_path, const std::vector<size_t>& partner_path)
        {
          for (size_t t = 1; t + 1 < source_path.size(); ++t)
          {
            append_source_unmatched(source_path[t]);
          }
          for (size_t t = 1; t + 1 < partner_path.size(); ++t)
          {
            append_partner_unmatched(partner_path[t]);
          }
        };

        for (size_t k = 0; k < m; ++k)
        {
          const size_t k1 = (k + 1) % m;
          const bool is_wrap_pair = (k + 1 == m);
          const size_t left_src = matched_src_vid[k];
          const size_t right_src = matched_src_vid[k1];
          const size_t left_partner = matched_partner_vid[k];
          const size_t right_partner = matched_partner_vid[k1];

          const std::vector<size_t> source_path
            = cycleArcPathInclusive(cycle, matched_cycle_idx[k], matched_cycle_idx[k1]);
          if (source_path.size() < 2 || isStrictlyLongerCycleArc(source_path.size(), cn))
          {
            continue;
          }

          std::vector<size_t> partner_path;
          std::vector<size_t> partner_triangles;
          const bool reciprocal = findReciprocalConsecutivePartnerPath(contribs[j].mesh, partner_info, i, left_partner,
            right_partner, source_path.size(), contribs, partner_path, partner_triangles);
          if (reciprocal)
          {
            if (is_wrap_pair && !close_full_cycle)
            {
              // Keep wrap pairs, but don't splice them onto an open chain (avoids false closed loops).
              flush_chain();
              append_matched(left_src, left_partner);
              append_arc_interiors(source_path, partner_path);
              append_matched(right_src, right_partner);
              flush_chain();
            }
            else
            {
              if (chain_len == 0)
              {
                append_matched(left_src, left_partner);
              }
              append_arc_interiors(source_path, partner_path);
              append_matched(right_src, right_partner);
            }
          }
          else
          {
            flush_chain();
            // Dump both arcs even when reciprocity fails (no unmatched-interior filter on partner).
            // Order unknown for pairing: list source interiors first, partner interiors second.
            // Indices matched to some contrib but not the consecutive endpoints get a trailing '?'.
            std::vector<size_t> partner_display_path;
            if (!findPartnerBoundaryArcPreferOpposite(partner_info, left_partner, right_partner, source_path.size(),
                  partner_display_path, /*require_unmatched_interiors=*/false))
            {
              partner_display_path = { left_partner, right_partner };
            }
            unmatched_pairs.push_back(formatGlueSeamPathScheme(
              source_path, partner_display_path, &source_info, &partner_info, /*mark_non_pair_matches=*/true));
          }
        }
        flush_chain();
      }

      if (out && (!match_chains.empty() || !unmatched_pairs.empty()))
      {
        out << "(" << label_i << ", " << label_j << ")";
        if (!match_chains.empty())
        {
          out << " matches:";
          for (size_t c = 0; c < match_chains.size(); ++c)
          {
            if (c > 0)
            {
              out << " |";
            }
            out << " " << match_chains[c];
          }
        }
        else
        {
          out << " matches: <none>";
        }
        out << "\n";
        if (!unmatched_pairs.empty())
        {
          out << "(" << label_i << ", " << label_j << ") unmatched consecutive pairs:";
          for (size_t u = 0; u < unmatched_pairs.size(); ++u)
          {
            if (u > 0)
            {
              out << ";";
            }
            out << " " << unmatched_pairs[u];
          }
          out << "\n";
        }
      }
    }
  }

  KINDS_INFO(oss.str());
  if (out)
  {
    out.close();
    KINDS_INFO("glue align: wrote boundary vertex coordinates to " << filepath.generic_string());
  }
}

void alignGlueEdgesOnContribCopies(const KineticDelaunay& kin_del, std::vector<SegmentAssemblyContrib>& contribs,
  size_t meshlet_index, int physics_segment_id)
{
  if (contribs.size() < 2)
  {
    return;
  }

  if (glueAlignLogEnabled(meshlet_index))
  {
    KINDS_INFO("glue align: entering computation for meshlet_index="
      << meshlet_index << " physics_segment_id=" << physics_segment_id << " contrib_count=" << contribs.size());
  }

  const size_t n = contribs.size();
  std::vector<MeshBoundaryInfo> boundaries(n);
  for (size_t i = 0; i < n; ++i)
  {
    boundaries[i] = buildMeshBoundaryInfo(contribs[i].mesh);
    initBoundaryMatchedMetadataFalse(contribs[i].mesh, boundaries[i], contribs[i]);
  }
  buildBoundaryMultiMatches(boundaries, contribs);
  logGlueBoundaryWalks(kin_del, meshlet_index, physics_segment_id, contribs, boundaries);

  // Propagate unmatched interiors for every consecutive matched pair with a reciprocal partner path.
  // Match tables / cycles stay frozen; only the local partner_path is updated after each split.
  // Each unordered contrib pair is processed once (source=i, partner=j with i < j) so roles never swap.
  std::vector<bool> contrib_topology_edited(n, false);
  std::vector<std::string> propagation_logs;
  for (size_t i = 0; i < n; ++i)
  {
    for (size_t j = i + 1; j < n; ++j)
    {
      if (processConsecutiveMatchedPairsTowardPartner(
            i, j, boundaries[i], boundaries[j], contribs, meshlet_index, physics_segment_id, &propagation_logs))
      {
        contrib_topology_edited[i] = true;
        contrib_topology_edited[j] = true;
      }
    }
  }

  if (glueAlignLogEnabled(meshlet_index) && !propagation_logs.empty())
  {
    const std::filesystem::path filepath = glueAlignBoundaryVertsTxtPath(kin_del, meshlet_index, physics_segment_id);
    std::ofstream out(filepath, std::ios::app);
    if (!out)
    {
      KINDS_WARNING("glue align: failed to append propagation paths to " << filepath.generic_string());
    }
    else
    {
      out << "\n# propagateUnmatchedInteriorsOntoPartner path pairs (includes propagated x0)\n";
      for (const std::string& line : propagation_logs)
      {
        out << line << "\n";
      }
      KINDS_INFO("glue align: appended " << propagation_logs.size() << " propagated path(s) to "
                                         << filepath.generic_string());
    }
  }

  for (size_t i = 0; i < n; ++i)
  {
    if (!contrib_topology_edited[i])
    {
      continue;
    }
    const NormalMode mode = contribs[i].mesh.getNormalMode();
    if (mode != NormalMode::NoNormals)
    {
      contribs[i].mesh.computeNormals(mode);
    }
  }
}
} // namespace

std::pair<std::vector<VoronoiMesh>, std::vector<std::vector<int>>> kinDS::SegmentBuilder::extractSegmentMeshlets(
  bool merge_by_segment) const
{
  if (!merge_by_segment)
  {
    std::vector<VoronoiMesh> raw_meshlets = meshes;
    std::vector<std::vector<int>> raw_neighbors;
    raw_neighbors.reserve(raw_meshlets.size());

    for (size_t meshlet_index = 0; meshlet_index < raw_meshlets.size(); ++meshlet_index)
    {
      raw_meshlets[meshlet_index].ensureFaceMetadataSize();
      std::vector<int> neighbors_for_meshlet;
      neighbors_for_meshlet.assign(raw_meshlets[meshlet_index].getTriangleCount(), -1);
      raw_neighbors.push_back(std::move(neighbors_for_meshlet));
    }

    return std::make_pair(raw_meshlets, raw_neighbors);
  }

  std::vector<VoronoiMesh> meshlets;
  std::vector<std::vector<int>> neighbor_segments; // accessed as [segment_id][triangle_index]

  // Same mapping as TreeMesher::mapMeshingToPhysicsSegmentIndices — meshlet_index == meshing segment id.
  std::vector<int> meshing_to_physics(segment_properties.size(), -1);
  {
    const auto& physics_map = kin_del.getStrandTree().getPhysicsStrandToSegmentIndices();
    const size_t strand_count = std::min(strand_to_segment_indices.size(), physics_map.size());
    for (size_t strand_id = 0; strand_id < strand_count; ++strand_id)
    {
      const size_t seg_count = std::min(strand_to_segment_indices[strand_id].size(), physics_map[strand_id].size());
      for (size_t segment_no = 0; segment_no < seg_count; ++segment_no)
      {
        const size_t meshing_segment_id = strand_to_segment_indices[strand_id][segment_no];
        if (meshing_segment_id < meshing_to_physics.size())
        {
          meshing_to_physics[meshing_segment_id] = physics_map[strand_id][segment_no];
        }
      }
    }
  }

  for (size_t segment_id = 0; segment_id < segment_properties.size(); ++segment_id)
  {
    // Match finalized meshlets (PerTriangleCorner). Empty segments must not stay NoNormals or
    // combined export can lock the wrong mode before the first geometric append.
    VoronoiMesh segment_mesh(MeshletExportMaterialNames, NormalMode::PerTriangleCorner);
    segment_mesh.setStoreMetadata(store_mesh_metadata);
    bool segment_mesh_initialized = false;
    std::vector<int> neighbor_segments_for_meshlet;
    const auto& properties = segment_properties[segment_id];
    double earliest_creation = std::numeric_limits<double>::quiet_NaN();

    std::vector<SegmentAssemblyContrib> contribs;
    contribs.reserve(properties.neighbor_count + 4);

    // Regular meshlets referenced through segment_properties.
    for (size_t neighbor_index = 0; neighbor_index < properties.neighbor_count; ++neighbor_index)
    {
      const size_t mesh_pair_index = properties.mesh_pair_indices[neighbor_index];
      SegmentAssemblyContrib contrib;
      contrib.mesh = meshes[mesh_pair_index];
      contrib.source_index = mesh_pair_index;
      contrib.is_boundary = false;
      const size_t seg0 = segment_mesh_pairs[mesh_pair_index].segment_index0;
      const size_t seg1 = segment_mesh_pairs[mesh_pair_index].segment_index1;
      if (seg0 != segment_id)
      {
        contrib.mesh.flipOrientation();
      }
      contrib.neighbor = (seg0 == segment_id) ? static_cast<int>(seg1) : static_cast<int>(seg0);
      const double mesh_ct = contrib.mesh.getCreationKineticTime();
      if (std::isfinite(mesh_ct) && (!std::isfinite(earliest_creation) || mesh_ct < earliest_creation))
      {
        earliest_creation = mesh_ct;
      }
      contribs.push_back(std::move(contrib));
    }

    // Boundary-interval meshlets: one owner segment, no orientation flip.
    for (size_t pair_idx = 0; pair_idx < intersection_segment_mesh_pairs.size(); ++pair_idx)
    {
      if (pair_idx >= intersection_meshes.size() || pair_idx >= intersection_mesh_pair_metadata.size())
      {
        continue;
      }
      const size_t owner_segment_id = intersection_mesh_pair_metadata[pair_idx].owner_segment_id;
      if (owner_segment_id == static_cast<size_t>(-1) || owner_segment_id != segment_id)
      {
        continue;
      }
      SegmentAssemblyContrib contrib;
      contrib.mesh = intersection_meshes[pair_idx];
      contrib.source_index = pair_idx;
      contrib.is_boundary = true;
      contrib.neighbor = -2;
      const double mesh_ct = contrib.mesh.getCreationKineticTime();
      if (std::isfinite(mesh_ct) && (!std::isfinite(earliest_creation) || mesh_ct < earliest_creation))
      {
        earliest_creation = mesh_ct;
      }
      contribs.push_back(std::move(contrib));
    }

    // Trace boundaries, match shared verts, align flexible interruptions on copies, then combine.
    if (align_flexible_glue_edges_postprocess_enabled)
    {
      const int physics_segment_id = (segment_id < meshing_to_physics.size()) ? meshing_to_physics[segment_id] : -1;
      // meshlet_index == meshing segment id == PerSegment / failed_meshlet OBJ index.
      alignGlueEdgesOnContribCopies(kin_del, contribs, segment_id, physics_segment_id);
    }

    for (SegmentAssemblyContrib& contrib : contribs)
    {
      // One OBJ group per contributor source (iN / bN) when this segment meshlet is exported.
      contrib.mesh.setGroupOffsets({ 0 });
      contrib.mesh.setGroupNames({ contrib.contributorLabel() });
      if (!segment_mesh_initialized)
      {
        segment_mesh = VoronoiMesh(MeshletExportMaterialNames, contrib.mesh.getNormalMode());
        segment_mesh.setStoreMetadata(store_mesh_metadata);
        segment_mesh_initialized = true;
      }
      neighbor_segments_for_meshlet.insert(
        neighbor_segments_for_meshlet.end(), contrib.mesh.getTriangleCount(), contrib.neighbor);
      segment_mesh += contrib.mesh;
    }

    if (std::isfinite(earliest_creation))
    {
      segment_mesh.setCreationKineticTime(earliest_creation);
    }
    neighbor_segments.push_back(neighbor_segments_for_meshlet);
    segment_mesh.mergeDuplicateVertices(0.0);
    postProcessFlexibleVerticesOnSegmentMeshlet(segment_mesh);
    segment_mesh.validateUVLayout("extractSegmentMeshlets merged segment mesh");
    segment_mesh.ensureFaceMetadataSize();
    meshlets.push_back(segment_mesh);

    if (diagnostics && segment_id == 0)
    {
      std::ostringstream oss;
      oss << "segment_id=0 verts=" << segment_mesh.getVertexCount() << " tris=" << segment_mesh.getTriangleCount()
          << " neighbor_count=" << properties.neighbor_count
          << " segment_mesh_initialized=" << (segment_mesh_initialized ? "true" : "false")
          << " contribs=" << contribs.size();
      strandInitDiagnosticLogLine(
        "extract_segment_meshlet", 0, std::numeric_limits<double>::quiet_NaN(), oss.str().c_str());
      if (!segment_mesh_initialized || segment_mesh.getVertexCount() == 0)
      {
        strandInitDiagnosticLogLine("extract_segment_meshlet_empty", 0, std::numeric_limits<double>::quiet_NaN(),
          "merged segment 0 meshlet is empty — check cap/strip wiring in accumulateSegmentProperties");
      }
    }
  }

  return std::make_pair(meshlets, neighbor_segments);
}

void SegmentBuilder::postProcessFlexibleVerticesOnSegmentMeshlet(VoronoiMesh& /*segment_mesh*/) const
{
  // Glue-edge flexible alignment runs during assembly in @ref extractSegmentMeshlets (on contrib copies).
  // Reserved for any remaining segment-meshlet flexible post-steps.
}

std::vector<std::string> kinDS::SegmentBuilder::extractSegmentMeshletExportSuffixes(bool merge_by_segment) const
{
  if (!merge_by_segment)
  {
    return meshlet_export_suffixes;
  }

  std::vector<std::string> suffixes;
  suffixes.reserve(segment_properties.size());

  std::vector<int> meshing_to_physics(segment_properties.size(), -1);
  std::vector<int> meshing_to_strand(segment_properties.size(), -1);
  for (size_t strand_id = 0; strand_id < strand_to_segment_indices.size(); ++strand_id)
  {
    for (size_t meshing_segment_id : strand_to_segment_indices[strand_id])
    {
      if (meshing_segment_id < meshing_to_strand.size())
      {
        meshing_to_strand[meshing_segment_id] = static_cast<int>(strand_id);
      }
    }
  }
  {
    const auto& physics_map = kin_del.getStrandTree().getPhysicsStrandToSegmentIndices();
    const size_t strand_count = std::min(strand_to_segment_indices.size(), physics_map.size());
    for (size_t strand_id = 0; strand_id < strand_count; ++strand_id)
    {
      const size_t seg_count = std::min(strand_to_segment_indices[strand_id].size(), physics_map[strand_id].size());
      for (size_t segment_no = 0; segment_no < seg_count; ++segment_no)
      {
        const size_t meshing_segment_id = strand_to_segment_indices[strand_id][segment_no];
        if (meshing_segment_id < meshing_to_physics.size())
        {
          meshing_to_physics[meshing_segment_id] = physics_map[strand_id][segment_no];
        }
      }
    }
  }

  for (size_t segment_id = 0; segment_id < segment_properties.size(); ++segment_id)
  {
    // meshlet{N}_segment{N}_strand{S}_phys{P}: N is meshing/meshlet index; S strand; P physics segment.
    std::string suffix = std::string("_segment") + std::to_string(segment_id);
    if (segment_id < meshing_to_strand.size() && meshing_to_strand[segment_id] >= 0)
    {
      suffix += "_strand" + std::to_string(meshing_to_strand[segment_id]);
    }
    if (segment_id < meshing_to_physics.size() && meshing_to_physics[segment_id] >= 0)
    {
      suffix += "_phys" + std::to_string(meshing_to_physics[segment_id]);
    }
    suffixes.push_back(std::move(suffix));
  }
  return suffixes;
}

std::vector<VoronoiMesh> kinDS::SegmentBuilder::extractBoundaryIntervalMeshlets() const { return intersection_meshes; }

std::vector<std::string> kinDS::SegmentBuilder::extractBoundaryIntervalMeshletExportSuffixes() const
{
  std::vector<std::string> out = intersection_meshlet_export_suffixes;
  if (out.size() < intersection_mesh_pair_metadata.size())
  {
    out.resize(intersection_mesh_pair_metadata.size());
  }
  const size_t count = std::min(out.size(), intersection_mesh_pair_metadata.size());
  for (size_t i = 0; i < count; ++i)
  {
    const auto& md = intersection_mesh_pair_metadata[i];
    std::string suffix = out[i];
    suffix += std::string("_cell") + std::to_string(md.voronoi_cell_id);
    if (md.start_delaunay_edge_id != static_cast<size_t>(-1))
    {
      suffix += std::string("_de0") + std::to_string(md.start_delaunay_edge_id);
    }
    if (md.end_delaunay_edge_id != static_cast<size_t>(-1))
    {
      suffix += std::string("_de1") + std::to_string(md.end_delaunay_edge_id);
    }
    out[i] = std::move(suffix);
  }
  return out;
}

const VoronoiMesh& kinDS::SegmentBuilder::getBoundaryMesh() const { return boundary_mesh; }

const std::vector<size_t>& kinDS::SegmentBuilder::getBoundaryVertexToStrandId() const
{
  return boundary_vertex_to_strand_id;
}

const std::vector<std::vector<size_t>>& kinDS::SegmentBuilder::getStrandToSegmentIndices() const
{
  return strand_to_segment_indices;
}

size_t kinDS::SegmentBuilder::strandIdForSegment(size_t segment_id) const
{
  for (size_t strand_id = 0; strand_id < strand_to_segment_indices.size(); ++strand_id)
  {
    for (size_t seg : strand_to_segment_indices[strand_id])
    {
      if (seg == segment_id)
      {
        return strand_id;
      }
    }
  }
  throw std::runtime_error("strandIdForSegment: unknown segment_id " + std::to_string(segment_id));
}

size_t kinDS::SegmentBuilder::strandIdForRawMeshlet(size_t meshlet_index) const
{
  if (meshlet_index >= segment_mesh_pairs.size())
  {
    throw std::runtime_error("strandIdForRawMeshlet: meshlet_index out of range.");
  }

  const auto& pair = segment_mesh_pairs[meshlet_index];
  size_t segment_id = pair.segment_index0;
  if (segment_id == static_cast<size_t>(-1))
  {
    segment_id = pair.segment_index1;
  }
  if (segment_id == static_cast<size_t>(-1))
  {
    throw std::runtime_error(
      "strandIdForRawMeshlet: meshlet " + std::to_string(meshlet_index) + " has no segment endpoint.");
  }
  return strandIdForSegment(segment_id);
}

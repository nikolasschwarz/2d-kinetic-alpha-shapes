#include "SegmentBuilder.hpp"

#include "KineticDelaunayEventPredicates.hpp"
#include "KineticDelaunayFlipEvent.hpp"
#include "KineticDelaunayRadiusEvent.hpp"
#include "Logger.hpp"
#include "SegmentBuilderCrossingCallback.hpp"
#include "SegmentBuilderFlipCallback.hpp"
#include "SegmentBuilderRadiusCallback.hpp"
#include "SegmentBuilderSectionCallback.hpp"
#include "SegmentBuilderSubdivisionCallback.hpp"
#include "SegmentBuilderSeparationCallback.hpp"
#include "Validator.hpp"

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
  o << std::setprecision(17) << std::showpoint << value;
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

bool triangulationPlaneXYEqual(const VoronoiMesh& mesh, size_t vertex_a, size_t vertex_b)
{
  const glm::dvec2 a = mesh.triangulationPlaneXY(vertex_a);
  const glm::dvec2 b = mesh.triangulationPlaneXY(vertex_b);
  return a.x == b.x && a.y == b.y;
}

std::filesystem::path makeTriangulateSimplePolygonDebugPath(const KineticDelaunay& kin_del, const VoronoiMesh& mesh,
  const char* tag, const char* extension)
{
  static size_t debug_counter = 0;
  ++debug_counter;

  const double kinetic_time = mesh.getCreationKineticTime();
  const std::string time_token
    = std::isfinite(kinetic_time) ? ("t" + numberLiteral(kinetic_time)) : std::string("t_unknown");
  const std::string filename = time_token + "_triangulateSimplePolygon_" + tag + "_"
    + std::to_string(debug_counter) + extension;

  std::filesystem::path filepath = std::filesystem::path("branch0") / filename;
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

void writeTriangulateSimplePolygonDebugTxt(const KineticDelaunay& kin_del, const VoronoiMesh& mesh,
  const std::filesystem::path& filepath, const char* tag,
  const std::vector<std::pair<std::string, std::vector<size_t>>>& rings)
{
  std::ofstream out(filepath);
  if (!out)
  {
    KINDS_WARNING("triangulateSimplePolygon: failed to open debug TXT " << filepath.generic_string());
    return;
  }

  const double kinetic_time = mesh.getCreationKineticTime();
  const auto& stored_vertices = mesh.getVertices();
  const auto& vertex_metadata = mesh.getVertexMetadata();

  out << "# tag=" << tag << '\n';
  out << "# kinetic_time=" << (std::isfinite(kinetic_time) ? numberLiteral(kinetic_time) : "unknown") << '\n';
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

    out << "index vertex_id profile_x profile_y object_x object_y object_z event_type\n";
    for (size_t i = 0; i < polygon_vertices.size(); ++i)
    {
      const size_t vertex_id = polygon_vertices[i];
      const glm::dvec2 profile = mesh.triangulationPlaneXY(vertex_id);
      const glm::dvec3 object = vertex_id < stored_vertices.size() ? stored_vertices[vertex_id] : glm::dvec3(0.0);
      out << i << ' ' << vertex_id << ' ' << numberLiteral(profile.x) << ' ' << numberLiteral(profile.y) << ' '
          << numberLiteral(object.x) << ' ' << numberLiteral(object.y) << ' ' << numberLiteral(object.z) << ' ';
      if (vertex_id < vertex_metadata.size())
      {
        const std::optional<std::string> event_type = metadataStringField(vertex_metadata[vertex_id], "event_type");
        if (event_type.has_value())
        {
          out << event_type.value();
        }
      }
      out << '\n';
    }
    out << '\n';
  }

  KINDS_WARNING("triangulateSimplePolygon: wrote debug TXT to " << filepath.generic_string());
}

std::vector<std::vector<size_t>> splitPolygonAtRepeatedVertices(const VoronoiMesh& mesh,
  const std::vector<size_t>& polygon)
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
        std::vector<size_t> sub1(ring.begin() + static_cast<std::ptrdiff_t>(start),
          ring.begin() + static_cast<std::ptrdiff_t>(end));

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

void writeTriangulateSimplePolygonFailSvg(const KineticDelaunay& kin_del, const VoronoiMesh& mesh,
  const std::vector<size_t>& polygon_vertices)
{
  if (polygon_vertices.size() < 3)
  {
    return;
  }

  const std::filesystem::path filepath
    = makeTriangulateSimplePolygonDebugPath(kin_del, mesh, "FAIL", ".svg");

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
  out << "<svg xmlns=\"http://www.w3.org/2000/svg\" viewBox=\"0 0 " << width << ' ' << height
      << "\" width=\"" << width * 100.0 << "\" height=\"" << height * 100.0 << "\">\n";
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
      const std::optional<std::string> event_type = metadataStringField(vertex_metadata[vertex_id], "event_type");
      if (event_type.has_value())
      {
        label << " " << event_type.value();
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

namespace
{
std::optional<glm::dvec3> meshVertexUv(const VoronoiMesh& mesh, size_t vertex_index)
{
  const std::vector<glm::dvec3>& uvs = mesh.getUVs();
  if (vertex_index < uvs.size())
  {
    return uvs[vertex_index];
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
  std::vector<glm::dvec3>& uvs = mesh.getUVs();
  if (vertex_index < uvs.size())
  {
    uvs[vertex_index] = uv;
  }
  for (size_t triangle_corner : mesh.findTriangleCorners(vertex_index))
  {
    if (mesh.hasValidUVIndex(triangle_corner))
    {
      mesh.setUV(uv, triangle_corner);
    }
  }
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

bool interpolateFlexibleVerticesAlongEdge(
  VoronoiMesh& mesh, std::vector<int>& flex, size_t anchor_old_vertex, size_t anchor_new_vertex)
{
  if (flex.empty())
  {
    return true;
  }
  auto& verts = mesh.getVertices();
  if (anchor_old_vertex >= verts.size() || anchor_new_vertex >= verts.size())
  {
    KINDS_WARNING("interpolateFlexibleVerticesAlongEdge: anchor index out of range (old=" << anchor_old_vertex
                                                                                         << " new="
                                                                                         << anchor_new_vertex
                                                                                         << " verts=" << verts.size()
                                                                                         << ").");
    return false;
  }
  const glm::dvec3 p0 = verts[anchor_old_vertex];
  const glm::dvec3 p1 = verts[anchor_new_vertex];
  if (!vertexPositionFinite(p0) || !vertexPositionFinite(p1))
  {
    KINDS_WARNING("interpolateFlexibleVerticesAlongEdge: non-finite anchor position (old=" << anchor_old_vertex
                                                                                           << " new="
                                                                                           << anchor_new_vertex
                                                                                           << ").");
    return false;
  }
  const double z0 = mesh.vertexKineticTime(anchor_old_vertex);
  const double z1 = mesh.vertexKineticTime(anchor_new_vertex);
  const double denom = z1 - z0;
  const std::optional<glm::dvec3> uv0 = meshVertexUv(mesh, anchor_old_vertex);
  const std::optional<glm::dvec3> uv1 = meshVertexUv(mesh, anchor_new_vertex);
  const bool interpolate_uv = uv0.has_value() && uv1.has_value();
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
      // Anchors share kinetic time: cannot parameterize by t; fall back to uniform spacing along the segment.
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
    if (interpolate_uv)
    {
      const glm::dvec3 uv_interp = *uv0 + s * (*uv1 - *uv0);
      setMeshVertexUv(mesh, fju, uv_interp);
    }
  }
  return true;
}

void snapFlexibleVerticesToAnchor(VoronoiMesh& mesh, const std::vector<int>& flex, size_t anchor_vertex)
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
  }
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
  const char* source,
  const std::optional<KineticDelaunay::CrossingData::EdgeIntersectionRef>& position_crossing,
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
    builder.addSize("delaunay_edge_id", geom_ref->delaunay_edge_id).addSize("voronoi_edge_id", geom_ref->voronoi_edge_id)
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

std::string makeBoundaryMeshFaceMetadataJson(double kinetic_time, const char* mesh_type, size_t half_edge_id,
  size_t delaunay_face_id, size_t input_branch_id)
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

bool diagnosticsNearMonitoredFlipTime(double t)
{
  return std::isfinite(t)
    && std::abs(t - SegmentBuilder::kDiagnosticsMonitoredFlipTime) <= SegmentBuilder::kDiagnosticsMonitoredTimeEpsilon;
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
    && delaunay_edge_id.value() == kDiagnosticsMonitoredDelaunayEdgeId;
  const bool pair_hit
    = mesh_pair_index.has_value() && mesh_pair_index.value() == kDiagnosticsMonitoredMeshPairId;
  if (!edge_hit && !pair_hit && !diagnosticsNearMonitoredFlipTime(t))
  {
    return;
  }
  logDiagnosticsMonitoredDelaunayEdgeState(t, event_context);
}

void kinDS::SegmentBuilder::logDiagnosticsMonitoredDelaunayEdgeState(double t, const char* event_context) const
{
  if (!diagnostics)
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
    oss << "{cell=" << md.voronoi_cell_id << " owner_seg=" << md.owner_segment_id << " start_d="
        << md.start_delaunay_edge_id << " end_d=" << md.end_delaunay_edge_id << "}";
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
      if (seg.start_crossing.has_value())
      {
        oss << " start_x{de=" << seg.start_crossing.value()->delaunay_edge_id
            << ",ve=" << seg.start_crossing.value()->voronoi_edge_id << "}";
      }
      if (seg.end_crossing.has_value())
      {
        oss << " end_x{de=" << seg.end_crossing.value()->delaunay_edge_id
            << ",ve=" << seg.end_crossing.value()->voronoi_edge_id << "}";
      }
      ++seg_i;
    }
    return oss.str();
  };

  const size_t d_edge = kDiagnosticsMonitoredDelaunayEdgeId;
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
  KINDS_INFO(header.str());

  if (he_odd >= graph.halfEdgeSlotCount())
  {
    KINDS_INFO("  edge slot out of bounds (he_slots=" << graph.halfEdgeSlotCount() << ")");
    return;
  }

  const bool on_boundary = kin_del.isOnComponentBoundary(he_even);
  const bool even_outside = on_boundary && kin_del.isOnComponentBoundaryOutside(he_even);
  KINDS_INFO("  alpha_boundary=" << (on_boundary ? "yes" : "no")
                                 << " even_is_outside=" << (even_outside ? "yes" : "no") << " even_origin="
                                 << graph.halfEdge(he_even).origin << " odd_origin=" << graph.halfEdge(he_odd).origin
                                 << " even_face=" << graph.halfEdge(he_even).face << " odd_face="
                                 << graph.halfEdge(he_odd).face);

  const auto& crossing_data = kin_del.getCrossingData();
  if (d_edge >= crossing_data.delaunay_edge_intersections.size())
  {
    KINDS_INFO("  no delaunay_edge_intersections slot");
    return;
  }

  const auto& refs = crossing_data.delaunay_edge_intersections[d_edge];
  KINDS_INFO("  crossing_count=" << refs.size());
  size_t list_idx = 0;
  for (const auto& ref : refs)
  {
    std::ostringstream line;
    line << "  [" << list_idx << "] ve=" << ref->voronoi_edge_id << " param=" << ref->delaunay_edge_param << " prev_pair="
         << formatMeshPairIndex(ref->prev_segment_mesh_pair_index) << " next_pair="
         << formatMeshPairIndex(ref->next_segment_mesh_pair_index);
    KINDS_INFO(line.str());
    if (ref->prev_segment_mesh_pair_index != static_cast<size_t>(-1))
    {
      KINDS_INFO("    prev_pair_meta " << format_intersection_mesh_pair_metadata(ref->prev_segment_mesh_pair_index)
                                       << " strip="
                                       << format_intersection_mesh_pair_strip_state(ref->prev_segment_mesh_pair_index));
    }
    if (ref->next_segment_mesh_pair_index != static_cast<size_t>(-1))
    {
      KINDS_INFO("    next_pair_meta " << format_intersection_mesh_pair_metadata(ref->next_segment_mesh_pair_index)
                                       << " strip="
                                       << format_intersection_mesh_pair_strip_state(ref->next_segment_mesh_pair_index));
    }
    if (ref->prev_segment_mesh_pair_index == kDiagnosticsMonitoredMeshPairId
      || ref->next_segment_mesh_pair_index == kDiagnosticsMonitoredMeshPairId)
    {
      KINDS_INFO("    ** references monitored mesh pair " << kDiagnosticsMonitoredMeshPairId << " **");
    }
    ++list_idx;
  }

  if (kDiagnosticsMonitoredMeshPairId < intersection_mesh_pair_metadata.size())
  {
    KINDS_INFO("  monitored_pair_" << kDiagnosticsMonitoredMeshPairId << "_meta="
                                     << format_intersection_mesh_pair_metadata(kDiagnosticsMonitoredMeshPairId));
    KINDS_INFO("  monitored_pair_" << kDiagnosticsMonitoredMeshPairId << "_strip="
                                     << format_intersection_mesh_pair_strip_state(kDiagnosticsMonitoredMeshPairId));
    if (kDiagnosticsMonitoredMeshPairId < intersection_meshes.size())
    {
      KINDS_INFO("  monitored_pair_" << kDiagnosticsMonitoredMeshPairId << "_verts="
                                       << intersection_meshes[kDiagnosticsMonitoredMeshPairId].getVertexCount()
                                       << " tris="
                                       << intersection_meshes[kDiagnosticsMonitoredMeshPairId].getTriangleCount()
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
  KINDS_INFO(oss.str());
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
    strandInitDiagnosticLogLine("init_summary_no_incident_strips", focus_strand, t,
      "no Voronoi strip meshlets indexed for incident half-edges");
  }

  const size_t section = static_cast<size_t>(t);
  const auto& branches_at_section = kin_del.getStrandTree().getStrandBranchesByHeight(section);
  for (size_t input_branch_id = 0; input_branch_id < branches_at_section.size(); ++input_branch_id)
  {
    const auto& branch_strands = branches_at_section[input_branch_id];
    const bool in_branch = std::find(branch_strands.begin(), branch_strands.end(), focus_strand) != branch_strands.end();
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
  if (nv > 0 && strips.empty())
  {
    KINDS_WARNING("meshlet_diag inconsistent after startNewMesh: dual_edge="
      << dual_edge << " t=" << t << " — mesh has " << nv << " vertices but strip list is empty.");
  }
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

glm::dvec3 SegmentBuilder::regularMeshStripIntervalEndpointPositionAt(const RegularMeshStripIntervalEndpoints& interval,
  bool at_start, size_t even_half_edge_id, size_t odd_half_edge_id, size_t voronoi_edge_id, double t,
  const RadiusBoundaryTransitionShiftContext* boundary_transition_shift) const
{
  const char* endpoint_label = at_start ? "start" : "end";
  if (at_start)
  {
    if (interval.start_open_voronoi_half_edge_id.has_value())
    {
      const glm::dvec3 p = computeVoronoiVertex(interval.start_open_voronoi_half_edge_id.value(), t);
      if (!std::isfinite(p.x) || !std::isfinite(p.y) || !std::isfinite(p.z) || (p.x == 0.0 && p.y == 0.0))
      {
        KINDS_WARNING("regularMeshStripIntervalEndpointPositionAt("
          << endpoint_label << "): degenerate circumcenter endpoint for voronoi_edge=" << voronoi_edge_id << " t=" << t
          << " -> (" << p.x << ", " << p.y << ", " << p.z << ").");
      }
      return p;
    }
    if (interval.start_crossing.has_value())
    {
      const KineticDelaunay::CrossingData::EdgeIntersectionRef orig_ref = interval.start_crossing.value();
      const RadiusBoundaryTransitionCrossingPlacement placement
        = resolveRadiusBoundaryTransitionCrossingPlacement(t, orig_ref, boundary_transition_shift);
      const glm::dvec3 old_chord_pos = crossingProfilePosition(t, placement.conceptual_intersection);
      const glm::dvec3 p = crossingProfilePositionFromPlacement(t, placement);
      if (!std::isfinite(p.x) || !std::isfinite(p.y) || !std::isfinite(p.z) || (p.x == 0.0 && p.y == 0.0))
      {
        KINDS_WARNING("regularMeshStripIntervalEndpointPositionAt("
          << endpoint_label << "): degenerate CrossingData endpoint for voronoi_edge=" << voronoi_edge_id << " t=" << t
          << " -> (" << p.x << ", " << p.y << ", " << p.z << ").");
      }
      if (placement.positionDiffersFromConceptual())
      {
        logRadiusBoundaryTransitionVertexShift("finishRegularMeshStripInterval_endpoint", t,
          placement.conceptual_intersection, placement.position_intersection, old_chord_pos, p);
      }
      return p;
    }
  }
  else
  {
    if (interval.end_open_voronoi_half_edge_id.has_value())
    {
      const glm::dvec3 p = computeVoronoiVertex(interval.end_open_voronoi_half_edge_id.value(), t);
      if (!std::isfinite(p.x) || !std::isfinite(p.y) || !std::isfinite(p.z) || (p.x == 0.0 && p.y == 0.0))
      {
        KINDS_WARNING("regularMeshStripIntervalEndpointPositionAt("
          << endpoint_label << "): degenerate circumcenter endpoint for voronoi_edge=" << voronoi_edge_id << " t=" << t
          << " -> (" << p.x << ", " << p.y << ", " << p.z << ").");
      }
      return p;
    }
    if (interval.end_crossing.has_value())
    {
      const KineticDelaunay::CrossingData::EdgeIntersectionRef orig_ref = interval.end_crossing.value();
      const RadiusBoundaryTransitionCrossingPlacement placement
        = resolveRadiusBoundaryTransitionCrossingPlacement(t, orig_ref, boundary_transition_shift);
      const glm::dvec3 old_chord_pos = crossingProfilePosition(t, placement.conceptual_intersection);
      const glm::dvec3 p = crossingProfilePositionFromPlacement(t, placement);
      if (placement.positionDiffersFromConceptual())
      {
        logRadiusBoundaryTransitionVertexShift("finishRegularMeshStripInterval_endpoint", t,
          placement.conceptual_intersection, placement.position_intersection, old_chord_pos, p);
      }
      return p;
    }
  }
  return glm::dvec3(0.0);
}

std::tuple<size_t, size_t> SegmentBuilder::finishRegularMeshStripInterval(VoronoiMesh& mesh,
  const std::vector<BoundaryPoint>& boundary_polygon, const glm::dvec2& centroid, size_t even_half_edge_id,
  size_t voronoi_edge_id, double t, size_t strand_vertex_id, int strand_even_origin_i, int strand_odd_origin_i,
  BoundaryEventType event_type, BoundarySegmentAction segment_action, const RegularMeshStripIntervalEndpoints& interval,
  size_t last_start_vertex_index, size_t last_end_vertex_index, const std::string& finish_face_metadata,
  const RadiusBoundaryTransitionShiftContext* boundary_transition_shift)
{
  const size_t odd_half_edge_id = even_half_edge_id + 1;
  const glm::dvec3 new_start_pos = regularMeshStripIntervalEndpointPositionAt(
    interval, true, even_half_edge_id, odd_half_edge_id, voronoi_edge_id, t, boundary_transition_shift);
  const glm::dvec3 new_end_pos = regularMeshStripIntervalEndpointPositionAt(
    interval, false, even_half_edge_id, odd_half_edge_id, voronoi_edge_id, t, boundary_transition_shift);

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
      const RadiusBoundaryTransitionCrossingPlacement placement
        = resolveRadiusBoundaryTransitionCrossingPlacement(t, crossing.value(), boundary_transition_shift);
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

  const size_t new_start_vertex_index
    = addMeshletVertex(mesh, boundary_polygon, centroid, new_start_pos, strand_vertex_id, t, false, start_vv, meta_start,
      std::nullopt,
      interval.start_crossing.has_value()
        ? MeshletVertexRuntimeInfo { false, false, interval.start_crossing, interval.start_crossing }
        : MeshletVertexRuntimeInfo {});
  const size_t new_end_vertex_index
    = addMeshletVertex(mesh, boundary_polygon, centroid, new_end_pos, strand_vertex_id, t, false, end_vv, meta_end,
      std::nullopt,
      interval.end_crossing.has_value()
        ? MeshletVertexRuntimeInfo { false, false, interval.end_crossing, interval.end_crossing }
        : MeshletVertexRuntimeInfo {});

  if (last_start_vertex_index == last_end_vertex_index)
  {
    addMeshletTriangle(mesh, new_start_vertex_index, last_end_vertex_index, new_end_vertex_index, finish_face_metadata);
  }
  else if (mesh.getVertices()[last_start_vertex_index][2] < mesh.getVertices()[last_end_vertex_index][2])
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
      KINDS_WARNING(
        "meshlet_diag finish_mesh: invalid pair index he=" << he_id << " t=" << t << " pair=" << segment_mesh_pair_index);
    }
    return operated_segments;
  }

  const bool retire_meshlet_at_subdivision
    = event_type == BoundaryEventType::Subdivision && segment_action == BoundarySegmentAction::SegmentCompleted;
  const auto mark_retired_if_subdivision = [&]()
  {
    if (retire_meshlet_at_subdivision)
    {
      markRegularMeshletCompleted(segment_mesh_pair_index);
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
      KINDS_WARNING("orientCrossingsAlongVoronoiEdge: skipping disconnected crossing at step " << step_index
        << " (voronoi_edge_id=" << voronoi_edge_id << ", delaunay_edge_id=" << d
        << ", current_face_id=" << current_face_id << ", face(he0)=" << graph.halfEdge(he0).face
        << ", face(he1)=" << graph.halfEdge(he1).face << ", voronoi_vertices=[" << left_voronoi_vertex_id << ","
        << right_voronoi_vertex_id << "])");
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
      ? intersectionCrossingVertexMetadata(
          MetadataBuilder()
            .addString("event_type", boundaryEventTypeToString(event_type))
            .addString("segment_action", boundarySegmentActionToString(segment_action))
            .build(),
          placement.conceptual_intersection, placement.position_intersection, "left",
          placement.explicit_profile_position.has_value())
      : std::string {};
    start_vertex_index = addMeshletVertex(mesh, boundary_polygon, centroid, pos,
      strand_id_for_inside_half_edge(start_half_edge_id), t, false, std::nullopt, meta_start, std::nullopt,
      MeshletVertexRuntimeInfo { false, placement.explicit_profile_position.has_value(),
        placement.position_intersection, placement.conceptual_intersection });
  }
  else if (interval.start_open_voronoi_half_edge_id.has_value())
  {
    const size_t open_he = interval.start_open_voronoi_half_edge_id.value();
    const glm::dvec3 pos = computeVoronoiVertex(open_he, t);
    const size_t voronoi_vertex_id = graph.halfEdge(open_he).face;
    const int origin = graph.halfEdge(open_he).origin;
    const std::string meta_start = composeRegularStripVertexMetadata(t, voronoi_edge_id, even_half_edge_id,
      strand_even_origin_i, strand_odd_origin_i, event_type, segment_action, std::nullopt, "left", nullptr,
      "Voronoi vertex");
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
      ? intersectionCrossingVertexMetadata(
          MetadataBuilder()
            .addString("event_type", boundaryEventTypeToString(event_type))
            .addString("segment_action", boundarySegmentActionToString(segment_action))
            .build(),
          placement.conceptual_intersection, placement.position_intersection, "right",
          placement.explicit_profile_position.has_value())
      : std::string {};
    end_vertex_index = addMeshletVertex(mesh, boundary_polygon, centroid, pos,
      strand_id_for_inside_half_edge(end_half_edge_id), t, false, std::nullopt, meta_end, std::nullopt,
      MeshletVertexRuntimeInfo { false, placement.explicit_profile_position.has_value(),
        placement.position_intersection, placement.conceptual_intersection });
  }
  else if (interval.end_open_voronoi_half_edge_id.has_value())
  {
    const size_t open_he = interval.end_open_voronoi_half_edge_id.value();
    const glm::dvec3 pos = computeVoronoiVertex(open_he, t);
    const size_t voronoi_vertex_id = graph.halfEdge(open_he).face;
    const int origin = graph.halfEdge(open_he).origin;
    const std::string meta_end = composeRegularStripVertexMetadata(t, voronoi_edge_id, even_half_edge_id,
      strand_even_origin_i, strand_odd_origin_i, event_type, segment_action, std::nullopt, "right", nullptr,
      "Voronoi vertex");
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
  const char* source,
  const std::optional<KineticDelaunay::CrossingData::EdgeIntersectionRef>& position_crossing,
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
  return makeBoundaryMeshFaceMetadataJson(
    kinetic_time, mesh_type, half_edge_id, delaunay_face_id, input_branch_id);
}

std::string SegmentBuilder::composeClosingMeshFaceMetadata(double kinetic_time, size_t strand_id) const
{
  if (!store_mesh_metadata)
  {
    return {};
  }
  return makeClosingMeshFaceMetadataJson(kinetic_time, strand_id);
}

void SegmentBuilder::configureMeshletStorage(VoronoiMesh& mesh) const
{
  mesh.setStoreMetadata(store_mesh_metadata);
}

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

void SegmentBuilder::markRegularMeshletCompleted(size_t meshlet_index)
{
  if (meshlet_index >= meshes.size())
  {
    return;
  }
  if (regular_meshlet_completed_.size() < meshes.size())
  {
    regular_meshlet_completed_.resize(meshes.size(), false);
  }
  regular_meshlet_completed_[meshlet_index] = true;
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

bool SegmentBuilder::warnAndSkipIfMeshletCompleted(
  const VoronoiMesh& mesh, const char* operation, double t) const
{
  if (const std::optional<size_t> regular_idx = regularMeshletIndexForMesh(mesh); regular_idx.has_value())
  {
    if (regular_idx.value() < regular_meshlet_completed_.size() && regular_meshlet_completed_[regular_idx.value()])
    {
      std::ostringstream oss;
      oss << "SegmentBuilder: skipping " << (operation != nullptr ? operation : "meshlet operation")
          << " on completed regular meshlet " << regular_idx.value();
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
    if (boundary_idx.value() < boundary_meshlet_completed_.size()
      && boundary_meshlet_completed_[boundary_idx.value()])
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
    KINDS_WARNING("neighborIntersectionOnTargetAlongVoronoiEdge: voronoi_ref stale for voronoi_edge="
      << voronoi_edge_id << ".");
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

std::optional<size_t> SegmentBuilder::radiusTransitionSharedCornerSiteVertex(
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
  const size_t he0 = 2 * s0;
  const size_t he1 = 2 * s1;
  if (he0 + 1 >= graph.halfEdgeSlotCount() || he1 + 1 >= graph.halfEdgeSlotCount())
  {
    return std::nullopt;
  }
  const int verts0[2] = { graph.halfEdge(he0).origin, graph.halfEdge(he0 + 1).origin };
  const int verts1[2] = { graph.halfEdge(he1).origin, graph.halfEdge(he1 + 1).origin };
  for (int v0 : verts0)
  {
    if (v0 < 0)
    {
      continue;
    }
    for (int v1 : verts1)
    {
      if (v1 == v0)
      {
        return static_cast<size_t>(v0);
      }
    }
  }
  return std::nullopt;
}

bool SegmentBuilder::isRadiusTransitionCornerAdjacentSourceIntersection(
  KineticDelaunay::CrossingData::EdgeIntersectionRef ref,
  const RadiusBoundaryTransitionShiftContext* boundary_transition_shift) const
{
  if (!radius_boundary_transition_shift_enabled || boundary_transition_shift == nullptr
    || !boundary_transition_shift->roles_valid
    || !boundary_transition_shift->interior_voronoi_vertex_id.has_value())
  {
    return false;
  }
  const size_t d_edge = ref->delaunay_edge_id;
  const size_t s0 = boundary_transition_shift->source_delaunay_edges[0];
  const size_t s1 = boundary_transition_shift->source_delaunay_edges[1];
  if (d_edge != s0 && d_edge != s1)
  {
    return false;
  }
  const auto shared_opt = radiusTransitionSharedCornerSiteVertex(boundary_transition_shift);
  if (!shared_opt.has_value())
  {
    return false;
  }
  const size_t shared_corner = shared_opt.value();
  const auto& graph = kin_del.getGraph();
  const size_t he_even = 2 * d_edge;
  if (he_even + 1 >= graph.halfEdgeSlotCount())
  {
    return false;
  }
  const std::vector<KineticDelaunay::CrossingData::EdgeIntersectionRef> refs
    = getBoundaryIntersectionsInBoundaryOrder(d_edge);
  if (refs.empty())
  {
    return false;
  }
  const int even_origin = graph.halfEdge(he_even).origin;
  const int odd_origin = graph.halfEdge(he_even + 1).origin;
  if (even_origin >= 0 && static_cast<size_t>(even_origin) == shared_corner)
  {
    return ref == refs.front();
  }
  if (odd_origin >= 0 && static_cast<size_t>(odd_origin) == shared_corner)
  {
    return ref == refs.back();
  }
  return false;
}

bool SegmentBuilder::voronoiEdgeHasEndpointFace(size_t voronoi_edge_id, size_t voronoi_vertex_id) const
{
  const auto& graph = kin_del.getGraph();
  const size_t he_even = 2 * voronoi_edge_id;
  const size_t he_odd = he_even + 1;
  if (he_odd >= graph.halfEdgeSlotCount())
  {
    return false;
  }
  return graph.halfEdge(he_even).face == voronoi_vertex_id
    || graph.halfEdge(he_odd).face == voronoi_vertex_id;
}

std::optional<KineticDelaunay::CrossingData::EdgeIntersectionRef>
SegmentBuilder::findInteriorVvAnchorCrossingOnTargetEdge(
  const RadiusBoundaryTransitionShiftContext* boundary_transition_shift) const
{
  if (!radius_boundary_transition_shift_enabled || boundary_transition_shift == nullptr
    || !boundary_transition_shift->roles_valid
    || !boundary_transition_shift->interior_voronoi_vertex_id.has_value())
  {
    return std::nullopt;
  }
  const size_t target_d_edge = boundary_transition_shift->target_delaunay_edge;
  const size_t interior_vv = boundary_transition_shift->interior_voronoi_vertex_id.value();
  const auto& crossing_data = kin_del.getCrossingData();
  if (target_d_edge >= crossing_data.delaunay_edge_intersections.size())
  {
    return std::nullopt;
  }

  std::optional<KineticDelaunay::CrossingData::EdgeIntersectionRef> anchor;
  for (const auto& ref : crossing_data.delaunay_edge_intersections[target_d_edge])
  {
    if (!voronoiEdgeHasEndpointFace(ref->voronoi_edge_id, interior_vv))
    {
      continue;
    }
    if (anchor.has_value())
    {
      KINDS_WARNING("findInteriorVvAnchorCrossingOnTargetEdge: multiple crossings on target delaunay_edge="
        << target_d_edge << " with Voronoi edge incident to interior_vv=" << interior_vv << ".");
    }
    anchor = ref;
  }
  return anchor;
}

std::optional<size_t> SegmentBuilder::interiorVvFictionalCornerSiteForTargetEdgeCrossing(
  KineticDelaunay::CrossingData::EdgeIntersectionRef ref,
  KineticDelaunay::CrossingData::EdgeIntersectionRef anchor_ref,
  const RadiusBoundaryTransitionShiftContext* boundary_transition_shift) const
{
  if (!radius_boundary_transition_shift_enabled || boundary_transition_shift == nullptr
    || !boundary_transition_shift->roles_valid
    || !boundary_transition_shift->interior_voronoi_vertex_id.has_value())
  {
    return std::nullopt;
  }
  const size_t target_d_edge = boundary_transition_shift->target_delaunay_edge;
  if (ref->delaunay_edge_id != target_d_edge || anchor_ref->delaunay_edge_id != target_d_edge || ref == anchor_ref)
  {
    return std::nullopt;
  }

  const std::vector<KineticDelaunay::CrossingData::EdgeIntersectionRef> refs
    = getBoundaryIntersectionsInBoundaryOrder(target_d_edge);
  size_t ref_idx = static_cast<size_t>(-1);
  size_t anchor_idx = static_cast<size_t>(-1);
  for (size_t i = 0; i < refs.size(); ++i)
  {
    if (refs[i] == ref)
    {
      ref_idx = i;
    }
    if (refs[i] == anchor_ref)
    {
      anchor_idx = i;
    }
  }
  if (ref_idx == static_cast<size_t>(-1) || anchor_idx == static_cast<size_t>(-1))
  {
    return std::nullopt;
  }

  const auto& graph = kin_del.getGraph();
  const size_t he_even = 2 * target_d_edge;
  if (he_even + 1 >= graph.halfEdgeSlotCount())
  {
    return std::nullopt;
  }
  const int even_origin = graph.halfEdge(he_even).origin;
  const int odd_origin = graph.halfEdge(he_even + 1).origin;
  if (even_origin < 0 || odd_origin < 0)
  {
    return std::nullopt;
  }

  // Boundary order is anchored at the even-half-edge origin (see corner wedge / one-null link logic).
  if (ref_idx < anchor_idx)
  {
    return static_cast<size_t>(even_origin);
  }
  return static_cast<size_t>(odd_origin);
}

std::optional<glm::dvec3> SegmentBuilder::interiorVvShiftAlongVoronoiEdgeToCornerLine(double t,
  KineticDelaunay::CrossingData::EdgeIntersectionRef ref, const glm::dvec2& vv_xy, size_t corner_site_id) const
{
  const glm::dvec2 corner_xy = kin_del.getPointAt(t, corner_site_id, false, false);
  const size_t v_edge = ref->voronoi_edge_id;
  const glm::dvec3 vor_left3 = kin_del.computeVoronoiVertexClampedInfinity(2 * v_edge, t, false, false);
  const glm::dvec3 vor_right3 = kin_del.computeVoronoiVertexClampedInfinity(2 * v_edge + 1, t, false, false);
  const glm::dvec2 vor_left(vor_left3.x, vor_left3.y);
  const glm::dvec2 vor_right(vor_right3.x, vor_right3.y);
  const glm::dvec2 out_xy = intersectAlongFirstSegmentWithLine2D(vor_left, vor_right, vv_xy, corner_xy);
  if (std::isfinite(out_xy.x) && std::isfinite(out_xy.y))
  {
    return glm::dvec3(out_xy, t);
  }
  return std::nullopt;
}

glm::dvec3 SegmentBuilder::crossingProfilePosition(double t,
  KineticDelaunay::CrossingData::EdgeIntersectionRef intersection_ref) const
{
  glm::dvec2 xy;
  if (tryComputeCrossingIntersectionPosition2D(kin_del, intersection_ref, t, xy, false, false))
  {
    return glm::dvec3(xy, t);
  }
  return closingMeshVoronoiDelaunayCrossingPosition(t, intersection_ref->voronoi_edge_id, intersection_ref->delaunay_edge_id);
}

glm::dvec3 SegmentBuilder::crossingProfilePositionFromPlacement(
  double t, const RadiusBoundaryTransitionCrossingPlacement& placement) const
{
  if (placement.explicit_profile_position.has_value())
  {
    return placement.explicit_profile_position.value();
  }
  return crossingProfilePosition(t, placement.position_intersection);
}

RadiusBoundaryTransitionCrossingPlacement SegmentBuilder::resolveRadiusBoundaryTransitionCrossingPlacement(double t,
  KineticDelaunay::CrossingData::EdgeIntersectionRef conceptual_ref,
  const RadiusBoundaryTransitionShiftContext* boundary_transition_shift) const
{
  RadiusBoundaryTransitionCrossingPlacement placement { conceptual_ref, conceptual_ref, std::nullopt };

  if (radius_boundary_transition_shift_enabled && boundary_transition_shift != nullptr
    && boundary_transition_shift->roles_valid
    && boundary_transition_shift->interior_voronoi_vertex_id.has_value())
  {
    const size_t interior_vv = boundary_transition_shift->interior_voronoi_vertex_id.value();
    const glm::dvec3 vv_pos = computeVoronoiVertex(kin_del.getGraph().face(interior_vv).half_edges[0], t);
    const glm::dvec2 vv_xy(vv_pos.x, vv_pos.y);
    const auto anchor_on_target_opt = findInteriorVvAnchorCrossingOnTargetEdge(boundary_transition_shift);

    if (isRadiusTransitionCornerAdjacentSourceIntersection(conceptual_ref, boundary_transition_shift))
    {
      KINDS_DEBUG("Radius boundary transition interior-vv corner crossing: d=" << conceptual_ref->delaunay_edge_id
        << " vv=" << interior_vv << " t=" << t << " out=(" << vv_pos.x << "," << vv_pos.y << ")");
      placement.explicit_profile_position = vv_pos;
      placement.snap_voronoi_vertex_id = interior_vv;
      return placement;
    }

    if (anchor_on_target_opt.has_value() && conceptual_ref == anchor_on_target_opt.value())
    {
      KINDS_DEBUG("Radius boundary transition interior-vv target anchor: d=" << conceptual_ref->delaunay_edge_id
        << " ve=" << conceptual_ref->voronoi_edge_id << " vv=" << interior_vv << " t=" << t << " out=(" << vv_pos.x
        << "," << vv_pos.y << ")");
      placement.explicit_profile_position = vv_pos;
      placement.snap_voronoi_vertex_id = interior_vv;
      return placement;
    }

    const size_t d_edge = conceptual_ref->delaunay_edge_id;
    const size_t s0 = boundary_transition_shift->source_delaunay_edges[0];
    const size_t s1 = boundary_transition_shift->source_delaunay_edges[1];
    const size_t target_d_edge = boundary_transition_shift->target_delaunay_edge;

    if (d_edge == s0 || d_edge == s1)
    {
      const auto shared_opt = radiusTransitionSharedCornerSiteVertex(boundary_transition_shift);
      if (shared_opt.has_value())
      {
        const auto& graph = kin_del.getGraph();
        if (auto other_opt = oppositeFiniteDelaunayVertexOnUndirectedEdge(graph, d_edge, shared_opt.value());
          other_opt.has_value())
        {
          if (auto shifted_opt
            = interiorVvShiftAlongVoronoiEdgeToCornerLine(t, conceptual_ref, vv_xy, other_opt.value());
            shifted_opt.has_value())
          {
            KINDS_DEBUG("Radius boundary transition interior-vv source crossing: d=" << d_edge << " vv=" << interior_vv
              << " ve=" << conceptual_ref->voronoi_edge_id << " corner_site=" << other_opt.value() << " t=" << t
              << " out=(" << shifted_opt->x << "," << shifted_opt->y << ")");
            placement.explicit_profile_position = shifted_opt.value();
            return placement;
          }
          const glm::dvec3 p_old = crossingProfilePosition(t, conceptual_ref);
          const glm::dvec2 other_xy = kin_del.getPointAt(t, other_opt.value(), false, false);
          const glm::dvec2 shared_xy = kin_del.getPointAt(t, shared_opt.value(), false, false);
          const glm::dvec2 old_xy(p_old.x, p_old.y);
          const double denom = glm::length(shared_xy - other_xy);
          if (denom > 1e-30)
          {
            const double alpha = glm::length(old_xy - other_xy) / denom;
            const glm::dvec2 fallback_xy = (1.0 - alpha) * other_xy + alpha * vv_xy;
            KINDS_DEBUG("Radius boundary transition interior-vv source crossing (parallel fallback): d=" << d_edge
              << " vv=" << interior_vv << " alpha=" << alpha << " t=" << t << " out=(" << fallback_xy.x << ","
              << fallback_xy.y << ")");
            placement.explicit_profile_position = glm::dvec3(fallback_xy, t);
            return placement;
          }
        }
      }
    }
    else if (d_edge == target_d_edge && anchor_on_target_opt.has_value())
    {
        if (auto corner_site_opt = interiorVvFictionalCornerSiteForTargetEdgeCrossing(
            conceptual_ref, anchor_on_target_opt.value(), boundary_transition_shift);
        corner_site_opt.has_value())
      {
        if (auto shifted_opt
          = interiorVvShiftAlongVoronoiEdgeToCornerLine(t, conceptual_ref, vv_xy, corner_site_opt.value());
          shifted_opt.has_value())
        {
          KINDS_DEBUG("Radius boundary transition interior-vv target crossing: d=" << d_edge << " vv=" << interior_vv
            << " ve=" << conceptual_ref->voronoi_edge_id << " corner_site=" << corner_site_opt.value() << " t=" << t
            << " out=(" << shifted_opt->x << "," << shifted_opt->y << ")");
          placement.explicit_profile_position = shifted_opt.value();
          return placement;
        }
      }
    }
    return placement;
  }

  if (radius_boundary_transition_shift_enabled && boundary_transition_shift != nullptr
    && boundary_transition_shift->roles_valid
    && !boundary_transition_shift->interior_voronoi_vertex_id.has_value())
  {
    const size_t d_edge = conceptual_ref->delaunay_edge_id;
    if (d_edge == boundary_transition_shift->source_delaunay_edges[0]
      || d_edge == boundary_transition_shift->source_delaunay_edges[1])
    {
      if (auto neighbor_opt
        = neighborIntersectionOnTargetAlongVoronoiEdge(conceptual_ref, boundary_transition_shift->target_delaunay_edge);
        neighbor_opt.has_value())
      {
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
    if (delaunayUndirectedEdgeHasVertex(graph, s0, site_vertex_id)
      && delaunayUndirectedEdgeHasVertex(graph, s1, site_vertex_id))
    {
      const size_t interior_vv = boundary_transition_shift->interior_voronoi_vertex_id.value();
      const glm::dvec3 vv_pos = computeVoronoiVertex(graph.face(interior_vv).half_edges[0], t);
      KINDS_DEBUG("Radius boundary transition interior-vv corner site: cell=" << site_vertex_id << " vv=" << interior_vv
        << " strip_d=" << strip_delaunay_edge_id << " t=" << t << " out=(" << vv_pos.x << "," << vv_pos.y << ")");
      return RadiusTransitionSitePlacement { vv_pos, interior_vv };
    }
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

  const auto anchor_old_and_new = [&](size_t delaunay_edge_id) -> std::optional<std::pair<glm::dvec3, glm::dvec3>>
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
      const glm::dvec3 p_old = crossingProfilePosition(t, ref.value());
      const glm::dvec3 p_new = crossingProfilePositionFromPlacement(
        t, resolveRadiusBoundaryTransitionCrossingPlacement(t, ref.value(), boundary_transition_shift));
      return std::make_pair(p_old, p_new);
    }
    if (!delaunayUndirectedEdgeHasVertex(graph, delaunay_edge_id, site_vertex_id))
    {
      return std::nullopt;
    }
    if (auto v_opp = oppositeFiniteDelaunayVertexOnUndirectedEdge(graph, delaunay_edge_id, site_vertex_id);
      v_opp.has_value())
    {
      const glm::dvec2 xy = kin_del.getPointAt(t, v_opp.value(), false, false);
      const glm::dvec3 p(xy, t);
      return std::make_pair(p, p);
    }
    return std::nullopt;
  };

  const auto a0 = anchor_old_and_new(d_strip);
  const auto a1 = anchor_old_and_new(d_other);
  if (!a0.has_value() || !a1.has_value())
  {
    return std::nullopt;
  }

  const glm::dvec3 p0_old = a0.value().first;
  const glm::dvec3 p1_old = a1.value().first;
  const glm::dvec2 site_xy = kin_del.getPointAt(t, site_vertex_id, false, false);
  const double d0 = glm::length(glm::dvec2(p0_old.x, p0_old.y) - site_xy);
  const double d1 = glm::length(glm::dvec2(p1_old.x, p1_old.y) - site_xy);
  const double sum = d0 + d1;
  if (!(sum > 1e-30))
  {
    return std::nullopt;
  }

  const glm::dvec3 p0_new = a0.value().second;
  const glm::dvec3 p1_new = a1.value().second;

  const double w0 = d1 / sum;
  const double w1 = d0 / sum;
  const glm::dvec3 out(w0 * p0_new.x + w1 * p1_new.x, w0 * p0_new.y + w1 * p1_new.y, t);
  KINDS_DEBUG("Radius boundary transition corner site: cell="
    << site_vertex_id << " strip_d=" << d_strip << " other_d=" << d_other << " w0=" << w0 << " w1=" << w1
    << " d0=" << d0 << " d1=" << d1 << " t=" << t << " out=(" << out.x << "," << out.y << ")");
  return RadiusTransitionSitePlacement { out, std::nullopt };
}

size_t SegmentBuilder::resolveIntersectionMeshPairIndex(size_t voronoi_cell_id,
  std::optional<KineticDelaunay::CrossingData::EdgeIntersectionRef> start_intersection,
  std::optional<KineticDelaunay::CrossingData::EdgeIntersectionRef> end_intersection, double event_time) const
{
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

    if (start_value->next_segment_mesh_pair_index == end_value->prev_segment_mesh_pair_index)
    {
      return start_value->next_segment_mesh_pair_index;
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
      return start_value->prev_segment_mesh_pair_index;
    }

    // One-sided link after radius boundary reset or stale crossing cleanup: trust the populated side.
    if (start_value->next_segment_mesh_pair_index != static_cast<size_t>(-1)
      && end_value->prev_segment_mesh_pair_index == static_cast<size_t>(-1))
    {
      KINDS_WARNING("resolveIntersectionMeshPairIndex: end prev unset; recovering from start next (start_next="
        << start_value->next_segment_mesh_pair_index << ", voronoi_cell_id=" << voronoi_cell_id
        << ", event_time=" << event_time << ").");
      return start_value->next_segment_mesh_pair_index;
    }
    if (end_value->prev_segment_mesh_pair_index != static_cast<size_t>(-1)
      && start_value->next_segment_mesh_pair_index == static_cast<size_t>(-1))
    {
      KINDS_WARNING("resolveIntersectionMeshPairIndex: start next unset; recovering from end prev (end_prev="
        << end_value->prev_segment_mesh_pair_index << ", voronoi_cell_id=" << voronoi_cell_id
        << ", event_time=" << event_time << ").");
      return end_value->prev_segment_mesh_pair_index;
    }

    std::ostringstream oss;
    oss << "resolveIntersectionMeshPairIndex: start/end intersection mesh pair index mismatch (start_next="
        << start_value->next_segment_mesh_pair_index << ", end_prev=" << end_value->prev_segment_mesh_pair_index
        << ", voronoi_cell_id=" << voronoi_cell_id << ", event_time=" << event_time << ").";
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

    return start_intersection.value()->next_segment_mesh_pair_index;
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

    return end_intersection.value()->prev_segment_mesh_pair_index;
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
    VoronoiMesh mesh_local(MeshletExportMaterialNames);
  configureMeshletStorage(mesh_local);
    VoronoiMesh& mesh = reuse_in_place ? intersection_meshes[intersection_pair_index] : mesh_local;

    const std::string base_boundary_meta = composeBoundaryMetadata(event_type, segment_action);

    struct CrossingOrSiteEndpoint
    {
      glm::dvec3 position {};
      std::optional<RadiusBoundaryTransitionCrossingPlacement> placement {};
      std::optional<size_t> snap_voronoi_vertex_id {};
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
        return { pos, placement, placement.snap_voronoi_vertex_id };
      }
      if (auto site_shifted
        = radiusTransitionInterpolatedSitePosition(t, voronoi_cell_id, delaunay_edge_id, boundary_transition_shift);
        site_shifted.has_value())
      {
        const glm::dvec2 p_site = kin_del.getPointAt(t, voronoi_cell_id, false, false);
        KINDS_DEBUG("Radius boundary transition [startNewMeshFromIntersections_site]: cell="
          << voronoi_cell_id << " strip_d=" << delaunay_edge_id << " t=" << t << " old=(" << p_site.x << "," << p_site.y
          << ") new=(" << site_shifted->position.x << "," << site_shifted->position.y << ")");
        return { site_shifted->position, std::nullopt, site_shifted->snap_voronoi_vertex_id };
      }
      const glm::dvec2 p = kin_del.getPointAt(t, voronoi_cell_id, false, false);
      return { glm::dvec3(p.x, p.y, t), std::nullopt, std::nullopt };
    };

    const CrossingOrSiteEndpoint start_endpoint = crossing_or_site_endpoint(interval_start_crossing);
    const CrossingOrSiteEndpoint end_endpoint = crossing_or_site_endpoint(interval_end_crossing);
    const glm::dvec3& start_pos = start_endpoint.position;
    const glm::dvec3& end_pos = end_endpoint.position;
    const auto& start_placement_opt = start_endpoint.placement;
    const auto& end_placement_opt = end_endpoint.placement;

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
      start_placement_opt.has_value()
        ? MeshletVertexRuntimeInfo { false, start_placement_opt.value().explicit_profile_position.has_value(),
            start_placement_opt.value().position_intersection, start_placement_opt.value().conceptual_intersection }
        : MeshletVertexRuntimeInfo {});
    const size_t end_vertex_index = same_endpoint
      ? start_vertex_index
      : addMeshletVertex(mesh, boundary_polygon, centroid, end_pos, voronoi_cell_id, t, false,
          end_endpoint.snap_voronoi_vertex_id, boundary_start_meta_right, glm::dvec3(0.0, 0.0, 1.0),
          end_placement_opt.has_value()
            ? MeshletVertexRuntimeInfo { false, end_placement_opt.value().explicit_profile_position.has_value(),
                end_placement_opt.value().position_intersection, end_placement_opt.value().conceptual_intersection }
            : MeshletVertexRuntimeInfo {});

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
      intersection_meshes.push_back(std::move(mesh_local));
      boundary_meshlet_completed_.push_back(false);
      intersection_meshlet_export_suffixes.push_back(std::string("_intersection_d") + std::to_string(delaunay_edge_id));
      intersection_mesh_pair_last_left_and_right_vertex.emplace_back(std::move(local_segments));
      intersection_meshes.back().setCreationKineticTime(t);
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
  };

  auto crossing_or_site_endpoint
    = [&](const std::optional<KineticDelaunay::CrossingData::EdgeIntersectionRef>& input_ref) -> CrossingOrSiteEndpoint
  {
    if (!input_ref.has_value())
    {
      if (auto site_shifted = radiusTransitionInterpolatedSitePosition(
            t, voronoi_cell_id, strip_delaunay_edge_id_for_site, boundary_transition_shift);
        site_shifted.has_value())
      {
        const glm::dvec2 p_site = kin_del.getPointAt(t, voronoi_cell_id, false, false);
        KINDS_DEBUG("Radius boundary transition [finishMeshFromIntersections_site]: cell="
          << voronoi_cell_id << " strip_d=" << strip_delaunay_edge_id_for_site << " t=" << t << " old=(" << p_site.x
          << "," << p_site.y << ") new=(" << site_shifted->position.x << "," << site_shifted->position.y << ")");
        return { site_shifted->position, std::nullopt, site_shifted->snap_voronoi_vertex_id };
      }
      const glm::dvec2 p = kin_del.getPointAt(t, voronoi_cell_id, false, false);
      return { glm::dvec3(p.x, p.y, t), std::nullopt, std::nullopt };
    }

    const RadiusBoundaryTransitionCrossingPlacement placement
      = resolveRadiusBoundaryTransitionCrossingPlacement(t, input_ref.value(), boundary_transition_shift);
    const glm::dvec3 old_pos = crossingProfilePosition(t, placement.conceptual_intersection);
    const glm::dvec3 pos = crossingProfilePositionFromPlacement(t, placement);
    if (placement.positionDiffersFromConceptual())
    {
      logRadiusBoundaryTransitionVertexShift("finishMeshFromIntersections_interval", t, placement.conceptual_intersection,
        placement.position_intersection, old_pos, pos);
    }
    return { pos, placement, placement.snap_voronoi_vertex_id };
  };

  const CrossingOrSiteEndpoint start_endpoint = crossing_or_site_endpoint(start_intersection);
  const CrossingOrSiteEndpoint end_endpoint = crossing_or_site_endpoint(end_intersection);
  const glm::dvec3 new_start_pos = start_endpoint.position;
  const glm::dvec3 new_end_pos = end_endpoint.position;
  start_placement_opt = start_endpoint.placement;
  end_placement_opt = end_endpoint.placement;
  start_snap_voronoi_vertex_id = start_endpoint.snap_voronoi_vertex_id;
  end_snap_voronoi_vertex_id = end_endpoint.snap_voronoi_vertex_id;

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
    start_placement_opt.has_value()
      ? MeshletVertexRuntimeInfo { false, start_placement_opt.value().explicit_profile_position.has_value(),
          start_placement_opt.value().position_intersection, start_placement_opt.value().conceptual_intersection }
      : MeshletVertexRuntimeInfo {});
  const size_t new_end_vertex_index = collapsed_finish_endpoints
    ? new_start_vertex_index
    : addMeshletVertex(mesh, boundary_polygon, centroid, new_end_pos, voronoi_cell_id, t, false,
        end_snap_voronoi_vertex_id, boundary_finish_meta_right, glm::dvec3(0.0, 0.0, 1.0),
        end_placement_opt.has_value()
          ? MeshletVertexRuntimeInfo { false, end_placement_opt.value().explicit_profile_position.has_value(),
              end_placement_opt.value().position_intersection, end_placement_opt.value().conceptual_intersection }
          : MeshletVertexRuntimeInfo {});

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
    KINDS_WARNING("finishMeshFromIntersections: left flexible interpolation failed for pair " << intersection_pair_index
                                                                                            << ".");
  }
  else
  {
    seg.flexible_left_vertex_ids.clear();
  }
  if (!interpolateFlexibleVerticesAlongEdge(mesh, seg.flexible_right_vertex_ids, static_cast<size_t>(old_fixed_end_id),
        uniform_finish_targets ? flex_interp_target : ordered_new_end_vertex_index))
  {
    KINDS_WARNING("finishMeshFromIntersections: right flexible interpolation failed for pair " << intersection_pair_index
                                                                                             << ".");
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
  maybeLogDiagnosticsMonitoredDelaunayEdgeTrigger(
    std::numeric_limits<double>::quiet_NaN(), "clearIntersectionMeshPairLinksOnDelaunayEdge", delaunay_edge_id,
    std::nullopt);
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

  const bool assigned_monitored = new_pair_index == kDiagnosticsMonitoredMeshPairId;
  const bool cleared_monitored
    = old_pair_index == kDiagnosticsMonitoredMeshPairId && new_pair_index != kDiagnosticsMonitoredMeshPairId;
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
  KINDS_INFO(oss.str());
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
      bool seen_start = false;
      bool seen_end = false;
      for (auto it = d_list.begin(); it != d_list.end(); ++it)
      {
        if (*it == start_intersection.value())
        {
          seen_start = true;
          break;
        }
        if (*it == end_intersection.value())
        {
          seen_end = true;
          break;
        }
      }
      if (seen_start)
      {
        assignIntersectionMeshPairLink(start_intersection.value(), false, intersection_pair_index,
          "writeIntersectionPairLinks:start_next", t);
        assignIntersectionMeshPairLink(end_intersection.value(), true, intersection_pair_index,
          "writeIntersectionPairLinks:end_prev", t);
        wrote_pair_index_to = "start_ref->next_segment_mesh_pair_index, end_ref->prev_segment_mesh_pair_index";
      }
      else if (seen_end)
      {
        assignIntersectionMeshPairLink(start_intersection.value(), true, intersection_pair_index,
          "writeIntersectionPairLinks:start_prev", t);
        assignIntersectionMeshPairLink(end_intersection.value(), false, intersection_pair_index,
          "writeIntersectionPairLinks:end_next", t);
        wrote_pair_index_to = "start_ref->prev_segment_mesh_pair_index, end_ref->next_segment_mesh_pair_index";
      }
      else
      {
        assignIntersectionMeshPairLink(start_intersection.value(), false, intersection_pair_index,
          "writeIntersectionPairLinks:start_next_fallback", t);
        assignIntersectionMeshPairLink(end_intersection.value(), true, intersection_pair_index,
          "writeIntersectionPairLinks:end_prev_fallback", t);
        wrote_pair_index_to
          = "start_ref->next_segment_mesh_pair_index, end_ref->prev_segment_mesh_pair_index (list lookup fallback)";
      }
    }
    else
    {
      assignIntersectionMeshPairLink(start_intersection.value(), false, intersection_pair_index,
        "writeIntersectionPairLinks:start_next_missing_list", t);
      assignIntersectionMeshPairLink(end_intersection.value(), true, intersection_pair_index,
        "writeIntersectionPairLinks:end_prev_missing_list", t);
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
  maybeLogDiagnosticsMonitoredDelaunayEdgeTrigger(t, "writeIntersectionPairLinks", trigger_d_edge, intersection_pair_index);
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
  regular_meshlet_completed_.push_back(false);
  meshlet_export_suffixes.push_back(std::move(suffix));
  return index;
}

void kinDS::SegmentBuilder::completeBoundaryMeshSection(size_t he_id, size_t new_left, size_t new_right, double t)
{
  const std::string face_metadata
    = composeBoundaryMeshFaceMetadata(t, "boundary_section", he_id);
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
      addBoundaryTriangle(last_left_and_right.second, new_right, half_edge_to_boundary_vertex_index[he_id], face_metadata);
      addBoundaryTriangle(new_left, last_left_and_right.first, half_edge_to_boundary_vertex_index[he_id], face_metadata);
      addBoundaryTriangle(new_left, half_edge_to_boundary_vertex_index[he_id], new_right, face_metadata);

      // reset the half-edge to boundary vertex index
      half_edge_to_boundary_vertex_index[he_id] = -1;
    }
  }
  else
  {
    assert(last_left_and_right.second == -1);
  }
}

size_t kinDS::SegmentBuilder::addBoundaryTriangle(size_t u, size_t v, size_t w, const std::string& metadata)
{
  // check bounds
  if (u >= boundary_mesh.getVertexCount() || v >= boundary_mesh.getVertexCount() || w >= boundary_mesh.getVertexCount())
  {
    KINDS_ERROR("Vertex index out of boundary mesh range.");
    return -1;
  }

  if (u >= boundary_mesh_raw_uvs.size() || v >= boundary_mesh_raw_uvs.size() || w >= boundary_mesh_raw_uvs.size())
  {
    KINDS_ERROR("Vertex index out of raw uv range.");
    return -1;
  }

  // get raw UVs
  glm::dvec3 uv_u = { boundary_mesh_raw_uvs[u][0], boundary_mesh_raw_uvs[u][1], 0.0 };
  glm::dvec3 uv_v = { boundary_mesh_raw_uvs[v][0], boundary_mesh_raw_uvs[v][1], 0.0 };
  glm::dvec3 uv_w = { boundary_mesh_raw_uvs[w][0], boundary_mesh_raw_uvs[w][1], 0.0 };

  // output UVs
  /*KINDS_DEBUG("Adding boundary triangle with raw UVs: u(" + std::to_string(uv_u[0]) + ", " + std::to_string(uv_u[1])
     +
                "), v(" + std::to_string(uv_v[0]) + ", " + std::to_string(uv_v[1]) + "), w(" + std::to_string(uv_w[0]) +
                ", " + std::to_string(uv_w[1]) + ")");*/

  // adjust UVs to avoid seams, first coordinate is the angle normalized to [-0.5, 0.5]
  // As a heuristic, we take the first angle and adjust the others such that they have less than 0.5 difference
  double base_angle = uv_u[0];
  double& angle_v = uv_v[0];
  double diff_v = angle_v - base_angle;
  double adjustment = std::round(diff_v);
  angle_v -= adjustment;

  double& angle_w = uv_w[0];
  double diff_w = angle_w - base_angle;
  adjustment = std::round(diff_w);
  angle_w -= adjustment;

  uv_u[0] *= uv_circum_factor;
  uv_v[0] *= uv_circum_factor;
  uv_w[0] *= uv_circum_factor;
  uv_u[1] *= uv_height_factor;
  uv_v[1] *= uv_height_factor;
  uv_w[1] *= uv_height_factor;

  // add adjusted UVs
  size_t uv_index_u = boundary_mesh.addUV(uv_u);
  size_t uv_index_v = boundary_mesh.addUV(uv_v);
  size_t uv_index_w = boundary_mesh.addUV(uv_w);

  /*KINDS_DEBUG("UVs after adjustment: u(" + std::to_string(uv_u[0]) + ", " + std::to_string(uv_u[1]) + "), v(" +
                std::to_string(uv_v[0]) + ", " + std::to_string(uv_v[1]) + "), w(" + std::to_string(uv_w[0]) + ", " +
                std::to_string(uv_w[1]) + ")");*/
  // warnIfTriangleKineticTimesNotInUnitSection(u, v, w, boundary_mesh.getVertices(), "boundary_mesh", 0);
  return boundary_mesh.addTriangle(u, v, w, uv_index_u, uv_index_v, uv_index_w, 0, metadata);
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
  double angle = std::atan2(centroid.y - profile_xy.y, centroid.x - profile_xy.x);

  glm::dvec2 raw_uv { angle / (2.0 * glm::pi<double>()), vertex[2] };

  const std::string metadata = store_mesh_metadata
    ? [&]()
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
  return mesh.addTriangle(u, v, w, u, v, w, material_id, metadata);
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
  if (inside_boundary_he_id < 0
    || static_cast<size_t>(inside_boundary_he_id) >= kin_del.getGraph().halfEdgeSlotCount())
  {
    return addMeshletTriangle(mesh, u, v, w, triangle_metadata, boundary_material_id);
  }

  // `inside_boundary_he_id` is the inside-directed boundary half-edge; its twin is the outside one on the same Delaunay
  // edge.
  const size_t outside_he = static_cast<size_t>(inside_boundary_he_id) ^ 1u;
  if (outside_he >= kin_del.getGraph().halfEdgeSlotCount())
  {
    return addMeshletTriangle(mesh, u, v, w, triangle_metadata, boundary_material_id);
  }

  if ((outside_he & 1u) != 0u)
  {
    std::swap(v, w);
  }
  return addMeshletTriangle(mesh, u, v, w, triangle_metadata, boundary_material_id);
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
  const char* context, size_t voronoi_vertex_id, const glm::dvec3& position) const
{
  const size_t containing_tri_id = kin_del.getCrossingDataContainingTriId(voronoi_vertex_id);
  if (kin_del.getFaceInside(containing_tri_id))
  {
    return;
  }
  KINDS_WARNING("SegmentBuilder: " << context << " - Voronoi vertex " << voronoi_vertex_id
                                   << " (containing Delaunay triangle " << containing_tri_id
                                   << ") is outside the alpha-shape; position (" << position.x << ", " << position.y
                                   << ", " << position.z << ").");
}

glm::dvec3 SegmentBuilder::transformFromInputBranchToObjectSpace(
  glm::dvec3 vertex, size_t strand_id, double t) const
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

glm::dvec3 SegmentBuilder::computeMeshSiteVertexPosition(
  glm::dvec3 profile_vertex, size_t strand_id, double t) const
{
  (void)profile_vertex;
  return getPointInMeshSpace(strand_id, t);
}

SegmentBuilder::MeshIntersectionObjectSpaceResult SegmentBuilder::computeMeshIntersectionObjectSpace(
  const MeshletVertexRuntimeInfo& runtime_info, glm::dvec3 fallback_profile_vertex, size_t fallback_strand_id,
  double t) const
{
  MeshIntersectionObjectSpaceResult result;
  result.position = computeMeshSiteVertexPosition(fallback_profile_vertex, fallback_strand_id, t);

  if (!runtime_info.position_intersection.has_value())
  {
    return result;
  }
  const auto ref = runtime_info.position_intersection.value();
  const size_t delaunay_edge_id = ref->delaunay_edge_id;

  const auto& graph = kin_del.getGraph();
  const size_t d_he0 = 2 * delaunay_edge_id;
  const size_t d_he1 = d_he0 + 1;
  if (d_he1 >= graph.halfEdgeSlotCount())
  {
    return result;
  }

  const int a = graph.halfEdge(d_he0).origin;
  const int b = graph.halfEdge(d_he1).origin;
  if (a < 0 || b < 0)
  {
    return result;
  }

  // Placement must be consistent with the live intersection ref updated by Crossing/Radius events.
  const double param = ref->delaunay_edge_param;
  if (!std::isfinite(param))
  {
    return result;
  }

  const size_t strand_a = static_cast<size_t>(a);
  const size_t strand_b = static_cast<size_t>(b);
  const glm::dvec2 profile_a = kin_del.getStrandTree().evaluate(strand_a, t);
  const glm::dvec2 profile_b = kin_del.getStrandTree().evaluate(strand_b, t);
  glm::dvec3 a_profile(profile_a.x, profile_a.y, t);
  glm::dvec3 b_profile(profile_b.x, profile_b.y, t);
  const glm::dvec3 a_mesh = computeMeshSiteVertexPosition(a_profile, strand_a, t);
  const glm::dvec3 b_mesh = computeMeshSiteVertexPosition(b_profile, strand_b, t);
  result.position = a_mesh * (1.0 - param) + b_mesh * param;
  result.mesh_interpolation = IntersectionInterpolationDebug { a_mesh, b_mesh, param };

  return result;
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

  // Containing Delaunay triangle (may differ from the dual triangle of this Voronoi vertex).
  const std::optional<size_t> containing_tri_opt
    = kin_del.getCrossingData().peekContainingTriId(voronoi_vertex_id);
  if (!containing_tri_opt.has_value())
  {
    KINDS_WARNING("computeMeshVoronoiVertexObjectSpace: Voronoi vertex " << voronoi_vertex_id
      << " has no registered containing triangle at t=" << t << "; falling back.");
    return result;
  }
  const size_t containing_tri_id = containing_tri_opt.value();
  result.containing_tri_id = containing_tri_id;
  if (containing_tri_id >= graph.faceSlotCount() || !graph.isLiveFace(containing_tri_id))
  {
    KINDS_WARNING("computeMeshVoronoiVertexObjectSpace: containing triangle " << containing_tri_id
      << " for Voronoi vertex " << voronoi_vertex_id << " is not live at t=" << t << "; falling back.");
    return result;
  }

  const std::array<int, 3> containing_vertices = graph.getTriangleVertexIndices(containing_tri_id);
  std::array<size_t, 3> site_ids {};
  for (size_t i = 0; i < containing_vertices.size(); ++i)
  {
    if (containing_vertices[i] < 0)
    {
      KINDS_WARNING("computeMeshVoronoiVertexObjectSpace: containing triangle " << containing_tri_id
        << " for Voronoi vertex " << voronoi_vertex_id << " is infinite at t=" << t << "; falling back.");
      return result;
    }
    site_ids[i] = static_cast<size_t>(containing_vertices[i]);
  }

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
    KINDS_WARNING("computeMeshVoronoiVertexObjectSpace: degenerate containing triangle " << containing_tri_id
      << " for Voronoi vertex " << voronoi_vertex_id << " at t=" << t << "; falling back.");
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
      KINDS_WARNING("addMeshletVertex(" << stage << "): non-finite vertex for strand " << strand_id << " at t=" << t
                                        << " -> (" << p.x << ", " << p.y << ", " << p.z << ").");
      return;
    }
    if (p.x == 0.0 && p.y == 0.0)
    {
      KINDS_WARNING("addMeshletVertex(" << stage << "): degenerate XY vertex (0,0,z) for strand " << strand_id
                                        << " at t=" << t << " -> (" << p.x << ", " << p.y << ", " << p.z << ").");
    }
  };

  const bool is_flexible_placeholder = runtime_info.is_flexible_placeholder;
  const bool radius_shift_explicit_profile_position = runtime_info.radius_shift_explicit_profile_position;
  if (!is_flexible_placeholder)
  {
    warn_degenerate_or_non_finite(vertex, "input");
  }
  glm::dvec2 profile_xy(vertex.x, vertex.y);
  // Never bake kinetic/mesh virtual separation offsets into stored mesh vertices.
  includes_virtual_shift = false;

  const bool is_intersection_vertex = runtime_info.isIntersectionVertex();
  const bool skip_world_intersection_interp = is_flexible_placeholder || radius_shift_explicit_profile_position;
  std::optional<MeshIntersectionObjectSpaceResult> intersection_object_space;
  if (is_intersection_vertex && !skip_world_intersection_interp)
  {
    intersection_object_space
      = computeMeshIntersectionObjectSpace(runtime_info, glm::dvec3(profile_xy.x, profile_xy.y, t), strand_id, t);
  }

  std::optional<MeshVoronoiVertexObjectSpaceResult> voronoi_object_space;
  if (!is_flexible_placeholder && is_intersection_vertex && !radius_shift_explicit_profile_position
    && intersection_object_space.has_value())
  {
    vertex = intersection_object_space.value().position;
    if (!vertexPositionFinite(vertex))
    {
      KINDS_WARNING("addMeshletVertex: non-finite mesh-space intersection vertex for strand " << strand_id << " at t="
        << t << "; falling back to site-equivalent profile placement.");
      vertex = computeMeshSiteVertexPosition(glm::dvec3(profile_xy.x, profile_xy.y, t), strand_id, t);
    }
  }
  else if (!is_flexible_placeholder && meshlet_voronoi_vertex_for_alpha_check.has_value())
  {
    // Prefer barycentric Voronoi placement (also used when radius transition snaps a site/crossing to a VV).
    voronoi_object_space = computeMeshVoronoiVertexObjectSpace(
      meshlet_voronoi_vertex_for_alpha_check.value(), glm::dvec3(profile_xy.x, profile_xy.y, t), strand_id, t);
    vertex = voronoi_object_space.value().position;
    if (!vertexPositionFinite(vertex))
    {
      KINDS_WARNING("addMeshletVertex: non-finite mesh-space Voronoi vertex for strand " << strand_id << " at t=" << t
        << "; falling back to site-equivalent profile placement.");
      vertex = computeMeshSiteVertexPosition(glm::dvec3(profile_xy.x, profile_xy.y, t), strand_id, t);
    }
  }
  else if (!is_flexible_placeholder && !is_intersection_vertex)
  {
    // True sites always use unshifted mesh placement (never keep a kinetically shifted profile,
    // including radius-transition synthetic site positions).
    vertex = computeMeshSiteVertexPosition(glm::dvec3(profile_xy.x, profile_xy.y, t), strand_id, t);
    profile_xy = glm::dvec2(vertex.x, vertex.y);
    includes_virtual_shift = false;
  }
  else if (create_transformed_mesh && !is_flexible_placeholder)
  {
    vertex = transformFromInputBranchToObjectSpace(vertex, strand_id, t);
    if (!vertexPositionFinite(vertex))
    {
      KINDS_WARNING("addMeshletVertex: non-finite transformed vertex for strand " << strand_id << " at t=" << t
        << "; falling back to input-branch transform.");
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
      glm::dvec3(profile_xy.x, profile_xy.y, t));
  }
  std::string vertex_metadata = metadata;
  if (store_mesh_metadata)
  {
    const std::optional<std::string> event_type_field = metadataStringField(metadata, "event_type");
    const std::string event_type = event_type_field.value_or("unknown_event");
    if (!event_type_field.has_value())
    {
      KINDS_WARNING("addMeshletVertex: missing metadata event_type for strand " << strand_id << " at t=" << t
                                                                               << "; using 'unknown_event'.");
    }
    else if (event_type == "unknown_event")
    {
      KINDS_WARNING("addMeshletVertex: unknown metadata event_type for strand " << strand_id << " at t=" << t << ".");
    }
    // Source follows how the vertex was placed. VV snap wins over intersection/site labels.
    const std::string source = meshlet_voronoi_vertex_for_alpha_check.has_value()
      ? "Voronoi vertex"
      : (is_intersection_vertex ? "intersection" : "site");
    const std::optional<std::string> source_field = metadataStringField(metadata, "source");
    if (source_field.has_value() && source_field.value() != source)
    {
      KINDS_WARNING("addMeshletVertex: caller metadata source '" << source_field.value()
                                                                << "' disagrees with placement source '" << source
                                                                << "' for strand " << strand_id << " at t=" << t
                                                                << "; using placement source.");
    }
    else if (!source_field.has_value())
    {
      KINDS_WARNING("addMeshletVertex: missing metadata source for strand " << strand_id << " at t=" << t
                                                                            << "; using placement source '" << source
                                                                            << "'.");
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
      }
    }
    else
    {
      builder.addSize("strand_id", strand_id);
    }

    if (!metadata_callback_phase_.empty())
    {
      builder.addString("callback", metadata_callback_phase_);
    }

    vertex_metadata = builder.addBool("shift", includes_virtual_shift)
                        .addDouble("x", profile_xy.x)
                        .addDouble("y", profile_xy.y)
                        .addDouble("t", t)
                        .build();
  }
  size_t index
    = mesh.addVertex(vertex, vertex_metadata, debug_color.has_value() ? debug_color.value() : glm::dvec3(1.0));
  mesh.setVertexKineticTime(index, t);
  if (create_transformed_mesh)
  {
    mesh.setProfilePlanePosition(index, profile_xy);
  }
  double rel_dist = relativeDistanceFromCenter(boundary_polygon, centroid, profile_xy);

  /*if (rel_dist > 1.0 + std::numeric_limits<double>::epsilon()) {
    KINDS_WARNING("Adding vertex that is too far outside, relative distance: " << rel_dist);
  }*/

  // TODO: this can be simplified to not use trigonometric functions
  double angle = std::atan2(centroid.y - profile_xy.y, centroid.x - profile_xy.x);
  double u = 0.5 + texture_diameter * rel_dist * 0.5 * std::cos(angle);
  double v = 0.5 + texture_diameter * rel_dist * 0.5 * std::sin(angle);
  mesh.addUV(u, v, vertex[2] * uv_height_factor);
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
  const size_t idx = addMeshletVertex(
    mesh, boundary_polygon, centroid, placeholder, strand_id, t, false, std::nullopt, vertex_meta, std::nullopt,
    MeshletVertexRuntimeInfo { true, false, std::nullopt, std::nullopt });
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
  if (neighbor_pair_idx >= intersection_meshes.size()
    || neighbor_pair_idx >= intersection_mesh_pair_last_left_and_right_vertex.size())
  {
    return;
  }

  const size_t shared_d_edge = shared_ref->delaunay_edge_id;
  const size_t shared_he_even = 2 * shared_d_edge;
  if (shared_he_even + 1 >= kin_del.getGraph().halfEdgeSlotCount()
    || !kin_del.isOnComponentBoundary(shared_he_even))
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
    logRadiusBoundaryTransitionVertexShift("extendIntersectionMeshAtSharedCrossing", t, placement.conceptual_intersection,
      placement.position_intersection, old_pos, crossing_pos);
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
      placement.conceptual_intersection });
  addBoundaryIntervalTriangleOriented(mesh, eff_l, eff_r, new_vid, inside_boundary_he_id, t, base_meta, neighbor_pair_idx);
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
    addFlexibleVertexToIntersectionMesh(mesh, seg, !fixed_start_side, boundary_polygon, centroid, strand_id, t,
      flexible_base_metadata);
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
      KINDS_WARNING("Unresolved flexible vertices in " << (context != nullptr ? context : "unknown")
                                                       << " without both fixed anchors: left="
                                                       << seg.flexible_left_vertex_ids.size()
                                                       << " right=" << seg.flexible_right_vertex_ids.size() << ".");
    }
    return;
  }

  const size_t start_anchor = static_cast<size_t>(seg.mesh_start_vertex_id);
  const size_t end_anchor = static_cast<size_t>(seg.mesh_end_vertex_id);

  if (!seg.flexible_left_vertex_ids.empty())
  {
    KINDS_WARNING("Resolving leftover left flexible vertices in " << (context != nullptr ? context : "unknown")
                                                                  << ": count="
                                                                  << seg.flexible_left_vertex_ids.size() << ".");
    if (!interpolateFlexibleVerticesAlongEdge(
          mesh, seg.flexible_left_vertex_ids, start_anchor, end_anchor))
    {
      snapFlexibleVerticesToAnchor(mesh, seg.flexible_left_vertex_ids, start_anchor);
    }
    seg.flexible_left_vertex_ids.clear();
  }
  if (!seg.flexible_right_vertex_ids.empty())
  {
    KINDS_WARNING("Resolving leftover right flexible vertices in " << (context != nullptr ? context : "unknown")
                                                                   << ": count="
                                                                   << seg.flexible_right_vertex_ids.size() << ".");
    if (!interpolateFlexibleVerticesAlongEdge(
          mesh, seg.flexible_right_vertex_ids, end_anchor, start_anchor))
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
    if (flex_index < mesh.getVertices().size()
      && !vertexPositionFinite(mesh.getVertices()[flex_index]))
    {
      KINDS_WARNING("resolveRemainingFlexibleVertices(" << (context != nullptr ? context : "unknown")
                                                       << "): anchor vertex " << flex_index << " remains non-finite.");
    }
  }
}

void SegmentBuilder::resolveAllIntersectionFlexibleVertices(const char* context)
{
  for (size_t mesh_id = 0;
    mesh_id < intersection_meshes.size() && mesh_id < intersection_mesh_pair_last_left_and_right_vertex.size(); ++mesh_id)
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

  std::vector<double> relative_center_distances(graph.getVertexCount(), 0.0);
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
    auto& boundary_polygon = kin_del.component_data.component_boundaries[component_index][0];
    auto& centroid = kin_del.component_data.component_centroids[component_index];

    const size_t mesh_vertex_index
      = addBoundaryVertex(glm::dvec3 { vertex[0], vertex[1], t + offset }, centroid, strand_id, t, false);

    strand_to_mesh_vertex[strand_id] = static_cast<int>(mesh_vertex_index);
    relative_center_distances[strand_id] = relativeDistanceFromCenter(boundary_polygon, centroid, vertex);
    ++strands_vertex_added;

    if (diagnostics && strand_id == 0)
    {
      std::ostringstream oss;
      oss << "input_branch_id=" << input_branch_id << " boundary_mesh_vertex=" << mesh_vertex_index << " offset=" << offset
          << " invert_orientation=" << (invert_orientation ? "true" : "false");
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

    size_t uv_indices[3];

    for (size_t i = 0; i < 3; i++)
    {
      const size_t mesh_vertex_index = static_cast<size_t>(strand_to_mesh_vertex[vertices[i]]);
      double rel_dist = relative_center_distances[vertices[i]];
      double angle = boundary_mesh_raw_uvs[mesh_vertex_index][0] * 2.0 * glm::pi<double>();
      double u = 0.5 + texture_diameter * rel_dist * 0.5 * std::cos(angle);
      double v = 0.5 + texture_diameter * rel_dist * 0.5 * std::sin(angle);
      uv_indices[i] = boundary_mesh.addUV(u, v, 0.0);
    }

    const size_t tri_v0 = static_cast<size_t>(strand_to_mesh_vertex[vertices[0]]);
    const size_t tri_v1 = static_cast<size_t>(strand_to_mesh_vertex[vertices[1]]);
    const size_t tri_v2 = static_cast<size_t>(strand_to_mesh_vertex[vertices[2]]);
    const std::string face_metadata = composeBoundaryMeshFaceMetadata(
      t, "boundary_delaunay", static_cast<size_t>(-1), face_index, input_branch_id);
    boundary_mesh.addTriangle(
      tri_v0, tri_v1, tri_v2, uv_indices[0], uv_indices[1], uv_indices[2], 1, face_metadata);
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
        << " boundary_tris_before=" << boundary_tris_before << " boundary_tris_after=" << boundary_mesh.getTriangleCount();
    strandInitDiagnosticLogLine("boundary_mesh_summary", branch_contains_strand_0 ? 0 : branch_strands.empty() ? 0
                                                                                                               : branch_strands.front(),
      t, oss.str().c_str());
    if (branch_contains_strand_0 && strands_vertex_added == 0)
    {
      strandInitDiagnosticLogLine("boundary_mesh_strand0_missing", 0, t,
        "strand 0 is in branch but no boundary mesh vertex was added");
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
  if (!isComponentLive(component_index))
  {
    return;
  }

  if (kin_del.component_data.component_last_updated[component_index] != t)
  {
    kin_del.component_data.component_boundaries[component_index]
      = kin_del.extractComponentBoundaries(kin_del.component_data.components[component_index], t, visited, false, false);
    kin_del.component_data.component_centroids[component_index]
      = polygonCentroid(kin_del.component_data.component_boundaries[component_index][0]);
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
    boundary_mesh_last_left_and_right_vertex.resize(he_slots, std::make_pair(static_cast<size_t>(-1), static_cast<size_t>(-1)));
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
      boundary_mesh_last_left_and_right_vertex[he_id] = std::make_pair(static_cast<size_t>(-1), static_cast<size_t>(-1));
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
      startNewMeshFromIntersections(
        first_cell, t, std::nullopt, refs.front(), false, BoundaryEventType::Section, BoundarySegmentAction::NewSegment);
    }
    for (size_t k = 0; k + 1 < refs.size(); ++k)
    {
      const size_t mid_cell
        = determineVoronoiCellForBoundaryIntersectionInterval(d_edge_id, refs[k], refs[k + 1]);
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

void kinDS::SegmentBuilder::onBeforeComponentGraphSplit(double /*t*/)
{
}

glm::dvec3 kinDS::SegmentBuilder::closingMeshVoronoiDelaunayCrossingPosition(
  double t, size_t voronoi_edge_id, size_t delaunay_edge_id) const
{
  auto& graph = kin_del.getGraph();
  const size_t voronoi_he0 = 2 * voronoi_edge_id;
  const size_t voronoi_he1 = voronoi_he0 + 1;

  const glm::dvec3 left_vertex = computeVoronoiVertex(voronoi_he0, t);
  const glm::dvec3 right_vertex = computeVoronoiVertex(voronoi_he1, t);
  const glm::dvec2 left2(left_vertex.x, left_vertex.y);
  const glm::dvec2 right2(right_vertex.x, right_vertex.y);

  const size_t d_he0 = 2 * delaunay_edge_id;
  const size_t d_he1 = d_he0 + 1;
  const size_t d_a = graph.halfEdge(d_he0).origin;
  const size_t d_b = graph.halfEdge(d_he1).origin;
  if (d_a == size_t(-1) || d_b == size_t(-1))
  {
    KINDS_WARNING("closingMeshVoronoiDelaunayCrossingPosition: degenerate Delaunay edge for voronoi_edge="
      << voronoi_edge_id << " delaunay_edge=" << delaunay_edge_id << " t=" << t);
    return glm::dvec3(std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN(), t);
  }

  const glm::dvec2 d0 = kin_del.getPointAt(t, d_a, false, false);
  const glm::dvec2 d1 = kin_del.getPointAt(t, d_b, false, false);
  const glm::dvec2 p = intersectSegments(left2, right2, d0, d1);
  return glm::dvec3(p, t);
}

auto kinDS::SegmentBuilder::extractSegmentsForVoronoiEdge(double t, int incident_edge_index, size_t voronoi_edge_id,
  const std::function<int(const glm::dvec3&, std::optional<size_t>, const std::string&)>& track_vertex, bool reverse)
  -> std::vector<MeshingData>
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
    out_segments.push_back(MeshingData { track_vertex(strand_cm_vertex, voronoi_vertex_id,
                                      makeClosingMeshVertexMetadata("Voronoi vertex", std::nullopt, voronoi_vertex_id)),
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
        makeClosingMeshVertexMetadata("Voronoi vertex", std::nullopt, voronoi_vertex_id));
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
      const glm::dvec3 crossing_pos
        = closingMeshVoronoiDelaunayCrossingPosition(t, voronoi_edge_id, crossing_ref->delaunay_edge_id);
      if (!std::isfinite(crossing_pos.x) || !std::isfinite(crossing_pos.y))
      {
        continue;
      }

      const int mesh_vertex_id = track_vertex(
        crossing_pos, std::nullopt, makeClosingMeshVertexMetadata("intersection", crossing_ref));

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
    out_segments.back().mesh_end_vertex_id = track_vertex(
      other_cm_vertex, voronoi_vertex_id, makeClosingMeshVertexMetadata("Voronoi vertex", std::nullopt, voronoi_vertex_id));
  }

  return out_segments;
}

auto kinDS::SegmentBuilder::closingMeshExtractRawSegmentsForVoronoiEdge(size_t strand_id, double t, VoronoiMesh& mesh,
  const std::vector<BoundaryPoint>& boundary_polygon, const glm::dvec2& centroid, int incident_edge_index,
  size_t incident_he,
  const std::function<int(const glm::dvec3&, std::optional<size_t>, const std::string&)>& track_vertex)
  -> std::vector<MeshingData>
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
    if (s->start_crossing.has_value()
      && tryComputeCrossingIntersectionPosition2D(kin_del, s->start_crossing, t, ip, false, false))
    {
      if (glm::distance(ip, ps) > k_closing_cap_geom_eps)
      {
        KINDS_ERROR("Closing mesh ordered_seg[" << si << "]: start_ref 2D position does not match mesh_start.");
      }
    }
    if (s->end_crossing.has_value()
      && tryComputeCrossingIntersectionPosition2D(kin_del, s->end_crossing, t, ip, false, false))
    {
      if (glm::distance(ip, pe) > k_closing_cap_geom_eps)
      {
        KINDS_ERROR("Closing mesh ordered_seg[" << si << "]: end_ref 2D position does not match mesh_end.");
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

  auto& cap_vertex_ids = result.mesh_vertex_ids;
  auto add_mesh_vertex = [&](const glm::dvec3& v, const std::string& metadata) -> int
  {
    const int id = closingMeshAppendVertex(mesh, boundary_polygon, centroid, strand_id, t, v, false, std::nullopt,
      metadata);
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

      polygon_push(polygon, static_cast<size_t>(current_segment->mesh_start_vertex_id));
      logClosingMeshStripEndpointVertex(
        kin_del, current_segment_index, current_segment->start_crossing, current_segment->closing_voronoi_edge_id);
      polygon_push(polygon, static_cast<size_t>(current_segment->mesh_end_vertex_id));
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
        const int vid = add_mesh_vertex(strand_pos, makeClosingMeshVertexMetadata("site", std::nullopt, std::nullopt,
                                                strand_id));
        polygon_push(polygon, static_cast<size_t>(vid));
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
              const glm::dvec3 inter_pos = closingMeshVoronoiDelaunayCrossingPosition(
                t, candidate_ref->voronoi_edge_id, candidate_ref->delaunay_edge_id);
              if (std::isfinite(inter_pos.x) && std::isfinite(inter_pos.y))
              {
                const int nv
                  = add_mesh_vertex(inter_pos, makeClosingMeshVertexMetadata("intersection", candidate_ref));
                polygon_push(polygon, static_cast<size_t>(nv));
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
              makeClosingMeshVertexMetadata("site", std::nullopt, std::nullopt, static_cast<size_t>(corner_vertex_id)));
            polygon_push(polygon, static_cast<size_t>(cv));
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
        glm::dvec2 xy(0.0, 0.0);
        const bool have_xy = tryComputeCrossingIntersectionPosition2D(kin_del, walk_closure_crossing_ref, t, xy, false, false);
        const auto& ir = *walk_closure_crossing_ref.value();
        const glm::dvec3 geom = closingMeshVoronoiDelaunayCrossingPosition(t, ir.voronoi_edge_id, ir.delaunay_edge_id);
        head << " Loop closed at crossing " << formatCrossingIntersectionForLog(kin_del, walk_closure_crossing_ref)
             << ".";
        if (have_xy)
        {
          head << " CrossingData 2D=(" << xy.x << "," << xy.y << ")";
        }
        else
        {
          head << " CrossingData 2D=(unavailable)";
        }
        if (std::isfinite(geom.x) && std::isfinite(geom.y))
        {
          head << "; Voronoi-Delaunay chord intersection=(" << geom.x << "," << geom.y << ").";
        }
        else
        {
          head << "; Voronoi-Delaunay chord intersection=(non-finite).";
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
      double area2 = 0.0;
      for (size_t i = 0; i < polygon.size(); ++i)
      {
        const glm::dvec3& p0 = mesh.getVertices()[polygon[i]];
        const glm::dvec3& p1 = mesh.getVertices()[polygon[(i + 1) % polygon.size()]];
        area2 += p0.x * p1.y - p1.x * p0.y;
      }
      if (area2 < 0.0)
      {
        std::reverse(polygon.begin() + 1, polygon.end());
      }
      result.polygons.push_back(std::move(polygon));
    }
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
  const std::string& metadata, int material_id, bool orient_upwards)
{
  constexpr double eps = 1e-12;
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
    writeTriangulateSimplePolygonDebugTxt(kin_del, mesh,
      makeTriangulateSimplePolygonDebugPath(kin_del, mesh, "SPLIT", ".txt"), "SPLIT", debug_rings);

    for (const std::vector<size_t>& sub_polygon : split_polygons)
    {
      triangulateSimplePolygon(mesh, sub_polygon, metadata, material_id, orient_upwards);
    }
    return;
  }
  vertices = split_polygons.front();

  auto cross_at = [&](size_t prev, size_t current, size_t next)
  {
    const glm::dvec2 a = mesh.triangulationPlaneXY(prev);
    const glm::dvec2 b = mesh.triangulationPlaneXY(current);
    const glm::dvec2 c = mesh.triangulationPlaneXY(next);
    return glm::cross(b - a, c - b);
  };

  bool removed_collinear = true;
  while (removed_collinear && vertices.size() > 3)
  {
    removed_collinear = false;
    for (size_t i = 0; i < vertices.size(); ++i)
    {
      const size_t prev = vertices[(i + vertices.size() - 1) % vertices.size()];
      const size_t current = vertices[i];
      const size_t next = vertices[(i + 1) % vertices.size()];
      if (std::abs(cross_at(prev, current, next)) <= eps)
      {
        vertices.erase(vertices.begin() + static_cast<std::ptrdiff_t>(i));
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
    for (size_t i = 0; i < vertices.size(); ++i)
    {
      const glm::dvec2 p0 = mesh.triangulationPlaneXY(vertices[i]);
      const glm::dvec2 p1 = mesh.triangulationPlaneXY(vertices[(i + 1) % vertices.size()]);
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
  }

  auto point_in_triangle = [&](size_t point_id, size_t a_id, size_t b_id, size_t c_id)
  {
    const glm::dvec2 p = mesh.triangulationPlaneXY(point_id);
    const glm::dvec2 a = mesh.triangulationPlaneXY(a_id);
    const glm::dvec2 b = mesh.triangulationPlaneXY(b_id);
    const glm::dvec2 c = mesh.triangulationPlaneXY(c_id);
    const double ab = glm::cross(b - a, p - a);
    const double bc = glm::cross(c - b, p - b);
    const double ca = glm::cross(a - c, p - c);
    return ab >= -eps && bc >= -eps && ca >= -eps;
  };

  while (vertices.size() > 3)
  {
    bool clipped_ear = false;
    for (size_t i = 0; i < vertices.size(); ++i)
    {
      const size_t prev = vertices[(i + vertices.size() - 1) % vertices.size()];
      const size_t current = vertices[i];
      const size_t next = vertices[(i + 1) % vertices.size()];
      if (cross_at(prev, current, next) <= eps)
      {
        continue;
      }

      bool contains_other_vertex = false;
      for (size_t candidate : vertices)
      {
        if (triangulationPlaneXYEqual(mesh, candidate, prev) || triangulationPlaneXYEqual(mesh, candidate, current)
          || triangulationPlaneXYEqual(mesh, candidate, next))
        {
          continue;
        }
        if (point_in_triangle(candidate, prev, current, next))
        {
          contains_other_vertex = true;
          break;
        }
      }
      if (contains_other_vertex)
      {
        continue;
      }

      if (orient_upwards)
      {
        addMeshletTriangle(mesh, prev, current, next, metadata, material_id);
      }
      else
      {
        addMeshletTriangle(mesh, prev, next, current, metadata, material_id);
      }
      vertices.erase(vertices.begin() + static_cast<std::ptrdiff_t>(i));
      clipped_ear = true;
      break;
    }

    if (!clipped_ear)
    {
      writeTriangulateSimplePolygonDebugTxt(kin_del, mesh,
        makeTriangulateSimplePolygonDebugPath(kin_del, mesh, "FAIL", ".txt"), "FAIL", {{"fail", vertices}});
      writeTriangulateSimplePolygonFailSvg(kin_del, mesh, vertices);
      throw std::runtime_error("triangulateSimplePolygon: failed to find an ear; polygon may be non-simple.");
    }
  }

  if (orient_upwards)
  {
    addMeshletTriangle(mesh, vertices[0], vertices[1], vertices[2], metadata, material_id);
  }
  else
  {
    addMeshletTriangle(mesh, vertices[0], vertices[2], vertices[1], metadata, material_id);
  }
}

void kinDS::SegmentBuilder::closingMeshTriangulatePolygons(
  VoronoiMesh& mesh, const std::vector<std::vector<size_t>>& polygons, double t, size_t strand_id)
{
  const std::string face_metadata = composeClosingMeshFaceMetadata(t, strand_id);
  for (const auto& polygon : polygons)
  {
    triangulateSimplePolygon(mesh, polygon, face_metadata);
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
  auto track_vertex = [&](const glm::dvec3& pos, std::optional<size_t> vv, const std::string& metadata) -> int
  {
    const int id = closingMeshAppendVertex(mesh, boundary_polygon, centroid, strand_id, t, pos, false, vv, metadata);
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
      strandInitDiagnosticLogLine("create_closing_mesh_empty", strand_id, t,
        "closing cap has no vertices/triangles after trace+triangulation");
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

  std::vector<bool> he_visited(kin_del.getGraph().halfEdgeSlotCount(), false);

  for (size_t i = 0; i < new_components.size(); i++)
  {
    size_t cid = component_ids[i];

    for (size_t v : new_components[i])
    {
      kin_del.component_data.component_map[v] = cid;
    }

    kin_del.component_data.components[cid] = new_components[i];
    kin_del.component_data.component_boundaries[cid]
      = kin_del.extractComponentBoundaries(new_components[i], t, he_visited, false, false);
    kin_del.component_data.component_centroids[cid]
      = polygonCentroid(kin_del.component_data.component_boundaries[cid][0]);
    kin_del.component_data.component_last_updated[cid] = t;
  }

  kin_del.notePendingBranchSplit(
    component_id, t, pre_split_parent_strands, new_components, component_ids);
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

  // Initialize the strand geometries at t = 0.0
  double t = 0.0; // TODO: might be customized later

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
          << " component_size=" << component.size() << " boundary_polygon_size=" << boundary_polygon.size();
      strandInitDiagnosticLogLine("init_cap_begin", strand_id, t, oss.str().c_str());
    }

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

  // now go through all half-edges and create a segment mesh pair
  for (size_t i = 0; i < half_edge_count; i += 2)
  {
    startNewMesh(i, t);
  }

  // Create boundary interval meshes at t=0 from precomputed crossing data.
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

  // Finalize closing meshes for strands still live in the graph (branches that ended earlier already got caps).
  for (size_t strand_id = 0; strand_id < graph.getVertexCount(); ++strand_id)
  {
    createClosingCapForStrand(strand_id, t);
  }

  accumulateSegmentProperties();

  resolveAllIntersectionFlexibleVertices("finalize intersection mesh");

  // compute normals
  for (auto& meshlet : meshes)
  {
    meshlet.ensureFaceMetadataSize();
    meshlet.computeNormals(NormalMode::PerTriangleCorner);
  }
  for (auto& meshlet : intersection_meshes)
  {
    meshlet.ensureFaceMetadataSize();
    meshlet.computeNormals(NormalMode::PerTriangleCorner);
  }

  auto remap1 = boundary_mesh.mergeDuplicateVertices();
  boundary_mesh.removeDegenerateTriangles();
  auto remap2 = boundary_mesh.removeIsolatedVertices();
  boundary_mesh.ensureFaceMetadataSize();
  boundary_mesh.computeNormals(NormalMode::PerTriangleCorner);

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
  for (size_t segment_id = 0; segment_id < segment_properties.size(); ++segment_id)
  {
    VoronoiMesh segment_mesh(MeshletExportMaterialNames);
    segment_mesh.setStoreMetadata(store_mesh_metadata);
    bool segment_mesh_initialized = false;
    std::vector<int> neighbor_segments_for_meshlet;
    const auto& properties = segment_properties[segment_id];
    double earliest_creation = std::numeric_limits<double>::quiet_NaN();
    auto append_oriented_mesh = [&](VoronoiMesh mesh, size_t seg0, size_t seg1)
    {
      const double mesh_ct = mesh.getCreationKineticTime();
      if (std::isfinite(mesh_ct))
      {
        if (!std::isfinite(earliest_creation) || mesh_ct < earliest_creation)
        {
          earliest_creation = mesh_ct;
        }
      }
      if (seg0 != segment_id)
      {
        mesh.flipOrientation();
      }
      const int neighbor = (seg0 == segment_id) ? static_cast<int>(seg1) : static_cast<int>(seg0);
      neighbor_segments_for_meshlet.insert(neighbor_segments_for_meshlet.end(), mesh.getTriangleCount(), neighbor);
      if (!segment_mesh_initialized)
      {
        segment_mesh = VoronoiMesh(MeshletExportMaterialNames, mesh.getNormalMode());
        segment_mesh.setStoreMetadata(store_mesh_metadata);
        segment_mesh_initialized = true;
      }
      segment_mesh += mesh;
    };

    // Regular meshlets referenced through segment_properties.
    for (size_t neighbor_index = 0; neighbor_index < properties.neighbor_count; ++neighbor_index)
    {
      size_t mesh_pair_index = properties.mesh_pair_indices[neighbor_index];
      VoronoiMesh mesh = meshes[mesh_pair_index];
      append_oriented_mesh(
        mesh, segment_mesh_pairs[mesh_pair_index].segment_index0, segment_mesh_pairs[mesh_pair_index].segment_index1);
    }

    // Boundary-interval meshlets are tracked in a separate pair/index space and must be merged explicitly.
    // Unlike regular Voronoi-edge strips, each boundary interval belongs to exactly one segment (its Voronoi cell),
    // and should not be orientation-flipped during merged export.
    for (size_t pair_idx = 0; pair_idx < intersection_segment_mesh_pairs.size(); ++pair_idx)
    {
      if (pair_idx >= intersection_meshes.size() || pair_idx >= intersection_mesh_pair_metadata.size())
      {
        continue;
      }

      const size_t owner_segment_id = intersection_mesh_pair_metadata[pair_idx].owner_segment_id;
      if (owner_segment_id == static_cast<size_t>(-1))
      {
        continue;
      }
      if (owner_segment_id != segment_id)
      {
        continue;
      }

      VoronoiMesh mesh = intersection_meshes[pair_idx];
      const double mesh_ct = mesh.getCreationKineticTime();
      if (std::isfinite(mesh_ct))
      {
        if (!std::isfinite(earliest_creation) || mesh_ct < earliest_creation)
        {
          earliest_creation = mesh_ct;
        }
      }
      neighbor_segments_for_meshlet.insert(neighbor_segments_for_meshlet.end(), mesh.getTriangleCount(), -2);
      if (!segment_mesh_initialized)
      {
        segment_mesh = VoronoiMesh(MeshletExportMaterialNames, mesh.getNormalMode());
        segment_mesh.setStoreMetadata(store_mesh_metadata);
        segment_mesh_initialized = true;
      }
      segment_mesh += mesh;
    }
    if (std::isfinite(earliest_creation))
    {
      segment_mesh.setCreationKineticTime(earliest_creation);
    }
    neighbor_segments.push_back(neighbor_segments_for_meshlet);
    segment_mesh.mergeDuplicateVertices(1e-4);
    segment_mesh.ensureFaceMetadataSize();
    meshlets.push_back(segment_mesh);

    if (diagnostics && segment_id == 0)
    {
      std::ostringstream oss;
      oss << "segment_id=0 verts=" << segment_mesh.getVertexCount() << " tris=" << segment_mesh.getTriangleCount()
          << " neighbor_count=" << properties.neighbor_count << " segment_mesh_initialized="
          << (segment_mesh_initialized ? "true" : "false");
      strandInitDiagnosticLogLine("extract_segment_meshlet", 0, std::numeric_limits<double>::quiet_NaN(),
        oss.str().c_str());
      if (!segment_mesh_initialized || segment_mesh.getVertexCount() == 0)
      {
        strandInitDiagnosticLogLine("extract_segment_meshlet_empty", 0, std::numeric_limits<double>::quiet_NaN(),
          "merged segment 0 meshlet is empty — check cap/strip wiring in accumulateSegmentProperties");
      }
    }
  }

  return std::make_pair(meshlets, neighbor_segments);
}

std::vector<std::string> kinDS::SegmentBuilder::extractSegmentMeshletExportSuffixes(bool merge_by_segment) const
{
  if (!merge_by_segment)
  {
    return meshlet_export_suffixes;
  }

  std::vector<std::string> suffixes;
  suffixes.reserve(segment_properties.size());
  for (size_t segment_id = 0; segment_id < segment_properties.size(); ++segment_id)
  {
    suffixes.push_back(std::string("_segment") + std::to_string(segment_id));
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
    throw std::runtime_error("strandIdForRawMeshlet: meshlet " + std::to_string(meshlet_index)
      + " has no segment endpoint.");
  }
  return strandIdForSegment(segment_id);
}

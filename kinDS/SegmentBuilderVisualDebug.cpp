#include "SegmentBuilderVisualDebug.hpp"

#include "DebugExportFormatting.hpp"
#include "HalfEdgeDelaunayGraphToSVG.hpp"
#include "KineticDelaunay.hpp"
#include "KineticDelaunayCrossingEvent.hpp"
#include "KineticDelaunayEventPredicates.hpp"
#include "Logger.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <exception>
#include <filesystem>
#include <limits>
#include <sstream>
#include <string>
#include <thread>
#include <unordered_map>
#include <unordered_set>

namespace kinDS
{
namespace
{
std::unordered_map<size_t, size_t> buildSiteInputBranchLabels(
  KineticDelaunay& kin_del, const HalfEdgeDelaunayGraph& graph, double occurrence_time)
{
  std::unordered_map<size_t, size_t> labels;
  const size_t vertex_count = graph.getVertexCount();
  const size_t section_index = kin_del.inputBranchSectionIndexAtIntervalUpperBound(
    eventIntervalUpperBound(occurrence_time));

  for (size_t strand_id = 0; strand_id < vertex_count; ++strand_id)
  {
    if (kin_del.isDummyBoundary(strand_id) || !kin_del.strandInActiveRuntimeScope(strand_id)
      || !kin_del.strandHasProfileAt(strand_id, occurrence_time))
    {
      continue;
    }

    labels[strand_id] = kin_del.getStrandTree().getBranchIndex(strand_id, section_index);
  }

  return labels;
}

void noteRuntimeBranchForStrand(
  std::unordered_set<size_t>& runtime_branches, const KineticDelaunay& kin_del, size_t strand_id,
  bool separate_pending_splits)
{
  if (kin_del.isDummyBoundary(strand_id))
  {
    return;
  }

  const auto& runtime_branch_map = kin_del.getRuntimeBranchMap();
  if (strand_id >= runtime_branch_map.size())
  {
    return;
  }

  const size_t branch_id = runtime_branch_map[strand_id];
  runtime_branches.insert(separate_pending_splits ? branch_id : kin_del.unsplitRuntimeBranchId(branch_id));
}

std::optional<size_t> uniqueBranchOrNull(const std::unordered_set<size_t>& branches)
{
  if (branches.size() == 1)
  {
    return *branches.begin();
  }
  return std::nullopt;
}

std::string formatBranchIdUnion(const std::unordered_set<size_t>& ids, bool none_is_unassigned)
{
  if (ids.empty())
  {
    return "none";
  }
  std::vector<size_t> sorted(ids.begin(), ids.end());
  std::sort(sorted.begin(), sorted.end());
  auto token = [&](size_t id) -> std::string
  {
    if (none_is_unassigned && id == KineticDelaunay::RuntimeBranchData::no_branch)
    {
      return "none";
    }
    return std::to_string(id);
  };
  if (sorted.size() == 1)
  {
    return token(sorted.front());
  }
  std::ostringstream oss;
  oss << "{";
  for (size_t i = 0; i < sorted.size(); ++i)
  {
    if (i != 0)
    {
      oss << ",";
    }
    oss << token(sorted[i]);
  }
  oss << "}";
  return oss.str();
}

void collectOmittedStrandBranchUnions(const KineticDelaunay& kin_del, const std::unordered_set<size_t>& omitted_ids,
  double t, bool separate_pending_splits, std::unordered_set<size_t>& runtime_branches,
  std::unordered_set<size_t>& input_branches)
{
  runtime_branches.clear();
  input_branches.clear();
  const auto& runtime_branch_map = kin_del.getRuntimeBranchMap();
  const size_t section = static_cast<size_t>(std::ceil(std::max(0.0, t)));
  for (size_t strand_id : omitted_ids)
  {
    if (kin_del.isDummyBoundary(strand_id))
    {
      continue;
    }
    if (strand_id < runtime_branch_map.size())
    {
      const size_t branch_id = runtime_branch_map[strand_id];
      runtime_branches.insert(separate_pending_splits ? branch_id : kin_del.unsplitRuntimeBranchId(branch_id));
    }
    else
    {
      runtime_branches.insert(KineticDelaunay::RuntimeBranchData::no_branch);
    }
    try
    {
      // Last sample when the strand already ended; that is the input branch that should have been retired.
      input_branches.insert(kin_del.getStrandTree().getBranchIndex(strand_id, section));
    }
    catch (const std::exception&)
    {
    }
  }
}

void warnOmittedStrandsForSvg(const KineticDelaunay& kin_del, const std::unordered_set<size_t>& considered_strands,
  EventTime t, const char* phase, bool separate_pending_splits, const std::string& svg_path)
{
  std::unordered_set<size_t> omitted;
  size_t considered = 0;
  for (size_t strand_id : considered_strands)
  {
    if (kin_del.isDummyBoundary(strand_id))
    {
      continue;
    }
    ++considered;
    if (!kin_del.strandHasProfileAt(strand_id, t))
    {
      omitted.insert(strand_id);
    }
  }
  if (omitted.empty() || considered == 0)
  {
    return;
  }

  std::unordered_set<size_t> runtime_branches;
  std::unordered_set<size_t> input_branches;
  collectOmittedStrandBranchUnions(kin_del, omitted, t, separate_pending_splits, runtime_branches, input_branches);
  std::filesystem::path path(svg_path);
  if (path.is_relative())
  {
    std::error_code ec;
    const std::filesystem::path cwd = std::filesystem::current_path(ec);
    if (!ec)
    {
      path = cwd / path;
    }
  }
  KINDS_WARNING("visual debug SVG omitted " << omitted.size() << "/" << considered
    << " strand(s) without a profile sample at t=" << t << " phase=" << (phase != nullptr ? phase : "?")
    << " runtime_branch=" << formatBranchIdUnion(runtime_branches, true)
    << " input_branch=" << formatBranchIdUnion(input_branches, false) << " file=" << path.generic_string());
}

std::optional<size_t> inferEventRuntimeBranchFromHighlight(
  const HalfEdgeDelaunayGraph& graph, const KineticDelaunay& kin_del, const VisualDebugHighlight& highlight,
  bool separate_pending_splits)
{
  std::unordered_set<size_t> runtime_branches;

  for (size_t strand_id : highlight.delaunay_vertices)
  {
    noteRuntimeBranchForStrand(runtime_branches, kin_del, strand_id, separate_pending_splits);
  }

  for (size_t he_id : highlight.directed_half_edges)
  {
    const int origin = graph.halfEdge(he_id).origin;
    if (origin >= 0)
    {
      noteRuntimeBranchForStrand(runtime_branches, kin_del, static_cast<size_t>(origin), separate_pending_splits);
    }
    const int dest = graph.destination(he_id);
    if (dest >= 0)
    {
      noteRuntimeBranchForStrand(runtime_branches, kin_del, static_cast<size_t>(dest), separate_pending_splits);
    }
  }

  for (size_t voronoi_vertex_id : highlight.voronoi_vertices)
  {
    if (voronoi_vertex_id >= graph.faceSlotCount() || !graph.isLiveFace(voronoi_vertex_id))
    {
      continue;
    }
    const auto verts = graph.getTriangleVertexIndices(voronoi_vertex_id);
    for (int v : verts)
    {
      if (v >= 0)
      {
        noteRuntimeBranchForStrand(runtime_branches, kin_del, static_cast<size_t>(v), separate_pending_splits);
      }
    }
  }

  for (size_t voronoi_edge_id : highlight.voronoi_edges)
  {
    const size_t he_even = 2 * voronoi_edge_id;
    if (he_even + 1 >= graph.halfEdgeSlotCount())
    {
      continue;
    }
    const int origin_even = graph.halfEdge(he_even).origin;
    const int origin_odd = graph.halfEdge(he_even + 1).origin;
    if (origin_even >= 0)
    {
      noteRuntimeBranchForStrand(runtime_branches, kin_del, static_cast<size_t>(origin_even), separate_pending_splits);
    }
    if (origin_odd >= 0)
    {
      noteRuntimeBranchForStrand(runtime_branches, kin_del, static_cast<size_t>(origin_odd), separate_pending_splits);
    }
  }

  return uniqueBranchOrNull(runtime_branches);
}

/// Chronologically sortable phase token: "!before", then "!error" / "!skipped" / "!failsafe", then "after".
std::string chronologicalPhaseToken(const char* phase)
{
  const std::string token(phase);
  if (token == "before")
  {
    return "!before";
  }
  // Sorts after !before and before "after": handleEvent check runs between those callbacks.
  if (token == "error")
  {
    return "!error";
  }
  if (token == "skipped")
  {
    return "!skipped";
  }
  if (token == "failsafe")
  {
    return "!failsafe";
  }
  return token;
}

bool delaunayEdgeTouchesStrands(
  const HalfEdgeDelaunayGraph& graph, size_t delaunay_edge_id, const std::unordered_set<size_t>* strands)
{
  if (strands == nullptr)
  {
    return true;
  }
  const size_t he0 = 2 * delaunay_edge_id;
  if (he0 + 1 >= graph.halfEdgeSlotCount())
  {
    return false;
  }
  const int origin_a = graph.halfEdge(he0).origin;
  const int origin_b = graph.halfEdge(he0 ^ 1).origin;
  if (origin_a < 0 || origin_b < 0)
  {
    const int finite = origin_a >= 0 ? origin_a : origin_b;
    return finite >= 0 && strands->count(static_cast<size_t>(finite)) != 0;
  }
  return strands->count(static_cast<size_t>(origin_a)) != 0 && strands->count(static_cast<size_t>(origin_b)) != 0;
}

std::string intersectionErrorFilenameSuffix(const VisualDebugHighlight& highlight, const HalfEdgeDelaunayGraph& graph,
  const std::unordered_set<size_t>* branch_strands)
{
  if (highlight.crossing_error_first_keys.empty() && highlight.crossing_error_repeat_keys.empty())
  {
    return {};
  }

  auto key_on_branch = [&](uint64_t packed) -> bool
  {
    return delaunayEdgeTouchesStrands(graph, static_cast<size_t>(packed >> 32), branch_strands);
  };

  bool has_first = false;
  bool has_repeat = false;
  for (uint64_t packed : highlight.crossing_error_first_keys)
  {
    if (key_on_branch(packed))
    {
      has_first = true;
      break;
    }
  }
  for (uint64_t packed : highlight.crossing_error_repeat_keys)
  {
    if (key_on_branch(packed))
    {
      has_repeat = true;
      break;
    }
  }

  std::string suffix;
  if (has_first)
  {
    suffix += "_error_first";
  }
  if (has_repeat)
  {
    suffix += "_error_recurring";
  }
  return suffix;
}

std::string currentDebugThreadIdToken()
{
  std::ostringstream os;
  os << std::this_thread::get_id();
  std::string token = os.str();
  for (char& c : token)
  {
    if (!std::isalnum(static_cast<unsigned char>(c)))
    {
      c = '_';
    }
  }
  return token;
}

std::string visualDebugSvgRelativePath(EventTime occurrence_time, const char* phase, const std::string& event_descriptor,
  std::optional<size_t> runtime_branch_id, const std::optional<std::filesystem::path>& output_root,
  std::optional<EventTime> creation_time, std::optional<uint64_t> event_id,
  const std::string& error_filename_suffix = {})
{
  std::string basename = formatDebugExportTimeToken(occurrence_time);
  if (event_id.has_value())
  {
    basename += "_eid" + std::to_string(*event_id);
  }
  basename += "_tid" + currentDebugThreadIdToken();
  basename += "_segmentbuilder_" + chronologicalPhaseToken(phase) + "_" + event_descriptor;
  if (creation_time.has_value())
  {
    basename += "_" + formatDebugExportTimeToken(*creation_time);
  }
  basename += error_filename_suffix;
  basename += ".svg";
  const std::string branch_folder = runtime_branch_id.has_value()
    ? ("branch" + std::to_string(runtime_branch_id.value()))
    : kVisualDebugUnresolvedBranchFolder;
  const std::string relative = branch_folder + "/" + basename;
  if (!output_root.has_value())
  {
    return relative;
  }
  return (*output_root / relative).generic_string();
}

/// Live graph vertices of the active scheduler; other forked schedulers' sites may not be valid at its time.
std::unordered_set<size_t> collectLiveStrandIds(const HalfEdgeDelaunayGraph& graph, const KineticDelaunay& kin_del)
{
  std::unordered_set<size_t> strand_ids;
  auto add_vertex = [&](int vertex)
  {
    if (vertex >= 0 && kin_del.strandInActiveRuntimeScope(static_cast<size_t>(vertex)))
    {
      strand_ids.insert(static_cast<size_t>(vertex));
    }
  };

  for (size_t he_id : graph.liveDelaunayEdges())
  {
    add_vertex(graph.halfEdge(he_id).origin);
    add_vertex(graph.destination(he_id));
  }

  return strand_ids;
}

std::unordered_set<size_t> collectActiveRuntimeBranches(
  const KineticDelaunay& kin_del, const std::unordered_set<size_t>& live_strand_ids, bool separate_pending_splits)
{
  std::unordered_set<size_t> runtime_branches;
  const auto& runtime_branch_map = kin_del.getRuntimeBranchMap();
  for (size_t strand_id : live_strand_ids)
  {
    if (strand_id < runtime_branch_map.size())
    {
      const size_t branch_id = runtime_branch_map[strand_id];
      runtime_branches.insert(
        separate_pending_splits ? branch_id : kin_del.unsplitRuntimeBranchId(branch_id));
    }
  }
  return runtime_branches;
}

std::unordered_set<size_t> collectStrandIdsForRuntimeBranch(
  const KineticDelaunay& kin_del, const std::unordered_set<size_t>& live_strand_ids, size_t runtime_branch_id,
  bool separate_pending_splits)
{
  std::unordered_set<size_t> strand_ids;
  if (separate_pending_splits)
  {
    // Exact runtime branch: retained parent or pending child strands only (already reassigned at note-pending).
    const auto& branches = kin_del.getRuntimeBranchData().branches;
    if (runtime_branch_id < branches.size())
    {
      for (size_t strand_id : branches[runtime_branch_id])
      {
        if (live_strand_ids.count(strand_id) != 0)
        {
          strand_ids.insert(strand_id);
        }
      }
    }
    return strand_ids;
  }

  // Default: unsplit branch id — gather its own strands plus those of any pending split-off children.
  for (size_t strand_id : kin_del.collectUnsplitRuntimeBranchStrands(runtime_branch_id))
  {
    if (live_strand_ids.count(strand_id) != 0)
    {
      strand_ids.insert(strand_id);
    }
  }
  return strand_ids;
}

std::unordered_set<size_t> collectStrandsNeededForBranchVoronoiGeometry(
  const HalfEdgeDelaunayGraph& graph, const KineticDelaunay& kin_del,
  const std::unordered_set<size_t>& branch_strands)
{
  // Start with the folder's own strands. For any finite Delaunay edge that will be drawn in this folder
  // (both primal origins in-set), also include finite sites of both adjacent faces. That covers Voronoi
  // edges that still span a pending sibling branch before the graph cut: the opposite-side site must be
  // positioned (in the folder frame) so circumcenters are not pulled toward (0,0) / a foreign frame.
  std::unordered_set<size_t> needed = branch_strands;
  for (size_t he_id : graph.liveDelaunayEdges())
  {
    const HalfEdgeDelaunayGraph::HalfEdge& he = graph.halfEdge(he_id);
    const HalfEdgeDelaunayGraph::HalfEdge& twin = graph.halfEdge(he_id ^ 1);
    if (he.origin < 0 || twin.origin < 0)
    {
      continue;
    }
    const size_t origin = static_cast<size_t>(he.origin);
    const size_t twin_origin = static_cast<size_t>(twin.origin);
    if (branch_strands.count(origin) == 0 || branch_strands.count(twin_origin) == 0)
    {
      continue;
    }

    for (size_t face_he_id : { he_id, he_id ^ 1 })
    {
      const int face_id = graph.halfEdge(face_he_id).face;
      if (face_id < 0)
      {
        continue;
      }
      for (int vertex : graph.getTriangleVertexIndices(static_cast<size_t>(face_id)))
      {
        if (vertex < 0 || kin_del.isDummyBoundary(static_cast<size_t>(vertex)))
        {
          continue;
        }
        needed.insert(static_cast<size_t>(vertex));
      }
    }
  }
  return needed;
}

std::pair<std::vector<glm::dvec2>, std::unordered_set<size_t>> buildVisualDebugStrandPositions(
  KineticDelaunay& kin_del, const HalfEdgeDelaunayGraph& graph, double occurrence_time,
  const std::unordered_set<size_t>& strands_to_position,
  const std::unordered_set<size_t>& frame_strand_ids, std::unordered_set<size_t>& omitted_strands)
{
  std::vector<glm::dvec2> points(graph.getVertexCount(), glm::dvec2(0.0));
  std::unordered_set<size_t> positioned_strands;
  positioned_strands.reserve(strands_to_position.size());

  // Prefer the folder's own strands for the shared input-branch frame so opposite-side sites are remapped
  // into this SVG's reference frame (not a min-id that might belong to the other pending branch).
  const std::unordered_set<size_t>& frame_source
    = frame_strand_ids.empty() ? strands_to_position : frame_strand_ids;
  std::vector<size_t> frame_strands;
  frame_strands.reserve(frame_source.size());
  for (size_t strand_id : frame_source)
  {
    if (kin_del.isDummyBoundary(strand_id))
    {
      continue;
    }
    if (!kin_del.strandHasProfileAt(strand_id, occurrence_time))
    {
      continue;
    }
    frame_strands.push_back(strand_id);
  }

  const size_t shared_reference_branch = frame_strands.empty()
    ? 0
    : kin_del.getSharedReferenceBranchForStrands(frame_strands, occurrence_time);

  for (size_t strand_id : strands_to_position)
  {
    if (kin_del.isDummyBoundary(strand_id))
    {
      continue;
    }
    if (!kin_del.strandHasProfileAt(strand_id, occurrence_time))
    {
      omitted_strands.insert(strand_id);
      continue;
    }
    try
    {
      points[strand_id] = kin_del.getPointInDelaunaySpace(strand_id, occurrence_time, shared_reference_branch);
      positioned_strands.insert(strand_id);
    }
    catch (const std::exception&)
    {
      omitted_strands.insert(strand_id);
    }
  }

  return { std::move(points), std::move(positioned_strands) };
}

std::unordered_map<size_t, glm::dvec3> buildSiteWorldPositions(
  KineticDelaunay& kin_del, double occurrence_time, const std::unordered_set<size_t>& strand_ids,
  std::unordered_set<size_t>& omitted_strands)
{
  std::unordered_map<size_t, glm::dvec3> world_positions;
  world_positions.reserve(strand_ids.size());

  for (size_t strand_id : strand_ids)
  {
    if (kin_del.isDummyBoundary(strand_id))
    {
      continue;
    }
    const size_t representative_strand_id = kin_del.representativeStrandIdForRuntimeBranch(strand_id);
    if (!kin_del.strandHasProfileAt(strand_id, occurrence_time)
      || !kin_del.strandHasProfileAt(representative_strand_id, occurrence_time))
    {
      omitted_strands.insert(strand_id);
      continue;
    }
    try
    {
      const glm::dvec2 profile_pos = kin_del.getPointInDelaunaySpace(strand_id, occurrence_time);
      glm::dvec3 profile_pos_3d(profile_pos, occurrence_time);
      world_positions[strand_id]
        = kin_del.getStrandTree().transformToObjectSpace(profile_pos_3d, representative_strand_id, occurrence_time);
    }
    catch (const std::exception&)
    {
      omitted_strands.insert(strand_id);
    }
  }

  return world_positions;
}

std::unordered_map<size_t, glm::dvec3> buildVoronoiVertexWorldPositions(
  KineticDelaunay& kin_del, const HalfEdgeDelaunayGraph& graph, double occurrence_time,
  std::unordered_set<size_t>& omitted_strands)
{
  std::unordered_map<size_t, glm::dvec3> world_positions;

  for (size_t voronoi_vertex_id : graph.liveFaces())
  {
    const auto& face = graph.face(voronoi_vertex_id);
    size_t reference_strand_id = static_cast<size_t>(-1);
    for (size_t he_id : face.half_edges)
    {
      const int origin = graph.halfEdge(he_id).origin;
      if (origin >= 0 && !kin_del.isDummyBoundary(static_cast<size_t>(origin)))
      {
        reference_strand_id = static_cast<size_t>(origin);
        break;
      }
    }
    if (reference_strand_id == static_cast<size_t>(-1) || !kin_del.strandInActiveRuntimeScope(reference_strand_id))
    {
      continue;
    }
    const size_t representative_strand_id = kin_del.representativeStrandIdForRuntimeBranch(reference_strand_id);
    if (!kin_del.strandHasProfileAt(reference_strand_id, occurrence_time)
      || !kin_del.strandHasProfileAt(representative_strand_id, occurrence_time))
    {
      omitted_strands.insert(reference_strand_id);
      continue;
    }
    try
    {
      const glm::dvec3 profile_pos
        = kin_del.computeVoronoiVertexClampedInfinity(face.half_edges[0], occurrence_time);
      glm::dvec3 world_pos_input = profile_pos;
      world_positions[voronoi_vertex_id]
        = kin_del.getStrandTree().transformToObjectSpace(world_pos_input, representative_strand_id, occurrence_time);
    }
    catch (const std::exception&)
    {
      omitted_strands.insert(reference_strand_id);
    }
  }

  return world_positions;
}

void writeVisualDebugSvgFile(const std::string& relative_path, const std::vector<glm::dvec2>& points,
  const HalfEdgeDelaunayGraph& graph, KineticDelaunay& kin_del,
  const SlotTable<size_t>& containing_tri_ids,
  const std::vector<HalfEdgeDelaunayGraphToSVG::IntersectionDebugInfo>& intersection_debug_data,
  const VisualDebugHighlight& highlight, const std::unordered_set<size_t>* highlighted_strands,
  const std::unordered_set<size_t>& positioned_strands,
  const std::unordered_map<size_t, size_t>& site_input_branch_labels,
  const std::unordered_map<size_t, glm::dvec3>& site_world_positions,
  const std::unordered_map<size_t, glm::dvec3>& voronoi_vertex_world_positions,
  const std::vector<HalfEdgeDelaunayGraphToSVG::SeparationOffsetSegment>* separation_offset_segments = nullptr,
  const std::vector<std::vector<glm::dvec2>>* seam_outlines = nullptr,
  const HalfEdgeDelaunayGraphToSVG::DebugOverlay* overlay = nullptr)
{
  const std::filesystem::path filepath(relative_path);
  if (filepath.has_parent_path())
  {
    std::filesystem::create_directories(filepath.parent_path());
  }
  const auto& face_inside_table = kin_del.getFacesInside();
  std::vector<bool> face_inside_flags(face_inside_table.size());
  for (size_t face_index = 0; face_index < face_inside_table.size(); ++face_index)
  {
    face_inside_flags[face_index] = face_inside_table[face_index] != 0;
  }
  const std::vector<size_t> containing_tri_copy(containing_tri_ids.begin(), containing_tri_ids.end());
  HalfEdgeDelaunayGraphToSVG::write(points, graph, relative_path, 0.1, &face_inside_flags, true,
    &containing_tri_copy, &intersection_debug_data, &highlight, highlighted_strands, &site_input_branch_labels,
    &positioned_strands, &site_world_positions, &voronoi_vertex_world_positions, separation_offset_segments,
    seam_outlines, overlay);
}

void writeVisualDebugSvgPair(EventTime occurrence_time, const char* phase, const std::string& event_descriptor,
  std::optional<size_t> runtime_branch_id, const std::optional<std::filesystem::path>& output_root,
  std::optional<EventTime> creation_time, std::optional<uint64_t> event_id, const std::vector<glm::dvec2>& points,
  const HalfEdgeDelaunayGraph& graph, KineticDelaunay& kin_del, const SlotTable<size_t>& containing_tri_ids,
  const std::vector<HalfEdgeDelaunayGraphToSVG::IntersectionDebugInfo>& intersection_debug_data,
  const VisualDebugHighlight& highlight_with_errors, const std::unordered_set<size_t>* highlighted_strands,
  const std::unordered_set<size_t>& positioned_strands,
  const std::unordered_map<size_t, size_t>& site_input_branch_labels,
  const std::unordered_map<size_t, glm::dvec3>& site_world_positions,
  const std::unordered_map<size_t, glm::dvec3>& voronoi_vertex_world_positions,
  const std::vector<HalfEdgeDelaunayGraphToSVG::SeparationOffsetSegment>* separation_offset_segments,
  const std::vector<std::vector<glm::dvec2>>* seam_outlines, const HalfEdgeDelaunayGraphToSVG::DebugOverlay* overlay)
{
  VisualDebugHighlight unmarked = highlight_with_errors;
  unmarked.crossing_error_first_keys.clear();
  unmarked.crossing_error_repeat_keys.clear();
  const std::string error_suffix
    = intersectionErrorFilenameSuffix(highlight_with_errors, graph, highlighted_strands);

  if (error_suffix.empty())
  {
    const std::string unmarked_path = visualDebugSvgRelativePath(
      occurrence_time, phase, event_descriptor, runtime_branch_id, output_root, creation_time, event_id, {});
    writeVisualDebugSvgFile(unmarked_path, points, graph, kin_del, containing_tri_ids, intersection_debug_data, unmarked,
      highlighted_strands, positioned_strands, site_input_branch_labels, site_world_positions,
      voronoi_vertex_world_positions, separation_offset_segments, seam_outlines, overlay);
    return;
  }

  const std::string marked_path = visualDebugSvgRelativePath(
    occurrence_time, phase, event_descriptor, runtime_branch_id, output_root, creation_time, event_id, error_suffix);
  writeVisualDebugSvgFile(marked_path, points, graph, kin_del, containing_tri_ids, intersection_debug_data,
    highlight_with_errors, highlighted_strands, positioned_strands, site_input_branch_labels, site_world_positions,
    voronoi_vertex_world_positions, separation_offset_segments, seam_outlines, overlay);
}
} // namespace

void writeSegmentBuilderVisualDebugSvg(bool visual_debug, KineticDelaunay& kin_del, const HalfEdgeDelaunayGraph& graph,
  EventTime occurrence_time, const char* phase, const std::string& event_descriptor,
  const VisualDebugHighlight& highlight, std::optional<size_t> event_runtime_branch_id,
  const std::vector<HalfEdgeDelaunayGraphToSVG::SeparationOffsetSegment>* separation_offset_segments,
  const std::vector<std::vector<glm::dvec2>>* seam_outlines, const std::vector<size_t>* explicit_runtime_branch_ids,
  std::optional<EventTime> creation_time, bool fan_out_active_runtime_branches, std::optional<uint64_t> event_id,
  const HalfEdgeDelaunayGraphToSVG::DebugOverlay* overlay)
{
  if (!visual_debug)
  {
    return;
  }
  try
  {
  const bool error_dump = !highlight.crossing_error_first_keys.empty()
    || !highlight.crossing_error_repeat_keys.empty();
  const bool force_failsafe = phase != nullptr && std::string(phase) == "failsafe";
  // Time window applies only to full --debug-files runs. Error dumps (intersection checks, overlays) always write.
  if (!error_dump && !force_failsafe && overlay == nullptr && kin_del.isVisualDebugEnabled()
    && !kin_del.shouldExportVisualDebugAt(occurrence_time.real_time))
  {
    return;
  }

  // Site / VV positions go through getPointInDelaunaySpace → separationOffsetAt, which reads
  // current_infinitesimal_t_. Bind it to this export's EventTime so before/after SVGs match the
  // filename's _i… token even if the event handler has already cleared the handle-path value.
  const KineticDelaunay::ScopedCurrentInfinitesimalTime infinitesimal_shift(kin_del, occurrence_time);

  const bool separate_pending_splits = kin_del.visualDebugSeparatePendingSplits();

  // By default, debug SVGs are grouped by unsplit branch (pending children collapse into the parent folder).
  // With --svg-separate-pending-splits, keep pending child ids so folders match future branches as soon as noted.
  if (event_runtime_branch_id.has_value() && !separate_pending_splits)
  {
    event_runtime_branch_id = kin_del.unsplitRuntimeBranchId(event_runtime_branch_id.value());
  }

  const auto& containing_tri_ids = kin_del.getCrossingData().getContainingTriIds();
  const auto intersection_debug_data = kin_del.getCrossingIntersectionDebugData(occurrence_time);
  VisualDebugHighlight svg_highlight = highlight;
  const std::unordered_map<size_t, size_t> site_input_branch_labels
    = buildSiteInputBranchLabels(kin_del, graph, occurrence_time);
  const std::unordered_set<size_t> live_strand_ids = collectLiveStrandIds(graph, kin_del);
  std::unordered_set<size_t> omitted_strands;
  const std::unordered_map<size_t, glm::dvec3> live_site_world_positions
    = buildSiteWorldPositions(kin_del, occurrence_time, live_strand_ids, omitted_strands);
  const std::unordered_map<size_t, glm::dvec3> voronoi_vertex_world_positions
    = buildVoronoiVertexWorldPositions(kin_del, graph, occurrence_time, omitted_strands);
  const std::optional<std::filesystem::path>& output_root = kin_del.getVisualDebugOutputRoot();
  const std::unordered_set<size_t> active_runtime_branches
    = collectActiveRuntimeBranches(kin_del, live_strand_ids, separate_pending_splits);

  auto try_write_for_runtime_branch = [&](size_t runtime_branch_id) -> bool
  {
    const size_t folder_branch_id
      = separate_pending_splits ? runtime_branch_id : kin_del.unsplitRuntimeBranchId(runtime_branch_id);
    const std::unordered_set<size_t> branch_strands
      = collectStrandIdsForRuntimeBranch(kin_del, live_strand_ids, folder_branch_id, separate_pending_splits);
    const std::unordered_set<size_t> geometry_strands
      = collectStrandsNeededForBranchVoronoiGeometry(graph, kin_del, branch_strands);
    const std::string error_suffix = intersectionErrorFilenameSuffix(svg_highlight, graph, &branch_strands);
    const std::string svg_path = visualDebugSvgRelativePath(occurrence_time, phase, event_descriptor, folder_branch_id,
      output_root, creation_time, event_id, error_suffix);
    warnOmittedStrandsForSvg(
      kin_del, branch_strands, occurrence_time, phase, separate_pending_splits, svg_path);
    const std::unordered_map<size_t, glm::dvec3> branch_site_world_positions
      = buildSiteWorldPositions(kin_del, occurrence_time, branch_strands, omitted_strands);
    const auto [points, positioned_strands]
      = buildVisualDebugStrandPositions(kin_del, graph, occurrence_time, geometry_strands, branch_strands,
        omitted_strands);
    if (positioned_strands.empty())
    {
      return false;
    }
    writeVisualDebugSvgPair(occurrence_time, phase, event_descriptor, folder_branch_id, output_root, creation_time,
      event_id, points, graph, kin_del, containing_tri_ids, intersection_debug_data, svg_highlight, &branch_strands,
      positioned_strands, site_input_branch_labels, branch_site_world_positions, voronoi_vertex_world_positions,
      separation_offset_segments, seam_outlines, overlay);
    return true;
  };

  auto resolve_branch_candidates = [&](std::optional<size_t> preferred_branch_id) -> std::vector<size_t>
  {
    std::vector<size_t> candidates;
    auto add_unique = [&](std::optional<size_t> branch_id)
    {
      if (!branch_id.has_value())
      {
        return;
      }
      const size_t resolved_branch_id = separate_pending_splits
        ? branch_id.value()
        : kin_del.unsplitRuntimeBranchId(branch_id.value());
      if (std::find(candidates.begin(), candidates.end(), resolved_branch_id) == candidates.end())
      {
        candidates.push_back(resolved_branch_id);
      }
    };

    add_unique(inferEventRuntimeBranchFromHighlight(graph, kin_del, highlight, separate_pending_splits));
    // Preferred id is second: callers sometimes pass a half-edge-derived fallback of 0 when the edge is
    // not live; trying that first would write a clipped SVG for the wrong branch and hide highlights.
    add_unique(preferred_branch_id);
    if (active_runtime_branches.size() == 1)
    {
      add_unique(*active_runtime_branches.begin());
    }
    return candidates;
  };

  const auto write_unresolved_branch_fallback = [&](const std::string& reason)
  {
    const std::string error_suffix = intersectionErrorFilenameSuffix(svg_highlight, graph, nullptr);
    const std::string svg_path = visualDebugSvgRelativePath(occurrence_time, phase, event_descriptor, std::nullopt,
      output_root, creation_time, event_id, error_suffix);
    warnOmittedStrandsForSvg(
      kin_del, live_strand_ids, occurrence_time, phase, separate_pending_splits, svg_path);
    const auto [points, positioned_strands]
      = buildVisualDebugStrandPositions(kin_del, graph, occurrence_time, live_strand_ids, live_strand_ids,
        omitted_strands);
    if (positioned_strands.empty())
    {
      KINDS_WARNING("writeSegmentBuilderVisualDebugSvg: branch resolution failed at t=" << occurrence_time
        << " phase=" << phase << " event=" << event_descriptor << " (" << reason
        << "); no live strands to export under " << kVisualDebugUnresolvedBranchFolder << ".");
      return;
    }
    KINDS_WARNING("writeSegmentBuilderVisualDebugSvg: branch resolution failed at t=" << occurrence_time
      << " phase=" << phase << " event=" << event_descriptor << " (" << reason << "); exporting all live strands to "
      << kVisualDebugUnresolvedBranchFolder << ".");
    writeVisualDebugSvgPair(occurrence_time, phase, event_descriptor, std::nullopt, output_root, creation_time, event_id,
      points, graph, kin_del, containing_tri_ids, intersection_debug_data, svg_highlight, nullptr, positioned_strands,
      site_input_branch_labels, live_site_world_positions, voronoi_vertex_world_positions, separation_offset_segments,
      seam_outlines, overlay);
  };

  // Explicit branch list: write only those folders (used by separation: parent while pending; parent+child after cut,
  // or parent+child while pending when --svg-separate-pending-splits is set; also crossing-param FAIL from edge ends).
  if (explicit_runtime_branch_ids != nullptr && !explicit_runtime_branch_ids->empty())
  {
    for (size_t branch_id : *explicit_runtime_branch_ids)
    {
      try_write_for_runtime_branch(branch_id);
    }
    return;
  }

  for (size_t branch_id : resolve_branch_candidates(event_runtime_branch_id))
  {
    if (try_write_for_runtime_branch(branch_id))
    {
      return;
    }
  }

  // Section events: no single owning branch — write one SVG per active runtime branch folder.
  if (fan_out_active_runtime_branches && !active_runtime_branches.empty())
  {
    bool wrote_any = false;
    for (size_t branch_id : active_runtime_branches)
    {
      wrote_any = try_write_for_runtime_branch(branch_id) || wrote_any;
    }
    if (!wrote_any)
    {
      write_unresolved_branch_fallback("active runtime branch fan-out had no positioned strands");
    }
    return;
  }

  // Known folder, but that branch has no positioned strands (e.g. section afterEvent after retire/vertex kill).
  // Do not dump another live branch's geometry into that event's filename / a shared folder.
  if (event_runtime_branch_id.has_value())
  {
    KINDS_WARNING("writeSegmentBuilderVisualDebugSvg: no positioned strands for runtime branch "
      << event_runtime_branch_id.value() << " at t=" << occurrence_time << " phase=" << phase
      << " event=" << event_descriptor << "; skipping (not exporting another branch).");
    return;
  }

  write_unresolved_branch_fallback("no unique runtime branch resolved");
  }
  catch (const std::exception& ex)
  {
    KINDS_WARNING("visual debug SVG skipped t=" << occurrence_time << " phase="
      << (phase != nullptr ? phase : "?") << " event=" << event_descriptor << ": " << ex.what());
  }
  catch (...)
  {
    KINDS_WARNING("visual debug SVG skipped t=" << occurrence_time << " phase="
      << (phase != nullptr ? phase : "?") << " event=" << event_descriptor << ": unknown error");
  }
}

} // namespace kinDS

#include "SegmentBuilderCrossingCallback.hpp"

#include "SegmentBuilder.hpp"
#include "KineticDelaunayCrossingEvent.hpp"
#include "SegmentBuilderVisualDebug.hpp"

#include <algorithm>
#include <array>
#include <optional>
#include <sstream>
#include <stdexcept>

namespace kinDS
{
namespace
{
void logCrossingMeshExtend(
  SegmentBuilder& sb, size_t voronoi_he_id, double t, VoronoiMesh& mesh, const char* branch, size_t tris_before)
{
  if (!sb.diagnostics)
  {
    return;
  }
  std::ostringstream o;
  o << "extend_crossing_" << branch << " d_tris=" << (mesh.getTriangleCount() - tris_before);
  sb.meshletDiagnosticLogLine("extend_mesh", voronoi_he_id, t, o.str().c_str());
}
} // namespace
void SegmentBuilderCrossingCallback::beforeEvent(KineticDelaunay::Event& e)
{
  auto* crossing = dynamic_cast<KineticDelaunay::CrossingEvent*>(&e);
  if (!crossing)
  {
    return;
  }
  SegmentBuilder::ScopedMetadataCallbackPhase callback_phase(segment_builder_, "before");

  auto& graph = segment_builder_.kin_del.getGraph();
  const size_t old_tri = graph.halfEdge(crossing->half_edge_id).face;
  const size_t new_tri = graph.halfEdge(crossing->half_edge_id ^ 1).face;
  int branch_vertex = graph.halfEdge(crossing->half_edge_id).origin;
  if (branch_vertex < 0)
  {
    branch_vertex = graph.destination(crossing->half_edge_id);
  }
  const size_t runtime_branch_id
    = segment_builder_.kin_del.getRuntimeBranchIdForStrand(static_cast<size_t>(branch_vertex));
  writeSegmentBuilderVisualDebugSvg(segment_builder_.visual_debug, segment_builder_.kin_del, graph,
    crossing->occurrence_time, "before",
    "crossing_v" + std::to_string(crossing->voronoi_vertex_id) + "_" + std::to_string(old_tri) + "_to_"
      + std::to_string(new_tri),
    VisualDebugHighlight::forCrossing(graph, crossing->half_edge_id, crossing->voronoi_vertex_id), runtime_branch_id);

  // Snapshot crossed-edge boundary interval links before CrossingData mutates them.
  crossing_edge_snapshot_.clear();
  crossing_edge_snapshot_delaunay_edge_id_ = crossing->half_edge_id / 2;
  if (!segment_builder_.kin_del.isOnComponentBoundary(crossing->half_edge_id))
  {
    return;
  }

  const auto& crossing_data = segment_builder_.kin_del.getCrossingData();
  if (crossing_edge_snapshot_delaunay_edge_id_ >= crossing_data.delaunay_edge_intersections.size())
  {
    return;
  }
  const auto& d_refs = crossing_data.delaunay_edge_intersections[crossing_edge_snapshot_delaunay_edge_id_];
  crossing_edge_snapshot_.reserve(d_refs.size());
  for (const auto& ref : d_refs)
  {
    CrossingEdgeSnapshotEntry s;
    s.voronoi_edge_id = ref->voronoi_edge_id;
    s.prev_pair_idx = ref->prev_segment_mesh_pair_index;
    s.next_pair_idx = ref->next_segment_mesh_pair_index;
    crossing_edge_snapshot_.push_back(s);
  }
  segment_builder_.maybeLogDiagnosticsMonitoredDelaunayEdgeTrigger(crossing->occurrence_time, "crossing_before_snapshot",
    crossing_edge_snapshot_delaunay_edge_id_, std::nullopt);
}

void SegmentBuilderCrossingCallback::afterEvent(KineticDelaunay::Event& e)
{
  auto* crossing = dynamic_cast<KineticDelaunay::CrossingEvent*>(&e);
  if (!crossing)
  {
    return;
  }
  SegmentBuilder::ScopedMetadataCallbackPhase callback_phase(segment_builder_, "after");
  const HalfEdgeDelaunayGraph& graph = segment_builder_.kin_del.getGraph();

  // `KineticDelaunay::CrossingEvent::handleEvent` already ran `updateAfterCrossingEvent`, which erases/inserts
  // `edge_intersections` records. Re-resolve iterators for strips incident to this Voronoi vertex before mesh edits.
  segment_builder_.refreshStripCrossingRefsIncidentToVoronoiVertex(crossing->voronoi_vertex_id);

  // Defer SVG until after SegmentBuilder updates crossing `prev`/`next` mesh-pair links (boundary merge/split) and
  // strip meshing so debug output matches runtime linkage.
  const auto write_crossing_visual_debug_svg = [&](size_t runtime_branch_id) {
    const size_t post_old_tri = graph.halfEdge(crossing->half_edge_id).face;
    const size_t post_new_tri = graph.halfEdge(crossing->half_edge_id ^ 1).face;
    writeSegmentBuilderVisualDebugSvg(segment_builder_.visual_debug, segment_builder_.kin_del, graph,
      crossing->occurrence_time, "after",
      "crossing_v" + std::to_string(crossing->voronoi_vertex_id) + "_" + std::to_string(post_old_tri) + "_to_"
        + std::to_string(post_new_tri),
      VisualDebugHighlight::forCrossing(graph, crossing->half_edge_id, crossing->voronoi_vertex_id), runtime_branch_id);
  };

  int branch_vertex = graph.halfEdge(crossing->half_edge_id).origin;
  if (branch_vertex < 0)
  {
    branch_vertex = graph.destination(crossing->half_edge_id);
  }
  const size_t runtime_branch_id
    = segment_builder_.kin_del.getRuntimeBranchIdForStrand(static_cast<size_t>(branch_vertex));

  const size_t voronoi_vertex_id = crossing->voronoi_vertex_id;
  const glm::dvec3 voronoi_vertex_position
    = segment_builder_.kin_del.computeVoronoiVertexClampedInfinity(
      graph.face(voronoi_vertex_id).half_edges[0], crossing->occurrence_time, false, false);
  const std::array<size_t, 3> half_edges = graph.face(voronoi_vertex_id).half_edges;
  if (!segment_builder_.kin_del.isOnComponentBoundary(crossing->half_edge_id))
  {
    write_crossing_visual_debug_svg(runtime_branch_id);
    return;
  }

  // Boundary-interval crossing handling on the crossed Delaunay edge.
  const size_t crossed_d_edge = crossing->half_edge_id / 2;
  const KineticDelaunay::CrossingData& crossing_data_after = segment_builder_.kin_del.getCrossingData();
  if (crossed_d_edge < crossing_data_after.delaunay_edge_intersections.size() && !crossing_edge_snapshot_.empty())
  {
    const KineticDelaunay::CrossingData::EdgeIntersectionRefListSlots::RefList& d_refs_after
      = crossing_data_after.delaunay_edge_intersections[crossed_d_edge];
    std::vector<size_t> after_voronoi_edge_ids;
    after_voronoi_edge_ids.reserve(d_refs_after.size());
    for (KineticDelaunay::CrossingData::EdgeIntersectionRef r : d_refs_after)
    {
      after_voronoi_edge_ids.push_back(r->voronoi_edge_id);
    }

    std::vector<size_t> old_voronoi_edge_ids;
    old_voronoi_edge_ids.reserve(crossing_edge_snapshot_.size());
    for (const CrossingEdgeSnapshotEntry& s : crossing_edge_snapshot_)
    {
      old_voronoi_edge_ids.push_back(s.voronoi_edge_id);
    }

    struct RemovedRef
    {
      size_t old_index = static_cast<size_t>(-1);
      CrossingEdgeSnapshotEntry snapshot;
    };
    std::vector<RemovedRef> removed;
    std::vector<std::pair<size_t, KineticDelaunay::CrossingData::EdgeIntersectionRef>> inserted;
    for (size_t i = 0; i < crossing_edge_snapshot_.size(); ++i)
    {
      const CrossingEdgeSnapshotEntry& old_ref = crossing_edge_snapshot_[i];
      if (std::find(after_voronoi_edge_ids.begin(), after_voronoi_edge_ids.end(), old_ref.voronoi_edge_id)
        == after_voronoi_edge_ids.end())
      {
        removed.push_back(RemovedRef { i, old_ref });
      }
    }
    for (size_t i = 0; i < d_refs_after.size(); ++i)
    {
      KineticDelaunay::CrossingData::EdgeIntersectionRefListSlots::RefList::const_iterator ref = d_refs_after.begin();
      std::advance(ref, static_cast<long long>(i));
      KineticDelaunay::CrossingData::EdgeIntersectionRef r = *ref;
      if (std::find(old_voronoi_edge_ids.begin(), old_voronoi_edge_ids.end(), r->voronoi_edge_id) == old_voronoi_edge_ids.end())
      {
        inserted.emplace_back(i, r);
      }
    }

    const glm::dvec3 event_pos = voronoi_vertex_position;
    const size_t component_id = segment_builder_.kin_del.component_data.component_map[graph.destination(crossing->half_edge_id)];
    std::vector<BoundaryPoint>& boundary_polygon
      = segment_builder_.kin_del.component_data.component_boundaries[component_id][0];
    const glm::dvec2 centroid_local = segment_builder_.kin_del.component_data.component_centroids[component_id];
    const bool boundary_even_he_is_outside = segment_builder_.kin_del.isOnComponentBoundaryOutside(2 * crossed_d_edge);
    const int inside_boundary_he_id
      = boundary_even_he_is_outside ? static_cast<int>(2 * crossed_d_edge + 1) : static_cast<int>(2 * crossed_d_edge);
    const std::string crossing_remove_meta = segment_builder_.composeBoundaryMetadata(
      SegmentBuilder::BoundaryEventType::Crossing, SegmentBuilder::BoundarySegmentAction::SegmentRemoved);
    const std::string crossing_update_meta = segment_builder_.composeBoundaryMetadata(
      SegmentBuilder::BoundaryEventType::Crossing, SegmentBuilder::BoundarySegmentAction::SegmentRemapped);
    auto with_crossing_case = [](const std::string& base_meta, const char* case_tag) -> std::string
    {
      return SegmentBuilder::MetadataBuilder::fromObject(base_meta)
        .ensureString("source", "Voronoi vertex")
        .addString("crossing_case", case_tag)
        .build();
    };

    auto update_pair_endpoint = [&](size_t pair_idx, bool update_start_endpoint,
                                  std::optional<KineticDelaunay::CrossingData::EdgeIntersectionRef> new_crossing,
                                  bool clear_segment_after_update, const std::string& metadata) {
      if (pair_idx == static_cast<size_t>(-1) || pair_idx >= segment_builder_.intersection_meshes.size()
        || pair_idx >= segment_builder_.intersection_mesh_pair_last_left_and_right_vertex.size())
      {
        return;
      }
      VoronoiMesh& mesh = segment_builder_.intersection_meshes[pair_idx];
      std::list<SegmentBuilder::MeshingData>& segs
        = segment_builder_.intersection_mesh_pair_last_left_and_right_vertex[pair_idx];
      if (segs.empty())
      {
        return;
      }
      SegmentBuilder::MeshingData& seg = segs.front();
      if (seg.mesh_start_vertex_id < 0 || seg.mesh_end_vertex_id < 0)
      {
        return;
      }
      bool resolved_update_start_endpoint = update_start_endpoint;
      // Keep endpoint updates consistent with the canonical crossed-edge direction:
      // CrossingData list order follows even half-edge direction, so pairs on the even-origin side
      // correspond to "prev" (update end), while odd-origin side corresponds to "next" (update start).
      const size_t he_even = 2 * crossed_d_edge;
      const size_t he_odd = he_even + 1;
      if (pair_idx < segment_builder_.intersection_mesh_pair_metadata.size() && he_odd < graph.halfEdgeSlotCount())
      {
        const MeshStructure::IntersectionMeshPairMetadata& pair_meta
          = segment_builder_.intersection_mesh_pair_metadata[pair_idx];
        const int even_origin = graph.halfEdge(he_even).origin;
        const int odd_origin = graph.halfEdge(he_odd).origin;
        if (even_origin >= 0 && odd_origin >= 0)
        {
          const size_t even_origin_cell = static_cast<size_t>(even_origin);
          const size_t odd_origin_cell = static_cast<size_t>(odd_origin);
          if (pair_meta.voronoi_cell_id == even_origin_cell)
          {
            resolved_update_start_endpoint = false; // prev-side => end endpoint.
          }
          else if (pair_meta.voronoi_cell_id == odd_origin_cell)
          {
            resolved_update_start_endpoint = true; // next-side => start endpoint.
          }
        }
      }
      // Closing the merged-away strip: one boundary vertex, not an interval endpoint — same convention as init `uniform`.
      const bool use_uniform_pos = clear_segment_after_update && !new_crossing.has_value();
      const std::string vertex_meta = segment_builder_.store_mesh_metadata
        ? SegmentBuilder::MetadataBuilder::fromObject(metadata)
            .ensureString("source", "Voronoi vertex")
            .addString("pos", use_uniform_pos ? "uniform" : (resolved_update_start_endpoint ? "left" : "right"))
            .build()
        : std::string {};
      const glm::dvec3 vertex_color = use_uniform_pos ? glm::dvec3(1.0, 0.0, 1.0)
                                                      : (resolved_update_start_endpoint ? glm::dvec3(1.0, 0.0, 0.0) : glm::dvec3(0.0, 0.0, 1.0));
      const size_t eff_l = segment_builder_.intersectionStripEffectiveVertexIndex(seg, true);
      const size_t eff_r = segment_builder_.intersectionStripEffectiveVertexIndex(seg, false);
      const size_t new_vid = segment_builder_.addMeshletVertex(mesh, boundary_polygon, centroid_local, event_pos,
        graph.destination(crossing->half_edge_id), crossing->occurrence_time, false, std::nullopt, vertex_meta, vertex_color);
      segment_builder_.addBoundaryIntervalTriangleOriented(
        mesh, eff_l, eff_r, new_vid, inside_boundary_he_id, crossing->occurrence_time, metadata, pair_idx);

      if (use_uniform_pos)
      {
        segment_builder_.applyIntersectionStripUniformClosureVertex(mesh, seg, new_vid);
      }
      else
      {
        segment_builder_.applyIntersectionStripOneSidedFixedVertex(mesh, seg, resolved_update_start_endpoint, new_vid,
          inside_boundary_he_id, new_crossing, boundary_polygon, centroid_local, graph.destination(crossing->half_edge_id),
          crossing->occurrence_time, !clear_segment_after_update);
      }

      if (clear_segment_after_update)
      {
        segs.clear();
        return;
      }
    };

    if (removed.size() == 2 && inserted.size() == 1)
    {
      KINDS_DEBUG("Boundary crossing case: merge_2_to_1 on de=" << crossed_d_edge << " t=" << crossing->occurrence_time);
      const std::string crossing_remove_meta_case = with_crossing_case(crossing_remove_meta, "merge_2_to_1");
      const std::string crossing_update_meta_case = with_crossing_case(crossing_update_meta, "merge_2_to_1");
      // Case 1 (merge): two removed, one inserted. Middle strip is whichever `next` of one equals the other's `prev`;
      // `removed` order does not matter.
      KineticDelaunay::CrossingData::EdgeIntersectionRef inserted_ref = inserted.front().second;
      const CrossingEdgeSnapshotEntry& s0 = removed[0].snapshot;
      const CrossingEdgeSnapshotEntry& s1 = removed[1].snapshot;

      // The completed strip is the mesh pair referenced from *both* removed crossings on the side between them.
      // The merged crossing must inherit only the outer neighbour links.
      size_t middle_pair = static_cast<size_t>(-1);
      size_t outer_prev_mesh = static_cast<size_t>(-1);
      size_t outer_next_mesh = static_cast<size_t>(-1);
      size_t update_end_of_strip_pair = static_cast<size_t>(-1);
      size_t update_start_of_strip_pair = static_cast<size_t>(-1);
      if (s0.next_pair_idx != static_cast<size_t>(-1) && s0.next_pair_idx == s1.prev_pair_idx)
      {
        middle_pair = s0.next_pair_idx;
        outer_prev_mesh = s0.prev_pair_idx;
        outer_next_mesh = s1.next_pair_idx;
        update_end_of_strip_pair = s0.prev_pair_idx;
        update_start_of_strip_pair = s1.next_pair_idx;
      }
      else if (s1.next_pair_idx != static_cast<size_t>(-1) && s1.next_pair_idx == s0.prev_pair_idx)
      {
        middle_pair = s1.next_pair_idx;
        outer_prev_mesh = s1.prev_pair_idx;
        outer_next_mesh = s0.next_pair_idx;
        update_end_of_strip_pair = s1.prev_pair_idx;
        update_start_of_strip_pair = s0.next_pair_idx;
      }
      else
      {
        std::ostringstream oss;
        oss << "Boundary crossing merge_2_to_1: could not match middle strip (expected removed[0].next==removed[1].prev "
               "or removed[1].next==removed[0].prev) on de="
            << crossed_d_edge << " t=" << crossing->occurrence_time << " s0.prev=" << s0.prev_pair_idx << " s0.next="
            << s0.next_pair_idx << " s1.prev=" << s1.prev_pair_idx << " s1.next=" << s1.next_pair_idx << ".";
        throw std::runtime_error(oss.str());
      }

      // The interval between the two removed intersections ends here.
      update_pair_endpoint(middle_pair, false, std::nullopt, true, crossing_remove_meta_case);

      // Adjacent intervals get a new vertex and become delimited by the newly inserted intersection.
      update_pair_endpoint(update_end_of_strip_pair, false, inserted_ref, false, crossing_update_meta_case);
      update_pair_endpoint(update_start_of_strip_pair, true, inserted_ref, false, crossing_update_meta_case);
      segment_builder_.assignIntersectionMeshPairLink(inserted_ref, true, outer_prev_mesh,
        "crossing_merge_2_to_1:inserted_prev", crossing->occurrence_time);
      segment_builder_.assignIntersectionMeshPairLink(inserted_ref, false, outer_next_mesh,
        "crossing_merge_2_to_1:inserted_next", crossing->occurrence_time);
    }
    else if (removed.size() == 1 && inserted.size() == 2)
    {
      KINDS_DEBUG("Boundary crossing case: split_1_to_2 on de=" << crossed_d_edge << " t=" << crossing->occurrence_time);
      const std::string crossing_update_meta_case = with_crossing_case(crossing_update_meta, "split_1_to_2");
      // Case 2 (split): one removed, two inserted — canonical (start,end) is adjacency in delaunay_edge_intersections.
      const KineticDelaunay::CrossingData::EdgeIntersectionRefListSlots::RefList& d_list_split
        = crossing_data_after.delaunay_edge_intersections[crossed_d_edge];
      const KineticDelaunay::CrossingData::EdgeIntersectionRef r0 = inserted[0].second;
      const KineticDelaunay::CrossingData::EdgeIntersectionRef r1 = inserted[1].second;
      const KineticDelaunay::CrossingData::EdgeIntersectionRefListSlots::RefList::const_iterator it0
        = std::find(d_list_split.begin(), d_list_split.end(), r0);
      const KineticDelaunay::CrossingData::EdgeIntersectionRefListSlots::RefList::const_iterator it1
        = std::find(d_list_split.begin(), d_list_split.end(), r1);
      if (it0 == d_list_split.end() || it1 == d_list_split.end())
      {
        std::ostringstream oss;
        oss << "Boundary crossing split_1_to_2: inserted crossing not found on delaunay_edge_intersections[" << crossed_d_edge
            << "] at t=" << crossing->occurrence_time << ".";
        throw std::runtime_error(oss.str());
      }
      KineticDelaunay::CrossingData::EdgeIntersectionRef start_ref;
      KineticDelaunay::CrossingData::EdgeIntersectionRef end_ref;
      if (std::next(it0) == it1)
      {
        start_ref = r0;
        end_ref = r1;
      }
      else if (std::next(it1) == it0)
      {
        start_ref = r1;
        end_ref = r0;
      }
      else
      {
        std::ostringstream oss;
        oss << "Boundary crossing split_1_to_2: inserted crossings are not adjacent in delaunay_edge_intersections["
            << crossed_d_edge << "] at t=" << crossing->occurrence_time << ".";
        throw std::runtime_error(oss.str());
      }

      const CrossingEdgeSnapshotEntry& old = removed[0].snapshot;

      // Outer topology: strips that met the old crossing keep the same mesh-pair ids on the open sides of the split.
      segment_builder_.assignIntersectionMeshPairLink(start_ref, true, old.prev_pair_idx,
        "crossing_split_1_to_2:start_prev", crossing->occurrence_time);
      segment_builder_.assignIntersectionMeshPairLink(end_ref, false, old.next_pair_idx,
        "crossing_split_1_to_2:end_next", crossing->occurrence_time);

      // Old adjacent intervals are advanced to event position and retargeted to the new delimiters.
      update_pair_endpoint(old.prev_pair_idx, false, start_ref, false, crossing_update_meta_case);
      update_pair_endpoint(old.next_pair_idx, true, end_ref, false, crossing_update_meta_case);

      // Middle strip between start and end; startNewMeshFromIntersections + writeIntersectionPairLinks set
      // start_ref->next_segment_mesh_pair_index and end_ref->prev_segment_mesh_pair_index to mid_pair.
      const size_t mid_cell
        = segment_builder_.determineVoronoiCellForBoundaryIntersectionInterval(crossed_d_edge, start_ref, end_ref);
      const size_t mid_pair = segment_builder_.startNewMeshFromIntersections(mid_cell, crossing->occurrence_time, start_ref,
        end_ref, false, SegmentBuilder::BoundaryEventType::Crossing, SegmentBuilder::BoundarySegmentAction::NewSegment,
        true);
      (void)mid_pair;
    }
    else
    {
      KINDS_WARNING("Boundary crossing case not handled yet on de=" << crossed_d_edge << " removed=" << removed.size()
                                                                    << " inserted=" << inserted.size() << " t="
                                                                    << crossing->occurrence_time);
    }
    segment_builder_.maybeLogDiagnosticsMonitoredDelaunayEdgeTrigger(
      crossing->occurrence_time, "crossing_after_boundary_interval", crossed_d_edge, std::nullopt);
  }

  size_t inside_he_id = crossing->half_edge_id;
  bool entering_boundary = false;
  if (segment_builder_.kin_del.isOnComponentBoundaryOutside(crossing->half_edge_id))
  {
    inside_he_id = crossing->half_edge_id ^ 1;
    entering_boundary = true;
  }

  size_t component_id = segment_builder_.kin_del.component_data.component_map[graph.destination(inside_he_id)];
  std::vector<BoundaryPoint>& boundary_polygon
    = segment_builder_.kin_del.component_data.component_boundaries[component_id][0];
  const glm::dvec2 centroid = segment_builder_.kin_del.component_data.component_centroids[component_id];

  for (size_t voronoi_he_id : half_edges)
  {
    const size_t strip_voronoi_edge_id = (voronoi_he_id & ~1) / 2;
    const size_t even_id = voronoi_he_id & ~static_cast<size_t>(1);
    const size_t odd_id = even_id + 1;
    const HalfEdgeDelaunayGraph::HalfEdge& he_even = graph.halfEdge(even_id);
    const HalfEdgeDelaunayGraph::HalfEdge& he_odd = graph.halfEdge(odd_id);
    const int strand_even_origin_i = static_cast<int>(he_even.origin);
    const int strand_odd_origin_i = static_cast<int>(he_odd.origin);

    const auto strip_vertex_meta = [&](const char* pos, const char* op,
                                      SegmentBuilder::BoundarySegmentAction segment_action) -> std::string
    {
      // These vertices are the event circumcenter (Voronoi vertex), not a Delaunay–Voronoi path intersection.
      // Boundary crossing refs are stored on MeshingData for strip topology only.
      return segment_builder_.composeRegularStripVertexMetadata(crossing->occurrence_time, strip_voronoi_edge_id, even_id,
        strand_even_origin_i, strand_odd_origin_i, SegmentBuilder::BoundaryEventType::Crossing, segment_action,
        std::nullopt, pos, op, "Voronoi vertex");
    };

    const auto strip_face_meta = [&](const char* op, SegmentBuilder::BoundarySegmentAction segment_action) -> std::string
    {
      return segment_builder_.composeRegularStripFaceMetadata(crossing->occurrence_time, strip_voronoi_edge_id, even_id,
        strand_even_origin_i, strand_odd_origin_i, SegmentBuilder::BoundaryEventType::Crossing, segment_action, op);
    };

    const auto strand_id_for_voronoi_strip = [&]() -> size_t
    {
      const int origin = graph.halfEdge(voronoi_he_id).origin;
      if (origin >= 0)
      {
        return static_cast<size_t>(origin);
      }
      if (inside_he_id >= 0)
      {
        const int inside_origin = graph.halfEdge(static_cast<size_t>(inside_he_id)).origin;
        if (inside_origin >= 0)
        {
          return static_cast<size_t>(inside_origin);
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
      throw std::runtime_error("SegmentBuilderCrossingCallback: no strand id for circumcenter vertex transform.");
    };
    const size_t strand_id_for_transform = strand_id_for_voronoi_strip();
    const std::optional<size_t> voronoi_vertex_for_alpha_check = voronoi_vertex_id;

    size_t& segment_mesh_pair_index = segment_builder_.half_edge_index_to_segment_mesh_pair_index[voronoi_he_id];
    VoronoiMesh& mesh = segment_builder_.meshes[segment_mesh_pair_index];

    std::list<SegmentBuilder::MeshingData>& segments
      = segment_builder_.segment_mesh_pair_last_left_and_right_vertex[segment_mesh_pair_index];
    std::list<SegmentBuilder::MeshingData>::iterator it = std::find_if(segments.begin(), segments.end(),
      [inside_he_id](const SegmentBuilder::MeshingData& data)
      {
        return data.end_half_edge_id == static_cast<int>(inside_he_id)
          || data.start_half_edge_id == static_cast<int>(inside_he_id);
      });
    if (it != segments.end())
    {
      if (!entering_boundary)
      {
        const std::string vertex_meta = strip_vertex_meta("circumcenter", "erase_strip",
          SegmentBuilder::BoundarySegmentAction::SegmentRemoved);
        const std::string face_meta = strip_face_meta("erase_strip", SegmentBuilder::BoundarySegmentAction::SegmentRemoved);
        size_t new_vertex_index = segment_builder_.addMeshletVertex(mesh, boundary_polygon, centroid,
          voronoi_vertex_position, strand_id_for_transform, crossing->occurrence_time, false, voronoi_vertex_for_alpha_check,
          vertex_meta);
        int last_left = it->mesh_start_vertex_id;
        int last_right = it->mesh_end_vertex_id;
        {
          const size_t tris_before = mesh.getTriangleCount();
          segment_builder_.addMeshletTriangle(mesh, last_left, last_right, new_vertex_index, face_meta);
          logCrossingMeshExtend(segment_builder_, voronoi_he_id, crossing->occurrence_time, mesh, "erase_strip", tris_before);
        }
        segments.erase(it);
      }
      else
      {
        if (it->start_half_edge_id == static_cast<int>(inside_he_id))
        {
          it->start_half_edge_id = -1;
          it->start_crossing.reset();
          const std::string vertex_meta = strip_vertex_meta("left", "enter_boundary_start",
            SegmentBuilder::BoundarySegmentAction::SegmentRemapped);
          const std::string face_meta
            = strip_face_meta("enter_boundary_start", SegmentBuilder::BoundarySegmentAction::SegmentRemapped);
          size_t new_vertex_index = segment_builder_.addMeshletVertex(mesh, boundary_polygon, centroid,
            voronoi_vertex_position, strand_id_for_transform, crossing->occurrence_time, false, voronoi_vertex_for_alpha_check,
            vertex_meta);
          int last_left = it->mesh_start_vertex_id;
          int last_right = it->mesh_end_vertex_id;
          {
            const size_t tris_before = mesh.getTriangleCount();
            segment_builder_.addMeshletTriangle(mesh, last_left, last_right, new_vertex_index, face_meta);
            logCrossingMeshExtend(
              segment_builder_, voronoi_he_id, crossing->occurrence_time, mesh, "enter_boundary_start", tris_before);
          }
          it->mesh_start_vertex_id = static_cast<int>(new_vertex_index);
        }
        else
        {
          assert(it->end_half_edge_id == static_cast<int>(inside_he_id));
          it->end_half_edge_id = -1;
          it->end_crossing.reset();
          const std::string vertex_meta = strip_vertex_meta("right", "enter_boundary_end",
            SegmentBuilder::BoundarySegmentAction::SegmentRemapped);
          const std::string face_meta
            = strip_face_meta("enter_boundary_end", SegmentBuilder::BoundarySegmentAction::SegmentRemapped);
          size_t new_vertex_index = segment_builder_.addMeshletVertex(mesh, boundary_polygon, centroid,
            voronoi_vertex_position, strand_id_for_transform, crossing->occurrence_time, false, voronoi_vertex_for_alpha_check,
            vertex_meta);
          int last_left = it->mesh_start_vertex_id;
          int last_right = it->mesh_end_vertex_id;
          {
            const size_t tris_before = mesh.getTriangleCount();
            segment_builder_.addMeshletTriangle(mesh, last_left, last_right, new_vertex_index, face_meta);
            logCrossingMeshExtend(
              segment_builder_, voronoi_he_id, crossing->occurrence_time, mesh, "enter_boundary_end", tris_before);
          }
          it->mesh_end_vertex_id = static_cast<int>(new_vertex_index);
        }
      }
    }
    else
    {
      const size_t end_voronoi_vertex_id = graph.halfEdge(odd_id).face;

      if (end_voronoi_vertex_id == voronoi_vertex_id)
      {
        if (entering_boundary)
        {
          const std::optional<KineticDelaunay::CrossingData::EdgeIntersectionRef> boundary_crossing
            = segment_builder_.closingMeshFindVoronoiEdgeIntersection(strip_voronoi_edge_id, inside_he_id);
          const std::string vertex_meta = strip_vertex_meta("left", "enter_boundary_tail_seed",
            SegmentBuilder::BoundarySegmentAction::SegmentRemapped);
          size_t new_vertex_index = segment_builder_.addMeshletVertex(mesh, boundary_polygon, centroid,
            voronoi_vertex_position, strand_id_for_transform, crossing->occurrence_time, false, voronoi_vertex_for_alpha_check,
            vertex_meta);
          SegmentBuilder::MeshingData seg { static_cast<int>(new_vertex_index), static_cast<int>(new_vertex_index),
            static_cast<int>(inside_he_id), -1 };
          seg.start_crossing = boundary_crossing;
          segments.emplace_back(std::move(seg));
        }
        else
        {
          assert(segments.back().end_half_edge_id == -1);
          const std::optional<KineticDelaunay::CrossingData::EdgeIntersectionRef> boundary_crossing
            = segment_builder_.closingMeshFindVoronoiEdgeIntersection(strip_voronoi_edge_id, inside_he_id);
          const std::string vertex_meta = strip_vertex_meta("right", "tail_close",
            SegmentBuilder::BoundarySegmentAction::SegmentRemapped);
          const std::string face_meta
            = strip_face_meta("tail_close", SegmentBuilder::BoundarySegmentAction::SegmentRemapped);
          size_t new_vertex_index = segment_builder_.addMeshletVertex(mesh, boundary_polygon, centroid,
            voronoi_vertex_position, strand_id_for_transform, crossing->occurrence_time, false, voronoi_vertex_for_alpha_check,
            vertex_meta);
          int last_left = segments.back().mesh_start_vertex_id;
          int last_right = segments.back().mesh_end_vertex_id;
          {
            const size_t tris_before = mesh.getTriangleCount();
            segment_builder_.addMeshletTriangle(mesh, last_left, last_right, new_vertex_index, face_meta);
            logCrossingMeshExtend(
              segment_builder_, voronoi_he_id, crossing->occurrence_time, mesh, "tail_close", tris_before);
          }
          segments.back().end_half_edge_id = static_cast<int>(inside_he_id);
          segments.back().end_crossing = boundary_crossing;
          segments.back().mesh_end_vertex_id = static_cast<int>(new_vertex_index);
        }
      }
      else
      {
        if (entering_boundary)
        {
          const std::optional<KineticDelaunay::CrossingData::EdgeIntersectionRef> boundary_crossing
            = segment_builder_.closingMeshFindVoronoiEdgeIntersection(strip_voronoi_edge_id, inside_he_id);
          const std::string vertex_meta = strip_vertex_meta("right", "enter_boundary_head_seed",
            SegmentBuilder::BoundarySegmentAction::SegmentRemapped);
          size_t new_vertex_index = segment_builder_.addMeshletVertex(mesh, boundary_polygon, centroid,
            voronoi_vertex_position, strand_id_for_transform, crossing->occurrence_time, false, voronoi_vertex_for_alpha_check,
            vertex_meta);
          SegmentBuilder::MeshingData seg { static_cast<int>(new_vertex_index), static_cast<int>(new_vertex_index), -1,
            static_cast<int>(inside_he_id) };
          seg.end_crossing = boundary_crossing;
          segments.emplace_front(std::move(seg));
        }
        else
        {
          assert(segments.front().start_half_edge_id == -1);
          const std::optional<KineticDelaunay::CrossingData::EdgeIntersectionRef> boundary_crossing
            = segment_builder_.closingMeshFindVoronoiEdgeIntersection(strip_voronoi_edge_id, inside_he_id);
          const std::string vertex_meta = strip_vertex_meta("left", "head_close",
            SegmentBuilder::BoundarySegmentAction::SegmentRemapped);
          const std::string face_meta
            = strip_face_meta("head_close", SegmentBuilder::BoundarySegmentAction::SegmentRemapped);
          size_t new_vertex_index = segment_builder_.addMeshletVertex(mesh, boundary_polygon, centroid,
            voronoi_vertex_position, strand_id_for_transform, crossing->occurrence_time, false, voronoi_vertex_for_alpha_check,
            vertex_meta);
          int last_left = segments.front().mesh_start_vertex_id;
          int last_right = segments.front().mesh_end_vertex_id;
          {
            const size_t tris_before = mesh.getTriangleCount();
            segment_builder_.addMeshletTriangle(mesh, last_left, last_right, new_vertex_index, face_meta);
            logCrossingMeshExtend(
              segment_builder_, voronoi_he_id, crossing->occurrence_time, mesh, "head_close", tris_before);
          }
          segments.front().start_half_edge_id = static_cast<int>(inside_he_id);
          segments.front().start_crossing = boundary_crossing;
          segments.front().mesh_start_vertex_id = static_cast<int>(new_vertex_index);
        }
      }
    }

    for (SegmentBuilder::MeshingData& seg : segments)
    {
      segment_builder_.refreshMeshingDataCrossingRefs(seg, strip_voronoi_edge_id);
    }
  }

  write_crossing_visual_debug_svg(runtime_branch_id);
}
} // namespace kinDS


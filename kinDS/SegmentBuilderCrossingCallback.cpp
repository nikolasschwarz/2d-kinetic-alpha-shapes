#include "SegmentBuilderCrossingCallback.hpp"

#include "SegmentBuilder.hpp"
#include "KineticDelaunayCrossingEvent.hpp"

#include <algorithm>
#include <optional>
#include <sstream>

namespace kinDS
{
namespace
{
void logCrossingMeshExtend(
  SegmentBuilder& sb, size_t voronoi_he_id, double t, VoronoiMesh& mesh, const char* branch, size_t tris_before)
{
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

  // Snapshot the crossed Delaunay edge intersection order + interval-mesh links before CrossingData mutates them.
  // We compare this against the post-event state in afterEvent to detect removed/inserted intersections.
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
}

void SegmentBuilderCrossingCallback::afterEvent(KineticDelaunay::Event& e)
{
  auto* crossing = dynamic_cast<KineticDelaunay::CrossingEvent*>(&e);
  if (!crossing)
  {
    return;
  }
  auto& graph = segment_builder_.kin_del.getGraph();

  // `KineticDelaunay::CrossingEvent::handleEvent` already ran `updateAfterCrossingEvent`, which erases/inserts
  // `edge_intersections` records. Re-resolve iterators for strips incident to this Voronoi vertex before mesh edits.
  segment_builder_.refreshStripCrossingRefsIncidentToVoronoiVertex(crossing->voronoi_vertex_id);

  if (segment_builder_.visual_debug)
  {
    std::vector<glm::dvec2> points = segment_builder_.kin_del.getPointsAt(crossing->occurrence_time);
    size_t old_tri = graph.getHalfEdges()[crossing->half_edge_id].face;
    size_t new_tri = graph.getHalfEdges()[crossing->half_edge_id ^ 1].face;
    std::string filename = "t" + std::to_string(crossing->occurrence_time) + "_segmentbuilder_after_crossing_v"
      + std::to_string(crossing->voronoi_vertex_id) + "_" + std::to_string(old_tri) + "_to_" + std::to_string(new_tri) + ".svg";
    const auto& containing_tri_ids = segment_builder_.kin_del.getCrossingData().getContainingTriIds();
    auto intersection_debug_data = segment_builder_.kin_del.getCrossingIntersectionDebugData();
    HalfEdgeDelaunayGraphToSVG::write(points, graph, filename, 0.1, &segment_builder_.kin_del.getFacesInside(), true,
      &containing_tri_ids, &intersection_debug_data);
    KINDS_INFO("SegmentBuilder wrote SVG: " << filename);
  }

  size_t voronoi_vertex_id = crossing->voronoi_vertex_id;
  glm::dvec3 voronoi_vertex_position = glm::dvec3(crossing->position, crossing->occurrence_time);
  auto half_edges = graph.getFaces()[voronoi_vertex_id].half_edges;
  if (!segment_builder_.kin_del.isOnComponentBoundary(crossing->half_edge_id))
  {
    return;
  }

  const size_t crossed_d_edge = crossing->half_edge_id / 2;
  const auto& crossing_data_after = segment_builder_.kin_del.getCrossingData();
  if (crossed_d_edge < crossing_data_after.delaunay_edge_intersections.size())
  {
    // Post-event intersection order along the crossed Delaunay edge.
    const auto& d_refs_after = crossing_data_after.delaunay_edge_intersections[crossed_d_edge];
    std::vector<size_t> after_voronoi_edge_ids;
    after_voronoi_edge_ids.reserve(d_refs_after.size());
    for (const auto& r : d_refs_after)
    {
      after_voronoi_edge_ids.push_back(r->voronoi_edge_id);
    }

    size_t component_id = segment_builder_.kin_del.component_data.component_map[graph.destination(crossing->half_edge_id)];
    auto& boundary_polygon = segment_builder_.kin_del.component_data.component_boundaries[component_id][0];
    auto centroid_local = segment_builder_.kin_del.component_data.component_centroids[component_id];
    const glm::dvec3 event_pos(crossing->position, crossing->occurrence_time);

    auto extend_pair_at_removed_intersection = [&](size_t pair_idx, bool remove_start_endpoint)
    {
      if (pair_idx == static_cast<size_t>(-1) || pair_idx >= segment_builder_.intersection_mesh_pair_metadata.size()
        || pair_idx >= segment_builder_.intersection_meshes.size()
        || pair_idx >= segment_builder_.intersection_mesh_pair_last_left_and_right_vertex.size())
      {
        return;
      }

      auto& mesh = segment_builder_.intersection_meshes[pair_idx];
      auto& segs = segment_builder_.intersection_mesh_pair_last_left_and_right_vertex[pair_idx];
      if (segs.empty())
      {
        return;
      }
      auto& seg = segs.front();
      const size_t new_vid = segment_builder_.addMeshletVertex(
        mesh, boundary_polygon, centroid_local, event_pos, graph.destination(crossing->half_edge_id), crossing->occurrence_time);
      const size_t last_left = static_cast<size_t>(seg.mesh_start_vertex_id);
      const size_t last_right = static_cast<size_t>(seg.mesh_end_vertex_id);
      // Finish/extend current strip piece up to the crossing event location.
      segment_builder_.addMeshletTriangle(mesh, last_left, last_right, new_vid);

      const auto& md = segment_builder_.intersection_mesh_pair_metadata[pair_idx];
      const bool is_outer = md.start_delaunay_edge_id == static_cast<size_t>(-1) || md.end_delaunay_edge_id == static_cast<size_t>(-1);
      if (is_outer)
      {
        // Outer interval terminates at this event (endpoint on the "null side"), so no active strip remains.
        segs.clear();
        return;
      }

      if (remove_start_endpoint)
      {
        // Inner interval: only one endpoint is replaced by the event vertex; the other endpoint keeps advancing.
        seg.mesh_start_vertex_id = static_cast<int>(new_vid);
        seg.start_half_edge_id = -1;
        seg.start_crossing.reset();
      }
      else
      {
        seg.mesh_end_vertex_id = static_cast<int>(new_vid);
        seg.end_half_edge_id = -1;
        seg.end_crossing.reset();
      }
    };

    for (const auto& old_ref : crossing_edge_snapshot_)
    {
      if (std::find(after_voronoi_edge_ids.begin(), after_voronoi_edge_ids.end(), old_ref.voronoi_edge_id)
        != after_voronoi_edge_ids.end())
      {
        continue;
      }
      if (old_ref.prev_pair_idx == static_cast<size_t>(-1) && old_ref.next_pair_idx == static_cast<size_t>(-1))
      {
        KINDS_WARNING("Boundary interval crossing removal has no linked mesh pair on de=" << crossed_d_edge
                                                                                           << " ve=" << old_ref.voronoi_edge_id
                                                                                           << " t=" << crossing->occurrence_time
                                                                                           << " (likely stale links after global recompute).");
      }
      // Along the crossed Delaunay edge: `prev` interval ends at removed intersection, `next` interval starts at it.
      extend_pair_at_removed_intersection(old_ref.prev_pair_idx, false);
      extend_pair_at_removed_intersection(old_ref.next_pair_idx, true);
    }

    std::vector<size_t> old_voronoi_edge_ids;
    old_voronoi_edge_ids.reserve(crossing_edge_snapshot_.size());
    for (const auto& s : crossing_edge_snapshot_)
    {
      old_voronoi_edge_ids.push_back(s.voronoi_edge_id);
    }

    for (const auto& new_ref : d_refs_after)
    {
      if (std::find(old_voronoi_edge_ids.begin(), old_voronoi_edge_ids.end(), new_ref->voronoi_edge_id) != old_voronoi_edge_ids.end())
      {
        continue;
      }
      const auto it_pos = std::find(after_voronoi_edge_ids.begin(), after_voronoi_edge_ids.end(), new_ref->voronoi_edge_id);
      if (it_pos == after_voronoi_edge_ids.end())
      {
        continue;
      }
      const size_t pos = static_cast<size_t>(std::distance(after_voronoi_edge_ids.begin(), it_pos));
      if (!(pos == 0 || pos + 1 == after_voronoi_edge_ids.size()))
      {
        KINDS_WARNING("Boundary interval crossing insertion is not at de-list endpoint on de=" << crossed_d_edge
                                                                                                << " ve="
                                                                                                << new_ref->voronoi_edge_id
                                                                                                << " pos=" << pos << "/"
                                                                                                << after_voronoi_edge_ids.size()
                                                                                                << " t=" << crossing->occurrence_time
                                                                                                << " (currently unsupported; kept for later incremental handling).");
        continue;
      }
      if (pos == 0 || pos + 1 == after_voronoi_edge_ids.size())
      {
        // Newly inserted boundary-end intersection creates a new outer interval [null,ref] or [ref,null].
        const bool at_front = pos == 0;
        const size_t cell = segment_builder_.determineVoronoiCellForBoundaryIntersectionInterval(
          crossed_d_edge, at_front ? std::nullopt : std::optional(new_ref), at_front ? std::optional(new_ref) : std::nullopt);
        segment_builder_.startNewMeshFromIntersections(
          cell, crossing->occurrence_time, at_front ? std::nullopt : std::optional(new_ref),
          at_front ? std::optional(new_ref) : std::nullopt, false);
      }
    }
  }

  size_t inside_he_id = crossing->half_edge_id;
  bool entering_boundary = false;
  if (segment_builder_.kin_del.isOnComponentBoundaryOutside(crossing->half_edge_id))
  {
    inside_he_id = crossing->half_edge_id ^ 1;
    entering_boundary = true;
  }

  size_t component_id = segment_builder_.kin_del.component_data.component_map[graph.destination(inside_he_id)];
  auto& boundary_polygon = segment_builder_.kin_del.component_data.component_boundaries[component_id][0];
  auto centroid = segment_builder_.kin_del.component_data.component_centroids[component_id];

  for (size_t voronoi_he_id : half_edges)
  {
    const size_t strip_voronoi_edge_id = (voronoi_he_id & ~1) / 2;
    auto& segment_mesh_pair_index = segment_builder_.half_edge_index_to_segment_mesh_pair_index[voronoi_he_id];
    VoronoiMesh& mesh = segment_builder_.meshes[segment_mesh_pair_index];

    auto& segments = segment_builder_.segment_mesh_pair_last_left_and_right_vertex[segment_mesh_pair_index];
    auto it = std::find_if(segments.begin(), segments.end(),
      [inside_he_id](const SegmentBuilder::MeshingData& data)
      {
        return data.end_half_edge_id == static_cast<int>(inside_he_id)
          || data.start_half_edge_id == static_cast<int>(inside_he_id);
      });
    if (it != segments.end())
    {
      if (!entering_boundary)
      {
        size_t new_vertex_index = segment_builder_.addMeshletVertex(mesh, boundary_polygon, centroid,
          voronoi_vertex_position, voronoi_vertex_id, crossing->occurrence_time,
          std::optional<size_t>(crossing->voronoi_vertex_id));
        int last_left = it->mesh_start_vertex_id;
        int last_right = it->mesh_end_vertex_id;
        {
          const size_t tris_before = mesh.getTriangleCount();
          segment_builder_.addMeshletTriangle(mesh, last_left, last_right, new_vertex_index);
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
          size_t new_vertex_index = segment_builder_.addMeshletVertex(mesh, boundary_polygon, centroid,
            voronoi_vertex_position, voronoi_vertex_id, crossing->occurrence_time,
            std::optional<size_t>(crossing->voronoi_vertex_id));
          int last_left = it->mesh_start_vertex_id;
          int last_right = it->mesh_end_vertex_id;
          {
            const size_t tris_before = mesh.getTriangleCount();
            segment_builder_.addMeshletTriangle(mesh, last_left, last_right, new_vertex_index);
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
          size_t new_vertex_index = segment_builder_.addMeshletVertex(mesh, boundary_polygon, centroid,
            voronoi_vertex_position, voronoi_vertex_id, crossing->occurrence_time,
            std::optional<size_t>(crossing->voronoi_vertex_id));
          int last_left = it->mesh_start_vertex_id;
          int last_right = it->mesh_end_vertex_id;
          {
            const size_t tris_before = mesh.getTriangleCount();
            segment_builder_.addMeshletTriangle(mesh, last_left, last_right, new_vertex_index);
            logCrossingMeshExtend(
              segment_builder_, voronoi_he_id, crossing->occurrence_time, mesh, "enter_boundary_end", tris_before);
          }
          it->mesh_end_vertex_id = static_cast<int>(new_vertex_index);
        }
      }
    }
    else
    {
      size_t even_id = voronoi_he_id & ~1;
      size_t odd_id = even_id + 1;
      size_t end_voronoi_vertex_id = graph.getHalfEdges()[odd_id].face;

      if (end_voronoi_vertex_id == voronoi_vertex_id)
      {
        if (entering_boundary)
        {
          size_t new_vertex_index = segment_builder_.addMeshletVertex(mesh, boundary_polygon, centroid,
            voronoi_vertex_position, voronoi_vertex_id, crossing->occurrence_time,
            std::optional<size_t>(crossing->voronoi_vertex_id));
          SegmentBuilder::MeshingData seg { static_cast<int>(new_vertex_index), static_cast<int>(new_vertex_index),
            static_cast<int>(inside_he_id), -1 };
          seg.start_crossing
            = segment_builder_.closingMeshFindVoronoiEdgeIntersection(strip_voronoi_edge_id, inside_he_id);
          segments.emplace_back(std::move(seg));
        }
        else
        {
          assert(segments.back().end_half_edge_id == -1);
          size_t new_vertex_index = segment_builder_.addMeshletVertex(mesh, boundary_polygon, centroid,
            voronoi_vertex_position, voronoi_vertex_id, crossing->occurrence_time,
            std::optional<size_t>(crossing->voronoi_vertex_id));
          int last_left = segments.back().mesh_start_vertex_id;
          int last_right = segments.back().mesh_end_vertex_id;
          {
            const size_t tris_before = mesh.getTriangleCount();
            segment_builder_.addMeshletTriangle(mesh, last_left, last_right, new_vertex_index);
            logCrossingMeshExtend(
              segment_builder_, voronoi_he_id, crossing->occurrence_time, mesh, "tail_close", tris_before);
          }
          segments.back().end_half_edge_id = static_cast<int>(inside_he_id);
          segments.back().end_crossing
            = segment_builder_.closingMeshFindVoronoiEdgeIntersection(strip_voronoi_edge_id, inside_he_id);
          segments.back().mesh_end_vertex_id = static_cast<int>(new_vertex_index);
        }
      }
      else
      {
        if (entering_boundary)
        {
          size_t new_vertex_index = segment_builder_.addMeshletVertex(mesh, boundary_polygon, centroid,
            voronoi_vertex_position, voronoi_vertex_id, crossing->occurrence_time,
            std::optional<size_t>(crossing->voronoi_vertex_id));
          SegmentBuilder::MeshingData seg { static_cast<int>(new_vertex_index), static_cast<int>(new_vertex_index), -1,
            static_cast<int>(inside_he_id) };
          seg.end_crossing
            = segment_builder_.closingMeshFindVoronoiEdgeIntersection(strip_voronoi_edge_id, inside_he_id);
          segments.emplace_front(std::move(seg));
        }
        else
        {
          assert(segments.front().start_half_edge_id == -1);
          size_t new_vertex_index = segment_builder_.addMeshletVertex(mesh, boundary_polygon, centroid,
            voronoi_vertex_position, voronoi_vertex_id, crossing->occurrence_time,
            std::optional<size_t>(crossing->voronoi_vertex_id));
          int last_left = segments.front().mesh_start_vertex_id;
          int last_right = segments.front().mesh_end_vertex_id;
          {
            const size_t tris_before = mesh.getTriangleCount();
            segment_builder_.addMeshletTriangle(mesh, last_left, last_right, new_vertex_index);
            logCrossingMeshExtend(
              segment_builder_, voronoi_he_id, crossing->occurrence_time, mesh, "head_close", tris_before);
          }
          segments.front().start_half_edge_id = static_cast<int>(inside_he_id);
          segments.front().start_crossing
            = segment_builder_.closingMeshFindVoronoiEdgeIntersection(strip_voronoi_edge_id, inside_he_id);
          segments.front().mesh_start_vertex_id = static_cast<int>(new_vertex_index);
        }
      }
    }

    for (auto& seg : segments)
    {
      segment_builder_.refreshMeshingDataCrossingRefs(seg, strip_voronoi_edge_id);
    }
  }
}
} // namespace kinDS


#include "VisualDebugHighlight.hpp"

#include <vector>

namespace kinDS
{
namespace
{
void processPendingUndirectedEdges(
  VisualDebugHighlight& highlight, const HalfEdgeDelaunayGraph& graph, std::vector<size_t>& edge_pending)
{
  while (!edge_pending.empty())
  {
    const size_t raw_he = edge_pending.back();
    edge_pending.pop_back();

    const size_t even_he = raw_he & ~static_cast<size_t>(1);
    const size_t undirected_edge_id = even_he / 2;
    if (!highlight.voronoi_edges.insert(undirected_edge_id).second)
    {
      continue;
    }

    highlight.directed_half_edges.insert(even_he);
    highlight.directed_half_edges.insert(even_he ^ 1);

    const auto& he = graph.getHalfEdges()[even_he];
    const auto& twin_he = graph.getHalfEdges()[even_he ^ 1];
    if (he.face != -1)
    {
      highlight.voronoi_vertices.insert(static_cast<size_t>(he.face));
    }
    if (twin_he.face != -1)
    {
      highlight.voronoi_vertices.insert(static_cast<size_t>(twin_he.face));
    }
    if (he.origin != -1)
    {
      highlight.delaunay_vertices.insert(static_cast<size_t>(he.origin));
    }
    if (twin_he.origin != -1)
    {
      highlight.delaunay_vertices.insert(static_cast<size_t>(twin_he.origin));
    }
  }
}

void flushPendingHighlightWork(VisualDebugHighlight& highlight, const HalfEdgeDelaunayGraph& graph,
  std::vector<size_t>& face_pending, std::vector<size_t>& edge_pending)
{
  while (!face_pending.empty() || !edge_pending.empty())
  {
    while (!face_pending.empty())
    {
      const size_t face_id = face_pending.back();
      face_pending.pop_back();
      if (!highlight.delaunay_faces.insert(face_id).second)
      {
        continue;
      }

      const std::array<int, 3> verts = graph.getTriangleVertexIndices(face_id);
      for (int v : verts)
      {
        if (v >= 0)
        {
          highlight.delaunay_vertices.insert(static_cast<size_t>(v));
        }
      }
      for (size_t he : graph.getFaces()[face_id].half_edges)
      {
        edge_pending.push_back(he);
      }
    }

    processPendingUndirectedEdges(highlight, graph, edge_pending);
  }
}
} // namespace

void VisualDebugHighlight::addUndirectedDelaunayEdge(const HalfEdgeDelaunayGraph& graph, size_t he_id)
{
  std::vector<size_t> edge_pending;
  edge_pending.push_back(he_id);
  std::vector<size_t> face_pending;
  flushPendingHighlightWork(*this, graph, face_pending, edge_pending);
}

void VisualDebugHighlight::addDelaunayTriangle(const HalfEdgeDelaunayGraph& graph, size_t face_id)
{
  std::vector<size_t> face_pending;
  face_pending.push_back(face_id);
  std::vector<size_t> edge_pending;
  flushPendingHighlightWork(*this, graph, face_pending, edge_pending);
}

bool VisualDebugHighlight::affectsDelaunayVertex(size_t vertex_id) const
{
  return delaunay_vertices.find(vertex_id) != delaunay_vertices.end();
}

bool VisualDebugHighlight::affectsDelaunayFace(size_t face_id) const
{
  return delaunay_faces.find(face_id) != delaunay_faces.end();
}

bool VisualDebugHighlight::affectsDirectedHalfEdge(size_t he_id) const
{
  return directed_half_edges.find(he_id) != directed_half_edges.end();
}

bool VisualDebugHighlight::affectsVoronoiVertex(size_t voronoi_vertex_id) const
{
  return voronoi_vertices.find(voronoi_vertex_id) != voronoi_vertices.end();
}

bool VisualDebugHighlight::affectsVoronoiEdge(size_t voronoi_edge_id) const
{
  return voronoi_edges.find(voronoi_edge_id) != voronoi_edges.end();
}

bool VisualDebugHighlight::affectsPrimaryVoronoiEdge(size_t voronoi_edge_id) const
{
  return primary_voronoi_edges.find(voronoi_edge_id) != primary_voronoi_edges.end();
}

bool VisualDebugHighlight::affectsCrossing(size_t delaunay_edge_id, size_t voronoi_edge_id) const
{
  return affectsVoronoiEdge(delaunay_edge_id) || affectsVoronoiEdge(voronoi_edge_id);
}

VisualDebugHighlight VisualDebugHighlight::forFlip(const HalfEdgeDelaunayGraph& graph, size_t flip_half_edge_id)
{
  VisualDebugHighlight highlight;
  std::vector<size_t> face_pending;
  std::vector<size_t> edge_pending;
  edge_pending.push_back(flip_half_edge_id);

  const auto& flip_he = graph.getHalfEdges()[flip_half_edge_id];
  const auto& flip_twin_he = graph.getHalfEdges()[flip_half_edge_id ^ 1];
  if (flip_he.face != -1)
  {
    face_pending.push_back(static_cast<size_t>(flip_he.face));
  }
  if (flip_twin_he.face != -1)
  {
    face_pending.push_back(static_cast<size_t>(flip_twin_he.face));
  }

  const std::array<size_t, 4> quad_hes = graph.getQuadBoundaryHalfEdgeIndices(flip_half_edge_id / 2);
  for (size_t quad_he : quad_hes)
  {
    edge_pending.push_back(quad_he);
  }

  flushPendingHighlightWork(highlight, graph, face_pending, edge_pending);
  highlight.primary_voronoi_edges.insert(flip_half_edge_id / 2);
  return highlight;
}

VisualDebugHighlight VisualDebugHighlight::forRadius(const HalfEdgeDelaunayGraph& graph, size_t radius_half_edge_id)
{
  VisualDebugHighlight highlight;
  const int triangle_face = graph.getHalfEdges()[radius_half_edge_id].face;
  if (triangle_face != -1)
  {
    std::vector<size_t> face_pending;
    face_pending.push_back(static_cast<size_t>(triangle_face));
    std::vector<size_t> edge_pending;
    flushPendingHighlightWork(highlight, graph, face_pending, edge_pending);
  }
  return highlight;
}

VisualDebugHighlight VisualDebugHighlight::forCrossing(
  const HalfEdgeDelaunayGraph& graph, size_t crossed_half_edge_id, size_t voronoi_vertex_id)
{
  VisualDebugHighlight highlight;
  highlight.voronoi_vertices.insert(voronoi_vertex_id);

  std::vector<size_t> face_pending;
  std::vector<size_t> edge_pending;
  edge_pending.push_back(crossed_half_edge_id);

  const int old_face = graph.getHalfEdges()[crossed_half_edge_id].face;
  const int new_face = graph.getHalfEdges()[crossed_half_edge_id ^ 1].face;
  if (old_face != -1)
  {
    face_pending.push_back(static_cast<size_t>(old_face));
  }
  if (new_face != -1)
  {
    face_pending.push_back(static_cast<size_t>(new_face));
  }

  for (size_t incident_he : graph.getFaces()[voronoi_vertex_id].half_edges)
  {
    edge_pending.push_back(incident_he);
  }

  flushPendingHighlightWork(highlight, graph, face_pending, edge_pending);
  return highlight;
}

VisualDebugHighlight VisualDebugHighlight::forSubdivisionStrand(const HalfEdgeDelaunayGraph& graph, size_t strand_vertex_id)
{
  VisualDebugHighlight highlight;
  std::vector<size_t> edge_pending;
  for (auto it = graph.incidentEdgesBegin(strand_vertex_id); it != graph.incidentEdgesEnd(strand_vertex_id); ++it)
  {
    edge_pending.push_back(*it);
  }
  highlight.delaunay_vertices.insert(strand_vertex_id);

  std::vector<size_t> face_pending;
  flushPendingHighlightWork(highlight, graph, face_pending, edge_pending);
  return highlight;
}

VisualDebugHighlight VisualDebugHighlight::forSectionBoundary(
  const HalfEdgeDelaunayGraph& graph, const std::function<bool(size_t even_half_edge_id)>& is_boundary_edge)
{
  VisualDebugHighlight highlight;
  std::vector<size_t> edge_pending;
  edge_pending.reserve(graph.getHalfEdges().size() / 6);
  for (size_t even_he = 0; even_he < graph.getHalfEdges().size(); even_he += 2)
  {
    if (is_boundary_edge(even_he))
    {
      edge_pending.push_back(even_he);
    }
  }

  std::vector<size_t> face_pending;
  flushPendingHighlightWork(highlight, graph, face_pending, edge_pending);
  return highlight;
}

} // namespace kinDS

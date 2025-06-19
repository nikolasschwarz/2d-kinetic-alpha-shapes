pragma once

class VoronoiGraph
{
public:
	VoronoiGraph();
	~VoronoiGraph();

	struct Site
	{
    float x, y; // Coordinates of the site
    int index; // Site number
    std::vector<size_t> edges; // Indices of edges connected to this site
	};

  struct VoronoiVertex
  {
    float x, y; // Coordinates of the vertex
    bool is_infinite; // Flag to indicate if the vertex is infinite
    int index; // Vertex number
    //std::vector<size_t> edges; // Indices of edges connected to this vertex
  };

  struct HalfEdge
  {
    size_t origin; // Index of the origin vertex
    size_t destination; // Index of the destination vertex
    size_t twin; // Index of the twin half-edge
    size_t next; // Index of the next half-edge in the face
    size_t site_index; // Index of the site this half-edge belongs to
  };

  std::vector<Site> sites;
  std::vector<VoronoiVertex> vertices;
  std::vector<HalfEdge> halfEdges;

private:

};

VoronoiGraph::VoronoiGraph()
{
}

VoronoiGraph::~VoronoiGraph()
{
}
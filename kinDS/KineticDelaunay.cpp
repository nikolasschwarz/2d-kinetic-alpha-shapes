#include "KineticDelaunay.hpp"

using namespace kinDS;

/* Compare to Leonidas Guibas and Jorge Stolfi. 1985. Primitives for the manipulation of general subdivisions and the
 * computation of Voronoi. ACM Trans. Graph. 4, 2 (April 1985), 74-123. https://doi.org/10.1145/282918.282923
 */
static Polynomial inCircle(const Polynomial& ax, const Polynomial& ay, const Polynomial& bx, const Polynomial& by,
  const Polynomial& cx, const Polynomial& cy, const Polynomial& px, const Polynomial& py)
{
  const Polynomial dx = ax - px;
  const Polynomial dy = ay - py;
  const Polynomial ex = bx - px;
  const Polynomial ey = by - py;
  const Polynomial fx = cx - px;
  const Polynomial fy = cy - py;

  const Polynomial ap = dx * dx + dy * dy;
  const Polynomial bp = ex * ex + ey * ey;
  const Polynomial cp = fx * fx + fy * fy;

  return (dx * (ey * cp - bp * fy) - dy * (ex * cp - bp * fx) + ap * (ex * fy - ey * fx));
}

static Polynomial ccw(const Polynomial& ax, const Polynomial& ay, const Polynomial& bx, const Polynomial& by,
  const Polynomial& cx, const Polynomial& cy)
{
  return (ax * by) + (bx * cy) + (cx * ay) - (ay * bx) - (by * cx) - (cy * ax);
}

// Polynomial that evaluates to zero iff the distance from A to the circumcenter equals the value r
static Polynomial circumradiusEquals(const Polynomial& ax, const Polynomial& ay, const Polynomial& bx,
  const Polynomial& by, const Polynomial& cx, const Polynomial& cy, double r)
{
  // We first do the same computations as for the circumcenter
  Polynomial D = (ax * (by - cy) + bx * (cy - ay) + cx * (ay - by)) * 2.0;

  // only compute the numerators
  Polynomial Nx = ((ax * ax + ay * ay) * (by - cy) + (bx * bx + by * by) * (cy - ay) + (cx * cx + cy * cy) * (ay - by));
  Polynomial Ny = ((ax * ax + ay * ay) * (cx - bx) + (bx * bx + by * by) * (ax - cx) + (cx * cx + cy * cy) * (bx - ax));

  // By taking the distance formula between the circumcenter and the point A and setting it equal to r, we get the
  // following after rearranging:
  Polynomial circumradius_eq = (Nx - ax * D) * (Nx - ax * D) + (Ny - ay * D) * (Ny - ay * D) - (D * D * (r * r));
  return circumradius_eq;
}

static double circumradius(const glm::dvec2& p0, const glm::dvec2& p1, const glm::dvec2& p2)
{
  const double x0 = p0[0], y0 = p0[1];
  const double x1 = p1[0], y1 = p1[1];
  const double x2 = p2[0], y2 = p2[1];

  // Side lengths
  const double a = std::hypot(x1 - x2, y1 - y2);
  const double b = std::hypot(x0 - x2, y0 - y2);
  const double c = std::hypot(x0 - x1, y0 - y1);

  // Twice the triangle area (cross product magnitude)
  const double area2 = std::abs((x1 - x0) * (y2 - y0) - (y1 - y0) * (x2 - x0));

  if (area2 == 0.0)
  {
    throw std::runtime_error("Degenerate triangle: circumradius undefined");
  }

  // R = (a * b * c) / (4 * A), and area2 = 2 * A
  return (a * b * c) / (2.0 * area2);
}

void KineticDelaunay::computeRadiusEvents(double t, size_t he_id)
{
  if (cutoff == std::numeric_limits<double>::infinity())
  {
    // no radius events wanted
    return;
  }

  const size_t section = static_cast<size_t>(t);
  const float fraction = t - section;

  size_t face_id = graph.getHalfEdges()[he_id].face;
  size_t u = graph.getHalfEdges()[he_id].origin;
  size_t v = graph.destination(he_id);
  size_t w = graph.triangleOppositeVertex(he_id);

  if (u == -1 || v == -1 || w == -1)
  {
    // one of the vertices is at infinity, no event possible
    return;
  }

  std::vector<Trajectory<2>> trajs;

  trajs.push_back(branch_trajs.getPiecePolynomial(u, section));
  trajs.push_back(branch_trajs.getPiecePolynomial(v, section));
  trajs.push_back(branch_trajs.getPiecePolynomial(w, section));

  Polynomial event_trigger
    = circumradiusEquals(trajs[0][0], trajs[0][1], trajs[1][0], trajs[1][1], trajs[2][0], trajs[2][1], cutoff);

  event_trigger.trim();
  auto zeros = event_trigger.realRoots();

  // print roots:
  for (const auto& root : zeros)
  {
    if (isnan(root))
    {
      continue; // Skip NaN roots
    }
    if (root > fraction && root <= 1)
    { // Check if the root is within the valid range
      double event_time = root + section;
      // std::cout << "Root found at t = " << event_time << std::endl;

      glm::dvec2 center {};

      for (const auto& traj : trajs)
      {
        center[0] += traj[0](root);
        center[1] += traj[1](root);
      }
      center[0] /= trajs.size();
      center[1] /= trajs.size();
      KINDS_DEBUG("Boundary Event at time " << event_time << " for half-edge ID " << he_id << " at center position"
                                            << glm::to_string(center));

      events.emplace(
        Event(event_time, he_id, t, center, Event::RADIUS)); // Store the event with the time and half-edge index
    }
  }
}

void KineticDelaunay::computeFlipEvents(double t, size_t quad_id)
{
  const size_t section = static_cast<size_t>(t);
  const float fraction = t - section;

  size_t he_id = quad_id * 2;
  Polynomial event_trigger;

  std::vector<Trajectory<2>> trajs;

  if (graph.isOnConvexBoundary(he_id) || graph.isOutsideConvexBoundary(he_id))
  {
    // boundary edges must be treated separately using ccw

    // need to get the inner half-edge so we have access to the triangle
    // in case both are outside, this swap does not matter, so we just let it happen
    if (graph.isOutsideConvexBoundary(he_id))
    {
      he_id = he_id ^ 1; // use the twin half-edge if the current one is on the boundary
    }

    // Depending on the half-edge, the infinite vertex could be in different places, so we just collect all and filter
    // it out
    int indices[4];
    indices[0] = graph.getHalfEdges()[he_id].origin; // First vertex
    indices[1] = graph.triangleOppositeVertex(he_id ^ 1); // Second vertex
    indices[2] = graph.getHalfEdges()[he_id ^ 1].origin; // Third vertex
    indices[3] = graph.triangleOppositeVertex(he_id); // Fourth vertex

    std::vector<int> filtered_indices;

    std::copy_if(indices, indices + 4, std::back_inserter(filtered_indices), [this](int index) { return index != -1; });

    int& a = filtered_indices[0]; // First vertex
    int& b = filtered_indices[1]; // Second vertex
    int& c = filtered_indices[2]; // Third vertex

    // print the triangle vertices:
    // std::cout << "Triangle vertices: " << a << ", " << b << ", " << c << std::endl;

    trajs.push_back(branch_trajs.getPiecePolynomial(a, section));
    trajs.push_back(branch_trajs.getPiecePolynomial(b, section));
    trajs.push_back(branch_trajs.getPiecePolynomial(c, section));

    event_trigger = ccw(trajs[0][0], trajs[0][1], trajs[1][0], trajs[1][1], trajs[2][0], trajs[2][1]);
  }
  else
  {
    int a = graph.getHalfEdges()[he_id].origin; // First vertex
    int b = graph.triangleOppositeVertex(he_id ^ 1); // Second vertex
    int c = graph.getHalfEdges()[he_id ^ 1].origin; // Third vertex
    int d = graph.triangleOppositeVertex(he_id); // Fourth vertex

    // print the quadrilateral vertices:
    // std::cout << "Quadrilateral vertices: " << a << ", " << b << ", " << c << ", " << d << std::endl;

    trajs.push_back(branch_trajs.getPiecePolynomial(a, section));
    trajs.push_back(branch_trajs.getPiecePolynomial(b, section));
    trajs.push_back(branch_trajs.getPiecePolynomial(c, section));
    trajs.push_back(branch_trajs.getPiecePolynomial(d, section));

    event_trigger = inCircle(
      trajs[0][0], trajs[0][1], trajs[1][0], trajs[1][1], trajs[2][0], trajs[2][1], trajs[3][0], trajs[3][1]);
  }
  event_trigger.trim();
  auto zeros = event_trigger.realRoots();

  for (const auto& root : zeros)
  {
    if (isnan(root))
    {
      continue; // Skip NaN roots
    }

    if (root > fraction && root <= 1)
    {
      // Check if the root is within the valid range
      double event_time = root + section;
      // std::cout << "Root found at t = " << event_time << std::endl;

      glm::dvec2 center {};

      for (const auto& traj : trajs)
      {
        center[0] += traj[0](root);
        center[1] += traj[1](root);
      }
      center[0] /= trajs.size();
      center[1] /= trajs.size();

      KINDS_DEBUG("Swap Event at time " << event_time << " for half-edge ID " << he_id << " at center position "
                                        << glm::to_string(center));

      events.emplace(
        Event(event_time, he_id, t, center, Event::FLIP)); // Store the event with the time and half-edge index
    }
  }
}

glm::dvec3 kinDS::KineticDelaunay::computeVoronoiVertexHomogenous(size_t voronoi_vertex_id, double t) const
{
  auto face = graph.getFaces()[voronoi_vertex_id];
  auto vertices = graph.adjacentTriangleVertices(face.half_edges[0]);

  std::vector<glm::dvec2> points;
  for (size_t i = 0; i < vertices.size(); i++)
  {
    if (vertices[i] != -1)
    {
      points.push_back(getPointAt(vertices[i], t));
    }
  }

  if (points.size() < 3)
  {
    // Voronoi vertex at infinity, return a point at infinity in homogeneous coordinates
    glm::dvec2 dir = glm::normalize(points[1] - points[0]);

    if(vertices[1] == -1)
    {
      dir = -dir;
    }
    return glm::dvec3(dir, 0.0);
  }

  auto circumcenter = graph.circumcenter(points[0], points[1], points[2]);
  return glm::dvec3(circumcenter, 1.0);
}

void kinDS::KineticDelaunay::computeCrossingEvents(double t, size_t voronoi_vertex_id)
{
  if (!on_the_fly_boundary)
  {
    return;
  }
  const size_t section = static_cast<size_t>(t);
  const float fraction = t - section;

  auto& dual_triangle = graph.getFaces()[voronoi_vertex_id];
  auto& containing_triangle = graph.getFaces()[crossing_data.getContainingTriId(voronoi_vertex_id)];

  // compute polynomials of two bisectors in homogeneous coordinates
  size_t v_i = graph.getHalfEdges()[containing_triangle.half_edges[0]].origin;
  size_t v_j = graph.getHalfEdges()[containing_triangle.half_edges[1]].origin;
  size_t v_k = graph.getHalfEdges()[containing_triangle.half_edges[2]].origin;

  // If a vertex is infinite, so is the Voronoi vertex and it cannot cross any edge, so we can skip this event.
  if (v_i == -1 || v_j == -1 || v_k == -1)
  {
    return;
  }

  Trajectory<2> traj_i = branch_trajs.getPiecePolynomial(v_i, section);
  Trajectory<2> traj_j = branch_trajs.getPiecePolynomial(v_j, section);
  Trajectory<2> traj_k = branch_trajs.getPiecePolynomial(v_k, section);

  Trajectory<3> bisector_ij;

  bisector_ij[0] = 2 * (traj_j[0] - traj_i[0]);
  bisector_ij[1] = 2 * (traj_j[1] - traj_i[1]);
  bisector_ij[2] = (traj_i[0] * traj_i[0] + traj_i[1] * traj_i[1]) - (traj_j[0] * traj_j[0] + traj_j[1] * traj_j[1]);

  Trajectory<3> bisector_ik;

  bisector_ik[0] = 2 * (traj_k[0] - traj_i[0]);
  bisector_ik[1] = 2 * (traj_k[1] - traj_i[1]);
  bisector_ik[2] = (traj_i[0] * traj_i[0] + traj_i[1] * traj_i[1]) - (traj_k[0] * traj_k[0] + traj_k[1] * traj_k[1]);

  // We only need the first event as any following events will be invalidated by the first crossing event.
  // TODO: The exception is the edge being crossed, but that would make this more complex. We can optimize this later if
  // needed.
  double event_time = std::numeric_limits<double>::infinity();
  size_t event_he_id = -1;

  // Construct polynomial predicates for each of the three edges of the containing triangle
  for (size_t edge_index = 0; edge_index < 3; edge_index++)
  {
    size_t he_id = containing_triangle.half_edges[edge_index];
    size_t a = graph.getHalfEdges()[he_id].origin;
    size_t b = graph.getHalfEdges()[he_id ^ 1].origin;

    Trajectory<3> line_ab;
    if (a != -1 && b != -1)
    {

      Trajectory<2> traj_a = branch_trajs.getPiecePolynomial(a, section);
      Trajectory<2> traj_b = branch_trajs.getPiecePolynomial(b, section);

      // line through a and b in homogeneous coordinates

      line_ab[0] = traj_a[1] - traj_b[1];
      line_ab[1] = traj_b[0] - traj_a[0];
      line_ab[2] = traj_a[0] * traj_b[1] - traj_a[1] * traj_b[0];
    }
    else
    {
      // Create a line through the finite vertex that is perpendicular to the line through the two neighboring vertices
      // on the convex hull.
      size_t finite_vertex = (a != -1) ? a : b;

      if (a == finite_vertex)
      {
        size_t prev_he_id = graph.prev(he_id);
        size_t next_he_id = graph.getHalfEdges()[he_id ^ 1].next;

        size_t c = graph.getHalfEdges()[prev_he_id].origin;
        size_t c_prime = graph.getHalfEdges()[next_he_id].origin;

        Trajectory<2> traj_a = branch_trajs.getPiecePolynomial(a, section);
        Trajectory<2> traj_c = branch_trajs.getPiecePolynomial(c, section);
        Trajectory<2> traj_c_prime = branch_trajs.getPiecePolynomial(c_prime, section);

        line_ab[0] = traj_c_prime[0] - traj_c[0];
        line_ab[1] = traj_c_prime[1] - traj_c[1];
        line_ab[2] = -(traj_c_prime[0] - traj_c[0]) * traj_a[0] - (traj_c_prime[1] - traj_c[1]) * traj_a[1];
      }
      else
      {
        size_t prev_he_id = graph.prev(he_id ^ 1);
        size_t next_he_id = graph.getHalfEdges()[he_id].next;

        size_t c_prime = graph.getHalfEdges()[prev_he_id].origin;
        size_t c = graph.getHalfEdges()[next_he_id].origin;

        Trajectory<2> traj_b = branch_trajs.getPiecePolynomial(b, section);
        Trajectory<2> traj_c = branch_trajs.getPiecePolynomial(c, section);
        Trajectory<2> traj_c_prime = branch_trajs.getPiecePolynomial(c_prime, section);

        line_ab[0] = traj_c_prime[0] - traj_c[0];
        line_ab[1] = traj_c_prime[1] - traj_c[1];
        line_ab[2] = -(traj_c_prime[0] - traj_c[0]) * traj_b[0] - (traj_c_prime[1] - traj_c[1]) * traj_b[1];
      }
    }

    // now compute the determinant of the matrix with bisector_ij, bisector_ik and line_ab as columns
    Polynomial event_trigger = bisector_ij[0] * bisector_ik[1] * line_ab[2]
      + bisector_ij[1] * bisector_ik[2] * line_ab[0] + bisector_ij[2] * bisector_ik[0] * line_ab[1]
      - bisector_ij[2] * bisector_ik[1] * line_ab[0] - bisector_ij[1] * bisector_ik[0] * line_ab[2]
      - bisector_ij[0] * bisector_ik[2] * line_ab[1];

    event_trigger.trim();
    auto zeros = event_trigger.realRoots();

    for (const auto& root : zeros)
    {
      if (isnan(root))
      {
        continue; // Skip NaN roots
      }

      if (root > fraction && root <= 1)
      {
        event_time = root + section;
        event_he_id = he_id;
      }
    }
  }

  if (event_time != std::numeric_limits<double>::infinity())
  {
    glm::dvec3 position_homogeneous;

    // use cross product to compute the intersection point of the two bisectors at the event time
    position_homogeneous[0] = bisector_ij[1](event_time) * bisector_ik[2](event_time)
      - bisector_ij[2](event_time) * bisector_ik[1](event_time);
    position_homogeneous[1] = bisector_ij[2](event_time) * bisector_ik[0](event_time)
      - bisector_ij[0](event_time) * bisector_ik[2](event_time);
    position_homogeneous[2] = bisector_ij[0](event_time) * bisector_ik[1](event_time)
      - bisector_ij[1](event_time) * bisector_ik[0](event_time);

    glm::dvec2 position(
      position_homogeneous.x / position_homogeneous.z, position_homogeneous.y / position_homogeneous.z);
    KINDS_DEBUG("Crossing Event at time " << event_time << " for Voronoi vertex ID " << voronoi_vertex_id
                                          << " crossing half-edge ID " << event_he_id << " at position "
                                          << glm::to_string(position));
    events.emplace(event_time, event_he_id, t, position, voronoi_vertex_id, Event::CROSSING);
  }
}

void KineticDelaunay::reassignVoronoiVerticesInQuadrilateral(size_t quad_index, double t)
{
  size_t he_id = quad_index * 2;
  size_t face_id0 = graph.getHalfEdges()[he_id].face;
  size_t face_id1 = graph.getHalfEdges()[he_id ^ 1].face;

  size_t u = graph.getHalfEdges()[he_id].origin;
  size_t v = graph.destination(he_id);

  glm::dvec2 edge_vector;
  glm::dvec2 pu, pv;
  // if a vertex is infinite, we need to compute a placeholder vector. Just as for the predicates we choose the vector
  // perpendicular to the line through the two neighboring vertices on the convex hull.
  if (u == -1 || v == -1)
  {
    size_t opposite0 = graph.triangleOppositeVertex(he_id);
    size_t opposite1 = graph.triangleOppositeVertex(he_id ^ 1);

    glm::vec2 p0 = getPointAt(opposite0, t);
    glm::vec2 p1 = getPointAt(opposite1, t);
    edge_vector = glm::normalize(glm::vec2(p1.y - p0.y, p0.x - p1.x));

    if (u == -1)
    {
      pv = getPointAt(v, t);
      pu = pv - edge_vector; // placeholder position for the infinite vertex
    }
    else
    {
      pu = getPointAt(u, t);
      pv = pu + edge_vector; // placeholder position for the infinite vertex
    }
  }
  else
  {
    pu = getPointAt(u, t);
    pv = getPointAt(v, t);

    edge_vector = pv - pu;
  }

  auto vertices0 = crossing_data.getVoronoiVerticesInTri(face_id0);
  auto vertices1 = crossing_data.getVoronoiVerticesInTri(face_id1);

  // TODO: what if the edge function evaluates to 0. Perhaps we should look at the derivative.
  for (size_t voronoi_vertex : vertices0)
  {
    // Compute the position of the Voronoi vertex at time t and check on which side of the edge it is.
    glm::dvec3 voronoi_pos = computeVoronoiVertexHomogenous(voronoi_vertex, t);
    if(voronoi_pos.z != 0)
    {
      glm::dvec2 event_pos = glm::dvec2(voronoi_pos.x / voronoi_pos.z, voronoi_pos.y / voronoi_pos.z);
      if (glm::cross(event_pos - pu, edge_vector) < 0)
      {
        crossing_data.moveVertex(voronoi_vertex, face_id1);
      }
    } else {
      // Must belong to the dual triangle
        crossing_data.moveVertex(voronoi_vertex, voronoi_vertex);
    }
  }

  for (size_t voronoi_vertex : vertices1)
  {
    // Compute the position of the Voronoi vertex at time t and check on which side of the edge it is.
    glm::dvec3 voronoi_pos = computeVoronoiVertexHomogenous(voronoi_vertex, t);
    if(voronoi_pos.z != 0)
    {
      glm::dvec2 event_pos = glm::dvec2(voronoi_pos.x / voronoi_pos.z, voronoi_pos.y / voronoi_pos.z);
      if (glm::cross(event_pos - pu, edge_vector) < 0)
      {
        crossing_data.moveVertex(voronoi_vertex, face_id1);
      }
    } else {
      // Must belong to the dual triangle
        crossing_data.moveVertex(voronoi_vertex, voronoi_vertex);
    }
  }

  // Recompute all crossing events
  for (size_t voronoi_vertex : vertices0)
  {
    computeCrossingEvents(t, voronoi_vertex);
  }
  for (size_t voronoi_vertex : vertices1)
  {
    computeCrossingEvents(t, voronoi_vertex);
  }
}

void KineticDelaunay::precomputeStep(double t)
{
  // TODO: Make sure it works where no change of sign occurs in the polynomial, i.e., roots that do not lead to a
  // change in the triangulation.
  size_t quad_count = graph.getHalfEdges().size() / 2;
  for (size_t i = 0; i < quad_count; i++)
  {
    computeFlipEvents(t, i);
  }

  size_t he_count = graph.getHalfEdges().size();
  for (size_t i = 0; i < face_inside.size(); i++)
  {
    size_t he_id = graph.getFaces()[i].half_edges[0];
    computeRadiusEvents(t, he_id);
  }

  for (size_t tri_id = 0; tri_id < graph.getFaces().size(); tri_id++)
  {
    computeCrossingEvents(t, tri_id);
  }
}

void KineticDelaunay::handleFlipEvent(EventHandler& event_handler, Event& event)
{
  // Check if the event is still valid
  if (event.creation_time < quadrilateral_last_updated[event.half_edge_id / 2])
  {
    // This event is outdated, skip it
    return;
  }

  // Process the event at the given time
  size_t face_id = graph.getHalfEdges()[event.half_edge_id].face;
  size_t twin_face_id = graph.getHalfEdges()[event.half_edge_id ^ 1].face;
  KINDS_DEBUG("Processing swap event at time " << event.time << " for half-edge ID " << event.half_edge_id
                                               << ". Faces inside " << face_inside[face_id] << " | "
                                               << face_inside[twin_face_id]);

  /*kinDS::HalfEdgeDelaunayGraphToSVG::write(
    getPointsAt(event.time), getGraph(), "test_" + std::to_string(event.time) + "_before.svg", 0.1);
  std::cout << "Wrote " << ("test_" + std::to_string(event.time) + "_before.svg") << std::endl;*/

  // Call the event handler if provided
  event_handler.beforeFlipEvent(event);

  // Faces swapped to the inside start out with an infinite circumradius, therefore their state depends on the cutoff
  if (graph.getHalfEdges()[event.half_edge_id].origin == -1)
  {
    //KINDS_DEBUG("Swapping face of half-edge " << event.half_edge_id << " to the inside at t = " << event.time);
    face_inside[twin_face_id] = (cutoff == std::numeric_limits<double>::infinity());
  }

  if (graph.getHalfEdges()[event.half_edge_id ^ 1].origin == -1)
  {
    //KINDS_DEBUG(
    //  "Swapping face of twin half-edge " << (event.half_edge_id ^ 1) << " to the inside at t = " << event.time);
    face_inside[face_id] = (cutoff == std::numeric_limits<double>::infinity());
  }

  //KINDS_DEBUG("Pre-flip: " << event.time << " for half-edge ID " << event.half_edge_id << ". Faces inside "
  //                         << face_inside[face_id] << " | " << face_inside[twin_face_id]);

  graph.flipEdge(event.half_edge_id);

  //KINDS_DEBUG("Post-flip:  " << event.time << " for half-edge ID " << event.half_edge_id << ". Faces inside "
  //                             << face_inside[face_id] << " | " << face_inside[twin_face_id]);

  // one of the triangles might have been swapped outside
  auto tri_verts1 = graph.adjacentTriangleVertices(event.half_edge_id);

  for (auto& v : tri_verts1)
  {
    if (v == -1)
    {
      size_t face_id = graph.getHalfEdges()[event.half_edge_id].face;
      //KINDS_DEBUG("Swapped face " << face_id << " of half-edge " << event.half_edge_id
      //                            << " to the outside at t = " << event.time);
      setFaceInside(face_id, false);
    }
  }

  auto tri_verts2 = graph.adjacentTriangleVertices(event.half_edge_id ^ 1);
  for (auto& v : tri_verts2)
  {
    if (v == -1)
    {
      size_t face_id = graph.getHalfEdges()[event.half_edge_id ^ 1].face;
      //KINDS_DEBUG("Swapped face " << face_id << " of half-edge " << (event.half_edge_id ^ 1)
      //                            << " to the outside at t = " << event.time);
      setFaceInside(face_id, false);
    }
  }

  //KINDS_DEBUG("Processed swap event at time " << event.time << " for half-edge ID " << event.half_edge_id
  //                                            << ". Faces inside " << face_inside[face_id] << " | "
  //                                            << face_inside[twin_face_id]);

  /*kinDS::HalfEdgeDelaunayGraphToSVG::write(
    getPointsAt(event.time), getGraph(), "test_" + std::to_string(event.time) + ".svg", 0.1, &face_inside);
  std::cout << "Wrote " << ("test_" + std::to_string(event.time) + ".svg") << std::endl;*/

  // After flipping the edge, we need to recompute the events for all surrounding half-edges
  size_t next1 = graph.getHalfEdges()[event.half_edge_id].next;
  size_t next2 = graph.getHalfEdges()[next1].next;

  size_t twin_next1 = graph.getHalfEdges()[event.half_edge_id ^ 1].next;
  size_t twin_next2 = graph.getHalfEdges()[twin_next1].next;

  computeFlipEvents(event.time, next1 / 2);
  quadrilateral_last_updated[next1 / 2] = event.time; // Update the last updated time for the quadrilateral

  computeFlipEvents(event.time, next2 / 2);
  quadrilateral_last_updated[next2 / 2] = event.time; // Update the last updated time for the quadrilateral

  computeFlipEvents(event.time, twin_next1 / 2);
  quadrilateral_last_updated[twin_next1 / 2] = event.time; // Update the last updated time for the quadrilateral

  computeFlipEvents(event.time, twin_next2 / 2);
  quadrilateral_last_updated[twin_next2 / 2] = event.time; // Update the last updated time for the quadrilateral

  // re-compute radius events for both triangles
  computeRadiusEvents(event.time, event.half_edge_id);
  face_last_updated[face_id] = event.time;

  computeRadiusEvents(event.time, event.half_edge_id ^ 1);
  face_last_updated[twin_face_id] = event.time;

  // trigger re-assignment of voronoi vertices needed for crossing events
  reassignVoronoiVerticesInQuadrilateral(event.half_edge_id / 2, event.time);

  event_handler.afterFlipEvent(event); // Call the event handler after processing the event
}

void KineticDelaunay::handleRadiusEvent(EventHandler& event_handler, Event& event)
{
  assert(event.type == Event::RADIUS);

  // Check if the event is still valid
  size_t face_id = graph.getHalfEdges()[event.half_edge_id].face;
  if (event.creation_time < face_last_updated[face_id])
  {
    // This event is outdated, skip it
    return;
  }

  // Process the event at the given time
  //KINDS_DEBUG("Processing radius event at time " << event.time << " for half-edge ID " << event.half_edge_id);
  /*kinDS::HalfEdgeDelaunayGraphToSVG::write(
    getPointsAt(event.time), getGraph(), "test_" + std::to_string(event.time) + "_before.svg", 0.1, &face_inside);
  std::cout << "Wrote " << ("test_" + std::to_string(event.time) + "_before.svg") << std::endl;*/
  // Call the event handler if provided
  // TODO: (probably in callback) Handle boundary.
  event_handler.beforeRadiusEvent(event);

  setFaceInside(face_id, !face_inside[face_id]);

  event_handler.afterRadiusEvent(event);
  /*kinDS::HalfEdgeDelaunayGraphToSVG::write(
    getPointsAt(event.time), getGraph(), "test_" + std::to_string(event.time) + ".svg", 0.1, &face_inside);
  std::cout << "Wrote " << ("test_" + std::to_string(event.time) + ".svg") << std::endl;*/
}
void KineticDelaunay::handleCrossingEvent(EventHandler& event_handler, Event& event)
{
  // Check if the event is still valid
  // TODO: I think this is actually redundant, perhaps remove it.
  if (event.creation_time < crossing_data.last_crossing[event.voronoi_vertex_id])
  {
    // This event is outdated, skip it
    return;
  }

  size_t containing_tri_id = crossing_data.getContainingTriId(event.voronoi_vertex_id);
  // The event is also outdated if the face has been updated in a flip event
  if (event.creation_time < face_last_updated[containing_tri_id])
  {
    return;
  }

  event_handler.beforeCrossingEvent(event);

  crossing_data.last_crossing[event.voronoi_vertex_id]
    = event.time; // Update the last crossing time for this Voronoi vertex

  // move to neighboring triangle
  crossing_data.moveVertex(event.voronoi_vertex_id, graph.getHalfEdges()[event.half_edge_id ^ 1].face);

  event_handler.afterCrossingEvent(event);

  // Re-compute crossing events for this Voronoi vertex
  computeCrossingEvents(event.time, event.voronoi_vertex_id);
}

void KineticDelaunay::handleEvents(EventHandler& event_handler)
{

  while (!events.empty())
  {
    Event event = events.top();
    events.pop();

    switch (event.type)
    {
    case Event::FLIP:
      handleFlipEvent(event_handler, event);
      break;
    case Event::RADIUS:
      handleRadiusEvent(event_handler, event);
      break;
    case Event::CROSSING:
      handleCrossingEvent(event_handler, event);
      break;
    }
  }
}

size_t KineticDelaunay::getBranchIndex(size_t strand_id, size_t t) const
{
  return branch_trajs.getBranchIndex(strand_id, t);
}

const std::vector<std::vector<size_t>>& KineticDelaunay::getBranches(size_t t) const
{
  return branch_trajs.getStrandsByBranchId()[t];
}

const std::vector<size_t>& KineticDelaunay::getBranchStrands(size_t t, size_t branch_id)
{
  return branch_trajs.getStrandsByBranchId()[t][branch_id];
}

KineticDelaunay::KineticDelaunay(const StrandTree& branch_trajs, double cutoff, bool add_dummy_splines)
  : branch_trajs(branch_trajs)
  , cutoff(cutoff)
  , add_dummy_boundary(add_dummy_splines)
{
  if (add_dummy_splines)
  {
    // first compute a bounding box:
    glm::dvec2 p_min { std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity() };
    glm::dvec2 p_max { -std::numeric_limits<double>::infinity(), -std::numeric_limits<double>::infinity() };

    for (const auto& points : branch_trajs.getPoints())
    {
      for (const auto& p : points)
      {
        for (int dim = 0; dim < 2; dim++)
        {
          if (p[dim] < p_min[dim])
          {
            p_min[dim] = p[dim];
          }
          if (p[dim] > p_max[dim])
          {
            p_max[dim] = p[dim];
          }
        }
      }
    }

    // We will need dummy points such that no voronoi vertices can slip outside
    double range = std::max(p_max[0] - p_min[0], p_max[1] - p_min[1]);
    double dist_from_bb = std::max(range, 2 * cutoff);

    dummy_boundary = {
      { p_min[0] - 0.75 * dist_from_bb, p_max[1] + 0.75 * dist_from_bb }, // corner_tl
      { p_min[0], p_max[1] + dist_from_bb }, // top_left
      { p_max[0], p_max[1] + dist_from_bb }, // top_right
      { p_max[0] + 0.75 * dist_from_bb, p_max[1] + 0.75 * dist_from_bb }, // corner_tr
      { p_max[0] + dist_from_bb, p_max[1] }, // right_top
      { p_max[0] + dist_from_bb, p_min[1] }, // right_bottom
      { p_max[0] + 0.75 * dist_from_bb, p_min[1] - 0.75 * dist_from_bb }, // corner_br
      { p_max[0], p_min[1] - dist_from_bb }, // bottom_right
      { p_min[0], p_min[1] - dist_from_bb }, // bottom_left
      { p_min[0] - 0.75 * dist_from_bb, p_min[1] - 0.75 * dist_from_bb }, // corner_bl
      { p_min[0] - dist_from_bb, p_min[1] }, // left_bottom
      { p_min[0] - dist_from_bb, p_max[1] } // left_top
    };

    size_t length = branch_trajs.getHeight() + 1;

    for (const auto& p : dummy_boundary)
    {
      std::vector<glm::dvec2> new_spline;
      for (size_t i = 0; i < length; i++)
      {
        new_spline.push_back(p);
      }
      this->branch_trajs.addTrajectory(new_spline);
    }
  }
}

bool KineticDelaunay::isDummyBoundary(size_t v)
{
  if (add_dummy_boundary)
  {
    return v >= branch_trajs.getPoints().size() - 12;
  }
  return false;
}

glm::dvec2 KineticDelaunay::getPointAt(size_t v, double t) const
{
  // get point transformed such that all points in the same component match
  size_t component_id = component_data.component_map[v];
  size_t representative_vertex = component_data.components[component_id].front();
  size_t reference_branch = branch_trajs.getBranchIndex(representative_vertex, std::ceil(t));

  return branch_trajs.evaluateTransformed(v, t, reference_branch);
}

glm::dvec2 KineticDelaunay::getPointAt(double t, size_t v) const { return getPointAt(v, t); }

std::vector<glm::dvec2> KineticDelaunay::getPointsAt(double t) const
{
  size_t vertex_count = graph.getVertexCount();
  std::vector<glm::dvec2> points;
  points.reserve(vertex_count);
  for (size_t v = 0; v < vertex_count; v++)
  {
    points.push_back(getPointAt(v, t));
  }
  return points;
}

glm::dvec3 KineticDelaunay::getPointInObjectSpace(size_t v, double t) const
{
  return branch_trajs.getPointInObjectSpace(v, t);
}

const StrandTree& KineticDelaunay::getStrandTree() const { return branch_trajs; }

void KineticDelaunay::computeComponentData(double t)
{
  auto& graph = getGraph();
  component_data.components = extractConnectedComponents();
  KINDS_DEBUG("Extracted " << component_data.components.size() << " components.");
  component_data.component_map = buildComponentMap(component_data.components, graph.getVertexCount());
  component_data.component_boundaries.resize(component_data.components.size());

  std::vector<bool> he_visited(graph.getHalfEdges().size(), false);

  for (size_t component_index = 0; component_index < component_data.components.size(); component_index++)
  {
    component_data.component_boundaries[component_index]
      = extractComponentBoundaries(component_data.components[component_index], t, he_visited);
  }

  component_data.component_centroids.resize(component_data.components.size());
  for (size_t component_index = 0; component_index < component_data.components.size(); component_index++)
  {
    if (!component_data.component_boundaries[component_index].empty())
    {
      component_data.component_centroids[component_index]
        = polygonCentroid(component_data.component_boundaries[component_index][0]);
    }
    else
    {
      // compute centroid from points in the component
      glm::dvec2 centroid { 0.0, 0.0 };
      for (auto& v : component_data.components[component_index])
      {
        glm::dvec2 p = getPointAt(t, v);
        centroid += p;
      }
      component_data.component_centroids[component_index]
        = centroid / double(component_data.components[component_index].size());
    }
  }

  component_data.component_last_updated.resize(component_data.components.size(), t);
}

const HalfEdgeDelaunayGraph& KineticDelaunay::init()
{
  graph.init(branch_trajs.getPoints());
  sections_advanced = 0; // Reset the section counter

  graph.printDebug();

  face_inside.clear();
  quadrilateral_last_updated.clear();
  face_last_updated.clear();

  crossing_data.init(graph.getFaces().size());

  face_inside.resize(graph.getFaces().size(), false);
  quadrilateral_last_updated.resize(graph.getHalfEdges().size() / 2, 0.0);
  face_last_updated.resize(graph.getFaces().size(), 0.0);

  for (size_t face_index = 0; face_index < graph.getFaces().size(); face_index++)
  {
    const HalfEdgeDelaunayGraph::Triangle& tri = graph.getFaces()[face_index];

    // compute circumradius at t = 0 and check if within cutoff
    auto vertices = graph.adjacentTriangleVertices(tri.half_edges[0]);
    std::vector<glm::dvec2> points;
    bool outer_face = false;
    for (const auto& v : vertices)
    {
      if (v == -1)
      {
        outer_face = true;
        break;
      }
      // assume that only one plane as frame of reference exists
      points.push_back(branch_trajs.evaluate(v, 0.0));
    }

    if (!outer_face)
    {

      // initialize face_inside based on the circumradius at t = 0
      double r = circumradius(points[0], points[1], points[2]);
      KINDS_DEBUG("Circumradius: " << r);
      if (r < cutoff)
      {
        setFaceInside(face_index, true);
      }
    }

    if (on_the_fly_boundary)
    {
      KINDS_DEBUG("Computing containing triangle for Voronoi vertex " << face_index);
      if (outer_face)
      {
        // The voronoi vertices dual to the outer face are always within it, so we just set it to itself, no events
        // necessary.
        crossing_data.setVoronoiVertexTriId(face_index, face_index);
        continue;
      }

      KINDS_DEBUG("Initial face has vertices at positions: " << glm::to_string(points[0]) << ", " << glm::to_string(points[1]) << ", "
                                 << glm::to_string(points[2]));
      KINDS_DEBUG("Vertex IDs: " << vertices[0] << ", " << vertices[1] << ", " << vertices[2]);
      // initialize voronoi_vertex_to_tri_id:
      // We can use the edge functions from Pineda's algorithm to find if the circumcenter is inside and if not which
      // edge must be crossed. First compute the circumcenter:
      glm::dvec2 circumcenter = HalfEdgeDelaunayGraph::circumcenter(points[0], points[1], points[2]);
      bool inside_triangle = false;

      // define a general edge function that takes the two points of the edge and the circumcenter as input. The sign of
      // the edge function determines on which side of the edge the circumcenter lies.
      auto edge_function = [](const glm::dvec2& a, const glm::dvec2& b, const glm::dvec2& c)
      { return -((b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x)); };

      auto edge_function_01 = edge_function(points[0], points[1], circumcenter);
      auto edge_function_12 = edge_function(points[1], points[2], circumcenter);
      auto edge_function_20 = edge_function(points[2], points[0], circumcenter);

      KINDS_DEBUG("Edge functions: " << edge_function_01 << ", " << edge_function_12 << ", " << edge_function_20);

      // note that only one of the edge functions can be negative by construction of the Delaunay triangulation
      glm::dvec2 midpoint;
      size_t next_crossed_edge_id = -1;
      if (edge_function_01 < 0)
      {
        next_crossed_edge_id = tri.half_edges[0];
        midpoint = (points[0] + points[1]) / 2.0;
      }
      else if (edge_function_12 < 0)
      {
        next_crossed_edge_id = tri.half_edges[1];
        midpoint = (points[1] + points[2]) / 2.0;
      }
      else if (edge_function_20 < 0)
      {
        next_crossed_edge_id = tri.half_edges[2];
        midpoint = (points[2] + points[0]) / 2.0;
      }
      else
      {
        inside_triangle = true;
      }

      while (!inside_triangle)
      {
        // TODO: Get next triangle index and check if circumcenter is inside. If not, compute which of the two edges in
        // the adjacent triangles are crossed using half-plane induced by bisector.
        size_t next_face_id = graph.getHalfEdges()[next_crossed_edge_id ^ 1].face;

        const HalfEdgeDelaunayGraph::Triangle& next_tri = graph.getFaces()[next_face_id];
        auto next_tri_half_edges = graph.getTriangleHalfEdgeIndices(next_crossed_edge_id ^ 1);

        // compute circumradius at t = 0 and check if within cutoff
        auto next_vertices = graph.adjacentTriangleVertices(next_crossed_edge_id ^ 1);
        std::vector<glm::dvec2> next_points;
        bool outer_face = false;
        for (int i = 0; i < next_vertices.size(); i++)
        {
          const auto& v = next_vertices[i];
          if (v == -1)
          {
            outer_face = true;
            break;
          }
          // assume that only one plane as frame of reference exists
          next_points.push_back(branch_trajs.evaluate(v, 0.0));
        }

        if (outer_face)
        {
          KINDS_DEBUG("Next face is an outer face, computing edge functions differently");
          // We need to use different edge functions in this case. For each finite vertex, we use the two neighbor
          // vertices that are also part of infinite triangles and choose a line perpendicular to the line through those
          // two vertices that passes through the vertex

          // TODO: This is only true if we first enter an outer face. If we come from another outer face, we must do
          // this differently.
          size_t finite_he_id = next_crossed_edge_id ^ 1;

          // first look at the origin and the incoming half-edge on the convex hull
          {
            size_t prev_he_id = graph.prevOnConvexBoundaryId(finite_he_id);
            size_t a = graph.getHalfEdges()[finite_he_id].origin;
            size_t c = graph.destination(finite_he_id);
            size_t c_prime = graph.getHalfEdges()[prev_he_id].origin;

            // compute points
            glm::dvec2 p_a = branch_trajs.evaluate(a, 0.0);
            glm::dvec2 p_c = branch_trajs.evaluate(c, 0.0);
            glm::dvec2 p_c_prime = branch_trajs.evaluate(c_prime, 0.0);

            KINDS_DEBUG("Considering points for incoming edge: " << glm::to_string(p_a) << ", " << glm::to_string(p_c) << " and " << glm::to_string(p_c_prime)
                                 << " for edge function computation");

            // compute edge in homogeneous coordinates
            glm::dvec3 edge { p_c_prime.x - p_c.x, p_c_prime.y - p_c.y,
              -(p_c_prime.x - p_c.x) * p_a.x - (p_c_prime.y - p_c.y) * p_a.y };

            double edge_function = glm::dot(edge, glm::dvec3(circumcenter, 1.0));
            KINDS_DEBUG("Edge function for incoming edge: " << edge_function);

            if (edge_function < 0)
            {
              // This is technically not the next crossed edge, but it gives us the only finite edge of the next
              // triangle. This means we do not check if we cross it, but that is not possible anyway because we can
              // never re-enter the convex hull after leaving it, so we can just continue.
              // We could optimize in the future by storing both, so we only need to check one edge in the next
              // iteration.
              next_crossed_edge_id = prev_he_id ^ 1;
              continue;
            }
          }
          // secondly look at the destination and the outgoing half-edge on the convex hull
          {
            size_t next_he_id = graph.nextOnConvexBoundaryId(finite_he_id);
            size_t a = graph.destination(finite_he_id);
            size_t c = graph.getHalfEdges()[finite_he_id].origin;
            size_t c_prime = graph.destination(next_he_id);

            // compute points
            glm::dvec2 p_a = branch_trajs.evaluate(a, 0.0);
            glm::dvec2 p_c = branch_trajs.evaluate(c, 0.0);
            glm::dvec2 p_c_prime = branch_trajs.evaluate(c_prime, 0.0);
            KINDS_DEBUG("Considering points for outgoing edge: " << glm::to_string(p_a) << ", " << glm::to_string(p_c) << " and " << glm::to_string(p_c_prime)
                                 << " for edge function computation");
            // compute edge in homogeneous coordinates
            glm::dvec3 edge { p_c_prime.x - p_c.x, p_c_prime.y - p_c.y,
              -(p_c_prime.x - p_c.x) * p_a.x - (p_c_prime.y - p_c.y) * p_a.y };

            double edge_function = glm::dot(edge, glm::dvec3(circumcenter, 1.0));
            KINDS_DEBUG("Edge function for outgoing edge: " << edge_function);

            if (edge_function < 0)
            {
              // Technically not the next crossed edge, see comment above.
              next_crossed_edge_id = next_he_id ^ 1;
              continue;
            }
          }

          inside_triangle = true;
          crossing_data.setVoronoiVertexTriId(face_index, next_face_id);
          break;
        }

        KINDS_DEBUG("Next face is an inner face, computing edge functions normally");
        KINDS_DEBUG("Next triangle vertices: " << next_vertices[0] << ", " << next_vertices[1] << ", " << next_vertices[2]);
        KINDS_DEBUG("Next triangle points: " << glm::to_string(next_points[0]) << ", " << glm::to_string(next_points[1]) << ", " << glm::to_string(next_points[2]));

        // We can skip edge 01 because that is where we came from and thus cannot be crossed again
        edge_function_12 = edge_function(next_points[1], next_points[2], circumcenter);
        // edge function for edge from points[2] to points[0]
        edge_function_20 = edge_function(next_points[2], next_points[0], circumcenter);
        KINDS_DEBUG("Edge functions for next triangle: " << edge_function_12 << ", " << edge_function_20);

        if (edge_function_12 >= 0 && edge_function_20 >= 0)
        {
          inside_triangle = true;
          crossing_data.setVoronoiVertexTriId(face_index, next_face_id);
        }
        else
        {
          if (edge_function_12 < 0 && edge_function_20 < 0)
          {
            // In this case, we need to check which edge intersects the bisector. The easiest way to do this is testing
            // the edge endpoints against the bisector and pick the one where the sign differs.
            auto check_v0 = edge_function(midpoint, circumcenter, next_points[0]);
            auto check_v1 = edge_function(midpoint, circumcenter, next_points[1]);
            auto check_v2 = edge_function(midpoint, circumcenter, next_points[2]);

            if (std::signbit(check_v0) != std::signbit(check_v1))
            {
              next_crossed_edge_id = next_tri_half_edges[0];
            }
            else if (std::signbit(check_v1) != std::signbit(check_v2))
            {
              next_crossed_edge_id = next_tri_half_edges[1];
            }
            else
            {
              throw std::runtime_error(
                "This should not happen, check_v0, check_v1 and check_v2 cannot all have the same sign");
            }
          }
          else if (edge_function_12 < 0)
          {
            next_crossed_edge_id = next_tri_half_edges[1];
          }
          else if (edge_function_20 < 0)
          {
            next_crossed_edge_id = next_tri_half_edges[2];
          }
        }
      }
    }
  }

  // initialize components
  computeComponentData(0.0);

  return graph;
}

const HalfEdgeDelaunayGraph& KineticDelaunay::advanceOneSection(EventHandler& event_handler)
{
  size_t section_count = branch_trajs.getHeight();
  assert(sections_advanced < section_count); // Ensure we do not exceed the number of sections
  KINDS_DEBUG("Advancing to section " << (sections_advanced + 1) << " of " << section_count);

  // update delaunay graph according to the components
  // For now we assume they can never be merged again
  if (component_data.components.size() > prev_component_count)
  {
    graph.update(branch_trajs.getPoints(), sections_advanced, component_data.components);
  }

  precomputeStep(static_cast<double>(sections_advanced));
  handleEvents(event_handler);
  sections_advanced++;

  return graph;
}

const HalfEdgeDelaunayGraph& KineticDelaunay::getGraph() const { return graph; }

size_t KineticDelaunay::getSectionCount() const { return branch_trajs.getHeight(); }

// Computes the Delaunay triangulation of the given splines
void KineticDelaunay::compute(EventHandler& event_handler)
{
  size_t section_count = getSectionCount(); // Assuming all splines have the same number of points

  ProgressBar progress_bar(
    0, branch_trajs.getHeight(), "Computing Kinetic Voronoi Sections", ProgressBar::Display::Absolute);
  for (size_t i = 0; i < section_count; ++i)
  {
    progress_bar.Update(i);

    assert(i == sections_advanced); // Ensure we are advancing one section at a time
    if (i != 0)
      event_handler.betweenSections(i); // Call the event handler for the section
    advanceOneSection(event_handler);
  }
  progress_bar.Finish();
}

std::vector<size_t> KineticDelaunay::extractConnectedComponent(size_t u, std::vector<bool>& visited) const
{
  std::vector<size_t> component;

  // Perform an iterative DFS with edges induced by inside faces

  std::vector<size_t> stack;
  stack.push_back(u);

  while (!stack.empty())
  {
    size_t v = stack.back();
    stack.pop_back();

    if (visited[v])
      continue;

    visited[v] = true;
    component.push_back(v);

    const auto nbrs = graph.inducedNeighbors(v, face_inside);

    // Push neighbors in reverse order, the same order as recursive DFS
    for (auto it = nbrs.rbegin(); it != nbrs.rend(); ++it)
    {
      size_t w = *it;
      if (!visited[w])
        stack.push_back(w);
    }
  }

  return component;
}

const std::vector<glm::dvec2>& KineticDelaunay::getDummyBoundary() const { return dummy_boundary; }

std::vector<std::vector<size_t>> KineticDelaunay::checkForSplit(const std::array<int, 3>& tri_vertices) const
{
  std::vector<std::vector<size_t>> components;
  std::vector<bool> visited(graph.getVertexCount(), false);

  size_t u = tri_vertices[0];

  std::vector<size_t> component;

  std::vector<size_t> queue;
  queue.push_back(u);
  visited[u] = true;

  size_t head = 0;

  while (head < queue.size())
  {
    size_t v = queue[head++];
    component.push_back(v);

    const auto nbrs = graph.inducedNeighbors(v, face_inside);

    for (size_t w : nbrs)
    {
      if (!visited[w])
      {
        visited[w] = true;

        // quit early if we found all triangle vertices
        if (visited[tri_vertices[1]] && visited[tri_vertices[2]])
        {
          return {}; // return empty to indicate no split
        }

        queue.push_back(w);
      }
    }
  }

  components.push_back(component);

  if (!visited[tri_vertices[1]])
  {
    auto component2 = extractConnectedComponent(tri_vertices[1], visited);
    components.push_back(component2);
  }

  if (!visited[tri_vertices[2]])
  {
    auto component3 = extractConnectedComponent(tri_vertices[2], visited);
    components.push_back(component3);
  }

  return components;
}

std::vector<std::vector<size_t>> KineticDelaunay::extractConnectedComponents() const
{
  std::vector<std::vector<size_t>> components;
  std::vector<bool> visited(graph.getVertexCount(), false);
  for (size_t u = 0; u < graph.getVertexCount(); u++)
  {
    if (visited[u])
    {
      continue;
    }

    auto component = extractConnectedComponent(u, visited);
    components.push_back(component);
  }

  return components;
}

std::vector<BoundaryPoint> KineticDelaunay::traverseBoundary(size_t start_he_id, double t) const
{
  // Walk the boundary to extract the boundary half-edges
  std::vector<BoundaryPoint> boundary_points;
  size_t he_id = start_he_id;
  do
  {
    size_t origin = graph.getHalfEdges()[he_id].origin;
    if (origin == -1)
    {
      KINDS_ERROR("Followed infinite edge.");
    }
    glm::dvec2 pos = getPointAt(origin, t);
    boundary_points.emplace_back(BoundaryPoint { origin, he_id, pos });
    he_id = nextOnComponentBoundaryId(he_id);
  } while (he_id != start_he_id);

  return boundary_points;
}

std::vector<std::vector<BoundaryPoint>> KineticDelaunay::extractComponentBoundaries(
  const std::vector<size_t>& component, double t, std::vector<bool>& he_visited) const
{
  KINDS_DEBUG("Extracting component boundaries at t = " << t);
  if (component.size() < 3)
  {
    return { {} };
  }

  std::vector<std::vector<BoundaryPoint>> boundaries;
  double min_x = std::numeric_limits<double>::infinity();
  // TODO: this is not perfectly safe if points of the outer and an inner boundary coincide at the minimum
  size_t min_x_id = 0;
  for (size_t i = 0; i < component.size(); i++)
  {
    const size_t& v = component[i];

    for (auto it = graph.incidentEdgesBegin(v); it != graph.incidentEdgesEnd(v); it++)
    {
      auto he_id = *it;

      if (he_visited[he_id] || !isOnComponentBoundaryOutside(he_id))
      {
        continue;
      }

      if (graph.destination(he_id) == -1)
      {
        KINDS_ERROR("Destination of half-edge is invalid");
      }
      if (graph.getHalfEdges()[he_id].origin == -1)
      {
        KINDS_ERROR("Origin of half-edge is invalid");
      }

      auto boundary_points = traverseBoundary(he_id, t);

      for (auto& bp : boundary_points)
      {
        he_visited[bp.he_id] = true;
        if (bp.p[0] < min_x)
        {
          min_x = bp.p[0];
          min_x_id = boundaries.size();
        }
      }

      boundaries.emplace_back(boundary_points);
    }
  }

  // swap the boundary with the minimum x to the front
  if (min_x_id != 0)
  {
    std::swap(boundaries[0], boundaries[min_x_id]);
  }

  return boundaries;
}

std::vector<BoundaryPoint> KineticDelaunay::extractComponentBoundary(
  const std::vector<size_t>& component, double t) const
{
  // Find an extreme point to start the boundary walk as it must be on the boundary
  // Note that merely being on the outside of the boundary is not sufficent as there can also be holes inside the
  // component

  size_t start_vertex_id = -1;
  double min_x = std::numeric_limits<double>::infinity();
  for (size_t i = 0; i < component.size(); i++)
  {
    const size_t& v = component[i];

    // Get position and check if it's the minimum x
    glm::dvec2 pos = getPointAt(v, t); // Evaluate at t=0 for starting point
    if (pos[0] < min_x)
    {
      min_x = pos[0];
      start_vertex_id = v;
    }
  }

  // From the starting vertex, find a half-edge that is on the boundary
  size_t start_he_id = -1;
  for (auto it = graph.incidentEdgesBegin(start_vertex_id); it != graph.incidentEdgesEnd(start_vertex_id); it++)
  {
    if (isOnComponentBoundaryOutside(*it))
    {
      start_he_id = *it;
      break;
    }
  }

  return traverseBoundary(start_he_id, t);
}

bool KineticDelaunay::getFaceInside(size_t face_index) const { return face_inside[face_index]; }

void KineticDelaunay::setFaceInside(size_t face_index, bool value)
{
  if (value)
  {
    auto tri_vertices = graph.adjacentTriangleVertices(graph.getFaces()[face_index].half_edges[0]);

    for (int& v : tri_vertices)
    {
      if (v == -1)
      {
        // cannot set face with infinite vertex to inside
        throw std::runtime_error("Cannot set face " + std::to_string(face_index) + " to inside!");
      }
    }
  }
  face_inside[face_index] = value;
}

bool KineticDelaunay::isOnComponentBoundary(size_t he_id) const
{
  size_t face_id = graph.getHalfEdges()[he_id].face;
  size_t twin_face_id = graph.getHalfEdges()[he_id ^ 1].face;
  return (face_inside[face_id] != face_inside[twin_face_id]);
}

bool KineticDelaunay::isOnComponentBoundaryOutside(size_t he_id) const
{
  size_t face_id = graph.getHalfEdges()[he_id].face;
  size_t twin_face_id = graph.getHalfEdges()[he_id ^ 1].face;
  return (!face_inside[face_id] && face_inside[twin_face_id]);
}

size_t KineticDelaunay::nextOnComponentBoundaryId(size_t he_id) const
{
  size_t next_he_id = graph.getHalfEdges()[he_id].next;

  while (!isOnComponentBoundaryOutside(next_he_id))
  {
    next_he_id = graph.twin(next_he_id);
    next_he_id = graph.getHalfEdges()[next_he_id].next;
  }

  return next_he_id;
}
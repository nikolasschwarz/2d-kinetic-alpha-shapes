#include "KineticDelaunay.hpp"
#include <glm/geometric.hpp>

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

// Numeric counterpart of the angular bisector predicate used for infinite / hull cases.
// Returns a value whose sign indicates on which side of the angular bisector at `a`
// (between rays a->c and a->c_prime) the query point `x` lies.
static double angularBisectorPredicate(
  const glm::dvec2& a, const glm::dvec2& c, const glm::dvec2& c_prime, const glm::dvec2& x)
{
  glm::dvec2 ac = c - a;
  glm::dvec2 ac_prime = c_prime - a;
  glm::dvec2 ax = x - a;

  double ac_sq = glm::dot(ac, ac);
  double ac_prime_sq = glm::dot(ac_prime, ac_prime);

  double term1 = glm::dot(ac, ax) * ac_prime_sq;
  double term2 = glm::dot(ac_prime, ax) * ac_sq;

  return term1 - term2;
}

static Polynomial angularBisectorHelper(Trajectory<3>& voronoi_homogeneous, Trajectory<2>& a, Trajectory<2>& c_i, Trajectory<2>& c_j){
  Trajectory<2> a_scaled = a * voronoi_homogeneous[2];
  Trajectory<2> voronoi_xy{voronoi_homogeneous[0], voronoi_homogeneous[1]};
  return (Trajectory<2>::dot(c_i - a, voronoi_xy - a_scaled) * (c_j - a).squaredNorm());
}
static Polynomial angularBisector(Trajectory<2>& a, Trajectory<2>& c, Trajectory<2>& c_prime, Trajectory<3> voronoi_homogeneous){
  return angularBisectorHelper(voronoi_homogeneous, a, c, c_prime) - angularBisectorHelper(voronoi_homogeneous, a, c_prime, c);
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

  auto event_times = findEvents(event_trigger, fraction);
  for (const auto& event_time : event_times)
  {    glm::dvec2 center {};
  
    for (const auto& traj : trajs)
    {
      center[0] += traj[0](event_time);
      center[1] += traj[1](event_time);
    }
    center[0] /= trajs.size();
    center[1] /= trajs.size();
    //KINDS_DEBUG("Boundary Event at time " << event_time + section << " for half-edge ID " << he_id << " at center position"
    //                                      << glm::to_string(center));

    events.emplace(
      Event(event_time + section, he_id, t, center, Event::RADIUS)); // Store the event with the time and half-edge index
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

  auto event_times = findEvents(event_trigger, fraction);
  for (const auto& event_time : event_times)
  {
    glm::dvec2 center {};

    for (const auto& traj : trajs)
    {
      center[0] += traj[0](event_time);
      center[1] += traj[1](event_time);
    }
    center[0] /= trajs.size();
    center[1] /= trajs.size();

    //KINDS_DEBUG("Event at time " << event_time + section << " for half-edge ID " << he_id << " at center position "
    //                                  << glm::to_string(center));

    events.emplace(
      Event(event_time + section, he_id, t, center, Event::FLIP)); // Store the event with the time and half-edge index
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

glm::dvec3 KineticDelaunay::computeVoronoiVertexClampedInfinity(size_t half_edge_id, double t) const
{
  const auto& graph = getGraph();
  const auto& half_edges = graph.getHalfEdges();
  const auto& he = half_edges[half_edge_id];
  const auto& twin_he = half_edges[half_edge_id ^ 1];

  // Compute the positions of the Voronoi vertices at time t.
  // First get the two adjacent triangles
  std::array<int, 3> triVertices = graph.adjacentTriangleVertices(half_edge_id);

  // now compute the circumcenters if the triangles are not infinite
  glm::dvec2 circumcenter;

  bool infinite = false;

  std::vector<glm::dvec2> points;
  std::vector<size_t> vertex_indices;

  size_t infinite_vertex_index = static_cast<size_t>(-1);

  for (size_t i = 0; i < 3; ++i)
  {
    if (triVertices[i] != -1)
    {
      points.push_back(getPointAt(t, triVertices[i]));
      vertex_indices.push_back(static_cast<size_t>(triVertices[i]));
    }
    else
    {
      infinite_vertex_index = i;
    }
  }

  if (points.size() == 3)
  {
    circumcenter = graph.circumcenter(points[0], points[1], points[2]);
    //KINDS_DEBUG("Computed circumcenter: " << glm::to_string(circumcenter));
    // circumcenter = (points[0] + points[1] + points[2]) / 3.0;
  }
  else
  {
    infinite = true;

    // get the triangle on the opposite side of the non-infinite edge
    size_t finite_he_id = half_edge_id;

    while (half_edges[finite_he_id].origin != -1)
    {
      finite_he_id = half_edges[finite_he_id].next;
    }
    finite_he_id = half_edges[finite_he_id].next;
    size_t inner_twin = graph.twin(finite_he_id);
    size_t opposite_vertex = graph.triangleOppositeVertex(inner_twin);
    glm::dvec2 opposite_point = getPointAt(t, opposite_vertex);

    glm::dvec2 neighboring_circumcenter = graph.circumcenter(points[0], points[1], opposite_point);
    // glm::dvec2 neighboring_circumcenter = (points[0] + points[1] + opposite_point) / 3.0;

    // make sure edge points in the correct direction
    if (triVertices[1] == -1)
    {
      std::swap(points[0], points[1]);
    }

    // For now just take the midpoint of the edge
    // circumcenter = (points[0] + points[1]) * 0.5;

    // move circumcenter far out in the direction perpendicular to the edge
    glm::dvec2 edge_dir = glm::normalize(points[1] - points[0]);
    glm::dvec2 perp_dir = glm::dvec2 { -edge_dir[1], edge_dir[0] };

    circumcenter = neighboring_circumcenter - perp_dir;
    //KINDS_DEBUG("Infinite case; replacement circumcenter: " << glm::to_string(circumcenter));
  }

  // place circumcenters into the mesh
  return glm::dvec3 { circumcenter[0], circumcenter[1], t };
}

size_t KineticDelaunay::getCrossingDataContainingTriId(size_t voronoi_vertex_id) const
{
  return crossing_data.getContainingTriId(voronoi_vertex_id);
}

std::vector<size_t> KineticDelaunay::getCrossingDataVoronoiVerticesInTri(size_t tri_id) const
{
  return crossing_data.getVoronoiVerticesInTri(tri_id);
}

glm::dvec3 KineticDelaunay::getVoronoiVertexHomogeneous(size_t voronoi_vertex_id, double t) const
{
  return computeVoronoiVertexHomogenous(voronoi_vertex_id, t);
}

std::vector<double> KineticDelaunay::findEvents(Polynomial& event_trigger, double min_fraction, bool only_positive_to_negative)
{
  if(event_trigger.degree() == -1){
    // No events possible, return empty vector
    return {};
  }

  event_trigger.trim();
  auto zeros = event_trigger.realRoots();
  std::vector<double> filtered_sorted_zeros;

  for (const auto& root : zeros)
  {
    if (isnan(root))
    {
      continue; // Skip NaN roots
    }

    if (root > min_fraction && root <= 1)
    {
      double event_time = root;
      filtered_sorted_zeros.emplace_back(event_time);
    }
  }

  if(filtered_sorted_zeros.empty()){
    // No valid events found, return empty vector
    return {};
  }

  // Sort events ascending by time
  std::sort(filtered_sorted_zeros.begin(), filtered_sorted_zeros.end());

  // Determine sign changes for each root
  std::vector<double> interval_signs(filtered_sorted_zeros.size() + 1);
  double test_point = (min_fraction + filtered_sorted_zeros[0]) / 2.0; // Start with a test point before the first root
  interval_signs[0] = event_trigger(test_point) > 0 ? 1 : -1;


  for (size_t i = 0; i < filtered_sorted_zeros.size(); ++i)
  {
    test_point = (filtered_sorted_zeros[i] + (i + 1 < filtered_sorted_zeros.size() ? filtered_sorted_zeros[i + 1] : 1.0)) / 2.0;
    interval_signs[i + 1] = event_trigger(test_point) > 0 ? 1 : -1;
  }

  std::vector<double> found_event_times;

  for (size_t i = 0; i < filtered_sorted_zeros.size(); ++i)
  {
    if (interval_signs[i] != interval_signs[i + 1])
    {
      if(only_positive_to_negative && interval_signs[i] < 0)
      {
        KINDS_DEBUG("Sign change from negative to positive at root " << filtered_sorted_zeros[i] << ", skipping event creation due to only_positive_to_negative flag.");
        continue; // Skip if we only want positive to negative sign changes
      }
      else
      {
      // Sign change detected, create an event
      double event_time = filtered_sorted_zeros[i];
      found_event_times.push_back(event_time);
      KINDS_DEBUG("Event found at time " << event_time << " with sign change from " << interval_signs[i] << " to " << interval_signs[i + 1]);
      }
    }
    else
    {
      KINDS_DEBUG("No sign change at root " << filtered_sorted_zeros[i] << ", skipping event creation.");
    }
  }

  return found_event_times;
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
  size_t v_i = graph.getHalfEdges()[dual_triangle.half_edges[0]].origin;
  size_t v_j = graph.getHalfEdges()[dual_triangle.half_edges[1]].origin;
  size_t v_k = graph.getHalfEdges()[dual_triangle.half_edges[2]].origin;

  // If a vertex is infinite, so is the Voronoi vertex and it cannot cross any edge, so we can skip this event.
  if (v_i == -1 || v_j == -1 || v_k == -1)
  {
    return;
  }


  // Check a special case: the containing triangle is infinite and adjacent to the dual triangle. In this case, we need a different predicate
  bool adjacent = false;
  size_t adjacent_edge_index = -1;
  size_t finite_he_id = -1;
  if(graph.isInfinite(containing_triangle.half_edges[0]) || graph.isInfinite(containing_triangle.half_edges[1]) || graph.isInfinite(containing_triangle.half_edges[2])){
    // One edge must be finite, find it
    
    size_t finite_edge_index = -1;
    for(size_t edge_index = 0; edge_index < 3; edge_index++){
      if(!graph.isInfinite(containing_triangle.half_edges[edge_index])){
        finite_he_id = containing_triangle.half_edges[edge_index];
        finite_edge_index = edge_index;
        break;
      }
    }

    // Now iterate over the edges of the dual triangle and check if any of them is the twin of the finite edge.
    for(size_t edge_index = 0; edge_index < 3; edge_index++){
      if(graph.twin(dual_triangle.half_edges[edge_index]) == finite_he_id){
        adjacent = true;
        adjacent_edge_index = edge_index;
        break;
      }
    }
  }

  if(adjacent){

    // re-assign vertices
    v_i = graph.triangleOppositeVertex(dual_triangle.half_edges[adjacent_edge_index]);
    v_j = graph.getHalfEdges()[dual_triangle.half_edges[adjacent_edge_index]].origin;
    v_k = graph.getHalfEdges()[dual_triangle.half_edges[adjacent_edge_index] ^ 1].origin;

    Trajectory<2> traj_i = branch_trajs.getPiecePolynomial(v_i, section);
    Trajectory<2> traj_j = branch_trajs.getPiecePolynomial(v_j, section);
    Trajectory<2> traj_k = branch_trajs.getPiecePolynomial(v_k, section);

    Trajectory<2> vector_ij;
    vector_ij[0] = traj_j[0] - traj_i[0];
    vector_ij[1] = traj_j[1] - traj_i[1];
    Trajectory<2> vector_ik;
    vector_ik[0] = traj_k[0] - traj_i[0];
    vector_ik[1] = traj_k[1] - traj_i[1];

    Polynomial event_trigger = -(vector_ij[0] * vector_ik[0] + vector_ij[1] * vector_ik[1]);
    auto fractional_event_times = findEvents(event_trigger, fraction, true);

    // Only need the first event as any following events will be invalidated by the first crossing event. TODO: The exception is the edge being crossed, but that would make this more complex. We can optimize this later if needed.
    if(!fractional_event_times.empty()){
      double fractional_event_time = fractional_event_times.front();
      double event_time = fractional_event_time + section;
      // Position must be the midpoint of the two vertices
      glm::dvec2 position = glm::vec2((traj_j[0](fractional_event_time) + traj_k[0](fractional_event_time)) / 2.0, (traj_j[1](fractional_event_time) + traj_k[1](fractional_event_time)) / 2.0);

      KINDS_DEBUG("Crossing (right angle) Event at time " << event_time << " for Voronoi vertex ID " << voronoi_vertex_id
                                            << " crossing half-edge ID " << finite_he_id << " at position "
                                            << glm::to_string(position));
      events.emplace(event_time, finite_he_id, t, position, voronoi_vertex_id, Event::CROSSING);
    }

  }
  else
  {
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
    Polynomial event_trigger;
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

        // now compute the determinant of the matrix with bisector_ij, bisector_ik and line_ab as columns
      event_trigger = bisector_ij[0] * bisector_ik[1] * line_ab[2]
      + bisector_ij[1] * bisector_ik[2] * line_ab[0] + bisector_ij[2] * bisector_ik[0] * line_ab[1]
      - bisector_ij[2] * bisector_ik[1] * line_ab[0] - bisector_ij[1] * bisector_ik[0] * line_ab[2]
      - bisector_ij[0] * bisector_ik[2] * line_ab[1];
      }
      else
      {
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
          Trajectory<3> voronoi_homogeneous = Trajectory<3>::cross(bisector_ij, bisector_ik);

          event_trigger = angularBisector(traj_a, traj_c, traj_c_prime, voronoi_homogeneous);

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
          Trajectory<3> voronoi_homogeneous = Trajectory<3>::cross(bisector_ij, bisector_ik);

          event_trigger = angularBisector(traj_b, traj_c, traj_c_prime, voronoi_homogeneous);
        }
      }

      auto fractional_event_times = findEvents(event_trigger, fraction, true);
      if (!fractional_event_times.empty())
      {
        double fractional_event_time = fractional_event_times.front();
        double candidate_event_time = fractional_event_time + section;
        if (candidate_event_time < event_time)
        {
          event_time = candidate_event_time;
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
}

void KineticDelaunay::reassignVoronoiVerticesOnBoundary(size_t he_id, double t){
  // Find which side of he_id is on the outside of the convex hull
  size_t face0 = graph.getHalfEdges()[he_id].face;
  size_t face1 = graph.getHalfEdges()[he_id ^ 1].face;

  // If either face has an infinite vertex, that face is the outside
  bool face0_is_outside = false;
  bool face1_is_outside = false;

  // A face is outside if any of its triangle vertices is -1
  auto isFaceOutside = [&](size_t face_id) -> bool {
    const auto& face = graph.getFaces()[face_id];
    for (size_t i = 0; i < 3; ++i) {
      if (graph.getHalfEdges()[face.half_edges[i]].origin == size_t(-1)) {
        return true;
      }
    }
    return false;
  };

  face0_is_outside = isFaceOutside(face0);
  face1_is_outside = isFaceOutside(face1);

  // Determine the outside face id
  size_t outside_face_id = face0_is_outside ? face0 : face1;

  // For each voronoi vertex in face0, reassign to outside face if not already there
  auto vertices0 = crossing_data.getVoronoiVerticesInTri(face0);
  for (size_t voronoi_vertex : vertices0) {
    if (crossing_data.getContainingTriId(voronoi_vertex) != outside_face_id) {
      crossing_data.moveVertex(voronoi_vertex, outside_face_id, t);
    }
  }

  // For each voronoi vertex in face1, reassign to outside face if not already there
  auto vertices1 = crossing_data.getVoronoiVerticesInTri(face1);
  for (size_t voronoi_vertex : vertices1) {
    if (crossing_data.getContainingTriId(voronoi_vertex) != outside_face_id) {
      crossing_data.moveVertex(voronoi_vertex, outside_face_id, t);
    }
  }

  // Secondly, we need to add all intersections to the new Delaunay edge (he_id/2).
  // To do this, find the "inner" edge (the one whose face is not outside) and get its triangle's two other half-edges.
  // Those two edges will have Voronoi–Delaunay intersections; we copy those intersections over to the new (now boundary) edge, adjusting appropriately.

  auto& boundary_edge_d_intersections = crossing_data.delaunay_edge_intersections[he_id / 2];
  auto it = boundary_edge_d_intersections.begin();
  while (it != boundary_edge_d_intersections.end())
  {
    // Capture next BEFORE erasing, because removeIntersection() erases from this list.
    auto next = std::next(it);
    crossing_data.removeIntersection(*it);
    it = next;
  }

  size_t inner_face_id = face0_is_outside ? face1 : face0;
  size_t inner_he_id = (graph.getHalfEdges()[he_id].face == inner_face_id) ? he_id : (he_id ^ 1);

  // Get the three half-edges for this triangle, and identify the other two that aren't inner_he_id
  const auto& inner_face = graph.getFaces()[inner_face_id];
  std::vector<size_t> triangle_half_edges;
  size_t next_he_id = graph.getHalfEdges()[inner_he_id].next;
  triangle_half_edges.push_back(next_he_id);
  next_he_id = graph.getHalfEdges()[next_he_id].next;
  triangle_half_edges.push_back(next_he_id);

  // Now, for each intersection on these two edges, copy the intersection (with changes) to he_id/2 boundary edge
  //for (size_t tri_he : triangle_half_edges) {
  auto copy_intersections = [&](size_t tri_he, bool backwards){
    auto& d_intersections = crossing_data.delaunay_edge_intersections[tri_he / 2];
    bool even = tri_he % 2 == 0;

    auto process_intersection = [&](CrossingData::VoronoiDelaunayEdgeIntersection intersection){
      intersection.delaunay_edge_id = inner_he_id / 2;

      auto& v_intersections = crossing_data.voronoi_edge_intersections[intersection.voronoi_edge_id];

      crossing_data.edge_intersections.emplace_back(intersection);
      auto intersection_it = std::prev(crossing_data.edge_intersections.end());

      // check if the start or end voronoi vertex lies in the outer triangle to determine on which side to insert
      size_t start_voronoi_vertex_id = graph.getHalfEdges()[intersection.voronoi_edge_id].face;
      size_t end_voronoi_vertex_id = graph.getHalfEdges()[intersection.voronoi_edge_id ^ 1].face;
      size_t start_containing_triangle_id = crossing_data.getContainingTriId(start_voronoi_vertex_id);
      size_t end_containing_triangle_id = crossing_data.getContainingTriId(end_voronoi_vertex_id);
      std::list<CrossingData::EdgeIntersectionRef>::iterator v_ref;
      if(start_containing_triangle_id == outside_face_id){
        v_ref = v_intersections.insert(v_intersections.begin(), intersection_it);
      }
      else if(end_containing_triangle_id == outside_face_id){
        v_ref = v_intersections.insert(v_intersections.end(), intersection_it);
      }
      else{
        throw std::runtime_error("Intersection is not on the boundary");
      }

      auto d_ref = boundary_edge_d_intersections.insert(boundary_edge_d_intersections.end(), intersection_it);
      intersection_it->voronoi_ref = v_ref;
      intersection_it->delaunay_ref = d_ref;
    };
    
    if(even != backwards){
      for (auto iter : d_intersections) {
        process_intersection(*iter);
      }
    }
    else
    {
      for(auto iter = d_intersections.rbegin(); iter != d_intersections.rend(); iter++){
        process_intersection(**iter);
      }
    }
  };

  if(inner_he_id % 2 == 0){
    copy_intersections(triangle_half_edges[1], true);
    copy_intersections(triangle_half_edges[0], true);
  }
  else
  {
    copy_intersections(triangle_half_edges[0], false);
    copy_intersections(triangle_half_edges[1], false);
  }

}

void KineticDelaunay::reassignVoronoiVerticesInQuadrilateral(size_t quad_index, double t,
  const std::map<size_t, size_t>& pre_flip_quad_faces)
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
  auto reassign_vertices = [&](const std::vector<size_t>& vertices, size_t target_face_id)
  {
    for (size_t voronoi_vertex : vertices)
    {
      // Compute the position of the Voronoi vertex at time t and check on which side of the edge it is.
      glm::dvec3 voronoi_pos = computeVoronoiVertexHomogenous(voronoi_vertex, t);
      if (voronoi_pos.z != 0)
      {
        glm::dvec2 event_pos = glm::dvec2(voronoi_pos.x / voronoi_pos.z, voronoi_pos.y / voronoi_pos.z);
        if (glm::cross(event_pos - pu, edge_vector) > 0)
        {
          crossing_data.moveVertex(voronoi_vertex, target_face_id, t);
        }
      }
      else
      {
        // Must belong to the dual triangle
        crossing_data.moveVertex(voronoi_vertex, voronoi_vertex, t);
      }

      
    }
  };

  reassign_vertices(vertices0, face_id1);
  reassign_vertices(vertices1, face_id0);

  // get updated vertices
  vertices0 = crossing_data.getVoronoiVerticesInTri(face_id0);
  vertices1 = crossing_data.getVoronoiVerticesInTri(face_id1);

  auto update_intersections = [&](const std::vector<size_t>& vertices)
  {
    // Find and handle all voronoi edges that lie entirely within the quadrilateral
    for (size_t voronoi_vertex : vertices)
    {
      size_t containing_tri_id = crossing_data.getContainingTriId(voronoi_vertex);
      // Detect for all Voronoi edges incident to the vertex if they cross the newly flipped edge
      auto adjacent_voronoi_edges = graph.getFaces()[voronoi_vertex].half_edges;
      for (size_t he_id : adjacent_voronoi_edges)
      {
        size_t other_voronoi_vertex = graph.getHalfEdges()[he_id ^ 1].face;

        size_t other_containing_tri_id = crossing_data.getContainingTriId(other_voronoi_vertex);
        if(other_containing_tri_id == face_id0 || other_containing_tri_id == face_id1)
        {
          // fully within the quadrilateral, might have one intersection removed or added with the flipped edge
          auto v_intersections = crossing_data.voronoi_edge_intersections[he_id/2];
          if(containing_tri_id == other_containing_tri_id)
          {
            // no intersection exists, loop through intersections and remove them
            
            for(auto intersection : v_intersections)
            {
              crossing_data.removeIntersection(intersection);
            }
          }
          else
          {
            // there is an intersection, either update the existing intersection or add a new one
            // Compute the parameter along the delaunay edge
            
            if(v_intersections.empty())
            {
              // add a new intersection
              crossing_data.edge_intersections.emplace_back();
              auto intersection = std::prev(crossing_data.edge_intersections.end());
              intersection->delaunay_edge_id = he_id/2;
              intersection->voronoi_edge_id = he_id/2;
              intersection->delaunay_edge_param = 0.0;
            }
            else
            {
              // update the existing intersection
              auto intersection = v_intersections.front();
              intersection->delaunay_edge_param = 0.0;
            }
          }
        }
      }
    }
  };

  update_intersections(vertices0);
  update_intersections(vertices1);

  // Now handle all voronoi edges that are partially outside the quadrilateral
  auto quad_he_ids = graph.getQuadBoundaryHalfEdgeIndices(quad_index);

  for(size_t he_id : quad_he_ids){
    auto& d_edge_intersections = crossing_data.delaunay_edge_intersections[he_id/2];
    

    for(CrossingData::EdgeIntersectionRef& intersection : d_edge_intersections){
      // get the next and previous intersections to check if they match with the flipped edge
      auto& v_edge_intersections = crossing_data.voronoi_edge_intersections[intersection->voronoi_edge_id];
      auto v_ref = intersection->voronoi_ref;
      std::list<CrossingData::EdgeIntersectionRef>::iterator v_next;
      std::list<CrossingData::EdgeIntersectionRef>::iterator v_prev;
      
      // check if next or previous matches the face
      size_t face_inside_old = pre_flip_quad_faces.at(he_id);
      size_t face_inside_new = graph.getHalfEdges()[he_id].face;

      v_next = std::next(v_ref);
      if(v_ref == v_edge_intersections.begin()){
        v_prev = v_edge_intersections.end();
      }
      else{
        v_prev = std::prev(v_ref);
      }

      bool use_prev = false;
      bool use_next = false;
      bool at_end = false;
      bool intersected_before = false;
      std::list<CrossingData::EdgeIntersectionRef>::iterator v_intersection = v_edge_intersections.end();

      if(v_prev != v_edge_intersections.end()){
        size_t d_edge_id = (*v_prev)->delaunay_edge_id;
        size_t prev_face0 = graph.getHalfEdges()[2 * d_edge_id].face;
        size_t prev_face1 = graph.getHalfEdges()[2 * d_edge_id + 1].face;

        use_prev = (face_inside_old == prev_face0) || (face_inside_old == prev_face1);

        if(use_prev){
          // check if it intersected before
          if((prev_face0 == face_id0 && prev_face1 == face_id1) || (prev_face0 == face_id1 && prev_face1 == face_id0))
          {
            intersected_before = true;
            v_intersection = v_prev;

            if(v_ref == v_edge_intersections.begin()){
              v_prev = v_edge_intersections.end();
            }
            else{
              v_prev = std::next(v_prev);
            }
          }
        }
      } else {
        size_t start_voronoi_vertex_id = graph.getHalfEdges()[2 * intersection->voronoi_edge_id].face;
        size_t containing_triangle_id = crossing_data.getContainingTriId(start_voronoi_vertex_id);

        // unlike the intersection data, this has already been updated, so we need to compare to the new face
        use_prev = containing_triangle_id == face_inside_new;
        if(use_prev){
          at_end = true;
        }
      }

      // now do the same for next
      if(!use_prev){
        if(v_next != v_edge_intersections.end()){
          size_t d_edge_id = (*v_next)->delaunay_edge_id;
          size_t next_face0 = graph.getHalfEdges()[2 * d_edge_id].face;
          size_t next_face1 = graph.getHalfEdges()[2 * d_edge_id + 1].face;

          use_next = (face_inside_new == next_face0) || (face_inside_new == next_face1);
          if(use_next){
            // check if it intersected before
            if((next_face0 == face_id0 && next_face1 == face_id1) || (next_face0 == face_id1 && next_face1 == face_id0))
            {
              intersected_before = true;
              v_intersection = v_next;
              v_next = std::next(v_next);
            }
          }
        } else {
          size_t end_voronoi_vertex_id = graph.getHalfEdges()[2 * intersection->voronoi_edge_id + 1].face;
          size_t containing_triangle_id = crossing_data.getContainingTriId(end_voronoi_vertex_id);

          use_next = containing_triangle_id == face_inside_old;
          if(use_next){
            at_end = true;
          }
        }
      }

      assert(use_next != use_prev);
      if(use_prev){
        v_next = v_prev;
      }

      if(v_next != v_edge_intersections.end()){ // Voronoi edge extends beyond quadrilateral
        auto quad_he_id_it = std::find(quad_he_ids.begin(), quad_he_ids.end(), (*v_next)->delaunay_edge_id * 2);
        if(quad_he_id_it == quad_he_ids.end())
        {
          // try again with odd id
          quad_he_id_it = std::find(quad_he_ids.begin(), quad_he_ids.end(), (*v_next)->delaunay_edge_id * 2 + 1);
        }

        size_t other_he_id = *quad_he_id_it;

        // now compare the face ids to determine if there is an intersection
        size_t face_id = graph.getHalfEdges()[he_id].face;
        size_t other_face_id = graph.getHalfEdges()[other_he_id].face;

        if(face_id != other_face_id){
          // new intersection
          if(!intersected_before){
            crossing_data.edge_intersections.emplace_back();

            CrossingData::EdgeIntersectionRef new_intersection = std::prev(crossing_data.edge_intersections.end());
            new_intersection->delaunay_edge_id = quad_index;
            new_intersection->voronoi_edge_id = intersection->voronoi_edge_id;
            
            auto new_intersection_params = delaunayVoronoiEdgeIntersection(quad_index, intersection->voronoi_edge_id, t);

            new_intersection->delaunay_edge_param = new_intersection_params.first;
            v_edge_intersections.insert(v_next, new_intersection);
            crossing_data.delaunay_edge_intersections[quad_index].push_back(new_intersection); // sort later
            
          }
        } else {
          if(intersected_before){
            crossing_data.removeIntersection(*v_intersection);
          }
        }
      }
      else// Voronoi edge is partially in quadrilateral, check which face its end belongs to.
      { 
        size_t containing_triangle_id;
        if(use_next){
          size_t end_voronoi_vertex_id = graph.getHalfEdges()[2 * intersection->voronoi_edge_id + 1].face;
          containing_triangle_id = crossing_data.getContainingTriId(end_voronoi_vertex_id);
        } else {
          size_t start_voronoi_vertex_id = graph.getHalfEdges()[2 * intersection->voronoi_edge_id].face;
          containing_triangle_id = crossing_data.getContainingTriId(start_voronoi_vertex_id);
        }

        if(containing_triangle_id == face_inside_new && intersected_before){
          // remove intersection
          crossing_data.removeIntersection(*v_intersection);
        }
        else if(containing_triangle_id != face_inside_new && !intersected_before)
        {
          // add new intersection
          crossing_data.edge_intersections.emplace_back();

            CrossingData::EdgeIntersectionRef new_intersection = std::prev(crossing_data.edge_intersections.end());
            new_intersection->delaunay_edge_id = quad_index;
            new_intersection->voronoi_edge_id = intersection->voronoi_edge_id;
            
            auto new_intersection_params = delaunayVoronoiEdgeIntersection(quad_index, intersection->voronoi_edge_id, t);

            new_intersection->delaunay_edge_param = new_intersection_params.first;
            v_edge_intersections.insert(v_next, new_intersection);
            crossing_data.delaunay_edge_intersections[quad_index].push_back(new_intersection); // sort later
        }
      }
    }
  }

  crossing_data.delaunay_edge_intersections[quad_index].sort([&](const CrossingData::EdgeIntersectionRef& a, const CrossingData::EdgeIntersectionRef& b)
      { return a->delaunay_edge_param < b->delaunay_edge_param; });

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

  // Before modifying the topology, store the face id for each half-edge in the quadrilateral
  // (three per triangle) so we can reason about pre-flip topology if needed.
  std::map<size_t, size_t> pre_flip_quad_faces;
  {
    size_t he0 = event.half_edge_id;
    size_t he1 = graph.getHalfEdges()[he0].next;
    size_t he2 = graph.getHalfEdges()[he1].next;
    size_t he3 = he0 ^ 1;
    size_t he4 = graph.getHalfEdges()[he3].next;
    size_t he5 = graph.getHalfEdges()[he4].next;

    pre_flip_quad_faces[he0] = graph.getHalfEdges()[he0].face;
    pre_flip_quad_faces[he1] = graph.getHalfEdges()[he1].face;
    pre_flip_quad_faces[he2] = graph.getHalfEdges()[he2].face;
    pre_flip_quad_faces[he3] = graph.getHalfEdges()[he3].face;
    pre_flip_quad_faces[he4] = graph.getHalfEdges()[he4].face;
    pre_flip_quad_faces[he5] = graph.getHalfEdges()[he5].face;
  }

  // Process the event at the given time
  size_t face_id = graph.getHalfEdges()[event.half_edge_id].face;
  size_t twin_face_id = graph.getHalfEdges()[event.half_edge_id ^ 1].face;
  KINDS_DEBUG("Processing flip event at time " << event.time << " for half-edge ID " << event.half_edge_id
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
  if(!graph.isOnConvexBoundary(event.half_edge_id)){
    reassignVoronoiVerticesInQuadrilateral(event.half_edge_id / 2, event.time, pre_flip_quad_faces);
  }
  else
  {
    reassignVoronoiVerticesOnBoundary(event.half_edge_id, event.time);
  }

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
  KINDS_DEBUG("Processing radius event at time " << event.time << " for half-edge ID " << event.half_edge_id);
  kinDS::HalfEdgeDelaunayGraphToSVG::write(
    getPointsAt(event.time), getGraph(), "test_" + std::to_string(event.time) + "_before.svg", 0.1, &face_inside);
  std::cout << "Wrote " << ("test_" + std::to_string(event.time) + "_before.svg") << std::endl;
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

  KINDS_DEBUG("Processing crossing event at time " << event.time << " for Voronoi vertex ID " << event.voronoi_vertex_id
                                               << " crossing half-edge ID " << event.half_edge_id);

  // move to neighboring triangle
  KINDS_DEBUG("Moving Voronoi vertex " << event.voronoi_vertex_id << " from triangle " << containing_tri_id << " to triangle " << graph.getHalfEdges()[event.half_edge_id ^ 1].face);
  crossing_data.moveVertex(event.voronoi_vertex_id, graph.getHalfEdges()[event.half_edge_id ^ 1].face, event.time);

  event_handler.afterCrossingEvent(event);

  // Update Voronoi–Delaunay edge intersections stored in crossing_data in response to this crossing.
  crossing_data.updateAfterCrossingEvent(*this, event);

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
  //KINDS_DEBUG("Extracted " << component_data.components.size() << " components.");
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

bool KineticDelaunay::computeBoundaryOnTheFly() const
{
  return on_the_fly_boundary;
}

glm::dvec3 lineToHomegeneous(const glm::dvec2& p, const glm::dvec2& dir){
  glm::dvec2 normal = glm::dvec2(-dir.y, dir.x);
  return glm::dvec3(normal, -glm::dot(normal, p));
}

std::pair<double, double> segmentIntersectionParameters(const glm::dvec2& p0, const glm::dvec2& p1, const glm::dvec2& q0, const glm::dvec2& q1){
  auto cross2D = [](const glm::dvec2& a, const glm::dvec2& b)
  {
    return a.x * b.y - a.y * b.x;
  };

  const glm::dvec2 p = p0;
  const glm::dvec2 r = p1 - p0;
  const glm::dvec2 q = q0;
  const glm::dvec2 s = q1 - q0;

  const double rxs = glm::cross(r, s);

  // Parallel (or nearly so): treat as no proper intersection
  if (std::abs(rxs) < 1e-12)
  {
    return std::pair<double, double>(std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity());
  }

  const double t = cross2D(q - p, s) / rxs;
  const double u = cross2D(q - p, r) / rxs;

  return std::pair<double, double>(t, u);
}

bool segmentRayIntersection(const glm::dvec2& segment_p0, const glm::dvec2& segment_p1, const glm::dvec2& ray_origin, const glm::dvec2& ray_dir)
{
  auto [t, u] = segmentIntersectionParameters(segment_p0, segment_p1, ray_origin, ray_origin + ray_dir);
  return t >= 0.0 && u >= 0.0 && u <= 1.0;
}

std::pair<glm::dvec2, glm::dvec2> KineticDelaunay::computeAngularBisector(size_t he_id, double t) const
{
  int a = graph.getHalfEdges()[he_id].origin;
  int b = graph.destination(he_id);
  size_t finite_he_id = graph.getHalfEdges()[he_id].next;
  size_t c, c_prime;

  if (a == -1)
  {
    a = b;
    c = graph.destination(finite_he_id);
    assert(c != size_t(-1));
    size_t prev_he_id = graph.prevOnConvexBoundaryId(finite_he_id);
    c_prime = graph.getHalfEdges()[prev_he_id].origin;
    assert(c_prime != size_t(-1));
  }
  else
  {
    finite_he_id = graph.getHalfEdges()[finite_he_id].next;
    c = graph.getHalfEdges()[finite_he_id].origin;
    assert(c != size_t(-1));
    size_t next_he_id = graph.nextOnConvexBoundaryId(finite_he_id);
    c_prime = graph.destination(next_he_id);
    assert(c_prime != size_t(-1));
  }

  glm::dvec2 p_a = branch_trajs.evaluate(a, t);
  glm::dvec2 p_c = branch_trajs.evaluate(c, t);
  glm::dvec2 p_c_prime = branch_trajs.evaluate(c_prime, t);

  glm::dvec2 angular_bisector_direction =
    glm::normalize(p_c_prime - p_a) + glm::normalize(p_c - p_a);

  // Fallback: If the direction is zero, we can just use the normal of the segment.
  if(glm::length(angular_bisector_direction) < 1e-12){
    glm::dvec2 tangent = glm::normalize(p_c_prime - p_c);
    angular_bisector_direction = glm::dvec2(-tangent.y, tangent.x);
    // TODO: I'm not sure if the sign is correct here.
  }

  // Since we are on the convex hull, the direction points inward, so we need to negate it.
  return std::pair(p_a, -angular_bisector_direction);
};

std::pair<double, double> KineticDelaunay::delaunayVoronoiEdgeIntersection(size_t delaunay_edge_id, size_t voronoi_edge_id, double t) const {

  glm::dvec2 start_point = computeVoronoiVertexClampedInfinity(voronoi_edge_id * 2, t);
  glm::dvec2 destination = computeVoronoiVertexClampedInfinity(voronoi_edge_id * 2 + 1, t);

  if(graph.isInfinite(delaunay_edge_id)){
    auto ray = computeAngularBisector(delaunay_edge_id, t);

    return segmentIntersectionParameters(ray.first, ray.first + ray.second, start_point, destination);
  }
  else{
    glm::dvec2 edge_start = branch_trajs.evaluate(graph.getHalfEdges()[delaunay_edge_id].origin, t);
    glm::dvec2 edge_end = branch_trajs.evaluate(graph.destination(delaunay_edge_id), t);
    return segmentIntersectionParameters(edge_start, edge_end, start_point, destination);
  }
}

  std::pair<std::vector<size_t>, std::vector<double>> KineticDelaunay::computeCrossedHalfEdges(
  size_t start_face_id, const glm::dvec2& destination, const glm::dvec2& start_point, double t) const
{
  std::vector<size_t> crossed_half_edge_ids;
  std::vector<double> crossed_half_edge_params;
    auto edge_function = [](const glm::dvec2& a, const glm::dvec2& b, const glm::dvec2& c)
  { return -((b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x)); };

  auto compute_angular_bisector_homogeneous = [&](size_t he_id)
  {
    auto ray = computeAngularBisector(he_id, t);
    return lineToHomegeneous(ray.first, ray.second);
  };

  auto compute_he_graph_edge_function = [&](size_t he_id, const glm::dvec2& query_point){
    size_t origin = graph.getHalfEdges()[he_id].origin;
    size_t dest = graph.destination(he_id);
    //KINDS_DEBUG("Computing edge function for half-edge " << he_id << " with origin " << origin << " and destination " << dest);
    if (origin != static_cast<size_t>(-1) && dest != static_cast<size_t>(-1))
    {
      glm::dvec2 p0 = branch_trajs.evaluate(origin, t);
      glm::dvec2 p1 = branch_trajs.evaluate(dest, t);
      return -((p1.x - p0.x) * (query_point.y - p0.y) - (p1.y - p0.y) * (query_point.x - p0.x));
    }
    else
    {
      glm::dvec3 edge = compute_angular_bisector_homogeneous(he_id);

      if(dest == -1){
        return -glm::dot(edge, glm::dvec3(query_point, 1.0));
      }
      else
      {
        return glm::dot(edge, glm::dvec3(query_point, 1.0));
      }
    }
  };

  int next_crossed_edge_id = -1;
  bool inside_triangle = false;
  size_t next_face_id = start_face_id;
  auto next_vertices = graph.getTriangleVertexIndices(start_face_id);
  auto next_tri_half_edges = graph.getFaces()[start_face_id].half_edges;

  //KINDS_DEBUG("Following line from " << glm::to_string(start_point) << " to " << glm::to_string(destination));

  while (!inside_triangle)
  {
    //KINDS_DEBUG("Next face ID: " << next_face_id);
    // check each edge for an intersection with the line we need to follow.
    inside_triangle = true;


    // first test if we are inside
    for(int edge_index = 0; edge_index < 3; edge_index++){
      size_t he_id = next_tri_half_edges[edge_index];
      double edge_function = compute_he_graph_edge_function(he_id, destination);
      if(edge_function < 0){
        inside_triangle = false;
        break;
      }
    }

    // Now determine direction of walk
    if(inside_triangle)
    {
      break;
    }

    size_t max_s_index = -1;
    double max_s = -1.0;
    double crossed_edge_param;

    for(int edge_index = 0; edge_index < 3; edge_index++){
      if(edge_index == 0 && next_face_id != start_face_id){
        // This is where we came from, so we don't need to check it as it will make us go backwards
        continue;
      }
      size_t he_id = next_tri_half_edges[edge_index];

      if(graph.isInfinite(he_id)){
        auto ray = computeAngularBisector(he_id, t);
        /*if(segmentRayIntersection(start_point, destination, ray.first, ray.second)){
          next_crossed_edge_id = he_id;
          inside_triangle = false;
          break;
        }*/
        auto [r, s] = segmentIntersectionParameters(ray.first, ray.first + ray.second, start_point, destination);
        if(r >= 0.0 && s <= 1.0){
          if(s > max_s){
            max_s = s;
            max_s_index = edge_index;
            crossed_edge_param = r;
          }
        }
      }
      else{
        glm::dvec2 edge_start = branch_trajs.evaluate(graph.getHalfEdges()[he_id].origin, t);
        glm::dvec2 edge_end = branch_trajs.evaluate(graph.destination(he_id), t);
        auto [r, s] = segmentIntersectionParameters(edge_start, edge_end, start_point, destination);
        if(r >= 0.0 && r <= 1.0 && s <= 1.0){
          if(s > max_s){
            max_s = s;
            max_s_index = edge_index;
            crossed_edge_param = r;
          }
        }
      }
    }
  
    next_crossed_edge_id = next_tri_half_edges[max_s_index];
    next_face_id = graph.getHalfEdges()[next_crossed_edge_id ^ 1].face;
    // Record the edge we are about to cross
    crossed_half_edge_ids.push_back(next_crossed_edge_id);
    crossed_half_edge_params.push_back(crossed_edge_param);
    next_vertices = graph.adjacentTriangleVertices(next_crossed_edge_id ^ 1);
    next_tri_half_edges = graph.getTriangleHalfEdgeIndices(next_crossed_edge_id ^ 1);
  }

  return std::make_pair(crossed_half_edge_ids, crossed_half_edge_params);
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
      //KINDS_DEBUG("Circumradius: " << r);
      if (r < cutoff)
      {
        setFaceInside(face_index, true);
      }
    }

    //KINDS_DEBUG("Computing containing triangle for Voronoi vertex " << face_index);
    if (outer_face)
    {
      // The voronoi vertices dual to the outer face are always within it, so we just set it to itself, no events
      // necessary.
      crossing_data.setVoronoiVertexTriId(face_index, face_index);
      continue;
    }

    //KINDS_DEBUG("Initial face has vertices at positions: " << glm::to_string(points[0]) << ", " << glm::to_string(points[1]) << ", "
    //                           << glm::to_string(points[2]));
    //KINDS_DEBUG("Vertex IDs: " << vertices[0] << ", " << vertices[1] << ", " << vertices[2]);
    // initialize voronoi_vertex_to_tri_id:
    // We can use the edge functions from Pineda's algorithm to find if the circumcenter is inside and if not which
    // edge must be crossed. First compute the circumcenter:
    glm::dvec2 circumcenter = HalfEdgeDelaunayGraph::circumcenter(points[0], points[1], points[2]);

    auto edge_function = [](const glm::dvec2& a, const glm::dvec2& b, const glm::dvec2& c)
    { return -((b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x)); };

    auto edge_function_01 = edge_function(points[0], points[1], circumcenter);
    auto edge_function_12 = edge_function(points[1], points[2], circumcenter);
    auto edge_function_20 = edge_function(points[2], points[0], circumcenter);

    bool inside_triangle = false;
    glm::dvec2 start_point;

    if (edge_function_01 < 0)
    {
      start_point = (points[0] + points[1]) / 2.0;
    }
    else if (edge_function_12 < 0)
    {
      start_point = (points[1] + points[2]) / 2.0;
    }
    else if (edge_function_20 < 0)
    {
      start_point = (points[2] + points[0]) / 2.0;
    }
    else
    {
      inside_triangle = true;
    }

    if (inside_triangle)
    {
      // Circumcenter lies within the initial triangle
      crossing_data.setVoronoiVertexTriId(face_index, face_index);
    }
    else
    {
      //KINDS_DEBUG("Determining containing triangle for voronoi vertex " << face_index);
      auto crossed_half_edges = computeCrossedHalfEdges(face_index, circumcenter, start_point, 0.0).first;

      if (crossed_half_edges.empty())
      {
        // Fallback: if no edges were reported as crossed, treat as inside start triangle
        crossing_data.setVoronoiVertexTriId(face_index, face_index);
      }
      else
      {
        size_t last_crossed_edge_id = crossed_half_edges.back();
        size_t containing_face_id = graph.getHalfEdges()[last_crossed_edge_id ^ 1].face;
        crossing_data.setVoronoiVertexTriId(face_index, containing_face_id);
      }
    }
  }

  // initialize components
  computeComponentData(0.0);

  // Precompute Voronoi–Delaunay edge intersections at t = 0 and store them in crossing_data.
  crossing_data.computeEdgeIntersections(*this, 0.0);

  return graph;
}

void KineticDelaunay::CrossingData::computeEdgeIntersections(const KineticDelaunay& kd, double t)
{
  // Clear any previous intersections
  for (auto& list : voronoi_edge_intersections)
  {
    list.clear();
  }
  edge_intersections.clear();

  size_t num_edges = kd.getGraph().getHalfEdges().size() / 2;
  voronoi_edge_intersections.resize(num_edges);
  delaunay_edge_intersections.resize(num_edges);

  const auto& graph = kd.getGraph();

  // Fill intersections for each Voronoi edge (dual to a Delaunay edge)
  for (size_t voronoi_edge_id = 0; voronoi_edge_id < num_edges; ++voronoi_edge_id)
  {
    size_t he_id0 = voronoi_edge_id * 2;
    size_t he_id1 = he_id0 + 1;

    glm::dvec3 left_pos = kd.computeVoronoiVertexClampedInfinity(he_id0, t);
    glm::dvec3 right_pos = kd.computeVoronoiVertexClampedInfinity(he_id1, t);

    glm::dvec2 left_vertex(left_pos.x, left_pos.y);
    glm::dvec2 right_vertex(right_pos.x, right_pos.y);

    size_t left_voronoi_vertex_id = graph.getHalfEdges()[he_id0].face;
    size_t left_containing_tri_id = kd.getCrossingDataContainingTriId(left_voronoi_vertex_id);

    auto crossed_half_edges_params
      = kd.computeCrossedHalfEdges(left_containing_tri_id, right_vertex, left_vertex, t);

    for (size_t i = 0; i < crossed_half_edges_params.first.size(); i++)
    {
      size_t delaunay_he_id = crossed_half_edges_params.first[i];
      double param = crossed_half_edges_params.second[i];

      edge_intersections.emplace_back();
      auto edge_itr = std::prev(edge_intersections.end());

      edge_itr->delaunay_edge_id = delaunay_he_id / 2;
      edge_itr->voronoi_edge_id = voronoi_edge_id;

      if (delaunay_he_id % 2 == 0)
      {
        edge_itr->delaunay_edge_param = param;
      }
      else
      {
        edge_itr->delaunay_edge_param = 1.0 - param;
      }

      auto voronoi_ref
        = voronoi_edge_intersections[voronoi_edge_id].emplace(voronoi_edge_intersections[voronoi_edge_id].end(), edge_itr);
      edge_itr->voronoi_ref = voronoi_ref;
    }
  }

  // Populate delaunay_edge_intersections and sort by parameter along the Delaunay edge.
  for (auto edge_itr = edge_intersections.begin(); edge_itr != edge_intersections.end(); ++edge_itr)
  {
    size_t delaunay_edge_id = edge_itr->delaunay_edge_id;
    if (delaunay_edge_id >= delaunay_edge_intersections.size())
    {
      KINDS_ERROR("Delaunay edge id out of bounds: " << delaunay_edge_id << " >= " << delaunay_edge_intersections.size());
      continue;
    }
    delaunay_edge_intersections[delaunay_edge_id].push_back(edge_itr);
  }

  for (auto& edge_list : delaunay_edge_intersections)
  {
    edge_list.sort([&](const EdgeIntersectionRef& a, const EdgeIntersectionRef& b)
      { return a->delaunay_edge_param < b->delaunay_edge_param; });

    for (auto it = edge_list.begin(); it != edge_list.end(); ++it)
    {
      (*it)->delaunay_ref = it;
    }
  }
}

void KineticDelaunay::CrossingData::removeIntersection(EdgeIntersectionRef intersection_ref)
{
  // Remove from Delaunay list if the cached iterator is valid.
  auto& d_list = delaunay_edge_intersections[intersection_ref->delaunay_edge_id];
  if (!d_list.empty())
  {
    d_list.erase(intersection_ref->delaunay_ref);
  }

  // Remove from Voronoi list if the cached iterator is valid.
  auto& v_list = voronoi_edge_intersections[intersection_ref->voronoi_edge_id];
  if (!v_list.empty())
  {
    v_list.erase(intersection_ref->voronoi_ref);
  }

  // Finally remove from the global list.
  edge_intersections.erase(intersection_ref);
}

void KineticDelaunay::CrossingData::updateAfterCrossingEvent(const KineticDelaunay& kd, const Event& e)
{
  auto& graph = kd.getGraph();
  size_t voronoi_vertex_id = e.voronoi_vertex_id;
  size_t crossed_delaunay_edge_id = e.half_edge_id / 2;
  auto& d_intersections = delaunay_edge_intersections[crossed_delaunay_edge_id];

  glm::dvec3 voronoi_vertex_position = glm::dvec3(e.position, e.time);
  auto half_edges = graph.getFaces()[voronoi_vertex_id].half_edges;

  bool erased[3] = { false, false, false };
  std::list<EdgeIntersectionRef>::iterator next_after_deletion;

  // First remove any intersection entries that involved this Voronoi vertex and the crossed Delaunay edge.
  for (size_t i = 0; i < half_edges.size(); i++)
  {
    size_t voronoi_he_id = half_edges[i];
    size_t voronoi_edge_id = voronoi_he_id / 2;
    auto& v_intersections = voronoi_edge_intersections[voronoi_edge_id]; // wrong as of now, first edge not deleted in preceding flip event

    auto is_matching = [&](EdgeIntersectionRef ref)
    {
      return ref->delaunay_edge_id == crossed_delaunay_edge_id && ref->voronoi_edge_id == voronoi_edge_id;
    };

    if (!v_intersections.empty() && is_matching(v_intersections.front()))
    {
      auto main_ref = v_intersections.front();
      next_after_deletion = std::next(main_ref->delaunay_ref);
      d_intersections.erase(main_ref->delaunay_ref);
      v_intersections.erase(main_ref->voronoi_ref);
      edge_intersections.erase(main_ref);
      erased[i] = true;
    }
    else if (!v_intersections.empty() && is_matching(v_intersections.back()))
    {
      auto main_ref = v_intersections.back();
      next_after_deletion = std::next(main_ref->delaunay_ref);
      d_intersections.erase(main_ref->delaunay_ref);
      v_intersections.erase(main_ref->voronoi_ref);
      edge_intersections.erase(main_ref);
      erased[i] = true;
    }
  }

  // Now insert new intersection entries for the updated configuration.
  std::list<EdgeIntersectionRef>::iterator inserted_it;
  bool inserted = false;

  for (size_t i = 0; i < half_edges.size(); i++)
  {
    size_t voronoi_he_id = half_edges[i];
    size_t voronoi_edge_id = voronoi_he_id / 2;
    auto& v_intersections = voronoi_edge_intersections[voronoi_edge_id];

    if (!erased[i])
    {
      VoronoiDelaunayEdgeIntersection new_int;
      new_int.delaunay_edge_id = crossed_delaunay_edge_id;
      new_int.voronoi_edge_id = voronoi_edge_id;
      // The exact parameter along the Delaunay edge can be computed if needed; for now, we keep 0.0 as before.
      new_int.delaunay_edge_param = 0.0;

      // Insert into main list
      edge_intersections.push_back(new_int);
      auto main_iter = std::prev(edge_intersections.end());

      // Insert into Delaunay list with appropriate ordering.
      if (!inserted)
      {
        inserted_it = d_intersections.insert(next_after_deletion, main_iter);
        main_iter->delaunay_ref = inserted_it;
        inserted = true;
      }
      else
      {
        bool insert_after = false;
        size_t cell_index_a = graph.getHalfEdges()[voronoi_he_id].origin;
        size_t cell_index_b = graph.destination(voronoi_he_id);
        if (next_after_deletion != d_intersections.end())
        {
          size_t next_voronoi_edge_id = (*next_after_deletion)->voronoi_edge_id;
          size_t next_cell_index_a = graph.getHalfEdges()[2 * next_voronoi_edge_id].origin;
          size_t next_cell_index_b = graph.destination(2 * next_voronoi_edge_id);

          if (cell_index_a == next_cell_index_a || cell_index_a == next_cell_index_b
            || cell_index_b == next_cell_index_a || cell_index_b == next_cell_index_b)
          {
            insert_after = true;
          }
        }
        else
        {
          size_t vertex_id = graph.destination(2 * crossed_delaunay_edge_id);
          if (cell_index_a == vertex_id || cell_index_b == vertex_id)
          {
            insert_after = true;
          }
        }

        if (insert_after)
        {
          inserted_it = d_intersections.insert(next_after_deletion, main_iter);
          main_iter->delaunay_ref = inserted_it;
        }
        else
        {
          inserted_it = d_intersections.insert(inserted_it, main_iter);
          main_iter->delaunay_ref = inserted_it;
        }
      }

      // Insert into Voronoi list (either at front or back depending on which end this Voronoi vertex corresponds to).
      if (e.voronoi_vertex_id == graph.getHalfEdges()[2 * voronoi_edge_id].face)
      {
        v_intersections.emplace_front(main_iter);
        main_iter->voronoi_ref = v_intersections.begin();
      }
      else
      {
        v_intersections.emplace_back(main_iter);
        main_iter->voronoi_ref = std::prev(v_intersections.end());
      }
    }
  }
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
  //KINDS_DEBUG("Extracting component boundaries at t = " << t);
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
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
      if (glm::cross(event_pos - pu, edge_vector) > 0)
      {
        crossing_data.moveVertex(voronoi_vertex, face_id1, t);
      }
    }
    else
    {
      // Must belong to the dual triangle
      crossing_data.moveVertex(voronoi_vertex, voronoi_vertex, t);
    }
  }

  for (size_t voronoi_vertex : vertices1)
  {
    // Compute the position of the Voronoi vertex at time t and check on which side of the edge it is.
    glm::dvec3 voronoi_pos = computeVoronoiVertexHomogenous(voronoi_vertex, t);
    if(voronoi_pos.z != 0)
    {
      glm::dvec2 event_pos = glm::dvec2(voronoi_pos.x / voronoi_pos.z, voronoi_pos.y / voronoi_pos.z);
      if (glm::cross(event_pos - pu, edge_vector) > 0)
      {
        crossing_data.moveVertex(voronoi_vertex, face_id0, t);
      }
    }
    else
    {
      // Must belong to the dual triangle
      crossing_data.moveVertex(voronoi_vertex, voronoi_vertex, t);
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

  KINDS_DEBUG("Processing crossing event at time " << event.time << " for Voronoi vertex ID " << event.voronoi_vertex_id
                                               << " crossing half-edge ID " << event.half_edge_id);

  // move to neighboring triangle
  KINDS_DEBUG("Moving Voronoi vertex " << event.voronoi_vertex_id << " from triangle " << containing_tri_id << " to triangle " << graph.getHalfEdges()[event.half_edge_id ^ 1].face);
  crossing_data.moveVertex(event.voronoi_vertex_id, graph.getHalfEdges()[event.half_edge_id ^ 1].face, event.time);

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

bool lineRayIntersection(
  const glm::dvec2& line_p0,
  const glm::dvec2& line_p1,
  const glm::dvec2& ray_origin,
  const glm::dvec2& ray_dir)
{
  auto cross2D = [](const glm::dvec2& a, const glm::dvec2& b)
  {
    return a.x * b.y - a.y * b.x;
  };

  const glm::dvec2 p = ray_origin;
  const glm::dvec2 r = ray_dir;
  const glm::dvec2 q = line_p0;
  const glm::dvec2 s = line_p1 - line_p0;

  const double rxs = cross2D(r, s);

  // Parallel (or nearly so): treat as no proper intersection
  if (std::abs(rxs) < 1e-12)
  {
    return false;
  }

  const double t_ray = cross2D(q - p, s) / rxs;

  // Intersection must lie on the ray, i.e. parameter >= 0
  return t_ray >= 0.0;
}

  std::vector<size_t> KineticDelaunay::computeCrossedHalfEdges(
  size_t start_face_id, const glm::dvec2& destination, const glm::dvec2& start_point, double t)
{
  // TODO: When I started implementing this, we didn't have an explicit start point and just relied on the start face ID. Perhaps we can simplify this.
  std::vector<size_t> crossed_half_edge_ids;

  auto edge_function = [](const glm::dvec2& a, const glm::dvec2& b, const glm::dvec2& c)
  { return -((b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x)); };

  auto compute_angular_bisector = [&](size_t he_id)
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

    return std::pair(p_a, angular_bisector_direction);
  };

  auto compute_angular_bisector_homogeneous = [&](size_t he_id)
  {
    auto ray = compute_angular_bisector(he_id);
    return lineToHomegeneous(ray.first, ray.second);
  };

  auto compute_he_graph_edge_function = [&](size_t he_id, const glm::dvec2& query_point){
    size_t origin = graph.getHalfEdges()[he_id].origin;
    size_t dest = graph.destination(he_id);
    if (origin != static_cast<size_t>(-1) && dest != static_cast<size_t>(-1))
    {
      glm::dvec2 p0 = branch_trajs.evaluate(origin, t);
      glm::dvec2 p1 = branch_trajs.evaluate(dest, t);
      return -((p1.x - p0.x) * (query_point.y - p0.y) - (p1.y - p0.y) * (query_point.x - p0.x));
    }
    else
    {
      glm::dvec3 edge = compute_angular_bisector_homogeneous(he_id);
      return glm::dot(edge, glm::dvec3(query_point, 1.0));
    }
  };

  int next_crossed_edge_id = -1;
  bool inside_triangle = false;
  size_t next_face_id = start_face_id;
  auto next_vertices = graph.getTriangleVertexIndices(start_face_id);
  auto next_tri_half_edges = graph.getFaces()[start_face_id].half_edges;

  double edge_functions[3] = {0.0, 0.0, 0.0};
  double& edge_function_01 = edge_functions[0];
  double& edge_function_12 = edge_functions[1];
  double& edge_function_20 = edge_functions[2];

  while (!inside_triangle)
  {
    
    std::vector<glm::dvec2> next_points;
    bool outer_face = false;
    bool infinite_vertex[3] = {false, false, false};
    for (int i = 0; i < next_vertices.size(); i++)
    {
      const auto& v = next_vertices[i];
      if (v == -1)
      {
        infinite_vertex[i] = true;
        next_points.push_back(glm::dvec2(0.0, 0.0));
      }
      else {
        next_points.push_back(branch_trajs.evaluate(v, t));
      }
    }

    edge_function_12 = compute_he_graph_edge_function(next_tri_half_edges[1], destination);
    edge_function_20 = compute_he_graph_edge_function(next_tri_half_edges[2], destination);

    if(next_face_id == start_face_id){
      edge_function_01 = compute_he_graph_edge_function(next_tri_half_edges[0], destination);
    }
    else {
      edge_function_01 = 1.0; // positive dummy
    }

    unsigned int violation_count = 0;
    for( double f : edge_functions){
      if(f < 0){
        violation_count++;
      }
    }

    if (violation_count == 0)
    {
      inside_triangle = true;
    }
    else if (violation_count == 2)
    {

      if(!outer_face){
        double check_v[3];

        for(int i = 0; i < 3; i++){
          check_v[i] = edge_function(start_point, destination, next_points[i]);
        }

        for(int i = 0; i < 3; i++){
          if (edge_functions[i] >= 0) // complement is easier to check
          {
            if (std::signbit(check_v[i]) != std::signbit(check_v[(i + 1) % 3]))
            {
              next_crossed_edge_id = next_tri_half_edges[i];
            }
            else if (std::signbit(check_v[(i + 1) % 3]) != std::signbit(check_v[(i + 2) % 3]))
            {
              next_crossed_edge_id = next_tri_half_edges[(i + 1) % 3];
            }
            else
            {
              throw std::runtime_error(
                "This should not happen, check_v0, check_v1 and check_v2 cannot all have the same sign");
            }
          }
        }
      }
      else
      {
        double check_v[3];

        for(int i = 0; i < 3; i++){
          if(!infinite_vertex[i]){
            check_v[i] = edge_function(start_point, destination, next_points[i]);
          }
        }
        // For the two infinite edges, we need to check if the line defined by start_point and destination intersects the edge.
        for(int i = 0; i < 3; i++)
        {
          if(i == 0 && next_face_id != start_face_id){
            // This is where we came from, so we shouldn't check it as it will make us go backwards
            // 
            continue;
          }
          // get edge
          size_t he_id = next_tri_half_edges[i];
          size_t origin = graph.getHalfEdges()[he_id].origin;
          size_t dest = graph.destination(he_id);
          if(origin != static_cast<size_t>(-1) && dest != static_cast<size_t>(-1)){
            // finite, check against the two vertices just as in the non-outer face case
            glm::dvec2 p0 = branch_trajs.evaluate(origin, t);
            glm::dvec2 p1 = branch_trajs.evaluate(dest, t);
            double check_v0 = edge_function(start_point, destination, p0);
            double check_v1 = edge_function(start_point, destination, p1);
           
            if(std::signbit(check_v0) != std::signbit(check_v1)){
              next_crossed_edge_id = he_id;
              break;
            }
          }
          else {
            // infinite, check against the angular bisector
            glm::dvec2 finite_point;
            int origin_index = -1;
            if(origin != -1){
              finite_point = branch_trajs.evaluate(origin, t);
              origin_index = i;
            }
            else
            {
              finite_point = branch_trajs.evaluate(dest, t);
              origin_index = (i + 1) % 3;
            }

            auto ray = compute_angular_bisector(he_id);
            KINDS_DEBUG("Checking ray intersection with finite point " << glm::to_string(finite_point) << " and ray direction " << glm::to_string(ray.second));
            bool ray_intersects = lineRayIntersection(start_point, destination, finite_point, ray.second);
            if(ray_intersects){
              next_crossed_edge_id = he_id;
              break;
            }
          }
        }
      }
    }
    else if (violation_count == 1)
    {
      // Here we can rely on sidedness only.
      if(edge_function_01 < 0)
      {
        next_crossed_edge_id = next_tri_half_edges[0];
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
    else
    {
       throw std::runtime_error(
            "This should not happen as it means that the edges do not form a triangle.");
    }

    if(!inside_triangle){
      next_face_id = graph.getHalfEdges()[next_crossed_edge_id ^ 1].face;
      // Record the edge we are about to cross
      crossed_half_edge_ids.push_back(next_crossed_edge_id);
      next_vertices = graph.adjacentTriangleVertices(next_crossed_edge_id ^ 1);
      next_tri_half_edges = graph.getTriangleHalfEdgeIndices(next_crossed_edge_id ^ 1);
    }
  }

  return crossed_half_edge_ids;
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
      auto crossed_half_edges = computeCrossedHalfEdges(face_index, circumcenter, start_point, 0.0);

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
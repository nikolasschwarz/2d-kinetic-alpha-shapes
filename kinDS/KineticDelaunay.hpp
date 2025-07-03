#pragma once
#include "CubicHermiteSpline.hpp"
#include "HalfEdgeDelaunayGraph.hpp"
#include "Polynomial.hpp"
#include <queue>

namespace kinDS
{

/**
 * \brief Class for computing the Delaunay triangulation of a set of cubic Hermite splines.
 *
 * This follows Guibas, L.J., Mitchell, J.S.B., Roos, T. (1992).
 * Voronoi diagrams of moving points in the plane. In:
 * Schmidt, G., Berghammer, R. (eds) Graph-Theoretic Concepts in Computer Science. WG 1991.
 * Lecture Notes in Computer Science, vol 570. Springer, Berlin, Heidelberg.
 * https://doi.org/10.1007/3-540-55121-2_11
 */
class KineticDelaunay
{
private:
    class Event
    {
    public:
        double time; // Time of the event
        size_t half_edge_id; // Half-edge index associated with the event
        double creation_time; // Time when the event was created, used do check validity after a quadrilateral is updated
        Event(double t, size_t he_id, double creation_time = 0.0)
            : time(t)
            , half_edge_id(he_id)
            , creation_time(creation_time)
        {
        }
        bool operator<(const Event& other) const
        {
            return time > other.time; // For priority queue, we want the earliest event first
        }
    };

    typedef std::priority_queue<Event> EventQueue;

    std::vector<CubicHermiteSpline<2>> splines;
    HalfEdgeDelaunayGraph graph;
    EventQueue events;
    size_t sections_advanced = 0; // Counter for the number of sections advanced

    /* Compare to Leonidas Guibas and Jorge Stolfi. 1985. Primitives for the manipulation of general subdivisions and the computation of Voronoi. ACM Trans. Graph. 4, 2 (April 1985), 74–123. https://doi.org/10.1145/282918.282923
     */
    static Polynomial inCircle(const Polynomial& ax, const Polynomial& ay, const Polynomial& bx, const Polynomial& by, const Polynomial& cx, const Polynomial& cy,
        const Polynomial& px, const Polynomial& py)
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

        // log all just computed polynomials with their variable names
        /* logger.log(DEBUG, "dx: %s", dx.to_string().c_str());
        logger.log(DEBUG, "dy: %s", dy.to_string().c_str());
        logger.log(DEBUG, "ex: %s", ex.to_string().c_str());
        logger.log(DEBUG, "ey: %s", ey.to_string().c_str());
        logger.log(DEBUG, "fx: %s", fx.to_string().c_str());
        logger.log(DEBUG, "fy: %s", fy.to_string().c_str());
        logger.log(DEBUG, "ap: %s", ap.to_string().c_str());
        logger.log(DEBUG, "bp: %s", bp.to_string().c_str());
        logger.log(DEBUG, "cp: %s", cp.to_string().c_str());*/

        return (dx * (ey * cp - bp * fy) - dy * (ex * cp - bp * fx) + ap * (ex * fy - ey * fx));
    }

    static Polynomial ccw(const Polynomial& ax, const Polynomial& ay, const Polynomial& bx, const Polynomial& by, const Polynomial& cx, const Polynomial& cy)
    {
        return (ax * by) + (bx * cy) + (cx * ay) - (ay * bx) - (by * cx) - (cy * ax);
    }

    static Point<2> circumcenter(const Point<2>& a, const Point<2>& b, const Point<2>& c)
    {
        // Calculate the circumcenter of the triangle formed by points a, b, c
        double D = 2 * (a[0] * (b[1] - c[1]) + b[0] * (c[1] - a[1]) + c[0] * (a[1] - b[1]));
        double Ux = ((a[0] * a[0] + a[1] * a[1]) * (b[1] - c[1]) + (b[0] * b[0] + b[1] * b[1]) * (c[1] - a[1]) + (c[0] * c[0] + c[1] * c[1]) * (a[1] - b[1])) / D;
        double Uy = ((a[0] * a[0] + a[1] * a[1]) * (c[0] - b[0]) + (b[0] * b[0] + b[1] * b[1]) * (a[0] - c[0]) + (c[0] * c[0] + c[1] * c[1]) * (b[0] - a[0])) / D;
        return { Ux, Uy };
    }

    /*static Trajectory<2> circumcenter(const Trajectory<2>& a, const Trajectory<2>& b, const Trajectory<2>& c)
    {
        // Calculate the circumcenter of the triangle formed by points a, b, c
        Polynomial D = 2 * (a[0] * (b[1] - c[1]) + b[0] * (c[1] - a[1]) + c[0] * (a[1] - b[1]));
        Polynomial Ux = ((a[0] * a[0] + a[1] * a[1]) * (b[1] - c[1]) + (b[0] * b[0] + b[1] * b[1]) * (c[1] - a[1]) + (c[0] * c[0] + c[1] * c[1]) * (a[1] - b[1])) / D;
        Polynomial Uy = ((a[0] * a[0] + a[1] * a[1]) * (c[0] - b[0]) + (b[0] * b[0] + b[1] * b[1]) * (a[0] - c[0]) + (c[0] * c[0] + c[1] * c[1]) * (b[0] - a[0])) / D;
        return { Ux, Uy };
    }*/

    void computeEvents(double t, size_t quad_id)
    {
        const size_t section = static_cast<size_t>(t);
        const float fraction = t - section;

        size_t he_id = quad_id * 2;
        if (graph.is_on_boundary(he_id) || graph.is_outside_boundary(he_id))
        {
            // boundary edges must be treated separately using ccw

            // need to get the inner half-edge so we have access to the triangle
            // in case both are outside, this swap does not matter, so we just let it happen
            if (graph.is_outside_boundary(he_id))
            {
                he_id = he_id ^ 1; // use the twin half-edge if the current one is on the boundary
            }

            // Depending on the half-edge, the infinite vertex could be in different places, so we just collect all and filter it out
            int indices[4];
            indices[0] = graph.get_half_edges()[he_id].origin; // First vertex
            indices[1] = graph.triangle_opposite_vertex(he_id ^ 1); // Second vertex
            indices[2] = graph.get_half_edges()[he_id ^ 1].origin; // Third vertex
            indices[3] = graph.triangle_opposite_vertex(he_id); // Fourth vertex

            std::vector<int> filtered_indices;

            std::copy_if(indices, indices + 4, std::back_inserter(filtered_indices),
                [this](int index)
                { return index != -1; });

            int& a = filtered_indices[0]; // First vertex
            int& b = filtered_indices[1]; // Second vertex
            int& c = filtered_indices[2]; // Third vertex

            // print the triangle vertices:
            // std::cout << "Triangle vertices: " << a << ", " << b << ", " << c << std::endl;

            Polynomial ax = splines[a].getPiecePolynomial(section)[0];
            Polynomial ay = splines[a].getPiecePolynomial(section)[1];
            Polynomial bx = splines[b].getPiecePolynomial(section)[0];
            Polynomial by = splines[b].getPiecePolynomial(section)[1];
            Polynomial cx = splines[c].getPiecePolynomial(section)[0];
            Polynomial cy = splines[c].getPiecePolynomial(section)[1];

            Polynomial ccw_func = ccw(ax, ay, bx, by, cx, cy);

            // print polynomial:
            // ccw_func.print();

            auto zeros = ccw_func.realRoots();

            // print roots:
            for (const auto& root : zeros)
            {
                if (root > fraction && root <= 1)
                { // Check if the root is within the valid range
                    std::cout << "Root found in triangle (" << a << ", " << b << ", " << c << ") at t = " << root + section << std::endl;
                    events.emplace(root, he_id, fraction); // Store the event with the time and half-edge index
                }
            }
        }
        else
        {
            int a = graph.get_half_edges()[he_id].origin; // First vertex
            int b = graph.triangle_opposite_vertex(he_id ^ 1); // Second vertex
            int c = graph.get_half_edges()[he_id ^ 1].origin; // Third vertex
            int d = graph.triangle_opposite_vertex(he_id); // Fourth vertex

            // print the quadrilateral vertices:
            std::cout << "Quadrilateral vertices: " << a << ", " << b << ", " << c << ", " << d << std::endl;

            Polynomial ax = splines[a].getPiecePolynomial(section)[0];
            Polynomial ay = splines[a].getPiecePolynomial(section)[1];
            Polynomial bx = splines[b].getPiecePolynomial(section)[0];
            Polynomial by = splines[b].getPiecePolynomial(section)[1];
            Polynomial cx = splines[c].getPiecePolynomial(section)[0];
            Polynomial cy = splines[c].getPiecePolynomial(section)[1];
            Polynomial dx = splines[d].getPiecePolynomial(section)[0];
            Polynomial dy = splines[d].getPiecePolynomial(section)[1];

            Polynomial in_circle = inCircle(ax, ay, bx, by, cx, cy, dx, dy);
            // print polynomial:
            // in_circle.print();

            auto zeros = in_circle.realRoots();

            // print roots:
            for (const auto& root : zeros)
            {
                if (root > fraction && root <= 1)
                { // Check if the root is within the valid range
                    std::cout << "Root found in quadrilateral (" << a << ", " << b << ", " << c << ", " << d << ") at t = " << root + section << std::endl;
                    events.emplace(root, he_id, fraction); // Store the event with the time and half-edge index
                }
            }
        }
    }

    void precomputeStep(double t)
    {
        // TODO: Make sure it works where no change of sign occurs in the polynomial, i.e., roots that do not lead to a change in the triangulation.
        size_t quad_count = graph.get_half_edges().size() / 2;
        for (size_t i = 0; i < quad_count; i++)
        {
            computeEvents(t, i);
        }
    }

    void handleEvents()
    {

        std::vector<double> quadrilateral_last_updated(graph.get_half_edges().size() / 2, 0.0);
        while (!events.empty())
        {
            auto event = events.top();
            events.pop();

            // Check if the event is still valid
            if (event.creation_time < quadrilateral_last_updated[event.half_edge_id / 2])
            {
                // This event is outdated, skip it
                continue;
            }

            // Process the event at the given time
            std::cout << "Processing event at time: " << event.time << " for half-edge ID: " << event.half_edge_id << std::endl;

            graph.flipEdge(event.half_edge_id);

            // After flipping the edge, we need to recompute the events for all surrounding half-edges
            size_t next1 = graph.get_half_edges()[event.half_edge_id].next;
            size_t next2 = graph.get_half_edges()[next1].next;

            size_t twin_next1 = graph.get_half_edges()[event.half_edge_id ^ 1].next;
            size_t twin_next2 = graph.get_half_edges()[twin_next1].next;

            computeEvents(event.time, next1 / 2);
            quadrilateral_last_updated[next1 / 2] = event.time; // Update the last updated time for the quadrilateral

            computeEvents(event.time, next2 / 2);
            quadrilateral_last_updated[next2 / 2] = event.time; // Update the last updated time for the quadrilateral

            computeEvents(event.time, twin_next1 / 2);
            quadrilateral_last_updated[twin_next1 / 2] = event.time; // Update the last updated time for the quadrilateral

            computeEvents(event.time, twin_next2 / 2);
            quadrilateral_last_updated[twin_next2 / 2] = event.time; // Update the last updated time for the quadrilateral
        }
    }

public:
    KineticDelaunay(const std::vector<CubicHermiteSpline<2>>& splines)
        : splines(splines)
    {
    }

    std::vector<Point<2>> getPointsAt(double t) const
    {
        std::vector<Point<2>> points;
        points.reserve(splines.size());
        for (const auto& spline : splines)
        {
            points.push_back(spline.evaluate(t)); // Get the first point of each spline
        }
        return points;
    }

    const HalfEdgeDelaunayGraph& init()
    {
        graph.init(splines);
        sections_advanced = 0; // Reset the section counter
        graph.print_debug();

        return graph;
    }

    const HalfEdgeDelaunayGraph& advanceOneSection()
    {
        size_t section_count = splines[0].pointCount() - 1;
        assert(sections_advanced < section_count); // Ensure we do not exceed the number of sections

        precomputeStep(static_cast<double>(sections_advanced));
        handleEvents();
        sections_advanced++;

        return graph;
    }

    const HalfEdgeDelaunayGraph& getGraph() const
    {
        return graph;
    }

    size_t getSectionCount() const
    {
        return splines[0].pointCount() - 1; // Assuming all splines have the same number of points
    }

    // Computes the Delaunay triangulation of the given splines
    void compute()
    {

        size_t section_count = getSectionCount(); // Assuming all splines have the same number of points

        for (size_t i = 0; i < section_count; ++i)
        {
            assert(i == sections_advanced); // Ensure we are advancing one section at a time
            advanceOneSection();
        }
    }
    // Other methods to manipulate and query the triangulation can be added here.
};
} // namespace kinDS

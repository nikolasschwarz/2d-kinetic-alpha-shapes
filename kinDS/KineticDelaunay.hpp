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
    std::vector<CubicHermiteSpline<2>> splines;
    HalfEdgeDelaunayGraph graph;

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
        logger.log(DEBUG, "dx: %s", dx.to_string().c_str());
        logger.log(DEBUG, "dy: %s", dy.to_string().c_str());
        logger.log(DEBUG, "ex: %s", ex.to_string().c_str());
        logger.log(DEBUG, "ey: %s", ey.to_string().c_str());
        logger.log(DEBUG, "fx: %s", fx.to_string().c_str());
        logger.log(DEBUG, "fy: %s", fy.to_string().c_str());
        logger.log(DEBUG, "ap: %s", ap.to_string().c_str());
        logger.log(DEBUG, "bp: %s", bp.to_string().c_str());
        logger.log(DEBUG, "cp: %s", cp.to_string().c_str());

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

    class Event
    {
    public:
        double time; // Time of the event
        size_t half_edge_id; // Half-edge index associated with the event
        Event(double t, size_t he_id)
            : time(t)
            , half_edge_id(he_id)
        {
        }
        bool operator<(const Event& other) const
        {
            return time > other.time; // For priority queue, we want the earliest event first
        }
    };

    typedef std::priority_queue<Event> EventQueue;

    EventQueue precomputeStep(double t)
    {
        // TODO: Make sure it works where no change of sign occurs in the polynomial, i.e., roots that do not lead to a change in the triangulation.
        EventQueue events; // Store events as pairs of (time, half-edge index)
        // Obtain all quadrilaterals, only even indices are considered to avoid duplicates from twin edges
        for (size_t i = 0; i < graph.get_half_edges().size(); i += 2)
        {
            size_t he_id = i;
            if (graph.is_on_boundary(i))
            {
                // boundary edges must be treated separately using ccw

                // need to get the inner half-edge so we have access to the triangle
                if (graph.get_half_edges()[i].face == -1)
                {
                    he_id = i ^ 1; // use the twin half-edge if the current one is on the boundary
                }

                size_t a = graph.get_half_edges()[he_id].origin; // First vertex
                size_t b = graph.destination(he_id); // Second vertex
                size_t c = graph.triangle_opposite_vertex(he_id); // Third vertex

                // print the triangle vertices:
                std::cout << "Triangle vertices: " << a << ", " << b << ", " << c << std::endl;

                Polynomial ax = splines[a].getPiecePolynomial(static_cast<size_t>(t))[0];
                Polynomial ay = splines[a].getPiecePolynomial(static_cast<size_t>(t))[1];
                Polynomial bx = splines[b].getPiecePolynomial(static_cast<size_t>(t))[0];
                Polynomial by = splines[b].getPiecePolynomial(static_cast<size_t>(t))[1];
                Polynomial cx = splines[c].getPiecePolynomial(static_cast<size_t>(t))[0];
                Polynomial cy = splines[c].getPiecePolynomial(static_cast<size_t>(t))[1];

                Polynomial ccw_func = ccw(ax, ay, bx, by, cx, cy);

                // print polynomial:
                ccw_func.print();

                auto zeros = ccw_func.realRoots();

                // print roots:
                for (const auto& root : zeros)
                {
                    if (root >= 0 && root < 1)
                    { // Check if the root is within the valid range
                        std::cout << "Root found in triangle (" << a << ", " << b << ", " << c << ") at t = " << root + t << std::endl;
                        events.emplace(root, he_id); // Store the event with the time and half-edge index
                    }
                }
            }
            else
            {
                size_t a = graph.get_half_edges()[he_id].origin; // First vertex
                size_t b = graph.triangle_opposite_vertex(he_id ^ 1); // Second vertex
                size_t c = graph.get_half_edges()[he_id ^ 1].origin; // Third vertex
                size_t d = graph.triangle_opposite_vertex(he_id); // Fourth vertex

                // print the quadrilateral vertices:
                std::cout << "Quadrilateral vertices: " << a << ", " << b << ", " << c << ", " << d << std::endl;

                Polynomial ax = splines[a].getPiecePolynomial(static_cast<size_t>(t))[0];
                Polynomial ay = splines[a].getPiecePolynomial(static_cast<size_t>(t))[1];
                Polynomial bx = splines[b].getPiecePolynomial(static_cast<size_t>(t))[0];
                Polynomial by = splines[b].getPiecePolynomial(static_cast<size_t>(t))[1];
                Polynomial cx = splines[c].getPiecePolynomial(static_cast<size_t>(t))[0];
                Polynomial cy = splines[c].getPiecePolynomial(static_cast<size_t>(t))[1];
                Polynomial dx = splines[d].getPiecePolynomial(static_cast<size_t>(t))[0];
                Polynomial dy = splines[d].getPiecePolynomial(static_cast<size_t>(t))[1];

                Polynomial in_circle = inCircle(ax, ay, bx, by, cx, cy, dx, dy);
                // print polynomial:
                in_circle.print();

                auto zeros = in_circle.realRoots();

                // print roots:
                for (const auto& root : zeros)
                {
                    if (root >= 0 && root < 1)
                    { // Check if the root is within the valid range
                        std::cout << "Root found in quadrilateral (" << a << ", " << b << ", " << c << ", " << d << ") at t = " << root + t << std::endl;
                        events.emplace(root, he_id); // Store the event with the time and half-edge index
                    }
                }
            }
        }

        return events;
    }

    void handleEvents(EventQueue& events)
    {
        while (!events.empty())
        {
            auto event = events.top();
            events.pop();
            // Process the event at the given time
            std::cout << "Processing event at time: " << event.time << " for half-edge ID: " << event.half_edge_id << std::endl;

            // flip if not on boundary
            graph.flipEdge(event.half_edge_id);
        }
    }

public:
    KineticDelaunay(const std::vector<CubicHermiteSpline<2>>& splines)
        : splines(splines)
    {
    }
    // Computes the Delaunay triangulation of the given splines
    void compute()
    {

        graph.init(splines);

        graph.print_debug();

        graph.flipEdge(6); // Example of flipping the first edge

        graph.print_debug(); // Print the graph after flipping the edge

        auto events = precomputeStep(0.0);
    }
    // Other methods to manipulate and query the triangulation can be added here.
};
} // namespace kinDS

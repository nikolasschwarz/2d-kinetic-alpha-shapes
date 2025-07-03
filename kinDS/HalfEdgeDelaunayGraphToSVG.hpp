#pragma once

#include "../simple_svg.hpp"
#include "HalfEdgeDelaunayGraph.hpp"

using namespace kinDS;

class HalfEdgeDelaunayGraphToSVG
{
public:
    HalfEdgeDelaunayGraphToSVG() = delete;
    ~HalfEdgeDelaunayGraphToSVG() = delete;

    /**
     * @brief Converts the half-edge Delaunay graph to an SVG representation.
     *
     * @param graph The half-edge Delaunay graph to convert.
     * @param filename The name of the output SVG file.
     */
    static void write(const std::vector<Point<2>> points, const HalfEdgeDelaunayGraph& graph, const std::string& filename)
    {
        // figure out the bounding box of the points
        double min_x = std::numeric_limits<double>::max();
        double min_y = std::numeric_limits<double>::max();
        double max_x = std::numeric_limits<double>::lowest();
        double max_y = std::numeric_limits<double>::lowest();

        for (const auto& point : points)
        {
            if (point[0] < min_x)
                min_x = point[0];
            if (point[1] < min_y)
                min_y = point[1];
            if (point[0] > max_x)
                max_x = point[0];
            if (point[1] > max_y)
                max_y = point[1];
        }

        // Create an SVG document with the bounding box
        double width = max_x - min_x;
        double height = max_y - min_y;

        svg::Dimensions dimensions(width, height);

        svg::Document doc(filename, svg::Layout(dimensions, svg::Layout::TopLeft, 1.0, svg::Point(-min_x, -min_y)));
        // Draw vertices
        for (auto& point : points)
        {
            doc << svg::Circle(svg::Point(point[0], point[1]), 3, svg::Fill(svg::Color::Red), svg::Stroke(1, svg::Color::Black));
        }
        // Draw edges
        for (size_t he_id = 0; he_id < graph.get_half_edges().size(); he_id += 2)
        {
            const HalfEdgeDelaunayGraph::HalfEdge& he = graph.get_half_edges()[he_id];
            if (he.next != -1 && he.face != -1) // Only draw edges that are not boundary edges
            {
                Point<2> start = points[graph.get_half_edges()[he_id].origin];
                Point<2> end = points[graph.get_half_edges()[he_id ^ 1].origin];
                doc << svg::Line(svg::Point(start[0], start[1]), svg::Point(end[0], end[1]),
                    svg::Stroke(0.5, svg::Color::Black));
            }
        }
        doc.save();
    }

private:
};

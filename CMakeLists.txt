cmake_minimum_required(VERSION 3.14)
project(kinDS-app VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(kinDS-app main.cpp
	"delaunay/Delaunator2D.cpp"
	
	"kinDS/Polynomial.cpp"
	"kinDS/HalfEdgeDelaunayGraph.cpp"
	"kinDS/VoronoiMesh.cpp"
	"kinDS/Point.cpp"
	"kinDS/SegmentBuilder.cpp"
)

include(CTest)
enable_testing()


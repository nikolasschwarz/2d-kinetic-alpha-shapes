cmake_minimum_required(VERSION 3.14)
project(kinDS-app VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_subdirectory(eigen)

set(KINDS_SOURCES
	"kinDS/Delaunator2D.cpp"
	"kinDS/HalfEdgeDelaunayGraph.cpp"
	"kinDS/Polynomial.cpp"
	"kinDS/SegmentBuilder.cpp"
	"kinDS/TreeMesher.cpp"
	"kinDS/VoronoiMesh.cpp"
 "kinDS/PlaneProjector.cpp")

add_executable(kinDS-app main.cpp
	${KINDS_SOURCES})

# glm
find_package(glm REQUIRED)
target_link_libraries(kinDS-app glm::glm)
add_definitions(-D GLM_ENABLE_EXPERIMENTAL)

# CGAL
find_package(CGAL REQUIRED)
target_link_libraries(kinDS-app CGAL::CGAL)

include(CTest)
enable_testing()


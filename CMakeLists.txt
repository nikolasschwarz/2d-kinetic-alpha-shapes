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
	"kinDS/PlaneProjector.cpp"
	"kinDS/MeshIntersection.cpp")

add_executable(kinDS-app main.cpp
	${KINDS_SOURCES})

# glm
find_package(glm REQUIRED)
target_link_libraries(kinDS-app glm::glm)
add_definitions(-D GLM_ENABLE_EXPERIMENTAL)

# CGAL
find_package(CGAL REQUIRED)
if(CGAL_FOUND)
	target_link_libraries(kinDS-app CGAL::CGAL)
	add_compile_definitions(USE_CGAL)
	MESSAGE(STATUS "EcoSysLab Plugin: CGAL found and linked.")
	else()
	MESSAGE(WARNING "EcoSysLab Plugin: CGAL not found. Some features may be disabled.")
endif()

include(CTest)
enable_testing()


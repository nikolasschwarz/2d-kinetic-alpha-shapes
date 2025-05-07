cmake_minimum_required(VERSION 3.14)
project(kinDS-app VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(kinDS-app main.cpp "voronoi/VoronoiDiagram2D.cpp")

find_package(Qhull CONFIG REQUIRED)
target_link_libraries(kinDS-app PRIVATE Qhull::qhull_r)

find_package(Qhull CONFIG REQUIRED)
target_link_libraries(kinDS-app PRIVATE Qhull::qhullcpp)

include(CTest)
enable_testing()


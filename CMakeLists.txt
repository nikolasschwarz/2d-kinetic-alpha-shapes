cmake_minimum_required(VERSION 3.14)
project(kinDS-app VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(BUILD_SHARED_LIBS OFF CACHE BOOL "Build Qhull as static library")

# Add Qhull
add_subdirectory(external/qhull EXCLUDE_FROM_ALL)

add_executable(kinDS-app main.cpp "voronoi/VoronoiDiagram2D.cpp")

# Include Qhull headers
target_include_directories(kinDS-app PRIVATE
    ${PROJECT_SOURCE_DIR}/external/qhull/src/libqhull_r
)

# Link to Qhull's reentrant library
target_link_libraries(kinDS-app PRIVATE qhullstatic_r)

include(CTest)
enable_testing()


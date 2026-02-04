cmake_minimum_required(VERSION 3.14)

# Detect if this is being used as a submodule
# If CMAKE_SOURCE_DIR != CMAKE_CURRENT_SOURCE_DIR, we're a submodule
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
    # This is the main project
    project(kinDS VERSION 0.1.0 LANGUAGES CXX)
    set(KINDS_BUILD_DEMO ON)
else()
    # This is being used as a submodule
    project(kinDS VERSION 0.1.0 LANGUAGES CXX)
    set(KINDS_BUILD_DEMO OFF)
endif()

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Add Eigen as subdirectory only if not already provided by parent project
# Check if Eigen3::Eigen target already exists (from parent project)
if(NOT TARGET Eigen3::Eigen)
    # Eigen not found, add our submodule version
    if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/eigen/CMakeLists.txt)
        add_subdirectory(eigen)
        message(STATUS "kinDS: Using bundled Eigen submodule.")
    else()
        # Try to find Eigen via find_package
        find_package(Eigen3 QUIET)
        if(NOT Eigen3_FOUND)
            message(FATAL_ERROR "kinDS: Eigen not found. Either provide Eigen3::Eigen target in parent project, or ensure the eigen submodule is initialized.")
        else()
            message(STATUS "kinDS: Using Eigen from find_package.")
        endif()
    endif()
else()
    message(STATUS "kinDS: Using Eigen provided by parent project.")
endif()

# ============================================================================
# kinDS Library
# ============================================================================

set(KINDS_SOURCES
	"kinDS/Delaunator2D.cpp"
	"kinDS/HalfEdgeDelaunayGraph.cpp"
	"kinDS/Polynomial.cpp"
	"kinDS/SegmentBuilder.cpp"
	"kinDS/TreeMesher.cpp"
	"kinDS/VoronoiMesh.cpp"
	"kinDS/PlaneProjector.cpp"
	"kinDS/MeshIntersection.cpp")

# Create the library target
add_library(kinDS ${KINDS_SOURCES})

# Set include directories - allows both #include "kinDS/..." and #include "..."
target_include_directories(kinDS
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
        $<INSTALL_INTERFACE:include>
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/kinDS
)

# Find and link dependencies
find_package(glm REQUIRED)
target_link_libraries(kinDS PUBLIC glm::glm)
target_compile_definitions(kinDS PUBLIC GLM_ENABLE_EXPERIMENTAL)

# CGAL (optional)
find_package(CGAL QUIET)
if(CGAL_FOUND)
	target_link_libraries(kinDS PUBLIC CGAL::CGAL)
	target_compile_definitions(kinDS PUBLIC USE_CGAL)
	message(STATUS "kinDS: CGAL found and linked.")
else()
	message(WARNING "kinDS: CGAL not found. Some features may be disabled.")
endif()

# Link Eigen (header-only, but we add it for consistency)
target_link_libraries(kinDS PUBLIC Eigen3::Eigen)

# Export the library for use in parent projects
set_target_properties(kinDS PROPERTIES
    VERSION ${PROJECT_VERSION}
    SOVERSION ${PROJECT_VERSION_MAJOR}
)

# ============================================================================
# Demo Executable (only build if this is the main project)
# ============================================================================

if(KINDS_BUILD_DEMO)
    add_executable(kinDS-demo main.cpp)
    target_link_libraries(kinDS-demo PRIVATE kinDS)
    
    # Set output directory for demo
    set_target_properties(kinDS-demo PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin
    )
endif()

# ============================================================================
# Testing
# ============================================================================

include(CTest)
enable_testing()


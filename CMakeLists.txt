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

# Report build type for visibility (important when used as submodule)
if(CMAKE_BUILD_TYPE)
    message(STATUS "kinDS: Build type = ${CMAKE_BUILD_TYPE}")
else()
    if(CMAKE_CONFIGURATION_TYPES)
        message(STATUS "kinDS: Multi-config generator (e.g. Visual Studio); build type chosen at build time (--config). Available: ${CMAKE_CONFIGURATION_TYPES}")
    else()
        message(STATUS "kinDS: Build type not set (CMAKE_BUILD_TYPE empty); default may apply depending on generator")
    endif()
endif()

if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
    message(STATUS "kinDS: Standalone = TRUE")
else()
    message(STATUS "kinDS: Standalone = FALSE")
endif()

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
	"kinDS/StrandTree.cpp"
	"kinDS/Delaunator2D.cpp"
	"kinDS/HalfEdgeDelaunayGraph.cpp"
    "kinDS/KineticDelaunay.cpp"
	"kinDS/Polynomial.cpp"
	"kinDS/SegmentBuilder.cpp"
	"kinDS/TreeMesher.cpp"
	"kinDS/VoronoiMesh.cpp"
	"kinDS/PlaneProjector.cpp"
	"kinDS/MeshIntersection.cpp")

# Create the library target
add_library(kinDS ${KINDS_SOURCES})

# Set include directories
# When used as submodule: allows #include "kinDS/Header.hpp" to work correctly
# When used standalone: allows #include "kinDS/Header.hpp" to work correctly
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
    # Standalone build - only need project root
    set(KINDS_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR})
    
    # For standalone builds: create directory structure to support "eigen/eigen/Eigen/..." includes
    # Create a directory junction/symlink in the build directory that maps eigen/eigen -> eigen
    # This allows the same include paths to work in both standalone and parent project builds
    if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/eigen/CMakeLists.txt)
        set(EIGEN_EIGEN_LINK_DIR ${CMAKE_BINARY_DIR}/eigen/eigen)
        set(EIGEN_TARGET_DIR ${CMAKE_CURRENT_SOURCE_DIR}/eigen)
        
        # Create parent directory
        file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/eigen)
        
        # Create symlink/directory junction (platform-specific)
        set(LINK_CREATED FALSE)
        if(WIN32)
            # On Windows, try to create a directory junction
            # Note: This may require admin privileges or developer mode
            # First, remove the directory if it exists (as a regular directory)
            if(EXISTS ${EIGEN_EIGEN_LINK_DIR})
                file(REMOVE_RECURSE ${EIGEN_EIGEN_LINK_DIR})
            endif()
            execute_process(
                COMMAND cmd /c mklink /J "${EIGEN_EIGEN_LINK_DIR}" "${EIGEN_TARGET_DIR}"
                ERROR_QUIET
                OUTPUT_QUIET
                RESULT_VARIABLE JUNCTION_RESULT
            )
            if(JUNCTION_RESULT EQUAL 0)
                set(LINK_CREATED TRUE)
                message(STATUS "kinDS: Created directory junction for eigen/eigen -> eigen")
            else()
                message(STATUS "kinDS: Could not create directory junction (may require admin/dev mode). Using alternative include path setup.")
            endif()
        else()
            # On Unix, create a symlink
            if(EXISTS ${EIGEN_EIGEN_LINK_DIR})
                file(REMOVE ${EIGEN_EIGEN_LINK_DIR})
            endif()
            execute_process(
                COMMAND ${CMAKE_COMMAND} -E create_symlink "${EIGEN_TARGET_DIR}" "${EIGEN_EIGEN_LINK_DIR}"
                ERROR_QUIET
                OUTPUT_QUIET
                RESULT_VARIABLE SYMLINK_RESULT
            )
            if(SYMLINK_RESULT EQUAL 0)
                set(LINK_CREATED TRUE)
                message(STATUS "kinDS: Created symlink for eigen/eigen -> eigen")
            else()
                message(STATUS "kinDS: Could not create symlink. Using alternative include path setup.")
            endif()
        endif()
        
        # Add build directory to include paths so eigen/eigen can be resolved
        # This allows "eigen/eigen/Eigen/..." to work in standalone builds
        list(APPEND KINDS_INCLUDE_DIRS ${CMAKE_BINARY_DIR})
        
        # If link creation failed at configure time, create a custom command to create it at build time
        if(NOT LINK_CREATED)
            if(WIN32)
                # Use a PowerShell script file to properly handle paths with spaces
                set(PS_SCRIPT ${CMAKE_BINARY_DIR}/create_eigen_link.ps1)
                file(TO_NATIVE_PATH "${EIGEN_EIGEN_LINK_DIR}" EIGEN_EIGEN_LINK_DIR_NATIVE)
                file(TO_NATIVE_PATH "${EIGEN_TARGET_DIR}" EIGEN_TARGET_DIR_NATIVE)
                file(WRITE ${PS_SCRIPT}
                    "if (Test-Path '${EIGEN_EIGEN_LINK_DIR_NATIVE}') {\n"
                    "    Remove-Item -Path '${EIGEN_EIGEN_LINK_DIR_NATIVE}' -Force -Recurse\n"
                    "}\n"
                    "New-Item -ItemType Junction -Path '${EIGEN_EIGEN_LINK_DIR_NATIVE}' -Target '${EIGEN_TARGET_DIR_NATIVE}'\n"
                )
                add_custom_command(
                    OUTPUT ${EIGEN_EIGEN_LINK_DIR}
                    COMMAND PowerShell -ExecutionPolicy Bypass -File "${PS_SCRIPT}"
                    COMMENT "Creating directory junction for eigen/eigen -> eigen"
                    VERBATIM
                )
            else()
                add_custom_command(
                    OUTPUT ${EIGEN_EIGEN_LINK_DIR}
                    COMMAND ${CMAKE_COMMAND} -E remove ${EIGEN_EIGEN_LINK_DIR}
                    COMMAND ${CMAKE_COMMAND} -E create_symlink "${EIGEN_TARGET_DIR}" "${EIGEN_EIGEN_LINK_DIR}"
                    COMMENT "Creating symlink for eigen/eigen -> eigen"
                    VERBATIM
                )
            endif()
            # Make the link a dependency of the kinDS target so it's created before compilation
            add_custom_target(create_eigen_link DEPENDS ${EIGEN_EIGEN_LINK_DIR})
            add_dependencies(kinDS create_eigen_link)
        endif()
    endif()
else()
    # Used as submodule - add both parent directory and project root
    # This allows #include "kinDS/Header.hpp" to work regardless of how
    # the parent project sets up include paths
    get_filename_component(KINDS_PARENT_DIR ${CMAKE_CURRENT_SOURCE_DIR} DIRECTORY)
    set(KINDS_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR} ${KINDS_PARENT_DIR})
endif()

target_include_directories(kinDS
    PUBLIC
        $<BUILD_INTERFACE:${KINDS_INCLUDE_DIRS}>
        $<INSTALL_INTERFACE:include>
    PRIVATE
        # For internal includes within the library
        ${CMAKE_CURRENT_SOURCE_DIR}/kinDS
)

# Find and link dependencies
# glm: allow parent project to provide glm::glm target if it already uses glm
if(NOT TARGET glm::glm)
    find_package(glm REQUIRED)
endif()
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

# Only build tests if this is the main project (not a submodule)
if(KINDS_BUILD_DEMO)
    # Use FetchContent to get Catch2 (header-only, no installation needed)
    include(FetchContent)
    FetchContent_Declare(
        Catch2
        GIT_REPOSITORY https://github.com/catchorg/Catch2.git
        GIT_TAG        v3.5.0
    )
    FetchContent_MakeAvailable(Catch2)
    
    # Create test executable
    file(GLOB TEST_SOURCES "tests/*.cpp")
    if(TEST_SOURCES)
        add_executable(kinDS-tests ${TEST_SOURCES})
        target_link_libraries(kinDS-tests PRIVATE kinDS Catch2::Catch2WithMain)
        target_include_directories(kinDS-tests PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
        
        # Register with CTest using Catch2's test discovery
        # Catch2 v3 provides the Catch.cmake module
        list(APPEND CMAKE_MODULE_PATH ${Catch2_SOURCE_DIR}/extras)
        include(Catch)
        
        # Try to use catch_discover_tests, fallback to manual registration if it fails
        if(COMMAND catch_discover_tests)
            catch_discover_tests(kinDS-tests)
        else()
            # Fallback: manually register the test executable
            add_test(NAME kinDS-tests COMMAND kinDS-tests)
            message(STATUS "kinDS: Using manual test registration (catch_discover_tests not available)")
        endif()
        
        message(STATUS "kinDS: Tests enabled. Run with: ctest (from build directory)")
    else()
        message(STATUS "kinDS: No test files found in tests/ directory")
    endif()
endif()


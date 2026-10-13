# Chained toolchain for vcpkg + Ninja on Windows with MSVC.
#
# vcpkg's CGAL 5.6.x does not compile with Clang (invalid operator bool in
# CGAL/boost/graph/iterator.h). Ninja presets require MSVC (cl.exe).

set(_kinds_vcpkg_root "")
if(DEFINED VCPKG_ROOT AND NOT VCPKG_ROOT STREQUAL "")
  set(_kinds_vcpkg_root "${VCPKG_ROOT}")
elseif(EXISTS "C:/src/vcpkg/installed/x64-windows/share/glm/glmConfig.cmake")
  set(_kinds_vcpkg_root "C:/src/vcpkg")
elseif(DEFINED ENV{VCPKG_ROOT} AND EXISTS "$ENV{VCPKG_ROOT}/installed/x64-windows/share/glm/glmConfig.cmake")
  set(_kinds_vcpkg_root "$ENV{VCPKG_ROOT}")
elseif(DEFINED ENV{VCPKG_ROOT})
  set(_kinds_vcpkg_root "$ENV{VCPKG_ROOT}")
endif()

if(_kinds_vcpkg_root STREQUAL "")
  message(FATAL_ERROR
    "kinDS: Could not locate a vcpkg installation with glm for x64-windows.\n"
    "Set VCPKG_ROOT in CMakePresets.json / CMakeUserPresets.json, or install glm:\n"
    "  vcpkg install glm:x64-windows cgal:x64-windows")
endif()

set(VCPKG_ROOT "${_kinds_vcpkg_root}" CACHE PATH "vcpkg root" FORCE)
set(_vcpkg_toolchain "${_kinds_vcpkg_root}/scripts/buildsystems/vcpkg.cmake")
if(NOT EXISTS "${_vcpkg_toolchain}")
  message(FATAL_ERROR "kinDS: vcpkg toolchain not found at ${_vcpkg_toolchain}")
endif()

message(STATUS "kinDS: Using vcpkg from ${VCPKG_ROOT}")
include("${_vcpkg_toolchain}")

if(NOT WIN32)
  return()
endif()

find_program(_kinds_msvc_cl NAMES cl cl.exe)
find_program(_kinds_clang NAMES clang++ clang-cl clang++)

if(_kinds_msvc_cl)
  set(CMAKE_C_COMPILER "${_kinds_msvc_cl}" CACHE FILEPATH "MSVC C compiler" FORCE)
  set(CMAKE_CXX_COMPILER "${_kinds_msvc_cl}" CACHE FILEPATH "MSVC C++ compiler" FORCE)
elseif(_kinds_clang)
  message(FATAL_ERROR
    "kinDS: Ninja builds on Windows require MSVC (cl.exe), not Clang.\n"
    "vcpkg CGAL headers fail to compile with Clang (CGAL/boost/graph/iterator.h).\n"
    "Open 'Developer PowerShell for VS 2022' (or run vcvars64.bat), then reconfigure.\n"
    "Alternatively use the Visual Studio generator preset: cmake --preset vcpkg-x64-vs")
else()
  message(FATAL_ERROR
    "kinDS: MSVC (cl.exe) not found on PATH.\n"
    "Open 'Developer PowerShell for VS 2022' or install the VS 2022 C++ workload.")
endif()

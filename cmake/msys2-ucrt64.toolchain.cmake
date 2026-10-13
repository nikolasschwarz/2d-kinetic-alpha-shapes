# MSYS2 UCRT64 GCC toolchain for headless debug runs (no MSVC assert dialogs).
# Uses glm from MSYS2 pacman; Eigen from the bundled submodule; CGAL is optional.

set(_msys2_ucrt64 "C:/msys64/ucrt64")

if(NOT EXISTS "${_msys2_ucrt64}/bin/g++.exe")
  message(FATAL_ERROR
    "kinDS: MSYS2 UCRT64 GCC not found at ${_msys2_ucrt64}/bin/g++.exe\n"
    "Install: pacman -S mingw-w64-ucrt-x86_64-gcc mingw-w64-ucrt-x86_64-cmake "
    "mingw-w64-ucrt-x86_64-ninja mingw-w64-ucrt-x86_64-glm")
endif()

set(CMAKE_C_COMPILER "${_msys2_ucrt64}/bin/gcc.exe" CACHE FILEPATH "MSYS2 UCRT64 GCC" FORCE)
set(CMAKE_CXX_COMPILER "${_msys2_ucrt64}/bin/g++.exe" CACHE FILEPATH "MSYS2 UCRT64 G++" FORCE)

list(PREPEND CMAKE_PREFIX_PATH "${_msys2_ucrt64}")
list(PREPEND CMAKE_PROGRAM_PATH "${_msys2_ucrt64}/bin")

message(STATUS "kinDS: MSYS2 UCRT64 toolchain (${CMAKE_CXX_COMPILER})")

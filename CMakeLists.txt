cmake_minimum_required(VERSION 3.10.0)
project(kinDS-app VERSION 0.1.0 LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 17)

add_executable(kinDS-app main.cpp)

include(CTest)
enable_testing()


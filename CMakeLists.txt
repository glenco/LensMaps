cmake_minimum_required(VERSION 3.0...3.8)

project(LensMaps)

add_executable(LensMaps LensMaps.cpp)

find_package(GLAMER NO_MODULE REQUIRED)
include_directories(${GLAMER_INCLUDE_DIRS})
target_link_libraries(LensMaps ${GLAMER_LIBRARIES})

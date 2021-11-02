cmake_minimum_required(VERSION 3.14)

project(fit_plane LANGUAGES CXX)

set(CMAKE_INCLUDE_CURRENT_DIR ON)

find_package(Eigen3 CONFIG REQUIRED)
find_package(fmt CONFIG REQUIRED)

set(HEADERS fit_plane.h)
set(SOURCES fit_plane.cpp main.cpp)

add_executable(fit_plane ${HEADERS} ${SOURCES})
target_link_libraries (fit_plane Eigen3::Eigen fmt)

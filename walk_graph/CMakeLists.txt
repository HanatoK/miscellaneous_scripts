cmake_minimum_required(VERSION 3.14)

project(walk_graph LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(TBB REQUIRED)

add_executable(test test.cpp)
target_link_libraries(test PRIVATE TBB::tbb)

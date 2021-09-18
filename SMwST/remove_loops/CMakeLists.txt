cmake_minimum_required(VERSION 3.10)
# set the project name and version
project(RemoveLoops VERSION 1.0)
project(RemoveLoops LANGUAGES CXX)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_INCLUDE_CURRENT_DIR ON)

set(HEADERS Graph.h  Helper.h  Pathway.h)
set(SOURCES Graph.cpp  Helper.cpp  main.cpp  Pathway.cpp)

add_executable(${PROJECT_NAME} ${SOURCES} ${HEADERS})
target_link_libraries(${PROJECT_NAME} fmt)

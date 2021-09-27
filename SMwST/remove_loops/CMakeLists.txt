cmake_minimum_required(VERSION 3.10)
# set the project name and version
project(RemoveLoops VERSION 1.0)
project(RemoveLoops LANGUAGES CXX)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_INCLUDE_CURRENT_DIR ON)

# in 2.6+, stubs are split off into a separate file
IF (${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION} GREATER 2.4)
  find_package(TclStub REQUIRED)
ELSE (${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION} GREATER 2.4)
  find_package(TCL REQUIRED)
ENDIF (${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION} GREATER 2.4)

find_package(fmt REQUIRED)

set(HEADERS Graph.h Helper.h Pathway.h)
set(SOURCES Graph.cpp Helper.cpp Pathway.cpp main.cpp)
set(LIB_HEADERS Matrix.h Spline.h Reparametrization.h Graph.h Helper.h Pathway.h)
set(LIB_SOURCES Matrix.cpp Spline.cpp Reparametrization.cpp Graph.cpp Helper.cpp Pathway.cpp SMwST_functions.cpp)

add_library(smwst SHARED ${LIB_HEADERS} ${LIB_SOURCES})
target_compile_definitions(smwst PRIVATE USE_TCL_STUBS)
target_link_libraries(smwst ${TCL_STUB_LIBRARY} fmt)
include_directories(${TCL_INCLUDE_PATH})

add_executable(${PROJECT_NAME} ${SOURCES} ${HEADERS})
target_link_libraries(${PROJECT_NAME} fmt)

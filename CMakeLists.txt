cmake_minimum_required (VERSION 3.10.2)

project(vspc_graph CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS_DEBUG "-D_DEBUG")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")
set(CMAKE_INSTALL_PREFIX ${CMAKE_SOURCE_DIR})

include_directories("include")

add_executable(main src/main.cpp)

install(TARGETS main DESTINATION bin)

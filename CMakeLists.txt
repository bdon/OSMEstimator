cmake_minimum_required (VERSION 3.5)
set(CMAKE_CXX_FLAGS_RELEASE "-O3")
set(CMAKE_CXX_FLAGS_DEBUG "-DDEBUG -g")
set(CMAKE_CXX_FLAGS "-Wno-deprecated")
set(CMAKE_CXX_FLAGS "-Wno-deprecated-declarations")
set(CMAKE_CXX_FLAGS "-pthread")

include_directories(../libosmium/include)
include_directories(../protozero/include)
include_directories(../cxxopts/include)
include_directories(/usr/local/include)

link_directories(/usr/local/lib)

add_executable(createtiff createtiff.cpp)
target_link_libraries(createtiff z bz2 expat geotiff tiff)
set_property(TARGET createtiff PROPERTY CXX_STANDARD 14)

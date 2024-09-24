cmake_minimum_required(VERSION 3.5.0)
project(ray VERSION 0.1.0 LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 20)
add_library(ray SHARED
Ray.hpp
Ray.cpp
Util.hpp
Geometry.cpp
Geometry.hpp
Window.cpp
Window.hpp
Render.cpp
Render.hpp
Material.hpp
Material.cpp
)

find_package(SDL2 REQUIRED)
target_link_libraries(ray ${SDL2_LIBRARIES})


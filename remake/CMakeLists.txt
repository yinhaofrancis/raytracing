cmake_minimum_required(VERSION 3.5.0)
project(rt VERSION 0.1.0 LANGUAGES C CXX)

add_executable(rt main.cpp)

add_subdirectory(ray)

target_link_libraries(rt ray)
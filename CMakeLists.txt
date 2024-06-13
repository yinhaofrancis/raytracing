cmake_minimum_required(VERSION 3.20.0)
project(glgo VERSION 0.1.0 LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 20)
add_executable(glgo main.cpp Screen.cpp Math.cpp Render.cpp)

find_package(glfw3 REQUIRED)
find_package(OpenGL 4.6 REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(SDL2 REQUIRED)
# message(${glfw3_LIBS})

if (SDL2_FOUND)
    message("SDL2")
    target_link_libraries(glgo SDL2)
endif(SDL2_FOUND)

if (OpenGL_FOUND)
    message("opengl found ")
    include_directories(glgo ${OPENGL_INCLUDE_DIR})
    target_link_libraries (glgo ${OPENGL_opengl_LIBRARY})
endif(OpenGL_FOUND)

if (glfw3_FOUND)
  target_link_libraries (glgo glfw)
endif (glfw3_FOUND)

if (Eigen3_FOUND)
  target_link_libraries(glgo Eigen3::Eigen)
endif(Eigen3_FOUND)




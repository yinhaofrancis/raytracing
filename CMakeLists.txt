cmake_minimum_required(VERSION 3.20.0)
project(glgo VERSION 0.1.0 LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 20)
add_executable(glgo main.cpp Screen.cpp Raytracing.cpp Texture.cpp Primitives.cpp Math.cpp Render.cpp)

find_package(Eigen3 REQUIRED)
find_package(SDL2 REQUIRED)



if (APPLE)
  if (SDL2_FOUND)
  include_directories(glgo ${SDL2_INCLUDE_DIRS})

  target_link_libraries(glgo ${SDL2_LIBRARIES})
  endif(SDL2_FOUND)

  if (Eigen3_FOUND)
    target_link_libraries(glgo Eigen3::Eigen)
    INCLUDE_DIRECTORIES (glgo /opt/homebrew/Cellar/eigen/3.4.0_1/include)
    INCLUDE_DIRECTORIES (glgo /usr/local/Cellar/eigen/3.4.0_1/include)
  endif(Eigen3_FOUND)
else()
  if (SDL2_FOUND)
    message("SDL2")
    target_link_libraries(glgo SDL2)
  endif(SDL2_FOUND)
  if (Eigen3_FOUND)
    target_link_libraries(glgo Eigen3::Eigen)
  endif(Eigen3_FOUND) 
endif()







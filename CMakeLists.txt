cmake_minimum_required(VERSION 3.27)
project(Scotty3D)

set(CMAKE_CXX_STANDARD 17)

include_directories(D:/CMU/2024 Spring/15662 Computer Graphic/Scotty3D/src)

add_executable(Scotty3D tests/a0/test.a0.task2.problems.cpp
        src/scene/transform.cpp
        tests/a1/test.a1.task1.cpp)

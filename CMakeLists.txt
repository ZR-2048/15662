cmake_minimum_required(VERSION 3.27)
project(Scotty3D)

set(CMAKE_CXX_STANDARD 17)

include_directories(Scotty3D src/test.h)

add_executable(Scotty3D tests/a0/test.a0.task2.problems.cpp
        src/scene/transform.cpp
        tests/a1/test.a1.task1.cpp)

cmake_minimum_required(VERSION 3.27)
project(Scotty3D)

set(CMAKE_CXX_STANDARD 17)

include_directories(Scotty3D src/test.h)

add_executable(Scotty3D tests/a0/test.a0.task2.problems.cpp
        src/scene/transform.cpp
        tests/a3/test.a3.task2.triangle.hit.cpp
        tests/a3/test.a3.task3.bvh.build.cpp)

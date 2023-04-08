cmake_minimum_required(VERSION 3.22)
project(b_tree)

set(CMAKE_CXX_STANDARD 14)

include_directories(.)
include_directories(stx-btree-0.9/include)

add_executable(b+_tree
        test2.cpp)

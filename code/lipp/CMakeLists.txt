cmake_minimum_required(VERSION 3.22)
project(lipp)

set(CMAKE_CXX_STANDARD 14)

include_directories(.)

add_executable(lipp
#        storage_management.h
        test.cpp
#        util.h
        )

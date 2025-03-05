cmake_minimum_required(VERSION 3.16.5)
project(ep_sda)

set(CMAKE_CXX_STANDARD 20)
# set(CMAKE_BUILD_TYPE "Release")

# Set CMAKE_BUILD_TYPE to Debug for development
set(CMAKE_BUILD_TYPE Debug)

# Add compiler flags for debugging
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g")

add_executable(ep_sda
        main.cpp
        main.h
        SDA.h
        SDA.cpp)

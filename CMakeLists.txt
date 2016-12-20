cmake_minimum_required(VERSION 3.5)
project(channel_polarization)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SOURCE_FILES main.cpp main.h Encoder.cpp)
add_executable(channel_polarization ${SOURCE_FILES})
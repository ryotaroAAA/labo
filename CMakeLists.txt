cmake_minimum_required(VERSION 3.5)
project(Polar_Code)

include_directories("lib")
set(sources GLOB "src/*.cpp")
message(STATUS ${sources})
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
add_executable(Polar_Code sources)
cmake_minimum_required(VERSION 3.5)
project(channel_polarization)

file(GLOB lib "lib/*.h")
file(GLOB sources "src/*.cpp")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} --std=c++11")

set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O2 -std=c++11")
add_executable(channel_polarization ${sources} ${lib})
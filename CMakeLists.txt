cmake_minimum_required(VERSION 3.4)
project(MP4 VERSION 0.1.0)

set(CMAKE_CXX_STANDARD 20)

add_executable(examples examples.cpp)
target_include_directories(examples PUBLIC C:/eigen-3.4.0)

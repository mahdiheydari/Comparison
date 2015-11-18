cmake_minimum_required(VERSION 2.6)
project(comparison)

add_executable(comparison comparison.cpp Alignment.cpp main.cpp)

install(TARGETS comparison RUNTIME DESTINATION bin)

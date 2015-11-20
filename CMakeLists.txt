cmake_minimum_required(VERSION 2.6)
project(comparison)

# set g++ specific flags
if (CMAKE_COMPILER_IS_GNUCXX OR CMAKE_COMPILER_IS_GNUCC)
        set(CMAKE_CXX_FLAGS "-Wno-deprecated -std=c++11 -static-libgcc -static-libstdc++ -fopenmp")
        set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g3 -Wall -pedantic -Wno-long-long")
        set(CMAKE_CXX_FLAGS_RELEASE "-O3")
        set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O2 -g3")
        set(CMAKE_C_FLAGS_DEBUG "-O0 -g3 -Wall -pedantic -Wno-long-long")
endif (CMAKE_COMPILER_IS_GNUCXX OR CMAKE_COMPILER_IS_GNUCC)

add_executable(comparison comparison.cpp Alignment.cpp main.cpp)

install(TARGETS comparison RUNTIME DESTINATION bin)

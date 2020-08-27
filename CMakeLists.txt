cmake_minimum_required(VERSION 3.17)
project(pendulum)

find_package(ifopt)

set(CMAKE_CXX_STANDARD 14)

add_executable(pendulum main.cpp p_variables.cpp p_variables.h p_constraint.cpp p_constraint.h p_cost.cpp p_cost.h)

target_link_libraries(pendulum PUBLIC ifopt::ifopt_ipopt)
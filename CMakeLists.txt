cmake_minimum_required(VERSION 3.16)
project(opt_lab)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(opt_lab
        main.cpp
        matrix.cpp
        ode_solver.cpp
        opt_alg.cpp
        solution.cpp
        user_funs.cpp
)

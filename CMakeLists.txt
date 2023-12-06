cmake_minimum_required(VERSION 3.25)
project(simplexoldClion)
#Eigen3 库配置
set(CMAKE_CXX_STANDARD 17)
add_executable(simplexoldClion main.cpp common/conmmon.h
        solver/MultiGrid.cpp
        solver/MultiGrid.h)
find_package(Eigen3 CONFIG REQUIRED)
target_link_libraries(simplexoldClion  Eigen3::Eigen)


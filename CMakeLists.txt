project(simplexoldClion)
cmake_minimum_required(VERSION 3.25)

set(CMAKE_CXX_STANDARD 17)
find_package(Eigen3 CONFIG REQUIRED)
find_package(yaml-cpp CONFIG REQUIRED)
# 源文件目录
set(SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR})

# 添加可执行文件
add_executable(simplexoldClion main.cpp)
# 添加 common 库
add_subdirectory(common)
#target_link_libraries(simplexoldClion common)
# 添加 solver 库
add_subdirectory(solver)
#target_link_libraries(simplexoldClion solver)

add_subdirectory(geometry)

add_subdirectory(optimization)

add_subdirectory(physics)
add_subdirectory(viewer)

#add_subdirectory(parallel)




# 链接 Eigen3 库和 solver 库
target_link_libraries(simplexoldClion Eigen3::Eigen common solver geometry
        optimization physics viewer )


target_link_libraries(simplexoldClion  yaml-cpp)
# CMakeLists.txt (common 子目录的 CMake 文件)

# 列出 common 文件夹下的所有源文件
file(GLOB PHYSICS_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp ${CMAKE_CURRENT_SOURCE_DIR}/*.h)

# 创建 common 库，包含所有源文件
add_library(physics ${PHYSICS_SOURCES}
        SoftBodyLinearFemGrid.cpp
        SoftBodyLinearFemGrid.h
        LinearFemFunc.cpp
        LinearFemFunc.h)

# 将 common 文件夹添加到项目中
target_include_directories(physics PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

# 链接 Eigen 库
target_link_libraries(physics Eigen3::Eigen)

target_link_libraries(physics common geometry optimization solver)
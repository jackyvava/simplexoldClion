# CMakeLists.txt (common 子目录的 CMake 文件)

# 列出 common 文件夹下的所有源文件
file(GLOB PARALLEL_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/*.cu ${CMAKE_CURRENT_SOURCE_DIR}/*.h)

# 创建 common 库，包含所有源文件
add_library(parallel ${PARALLEL_SOURCES}
        AuxFuncCuda.cu
        AuxFuncCuda.cuh

)

set_target_properties(parallel PROPERTIES CUDA_ARCHITECTURES 75)  # 例如，这里指定了 GPU 架构为 75

# 将 common 文件夹添加到项目中
target_include_directories(parallel PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
# 链接 Eigen 库
#target_link_libraries(parallel Eigen3::Eigen)
#target_link_libraries(parallel  yaml-cpp)

target_link_libraries(parallel  common)
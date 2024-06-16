# CMakeList.txt: SurfaceParameterization 的 CMake 项目，在此处包括源代码并定义
# 项目特定的逻辑。
#
cmake_minimum_required (VERSION 3.25)

# Enable Hot Reload for MSVC compilers if supported.
if (POLICY CMP0141)
  cmake_policy(SET CMP0141 NEW)
  set(CMAKE_MSVC_DEBUG_INFORMATION_FORMAT "$<IF:$<AND:$<C_COMPILER_ID:MSVC>,$<CXX_COMPILER_ID:MSVC>>,$<$<CONFIG:Debug,RelWithDebInfo>:EditAndContinue>,$<$<CONFIG:Debug,RelWithDebInfo>:ProgramDatabase>>")
endif()

project ("SurfaceParameterization")

OPTION(ENABLE_OPENMESH "Enable OpenMesh" True)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

include_directories(include)
link_directories(lib)

add_compile_definitions(_USE_MATH_DEFINES)

# 定义一个字符串变量
set(OFF_DATA_WRITE_PATH ${CMAKE_SOURCE_DIR}/data/)

# 将该变量传递给编译器
add_definitions(-DOFF_DATA_WRITE_PATH="${OFF_DATA_WRITE_PATH}")

# 将源代码添加到此项目的可执行文件。
#add_executable (${PROJECT_NAME} "SurfaceParameterization.cpp" "SurfaceParameterization.h")
add_executable (${PROJECT_NAME} "main.cpp")

set (CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/cmake)

############ for test: libigl ############
# list(PREPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/tests/cmake)
# Libigl
include(libigl)
# Enable the target igl::glfw
igl_include(glfw)
igl_include(imgui)
##########################################

# TODO: 如有需要，请添加测试并安装目标。
target_link_libraries(${PROJECT_NAME} PUBLIC igl::glfw igl::imgui )



if(ENABLE_OPENMESH)
    # 定义了 USE_OpenMesh 宏
    add_definitions(-D USE_OpenMesh)
    target_link_libraries(${PROJECT_NAME} PUBLIC OpenMeshCored OpenMeshToolsd)
    message("OpenMesh is enabled.")
else()
    # 没有定义 USE_OpenMesh 宏
    message("OpenMesh is disabled.")
endif()


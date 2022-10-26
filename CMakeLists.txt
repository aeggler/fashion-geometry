cmake_minimum_required(VERSION 3.16)
project(example)

list(PREPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# Libigl
include(libigl)

# Enable the target igl::glfw
igl_include(glfw)
igl_include(imgui)
#set(PMP_VERSION 1.2.0)
#include_directories("ext/pmp-library/external/imgui")

# which directories to process


# Add your project files
file(GLOB SRC_FILES main.cpp
        toolbox/adjacency.cpp               toolbox/adjacency.h
        toolbox/Timer.cpp                   toolbox/Timer.h
        toolbox/PositionBasedDynamics.cpp   toolbox/PositionBasedDynamics.h
        toolbox/MathFunctions.cpp           toolbox/MathFunctions.h
        toolbox/constraint_utils.cpp        toolbox/constraint_utils.h
        toolbox/body_interpolation.cpp      toolbox/body_interpolation.h
#

        )
add_executable(${PROJECT_NAME} ${SRC_FILES})
target_link_libraries(${PROJECT_NAME} PUBLIC igl::glfw igl::imgui)

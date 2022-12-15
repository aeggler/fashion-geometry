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
        toolbox/garment_adaption.cpp        toolbox/garment_adaption.h
        toolbox/seam.cpp                    toolbox/seam.h
        toolbox/patternAdaption.cpp         toolbox/patternAdaption.h
#

        )
add_executable(${PROJECT_NAME} ${SRC_FILES})
target_link_libraries(${PROJECT_NAME} PUBLIC igl::glfw igl::imgui)
#
#list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)
#set(CMAKE_FIND_LIBRARY_SUFFIXES ".so" ".a" ".dylib")
find_package(GUROBI REQUIRED)
#


#In my Mac GUROBI_INCLUDE_DIRS points to "/Library/gurobi951/macos_universal2/include"

#target_link_libraries(${CMAKE_PROJECT_NAME} ${GUROBI_CXX_LIBRARY} ${GUROBI_LIBRARY})
#which links the compiled gurobi libraries.
#In my Mac GUROBI_CXX_LIBRARY points to "/Library/gurobi951/macos_universal2/lib/libgurobi_c++.a"
#and GUROBI_LIBRARY points to "/Library/gurobi951/macos_universal2/lib/libgurobi95.dylib"


target_include_directories (${PROJECT_NAME} PUBLIC ${GUROBI_INCLUDE_DIRS})
#
#if(CXX)
#    set(CMAKE_CXX_STANDARD 17)
#    target_link_libraries(example PUBLIC )
#endif()
SET(GUROBI_CXX_LIBRARY "/Library/gurobi1000/macos_universal2/lib/libgurobi_c++.a" )
target_link_libraries(example PUBLIC ${GUROBI_CXX_LIBRARY} ${GUROBI_LIBRARY})

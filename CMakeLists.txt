cmake_minimum_required(VERSION 3.20)
project(MWMuscle LANGUAGES CXX C)

set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTOUIC ON)
set(CMAKE_AUTORCC ON)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(Qt5 REQUIRED COMPONENTS Widgets)
find_package(VTK REQUIRED)
find_package(Eigen3 REQUIRED)

#file(GLOB_RECURSE MODEL_SRC ${CMAKE_CURRENT_SOURCE_DIR}/src/model/*.cpp ${CMAKE_CURRENT_SOURCE_DIR}/src/model/*.h)
#file(GLOB_RECURSE RENDERER_SRC ${CMAKE_CURRENT_SOURCE_DIR}/src/renderer/*.cpp ${CMAKE_CURRENT_SOURCE_DIR}/src/renderer/*.h)
#file(GLOB_RECURSE VIEWPORT_SRC ${CMAKE_CURRENT_SOURCE_DIR}/src/viewport/*.cpp ${CMAKE_CURRENT_SOURCE_DIR}/src/viewport/*.h)
#file(GLOB_RECURSE GUI_SRC ${CMAKE_CURRENT_SOURCE_DIR}/src/gui/*.cpp ${CMAKE_CURRENT_SOURCE_DIR}/src/gui/*.h)
#file(GLOB_RECURSE APP_SRC ${CMAKE_CURRENT_SOURCE_DIR}/src/app/*.cpp ${CMAKE_CURRENT_SOURCE_DIR}/src/app/*.h)
file(GLOB_RECURSE SIMPLESIM_SRC ${CMAKE_CURRENT_SOURCE_DIR}/src/simplesim/*.cpp ${CMAKE_CURRENT_SOURCE_DIR}/src/simplesim/*.h)

add_executable(MWMuscle
    src/main.cpp
    src/utils/rawSimulation.cpp
    src/simpleSimulation/casadiSystem.cpp
    src/simpleSimulation/SSMeshes.cpp
    src/simpleSimulation/SSMuscle.cpp
    src/simpleSimulation/VTKSimViewerSimple.cpp
    src/simpleSimulation/SSBody.cpp
    src/simpleSimulation/SSTissue.cpp
    #${MODEL_SRC}
    #${RENDERER_SRC}
    #${VIEWPORT_SRC}
    #${GUI_SRC}
    #${APP_SRC}
    ${SIMPLESIM_SRC}
)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src)

target_include_directories(MWMuscle PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/src
    ${EIGEN3_INCLUDE_DIRS}   # <- Eigen Header-only
)

target_link_libraries(MWMuscle
    Qt5::Widgets
    ${VTK_LIBRARIES}         # VTK muss hier bleiben
    casadi
)

vtk_module_autoinit(
    TARGETS MWMuscle
    MODULES ${VTK_LIBRARIES}
)

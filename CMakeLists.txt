cmake_minimum_required(VERSION 3.23)
project(implementacja)

set(CMAKE_CXX_STANDARD 14)
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/capdMake/cmake/")
include_directories(
        "/home/amrowca/CAPD/capdDynSys/include"
        "/home/amrowca/CAPD/capdAlg/include"
        "/home/amrowca/CAPD/capdAux/include"
        "/home/amrowca/CAPD/capdExt/include"
        "/home/amrowca/CAPD/capdMake/libcapd/include"
        "/home/amrowca/CAPD/capdMake/bin/include"
)

link_directories(
        "/home/amrowca/CAPD/capdDynSys/lib"
        "/home/amrowca/CAPD/capdAlg/lib"
        "/home/amrowca/CAPD/capdAux/lib"
        "/home/amrowca/CAPD/capdExt/lib"
        "/home/amrowca/CAPD/capdMake/libcapd/lib"
        "/home/amrowca/CAPD/capdMake/bin/lib"
)
find_package(CAPD REQUIRED PATHS /home/amrowca/CAPD/build)
add_executable(implementacja main.cpp)
target_link_libraries(implementacja /home/amrowca/CAPD/build/libcapd.a /home/amrowca/CAPD/build/capdExt/filibsrc/libfilib.a)



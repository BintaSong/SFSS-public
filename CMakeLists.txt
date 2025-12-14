#############################################
#               Build SFSS                  #
#############################################
cmake_minimum_required(VERSION 3.12)
project(SFSS C CXX)

find_path(CMAKE_FOLDER NAMES cmake/emp-tool-config.cmake)
include(${CMAKE_FOLDER}/cmake/emp-base.cmake)

find_package(emp-tool REQUIRED)
include_directories(${EMP-TOOL_INCLUDE_DIRS})

include_directories(/opt/homebrew/opt/gmp/include)
link_directories(/opt/homebrew/opt/gmp/lib)

# find_package(OpenMP REQUIRED)
# message(STATUS "OpenMP_CXX_FOUND: ${OpenMP_CXX_FOUND}")
# if(OpenMP_CXX_FOUND)
#     target_link_libraries(your_target_name PUBLIC OpenMP::OpenMP_CXX)
# endif()

add_executable(sfss_main "sfss_main.cpp")
target_link_libraries(sfss_main ${EMP-TOOL_LIBRARIES} gmp gmpxx)

add_executable(field "field.cpp")
target_link_libraries(field ${EMP-TOOL_LIBRARIES} gmp gmpxx)

add_executable(client_runner_main "client_runner_main.cpp")
target_link_libraries(client_runner_main ${EMP-TOOL_LIBRARIES} gmp gmpxx)

add_executable(server_runner_main "server_runner_main.cpp")
target_link_libraries(server_runner_main ${EMP-TOOL_LIBRARIES} gmp gmpxx)

set(CMAKE_CXX_FLAGS  "${CMAKE_C_FLAGS} -O3 -std=c++20") # MUST SET! Otherwise you will get -std=c++11 and cannot compile

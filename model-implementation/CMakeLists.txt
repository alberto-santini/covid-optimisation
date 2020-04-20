cmake_minimum_required(VERSION 3.15)
project(covid)

set(CMAKE_CXX_STANDARD 17)
list(INSERT CMAKE_MODULE_PATH 0 ${CMAKE_CURRENT_SOURCE_DIR})

add_executable(covid
        src/main.cpp
        src/TestLabAssignmentInstance.h
        src/TestLabAssignmentInstance.cpp
        src/TestLabAssignmentSolver.h
        src/TestLabAssignmentSolver.cpp
        src/TestLabAssignmentSolution.h
        src/TestLabAssignmentSolution.cpp)

set(WARNING_OPTIONS             -Wall -Wextra -Werror)
set(GCC_WARNING_OPTIONS         -Wduplicated-cond -Wlogical-op -Wuseless-cast -Wrestrict)
set(ARCH_OPTIONS                -march=native)
set(RELEASE_OPTIONS             -O3)
set(DEBUG_OPTIONS               -O0 -g)
set(GDB_DEBUG_OPTIONS           -ggdb)
set(LINKER_OPTIONS              -flto)

find_package(GUROBI REQUIRED)

target_include_directories( covid SYSTEM PRIVATE "src/utils")
target_compile_options(     covid PRIVATE ${WARNING_OPTIONS} ${ARCH_OPTIONS})
target_compile_options(     covid PRIVATE $<$<CXX_COMPILER_ID:GNU>:${GCC_WARNING_OPTIONS}>)
target_compile_options(     covid PRIVATE $<$<CONFIG:RELEASE>:${RELEASE_OPTIONS}>)
target_compile_options(     covid PRIVATE $<$<CONFIG:DEBUG>:${DEBUG_OPTIONS}>)
target_compile_options(     covid PRIVATE $<$<AND:$<CONFIG:DEBUG>,$<CXX_COMPILER_ID:GNU>>:${GDB_DEBUG_OPTIONS}>)
target_link_libraries(      covid PRIVATE ${LINKER_OPTIONS})
target_link_libraries(      covid PRIVATE gurobi gurobi_c++)
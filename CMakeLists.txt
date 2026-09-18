cmake_minimum_required(VERSION 3.1)

if(APPLE)
  # AppleClang's driver searches /usr/local/include before any -isystem
  # path, so a stray system-wide Eigen there would otherwise shadow
  # libigl's bundled Eigen and break its unsupported/ headers. This must
  # be set before project()/enable_language(), so it is picked up during
  # compiler detection and recorded in CMAKE_CXX_IMPLICIT_INCLUDE_DIRECTORIES;
  # CMake then skips re-adding this directory as -isystem when igl::core
  # later requests it, leaving only this higher-precedence plain -I. Not
  # gated on the directory existing yet: libigl only clones it lazily,
  # later, inside find_package(LIBIGL) below, so on a fresh checkout it
  # isn't there yet at this point -- but it exists by the time anything
  # actually gets compiled, and an extra nonexistent -I is harmless.
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -I${CMAKE_CURRENT_SOURCE_DIR}/libigl/external/eigen")
endif()

project(sparse_solver_benchmark)

option(IGL_WITH_CHOLMOD   "Use Cholmod" ON)
option(IGL_WITH_MKL       "Use Intel MKL" ON)
option(IGL_WITH_GPL       "Use GPL features in Eigen" ON)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# catamari requires C++14
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
# https://stackoverflow.com/q/46724267/148668
include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-march=native" COMPILER_SUPPORTS_MARCH_NATIVE)
if(COMPILER_SUPPORTS_MARCH_NATIVE)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
endif()

# openmp
if(APPLE AND NOT OpenMP_C_FLAGS AND EXISTS "/opt/homebrew/opt/libomp")
  # AppleClang ships no OpenMP runtime, and find_package(OpenMP) cannot
  # locate Homebrew's keg-only libomp on its own, so point it there directly.
  set(OpenMP_C_FLAGS "-Xpreprocessor -fopenmp -I/opt/homebrew/opt/libomp/include")
  set(OpenMP_C_LIB_NAMES "omp")
  set(OpenMP_CXX_FLAGS "-Xpreprocessor -fopenmp -I/opt/homebrew/opt/libomp/include")
  set(OpenMP_CXX_LIB_NAMES "omp")
  set(OpenMP_omp_LIBRARY "/opt/homebrew/opt/libomp/lib/libomp.dylib")
endif()
find_package(OPENMP REQUIRED QUIET)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")

# libigl
find_package(LIBIGL REQUIRED QUIET)

# Add your project files
if(LIBIGL_USE_STATIC_LIBRARY)
  file(GLOB SRCFILES include/igl/*.cpp)
endif()
list(APPEND SRCFILES main.cpp)
add_executable(${PROJECT_NAME} ${SRCFILES})

target_link_libraries(${PROJECT_NAME} igl::core OpenMP::OpenMP_CXX)

if(IGL_WITH_MKL)
  find_package(MKL REQUIRED)
  target_compile_definitions(${PROJECT_NAME} PUBLIC -DIGL_WITH_MKL)
  target_compile_definitions(${PROJECT_NAME} PUBLIC -DEIGEN_USE_MKL_ALL)
  include_directories(${MKL_INCLUDE_DIR})
  target_link_libraries(${PROJECT_NAME} ${MKL_LIBRARIES})
endif()

if(IGL_WITH_GPL)
  if(IGL_WITH_CHOLMOD)
    message(STATUS "Using fastest solver: Cholmod")
    target_compile_definitions(${PROJECT_NAME} PUBLIC -DIGL_WITH_CHOLMOD)
    # heavy handed way of including and linking to SuiteSparse, surely there's a
    # more modern cmake way. This still builds umfpack, spqr and klu though those
    # are not used.
    set(WITH_FORTRAN OFF CACHE BOOL "")
    set(WITH_CUDA    OFF CACHE BOOL "")
    set(WITH_OPENMP  OFF CACHE BOOL "")
    set(WITH_PARTITION OFF CACHE BOOL "")
    set(WITH_DEMOS OFF CACHE BOOL "")
    set(BUILD_CXSPARSE OFF CACHE BOOL "")
    add_subdirectory("SuiteSparse/")
    target_link_libraries(${PROJECT_NAME} cholmod umfpack)
  else()
    if(NOT DEFINED IGL_WITH_MKL)
      message(WARNING "Using slower solver: SimplicialLLT")
    endif()
    target_compile_definitions(${PROJECT_NAME} PUBLIC -DIGL_WITH_GPL)
  endif()
else()
  message(WARNING "Using slowest solver: SparseLU")
endif()


target_include_directories(${PROJECT_NAME} PUBLIC catamari/include)
target_include_directories(${PROJECT_NAME} PUBLIC mantis/include)
target_include_directories(${PROJECT_NAME} PUBLIC quotient/include)
target_compile_definitions(${PROJECT_NAME} PUBLIC -DCATAMARI_HAVE_OPENBLAS)
function(add_framework target framework)
    if(APPLE)
        find_library(found ${framework})

        if(${found} STREQUAL "${framework}-NOTFOUND")
            message(FATAL_ERROR "ERROR: ${framework} not found (${CMAKE_OSX_SYSROOT})")
        endif()

        #message("Found framework ${found}")

        target_link_libraries(${target} ${found})
    endif()
endfunction()
add_framework(${PROJECT_NAME} Accelerate)

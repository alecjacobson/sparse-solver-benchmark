cmake_minimum_required(VERSION 3.1)
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

# Find OpenMP
if(APPLE)
    #    SET(CMAKE_C_COMPILER "/usr/local/opt/llvm/bin/clang")
    #    SET(CMAKE_CXX_COMPILER "/usr/local/opt/llvm/bin/clang++")
    #    SET(OPENMP_LIBRARIES "/usr/local/opt/llvm/lib")
    #    SET(OPENMP_INCLUDES "/usr/local/opt/llvm/include")
    if(CMAKE_C_COMPILER_ID MATCHES "Clang")
        set(OpenMP_C "${CMAKE_C_COMPILER}")
        set(OpenMP_C_FLAGS "-fopenmp=libomp -Wno-unused-command-line-argument")
        set(OpenMP_C_LIB_NAMES "libomp" "libgomp" "libiomp5")
        set(OpenMP_libomp_LIBRARY ${OpenMP_C_LIB_NAMES})
        set(OpenMP_libgomp_LIBRARY ${OpenMP_C_LIB_NAMES})
        set(OpenMP_libiomp5_LIBRARY ${OpenMP_C_LIB_NAMES})
    endif()
    if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
        set(OpenMP_CXX "${CMAKE_CXX_COMPILER}")
        set(OpenMP_CXX_FLAGS "-fopenmp=libomp -Wno-unused-command-line-argument")
        set(OpenMP_CXX_LIB_NAMES "libomp" "libgomp" "libiomp5")
        set(OpenMP_libomp_LIBRARY ${OpenMP_CXX_LIB_NAMES})
        set(OpenMP_libgomp_LIBRARY ${OpenMP_CXX_LIB_NAMES})
        set(OpenMP_libiomp5_LIBRARY ${OpenMP_CXX_LIB_NAMES})
    endif()
endif()

#find_package(OPENMP)
find_package(OpenMP)

if (OPENMP_FOUND)
    #include_directories("${OPENMP_INCLUDES}")
    #include_directories("/usr/local/include" "/usr/local/opt/llvm/include")
    link_directories("${OPENMP_LIBRARIES}")
    #link_directories("/usr/local/lib" "/usr/local/opt/llvm/lib")
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    # set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()


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
    target_link_libraries(${PROJECT_NAME} cholmod)
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
target_include_directories(${PROJECT_NAME} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/sympiler-eigen/includes/)
target_include_directories(${PROJECT_NAME} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/sympiler-eigen//sympiler/sparse_blas/includes/parsy)
link_directories( ${CMAKE_CURRENT_SOURCE_DIR}/sympiler_eigen//sympiler/lib)
target_link_libraries(${PROJECT_NAME}
    lbc
    sym_sparse_blas
    sym_dense_blas
    sympiler_utils
        ${METIS_LIBRARY} ${METIS_LIBRARIES}
        )

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

add_subdirectory(sympiler-eigen)

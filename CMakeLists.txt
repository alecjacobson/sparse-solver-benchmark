cmake_minimum_required(VERSION 3.16)
project(sparse_solver_benchmark)

option(IGL_WITH_CHOLMOD   "Use Cholmod" ON)
option(IGL_WITH_MKL       "Use Intel MKL" ON)
option(IGL_WITH_GPL       "Use GPL features in Eigen" ON)
option(IGL_WITH_CUDSS     "Use NVIDIA cuDSS/cuSOLVER (requires the CUDA toolkit)" ON)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# libigl requires C++17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
# https://stackoverflow.com/q/46724267/148668
include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-march=native" COMPILER_SUPPORTS_MARCH_NATIVE)
if(COMPILER_SUPPORTS_MARCH_NATIVE)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
endif()

# openmp (optional: only used for one diagnostic thread-count print in
# main.cpp, guarded by #ifdef _OPENMP. Not REQUIRED since AppleClang on
# macOS has no OpenMP support without Homebrew's libomp + extra flags, which
# CI's macOS job intentionally doesn't install to keep that job minimal.)
find_package(OpenMP QUIET)

# libigl (header-only igl::core target; add_subdirectory instead of the old
# find_package(LIBIGL) shim which predates libigl's target-based CMake build)
add_subdirectory(libigl)

add_executable(${PROJECT_NAME} main.cpp)

target_link_libraries(${PROJECT_NAME} igl::core)
if(OpenMP_CXX_FOUND)
  target_link_libraries(${PROJECT_NAME} OpenMP::OpenMP_CXX)
endif()

if(IGL_WITH_MKL)
  # Ubuntu's intel-mkl/libmkl-full-dev package predates Intel's oneAPI
  # MKLConfig.cmake, so it's still found via our own hand-rolled finder.
  find_package(MKL REQUIRED)
  target_compile_definitions(${PROJECT_NAME} PUBLIC -DIGL_WITH_MKL)
  # Deliberately NOT defining EIGEN_USE_MKL_ALL: it routes every Eigen dense
  # matrix product through MKL's BLAS, including the many small internal
  # panel updates Eigen::SparseLU performs during factorization. Verified via
  # gdb that MKL's dgemm_direct path calls pthread_create() per small update
  # in that hot loop, turning a multi-second sparse LU into an 800+ second
  # one (~100x) from pure thread-spawn overhead. Eigen::PardisoLLT doesn't
  # need this define — it links MKL directly, independent of Eigen's own
  # dense dispatch.
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
    set(WITH_DEMOS OFF CACHE BOOL "")
    set(BUILD_CXSPARSE OFF CACHE BOOL "")
    if(IGL_WITH_MKL)
      # SuiteSparse's own find_package(BLAS) is independent of our FindMKL.cmake
      # above and would otherwise pick MKL's Intel/libiomp5 threading layer while
      # we link MKL's GNU/libgomp threading layer — two OpenMP runtimes fighting
      # in one process causes ~100x slowdowns from thread-pool contention (not a
      # measurement artifact: verified via a run where a single-threaded solver
      # got 14% CPU on one core while the process burned 350%+ across 180+
      # threads). Force sequential (non-threaded) MKL for SuiteSparse's BLAS so
      # only our own GNU-threaded MKL linkage spins up worker threads.
      set(BLA_VENDOR Intel10_64lp_seq CACHE STRING "")
    endif()
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

if(IGL_WITH_CUDSS)
  # nvcc typically isn't on PATH even when the CUDA toolkit is installed;
  # look in the standard toolkit locations before giving up.
  find_program(CMAKE_CUDA_COMPILER nvcc
    HINTS /usr/local/cuda/bin $ENV{CUDA_HOME}/bin $ENV{CUDA_PATH}/bin
    PATHS /usr/local/cuda-*/bin)
  include(CheckLanguage)
  check_language(CUDA)
  if(CMAKE_CUDA_COMPILER)
    enable_language(CUDA)
    find_package(CUDAToolkit REQUIRED)
    # The libcudss0-dev-cuda-12 apt package's own cudss-config.cmake assumes
    # headers live under the same prefix as the libraries; Debian/Ubuntu's
    # packaging splits them (/usr/include/libcudss/<ver> vs.
    # /usr/lib/<arch>/libcudss/<ver>), which trips up its relative-path
    # search. Declare the imported target ourselves instead.
    find_path(CUDSS_INCLUDE_DIR cudss.h PATH_SUFFIXES libcudss/12 libcudss/11)
    find_library(CUDSS_LIBRARY NAMES cudss PATH_SUFFIXES libcudss/12 libcudss/11)
    if(CUDSS_INCLUDE_DIR AND CUDSS_LIBRARY)
      message(STATUS "Using GPU solvers: cuDSS + cuSOLVER")
      add_library(cudss SHARED IMPORTED)
      set_target_properties(cudss PROPERTIES
        IMPORTED_LOCATION "${CUDSS_LIBRARY}"
        INTERFACE_INCLUDE_DIRECTORIES "${CUDSS_INCLUDE_DIR}")
      target_compile_definitions(${PROJECT_NAME} PUBLIC -DIGL_WITH_CUDSS)
      target_link_libraries(${PROJECT_NAME} cudss CUDA::cusolver CUDA::cusparse CUDA::cudart)
    else()
      message(WARNING "cuDSS not found (install libcudss0-dev-cuda-12); disabling GPU solvers.")
      set(IGL_WITH_CUDSS OFF)
    endif()
  else()
    message(WARNING "No CUDA compiler found; disabling GPU solvers (IGL_WITH_CUDSS).")
    set(IGL_WITH_CUDSS OFF)
  endif()
endif()


# Patch an upstream mantis bug: MultiplyAdd(float) is missing a semicolon on
# its FP_FAST_FMAF branch, which only fails to compile when FP_FAST_FMAF is
# actually defined (e.g. via -march=native enabling FMA on this machine) — so
# it's easy to pull a submodule pin that silently doesn't build here. Patched
# at configure time (idempotent) rather than by hand-editing the submodule
# checkout, so a fresh `git submodule update` doesn't reintroduce the failure.
set(_mantis_fmaf_bug "${CMAKE_CURRENT_SOURCE_DIR}/mantis/include/mantis/util-impl.hpp")
if(EXISTS "${_mantis_fmaf_bug}")
  file(READ "${_mantis_fmaf_bug}" _mantis_fmaf_bug_content)
  string(REPLACE "return std::fmaf(x, y, z)\n#elif defined(MANTIS_XLC_FUSED_MULTIPLY_ADD)"
                 "return std::fmaf(x, y, z);\n#elif defined(MANTIS_XLC_FUSED_MULTIPLY_ADD)"
                 _mantis_fmaf_bug_patched "${_mantis_fmaf_bug_content}")
  if(NOT _mantis_fmaf_bug_patched STREQUAL _mantis_fmaf_bug_content)
    file(WRITE "${_mantis_fmaf_bug}" "${_mantis_fmaf_bug_patched}")
    message(STATUS "Patched upstream mantis bug: missing semicolon in util-impl.hpp (FP_FAST_FMAF branch)")
  endif()
endif()

target_include_directories(${PROJECT_NAME} PUBLIC catamari/include)
target_include_directories(${PROJECT_NAME} PUBLIC mantis/include)
target_include_directories(${PROJECT_NAME} PUBLIC quotient/include)
if(IGL_WITH_MKL)
  # catamari's BLAS-accelerated kernels; -DCATAMARI_HAVE_OPENBLAS just picks
  # the OpenBLAS-compatible cblas_*/*_ symbol names, which MKL also provides
  # (MKL_LIBRARIES is already linked above). Without MKL, catamari falls back
  # to its portable non-BLAS kernels rather than leave unresolved symbols.
  target_compile_definitions(${PROJECT_NAME} PUBLIC -DCATAMARI_HAVE_OPENBLAS)
endif()
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

# Correctness regression tests: run on a small synthetic grid mesh (not the
# 720K-vertex dragon) and assert every compiled-in solver's L∞ residual
# stays under tolerance, for each of k=1,2,3.
enable_testing()
add_test(NAME harmonic_residuals
  COMMAND ${PROJECT_NAME} --grid 20 --check)

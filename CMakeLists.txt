cmake_minimum_required(VERSION 3.16)
project(sparse_solver_benchmark)

option(IGL_WITH_CHOLMOD   "Use Cholmod" ON)
option(IGL_WITH_MKL       "Use Intel MKL" ON)
option(IGL_WITH_GPL       "Use GPL features in Eigen" ON)
option(IGL_WITH_CUDSS     "Use NVIDIA cuDSS/cuSOLVER (requires the CUDA toolkit)" ON)
option(IGL_WITH_NASOQ     "Use NASOQ's LBL (parallel symmetric-indefinite) solver" ON)
option(IGL_WITH_CATAMARI  "Use catamari's SparseLDL solver (needs the mantis/quotient submodules too)" ON)
option(IGL_WITH_MA57      "Use MA57 (symla), a from-scratch header-only symmetric-indefinite solver" ON)
option(IGL_WITH_ACCELERATE_SPARSE "Use Apple's Accelerate Sparse Solvers library (macOS-only)" ON)

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
  # MKLConfig.cmake, so it's still found via our own hand-rolled finder. Not
  # REQUIRED: autoconfigure to whatever's actually installed rather than
  # hard-failing cmake configure on a machine without MKL (PardisoLLT/LDLT,
  # NASOQ, and catamari's BLAS-accelerated kernels are all gated behind
  # IGL_WITH_MKL and simply don't appear on the leaderboard when it's off).
  find_package(MKL)
  if(NOT MKL_FOUND)
    message(WARNING "MKL not found; disabling IGL_WITH_MKL (skips PardisoLLT/LDLT, NASOQ, and catamari's BLAS kernels).")
    set(IGL_WITH_MKL OFF)
  endif()
endif()
if(IGL_WITH_MKL)
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

if(IGL_WITH_CHOLMOD AND NOT EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/SuiteSparse/CMakeLists.txt")
  message(WARNING "SuiteSparse/ submodule not checked out; disabling IGL_WITH_CHOLMOD "
    "(run `git submodule update --init SuiteSparse` to enable it).")
  set(IGL_WITH_CHOLMOD OFF)
endif()

# --- CUDA toolkit detection, shared by CHOLMOD's own GPU-accelerated
# factorization (SuiteSparse's WITH_CUDA, below) and IGL_WITH_CUDSS's
# cuDSS/cuSOLVER (further down) -- both need CMAKE_CUDA_COMPILER set BEFORE
# add_subdirectory("SuiteSparse/") runs, since SuiteSparse's own CMakeLists
# calls enable_language(CUDA) itself based on whether CMAKE_CUDA_COMPILER is
# already found at that point.
#
# nvcc typically isn't on PATH even when a CUDA toolkit is installed. Besides
# the standard system locations, also probe two conventions used on GPU dev
# machines that don't have a system/apt CUDA install: a toolkit tree
# assembled under $HOME (e.g. for Warp's native builds -- see CLAUDE.md's
# WARP_CUDA_PATH note) and the CUDA_ROOT env var some sibling projects on
# such machines already use to point at exactly that. All of these are
# no-ops (glob to nothing, hint unused) on a machine that doesn't have them.
set(IGL_HAVE_CUDA_COMPILER OFF)
if(IGL_WITH_CUDSS OR IGL_WITH_CHOLMOD)
  file(GLOB _home_cuda_bin_dirs "$ENV{HOME}/warp-cuda/*/bin" "$ENV{HOME}/projects/cuda*/bin")
  find_program(CMAKE_CUDA_COMPILER nvcc
    HINTS /usr/local/cuda/bin $ENV{CUDA_HOME}/bin $ENV{CUDA_PATH}/bin $ENV{CUDA_ROOT}/bin
    PATHS /usr/local/cuda-*/bin ${_home_cuda_bin_dirs})
  include(CheckLanguage)
  check_language(CUDA)
  if(CMAKE_CUDA_COMPILER)
    enable_language(CUDA)
    # enable_language(CUDA) leaves CMAKE_CUDA_HOST_COMPILER defined-but-empty
    # in this project's scope; SuiteSparse's own CMakeLists checks
    # `if(DEFINED CMAKE_CUDA_HOST_COMPILER)` before falling back to
    # CMAKE_CXX_COMPILER, so it inherits that empty value instead of our
    # actual g++ and passes nvcc a blank `--compiler-bindir=`, which fails
    # with "nvcc fatal: Failed to preprocess host compiler properties." Set
    # it explicitly so both our own CUDA usage and SuiteSparse's inherit a
    # real host compiler.
    set(CMAKE_CUDA_HOST_COMPILER "${CMAKE_CXX_COMPILER}")
    find_package(CUDAToolkit REQUIRED)
    set(IGL_HAVE_CUDA_COMPILER ON)

    # On a machine whose CUDA math libraries come from pip (nvidia-cusolver-cu12/
    # nvidia-cusparse-cu12/nvidia-cublas-cu12/nvidia-cudss-cu12 -- common as
    # transitive PyTorch/Warp dependencies, e.g. via `uv`/`pip install torch`)
    # rather than a system/apt CUDA install, `find_package(CUDAToolkit)` above
    # only gives us what's actually inside the nvcc toolkit tree (nvcc itself,
    # cudart) -- CUDA::cusolver/CUDA::cusparse/CUDA::cublas targets won't
    # exist in that case, since those libraries live in a completely separate
    # pip site-packages tree, not next to nvcc. Ask python where that tree is
    # and fall back to it. Those wheels also ship only version-suffixed .so
    # files (e.g. libcusolver.so.11, no unversioned libcusolver.so symlink),
    # which find_library's normal lib<name>.so naming convention won't match
    # even with the right search path -- file(GLOB) sidesteps that.
    set(_pip_nvidia_dir "")
    find_package(Python3 QUIET COMPONENTS Interpreter)
    if(Python3_EXECUTABLE)
      execute_process(
        COMMAND ${Python3_EXECUTABLE} -c "import nvidia, os; print(os.path.dirname(nvidia.__file__))"
        OUTPUT_VARIABLE _pip_nvidia_dir OUTPUT_STRIP_TRAILING_WHITESPACE
        ERROR_QUIET)
    endif()

    if(TARGET CUDA::cusolver)
      set(CUSOLVER_TARGET CUDA::cusolver)
    elseif(_pip_nvidia_dir)
      find_path(CUSOLVER_INCLUDE_DIR cusolverSp.h HINTS "${_pip_nvidia_dir}/cusolver/include")
      file(GLOB CUSOLVER_LIBRARY "${_pip_nvidia_dir}/cusolver/lib/libcusolver.so*")
      if(CUSOLVER_INCLUDE_DIR AND CUSOLVER_LIBRARY)
        list(GET CUSOLVER_LIBRARY 0 CUSOLVER_LIBRARY)
        add_library(pip_cusolver SHARED IMPORTED)
        set_target_properties(pip_cusolver PROPERTIES
          IMPORTED_LOCATION "${CUSOLVER_LIBRARY}"
          INTERFACE_INCLUDE_DIRECTORIES "${CUSOLVER_INCLUDE_DIR}")
        set(CUSOLVER_TARGET pip_cusolver)
      endif()
    endif()

    if(TARGET CUDA::cusparse)
      set(CUSPARSE_TARGET CUDA::cusparse)
    elseif(_pip_nvidia_dir)
      find_path(CUSPARSE_INCLUDE_DIR cusparse.h HINTS "${_pip_nvidia_dir}/cusparse/include")
      file(GLOB CUSPARSE_LIBRARY "${_pip_nvidia_dir}/cusparse/lib/libcusparse.so*")
      if(CUSPARSE_INCLUDE_DIR AND CUSPARSE_LIBRARY)
        list(GET CUSPARSE_LIBRARY 0 CUSPARSE_LIBRARY)
        add_library(pip_cusparse SHARED IMPORTED)
        set_target_properties(pip_cusparse PROPERTIES
          IMPORTED_LOCATION "${CUSPARSE_LIBRARY}"
          INTERFACE_INCLUDE_DIRECTORIES "${CUSPARSE_INCLUDE_DIR}")
        set(CUSPARSE_TARGET pip_cusparse)
      endif()
    endif()

    # SuiteSparse's own CMakeLists hardcodes the literal target name
    # CUDA::cublas (target_link_libraries/$<TARGET_PROPERTY:CUDA::cublas,...>
    # in its CHOLMOD GPU block), so -- unlike cusolver/cusparse above, which
    # we reference via our own CUSOLVER_TARGET/CUSPARSE_TARGET indirection --
    # this one must be created under that exact name if missing.
    if(NOT TARGET CUDA::cublas AND _pip_nvidia_dir)
      find_path(CUBLAS_INCLUDE_DIR cublas_v2.h HINTS "${_pip_nvidia_dir}/cublas/include")
      file(GLOB CUBLAS_LIBRARY "${_pip_nvidia_dir}/cublas/lib/libcublas.so*")
      if(CUBLAS_INCLUDE_DIR AND CUBLAS_LIBRARY)
        list(GET CUBLAS_LIBRARY 0 CUBLAS_LIBRARY)
        add_library(CUDA::cublas SHARED IMPORTED GLOBAL)
        set_target_properties(CUDA::cublas PROPERTIES
          IMPORTED_LOCATION "${CUBLAS_LIBRARY}"
          INTERFACE_INCLUDE_DIRECTORIES "${CUBLAS_INCLUDE_DIR}")
      endif()
    endif()
  endif()
endif()

if(IGL_WITH_GPL)
  if(IGL_WITH_CHOLMOD)
    message(STATUS "Using fastest solver: Cholmod")
    target_compile_definitions(${PROJECT_NAME} PUBLIC -DIGL_WITH_CHOLMOD)
    # heavy handed way of including and linking to SuiteSparse, surely there's a
    # more modern cmake way. This still builds umfpack, spqr and klu though those
    # are not used.
    set(WITH_FORTRAN OFF CACHE BOOL "")
    # GPU-accelerated supernodal factorization (CHOLMOD's own CUDA path, not
    # to be confused with cuDSS/cuSOLVER below -- this is literally CHOLMOD
    # offloading its own dense supernode BLAS3 updates to the GPU via cuBLAS).
    # Needs both a CUDA compiler and CUDA::cublas (see the detection block
    # above); gracefully stays on CPU-only CHOLMOD otherwise.
    if(IGL_HAVE_CUDA_COMPILER AND TARGET CUDA::cublas)
      set(WITH_CUDA ON CACHE BOOL "")
    else()
      set(WITH_CUDA OFF CACHE BOOL "")
    endif()
    set(WITH_OPENMP  OFF CACHE BOOL "")
    set(WITH_DEMOS OFF CACHE BOOL "")
    set(BUILD_CXSPARSE OFF CACHE BOOL "")
    # This SuiteSparse fork's CHOLMOD renamed the old WITH_PARTITION toggle to
    # WITH_METIS (default ON), which this project's CMakeLists never updated
    # to match, so METIS-based nested-dissection ordering was silently
    # enabled. A/B timing on this machine (full dragon mesh, k=1..3) showed it
    # makes no measurable difference for k=1,2 (within run-to-run noise) but a
    # ~4x factorization slowdown for k=3 (1.1s vs 4.4s) -- CHOLMOD's default
    # AMD-based ordering is a much better fit than METIS for the heavy fill-in
    # the triharmonic system produces. Disable it.
    set(WITH_METIS OFF CACHE BOOL "")
    if(IGL_WITH_MKL)
      # SuiteSparse's own find_package(BLAS) is independent of our FindMKL.cmake
      # above and could in principle pick a different MKL threading layer than
      # our own linkage, which was once diagnosed (via gdb) as a ~100x
      # thread-pool-contention slowdown and "fixed" by forcing SuiteSparse's
      # BLAS to Intel10_64lp_seq (sequential/non-threaded). That fix was a
      # bigger hammer than needed: it also serializes CHOLMOD/UMFPACK's own
      # internal BLAS3 calls, and A/B timing on this machine (full dragon
      # mesh, k=1..3) showed plain threaded MKL BLAS is faster across the
      # board (5.7/25/38s vs. 6.5/31/56s factor time) with no reappearance of
      # the original contention. Leave BLAS threaded; if the contention bug
      # resurfaces, this is the first place to look.
      set(BLA_VENDOR Intel10_64lp CACHE STRING "")
    endif()
    add_subdirectory("SuiteSparse/")
    target_link_libraries(${PROJECT_NAME} cholmod umfpack)
    # WITH_CUDA may have been forced back OFF inside SuiteSparse's own
    # cmake_dependent_option handling (e.g. WITH_LICENSE mismatch) even
    # though we requested it above -- check the actual outcome, not our
    # request, before telling main.cpp the GPU CHOLMOD path exists.
    if(WITH_CUDA)
      message(STATUS "CHOLMOD built with its own CUDA-accelerated supernodal factorization")
      target_compile_definitions(${PROJECT_NAME} PUBLIC -DIGL_WITH_CHOLMOD_CUDA)
    endif()
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
  # CMAKE_CUDA_COMPILER/CUDAToolkit/CUSOLVER_TARGET/CUSPARSE_TARGET/
  # _pip_nvidia_dir were already resolved by the shared CUDA-detection block
  # above (before SuiteSparse's add_subdirectory, which needs them earlier).
  if(IGL_HAVE_CUDA_COMPILER)
    # The libcudss0-dev-cuda-12 apt package's own cudss-config.cmake assumes
    # headers live under the same prefix as the libraries; Debian/Ubuntu's
    # packaging splits them (/usr/include/libcudss/<ver> vs.
    # /usr/lib/<arch>/libcudss/<ver>), which trips up its relative-path
    # search. Declare the imported target ourselves instead -- also checking
    # the pip nvidia-cudss-cu12 tree (versioned .so, same GLOB workaround as
    # cusolver/cusparse above) alongside the apt/system locations.
    find_path(CUDSS_INCLUDE_DIR cudss.h PATH_SUFFIXES libcudss/12 libcudss/11
      HINTS "${_pip_nvidia_dir}/cu12/include" "${_pip_nvidia_dir}/cu13/include")
    find_library(CUDSS_LIBRARY NAMES cudss PATH_SUFFIXES libcudss/12 libcudss/11)
    if(NOT CUDSS_LIBRARY AND _pip_nvidia_dir)
      file(GLOB CUDSS_LIBRARY "${_pip_nvidia_dir}/cu12/lib/libcudss.so*" "${_pip_nvidia_dir}/cu13/lib/libcudss.so*")
      list(GET CUDSS_LIBRARY 0 CUDSS_LIBRARY)
    endif()
    if(CUDSS_INCLUDE_DIR AND CUDSS_LIBRARY AND CUSOLVER_TARGET AND CUSPARSE_TARGET)
      message(STATUS "Using GPU solvers: cuDSS + cuSOLVER")
      add_library(cudss SHARED IMPORTED)
      set_target_properties(cudss PROPERTIES
        IMPORTED_LOCATION "${CUDSS_LIBRARY}"
        INTERFACE_INCLUDE_DIRECTORIES "${CUDSS_INCLUDE_DIR}")
      target_compile_definitions(${PROJECT_NAME} PUBLIC -DIGL_WITH_CUDSS)
      target_link_libraries(${PROJECT_NAME} cudss ${CUSOLVER_TARGET} ${CUSPARSE_TARGET} CUDA::cudart)
    else()
      message(WARNING "cuDSS/cuSOLVER/cuSPARSE not found (install libcudss0-dev-cuda-12, or "
        "pip install nvidia-cudss-cu12 nvidia-cusolver-cu12 nvidia-cusparse-cu12); "
        "disabling GPU solvers.")
      set(IGL_WITH_CUDSS OFF)
    endif()
  else()
    message(WARNING "No CUDA compiler found; disabling GPU solvers (IGL_WITH_CUDSS).")
    set(IGL_WITH_CUDSS OFF)
  endif()
endif()

if(IGL_WITH_NASOQ AND NOT EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/nasoq/CMakeLists.txt")
  message(WARNING "nasoq/ submodule not checked out; disabling IGL_WITH_NASOQ "
    "(run `git submodule update --init nasoq` to enable it).")
  set(IGL_WITH_NASOQ OFF)
endif()

if(IGL_WITH_NASOQ)
  if(IGL_WITH_MKL)
    # NASOQ's own cmake/FindMKL.cmake has the same missing-header-search-hint
    # bug our original cmake/FindMKL.cmake had (FIND_PATH with no HINTS, so
    # it never looks in /usr/include/mkl); pre-seed the cache var so its
    # find_path is a no-op instead of failing.
    set(MKL_INCLUDE_DIR /usr/include/mkl CACHE PATH "")
    set(NASOQ_BLAS_BACKEND MKL CACHE STRING "")
    set(NASOQ_BUILD_CLI OFF CACHE BOOL "")
    set(NASOQ_BUILD_EXAMPLES OFF CACHE BOOL "")
    set(NASOQ_BUILD_TESTS OFF CACHE BOOL "")
    set(NASOQ_BUILD_DOCS OFF CACHE BOOL "")
    set(NASOQ_WITH_EIGEN ON CACHE BOOL "")
    set(NASOQ_WITH_MATLAB OFF CACHE BOOL "")
    set(NASOQ_USE_CLAPACK OFF CACHE BOOL "")
    add_subdirectory(nasoq)
    target_compile_definitions(${PROJECT_NAME} PUBLIC -DIGL_WITH_NASOQ)
    target_link_libraries(${PROJECT_NAME} nasoq::eigen_interface)
  else()
    message(WARNING "NASOQ LBL requires IGL_WITH_MKL for its BLAS backend; disabling.")
    set(IGL_WITH_NASOQ OFF)
  endif()
endif()

if(IGL_WITH_MA57 AND NOT EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/ma57/CMakeLists.txt")
  message(WARNING "ma57/ submodule not checked out; disabling IGL_WITH_MA57 "
    "(run `git submodule update --init ma57` to enable it).")
  set(IGL_WITH_MA57 OFF)
endif()

if(IGL_WITH_MA57)
  # symla (ma57/) is a header-only INTERFACE target with its own CMakeLists
  # that hard find_package(Eigen3 ... REQUIRED)s and FATAL_ERRORs if
  # SuiteSparse AMD isn't found system-wide -- neither of which is otherwise
  # a dependency of this project (libigl fetches its own Eigen; our own
  # CHOLMOD/UMFPACK come from the SuiteSparse/ submodule, not system
  # libsuitesparse-dev). Probe for both ourselves first so a machine without
  # them just skips MA57 instead of hard-failing cmake configure.
  find_package(Eigen3 3.4 QUIET NO_MODULE)
  find_path(MA57_SUITESPARSE_INCLUDE_DIR amd.h PATH_SUFFIXES suitesparse)
  find_library(MA57_AMD_LIBRARY NAMES amd)
  if(NOT Eigen3_FOUND OR NOT MA57_SUITESPARSE_INCLUDE_DIR OR NOT MA57_AMD_LIBRARY)
    message(WARNING "MA57 (symla) requires system Eigen3 (libeigen3-dev) and SuiteSparse AMD "
      "(libsuitesparse-dev); not found -- disabling IGL_WITH_MA57.")
    set(IGL_WITH_MA57 OFF)
  endif()
endif()

if(IGL_WITH_MA57)
  set(SYMLA_BUILD_TESTS OFF CACHE BOOL "")
  set(SYMLA_BUILD_BENCH OFF CACHE BOOL "")
  set(SYMLA_WITH_SUITESPARSE_BENCH OFF CACHE BOOL "")
  add_subdirectory(ma57)
  target_compile_definitions(${PROJECT_NAME} PUBLIC -DIGL_WITH_MA57)
  target_link_libraries(${PROJECT_NAME} symla::symla)
endif()


if(IGL_WITH_CATAMARI)
  foreach(_dep catamari mantis quotient)
    if(NOT EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${_dep}/include")
      message(WARNING "${_dep}/ submodule not checked out; disabling IGL_WITH_CATAMARI "
        "(run `git submodule update --init catamari mantis quotient` to enable it).")
      set(IGL_WITH_CATAMARI OFF)
    endif()
  endforeach()
endif()

if(IGL_WITH_CATAMARI)
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

  target_compile_definitions(${PROJECT_NAME} PUBLIC -DIGL_WITH_CATAMARI)
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

if(IGL_WITH_ACCELERATE_SPARSE AND NOT APPLE)
  message(WARNING "Accelerate's Sparse Solvers library is macOS-only; disabling IGL_WITH_ACCELERATE_SPARSE.")
  set(IGL_WITH_ACCELERATE_SPARSE OFF)
endif()
if(IGL_WITH_ACCELERATE_SPARSE)
  # Header-only (Accelerate/Sparse/Solve.h, macOS 10.13+); Accelerate itself
  # is already linked above via add_framework, so nothing further to link.
  target_compile_definitions(${PROJECT_NAME} PUBLIC -DIGL_WITH_ACCELERATE_SPARSE)
endif()

if(IGL_WITH_NASOQ AND OpenMP_CXX_FOUND AND UNIX AND NOT APPLE AND CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  # Force GNU libgomp to win OpenMP symbol resolution for this executable,
  # even though NASOQ's own MKL discovery (nasoq/cmake/third_party/mkl.cmake,
  # via find_package(MKL) with the default MKL_THREADING=OMP) also links
  # libmkl_intel_thread.so, which needs Intel's OpenMP runtime (libiomp5.so
  # -- on Debian/Ubuntu a symlink straight to LLVM's libomp.so.5, confirm
  # with `readlink -f /usr/lib/x86_64-linux-gnu/libiomp5.so`). Without this,
  # `ld`'s default --as-needed drops the earlier-declared OpenMP::OpenMP_CXX
  # (libgomp) from this executable's DT_NEEDED entirely, once libiomp5/
  # libomp.so.5 (pulled in later, transitively, by nasoq::eigen_interface)
  # has already satisfied every GOMP_*/omp_* symbol GCC's own -fopenmp code
  # needs -- confirmed directly with `ldd`: without this fix, the built
  # binary depends on libomp.so.5 only, not libgomp at all, even though
  # most of this binary's own OpenMP-using code (symla/MA57, and MKL's own
  # mkl_gnu_thread.so) was compiled/intended for libgomp's ABI.
  #
  # That matters far beyond just "which .so provides the symbol": on the
  # L40/dual-Xeon dev machine this was diagnosed on, LLVM's libomp.so.5 was
  # measured to be catastrophically slower than GNU's libgomp specifically
  # for CPU-bound work executed while many sibling threads sit idle inside
  # a large active OpenMP team (e.g. symla's multifrontal factorize()
  # parallel driver) -- MA57's factor time dropped from 39-56s to 3-6s on
  # k=1 after this fix. A minimal, solver-independent repro (57,741 small
  # Eigen matrix ops, one thread busy inside a 128-thread team) measured
  # 0.7s under libgomp vs. 16s under libomp.so.5, scaling worse with team
  # size; not explained by allocator contention, and unaffected by
  # OMP_WAIT_POLICY=PASSIVE or KMP_BLOCKTIME=0 -- looks like a genuine
  # LLVM-libomp idle-thread-in-a-team pathology on that specific machine.
  # Only observed/tested there; scoped narrowly here (Linux + GCC + NASOQ,
  # the only combination that actually triggers the underlying link-order
  # issue in this project's CI matrix) since forcing `-lgomp` unconditionally
  # would break Clang/macOS/Windows builds where no such library exists.
  #
  # `--no-as-needed`/`--as-needed` bracket only `gomp` itself (not everything
  # linked after it), so this doesn't change --as-needed pruning behavior
  # for anything else in the link. Both runtimes still coexist in the
  # process afterward -- this only changes which one wins the *shared*
  # (GOMP_*/omp_*) symbol names; libiomp5/libomp.so.5 remains available
  # (transitively, via libmkl_intel_thread.so's own DT_NEEDED) for MKL's own
  # Intel-ABI (__kmpc_*) calls, which libgomp cannot satisfy.
  target_link_libraries(${PROJECT_NAME} -Wl,--no-as-needed gomp -Wl,--as-needed)
endif()

# Correctness regression tests: run on a small synthetic grid mesh (not the
# 720K-vertex dragon) and assert every compiled-in solver's L∞ residual
# stays under tolerance, for each of k=1,2,3.
enable_testing()
add_test(NAME harmonic_residuals
  COMMAND ${PROJECT_NAME} --grid 20 --check)

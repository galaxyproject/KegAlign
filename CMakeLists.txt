# 3.30, not 3.10: CMP0167 below is a 3.30 policy and cmake_policy(SET ...) is an ERROR on a CMake
# that does not know it, so 3.30 has been the real floor for as long as that line has been here.
cmake_minimum_required (VERSION 3.30)
cmake_policy(SET CMP0104 NEW)   # CMake 3.18
cmake_policy(SET CMP0128 NEW)   # CMake 3.22
cmake_policy(SET CMP0146 NEW)   # CMake 3.27
cmake_policy(SET CMP0167 NEW)   # CMake 3.30

execute_process(COMMAND bash "-c" "${CMAKE_SOURCE_DIR}/scripts/get-cuda-arches.bash"
                OUTPUT_VARIABLE KEGALIGN_CUDA_ARCHITECTURES
                RESULT_VARIABLE KEGALIGN_ARCH_SCRIPT_RESULT
                ERROR_VARIABLE  KEGALIGN_ARCH_SCRIPT_STDERR)
set(CMAKE_CUDA_ARCHITECTURES "${KEGALIGN_CUDA_ARCHITECTURES}")

# ---- TEMPORARY DIAGNOSTIC, NOT FOR MERGE -------------------------------------
# Which architectures actually reach nvcc? The shipped 0.1.2.10 and 0.1.2.11
# binaries both record 14 in __CUDA_ARCH_LIST__, but this script emits 19. Print
# every candidate source so we can see who drops the other five.
message(STATUS "KEGALIGN_DIAG script_exit=[${KEGALIGN_ARCH_SCRIPT_RESULT}]")
message(STATUS "KEGALIGN_DIAG script_stderr=[${KEGALIGN_ARCH_SCRIPT_STDERR}]")
message(STATUS "KEGALIGN_DIAG script_stdout=[${KEGALIGN_CUDA_ARCHITECTURES}]")
message(STATUS "KEGALIGN_DIAG env_CUDAARCHS=[$ENV{CUDAARCHS}]")
message(STATUS "KEGALIGN_DIAG env_cuda_compiler_version=[$ENV{cuda_compiler_version}]")
message(STATUS "KEGALIGN_DIAG before_project=[${CMAKE_CUDA_ARCHITECTURES}]")
# ------------------------------------------------------------------------------

project(kegalign LANGUAGES C CXX CUDA)

# ---- TEMPORARY DIAGNOSTIC, NOT FOR MERGE -------------------------------------
message(STATUS "KEGALIGN_DIAG after_project=[${CMAKE_CUDA_ARCHITECTURES}]")
list(LENGTH CMAKE_CUDA_ARCHITECTURES KEGALIGN_DIAG_N)
message(STATUS "KEGALIGN_DIAG after_project_count=[${KEGALIGN_DIAG_N}]")
message(STATUS "KEGALIGN_DIAG cmake_version=[${CMAKE_VERSION}] policy_floor=[${CMAKE_MINIMUM_REQUIRED_VERSION}]")
# ------------------------------------------------------------------------------

set(CMAKE_CUDA_COMPILER nvcc)

set(CMAKE_CUDA_COMPILER_ENV_VAR "CUDACXX")
set(CMAKE_CUDA_HOST_COMPILER_ENV_VAR "CUDAHOSTCXX")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O4 -std=c++17 -l c -lboost_program_options")

find_package(TBB REQUIRED tbbmalloc tbbmalloc_proxy tbb_preview)

find_package(ZLIB REQUIRED)
find_path(ZLIB_INCLUDE_DIR zlib.h)
include_directories(${ZLIB_INCLUDE_DIR})

find_package(Boost COMPONENTS program_options REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

add_executable(kegalign
    common/DRAM.cpp
    common/ntcoding.cpp
    common/seed_filter_interface.cu
    common/seed_pos_table.cu
    common/utilities.c
    common/dna_utilities.c
    common/scoring.c
    src/seed_filter.cu
    src/seeder.cpp
    src/segment_printer.cpp
    src/main.cpp
    )


set(LIBS rt stdc++ ${TBB_IMPORTED_TARGETS} ${Boost_LIBRARIES} ${CUDA_LIBRARIES} ${ZLIB_LIBRARIES})

set_target_properties(
    kegalign 
    PROPERTIES 
      CUDA_SEPARABLE_COMPILATION ON
      # CUDA 13 bundles a CCCL whose Thrust, CUB and libcu++ headers all refuse to compile below
      # C++17. REQUIRED so a toolchain that cannot provide it fails here rather than quietly
      # dropping to an older standard and failing later inside those headers.
      CUDA_STANDARD 17
      CUDA_STANDARD_REQUIRED ON
      )

target_compile_options(kegalign PRIVATE -DTBB_SUPPRESS_DEPRECATED_MESSAGES)

target_link_libraries(kegalign PRIVATE ${LIBS})

target_include_directories(kegalign PRIVATE ${PROJECT_SOURCE_DIR}/common/ ${PROJECT_SOURCE_DIR}/src/)

# ---- TEMPORARY DIAGNOSTIC, NOT FOR MERGE -------------------------------------
# The target property is what actually becomes -gencode; it can differ from the
# directory variable above.
get_target_property(KEGALIGN_DIAG_TGT kegalign CUDA_ARCHITECTURES)
message(STATUS "KEGALIGN_DIAG target_property=[${KEGALIGN_DIAG_TGT}]")
list(LENGTH KEGALIGN_DIAG_TGT KEGALIGN_DIAG_TGT_N)
message(STATUS "KEGALIGN_DIAG target_property_count=[${KEGALIGN_DIAG_TGT_N}]")
# ------------------------------------------------------------------------------



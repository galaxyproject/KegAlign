cmake_minimum_required (VERSION 3.10)

# if(POLICY ...) around each of these, because cmake_policy(SET ...) on a policy
# the running CMake has never heard of is an error. Guarding them is what lets
# this file name a low floor and still ask for NEW behaviour where it is
# available -- raising cmake_minimum_required is not the way to do that, and
# doing so once changed which CUDA architectures this build used.
if(POLICY CMP0104)
    cmake_policy(SET CMP0104 NEW)   # CMake 3.18
endif()
if(POLICY CMP0128)
    cmake_policy(SET CMP0128 NEW)   # CMake 3.22
endif()
if(POLICY CMP0146)
    cmake_policy(SET CMP0146 NEW)   # CMake 3.27
endif()
if(POLICY CMP0167)
    cmake_policy(SET CMP0167 NEW)   # CMake 3.30
endif()

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")
include(SelectCudaArchitectures)
kegalign_select_cuda_architectures(CMAKE_CUDA_ARCHITECTURES)

project(kegalign LANGUAGES C CXX CUDA)

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


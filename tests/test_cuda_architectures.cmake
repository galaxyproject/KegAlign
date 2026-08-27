# Tests for kegalign_select_cuda_architectures().
#
# Run with:  cmake -P tests/test_cuda_architectures.cmake
#
# No compiler, no CUDA toolkit and no GPU are needed: the function is pure list
# handling plus one execute_process, and the tests substitute a stub for that.

cmake_minimum_required(VERSION 3.10)

get_filename_component(_repo "${CMAKE_CURRENT_LIST_DIR}/.." ABSOLUTE)
include("${_repo}/cmake/SelectCudaArchitectures.cmake")

set(_failures 0)

function(expect_eq label actual expected)
  if(NOT "${actual}" STREQUAL "${expected}")
    message("not ok - ${label}\n     got: [${actual}]\nexpected: [${expected}]")
    math(EXPR _failures "${_failures} + 1")
    set(_failures "${_failures}" PARENT_SCOPE)
  else()
    message("ok - ${label}")
  endif()
endfunction()

# A stub standing in for get-cuda-arches.bash, so the tests never invoke nvcc.
set(_stub "${CMAKE_CURRENT_LIST_DIR}/stub-arches.bash")

# ---------------------------------------------------------------------------
# The bug this all exists for: conda-forge exports CUDAARCHS, and KegAlign used
# to overwrite it unconditionally with its own longer list. That list then blew
# past what CCCL's namespace macro can expand, and Thrust failed to compile.
# An architecture list supplied by the environment must win.
# ---------------------------------------------------------------------------
set(ENV{CUDAARCHS} "50-real;52-real;90a-real;121-virtual")
unset(CMAKE_CUDA_ARCHITECTURES)
kegalign_select_cuda_architectures(got SCRIPT "${_stub}")
expect_eq("CUDAARCHS from the environment is used verbatim"
          "${got}" "50-real;52-real;90a-real;121-virtual")

# An explicit caller choice (-DCMAKE_CUDA_ARCHITECTURES=... or a toolchain
# file) outranks even CUDAARCHS.
set(ENV{CUDAARCHS} "50-real;52-real")
set(CMAKE_CUDA_ARCHITECTURES "80-real;90")
kegalign_select_cuda_architectures(got SCRIPT "${_stub}")
expect_eq("an explicit CMAKE_CUDA_ARCHITECTURES outranks CUDAARCHS"
          "${got}" "80-real;90")
unset(CMAKE_CUDA_ARCHITECTURES)

# With nothing supplied, fall back to asking nvcc.
unset(ENV{CUDAARCHS})
kegalign_select_cuda_architectures(got SCRIPT "${_stub}")
expect_eq("falls back to the script when nothing is supplied"
          "${got}" "75-real;80-real;86")

# An empty CUDAARCHS is not a choice; it must not shadow the fallback.
set(ENV{CUDAARCHS} "")
kegalign_select_cuda_architectures(got SCRIPT "${_stub}")
expect_eq("an empty CUDAARCHS falls through to the script"
          "${got}" "75-real;80-real;86")
unset(ENV{CUDAARCHS})

# ---------------------------------------------------------------------------
# Counting, which is what actually decides whether the build survives. CCCL
# expands THRUST_NAMESPACE_BEGIN as
#     _CCCL_PP_SPLICE_WITH(_, THRUST, THRUST_VERSION, SM, __CUDA_ARCH_LIST__, NS)
# so the argument count is 4 + (distinct architectures). Measured against the
# 12.9.27 headers with cpp alone: <=15 architectures expand correctly, 16 is
# silently garbled, 17-26 is a hard error, >=27 is silently garbled again.
# -real/-virtual suffixes and the a/f variants all collapse to one architecture.
# ---------------------------------------------------------------------------
kegalign_count_cuda_architectures(n "50-real;52-real;90a-real;90-real;121-virtual;121")
expect_eq("suffixes and a/f variants collapse to distinct architectures" "${n}" "4")

kegalign_count_cuda_architectures(n "100f-real;100a-real;100")
expect_eq("f and a variants of one architecture count once" "${n}" "1")

# conda-forge's own CUDA 12.9 list -- 14, which is why those builds were green.
kegalign_count_cuda_architectures(n
  "50-real;52-real;60-real;61-real;70-real;75-real;80-real;86-real;89-real;90a-real;100f-real;103f-real;120a-real;121f-real;121-virtual")
expect_eq("conda-forge's 12.9 list is 14 architectures" "${n}" "14")

# What get-cuda-arches.bash produces for 12.9 -- 19, over the cliff.
kegalign_count_cuda_architectures(n
  "50-real;52-real;53-real;60-real;61-real;62-real;70-real;72-real;75-real;80-real;86-real;87-real;89-real;90-real;100-real;101-real;103-real;120-real;121")
expect_eq("get-cuda-arches.bash's 12.9 list is 19 architectures" "${n}" "19")

# The limit itself, so a future toolkit that adds architectures trips this test
# rather than the CI build.
expect_eq("the documented safe limit is 15" "${KEGALIGN_CCCL_MAX_ARCHITECTURES}" "15")

# ---------------------------------------------------------------------------
# The guard must actually stop the build. Run it in a child cmake, since it
# reports the problem with FATAL_ERROR.
# ---------------------------------------------------------------------------
set(_over "${CMAKE_CURRENT_LIST_DIR}/over-limit.cmake")
file(WRITE "${_over}"
"include(\"${_repo}/cmake/SelectCudaArchitectures.cmake\")\n"
"set(CMAKE_CUDA_ARCHITECTURES \"50;52;53;60;61;62;70;72;75;80;86;87;89;90;100;101;103;120;121\")\n"
"kegalign_select_cuda_architectures(got)\n")
execute_process(COMMAND "${CMAKE_COMMAND}" -P "${_over}"
                RESULT_VARIABLE _rc ERROR_VARIABLE _err OUTPUT_QUIET)
file(REMOVE "${_over}")
if(_rc EQUAL 0)
  message("not ok - 19 architectures must be refused, but configure succeeded")
  math(EXPR _failures "${_failures} + 1")
elseif(NOT _err MATCHES "19 CUDA architectures")
  message("not ok - refusal did not name the count\n  got: ${_err}")
  math(EXPR _failures "${_failures} + 1")
else()
  message("ok - 19 architectures is refused with a legible message")
endif()

if(_failures GREATER 0)
  message(FATAL_ERROR "${_failures} test(s) failed")
endif()
message("\nall tests passed")

# Choosing the CUDA architectures to build for.
#
# WHY THIS IS NOT JUST `set(CMAKE_CUDA_ARCHITECTURES <our list>)`.
#
# Packagers already know which architectures they want. conda-forge exports
# CUDAARCHS, which CMake reads, and its list is maintained by the people who run
# the CUDA migrations -- it uses the a/f variants where they matter (90a, 100f,
# 103f, 120a, 121f) and ends in 121-virtual so PTX is embedded for future GPUs.
# Overwriting that with a list of our own is both presumptuous and, as it turned
# out, dangerous: the longer list broke the build. So an architecture list that
# reaches us from outside always wins, and get-cuda-arches.bash is the fallback
# for people building this by hand.
#
# THE LIMIT. CCCL opens Thrust's ABI namespace with
#
#     _CCCL_PP_SPLICE_WITH(_, THRUST, THRUST_VERSION, SM, __CUDA_ARCH_LIST__, NS)
#
# whose argument count is 4 + (distinct architectures). CCCL's argument counter
# does not survive arbitrary widths. Measured directly against the CUDA 12.9.27
# headers, with the preprocessor alone:
#
#     <= 15 architectures   expands correctly
#        16 architectures   expands to the WRONG namespace, silently
#     17 - 26 architectures   hard error inside Thrust and CUB
#     >= 27 architectures   silently wrong again
#
# 15 is therefore the last width known to be correct, and going over it is
# checked here so the failure is one legible message rather than a thousand
# lines of macro errors -- or, at exactly 16, no error at all and a binary whose
# Thrust symbols are in a namespace nobody meant.
set(KEGALIGN_CCCL_MAX_ARCHITECTURES 15)

# Captured at include time. CMAKE_CURRENT_FUNCTION_LIST_DIR would be the natural
# thing to use inside the function below, but it only exists from CMake 3.17 and
# this file must not quietly depend on more than cmake_minimum_required claims.
set(KEGALIGN_CMAKE_MODULE_DIR "${CMAKE_CURRENT_LIST_DIR}")

# Count the distinct architectures in a CMAKE_CUDA_ARCHITECTURES value.
# `-real` and `-virtual` suffixes, and the `a`/`f` variants, all name the same
# architecture as far as __CUDA_ARCH_LIST__ is concerned: 90, 90a and 90-real
# are one architecture, not three.
function(kegalign_count_cuda_architectures out_var architectures)
    set(_numbers "")
    foreach(_entry IN LISTS architectures)
        string(REGEX REPLACE "-(real|virtual)$" "" _entry "${_entry}")
        string(REGEX REPLACE "[af]$" "" _entry "${_entry}")
        if(_entry MATCHES "^[0-9]+$")
            list(APPEND _numbers "${_entry}")
        endif()
    endforeach()
    if(_numbers)
        list(REMOVE_DUPLICATES _numbers)
    endif()
    list(LENGTH _numbers _count)
    set(${out_var} "${_count}" PARENT_SCOPE)
endfunction()

# kegalign_select_cuda_architectures(<out-var> [SCRIPT <path>])
#
# Sets <out-var> to the architectures to build for, in precedence order:
#   1. CMAKE_CUDA_ARCHITECTURES already set by the caller (-D, or a toolchain file)
#   2. the CUDAARCHS environment variable
#   3. scripts/get-cuda-arches.bash, which asks nvcc
function(kegalign_select_cuda_architectures out_var)
    cmake_parse_arguments(PARSE_ARGV 1 _arg "" "SCRIPT" "")

    if(NOT _arg_SCRIPT)
        set(_arg_SCRIPT "${KEGALIGN_CMAKE_MODULE_DIR}/../scripts/get-cuda-arches.bash")
    endif()

    if(DEFINED CMAKE_CUDA_ARCHITECTURES AND NOT CMAKE_CUDA_ARCHITECTURES STREQUAL "")
        set(_architectures "${CMAKE_CUDA_ARCHITECTURES}")
        set(_source "the caller")
    elseif(DEFINED ENV{CUDAARCHS} AND NOT "$ENV{CUDAARCHS}" STREQUAL "")
        set(_architectures "$ENV{CUDAARCHS}")
        set(_source "the CUDAARCHS environment variable")
    else()
        execute_process(
            COMMAND bash "-c" "${_arg_SCRIPT}"
            OUTPUT_VARIABLE _architectures
            RESULT_VARIABLE _result
            ERROR_VARIABLE _stderr)
        if(NOT _result EQUAL 0)
            message(FATAL_ERROR
                "${_arg_SCRIPT} failed (exit ${_result}): ${_stderr}")
        endif()
        string(STRIP "${_architectures}" _architectures)
        if(_architectures STREQUAL "")
            message(FATAL_ERROR
                "${_arg_SCRIPT} produced no architectures. Set CUDAARCHS or pass "
                "-DCMAKE_CUDA_ARCHITECTURES=... to choose them explicitly.")
        endif()
        set(_source "get-cuda-arches.bash")
    endif()

    kegalign_count_cuda_architectures(_count "${_architectures}")
    if(_count GREATER KEGALIGN_CCCL_MAX_ARCHITECTURES)
        message(FATAL_ERROR
            "Building for ${_count} CUDA architectures, from ${_source}, but CCCL's "
            "Thrust/CUB namespace macro only expands correctly for up to "
            "${KEGALIGN_CCCL_MAX_ARCHITECTURES}. Above that it either fails deep "
            "inside the Thrust headers or, at exactly 16, silently produces the "
            "wrong namespace. Narrow the list with CUDAARCHS or "
            "-DCMAKE_CUDA_ARCHITECTURES=...\n"
            "  architectures: ${_architectures}")
    endif()

    message(STATUS
        "KegAlign: building for ${_count} CUDA architectures, from ${_source}")
    set(${out_var} "${_architectures}" PARENT_SCOPE)
endfunction()

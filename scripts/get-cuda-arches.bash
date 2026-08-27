#!/usr/bin/env bash
#
# Emit a CMAKE_CUDA_ARCHITECTURES value covering every GPU this nvcc can target.
#
# ASK NVCC, DO NOT MAINTAIN A TABLE. `nvcc --list-gpu-code` reports the architectures the compiler
# in front of us actually supports, so the list cannot go stale, and it is right on the day a new
# CUDA release lands rather than whenever somebody notices.
#
# The table this replaced was transcribed from NVIDIA's documentation, and it drifted exactly the way
# transcriptions do. Its CUDA 13.0 row was the 12.9 row with the pre-Turing entries trimmed off the
# front: it still listed 101/101a/101f, which 13.0 REMOVED, and never gained 88 or 110/110a/110f,
# which 13.0 ADDED. A build against 13.0 therefore asked for an architecture that no longer exists.
# The failure surfaced only on the next CUDA upgrade, which is the worst time to discover it.
#
# ⚠ ARCHITECTURE-SPECIFIC TARGETS ARE DELIBERATELY NOT BUILT. `--list-gpu-code` reports only the
# "non-architecture-specific" architectures -- sm_100, never sm_100a or sm_100f. That is what we
# want: the `a` and `f` variants exist for code using architecture-specific instructions, and
# KegAlign's kernels use none (no __CUDA_ARCH__ guards, no wgmma/tcgen05, no cluster APIs). Building
# them would triple the compile time and the fatbinary for identical generic code and zero extra
# device coverage, since sm_100 runs on every SM 10.0 device. If a kernel ever does use such an
# instruction, add that one target here for that stated reason -- nvcc will not enumerate them.
#
# The highest architecture is emitted WITHOUT `-real`, so its PTX is embedded and can be JIT-compiled
# for GPUs newer than this toolkit knows about. Every other architecture is `-real`: SASS only, no
# redundant PTX.

set -o errexit
set -o nounset
set -o pipefail

# in case cuda_compiler_version isn't already set
if [ -z ${cuda_compiler_version+x} ]; then
    type -p nvcc &> /dev/null || {
        >&2 echo "error: unable to find nvcc command"
        exit 1
    }

    cuda_compiler_version=$(nvcc --version | sed -n 's/^.*release \([0-9]\+\.[0-9]\+\).*$/\1/p')
fi

declare -a CUDA_CONFIG_ARGS
if [ "${cuda_compiler_version}" != "None" ]; then
    type -p nvcc &> /dev/null || {
        >&2 echo "error: unable to find nvcc command"
        exit 1
    }

    # sm_75 sm_80 sm_86 ... -> 75 80 86 ...
    declare -a ARCHES
    mapfile -t ARCHES < <(nvcc --list-gpu-code | sed -n 's/^sm_\([0-9]\+\)$/\1/p' | sort -n -u)

    if [ ${#ARCHES[@]} -eq 0 ]; then
        >&2 echo "error: nvcc --list-gpu-code reported no usable architectures"
        >&2 echo "       cuda_compiler_version=${cuda_compiler_version}"
        exit 1
    fi

    LATEST_ARCH="${ARCHES[-1]}"
    unset "ARCHES[${#ARCHES[@]}-1]"

    CMAKE_CUDA_ARCHS=""
    for arch in ${ARCHES[@]+"${ARCHES[@]}"}; do
        CMAKE_CUDA_ARCHS="${CMAKE_CUDA_ARCHS:+${CMAKE_CUDA_ARCHS};}${arch}-real"
    done

    CMAKE_CUDA_ARCHS="${CMAKE_CUDA_ARCHS:+${CMAKE_CUDA_ARCHS};}${LATEST_ARCH}"

    CUDA_CONFIG_ARGS+=(
        "${CMAKE_CUDA_ARCHS}"
    )
fi

echo -n ${CUDA_CONFIG_ARGS+"${CUDA_CONFIG_ARGS[@]}"}

#!/usr/bin/env bash

set -o errexit
set -o nounset

# in case cuda_compiler_version isn't already set
if [ -z ${cuda_compiler_version+x} ]; then
    type -p nvcc &> /dev/null || {
        >&2 echo "error: unable to find nvcc command"
        exit 1
    }

    cuda_compiler_version=$(nvcc --version | sed -n 's/^.*release \([0-9]\+\.[0-9]\+\).*$/\1/p')
fi

# function to facilitate version comparison; cf. https://stackoverflow.com/a/37939589
version2int () { echo "$@" | awk -F. '{ printf("%d%02d\n", $1, $2); }'; }

declare -a CUDA_CONFIG_ARGS
if [ "${cuda_compiler_version}" != "None" ]; then
    cuda_compiler_version_int=$(version2int "$cuda_compiler_version")

    ARCHES=()
    if   [ $cuda_compiler_version_int -ge $(version2int "13.0") ]; then # 2025-08
        ARCHES=(75 80 86 87 89 90 90a 100 100a 100f 101 101a 101f 103 103a 103f 120 120a 120f 121 121a 121f)
    elif [ $cuda_compiler_version_int -ge $(version2int "12.9") ]; then # 2025-05
        ARCHES=(50 52 53 60 61 62 70 72 75 80 86 87 89 90 90a 100 100a 100f 101 101a 101f 103 103a 103f 120 120a 120f 121 121a 121f)
    elif [ $cuda_compiler_version_int -ge $(version2int "12.8") ]; then # 2025-01
        ARCHES=(50 52 53 60 61 62 70 72 75 80 86 87 89 90 90a 100 100a 101 101a 120 120a)
    elif [ $cuda_compiler_version_int -ge $(version2int "12.6") ]; then # 2024-08
        ARCHES=(50 52 53 60 61 62 70 72 75 80 86 87 89 90 90a)
    elif [ $cuda_compiler_version_int -ge $(version2int "11.8") ]; then # 2022-10
        ARCHES=(35 37 50 52 53 60 61 62 70 72 75 80 86 87 89 90)
    elif [ $cuda_compiler_version_int -ge $(version2int "11.7") ]; then # 2022-05
        ARCHES=(35 37 50 52 53 60 61 62 70 72 75 80 86 87)
    elif [ $cuda_compiler_version_int -ge $(version2int "11.4") ]; then # 2021-06
        ARCHES=(35 37 50 52 53 60 61 62 70 72 75 80 86)
    elif [ $cuda_compiler_version_int -ge $(version2int "11.1") ]; then # 2020-09
        ARCHES=(35 37 50 52 53 60 61 62 70 72 75 80)
    elif [ $cuda_compiler_version_int -ge $(version2int "10.0") ]; then # 2018-09
        ARCHES=(30 32 35 50 52 53 60 61 62 70 72 75)
    elif [ $cuda_compiler_version_int -ge $(version2int "9.0") ]; then # 2017-09
        ARCHES=(30 32 35 50 52 53 60 61 62 70)
    elif [ $cuda_compiler_version_int -ge $(version2int "8.0") ]; then # 2017-02
        ARCHES=(20 30 32 35 50 52 53)
    fi

    LATEST_ARCH="${ARCHES[-1]}"
    unset "ARCHES[${#ARCHES[@]}-1]"

    for arch in "${ARCHES[@]}"; do
        CMAKE_CUDA_ARCHS="${CMAKE_CUDA_ARCHS+${CMAKE_CUDA_ARCHS};}${arch}-real"
    done

    CMAKE_CUDA_ARCHS="${CMAKE_CUDA_ARCHS+${CMAKE_CUDA_ARCHS};}${LATEST_ARCH}"

    CUDA_CONFIG_ARGS+=(
        "${CMAKE_CUDA_ARCHS}"
    )
fi

echo -n ${CUDA_CONFIG_ARGS+"${CUDA_CONFIG_ARGS[@]}"}

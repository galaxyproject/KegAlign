#!/usr/bin/env bash

set -o errexit
set -o nounset
set -o pipefail

project_root="."
conda_env_file="conda-env.bash"
conda_env_dev_file="conda-env-dev.bash"
conda_root="$project_root/.conda"

test -e "$conda_root" || mkdir -p "$conda_root"

dev=0
if [ $# -eq 1 ]; then
    if [ "$1" = "dev" -o  "$1" = "-dev" -o "$1" = "--dev" ]; then
        dev=1
    fi
fi

##
## download the latest Miniforge
##

miniforge_root="$conda_root/miniforge3"
if [ ! -e "$miniforge_root" ]; then
    repo="conda-forge/miniforge"
    tag=$(curl --silent "https://api.github.com/repos/$repo/releases/latest" | jq -r ".tag_name // empty")
    if [ -z "$tag" ]; then
        echo "unable to get latest tag for $repo"
        exit 1
    fi
    filename="Miniforge3-${tag}-$(uname -s)-$(uname -m).sh"
    url="https://github.com/$repo/releases/download/$tag/$filename"

    curl --silent --location -o "$conda_root/$filename" "$url"
    curl --silent --location -o "$conda_root/${filename}.sha256" "${url}.sha256"
    expected_checksum=$(awk '{print $1}' "$conda_root/${filename}.sha256")
    actual_checksum=$(openssl dgst -sha256 "$conda_root/$filename" | awk '{print $2}')

    if [ "$expected_checksum" != "$actual_checksum" ]; then
        echo "Checksum mismatch: $conda_root/$filename"
        exit 1
    fi

    ##
    ## install Miniforge
    ##

    bash "$conda_root/$filename" -b -p "$miniforge_root"
    rm -f "$conda_root/$filename" "$conda_root/${filename}.sha256"
fi

if [ $dev -eq 0 ]; then
    echo "source \"${miniforge_root}/etc/profile.d/conda.sh\"" > "$conda_env_file"
    echo "source \"${miniforge_root}/etc/profile.d/mamba.sh\"" >> "$conda_env_file"
    source "$conda_env_file"
    echo "conda activate kegalign" >> "$conda_env_file"

    installed=$(conda env list | grep -E "^kegalign\b" || true)
    if [ -z "$installed" ]; then
        mamba create \
            --name kegalign \
            --channel conda-forge \
            --channel bioconda \
            --channel defaults \
            --override-channels \
            --strict-channel-priority \
            --yes \
            "kegalign-full"
    fi
elif [ $dev -eq 1 ]; then
    echo "source \"${miniforge_root}/etc/profile.d/conda.sh\"" > "$conda_env_dev_file"
    echo "source \"${miniforge_root}/etc/profile.d/mamba.sh\"" >> "$conda_env_dev_file"
    source "$conda_env_dev_file"
    echo "conda activate kegalign-dev" >> "$conda_env_dev_file"

    installed=$(conda env list | grep -E "^kegalign-dev\b" || true)
    if [ -z "$installed" ]; then
        mamba create \
            --name kegalign-dev \
            --channel conda-forge \
            --channel bioconda \
            --channel defaults \
            --override-channels \
            --strict-channel-priority \
            --yes \
            "black" \
            "cmake" \
            "flake8" \
            "gxx=11" \
            "libboost-devel>=1.70" \
            "mypy" \
            "nvidia-ml-py" \
            "python=3.12" \
            "tbb-devel=2020.2.*" \
            "zlib"
    fi
fi


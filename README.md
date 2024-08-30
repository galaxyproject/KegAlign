[license-badge]: https://img.shields.io/badge/License-MIT-yellow.svg 
[license-link]: https://opensource.org/licenses/MIT

[![License][license-badge]][license-link]
[![Published in SC20](https://img.shields.io/badge/published%20in-SC20-blue.svg)](https://doi.ieeecomputersociety.org/10.1109/SC41405.2020.00043)
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/kegalign-full/README.html)

<img src="kegalign_logo.webp" width="300">

This is a [@galaxyproject](https://github.com/galaxyproject)'s modified fork of the original [SegAlign](https://github.com/gsneha26/SegAlign). 

## Table of Contents

- [Overview](#overview)
- [Installation](#installation)
  - [Dependencies](#dependencies)
- [Usage](#usage)
  - [Alignment](#alignment)
  - [Scoring Options](#scoring)
  - [Output Options](#output)
- [Citing KegAlign](#cite_kegalign)

## <a name="overview"></a> Overview

Precise genome aligner efficiently leveraging GPUs.

### <a name="changes"></a> Changes from the original implementation

- Added advanced runner script allowing the usage of MIG and/or MPS for better GPU utilization
- Updated to compile with TBB (Threading Building Blocks) [version 2020.2](https://github.com/oneapi-src/oneTBB/releases/tag/v2020.2)
- Fixed the **--scoring** option. It can now read and use the substitution matrix from a LASTZ [Scoring File](https://lastz.github.io/lastz/#fmt_scoring)
- Added **--num_threads** option to limit the number of threads used
- Added **--segment_size** option to limit maximum number of HSPs per segment file for CPU load balancing
- Cleaned up build files and addressed compiler warnings

## <a name="installation"></a> Installation

For standalone installation use Conda: `conda install conda-forge::kegalign`
For standalone installation with additional tools use Bioconda: `conda install bioconda::kegalign-full`

For installation in Galaxy we currently use the wrappers [`richard-burhans:kegalign`](https://toolshed.g2.bx.psu.edu/view/richard-burhans/kegalign
) and [`richard-burhans:batched_lastz`](https://toolshed.g2.bx.psu.edu/view/richard-burhans/batched_lastz
) from the [Main Tool Shed](https://toolshed.g2.bx.psu.edu/).
Try the tools at usegalaxy.org: [kegalign](https://usegalaxy.org/root?tool_id=toolshed.g2.bx.psu.edu/repos/richard-burhans/kegalign/kegalign), [batched_lastz](https://usegalaxy.org/root?tool_id=toolshed.g2.bx.psu.edu/repos/richard-burhans/batched_lastz/batched_lastz)

* Script to create conda environment

```bash
git clone https://github.com/galaxyproject/KegAlign.git
cd KegAlign
./scripts/make-conda-env.bash
source ./conda-env.bash
```

* Script to install development enviroment

```bash
git clone https://github.com/galaxyproject/KegAlign.git
cd KegAlign
./scripts/make-conda-env.bash -dev
source ./conda-env-dev.bash

mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make
```
## <a name="dependencies"></a> Dependencies
The following dependencies are required by KegAlign:

  * [CMake](https://cmake.org/) >= 3.8
  * oneAPI Threading Building Blocks ([oneTBB](https://oneapi-src.github.io/oneTBB/)) [2020.2](https://github.com/oneapi-src/oneTBB/releases/tag/v2020.2)
  * [Boost C++ Libraries](https://www.boost.org/) >= 1.70
  * [LASTZ](https://github.com/lastz/lastz) 1.04.22
  * faToTwoBit (from [UCSC Genome Browser source](https://github.com/ucscGenomeBrowser/kent))

## <a name="usage"></a> Usage

## <a name="alignment"></a> Alignment

#### Running a Sample Alignment
```bash
# install kegalign
git clone https://github.com/galaxyproject/KegAlign.git
cd KegAlign
./scripts/make-conda-env.bash
source ./conda-env.bash

# convert target (ref) and query to 2bit
mkdir work
faToTwoBit <(gzip -cdfq ./test-data/apple.fasta.gz) work/ref.2bit
faToTwoBit <(gzip -cdfq ./test-data/orange.fasta.gz) work/query.2bit

# generate LASTZ keg
python ./scripts/runner.py --diagonal-partition --format maf- --num-cpu 16 --num-gpu 1 --output-file data_package.tgz --output-type tarball --tool_directory ./scripts test-data/apple.fasta.gz test-data/orange.fasta.gz
python ./scripts/package_output.py --format_selector maf --tool_directory ./scripts

# run LASTZ keg
python ./scripts/run_lastz_tarball.py --input=data_package.tgz --output=apple_orange.maf --parallel=16

# check output
diff apple_orange.maf <(gzip -cdfq ./test-data/apple_orange.maf.gz)

# command-line kegalign
kegalign test-data/apple.fasta.gz test-data/orange.fasta.gz work/ --num_gpu 1 --num_threads 16 > lastz-commands.txt
bash lastz-commands.txt
(echo "##maf version=1"; cat *.maf-) > apple_orange.maf
```

#### Running with MIG/MPS
GPU utilization can be increased by using MIG and/or MPS, leading up to 20% faster alignments.

* Preparing inputs

```bash
mkdir query_split target_split
./scripts/mps-mig/split_input.py --input ./test-data/apple.fasta.gz --out query_split --to_2bit --goal_bp 20000000
./scripts/mps-mig/split_input.py --input ./test-data/orange.fasta.gz --out target_split --to_2bit --goal_bp 20000000
mkdir tmp
```

* Select GPU UUIDs to run on using

```bash
nvidia-smi -L
```

* run on two GPUs with 4 MPS processes per GPU (replace [GPU-UUID#] with outputs from above command)

```bash
python ./scripts/mps-mig/run_mig.py [GPU-UUID1],[GPU-UUID2] --MPS 4 --target ./target_split --query ./query_split  --tmp_dir ./tmp/ --mps_pipe_dir ./tmp/ --output ./apples_oranges.maf --num_threads 64
```

### <a name="scoring"></a>Scoring Options

By default the HOXD70 substitution scores are used (from [Chiaromonte et al. 2002](https://doi.org/10.1142/9789812799623_0012))

    bad_score          = X:-1000  # used for sub['X'][*] and sub[*]['X']
    fill_score         = -100     # used when sub[*][*] is not defined
    gap_open_penalty   =  400
    gap_extend_penalty =   30

         A     C     G     T
    A   91  -114   -31  -123
    C -114   100  -125   -31
    G  -31  -125   100  -114
    T -123   -31  -114    91

Matrix can be supplied as an input to **--scoring** parameter. Substitution matrix can be inferred from your data using another LASTZ-based tool (LASTZ_D: Infer substitution scores).

### <a name="output"></a>Output Options

The default output is a MAF alignment file. Other formats can be selected with the **--format** parameter.  See [LASTZ manual](https://lastz.github.io/lastz/#formats) for description of possible formats.

## <a name="cite_kegalign"></a> Citing KegAlign


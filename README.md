[license-badge]: https://img.shields.io/badge/License-MIT-yellow.svg 
[license-link]: https://opensource.org/licenses/MIT

[![License][license-badge]][license-link]
[![Published in SC20](https://img.shields.io/badge/published%20in-SC20-blue.svg)](https://doi.ieeecomputersociety.org/10.1109/SC41405.2020.00043)

<img src="logo.png" width="300">

This is a [@galaxyproject](https://github.com/galaxyproject)'s modified fork of the original [SegAlign](https://github.com/gsneha26/SegAlign). 

## Table of Contents

- [Overview](#overview)
- [Installation](#installation)
- [Usage](#usage)
  - [Alignment](#alignment)
  - [Scoring Options](#scoring)
  - [Output Options](#output)
- [Citing SegAlign](#cite_segalign)

## <a name="overview"></a> Overview

Precise genome aligner efficiently leveraging GPUs.

### <a name="what"></a> What it does

SegAlign processes **Target** and **Query** sequences to identify highly similar regions where gapped extension will be performed to create actual alignments. 

### <a name="changes"></a> Changes from the original implementation

- Cleaned up build files and addressed compiler warnings
- Updated to compile with TBB (Threading Building Blocks) [version 2020.2](https://github.com/oneapi-src/oneTBB/releases/tag/v2020.2)
- Fixed the --scoring option. It can now read and use the substitution matrix from a LASTZ [Scoring File](https://lastz.github.io/lastz/#fmt_scoring)
- Added --num_threads option to limit the number of threads used
- Added --segment_size option to limit maximum number of HSPs per segment file for CPU load balancing.
- Added optional runner script using MIG and/or MPS for better GPU utilization.

## <a name="usage"></a> Usage

### <a name="alignment"></a>Alignment

#### Running a Sample Alignment
```
cd $PROJECT_DIR
mkdir test
cd test
wget https://usegalaxy.org/api/datasets/f9cad7b01a47213501e23cde09bc3eb2/display?to_ext=fasta 
wget https://usegalaxy.org/api/datasets/f9cad7b01a4721352af44e7304057a1c/display?to_ext=fasta 
run_segalign ce11.fa cb4.fa --output=ce11.cb4.maf --num_gpus 1 --num_threads 32 --segment_size -1
```
* For a list of options run
```
run_segalign --help
```
#### Running with MIG/MPS 
GPU utilization can be increased by using MIG and/or MPS, leading up to 20% faster alignments.
* Preparing inputs
```
wget https://usegalaxy.org/api/datasets/f9cad7b01a47213501e23cde09bc3eb2/display?to_ext=fasta 
wget https://usegalaxy.org/api/datasets/f9cad7b01a4721352af44e7304057a1c/display?to_ext=fasta
mkdir query_split
mkdir target_split
split_input.py --input apple.fasta --out ./query_split --to_2bit true --goal_bp 20000000
split_input.py --input orange.fasta --out ./target_split --to_2bit true --goal_bp 20000000
mkdir tmp
```
* Select GPU UUIDs to run on using
```
nvidia-smi -L
```
* run on two GPUs with 4 MPS processes per GPU (replace [GPU-UUID#] with outputs from above command)
```
run_mig.py [GPU-UUID1],[GPU-UUID2] --MPS 4 --target ./target_split --query ./query_split  --tmp_dir ./tmp/ --mps_pipe_dir ./tmp/ --output ./apples_oranges.maf --num_threads 64
```


### <a name="scoring"></a>Scoring Options

By default the HOXD70 substitution scores are used (from `Chiaromonte et al. 2002 <https://www.ncbi.nlm.nih.gov/pubmed/11928468>`_)::

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

The default output is a MAF alignment file. Other formats can be selected with the *--format* parameter.  See LASTZ manual <https://lastz.github.io/lastz/#formats>`_ for description of possible formats.

## <a name="installation"></a> Installation

For standalone installation use Conda: `conda install conda-forge::segalign`

For installation in Galaxy we currently use the wrappers `richard-burhans:segalign` and `richard-burhans:batched_lastz` from the [Main Tool Shed](https://toolshed.g2.bx.psu.edu/).
Try the tools at usegalaxy.org: [segalign](https://usegalaxy.org/root?tool_id=toolshed.g2.bx.psu.edu/repos/richard-burhans/segalign/segalign/), [batched_lastz](https://usegalaxy.org/root?tool_id=toolshed.g2.bx.psu.edu/repos/richard-burhans/batched_lastz/batched_lastz/)

## <a name="cite_segalign"></a> Citing SegAlign

S. Goenka, Y. Turakhia, B. Paten and M. Horowitz,  "SegAlign: A Scalable GPU-Based Whole Genome Aligner," in 2020 SC20: International Conference for High Performance Computing, Networking, Storage and Analysis (SC), Atlanta, GA, US, 2020 pp. 540-552. doi: 10.1109/SC41405.2020.00043

[license-badge]: https://img.shields.io/badge/License-MIT-yellow.svg 
[license-link]: https://opensource.org/licenses/MIT

[![License][license-badge]][license-link]
[![Published in SC20](https://img.shields.io/badge/published%20in-SC20-blue.svg)](https://doi.ieeecomputersociety.org/10.1109/SC41405.2020.00043)

<img src="kegalign_logo.webp" width="300">

This is a [@galaxyproject](https://github.com/galaxyproject)'s modified fork of the original [SegAlign](https://github.com/gsneha26/SegAlign). 

## Table of Contents

- [Overview](#overview)
- [Dependencies](#dependencies)
- [Installation](#installation)
- [Usage](#usage)
  - [Alignment](#alignment)
  - [Scoring Options](#scoring)
  - [Output Options](#output)
- [Citing SegAlign](#cite_segalign)

## <a name="overview"></a> Overview

Precise genome aligner efficiently leveraging GPUs.

## <a name="dependencies"></a> Dependencies
The following dependencies are required by SegAlign:

  * [CMake](https://cmake.org/) >= 3.8
  * oneAPI Threading Building Blocks ([oneTBB](https://oneapi-src.github.io/oneTBB/)) [2020.2](https://github.com/oneapi-src/oneTBB/releases/tag/v2020.2)
  * [Boost C++ Libraries](https://www.boost.org/) >= 1.70
  * [LASTZ](https://github.com/lastz/lastz) 1.04.22
  * faToTwoBit (from [UCSC Genome Browser source](https://github.com/ucscGenomeBrowser/kent))

### <a name="changes"></a> Changes from the original implementation

- Cleaned up build files and addressed compiler warnings
- Updated to compile with TBB (Threading Building Blocks) [version 2020.2](https://github.com/oneapi-src/oneTBB/releases/tag/v2020.2)
- Fixed the **--scoring** option. It can now read and use the substitution matrix from a LASTZ [Scoring File](https://lastz.github.io/lastz/#fmt_scoring)
- Added **--num_threads** option to limit the number of threads used
- Added **--segment_size** option to limit maximum number of HSPs per segment file for CPU load balancing
- Added optional runner script using MIG and/or MPS for better GPU utilization

## <a name="installation"></a> Installation

For standalone installation use Conda: `conda install conda-forge::segalign`

For installation in Galaxy we currently use the wrappers `richard-burhans:segalign` and `richard-burhans:batched_lastz` from the [Main Tool Shed](https://toolshed.g2.bx.psu.edu/).
Try the tools at usegalaxy.org: [segalign](https://usegalaxy.org/root?tool_id=toolshed.g2.bx.psu.edu/repos/richard-burhans/segalign/segalign/), [batched_lastz](https://usegalaxy.org/root?tool_id=toolshed.g2.bx.psu.edu/repos/richard-burhans/batched_lastz/batched_lastz/)

## <a name="usage"></a> Usage

## <a name="alignment"></a> Alignment

#### Running a Sample Alignment

#### Running with MIG/MPS

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

## <a name="cite_segalign"></a> Citing SegAlign

S. Goenka, Y. Turakhia, B. Paten and M. Horowitz,  "SegAlign: A Scalable GPU-Based Whole Genome Aligner," in 2020 SC20: International Conference for High Performance Computing, Networking, Storage and Analysis (SC), Atlanta, GA, US, 2020 pp. 540-552. doi: 10.1109/SC41405.2020.00043

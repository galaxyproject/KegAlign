[license-badge]: https://img.shields.io/badge/License-MIT-yellow.svg 
[license-link]: https://opensource.org/licenses/MIT

[![License][license-badge]][license-link]
[![Published in SC20](https://img.shields.io/badge/published%20in-SC20-blue.svg)](https://doi.ieeecomputersociety.org/10.1109/SC41405.2020.00043)

<img src="logo.png" width="300">

This is a [@galaxyproject](https://github.com/galaxyproject)'s modified fork of the original [SegAlign](https://github.com/gsneha26/SegAlign). 

## Table of Contents

- [Overview](#overview)
- [Installation](#installation)
- [Citing SegAlign](#cite_segalign)

## <a name="overview"></a> Overview

Precise genome aligner efficiently leveraging GPUs.

## <a name="installation"></a> Installation

For standalone installation use Conda: `conda install conda-forge::segalign`

For installation in Galaxy we currently use the wrappers `richard-burhans:segalign` and `richard-burhans:batched_lastz` from the [Main Tool Shed](https://toolshed.g2.bx.psu.edu/).
Try the tools at usegalaxy.org: [segalign](https://usegalaxy.org/root?tool_id=toolshed.g2.bx.psu.edu/repos/richard-burhans/segalign/segalign/), [batched_lastz](https://usegalaxy.org/root?tool_id=toolshed.g2.bx.psu.edu/repos/richard-burhans/batched_lastz/batched_lastz/)

## <a name="cite_segalign"></a> Citing SegAlign

S. Goenka, Y. Turakhia, B. Paten and M. Horowitz,  "SegAlign: A Scalable GPU-Based Whole Genome Aligner," in 2020 SC20: International Conference for High Performance Computing, Networking, Storage and Analysis (SC), Atlanta, GA, US, 2020 pp. 540-552. doi: 10.1109/SC41405.2020.00043

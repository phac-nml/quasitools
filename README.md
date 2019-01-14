# Introduction

[![License](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)  [![Anaconda-Server Badge](https://anaconda.org/bioconda/quasitools/badges/installer/conda.svg)](https://anaconda.org/bioconda/quasitools)

Quasitools is a collection of tools for analysing viral quasispecies data. The following tools are currently available in quasitools:

* **aacoverage**: builds an amino acid consensus and returns its coverage
* **call**: contains tools for identifying variants between a NGS dataset a fasta reference sequence
	* **ntvar**: call nucleotide variants from a BAM file
	* **aavar**: call amino acid variants from a BAM file
	* **codonvar**: call codon variants from a BAM file
* **complexity**: reports the complexity of a quasispecies using several measures
* **consensus**: generate a consensus sequence from a BAM file
* **distance**: measures the distance between quasispecies using angular cosine distance
* **dnds**: calculates the dn/ds value for each region in a bed file
* **drmutations**: identifies amino acid mutations
* **hydra**: identify HIV drug resistance mutations in an NGS dataset
* **quality**: perform quality control on FASTQ reads

# Resources

* __Website__: https://phac-nml.github.io/quasitools/
* __Installation__: https://phac-nml.github.io/quasitools/installation/

# Release

__quasitools 0.6.0__

This release improves the documentation, integrates PyAAVF, and fixes a bug in quasitools distance that sometimes caused the software to crash.

# Installation

It is strongly recommended you refer to the documentation for full installation instructions. quasitools may be installed on any 64-bit Linux system from [Bioconda](https://bioconda.github.io/) with [Conda](https://conda.io/docs/) ([installation instructions](https://bioconda.github.io/#install-conda)):

1. Install [Bioconda](https://bioconda.github.io/)
2. Install the "quasitools" Bioconda package (conda install quasitools)

# Contact

**Eric Enns**: eric.enns@canada.ca

# Legal

Copyright Government of Canada 2017-2018

Written by: National Microbiology Laboratory, Public Health Agency of Canada

Licensed under the Apache License, Version 2.0 (the "License"); you may not use
this work except in compliance with the License. You may obtain a copy of the
License at:

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software distributed
under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR
CONDITIONS OF ANY KIND, either express or implied. See the License for the
specific language governing permissions and limitations under the License.



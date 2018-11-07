quasitools
==========

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

## Requirements and Dependencies

The Hydra pipeline requires:

* pysam >= 0.8.1
* SAMTools v1.3
* Bowtie 2 v2.2.6

## Installtion

### With [Conda](https://conda.io/docs/)

Install ``quasitools`` from [Bioconda](https://bioconda.github.io/) with [Conda](https://conda.io/docs/) ([Conda installation instructions](https://bioconda.github.io/#install-conda)):

```bash
# setup Conda channels for Bioconda and Conda-Forge (https://bioconda.github.io/#set-up-channels)
conda config --add channels defaults
conda config --add channels conda-forge
conda config --add channels bioconda
# install quasitools
conda install quasitools
```

## Legal

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

Contact
=======

**Eric Enns**: eric.enns@canada.ca  

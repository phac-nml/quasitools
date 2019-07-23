# Introduction

[![License](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)  [![Anaconda-Server Badge](https://anaconda.org/bioconda/quasitools/badges/installer/conda.svg)](https://anaconda.org/bioconda/quasitools)

Quasitools is a collection of tools for analysing viral quasispecies data. The following tools are currently available in quasitools:

* **[aacoverage](aacoverage.md)**: builds an amino acid consensus and returns its coverage
* **[call aavr](aavar.md)**: call amino acid mutations for a BAM file and a supplied reference file
* **[call codonvar](codonvar.md)**: call codon variants for a BAM file and a supplied reference file
* **[call ntvar](ntvar.md)**: call nucleotide variants for a BAM file and a supplied reference file
* **[complexity](complexity.md)**: reports the complexity of a quasispecies using several measures
* **[consensus](consensus.md)**: generate a consensus sequence from a BAM file
* **[distance](distance.md)**: measures the distance between quasispecies using angular cosine distance
* **[dnds](dnds.md)**: calculates the dn/ds value for each region in a bed file
* **[drmutations](drmutations.md)**: identifies amino acid mutations
* **[hydra](hydra.md)**: identify HIV drug resistance mutations in an NGS dataset
* **[quality](quality.md)**: perform quality control on FASTQ reads

# Release

__quasitools 0.7.0__

This release provides greater functionality for the complexity tool, including the ability to run on aligned read data and filter out low frequency haplotypes.

# Resources

* __Website__: [https://phac-nml.github.io/quasitools/](https://phac-nml.github.io/quasitools/)
* __Installation__: [https://phac-nml.github.io/quasitools/installation/](https://phac-nml.github.io/quasitools/installation/)
* __GitHub__: [https://github.com/phac-nml/quasitools](https://github.com/phac-nml/quasitools)

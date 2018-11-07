# quasitools

[![License](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)  [![Anaconda-Server Badge](https://anaconda.org/bioconda/quasitools/badges/installer/conda.svg)](https://anaconda.org/bioconda/quasitools)

Quasitools is a collection of tools for analysing viral quasispecies data. The following tools are currently available in quasitools:

* **aacoverage**: builds an amino acid consensus and returns its coverage
* **call**: contains tools for identifying variants between a NGS dataset and a fasta reference sequence
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

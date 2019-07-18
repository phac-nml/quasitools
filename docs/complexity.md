# complexity

Reports various measures of viral quasispecies complexity.

## Basic Usage

FASTA Input:

```bash
quasitools complexity fasta [OPTIONS] <FASTA READS FILE>
```

BAM Input:

```bash
quasitools complexity bam <FASTA REFERENCE FILE> <BAM FILE> <K-MER SIZE> [OPTIONS] 
```

## Arguments

### FASTA Reads File

This input file is only necessary when running the tool in FASTA mode.

An aligned FASTA file containing multiple aligned sequences, representing haplotypes of a genomic region from the mutant spectrum of interest. This would FASTA file would likely be created a multiple sequence alignment from aligning amplicon sequencing data.

### FASTA Reference File

This input file is only necessary when running the tool in BAM mode.

### BAM File

This input file is only necessary when running the tool in BAM mode.

A BAM file describing the alignments of reads to the same reference provided as input. These reads should be derived from a quasispecies mutant spectrum. This BAM file would likely be created using a read aligner which aligns FASTQ reads to a FASTA reference.

### K-Mer SIZE

This input is only necessary when running the tool in BAM mode.

The K-Mer size provides the sequence length for reads from a given starting position.


## Options

### FILTER

###### FLAG
```
-f [INTEGER]
```

This option is only available when running the tool in BAM mode.

Using this option allows for a user defined filter size between 0 and 100, haplotypes under the filter size will be removed from the positional list produced by running the command.

### OUTPUT FILE

###### FLAG
```
-o [USER-DEFINED-FILE-NAME.CSV]
```

This option is availble when running the tool in both BAM and FASTA mode.

Using this option allows users to define an output file in *CSV format* where the program output will be written.

## Output

The quasispecies complexity measures are taken from *Gregori, Josep, et al. 2016*. These include the following various indices:

Incidence (Entity Level):

* Number of haplotypes
* Number of polymorphic sites
* Number of mutations

Abundance (Molecular Level):

* Shannon entropy
* Simpson index
* Gini-Simpson index
* Hill numbers

Functional (Incidence):

* Minimum mutation frequency (Mf min)
* Mutation frequency (Mfe)
* FAD
* Sample nucleotide diversity

Functional (Abundance):

* Maximum mutation frequency (Mf max)
* Population nucleotide diversity

## Applications

* Assessing the quasispecies complexity of a genomic region.
* Comparing the quasispecies complexity of multiple genomic regions from the same mutant spectrum.

## Example: FASTA Reads File

### Data

The following example data may be used to run the tool:

* [aligned.fasta](data/aligned.fasta)

### Command

```bash
quasitools complexity fasta aligned.fasta
```

### Output

```text
Position,Number of Haplotypes,Haplotype Population,Number of Polymorphic Sites,Number of Mutations,Shannon Entropy,Shannon Entropy Normalized to N,Shannon Entropy Normalized to H,Simpson Index,Gini-Simpson Index,Hill Number #0,HIll Number #1,Hill Number #2,Hill Number #3,Minimum Mutation Frequency,Mutation Frequency,Functional Attribute Diversity,Sample Nucleotide Diversity (Entity),Maximum Mutation Frequency,Population Nucleotide Diversity,Sample Nucleotide Diversity
0,9,30,38,40,1.8774672554524843,0.5520018525167073,0.8544721713101401,0.19111111111111112,0.8088888888888889,9.0,6.536927510444632,5.232558139534883,4.543368996115371,0.013333333333333334,0.05555555555555555,7.379999999999999,0.10249999999999998,0.03866666666666667,0.06682222222222223,0.06912643678160921
```

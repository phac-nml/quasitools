# complexity

Reports various measures of viral quasispecies complexity.

## Basic Usage

```bash
quasitools complexity [options] <FASTA input>
```

## Arguments

### FASTA

An aligned FASTA file containing multiple aligned sequences, representing haplotypes of a genomic region from the mutant spectrum of interest. This would FASTA file would likely be created a multiple sequence alignment from aligning amplicon sequencing data.

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

## Example

### Data

The following example data may be used to run the tool:

* [aligned.fasta](data/aligned.fasta)

### Command

```bash
quasitools complexity aligned.fasta
```

### Output

```text
Starting...

Calculating the complexity from file: aligned.fasta


Incidence - Entity Level
-------------------------
Number of haplotypes (H) : 9
Number of polymorphic sites (P) : 38
Number of unique mutations (M) : 40


Abundance - Molecular Level
---------------------------
Shannon Entropy (Hs) : 1.87746725545
Shannon Entropy, normalized to log(N) (Hsn) : 0.552001852517
Shannon Entropy, normalized to log(H) (Hsh) : 0.85447217131
Simpson Index (Hsi) : 0.191111111111
Gini-Simpson Index (Hgs) : 0.808888888889
Hill numbers
  q = 0 : 9.0
  q = 1 : 6.53692751044
  q = 2 : 5.23255813953
  q = 3 : 4.54336899612


Functional, Indidence - Entity Level
-------------------------------------
Minimum Mutation Frequency (Mf min) : 0.0133333333333
Mutation Frequency (Mfe) : 0.0555555555556
Functional Attribute Diversity (FAD) : 7.379999999999999
Sample Nucleotide Diversity, Entity Level (^PIe) : 0.1025


Functional, Abundance - Molecular Level
----------------------------------------
Maximum Mutation Frequency (Mf max) : 0.03866666666666667
Population Nucleotide Diversity (PI) : 0.06682222222222223


Other
------
Sample Nucleotide Diversity (^PI) : 0.06912643678160921

Complete!

```

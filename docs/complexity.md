# Complexity

Reports various measures of viral quasispecies complexity. The input data is a single, aligned FASTA file, representing a genomic region from the mutant spectrum of interest.

## Basic Usage

```
quasitools complexity [options] <FASTA input>
```

## Output

The quasispecies complexity measures are taken from Gregori, Josep, et al. 2016. These include various indices:

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
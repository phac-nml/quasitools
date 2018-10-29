# Distance

Measures and reports the distance between multiple read pileups from multiple quasispecies samples, aligned to the same genomic region. This tool determines the cosine relatedness between viral quasispecies, reporting either angular cosine distance or cosine similarity as measures of relatedness. These measures of relatedness should be understood as approximations for evolutionary distance.

The software represents quasispecies pileup data as vectors and measures the cosine angle between every pair of vectors. This has the benefit of comparing the relative composition of the quasispecies and is robust against widely varying degrees of coverage. This method does not reduce the quasipecies data into a consensus vector and therefore can capture more nuanced differences between two quasispecies pileups. The software can also normalize the pileup vectors to prevent bias of relative coverage within a single pileup.

## Basic Usage

```
quasitools distance [options] <reference> (<BAM inputs>)+
```

## Output

A distance matrix with the distances between all pairs of quasispecies pileups will be written to standard output. The user may also have the matrix written to a file by specifying the location with the `-o` command.

## Applications

* Generating a distance matrix between inputs for the purpose of clustering.
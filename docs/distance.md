# Distance

Measures and reports the distance between multiple read pileups from multiple quasispecies samples, aligned to the same genomic region. This tool determines the cosine relatedness between viral quasispecies, reporting either angular cosine distance or cosine similarity as measures of relatedness. These measures of relatedness should be understood as approximations for evolutionary distance.

The software represents quasispecies pileup data as vectors and measures the cosine angle between every pair of vectors. This has the benefit of comparing the relative composition of the quasispecies and is robust against widely varying degrees of coverage. This method does not reduce the quasipecies data into a consensus vector and therefore can capture more nuanced differences between two quasispecies pileups. The software can also normalize the pileup vectors to prevent bias of relative coverage within a single pileup.

## Basic Usage

```
quasitools distance [options] <reference> (<BAM inputs>)+
```

## Options

### Normalize

```text
-n, --normalize / -dn, --dont_normalize
```

Whether or not to normalize read count data so that counts at each position sum to one. This is useful if one region of the reference has more coverage depth than another region. Normalization will ensure each base position will contribute equally to the distance measurement. Normalization is done by dividing base read counts (A, C, T, G) inside every 4-tuple by the sum of the read counts inside the same tuple.

### Output Distance

```text
-od, --output_distance / -os, --output_similarity
```

Either output an angular distance matrix (by default) or output a cosine similarity matrix. Please be aware that cosine similarity is not a metric.

### Start Position

```text
-s, --startpos INTEGER
```

Sets the start base position of the reference to use in the distance or similarity calculation. The start position is one-indexed [1, n].

### End Position

```text
-e, --endpos INTEGER
```

Sets the end base position of the reference to use in the distance or similarity calculation. The end position is one-indexed [1, n].

### Output

```text
-o, --output FILENAME
```

The file output location of the quasispecies distance or similarity matrix in CSV format.

### Truncate

```text
-t, --truncate
```

Ignores congiguous start and end pileup regions with no coverage.

### Remove No Coverage

```text
-r, --remove_no_coverage
```

Remove all regions of the pileup with no coverage from the distance calculations.

### Keep No Coverage

```text
-k, --keep_no_coverage
```

Do not remove regions of the pileup with no coverage from the distance calculations.

### Help

```text
--help
```

Show a help message and exit.


## Output

A distance matrix with the distances between all pairs of quasispecies pileups will be written to standard output. The user may also have the matrix written to a file by specifying the location with the `-o` command.

## Applications

* Generating a distance matrix between inputs for the purpose of clustering.

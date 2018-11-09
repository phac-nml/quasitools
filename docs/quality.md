# Quality

Performs quality control measures of FASTQ reads and outputs filtered FASTQ reads. The tool may filter reads based on average quality score, median quality score, mask low quality bases, perform iterative read trimming, and filter based on read length.

## Basic Usage

```
quasitools quality [options] <FASTQ forward> [<FASTQ reverse>] -o <output directory>
```

## Options

### Output Directory

```text
-o, --output_dir DIRECTORY
```

The location of the output directory.

### Trim Reads

```text
-tr, --trim_reads
```

If enabled, the software will iteratively trim the ends of reads until they either meet filter values or become to short. If trimmed reads become too short, they will be discarded.

If this option is disabled, the software will instead discard reads that do not meet filter values.

### Mask Reads

```text
-mr, --mask_reads
```

This option will mask low quality regions with "N" if they are below the minimum read quality threshold. This option and N-filtering (--ns) cannot be enabled simultaneously.

### Minimum Read Quality

```text
-rq, --min_read_qual INTEGER
```

Applies when read masking is enabled. The minimum quality for positions in the read before becoming masked with an "N".


### Length Cutoff

```text
-lc, --length_cutoff INTEGER
```

Filters out any reads shorter than this length.

### Score Cutoff

```text
-sc, --score_cutoff INTEGER
```

This software will filter out reads that have a median (default) or mean score less than this value.

### Median or Mean Score Filtering

```text
-me, --median / -mn, --mean
```

This determines whether the software will use median (default) or mean score as a cutoff for read filtering.

### Filtering Ns

```text
-n, --ns
```

If enabled, the software will discard any reads that contain N characters.

### Help

```text
--help
```

Displays the help message and exits.

## Output

The filtered reads will be writen to the output directory.

## Applications

* Read quality control before other performing other analyses.

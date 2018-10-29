# Quality

Performs quality control measures of FASTQ reads and outputs filtered FASTQ reads. The tool may filter reads based on average quality score, median quality score, mask low quality bases, perform iterative read trimming, and filter based on read length.

## Basic Usage

```
quasitools quality [options] <FASTQ forward> [<FASTQ reverse>] -o <output directory>
```

## Output

The filtered reads will be writen to the output directory.

## Applications

* Read quality control before other performing other analyses.
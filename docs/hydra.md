# HyDRA

Identify HIV drug resistance in a next generation sequencing dataset.

## Basic Usage

```
quasitools hydra [OPTIONS] <FORWARD READS> [REVERSE READS]
```

## Options

### Output

```text
-o, --output DIRECTORY
```

This is used to redirect from standard output to a file.

### Mutation Database

```text
-m, --mutation_db FILE
```

### Reporting Threshold

```text
-rt, --reporting_threshold INTEGER
```

Minimum mutation frequency percent to report.

### Generate Consensus

```text
-gc, --generate_consensus
```

Generate a mixed base consensus sequence.

### Consensus Percentage

```text
-cp, --consensus_pct INTEGER
```

The minimum percentage a base needs to be incorporated into the consensus sequence. This must be used with the `-gc/--generate_consensus` flag turned on.

### Quiet

```text
-q, --quiet
```

This is used to suppress all normal output throughout the pipeline. This does not affect any file generation.

### Trim Reads

```text
-tr, --trim_reads
```

If enabled, the pipeline will iteratively trim the ends of reads until they either meet filter values or become too short. If trimmed reads become too short, they will be discarded.
If not enabled, the pipeline will remove reads which do not meet filter values.

### Mask Reads

```text
-mr, --mask_reads
```

Mask low quality regions in reads with an 'N', if below the minimum read quality. This option and `-n/--ns` cannot be enabled simultaneously.

### Minimum Read Quality

```text
-rq, --min_read_qual INTEGER
```

The minimum quality that a position in a read must have. If below this threshold, the position will be masked as an 'N'.

### Length Cutoff

```text
-lc, --length_cutoff INTEGER
```

Reads which fall short of the specified length will be filtered out.

### Score Cutoff

```text
-sc, --score_cutoff INTEGER
```

Reads that have a median or mean quality score (depending on the score type specified) less than the score cutoff value will be filtered out.

### Mean or Median Score

```text
-me, --median / -mn, --mean
```

This determines whether the pipeline will use a mean or median (default) score as a cutoff for read filtering.

### Filter N's

```text
-n, --ns
```

If enabled, the pipeline will discard any reads that contain N characters.

### Error Rate

```text
-e, --error_rate FLOAT
```

Error rate for the sequencing platform.

### Minimum Variant Quality

```text
-vq, --minimum_variant_qual INTEGER
```

Minimum quality for variant to be included for processing in the pipeline.

### Minimum Depth

```text
-md, --min_dp INTEGER
```

Minimum required read depth for variant included for processing in the pipeline.

### Minimum Allele Count

```text
-ma, --min_ac INTEGER
```

The minimum required allele count for variant to be included for processing in the pipeline.

### Minimum Frequency

```text
-mf, --min_freq FLOAT
```

The minimum required frequency for mutation to be included for processing in the pipeline.

### ID

```text
-i, --id TEXT
```

This is used to specifiy a FASTA sequence identifier to be used in the consensus report.

## Output

The output directory location will default to the current directory and will be called `data/`. It will include the following output files:

* align.bam
* align.bam.bai
* consensus.fasta
* coverage_file.csv
* dr_report.csv
* filtered.fastq
* hydra.vcf
* mutation_report.aavf
* stats.txt


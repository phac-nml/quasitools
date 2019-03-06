# Quality

Performs quality control measures of FASTQ reads and outputs filtered FASTQ reads. The tool may filter reads based on average quality score, median quality score, mask low quality bases, perform iterative read trimming, and filter based on read length.

## Basic Usage

```bash
quasitools quality [options] <FASTQ forward> [<FASTQ reverse>] -o <output directory>
```

## Arguments

### FASTQ Forward

The forward or single-end FASTQ-format reads to be correted.

### FASTQ Reverse

The reverse reads, associated with paired-end data, to be corrected. This argument is not used if the tool is provided with only single-end reads.

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

## Example

This example masks all positions with quality scores less than 30 as Ns.

### Data

* [reads.fastq](data/reads.fastq)

### Command

```bash
quasitools quality -mr -rq 30 reads.fastq -o output
```

### Output

```text
@AF033819.3-1820
AACAAACTTGGCNATGAAAGCAACACTTTTTACAATANCAATTGGTACAAGCAGTNTTAGGCNGACNTCCTNGNTGCTTNNAGGGCTCTAGTCTAGNANC
+
CCCFFAFFHHHF)IGCFIFJIJIJIJJIIE@HBJGIJ(IIHCJHIBIJIJIBHDA/CJGI@H=DFE3EDE@'G;DCA?D0:AADDCDDCADDEDDC9E+D
@AF033819.3-1819
TAATAAGACGNTCAATGNAACAGGACCATGTACAAANGTNAGCACAGTANAATGTACACATGGNATTAGGCCAGTAGTANCANCTCNNCTNCTGTTAANN
+
?@CFFFBBHF2DDIJJG<JJJIJJEJJIGFJHIIJB=@G)GGBIJ@G?C6JJHJHIHC@JFJC-DC@EDDDHIEDDFFF:HE(DDD>>DB(BABCDDD>(
@AF033819.3-1818
TAATNCNGACGCTCTCGGANCCATCTCNCTCCTTCTNGCCNNCGCNAGTCAAAATTTTTGNCGTACTCACNAGTNNCCNCCNCTCGCCTCTTGCCGTGNG
+
@C?F1B;DHGGHGGJIIGI<JJIBGIG<JIJCGJIG3JJJ;<GHG*JJDFJJAG?GJGHI8CJE@HJADF2DFD(;CH9EA9DCDEEDCCBCDBDEFD>C
@AF033819.3-1817
CNGCCATTGTCANTATGTATTGTTTTTANTGGCCATNNTCCNGCTAATTTTAAAAGAAAATATGCTGTTTCCNNCCCTNTTTNTGCTGNAATNACTTCTN
+
C=?FDFDDHDHB0@JJIEJICDIHJDIF1IFJIGJJ):?IE>IFCHDIIJ@J@JCJJJ?JDEIJIFIHGJCJ95B?CD,@D?2DFBDE;DEC>DCEBBD3

(output truncated)
```

# Amino Acid Coverage  

Builds an amino acid census and returns its coverage.

## Basic Usage  

```bash
quasitools aacoverage [options] <BAM file> <reference file> <bed file>  
```

## Options

### BAM File

A [BAM file](https://samtools.github.io/hts-specs/SAMv1.pdf) (.bam) of sequences aligned to a related reference. A BAM index file (.bai) is also required and should be named the same as the BAM file, with the extension instead changed from ".bam" to ".bai".

### Reference File

A reference file related to the aligned BAM sequences. The provided reference file must be the same reference file used when producing the BAM and BAM index files.

### BED File

A [BED file](https://bedtools.readthedocs.io/en/latest/content/general-usage.html) that specifies the coordinates of genes, with repsect to the provided reference. This BED file must be a BED4+ file. That is, the BED file must contain at least the first 4 BED file columns.

### Output

```text
-o, --output FILENAME
```

The file output location of the amino acid coverage.

## Output  

The amino acid coverage will be output in CSV format. The output will have one entry per line, with the amino acid position and the coverage at that position. By default, this will be printed to standard output. The user may direct the output to a file by specifying a file name with the `-o/--output` option. 

## Example

### Data

The following example data may be used to run the tool:

* [hiv.fasta](data/hiv.fasta)
* [hiv.bed](data/hiv.bed)
* [variant.bam](data/variant.bam)
* [variant.bai](data/variant.bai)

### Command

```bash
quasitools aacoverage variant.bam hiv.fasta hiv.bed
```

### Output

```text
frame: 0
1,0
2,4
3,5
4,7
5,11
6,12
7,16
8,17
9,23
10,27
11,33
12,36
13,40
14,46
15,49
16,52
17,56
18,57
19,60
20,61
21,63
22,66
23,68

...
```


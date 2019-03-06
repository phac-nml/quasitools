# Amino Acid Coverage  

Builds an amino acid census and returns its coverage.

## Basic Usage  

```bash
quasitools aacoverage [options] <BAM file> <reference file> <bed file>  
```

## Arguments

### BAM File

A [BAM](https://samtools.github.io/hts-specs/SAMv1.pdf) file (.bam) of sequences aligned to a related reference. A BAM index file (.bai) is also required and should be named the same as the BAM file, with the extension instead changed from ".bam" to ".bai".

### Reference File

A reference file related to the aligned BAM sequences. The provided reference file must be the same reference file used when producing the BAM and BAM index files.

### BED File

A [BED file](https://bedtools.readthedocs.io/en/latest/content/general-usage.html) that specifies the coordinates of genes, with repsect to the provided reference. This BED file must be a BED4+ file. That is, the BED file must contain at least the first 4 BED file columns.

## Options

### Output

```text
-o, --output FILENAME
```

The file output location of the amino acid coverage.

## Output  

The amino acid coverage will be output in CSV format. The output will have one entry per line, with the amino acid position and the coverage at that position. By default, this will be printed to standard output. The user may direct the output to a file by specifying a file name with the `-o/--output` option.

## Applications

* Generating a report of the amino acid coverage, with respect to a referemce, from a BAM alignment file.

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
4,9
5,11
6,13
7,15
8,17
9,21
10,22
11,26
12,26
13,28
14,32
15,33
16,37
17,38
18,43
19,44
20,45
21,50
22,54
23,57

(output truncated)
```


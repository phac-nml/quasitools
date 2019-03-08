# Data Formats

This resource provides a detailed description of the various data formats used by quasitools, either as input or output, and their relationship to each other. In general, the data inputs usually must be consistent with each other. If you change one of the inputs, it is possible that some of the other inputs will need to change as well.

## Overview

The following is a summary of the formats used by quasitools:

| Format | Description |
| --- | --- |
| Reads | A FASTQ file containing sequencing reads. |
| Reference | A FASTA reference of either a gene, chromosome, or genome. |
| BAM | Specifies sequence alignments between the FASTA reference and FASTQ reads. |
| BAI | Indexes the BAM file for faster processing. |
| BED4 | Specifies the coordinates of coding sequences in the reference. |
| Mutation Database | Specifies meaningful amino acid mutations within coding sequences. |
| Codon Variants CSV | Specifies nucleotide variants within codons and resulting amino acid mutations. | 
| VCF | Specifies observed nucleotide variants and related information. |
| AAVF | Specifies observed amino acid variants and related information. |

## Relationships

The following is a summary of the relationships between quasitools inputs:

| Input | Input | Relationship |
| --- | --- | --- |
| BAM | Reads | The BAM file describes the alignment of reads to a reference. |
| BAM | Reference | The BAM file describes the alignment of reads to a reference. |
| BAM | BAI | The BAM index (BAI) file must be created from the BAM file. |
| BED4 | Reference | The BED4 file describes coding sequences on the reference. |
| BED | Mutation Database | The BED "name" column must be consistent with the MutationDB "coding sequences" column. |

## Reads

Most tools in quasitools do not operated directly on reads and instead require BAM alignment files, which must be created manually by the user. Please see the BAM format section for more information.

## Reference

The reference file must be in [FASTA format](https://en.wikipedia.org/wiki/FASTA_format). When providing a BAM file, BAM index file (BAI), and a reference file together to the same program within quasitools, the BAM file must have been generated from an alignment with the reference and the BAM index file generated from the BAM file.

## BAM

A [BAM](https://samtools.github.io/hts-specs/SAMv1.pdf) file (.bam) is the binary version of a [SAM](https://samtools.github.io/hts-specs/SAMv1.pdf) file. A SAM (Sequence Alignment/Map) file is a tab-delimited text file that specifies sequence alignments. SAM files are often converted into BAM files to reduce storage space and allow faster processing.

Within quasitools, BAM files are used to provide information about how sequencing reads aligned to a reference file. A BAM index file (.bai) is also required by quasitools and should be named the same as the BAM file, with the extension instead changed from ".bam" to ".bai". When running quasitools, it is important to ensure that the BAM file, BAM index file, and reference file are consistent with each other. The BAM file should be generated from an alignment against the reference and the BAM index file should be generated from the BAM file.

BAM alignment files may be generated using read or sequence alignment software, such as [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml). If your sequence alignment software does not output a BAI file, you may use [samtools](http://www.htslib.org/doc/) to index your BAM file.

## BED4

A [BED file](https://bedtools.readthedocs.io/en/latest/content/general-usage.html) specifies the coordinates of coding sequences, with repsect to sequences within a FASTA reference. The FASTA sequences are larger sequences, such a chromosomes or genes. The specified coding sequences are smaller sequences within the FASTA sequences. These small sequences within FASTA sequences may be genes or gene products. However, the should not contain non-coding sequences. The coordinates of the coding sequences must be specified in 0-based nucleotide coordinates and the length of each coding sequence must be divisible by three.

The BED files used with quasitools must be BED4+ files and therefore must contain at least the first 4 BED file columns. It is very likely that these BED4 files will need to be manually created with knowledge of the locations of coding sequences in the reference. It is very important that the ```names``` of these coding sequences are the same names used in the ```coding sequences``` column of the mutation database, if the user is providing both files as inputs to a tool. The format of the BED files expected by quasitools is as follows:

```
[identifier] [start] [end] [name]
```

There must only be one record per line in the BED4 file. The items in each record are specified as follows:

| Item | Description | Values |
| --- | --- | --- |
| identifier | The sequence identifier of the FASTA sequence within the FASTA reference file. | string |
| start | The starting nucleotide coordinate of the coding sequence. | integer |
| end | The ending nucleotide coordinate of the coding sequence. | integer |
| name | The name of the coding sequence. | string |

### Example: Genes Witin a Chromsome 

**Reference File (FASTA)**

```text
>chromosome1
ACGTACGT ...
>chromesome2
GGAATTCC ...
```

**BED4 File**

```text
chromosome1    300    599     gene1
chromosome1    900    1199    gene2
chromosome2    300    599     gene3
```

Observe that ```chromosome1``` and ```chromosome2``` are the names of contigs in a reference FASTA file. Our BED4 file specifies three genes within the reference file. ```gene1``` and ```gene2``` reside on the ```chromosome1``` contig and ```gene3``` resides on the ```chromosome2``` contig.

### Example: Gene Products Witin a Gene

**Reference File (FASTA)**

```text
>gene1
ACGTACGT ...
```

**BED4 File**

```text
gene1    300    599     product1
gene1    900    1199    product2
```

Observe that our reference file specifics one contig, ```gene1```, which corresponds to a gene. The coding sequences specific in the BED4 file are instead gene products: ```product1``` and ```product2```.

### Example: Gene Products Witin a Chromosome

**Reference File (FASTA)**

```text
>chromosome1
ACGTACGT ...
>chromesome2
GGAATTCC ...
```

**BED4 File**

```text
chromosome1    300    399    product1
chromosome1    900    999    product2
chromosome2    300    399    product3
```

Observe that ```chromosome1``` and ```chromosome2``` are the names of contigs in a reference FASTA file. Our BED4 file specifies three coding sequences which are gene products (assumed to be within genes). ```product1``` and ```product2``` reside on ```chromosome1``` and ```product3``` resides on ```chromosome2```.

## Mutation Database

The mutation database describes specific mutations within the named coding sequences specified previously in the BED4 file. It is used to report animo acid mutations that the user has identified. The mutation database is represented as a specifically formated TSV (tab separated values) file. The mutation database format is as follows:

```text
[coding sequence] [wildtype] [position] [mutation] [category] [surveillance] [comment]
```

| Item | Description | Values |
| --- | --- | --- |
| coding sequence | The name of the coding sequence (CDS or gene product). | string |
| wildtype   | The amino acid of wildtype. | character |
| position | The position of the amino acid within the coding sequence. | integer |
| mutation | The amino acid of the muation. | character |
| category | A categorization of the mutation. | string |
| surveillance | Whether or not the mutation is part of a surveillance program. | string |
| comment | A user-provided comment about the mutation. | string |

In this format, there is one mutation entry per line in the file. Each item of every entry is separated by tabs. A line may be a comment if the very first character of the line is a ```#``` character. Please note that the ```position``` is the amino acid position within that particular coding sequence, in amino acid coordinates. It is very important that the values of the "coding sequence" column are the same names used in the "name" column of the BED4 file, if the user is providing both files as inputs to a tool. Every name that appears in the ```coding sequence``` column must also appear in the ```name``` column of the BED4 file. However, not every ```name``` in the BED4 file must appear under the ```coding sequence``` column in the mutation database.

It is very likely that these mutation database files will need to be manually created with knowledge of the locations of mutations within the coding sequences specified in the BED file.

### Example: Mutations Within a Gene

**Mutation Database**

```text
#genetic region wildtype	position	mutation	category	surveillance	comment
gene1	M	10	F	major	Yes	comment1
gene1	K	20	I	minor	No	comment2
gene2	W	5	V	major	No	comment3
```

**BED4 File**

```text
chromosome1    300    599     gene1
chromosome1    900    1199    gene2
chromosome2    300    599     gene3
```

Observe that the all names in the first column of the mutation database (```gene1```, ```gene2```) appear in the last column of the BED4 file. However, not every named coding sequence in the BED4 (```gene3```) has to appear in the mutation database ```coding sequence``` column.

## Codon Variants CSV

The codon variants CSV file is consistent with the standard CSV format. The file describes nucleotide variants within codons and their resulting amino acid variants. It also clarifies whether the nucleotide variants are synonymous or a non-synonymous mutations. The codon variants CSV file has the following comment header:

```text
#gene,nt position (gene),nt start position,nt end position,ref codon,mutant codon,ref AA,mutant AA,coverage,mutant frequency,mutant type,NS count,S count
```

The items on each line correspond to the following information:

| Position | Name | Description |
| --- | --- | --- |
| 0 | gene | The name of the coding region. |
| 1 | nt position (gene) | The start and end nucleotide positions of the gene which contains the codon. |
| 2 | nt start position | The start nucleotide position of the codon. |
| 3 | nt end position | The end nucleotide position of the codon. |
| 4 | ref codon | The codon in the reference. |
| 5 | mutant codon | The mutant codon in the data. |
| 6 | ref AA | The corresponding amino acid in the reference. |
| 7 | mutant AA | The corresponding amino acid in the mutant. |
| 8 | coverage | The coverage of the codon. |
| 9 | mutant frequency | The observed frequency of the mutant codon. |
| 10 | mutant type | Whether or not the mutation is synonymous (S) or nonsynonymous (NS). |
| 11 | NS count | The expected number of nonsynonymous sites in the codon. A number between [0, 3]. |
| 12 | S count | The expected number of synonymous sites in the codon. A number between [0, 3]. |

Please refer to [Nei and Gojobori 1986](https://academic.oup.com/mbe/article/3/5/418/988012) for more information about how to calculate ```NS count``` and ```S count```.

## VCF

The VCF format used within quasitools is consistent with the standard VCF format. However, quasitools uses custom values in the ```FILTER``` and ```INFO``` columns.

### FILTER

| Name | Meaning |
| --- | --- |
| dp100 | This variant was filtered because the coverage depth was less than 100. |
| q30 | This variant was filtered because the quality of the variant was less than 30. |
| ac5 | This variant was filtered because the variant was observed less than 5 times. |

Of particular interest is the ```q30``` flag. This flag is set when the estimated quality of the variant is less than 30. quasitools calculates the probability of a variant being legitimate using the Poisson cumulative distribution function. In this framework, λ is the expected number of errors at a particular position (the coverage depth of that position multiplied by the error rate).

In order for the variant to accepted as a real mutation, the probability of the observed variant being caused entirely by errors must be sufficiently low. In other words, at least some of the observed variant was probably caused by at least one real mutation. In order for the probability to be sufficiently low, it must be less than Q30 (1 in 1000 chance). When performing this probability calculation, we assume the worst case scenario: all expected substution errors at a particular position are the same nucleotide as the variant being tested, rather than being evenly distributed evenly across all possible substitions.

### INFO

| Name | Meaning |
| --- | --- |
| DP | The total coverage depth of the pileup at this position. |
| AC | The number of times this particular variants was observed in the pileup at this position. |
| AF | The frequency of this particular variants was observed in the pileup at this position. |

## AAVF

AAVF is a text file format, inspired by the Variant Call Format (VCF). It contains meta-information lines, a head line, and then data lines each containing information about a position in a gene within a genome. Please refer to the [AAVF documentation](https://github.com/winhiv/aavf-spec/raw/master/AAVFv1.0.pdf) for more information.


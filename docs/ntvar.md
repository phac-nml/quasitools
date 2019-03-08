# ntvar

Calls nucleotide variants observed witin an aligned BAM file when compared against a supplied reference file.

## Basic Usage  

```bash
quasitools ntvar [options] <BAM file> <reference file>
```

## Arguments

### BAM File

A [BAM](https://samtools.github.io/hts-specs/SAMv1.pdf) file (.bam) of sequences aligned to a related reference. A BAM index file (.bai) is also required and should be named the same as the BAM file, with the extension instead changed from ".bam" to ".bai".

### Reference File

A reference file related to the aligned BAM sequences. The provided reference file must be the same reference file used when producing the BAM and BAM index files.

## Options

### Error Rate

```text
-e, --error_rate FLOAT
```

This is the expected substitution sequencing error rate. The default value is 0.0021 substitutions per sequenced base.

### Output

```text
-o, --output FILENAME
```

The file output location to write the identified nucleotide variants in VCF format. Otherwise, the results are printed to standard output.

## Output

The output of the tool is a list of the variants observed witin the aligned BAM file when compared against the supplied reference file. The output is in VCF format and the quasitools provides the following additional information in the ```FILTER``` and ```INFO``` columns of the VCF file.

Please see [Data Formats](../formats) for more information about the custom fields used by quasitools. However, a short description is provided below.

### FILTER

| Name | Meaning |
| --- | --- |
| dp100 | This variant was filtered because the coverage depth was less than 100. |
| q30 | This variant was filtered because the quality of the variant was less than 30. |
| ac5 | This variant was filtered because the variant was observed less than 5 times. |

### INFO

| Name | Meaning |
| --- | --- |
| DP | The total coverage depth of the pileup at this position. |
| AC | The number of times this particular variants was observed in the pileup at this position. |
| AF | The frequency of this particular variants was observed in the pileup at this position. |

## Example

### Data

The following example data may be used to run the tool:

* [hiv.fasta](data/hiv.fasta)
* [variant.bam](data/variant.bam)
* [variant.bai](data/variant.bai)

### Command

```bash
quasitools call ntvar variant.bam hiv.fasta
```

### Output

```text
##fileformat=VCFv4.2
##fileDate=20190206
##source=quasitools
##contig=<ID=AF033819.3,length=9181>
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele Count">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
##FILTER=<ID=dp100,Description="Set if True; DP<100">
##FILTER=<ID=q30,Description="Set if True; QUAL<30">
##FILTER=<ID=ac5,Description="Set if True; AC<5">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
AF033819.3	153	.	t	a	100	PASS	DP=129;AC=129;AF=1.0000
AF033819.3	342	.	g	c	100	PASS	DP=141;AC=141;AF=1.0000
AF033819.3	719	.	c	g	100	PASS	DP=132;AC=132;AF=1.0000
AF033819.3	1025	.	a	t	100	PASS	DP=129;AC=129;AF=1.0000
AF033819.3	1351	.	c	a	100	PASS	DP=135;AC=135;AF=1.0000
AF033819.3	1657	.	a	c	100	PASS	DP=111;AC=111;AF=1.0000
AF033819.3	1917	.	g	a	100	PASS	DP=128;AC=128;AF=1.0000
AF033819.3	2052	.	t	g	100	PASS	DP=147;AC=147;AF=1.0000
AF033819.3	2368	.	a	t	100	PASS	DP=140;AC=140;AF=1.0000
AF033819.3	2422	.	g	c	100	PASS	DP=146;AC=146;AF=1.0000
AF033819.3	2989	.	a	g	100	PASS	DP=108;AC=108;AF=1.0000
AF033819.3	5707	.	t	a	100	PASS	DP=119;AC=119;AF=1.0000
AF033819.3	5970	.	a	c	100	PASS	DP=117;AC=117;AF=1.0000
AF033819.3	6139	.	c	g	100	PASS	DP=138;AC=138;AF=1.0000
AF033819.3	6674	.	a	t	100	PASS	DP=142;AC=142;AF=1.0000
AF033819.3	7366	.	a	c	100	PASS	DP=129;AC=129;AF=1.0000
AF033819.3	8631	.	c	t	100	PASS	DP=144;AC=144;AF=1.0000
```

# codonvar  

Call codon variants for a given [BAM](https://samtools.github.io/hts-specs/SAMv1.pdf) alignment file. A report is generated that details nucleotide variants within codons and the resulting amino acid variants. The report indicates whether the nucleotide variants correspond to a synonymous or a non-synonymous mutation.

## Basic Usage

```bash
quasitools codonvar [options] <BAM file> <reference file> <offset> <bed file>
```

## Arguments

### BAM File

A [BAM](https://samtools.github.io/hts-specs/SAMv1.pdf) file (.bam) of sequences aligned to a related reference. A BAM index file (.bai) is also required and should be named the same as the BAM file, with the extension instead changed from ".bam" to ".bai".

### Reference File

A reference file related to the aligned BAM sequences. The provided reference file must be the same reference file used when producing the BAM and BAM index files.

### Offset

An integer to offset the reported positions in the output. This does not change the frame or coordinates of coding sequences in the BED file.

It may be useful to provide an offset when providing a gene as a reference and gene products in the BED file, but want to report codon variants with respect to the entire chromosome. In this circumstance, the offset would be the position of the reference gene, with respect to the chromosome on which it resides.

### BED File

A [BED file](https://bedtools.readthedocs.io/en/latest/content/general-usage.html) that specifies the coordinates of genes, with repsect to the provided reference. This BED file must be a BED4+ file and therefore contain at least the first 4 BED file columns. The "names" of these genetic regions in the BED4 file must be the same names used in the "genetic regions" column of the mutation database. Please refer to [Data Formats](../formats) for more information.

### Mutation Database

The mutation database describes specific mutations within the named genetic regions specified previously in the BED4 file. The entries in the "genetic regions" colummn of this database must match the "names" column of the provided BED4 file. Please refer to [Data Formats](../formats) for more information.

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

The file output location to write the identified amino acid mutations.
## Example

### Data

* [variant.bam](data/variant.bam)
* [variant.bai](data/variant.bai)
* [hiv.fasta](data/hiv.fasta)
* [hiv.bed](data/hiv.bed)

### Command

```bash
quasitools call codonvar variant.bam hiv.fasta 0 hiv.bed
```

### Output

```text
#gene,nt position (gene),nt start position,nt end position,ref codon,mutant codon,ref AA,mutant AA,coverage,mutant frequency,mutant type,NS count,S count
vpu,5607-5855,5706,5708,tta,Ata,L,I,118,100.00,NS,1.0000,0.0000
nef,8342-8713,8630,8632,ctg,Ttg,L,L,143,100.00,S,0.0000,1.0000
pol,1630-4641,2050,2052,att,aGt,I,S,145,100.00,NS,1.0000,0.0000
pol,1630-4641,2986,2988,gaa,gaG,E,E,107,100.00,S,0.0000,1.0000
pol,1630-4641,2419,2421,ctg,ctC,L,L,143,100.00,S,0.0000,1.0000
pol,1630-4641,1654,1656,cta,ctC,L,L,110,100.00,S,0.0000,1.0000
pol,1630-4641,1915,1917,gga,gAa,G,E,128,100.00,NS,1.0000,0.0000
pol,1630-4641,2365,2367,caa,caT,Q,H,139,100.00,NS,1.0000,0.0000
env,5770-8340,5968,5970,aat,aCt,N,T,116,100.00,NS,1.0000,0.0000
env,5770-8340,6136,6138,acc,acG,T,T,134,100.00,S,0.0000,1.0000
env,5770-8340,7363,7365,gca,gcC,A,A,128,100.00,S,0.0000,1.0000
env,5770-8340,6673,6675,aat,Tat,N,Y,138,100.00,NS,1.0000,0.0000
gag,335-1837,1655,1657,tac,tCc,Y,S,109,100.00,NS,1.0000,0.0000
gag,335-1837,716,718,gtc,gtG,V,V,131,100.00,S,0.0000,1.0000
gag,335-1837,341,343,gcg,Ccg,A,P,140,100.00,NS,1.0000,0.0000
gag,335-1837,1022,1024,gaa,gaT,E,D,126,100.00,NS,1.0000,0.0000
gag,335-1837,1349,1351,cca,cAa,P,Q,133,100.00,NS,1.0000,0.0000
```

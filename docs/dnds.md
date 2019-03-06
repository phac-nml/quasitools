# dNdS

Determines the dN/dS ratio for each codon variant in a supplied csv file (codon variants). The dN/dS ratio is the ratio between the number of nonsynonymous substitutions per non-synonymous site to the number of synonymous substitutions per synonymous site. 

The reported ratio is calculated by observing by first observing what proportion of nonsynonymous and synonymous sites have mutations, pN and pS. However, pN and pS are calculated with consideration for many reads are covering the same mutation sites. pN and pS are then converted to dN and dS using the Jukes-Cantor correction, which estimates the expected proportion of changes between two sequences (unobservable), based on our observed proportions (pN and pS). The corrected Jukes-Cantor rates will be higher, because the expectation is that over time, multiple mutations will have happened at the same site and this information is lost in our observed proportions; a slightly higher mutation rate (dN and dS) would have been required to observe our proportions (pN and pS).

This dN/dS ratio can be used as an indicator of evolutionary pressure acting on a codon. For more information about the dN/dS calculation and how to interpret the ratio, please view [Morelli et al. 2013](https://veterinaryresearch.biomedcentral.com/articles/10.1186/1297-9716-44-12) and [Nei and Gojobori 1986](https://academic.oup.com/mbe/article/3/5/418/988012).

## Basic Usage 

```bash
quasitools dnds [options] <codon variants csv> <reference file> <offset>
```

## Arguments

### Codon Variants CSV

The codon variants CSV file should be taken directly from the output of [call codonvar](../codonvar) and provided as input when running this tool. For more information about this format, please refer to [Data Formats](../formats).

### Reference File

A reference file to compare the codon variants against. The provided reference file must be the same reference file used when producing the codon variants CSV file.

### Offset

An integer to offset the reported positions in the output. This does not change the frame or coordinates of coding sequences in the BED file.

It may be useful to provide an offset when providing a gene as a reference and gene products in the BED file, but want to report codon variants with respect to the entire chromosome. In this circumstance, the offset would be the position of the reference gene, with respect to the chromosome on which it resides.

## Options

### Output

```text
-o, --output FILENAME
```

The file output location to write the dN/dS ratios.

## Example

### Data

* [hiv.fasta](data/hiv.fasta)
* [variant_codonvar.csv](data/variant_codonvar.csv)


### Command

```bash
quasitools dnds variant_codonvar.csv hiv.fasta 0
```

### Output

```text
#gene,pn,ps,pn_sites,ps_sites,dn/ds
pol,0.4345,1.5000,3,3,N/A
env,0.3750,1.0000,2,2,N/A
gag,0.4375,1.0000,4,1,N/A
```

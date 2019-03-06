# Drug Resistance Mutations

Generates a report detailing the drug resistant mutations found. This tool identifies which nucleotide mutations, identified in the VCF file, have resulted in a noteworthy amino acid mutation, as specified by the mutation database. Please refer to [Data Formats](../formats) for detailed information about the the expected input formats for this tool.

## Basic Usage

```bash
quasitools drmutations [options] <BAM file> <reference file> <variants file> <bed file> <mutation db>
```

## Arguments

### BAM File

A [BAM](https://samtools.github.io/hts-specs/SAMv1.pdf) file (.bam) of sequences aligned to a related reference. A BAM index file (.bai) is also required and should be named the same as the BAM file, with the extension instead changed from ".bam" to ".bai".

### Reference File

A reference file related to the aligned sequences in the BAM file. The provided reference file must be the same reference file used when producing the BAM and BAM index files.

### Variants File

A VCF-formated file containing information about nucleotide variants in the reads with respect to an aligned to the provided reference. The read alignment information is taken from the provided BAM file. Only mutations which ```PASS``` the VCF filters will be considered when running the tool. Additionally, a VCF file may be created using [quasitools call ntvar](../ntvar).

### BED File

A [BED file](https://bedtools.readthedocs.io/en/latest/content/general-usage.html) that specifies the coordinates of genes, with repsect to the provided reference. This BED file must be a BED4+ file and therefore contain at least the first 4 BED file columns. The "names" of these genetic regions in the BED4 file must be the same names used in the "genetic regions" column of the mutation database. Please refer to [Data Formats](../formats) for more information.

### Mutation Database

The mutation database describes specific mutations within the named genetic regions specified previously in the BED4 file. When provided to the tool. The entries in the "genetic regions" colummn of this database must match the "names" column of the provided BED4 file. Please refer to [Data Formats](../formats) for more information.

## Options  

### Reporting Threshold

```text
-t, --reporting_threshold INTEGER  
```

The minimum number of observations required in the read data (BAM file) for an entry to be reported in the drug resistance report. Mutations with a number of observations less than this will not be reported. The default value is 1. 

### Output

```text
-o, --output FILENAME
```

The file output location to write the identified drug resistant mutations.


## Example

### Data

* [variant.bam](data/variant.bam)
* [variant.bai](data/variant.bai)
* [hiv.fasta](data/hiv.fasta)
* [variant.vcf](data/variant.vcf)
* [hiv.bed](data/hiv.bed)
* [hiv_db.tsv](data/hiv_db.tsv)

### Command

```bash
quasitools drmutations variant.bam hiv.fasta variant.vcf hiv.bed hiv_db.tsv
```

### Output

```text
Chromosome,Gene,Category,Surveillance,Wildtype,Position,Mutation,Mutation Frequency,Coverage
AF033819.3,gag,major,Yes,A,3,P,100.00,140
AF033819.3,gag,minor,No,Y,441,S,100.00,109
AF033819.3,pol,major,Yes,G,96,E,100.00,128
AF033819.3,pol,minor,No,Q,246,H,100.00,139
AF033819.3,vpu,major,Yes,L,34,I,100.00,118
AF033819.3,env,minor,No,N,67,T,100.00,116
```


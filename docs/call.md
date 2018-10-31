# Call

# ntvar

Call nucleotide variants from a BAM file and a supplied reference file.

## Basic Usage  

```
quasitools ntvar [options] <BAM file> <reference file>
```

---

# aavar  

Identifies amino acid mutations.

## Basic Usage  

```
quasitools aavar [options] <BAM file> <reference file> <variants file> <genes file> <mutation db>
```

---

# codonvar  

Call codon variants for a given BAM file. A report is generated that details nucleotide variants within a codon and the resulting AA variants. The report indicates whether the nucleotide variants correspond to a synonymous or a non-synonymous mutation.  

## Basic Usage

```
quasitools codonvar [options] <BAM file> <reference file> <offset> <genes file>
```

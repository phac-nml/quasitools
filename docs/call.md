# Call

Contains tools for identifying variants between a NGS dataset and a fasta reference sequence.

## Basic Usage  

```
quasitools call SUBCOMMAND
```

## Subcommands

## ntvar

Call nucleotide variants from a BAM file and a supplied reference file.

### Basic Usage  

```
quasitools ntvar [options] <BAM file> <reference file>
```

### Options  

**Error Rate**
```text
-e, --error_rate FLOAT
```

This is the expected sequencing error rate. Defaults to 0.0021.

**Output**
```text
-o, --output FILENAME
```

This is used to redirect from standard output to a file.  

---

## aavar  

Identifies amino acid mutations for a BAM file.

### Basic Usage  

```
quasitools aavar [options] <BAM file> <reference file> <variants file> <genes file> <mutation db>
```

### Output  

**Minimum Frequency**  

```text
-f, --min_freq FLOAT
```

The minimum required frequency for variant to be included. Defaults to 0.01

**Error Rate**  

```text
-e, --error_rate FLOAT
```

This is the expected sequencing error rate. Defaults to 0.0021.  

**Output**  

```text
-o, --output FILENAME
```  

This is used to redirect from standard output to a file.  

---

## codonvar  

Call codon variants for a given BAM file. A report is generated that details nucleotide variants within a codon and the resulting AA variants. The report indicates whether the nucleotide variants correspond to a synonymous or a non-synonymous mutation.  

### Basic Usage

```
quasitools codonvar [options] <BAM file> <reference file> <offset> <genes file>
```

### Options  

**Error Rate**  

```text
-e, --error_rate FLOAT
```

This is the expected sequencing error rate. Defaults to 0.0021.  

**Output**  

```text
-o, --output FILENAME
```

This is used to redirect from standard output to a file.    

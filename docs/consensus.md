# Consensus  

Generate a consensus sequence from a BAM file.

## Basic Usage  

```
quasitools consensus [options] <BAM file> <reference file>
```

## Options  

### Percentage  

```text
-p, --percentage INTEGER
```

This is the percentage threshold for the incorporation of an 'n' base.  
When percentage is set to 100, the most frequent base will be incorporated (note: in the case of a tie, the base will be chosen in reverse alphabetical order). Insertions that are at least a multiple of 3 will be incorporated (i.e. codon length). When percentage is less than 100, the base that has frequency greater than or equal to percentage set will be incorporated, with tie-breaking the same as above. When there is zero coverage, or no bases meet the percentage threshold (only when percentage is < 100) a 'n' will be incorporated.  

### ID  

```text
-i, --id TEXT
```

Specify the default FASTA sequence identifier to be used for sequences without an RG tag.

### Output

```text
-o, --ouput FILENAME
```

This is used to redirect from standard output to a file.

## Output  

By default, the consensus sequence will print to standard out. The user may direct the output to a file by specifying a file name with the `-o/--output` option.

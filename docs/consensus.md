# Consensus  

Generate a consensus sequence from a BAM file.

## Basic Usage  

```
quasitools consensus [options] <BAM file> <reference file>
```

Options: 

* -p, --percentage [INTEGER]:  
When percentage is set to 100, the most frequent base will be incorporated (note: in the case of a tie, the base will be chosen in reverse alphabetical order). Insertions that are at least a multiple of 3 will be incorporated (i.e. codon length). When percentage is less than 100, the base that has frequency greater than or equal to percentage set will be incorporated, with tie-breaking the same as above. When there is zero coverage, or no bases meet the percentage threshold (only when percentage is < 100) a 'n' will be incorporated.  

* -i, --id [TEXT]:  
Specify the default FASTA sequence identifier to be used for sequences without an RG tag.  

* -o, --output [FILENAME]:  
Specify the name of an output file.

# Amino Acid Coverage  

Builds an amino acid census and returns its coverage.

## Basic Usage  

```
quasitools aacoverage [options] <BAM file> <reference file> <genes file>  
```

## Options

### Output  

```text
-o, --output FILENAME
```

## Output  

The output will have one entry per line, with the AA position and the coverage at that position. By default, this will be printed to standard output. The user may direct the output to a file by specifying a file name with the `-o/--output` option. 


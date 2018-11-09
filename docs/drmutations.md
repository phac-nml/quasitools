# Drug Resistance Mutations

Generates a report detailing the drug resistant mutations found.  

## Basic Usage

```
quasitools drmutations [options] <BAM file> <reference file> <variants file> <genes file> <mutation db file>
```

## Options  

```text
-f, --min_freq FLOAT  
```

The minimum required frequency. Defaults to 0.01.  

```text
-t, --reporting_threshold INTEGER  
```

The minimum percentage required for an entry in the drug resistance report. Default: 1. 

```text
-o, --output FILENAME
```

This is used to direct from standard output to a file.



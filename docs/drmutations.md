# Drug Resistance Mutations

Generates a report detailing the drug resistant mutations found, above the reporting threshold (default: 1%).  

## Basic Usage

```
quasitools drmutations [options] <BAM file> <reference file> <variants file> <genes file> <mutation db file>
```

Options:  

* -f, --min_freq [FLOAT]  
The minimum required frequency. Default: 0.01.  

* -t, --reporting_threshold [INTEGER]  
The minimum percentage required for an entry in the drug resistance report. Default: 1. 

* -o, --output [FILENAME]  


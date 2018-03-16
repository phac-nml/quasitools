# Change Log

All notable changes to Quasitools will be documented in this file.

## 0.2.3 ##

2018-03-16

### Changed ###

- Fixed FASTA identifier in consensus output for the consensus and hydra commands
- User can pass in a default id as an option
- RG tag is used if present in the supplied .bam file

## 0.2.2 ##

2017-11-27

### Changed ###

 - Updated samtools dependency to v1.3

## 0.2.1 ##

2017-11-22

### Changed ###

 - fixed call codonvar command to properly output to either file or stdout depending on flag

## 0.2.0 ##

2017-11-03

### Added ###

 - commands:
   - hydra: which identifies HIV Drug Resistance in a NGS sample
   - call: new base command for nt, aa, and codon variant calling
     - ntvar: formally just callntvar
     - aavar: formally aavariants
     - codonvar: reports codon variants for use with dnds command
   - dnds: calculates dnds values per coding region using codonvar output

### Changed ###

 - Renamed aa_census command to aa_coverage
 - Hardended command options and parameters

## 0.1.0 ##

2017-09-06

This is the initial release of Quasitools

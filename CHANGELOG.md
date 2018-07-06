# Change Log

All notable changes to Quasitools will be documented in this file.

## 0.4.0 ##

### Added ###

 - commands:
  - quality: performs quality control on FASTQ reads
 - modules:
  - quality_control.py: contains QualityControl class

### Changed ###

- Updated hydra command to use new quality control class, masking and/or trimming reads based on command-line options that user has specified.
- Fixed bug in cmd_consensus.py which occured when there are multiple records in the reference file and the user is writing to a file.

## 0.3.1 ##

2018-04-13

### Added ###

- Updated hydra command to accept an id to be used as the sequence identifier
  in the consensus report, and as the RG-id in the bam file

## 0.3.0 ##

2018-04-10

### Added ###

 - commands:
   - distance: returns the evolutionary distances between viral quasispecies as a distance matrix
 - modules:
   - pileup.py: contains Pileup and Pileup_List classes
   - distance.py: contains Distance_Matrix class

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

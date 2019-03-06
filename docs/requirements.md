# Requirements

If quasitools is installed through Conda or Galaxy, then dependencies should be installed automatically. However, several tools within quasitools require BAM and BAM index files as input. These will need to be created using alignment software, such as [Bowtie 2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) to generate BAM alignment files and [SAMTools](http://www.htslib.org/) to generate BAM index files. Please see the [Data Formats](../formats) resource for more information about how to prepare input to quasitools.

## HyDRA Pipeline

When running quasitools directly from source, the HyDRA pipeline requires the follow dependencies:

* [pysam](https://pysam.readthedocs.io/en/latest/api.html) >= 0.8.1
* [SAMTools](http://www.htslib.org/) v1.3
* [Bowtie 2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) v2.2.6


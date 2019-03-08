# aavar  

Call amino acid mutations from a [BAM](https://samtools.github.io/hts-specs/SAMv1.pdf) alignment file and a supplied reference file. Please refer to [Data Formats](../formats) for detailed information about the the expected input formats for this tool.

## Basic Usage  

```bash
quasitools aavar [options] <BAM file> <reference file> <bed file> [variants file] [mutation db]
```

## Arguments

### BAM File

A [BAM](https://samtools.github.io/hts-specs/SAMv1.pdf) file (.bam) of sequences aligned to a related reference. A BAM index file (.bai) is also required and should be named the same as the BAM file, with the extension instead changed from ".bam" to ".bai".

### Reference File

A reference file related to the aligned sequences in the BAM file. The provided reference file must be the same reference file used when producing the BAM and BAM index files.

### BED File

A [BED file](https://bedtools.readthedocs.io/en/latest/content/general-usage.html) that specifies the coordinates of genes, with repsect to the provided reference. This BED file must be a BED4+ file and therefore contain at least the first 4 BED file columns. The "names" of these genetic regions in the BED4 file must be the same names used in the "genetic regions" column of the mutation database. Please refer to [Data Formats](../formats) for more information.

### Variants File

A [VCF (Variant Call Format)](https://en.wikipedia.org/wiki/Variant_Call_Format) file format specifying identified variants in the BAM file, with respect to the passed reference file. When this file is provided, the computational running time of the program is improved. This variants file should be generated using the same BAM and reference files passed as parameters to program. The VCF output of [call ntvar](../ntvar) may be used as input to this program.

### Mutation Database

The mutation database describes specific mutations within the named genetic regions specified previously in the BED4 file. When provided to the tool, the surveillence and drug resistance category information are included in the output mutation annotation.

The entries in the "genetic regions" colummn of this database must match the "names" column of the provided BED4 file. Please refer to [Data Formats](../formats) for more information.

## Options

### Minimum Frequency

```text
-f, --min_freq FLOAT
```

The minimum observed frequency for a variant to reported. The default frequency is 0.01.

### Error Rate

```text
-e, --error_rate FLOAT
```

This is the expected substitution sequencing error rate. The default value is 0.0021 substitutions per sequenced base.

### Output

```text
-o, --output FILENAME
```

The file output location to write the identified amino acid mutations.

## Output

The amino acid mutations that exceed the minimum frequency threshold will be output in [AAVF](https://github.com/winhiv/aavf-spec/blob/master/AAVFv1.0.pdf) format. By default, the results will be printed to standard output. The user may direct the output to a file by specifying a file name with the ```-o / --output``` option.

## Examples

### Data

The following example data may be used to run the tool:

* [hiv.fasta](data/hiv.fasta)
* [variant.bam](data/variant.bam)
* [variant.bai](data/variant.bai)
* [variant.vcf](data/variant.vcf)
* [hiv_db.tsv](data/hiv_db.tsv)

### Example: No Database

#### Command

```
quasitools call aavar variant.bam hiv.fasta hiv.bed
```

#### Output

```text
##reference=hiv.fasta
##source=quasitools:aavar
##fileformat=AAVFv1.0
##fileDate=20190206
##INFO=<ID=SRVL,Number=.,Type=String,Description="Drug Resistance Surveillance">
##INFO=<ID=AC,Number=.,Type=String,Description="Alternate Codon">
##INFO=<ID=CAT,Number=.,Type=String,Description="Drug Resistance Category">
##INFO=<ID=ACF,Number=.,Type=Float,Description="Alternate Codon Frequency,for each Alternate Codon,in the same order aslisted.">
##INFO=<ID=RC,Number=1,Type=String,Description="Reference Codon">
##FILTER=<ID=af0.01,Description="Set if True; alt_freq<0.01">
#CHROM	GENE	POS	REF	ALT	FILTER	ALT_FREQ	COVERAGE	INFO
AF033819.3	gag	339	P	Q	PASS	1.0000	133	SRVL=.;AC=cAa;CAT=.;ACF=1.0000;RC=cca
AF033819.3	env	67	N	T	PASS	1.0000	116	SRVL=.;AC=aCt;CAT=.;ACF=1.0000;RC=aat
AF033819.3	gag	441	Y	S	PASS	1.0000	109	SRVL=.;AC=tCc;CAT=.;ACF=1.0000;RC=tac
AF033819.3	pol	141	I	S	PASS	1.0000	145	SRVL=.;AC=aGt;CAT=.;ACF=1.0000;RC=att
AF033819.3	vpu	34	L	I	PASS	1.0000	118	SRVL=.;AC=Ata;CAT=.;ACF=1.0000;RC=tta
AF033819.3	env	302	N	Y	PASS	1.0000	138	SRVL=.;AC=Tat;CAT=.;ACF=1.0000;RC=aat
AF033819.3	gag	3	A	P	PASS	1.0000	140	SRVL=.;AC=Ccg;CAT=.;ACF=1.0000;RC=gcg
AF033819.3	pol	246	Q	H	PASS	1.0000	139	SRVL=.;AC=caT;CAT=.;ACF=1.0000;RC=caa
AF033819.3	gag	230	E	D	PASS	1.0000	126	SRVL=.;AC=gaT;CAT=.;ACF=1.0000;RC=gaa
AF033819.3	pol	96	G	E	PASS	1.0000	128	SRVL=.;AC=gAa;CAT=.;ACF=1.0000;RC=gga
```

Observe that there is no information under ```INFO``` column of the output for ```SRVL``` (surveillance) and ```CAT``` (drug resistance category). This is because a mutation database was not provided.

### Example: With Database

#### Command

```bash
quasitools call aavar variant.bam hiv.fasta hiv.bed variant.vcf hiv_db.tsv
```

#### Output 

```text
##reference=hiv.fasta
##source=quasitools:aavar
##fileformat=AAVFv1.0
##fileDate=20190206
##INFO=<ID=SRVL,Number=.,Type=String,Description="Drug Resistance Surveillance">
##INFO=<ID=AC,Number=.,Type=String,Description="Alternate Codon">
##INFO=<ID=CAT,Number=.,Type=String,Description="Drug Resistance Category">
##INFO=<ID=ACF,Number=.,Type=Float,Description="Alternate Codon Frequency,for each Alternate Codon,in the same order aslisted.">
##INFO=<ID=RC,Number=1,Type=String,Description="Reference Codon">
##FILTER=<ID=af0.01,Description="Set if True; alt_freq<0.01">
#CHROM	GENE	POS	REF	ALT	FILTER	ALT_FREQ	COVERAGE	INFO
AF033819.3	gag	339	P	Q	PASS	1.0000	133	SRVL=.;AC=cAa;CAT=.;ACF=1.0000;RC=cca
AF033819.3	env	67	N	T	PASS	1.0000	116	SRVL=No;AC=aCt;CAT=minor;ACF=1.0000;RC=aat
AF033819.3	gag	441	Y	S	PASS	1.0000	109	SRVL=No;AC=tCc;CAT=minor;ACF=1.0000;RC=tac
AF033819.3	pol	141	I	S	PASS	1.0000	145	SRVL=.;AC=aGt;CAT=.;ACF=1.0000;RC=att
AF033819.3	vpu	34	L	I	PASS	1.0000	118	SRVL=Yes;AC=Ata;CAT=major;ACF=1.0000;RC=tta
AF033819.3	env	302	N	Y	PASS	1.0000	138	SRVL=.;AC=Tat;CAT=.;ACF=1.0000;RC=aat
AF033819.3	gag	3	A	P	PASS	1.0000	140	SRVL=Yes;AC=Ccg;CAT=major;ACF=1.0000;RC=gcg
AF033819.3	pol	246	Q	H	PASS	1.0000	139	SRVL=No;AC=caT;CAT=minor;ACF=1.0000;RC=caa
AF033819.3	gag	230	E	D	PASS	1.0000	126	SRVL=.;AC=gaT;CAT=.;ACF=1.0000;RC=gaa
AF033819.3	pol	96	G	E	PASS	1.0000	128	SRVL=Yes;AC=gAa;CAT=major;ACF=1.0000;RC=gga
```

Observe that under the ```INFO``` column of the output, that the ```SRVL``` and ```CAT``` information is included in the annotation, if the mutation was specified in the mutation database.



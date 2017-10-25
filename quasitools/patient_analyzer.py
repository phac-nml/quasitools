"""
Copyright Government of Canada 2017

Written by: Camy Tran, National Microbiology Laboratory,
            Public Health Agency of Canada

Licensed under the Apache License, Version 2.0 (the "License"); you may not use
this work except in compliance with the License. You may obtain a copy of the
License at:

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software distributed
under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR
CONDITIONS OF ANY KIND, either express or implied. See the License for the
specific language governing permissions and limitations under the License.
"""

import os
from quasitools.parsers.genes_file_parser import parse_genes_file
from quasitools.parsers.reference_parser import parse_references_from_fasta
from quasitools.parsers.mapped_read_parser import parse_mapped_reads_from_bam
from quasitools.nt_variant import NTVariantCollection
from quasitools.aa_variant import AAVariantCollection
from quasitools.mutations import MutationDB
from quasitools.aa_census import AACensus, CONFIDENT
import Bio.SeqIO


class PatientAnalyzer():
    def __init__(self, id, output_dir, reads, reference,
                 genes_file, mutation_db, quiet, consensus_pct):
        self.id = id
        self.output_dir = output_dir
        self.reads = reads
        self.reference = reference
        self.mutation_db = mutation_db
        self.genes_file = genes_file

        self.quiet = quiet
        self.consensus_pct = consensus_pct

        self.filtered = {}
        self.filtered["status"] = 0
        self.filtered["length"] = 0
        self.filtered["score"] = 0
        self.filtered["ns"] = 0

        self.downsample = {}
        self.downsample["status"] = 0
        self.downsample["size"] = 0

        self.input_size = 0
        self.determine_input_size()

        self.references = parse_references_from_fasta(self.reference)
        self.genes = parse_genes_file(genes_file, self.references[0].name)

        self.filtered_reads = "%s/filtered.fastq" % output_dir
        self.downsampled_reads = "%s/downsampled.fas" % output_dir

        os.mkdir(output_dir)

    def determine_input_size(self):
        sequences = Bio.SeqIO.parse(self.reads, "fastq")

        for seq in sequences:
            self.input_size += 1

    def filter_reads(self, filters):
        if not self.quiet:
            print("# Filtering reads...")

        filtered_reads_file = open(self.filtered_reads, "w+")

        seq_rec_obj = Bio.SeqIO.parse(self.reads, "fastq")

        for seq in seq_rec_obj:

            avg_score = (sum(seq.letter_annotations['phred_quality']) /
                         len(seq.letter_annotations['phred_quality']))

            length = len(seq.seq)

            if length <= filters["length_cutoff"]:
                self.filtered["length"] += 1
            elif avg_score <= filters["score_cutoff"]:
                self.filtered["score"] += 1
            elif filters['ns'] and 'n' in seq.seq.lower():
                self.filtered['ns'] += 1
            elif (length > filters["length_cutoff"] and
                  avg_score > filters["score_cutoff"]):
                Bio.SeqIO.write(seq, filtered_reads_file, "fastq")

        self.filtered["status"] = 1
        filtered_reads_file.close()

    def downsample_reads(self, target_coverage):
        if not self.quiet:
            print("\n# Downsampling reads...")

        seq_rec_obj = Bio.SeqIO.parse(self.filtered_reads, "fastq")

        total = 0
        for seq in seq_rec_obj:
            total += len(seq.seq)

        length = len(self.references[0].seq)

        raw_coverage = total / length
        total_reads = (self.input_size - self.filtered["length"] -
                       self.filtered["score"] - self.filtered["ns"])

        if raw_coverage > target_coverage:
            subsample_pct = float(target_coverage) / float(raw_coverage)
            subsample_size = subsample_pct * total_reads

            command = "seqtk sample -s 42 %s %i > %s" % (
                self.filtered_reads, subsample_size,
                self.downsampled_reads)

            os.system(command)

            if not self.quiet:
                print("\n# Reads downsampled to %0.2f%% to achieve "
                      "%i coverage..." % (subsample_pct*100, target_coverage))

            self.downsample["status"] = 1
            self.downsample["size"] = total_reads - subsample_size
        else:
            if not self.quiet:
                print("\n# Reads were not downsampled to %i coverage,"
                      "as raw coverage is %i..." % (
                          target_coverage, raw_coverage))

    def analyze_reads(self, filters, reporting_threshold, generate_consensus):
        # Map reads against reference using bowtietwo
        if not self.quiet:
            print("\n# Mapping reads...")

        bam = self.generate_bam()

        if not self.quiet:
            print("\n# Loading read mappings...")

        # cmd_consensus
        if generate_consensus:
            cons_seq_file = open("%s/consensus.fasta" % self.output_dir, "w+")

        mapped_read_collection_arr = []
        for r in self.references:
            mrc = parse_mapped_reads_from_bam(r, bam)
            mapped_read_collection_arr.append(mrc)
            if generate_consensus:
                cons_seq_file.write('>{0}_{1}_{2}\n{3}'.format(
                    'blah', reporting_threshold, r.name,
                    mrc.to_consensus(self.consensus_pct)))

        if generate_consensus:
            cons_seq_file.close()

        # cmd_callntvar
        if not self.quiet:
            print("\n# Identifying variants...")

        variants = NTVariantCollection.from_mapped_read_collections(
            filters["error_rate"], self.references,
            *mapped_read_collection_arr)

        variants.filter('q%s' % filters["min_qual"],
                        'QUAL<%s' % filters["min_qual"], True)
        variants.filter('ac%s' % filters["min_ac"],
                        'AC<%s' % filters["min_ac"], True)
        variants.filter('dp%s' % filters["min_dp"],
                        'DP<%s' % filters["min_dp"], True)

        vcf_file = open("%s/hydra.vcf" % self.output_dir, "w+")
        vcf_file.write(variants.to_vcf_file())
        vcf_file.close()

        # cmd_aa_census
        if not self.quiet:
            print("\n# Masking filtered variants...")

        for mrc in mapped_read_collection_arr:
            mrc.mask_unconfident_differences(variants)

        if not self.quiet:
            print("\n# Building amino acid census...")

        # Determine which frames our genes are in
        frames = set()

        for gene in self.genes:
            frames.add(self.genes[gene]['frame'])

        aa_census = AACensus(self.reference, mapped_read_collection_arr,
                             self.genes, frames)

        coverage_file = open("%s/coverage_file.csv" % self.output_dir, "w+")
        coverage_file.write(aa_census.coverage(next(iter(frames))))
        coverage_file.close()

        # cmd_aavariants
        if not self.quiet:
            print("\n# Finding amino acid mutations...")

        # Create AAVar collection and print the hmcf file
        aa_vars = AAVariantCollection.from_aacensus(
            aa_census, next(iter(frames)))

        # Filter for mutant frequency
        aa_vars.filter('mf%s' % filters['min_freq'],
                       'freq<%s' % filters['min_freq'], True)

        # Build the mutation database and update collection
        if self.mutation_db is not None:
            mutation_db = MutationDB(self.mutation_db, self.genes)
            aa_vars.apply_mutation_db(mutation_db)

        mut_report = open("%s/mutation_report.hmcf" % self.output_dir, "w+")
        mut_report.write(aa_vars.to_hmcf_file(CONFIDENT))
        mut_report.close()

        # cmd_drmutations
        if not self.quiet:
            print("\n# Writing drug resistant mutation report...")

        dr_report = open("%s/dr_report.csv" % self.output_dir, "w+")
        dr_report.write(aa_vars.report_dr_mutations(mutation_db,
                                                    reporting_threshold))
        dr_report.close()

        self.output_stats(mapped_read_collection_arr)

    def generate_bam(self):
        """ Runs bowtietwo local alignment on self.reads
            to generate a bam file """

        reads = self.filtered_reads

        if self.downsample["status"] == 1:
            reads = self.downsampled_reads

        sorted_bam_fn = "%s/align.bam" % self.output_dir
        bowtietwo_bam_output = sorted_bam_fn[0:sorted_bam_fn.index(".")]
        bam_fn = "%s/tmp.bam" % self.output_dir
        sam_fn = "%s/tmp.sam" % self.output_dir

        # create the files
        bam_fh = open(bam_fn, "w+")
        sam_fh = open(sam_fn, "w+")

        bowtietwo_index = self.reference[0:self.reference.index(".")]

        bowtietwo_cmd = (("bowtie2 --local --rdg '8,3' --rfg '8,3' "
                          "--ma 1 --mp '2,2' -S %s -x %s -U %s") %
                         (sam_fn, bowtietwo_index, reads))

        os.system(bowtietwo_cmd)

        # Convert sam output to bam output
        sam_to_bam_cmd = "samtools view -bt %s.fai -o %s %s" % (self.reference,
                                                                bam_fn, sam_fn)

        os.system(sam_to_bam_cmd)

        # Sort bam output
        sort_bam_cmd = "samtools sort %s %s" % (bam_fn, bowtietwo_bam_output)

        os.system(sort_bam_cmd)

        # Index bam output
        index_bam_cmd = "samtools index %s" % sorted_bam_fn

        os.system(index_bam_cmd)

        bam_fh.close()
        sam_fh.close()

        os.unlink(bam_fn)
        os.unlink(sam_fn)

        return sorted_bam_fn

    def output_stats(self, mapped_read_collection_arr):
        mr_len = len(mapped_read_collection_arr[0].mapped_reads)

        stats_report = open("%s/stats.txt" % self.output_dir, "w+")

        stats_report.write("Input Size: %i\n" % self.input_size)
        stats_report.write("Number of reads filtered due to length: %i\n" %
                           self.filtered["length"])
        stats_report.write(("Number of reads filtered due to average "
                            "quality score: %i\n") % self.filtered["score"])
        stats_report.write(("Number of reads filtered due to presence "
                            "of Ns: %i\n") % self.filtered["ns"])
        if self.downsample["status"]:
            stats_report.write(("Number of reads filtered due to excess "
                                "coverage: %i\n") % self.downsample["size"])
            stats_report.write(("Number of reads filtered due to poor "
                                "mapping: %i\n") %
                               (self.input_size - self.filtered["length"] -
                                self.filtered["score"] - self.filtered["ns"] -
                                self.downsample["size"] - mr_len))
        else:
            stats_report.write("Number of reads filtered due to excess "
                               "coverage: 0\n")
            stats_report.write(("Number of reads filtered due to poor "
                                "mapping: %i\n") %
                               (self.input_size - self.filtered["length"] -
                                self.filtered["score"] - self.filtered["ns"] -
                                mr_len))
        stats_report.write("Percentage of reads filtered: %0.2f" %
                           (float(self.input_size - mr_len) /
                            self.input_size * 100))

        stats_report.close()

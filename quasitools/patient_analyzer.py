"""
Copyright Government of Canada 2017 - 2018

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
from quasitools.quality_control import QualityControl
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

        self.input_size = 0
        self.determine_input_size()

        self.references = parse_references_from_fasta(self.reference)
        self.genes = parse_genes_file(genes_file, self.references[0].name)

        self.quality = QualityControl()

        if not os.path.isdir(output_dir):
            os.mkdir(output_dir)

        self.filtered_reads_dir = "%s/filtered.fastq" % output_dir

    def determine_input_size(self):
        sequences = Bio.SeqIO.parse(self.reads, "fastq")

        for seq in sequences:
            self.input_size += 1

    def analyze_reads(self, fasta_id, quality_filters, variant_filters,
                      reporting_threshold, generate_consensus):
        # Calls quality_control function
        if not self.quiet:
            print("# Performing quality control on reads...")

        self.quality.filter_reads(self.reads,
                                  self.filtered_reads_dir,
                                  quality_filters)

        # Map reads against reference using bowtietwo
        if not self.quiet:
            print("# Mapping reads...")

        bam = self.generate_bam(fasta_id)

        if not self.quiet:
            print("# Loading read mappings...")

        # cmd_consensus
        if generate_consensus:
            cons_seq_file = open("%s/consensus.fasta" % self.output_dir, "w+")

        mapped_read_collection_arr = []
        for r in self.references:
            mrc = parse_mapped_reads_from_bam(r, bam)
            mapped_read_collection_arr.append(mrc)
            if generate_consensus:
                cons_seq_file.write('>{0}_{1}_{2}\n{3}'.format(
                    fasta_id, reporting_threshold, r.name,
                    mrc.to_consensus(self.consensus_pct)))

        if generate_consensus:
            cons_seq_file.close()

        # cmd_callntvar
        if not self.quiet:
            print("# Identifying variants...")

        variants = NTVariantCollection.from_mapped_read_collections(
            variant_filters["error_rate"], self.references,
            *mapped_read_collection_arr)

        variants.filter('q%s' % variant_filters["min_qual"],
                        'QUAL<%s' % variant_filters["min_qual"], True)
        variants.filter('ac%s' % variant_filters["min_ac"],
                        'AC<%s' % variant_filters["min_ac"], True)
        variants.filter('dp%s' % variant_filters["min_dp"],
                        'DP<%s' % variant_filters["min_dp"], True)

        vcf_file = open("%s/hydra.vcf" % self.output_dir, "w+")
        vcf_file.write(variants.to_vcf_file())
        vcf_file.close()

        # cmd_aa_census
        if not self.quiet:
            print("# Masking filtered variants...")

        for mrc in mapped_read_collection_arr:
            mrc.mask_unconfident_differences(variants)

        if not self.quiet:
            print("# Building amino acid census...")

        # Determine which frames our genes are in
        frames = set()

        for gene in self.genes:
            frames.add(self.genes[gene]['frame'])

        aa_census = AACensus(self.reference, mapped_read_collection_arr,
                             self.genes, frames)

        coverage_file = open("%s/coverage_file.csv" % self.output_dir, "w+")
        coverage_file.write(aa_census.coverage(frames))
        coverage_file.close()

        # cmd_aavariants
        if not self.quiet:
            print("# Finding amino acid mutations...")

        # Create AAVar collection and print the hmcf file
        aa_vars = AAVariantCollection.from_aacensus(aa_census)

        # Filter for mutant frequency
        aa_vars.filter('mf%s' % variant_filters['min_freq'],
                       'freq<%s' % variant_filters['min_freq'], True)

        # Build the mutation database and update collection
        if self.mutation_db is not None:
            mutation_db = MutationDB(self.mutation_db, self.genes)
            aa_vars.apply_mutation_db(mutation_db)

        mut_report = open("%s/mutation_report.hmcf" % self.output_dir, "w+")
        mut_report.write(aa_vars.to_hmcf_file(CONFIDENT))
        mut_report.close()

        # cmd_drmutations
        if not self.quiet:
            print("# Writing drug resistant mutation report...")

        dr_report = open("%s/dr_report.csv" % self.output_dir, "w+")
        dr_report.write(aa_vars.report_dr_mutations(mutation_db,
                                                    reporting_threshold))
        dr_report.close()

        self.output_stats(mapped_read_collection_arr)

    # This is a helper method that generates the bam file.
    # It takes as an argument the fasta_id, which is used by bowtie2 as the
    # RG-ID in the output bam file.
    def generate_bam(self, fasta_id):
        """ Runs bowtietwo local alignment on self.reads
            to generate a bam file """

        sorted_bam_fn = "%s/align.bam" % self.output_dir
        bowtietwo_bam_output = sorted_bam_fn[0:sorted_bam_fn.rindex(".")]
        bam_fn = "%s/tmp.bam" % self.output_dir
        sam_fn = "%s/tmp.sam" % self.output_dir

        # create the files
        bam_fh = open(bam_fn, "w+")
        sam_fh = open(sam_fn, "w+")

        bowtietwo_index = self.reference[0:self.reference.rindex(".")]

        bowtietwo_cmd = (("bowtie2 --local --rdg '8,3' --rfg '8,3' "
                          "--rg-id %s --ma 1 --mp '2,2' -S %s -x %s "
                          "-U %s") % (fasta_id, sam_fn, bowtietwo_index,
                                      self.filtered_reads_dir))

        os.system(bowtietwo_cmd)

        # Convert sam output to bam output
        sam_to_bam_cmd = "samtools view -bt %s.fai -o %s %s" % (self.reference,
                                                                bam_fn, sam_fn)

        os.system(sam_to_bam_cmd)

        # Sort bam output
        sort_bam_cmd = "samtools sort %s -T %s -o %s" % (bam_fn,
                                                         bowtietwo_bam_output,
                                                         sorted_bam_fn)

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
        self.amount_filtered = self.quality.get_amount_filtered()
        mr_len = len(mapped_read_collection_arr[0].mapped_reads)

        stats_report = open("%s/stats.txt" % self.output_dir, "w+")

        stats_report.write("Input Size: %i\n" % self.input_size)
        stats_report.write("Number of reads filtered due to length: %i\n" %
                           self.amount_filtered["length"])
        stats_report.write(("Number of reads filtered due to average "
                            "quality score: %i\n")
                           % self.amount_filtered["score"])
        stats_report.write(("Number of reads filtered due to presence "
                            "of Ns: %i\n") % self.amount_filtered["ns"])
        stats_report.write("Number of reads filtered due to excess "
                           "coverage: 0\n")
        stats_report.write(("Number of reads filtered due to poor "
                            "mapping: %i\n") %
                           (self.input_size - self.amount_filtered["length"] -
                            self.amount_filtered["score"] -
                            self.amount_filtered["ns"] -
                            mr_len))
        stats_report.write("Percentage of reads filtered: %0.2f" %
                           (float(self.input_size - mr_len) /
                            self.input_size * 100))

        stats_report.close()

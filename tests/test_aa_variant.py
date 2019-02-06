"""
Copyright Government of Canada 2017

Written by: Cole Peters, Eric Chubaty, National Microbiology Laboratory, Public Health Agency of Canada

Licensed under the Apache License, Version 2.0 (the "License"); you may not use
this work except in compliance with the License. You may obtain a copy of the
License at:

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software distributed
under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR
CONDITIONS OF ANY KIND, either express or implied. See the License for the
specific language governing permissions and limitations under the License.
"""

import pytest
import os
import re
import string
import PyAAVF.parser as parser
from quasitools.aa_census import AACensus, CONFIDENT
from quasitools.aa_variant import AAVariantCollection
from quasitools.nt_variant import NTVariantCollection
from quasitools.mutations import MutationDB
from quasitools.parsers.mapped_read_parser import parse_mapped_reads_from_bam
from quasitools.parsers.genes_file_parser import parse_BED4_file
from quasitools.parsers.reference_parser import parse_references_from_fasta
from quasitools.parsers.nt_variant_file_parser \
    import parse_nt_variants_from_vcf

TEST_PATH = os.path.dirname(os.path.abspath(__file__))
VARIANTS_FILE = TEST_PATH + "/data/output/nt_variants.vcf"
VALID_AA_VARIANTS_AAVF = TEST_PATH + "/data/output/mutation_report.aavf"
VALID_DR_MUTATIONS_CSV = TEST_PATH + "/data/output/dr_report.csv"


class TestAAVariant:

    @classmethod
    def setup(self):
        self.reference = TEST_PATH + "/data/hxb2_pol.fas"
        bam = TEST_PATH + "/data/align.bam"
        BED4_file = TEST_PATH + "/data/hxb2_pol.bed"
        mutation_db = TEST_PATH + "/data/mutation_db.tsv"
        min_freq = 0.01

        rs = parse_references_from_fasta(self.reference)

        mapped_read_collection_arr = []
        for r in rs:
            # Create a MappedReadCollection object
            mapped_read_collection_arr.append(
                parse_mapped_reads_from_bam(r, bam))

        variants_obj = parse_nt_variants_from_vcf(VARIANTS_FILE, rs)

        # Mask the unconfident differences
        for mrc in mapped_read_collection_arr:
            mrc.mask_unconfident_differences(variants_obj)

        # Parse the genes from the gene file
        genes = parse_BED4_file(BED4_file, rs[0].name)

        # Determine which frames our genes are in
        frames = set()

        for gene in genes:
            frames.add(genes[gene]['frame'])

        # Create an AACensus object
        aa_census = AACensus(self.reference, mapped_read_collection_arr, genes,
                             frames)

        # Find the AA mutations
        self.aa_collection = AAVariantCollection.from_aacensus(
            aa_census)

        # Build the mutation database
        self.mutation_db = MutationDB(mutation_db, genes)

    def test_dr_mutations(self):
        reporting_threshold = 1

        with open(VALID_DR_MUTATIONS_CSV, "r") as input:
            # Read both sample and collection, remove nl, and compare
            dr_mutations = [x.replace("\n", "") for x in input.readlines()]

            # Apply filter to the collection
            self.aa_collection.filter('mf0.01', 'freq<0.01', True)

            # Get the report and compare outputs
            collection = self.aa_collection.report_dr_mutations(
                self.mutation_db, reporting_threshold
            ).split("\n")

            assert dr_mutations == collection

    def test_aa_variants(self):
        aavf_path = TEST_PATH + "/data/output/temp.aavf"
        aavf_out = open(aavf_path, "w")

        # Read from file and make sure there are no empty lines
        with open(VALID_AA_VARIANTS_AAVF, "r") as input:
            valid_variants = input.read()

        # Sort and filter so comparison order will be fine afterwards
        valid_variants_lines = sorted(filter(None, valid_variants.split("\n")))

        # Apply the filter to the collection
        self.aa_collection.filter('af0.01', 'freq<0.01', True)

        # Do the thing with the mutation_db
        self.aa_collection.apply_mutation_db(self.mutation_db)

        aavf_obj = self.aa_collection.to_aavf_obj(
            "test", [os.path.basename(self.reference)], CONFIDENT)
        records = list(aavf_obj)

        writer = parser.Writer(aavf_out, aavf_obj)

        for record in records:
            writer.write_record(record)

        aavf_out.close()

        with open(aavf_path, "r") as input:
            aa_variants = input.read()

        # Make sure it's sorted and has no empty strings
        aa_variants_lines = sorted(filter(None, aa_variants.split("\n")))

        # Check the length
        assert len(valid_variants_lines) == len(aa_variants_lines)

        # Make sure all the tokens that need to be there are there
        for pos in range(0, len(valid_variants_lines)):
            if valid_variants_lines[pos][0:1] != "#":
                valid_variants_tokens = \
                    re.split("[,=;\t]", valid_variants_lines[pos].rstrip())

                aa_variants_tokens = re.split(
                    "[,=;\t]", aa_variants_lines[pos])

                for token in aa_variants_tokens:
                    assert token in valid_variants_tokens

        writer.close()
        # remove temp test file

    def test_aa_variants_nodb(self):
        # Same as previous test but for no mutation db
        aavf_path = TEST_PATH + "/data/output/temp.aavf"
        aavf_out = open(aavf_path, "w")

        # Read from file and make sure there are no empty lines
        with open(VALID_AA_VARIANTS_AAVF, "r") as input:
            valid_variants = input.read()

        valid_variants_lines = sorted(filter(None, valid_variants.split("\n")))

        # Replace category and surveillance with "."s
        # Okay because comparisons only done on non "#" lines
        for i, x in enumerate(valid_variants_lines):
            tokens = x.split(";")

            # Change result to be what it would be without a db
            if len(tokens) > 2:
                x = x.replace(tokens[-2], "CAT=.", 1)
                x = x.replace(tokens[-1], "SRVL=.", 1)
                valid_variants_lines[i] = x

        # Apply the filter to the collection
        self.aa_collection.filter('af0.01', 'freq<0.01', True)

        aavf_obj = self.aa_collection.to_aavf_obj(
            "test", [os.path.basename(self.reference)], CONFIDENT)
        records = list(aavf_obj)

        writer = parser.Writer(aavf_out, aavf_obj)

        for record in records:
            writer.write_record(record)

        aavf_out.close()

        with open(aavf_path, "r") as input:
            aa_variants = input.read()

        # Make sure it's sorted and has no empty strings
        aa_variants_lines = sorted(filter(None, aa_variants.split("\n")))

        # Check the length
        assert len(valid_variants_lines) == len(aa_variants_lines)

        # Make sure all the tokens that need to be there are there
        for pos in range(0, len(valid_variants_lines)):
            if valid_variants_lines[pos][0:1] != "#":
                valid_variants_tokens = \
                    re.split("[,=;\t]", valid_variants_lines[pos].rstrip())

                aa_variants_tokens = re.split(
                    "[,=;\t]", aa_variants_lines[pos])

                for token in aa_variants_tokens:
                    assert token in valid_variants_tokens

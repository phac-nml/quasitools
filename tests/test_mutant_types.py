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
import pytest
from quasitools.codon_variant import CodonVariantCollection
from quasitools.nt_variant import NTVariantCollection
from quasitools.aa_census import AACensus
from quasitools.parsers.mapped_read_parser import \
    parse_mapped_reads_from_bam
from quasitools.parsers.reference_parser import \
    parse_references_from_fasta
from quasitools.parsers.genes_file_parser import parse_genes_file


TEST_PATH = os.path.dirname(os.path.abspath(__file__))
VALID_MUTANT_TYPES_CSV = TEST_PATH + "/data/output/mutant_types.csv"

class TestMutantTypes:

    @classmethod
    def setup(self):
        bam = TEST_PATH + "/data/align.bam"
        reference = TEST_PATH + "/data/hxb2_pol.fas"
        genes_file = TEST_PATH + "/data/hxb2_pol.bed"
        error_rate = 0.0038

        rs = parse_references_from_fasta(reference)
        mapped_read_collection_arr = []

        # Create a MappedReadCollection object
        for r in rs:
            mapped_read_collection_arr.append(
                parse_mapped_reads_from_bam(r, bam))

        variants = NTVariantCollection.from_mapped_read_collections(
                        error_rate, rs, *mapped_read_collection_arr)
        variants.filter('q30', 'QUAL<30', True)
        variants.filter('ac5', 'AC<5', True)
        variants.filter('dp100', 'DP<100', True)

        # Mask the unconfident differences
        for mrc in mapped_read_collection_arr:
            mrc.mask_unconfident_differences(variants)

        # Parse the genes from the gene file
        genes = parse_genes_file(genes_file, rs[0].name)

        # Determine which frames our genes are in
        frames = set()

        for gene in genes:
            frames.add(genes[gene]['frame'])

        aa_census = AACensus(reference,
                             mapped_read_collection_arr,
                             genes,
                             frames)

        self.codon_variants = CodonVariantCollection.from_aacensus(
                            aa_census, next(iter(frames)))

    def test_mutant_types(self): 
        offset = 1269

        # Read from file and make sure there are no empty lines
        with open(VALID_MUTANT_TYPES_CSV, "r") as input:
            valid_mutant_types = input.read()

        # Sort and filter for comparison
        valid_mutant_types_lines = sorted(filter(None, 
            valid_mutant_types.split("\n")))

        # Create the report string
        mutant_types = self.codon_variants.to_csv_file(offset)

        # Split into lines and sort
        mutant_types_lines = sorted(mutant_types.split("\n"))

        assert len(valid_mutant_types_lines) == len(mutant_types_lines)

        # Compare each line in the test report to the valid report
        for pos in range(0, len(valid_mutant_types_lines)):
            if valid_mutant_types_lines[pos][0:1] != "#":
                assert valid_mutant_types_lines[pos] == \
                    mutant_types_lines[pos]

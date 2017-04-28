"""
Copyright Government of Canada 2017

Written by: Cole Peters, National Microbiology Laboratory, Public Health Agency of Canada

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
from quasitools.aa_census import AACensus
from quasitools.parsers.mapped_read_parser import parse_mapped_reads_from_bam
from quasitools.parsers.genes_file_parser import parse_genes_file
from quasitools.parsers.reference_parser import parse_references_from_fasta


TEST_PATH = os.path.dirname(os.path.abspath(__file__))
VALID_COVERAGE_CSV = TEST_PATH + "/data/output/coverage_file.csv"

class TestMappedRead:
    @classmethod
    def setup_class(self):
        reference = TEST_PATH + "/data/hxb2_pol.fas"
        bam = TEST_PATH + "/data/align.bam"
        genes_file = TEST_PATH + "/data/hxb2_pol.bed"

        rs = parse_references_from_fasta(reference)

        mapped_read_collection_arr = []
        for r in rs:
            # create MappedReadCollection object
            mapped_read_collection_arr.append(
                parse_mapped_reads_from_bam(r, bam))

        genes = parse_genes_file(genes_file, rs[0].name)

        # Determine which frames our genes are in
        self.frames = set()

        for gene in genes:
            self.frames.add(genes[gene]["frame"])

        self.aa_census = AACensus(reference, mapped_read_collection_arr,
                                  genes, self.frames)

    def test_coverage(self):

        with open(VALID_COVERAGE_CSV, "r") as input:
            coverage = input.read()

        assert self.aa_census.coverage(next(iter(self.frames))) == coverage

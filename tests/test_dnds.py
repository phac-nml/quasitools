"""
Copyright Government of Canada 2015

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

import pdb
import pytest
import os
from quasitools.parsers.reference_parser import parse_references_from_fasta
from quasitools.parsers.codon_variant_file_parser import parse_codon_variants


TEST_PATH = os.path.dirname(os.path.abspath(__file__))
VALID_DNDS_REPORT = TEST_PATH + "/data/output/dnds_report.csv"

class TestDnds:
    @classmethod
    def setup(self):
        csv = TEST_PATH + "/data/output/mutant_types.csv"
        reference = TEST_PATH + "/data/hxb2_pol.fas"
        self.offset = 1269

        rs = parse_references_from_fasta(reference)
        self.ref_seq = rs[0].seq

        self.codon_variants = parse_codon_variants(csv, rs)


    def test_dnds(self):
        # Read from file and make sure there are no empty lines
        with open(VALID_DNDS_REPORT, "r") as input:
            valid_report = input.read()

        # Sort and filter for comparison
        valid_dnds_values = sorted(filter(None, 
            valid_report.split("\n")))

        # Create the report string
        test_report = self.codon_variants.report_dnds_values(self.ref_seq, self.offset)

        # Split into lines and sort
        test_values = sorted(test_report.split("\n"))

        assert len(valid_dnds_values) == len(test_values)
        
        # Compare each line in the test report to the valid report
        for pos in range(0, len(valid_dnds_values)):
            if valid_dnds_values[pos][0:1] != "#":
                assert valid_dnds_values[pos] == \
                    test_values[pos]


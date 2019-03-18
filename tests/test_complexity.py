"""
# =============================================================================
Copyright Government of Canada 2019

Written by: Ahmed Kidwai, Public Health Agency of Canada,
    National Microbiology Laboratory

Funded by the National Micriobiology Laboratory and the Genome Canada / Alberta
    Innovates Bio Solutions project "Listeria Detection and Surveillance
    using Next Generation Genomics"

Licensed under the Apache License, Version 2.0 (the "License"); you may not use
this file except in compliance with the License. You may obtain a copy of the
License at:

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software distributed
under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR
CONDITIONS OF ANY KIND, either express or implied. See the License for the
specific language governing permissions and limitations under the License.

# =============================================================================
"""
import os
import pytest
import quasitools.commands.cmd_complexity as complexity
import quasitools.haplotype as haplotype
from click.testing import CliRunner
import click

front = ["AAGAGTTGTGAGGAGTACTCACCTCAGTAGACAAGGAGAGCTATATCGCGCAGCATCGCTTCGATA"+
          "TTTTGCTCCTACGCATCCACACGTTGAAAGGGCC", 
         "AAGAGTTGTGAGGAGTACTCACCTCAGTAGACAAGGAGAGCTATATCGCGCGGCATCGCTTCGATA"+
          "TTTTGCTCCTACGCATCCACACGTTGAAAGGGCA",
         "TAGAGTTGTGAGGAGTACTCACCTCAGTAGACAAGGAGAGCTATATCGCGCTGCATCGCTTCGATA"+
         "TTTTGCTCCTACGCATCCACACGTTGAAAGGGCA",
         "AAGAGTTGTGAGGAGTACTCACCTCAGTAGACAAGGAGAGCTATATCGCGCAGCATCGCTTCGATA"+
         "TTTTGCTCCTACGCATCCACACGTTGAAAGGGCC",
         "TAGAGTTGTGAGGAGTACTCACCTCAGTAGACAAGGAGAGCTATATCGCGCAGCATCGCTTCGATA"+
         "TTTTGCTCCTACGCATCCACACGTTGAAAGGGCG",
         "AAGAGTTGTGAGGAGTACTCACCTCAGTAGACAAGGAGAGCTATATCGCGCAGCATCGCTTCGATA"+
         "TTTTGCTCCTACGCATCCACACGTTGAAAGGGCC",
         "AAGAGTTGTGAGGAGTACTCACCTCAGTAGACAAGGAGAGCTATATCGCGCAGCATCGCTTCGATA"+
         "TTTTGCTCCTACGCATCCACACGTTGAAAGGGCC",
         "AAGAGTTGTGAGGAGTACTCACCTCAGTAGACAAGGAGAGCTATATCGCGCAGCATCGCTTCGATA"+
         "TTTTGCTCCTACGCATCCACACGTTGAAAGGGCC",
         "AAGAGTTGTGAGGAGTACTCACCTCAGTAGACAAGGAGAGCTATATCGCGCAGCATCGCTTCGATA"+
         "TTTTGCTCCTACGCATCCACACGTTGAAAGGGCC"
         "CAGAGTTGTGAGGAGTACTCACCTCAGTAGACAAGGAGAGCTATATCGCGCTGCATCGCTTCGATA"+
         "TTTTGCTCCTACGCATCCACACGTTGAAAGGGCA"]

back = ["AAATGATGGTCGCTCCACTTTCTCTCTAGAACGATCGTTATGTCGAAGTCCTAGCACTACGGGGTAC"+
        "TACAAGTATACCACGAGCCTACTTACCTCCAAA",
        "AAATGATGGTCGCTCCACTTTCTCTCTAGAACGATCGTTATGTCGAAGTCCAAGCACTACGGGGTAC \
         TACAAGTATACCACGAGCCTACTTACCTCCAAA",
        "AAATGATGGTCGCTCCACTTTCTCTCTAGAACGATCGTTATGTCGAAGTCCAAGCACTACGGGGTAC \
         TACAAGTATACCACGAGCCTACTTACCTCCAAA",
        "GAATGATGGTCGCTCCACTTTCTCTCTAGAACGATCGTTATGTCGAAGTCCAAGCACTACGGGGTAC \
         TACAAGTATACCACGAGCCTACTTACCTCCAAA"
        "AAATGATGGTCGCTCCACTTTCTCTCTAGAACGATCGTTATGTCGAAGTCCGAGCACTACGGGGTAC \
         TACAAGTATACCACGAGCCTACTTACCTCCAAT",
        "CAATGATGGTCGCTCCACTTTCTCTCTAGAACGATCGTTATGTCGAAGTCCGAGCACTACGGGGTAC \
         TACAAGTATACCACGAGCCTACTTACCTCCAAT",
        "AAATGATGGTCGCTCCACTTTCTCTCTAGAACGATCGTTATGTCGAAGTCCAAGCACTACGGGGTAC \
         TACAAGTATACCACGAGCCTACTTACCTCCAAA",
        "CAATGATGGTCGCTCCACTTTCTCTCTAGAACGATCGTTATGTCGAAGTCCAAGCACTACGGGGTAC \
         TACAAGTATACCACGAGCCTACTTACCTCCAAT",
        "CAATGATGGTCGCTCCACTTTCTCTCTAGAACGATCGTTATGTCGAAGTCCAAGCACTACGGGGTAC \
         TACAAGTATACCACGAGCCTACTTACCTCCAAC",
        "GAATGATGGTCGCTCCACTTTCTCTCTAGAACGATCGTTATGTCGAAGTCCTAGCACTACGGGGTAC \
         TACAAGTATACCACGAGCCTACTTACCTCCAAT"]

TEST_PATH = os.path.dirname(os.path.abspath(__file__))



class Test_CONSENSUS:

    def test_consensus_start(build_consensus):

        haplotype_list = []

        for i in range(len(front)):

            if front[i] not in haplotype_list:
                haplotype_list.append(haplotype.Haplotype(front[i]))


        expected_consensus = "AAGAGTTGTGAGGAGTACTCACCTCAGTAGACAAGGAGAGCTATATCGCGCAGCATCGCTTCGATATTTTGCTCCTACGCATCCACACGTTGAAAGGGCC"

        result = haplotype.build_consensus_from_haplotypes( \
                 haplotype_list)
        assert expected_consensus == result
        


# Test to see if method runs.
# NOTE: will not produce an output file.
class Test_BAM_Complexity:
    @classmethod
    def setup(self):
        self.bam_location = TEST_PATH + '/data/complexity.bam'
        self.reference_location =  TEST_PATH + '/data/complexity_reference.fasta'
        self.output_location =  TEST_PATH + '/data/output.csv'

    def test_complexity_bam(self):
        
        runner = CliRunner()
        result = runner.invoke(complexity.bam, [self.reference_location,\
                self.bam_location, "1",'--output_location', \
                self.output_location])
        
        # If method ran successfully the exit code is 0.
        assert result.exit_code == 0
        # output checks for print messages at the end of method.
        # the bam method in complexity has no print message.
        assert result.output == ""

class Test_FASTA_Complexity:
    @classmethod
    def setup(self):
        self.fasta_location = TEST_PATH + '/data/complexity_fasta.fasta'
        self.output_location =  TEST_PATH + '/data/output.csv'


    def test_complexity_fasta(self):
        
        runner = CliRunner()
        result = runner.invoke(complexity.fasta, [self.fasta_location, "--output_location", self.output_location])
        
        # If method ran successfully the exit code is 0.
        assert result.exit_code == 0
        # output checks for print messages at the end of method.
        # the bam method in complexity has no print message.
        assert result.output == ""

            

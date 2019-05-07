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
import quasitools.calculate as calculate
import quasitools.pileup as pileup
from click.testing import CliRunner
import click
import csv

front = ["AAGAGTTGTGAGGAGTACTCACCTCAGTAGACAAGGAGAGCTATATCGCGCAGCATCGCTTCGATATTTTGCTCCTACGCATCCACACGTTGAAAGGGCC",
         "AAGAGTTGTGAGGAGTACTCACCTCAGTAGACAAGGAGAGCTATATCGCGCGGCATCGCTTCGATATTTTGCTCCTACGCATCCACACGTTGAAAGGGCA",
         "TAGAGTTGTGAGGAGTACTCACCTCAGTAGACAAGGAGAGCTATATCGCGCTGCATCGCTTCGATATTTTGCTCCTACGCATCCACACGTTGAAAGGGCA",
         "AAGAGTTGTGAGGAGTACTCACCTCAGTAGACAAGGAGAGCTATATCGCGCAGCATCGCTTCGATATTTTGCTCCTACGCATCCACACGTTGAAAGGGCC",
         "TAGAGTTGTGAGGAGTACTCACCTCAGTAGACAAGGAGAGCTATATCGCGCAGCATCGCTTCGATATTTTGCTCCTACGCATCCACACGTTGAAAGGGCG",
         "AAGAGTTGTGAGGAGTACTCACCTCAGTAGACAAGGAGAGCTATATCGCGCAGCATCGCTTCGATATTTTGCTCCTACGCATCCACACGTTGAAAGGGCC",
         "AAGAGTTGTGAGGAGTACTCACCTCAGTAGACAAGGAGAGCTATATCGCGCAGCATCGCTTCGATATTTTGCTCCTACGCATCCACACGTTGAAAGGGCC",
         "AAGAGTTGTGAGGAGTACTCACCTCAGTAGACAAGGAGAGCTATATCGCGCAGCATCGCTTCGATATTTTGCTCCTACGCATCCACACGTTGAAAGGGCC",
         "AAGAGTTGTGAGGAGTACTCACCTCAGTAGACAAGGAGAGCTATATCGCGCAGCATCGCTTCGATATTTTGCTCCTACGCATCCACACGTTGAAAGGGCC",
         "CAGAGTTGTGAGGAGTACTCACCTCAGTAGACAAGGAGAGCTATATCGCGCTGCATCGCTTCGATATTTTGCTCCTACGCATCCACACGTTGAAAGGGCA"]

haplotypes_list_front = [
        haplotype.Haplotype("AAGAGTTGTGAGGAGTACTCACCTCAGTAGACAAGGAGAGCTATATCGCGCAGCATCGCTTCGATATTTTGCTCCTACGCATCCACACGTTGAAAGGGCC", 6),
        haplotype.Haplotype("AAGAGTTGTGAGGAGTACTCACCTCAGTAGACAAGGAGAGCTATATCGCGCGGCATCGCTTCGATATTTTGCTCCTACGCATCCACACGTTGAAAGGGCA", 1),
        haplotype.Haplotype("TAGAGTTGTGAGGAGTACTCACCTCAGTAGACAAGGAGAGCTATATCGCGCTGCATCGCTTCGATATTTTGCTCCTACGCATCCACACGTTGAAAGGGCA", 1),
        haplotype.Haplotype("TAGAGTTGTGAGGAGTACTCACCTCAGTAGACAAGGAGAGCTATATCGCGCAGCATCGCTTCGATATTTTGCTCCTACGCATCCACACGTTGAAAGGGCG", 1),
        haplotype.Haplotype("CAGAGTTGTGAGGAGTACTCACCTCAGTAGACAAGGAGAGCTATATCGCGCTGCATCGCTTCGATATTTTGCTCCTACGCATCCACACGTTGAAAGGGCA", 1)
        ]

back = ["AAATGATGGTCGCTCCACTTTCTCTCTAGAACGATCGTTATGTCGAAGTCCTAGCACTACGGGGTAC"+
         "TACAAGTATACCACGAGCCTACTTACCTCCAAA",
        "AAATGATGGTCGCTCCACTTTCTCTCTAGAACGATCGTTATGTCGAAGTCCAAGCACTACGGGGTAC"+
         "TACAAGTATACCACGAGCCTACTTACCTCCAAA",
        "AAATGATGGTCGCTCCACTTTCTCTCTAGAACGATCGTTATGTCGAAGTCCAAGCACTACGGGGTAC"+
         "TACAAGTATACCACGAGCCTACTTACCTCCAAA",
        "GAATGATGGTCGCTCCACTTTCTCTCTAGAACGATCGTTATGTCGAAGTCCAAGCACTACGGGGTAC"+
         "TACAAGTATACCACGAGCCTACTTACCTCCAAA"
        "AAATGATGGTCGCTCCACTTTCTCTCTAGAACGATCGTTATGTCGAAGTCCGAGCACTACGGGGTAC"+
         "TACAAGTATACCACGAGCCTACTTACCTCCAAT",
        "CAATGATGGTCGCTCCACTTTCTCTCTAGAACGATCGTTATGTCGAAGTCCGAGCACTACGGGGTAC"+
         "TACAAGTATACCACGAGCCTACTTACCTCCAAT",
        "AAATGATGGTCGCTCCACTTTCTCTCTAGAACGATCGTTATGTCGAAGTCCAAGCACTACGGGGTAC"+
         "TACAAGTATACCACGAGCCTACTTACCTCCAAA",
        "CAATGATGGTCGCTCCACTTTCTCTCTAGAACGATCGTTATGTCGAAGTCCAAGCACTACGGGGTAC"+
         "TACAAGTATACCACGAGCCTACTTACCTCCAAT",
        "CAATGATGGTCGCTCCACTTTCTCTCTAGAACGATCGTTATGTCGAAGTCCAAGCACTACGGGGTAC"+
         "TACAAGTATACCACGAGCCTACTTACCTCCAAC",
        "GAATGATGGTCGCTCCACTTTCTCTCTAGAACGATCGTTATGTCGAAGTCCTAGCACTACGGGGTAC"+
         "TACAAGTATACCACGAGCCTACTTACCTCCAAT"]

TEST_PATH = os.path.dirname(os.path.abspath(__file__))


# Test Consensus made by two expected haplotype lists.
class Test_Consensus:
    
    def test_consensus_start(self):
        
        haplotype_list = []
 
        for i in range(len(front)):
 
             if front[i] not in haplotype_list:
                 haplotype_list.append(haplotype.Haplotype(front[i])) 

        expected_consensus =  "AAGAGTTGTGAGGAGTACTCACCTCAGTAGACAAGGAGAGCTATAT"+\
                    "CGCGCAGCATCGCTTCGATATTTTGCTCCTACGCATCCACACGTTGAAAGGGCA"

        
        result = haplotype.build_consensus_from_haplotypes( \
                 haplotypes_list_front)
        assert expected_consensus == result


    def test_consensus_end(self):
        haplotype_list = []
 
        for i in range(len(back)):
 
             if front[i] not in haplotype_list:
                 haplotype_list.append(haplotype.Haplotype(back[i]))

        expected_consensus =  "AAATGATGGTCGCTCCACTTTCTCTCTAGAACGATCGTTATGTCG"+\
                    "AAGTCCAAGCACTACGGGGTACTACAAGTATACCACGAGCCTACTTACCTCCAAA"
        result = haplotype.build_consensus_from_haplotypes( \
                 haplotype_list)
        assert expected_consensus == result

      
# Test to see if bam subcommand runs.
# NOTE : will not produce an output file.
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

# Test to see if fasta subcommand runs.
# NOTE : will not produce an output file.
class Test_FASTA_Complexity:
    @classmethod
    def setup(self):
        self.fasta_location = TEST_PATH + '/data/complexity.fasta'
        self.output_location =  TEST_PATH + '/data/output.csv'


    def test_complexity_fasta(self):
        
        runner = CliRunner()
        result = runner.invoke(complexity.fasta, [self.fasta_location, "--output_location", self.output_location], catch_exceptions=False)
        
        # If method ran successfully the exit code is 0.
        assert result.exit_code == 0
        # output checks for print messages at the end of method.
        # the bam method in complexity has no print message.
        assert result.output == ""

# Test each measurment
class Test_Measurements():
    
    @classmethod
    def setup(self):
        
        self.consensus =  "AAGAGTTGTGAGGAGTACTCACCTCAGTAGACAAGGAGAGCTATAT"+\
                    "CGCGCAGCATCGCTTCGATATTTTGCTCCTACGCATCCACACGTTGAAAGGGCA"

        self.sorted_hap = haplotype.sort_haplotypes(haplotypes_list_front, self.consensus)
        
        self.pileup = haplotype.build_pileup_from_haplotypes(self.sorted_hap)

        self.frequencies =  haplotype.build_frequencies(self.sorted_hap)

        self.distance_matrix = haplotype.build_distiance_matrix(self.sorted_hap)
        
        self.counts = haplotype.build_counts(self.sorted_hap)

    def test_number_of_mutations(self):

       unique_mutations = self.pileup.count_unique_mutations()

       assert unique_mutations == 6

    
    def test_number_of_polymorphic_sites(self):

        number_of_polymorphic_sites = self.pileup.count_polymorphic_sites()

        assert number_of_polymorphic_sites == 3

    def test_shannon_entropy(self):
        
        shannon_entropy = calculate.shannon_entropy(self.frequencies)

        assert float("{0:.2f}".format(shannon_entropy))== 1.23

    def test_shannon_entropy_normalized_to_n(self):

        Hs = 1.23

        Hsn = complexity.get_shannon_entropy_normalized_to_n(self.sorted_hap, Hs)

        assert float("{0:.2f}".format(Hsn)) == 0.53

    def test_shannon_entropy_normalized_to_h(self):

        Hs = 1.23

        
        Hsn = complexity.get_shannon_entropy_normalized_to_h(self.sorted_hap, Hs)

        assert float("{0:.2f}".format(Hsn)) == 0.76 

    
    def test_minimum_muttion_frequency(self):

        minimum_mutation_frequency = complexity.get_minimum_mutation_frequency(self.sorted_hap,self.pileup)
        
        assert minimum_mutation_frequency == 0.006

    def test_mutation_frequency(self):

        mutation_frequency = complexity.get_mutation_frequency(self.distance_matrix)

        assert mutation_frequency == 0.02


    def test_functional_attribute_diversity(self):

        FAD = complexity.get_FAD(self.distance_matrix)
        assert float("{0:.2f}".format(FAD)) == 0.46

    def test_sample_nucleotide_diversity_entity(self):

        SND = complexity.get_sample_nucleotide_diversity(self.distance_matrix, self.frequencies, self.sorted_hap)

        assert float("{0:.2f}".format(SND)) == 0.02

    def test_maximum_mutation(self):

        maximum_mutation_frequency =  complexity.get_maximum_mutation_frequency(self.counts, self.distance_matrix, self.frequencies)
 
        assert float("{0:.2f}".format(maximum_mutation_frequency)) == 0.01

    def test_population_nucleotide_diversity(self):
        
        PND = complexity.get_population_nucleotide_diversity(self.distance_matrix, self.frequencies)

        assert float("{0:.2f}".format(PND)) == 0.01

    def test_sample_nucleotide_diversity_entity(self):

        snde = complexity.get_sample_nucleotide_diversity_entity(self.distance_matrix, self.frequencies)

        assert float("{0:.2f}".format(snde)) == 0.02
       

class Test_CSV_Building:

    @classmethod
    def setup(self):
        self.measurements = []
        self.measurements.append(['0'])
            
    def test_measurements_to_csv(self):

            with open('output_location.csv', 'w')  as \
            complexity_file:
                # Assert nothing because we should return a none type if method runs correctly
                assert None ==  complexity.measurement_to_csv(self.measurements, complexity_file)

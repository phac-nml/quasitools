
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


TEST_PATH = os.path.dirname(os.path.abspath(__file__))


# Test Consensus made by two expected haplotype lists.
class Test_Consensus:
    
    def test_consensus_simple(self):
        
        haplotypes_list = [
                haplotype.Haplotype("AAGAGTTGTGAGGAGTACTCACCTCAGTAGACAAGGAGAGCTA", 6),
                haplotype.Haplotype("AAGAGTTGTGAGGAGTACTCACCTCAGTAGACAAGGAGAGCTA", 1),
                haplotype.Haplotype("TAGAGTTGTGAGGAGTACTCACCTCAGTAGACAAGGAGAGCTA", 1),
                haplotype.Haplotype("TAGAGTTGTGAGGAGTACTCACCTCAGTAGACAAGGAGAGCTA", 1),
                haplotype.Haplotype("CAGAGTTGTGAGGAGTACTCACCTCAGTAGACAAGGAGAGCTA", 1)
                ]

        expected_consensus =  "AAGAGTTGTGAGGAGTACTCACCTCAGTAGACAAGGAGAGCTA"
        result = haplotype.build_consensus_from_haplotypes( \
                 haplotypes_list)

        assert expected_consensus == result

    def test_consensus_degenerate(self):
        
        haplotypes_list = [
                haplotype.Haplotype("AAGAGTTGTGAGGAGTACTCACCTCAGTAGACAAGGAGAGCTA", 1),
                haplotype.Haplotype("AAGAGTTGTGAGGAGTACTCACCTCAGTAGACAAGGAGAGCTA", 1),
                haplotype.Haplotype("AAGAGTTGTGAGGAGTACTCACCTCAGTAGACAAGGAGAGCTA", 1),
                haplotype.Haplotype("TAGAGTTGTGAGGAGTACTCACCTCAGTAGACAAGGAGAGCTA", 1),
                haplotype.Haplotype("CAGAGTTGTGAGGAGTACTCACCTCAGTAGACAAGGAGAGCTA", 1)
                ]

        
        expected_consensus =  "AAGAGTTGTGAGGAGTACTCACCTCAGTAGACAAGGAGAGCTA"
        result = haplotype.build_consensus_from_haplotypes( \
                 haplotypes_list)

        assert expected_consensus == result

    def test_consensus_tie(self):
        
        haplotypes_list = [
                haplotype.Haplotype("AAGAGTTGTGAGGAGTACTCACCTCAGTAGACAAGGAGAGCTA", 6),
                haplotype.Haplotype("AAGAGTTGTGAGGAGTACTCACCTCAGTAGACAAGGAGAGCTA", 1),
                haplotype.Haplotype("TAGAGTTGTGAGGAGTACTCACCTCAGTAGACAAGGAGAGCTA", 1),
                haplotype.Haplotype("TAGAGTTGTGAGGAGTACTCACCTCAGTAGACAAGGAGAGCTA", 1),
                haplotype.Haplotype("CAGAGTTGTGAGGAGTACTCACCTCAGTAGACAAGGAGAGCTA", 6)
                ]
        expected_consensus =  "AAGAGTTGTGAGGAGTACTCACCTCAGTAGACAAGGAGAGCTA"        
        result = haplotype.build_consensus_from_haplotypes( \
                 haplotypes_list)

        assert expected_consensus == result

    def test_consensus_missing(self):
        
        haplotypes_list = []

        result = haplotype.build_consensus_from_haplotypes( \
                 haplotypes_list)

        expected_consensus = "" 

        assert expected_consensus == result

      
# Test to see if bam subcommand runs.

class Test_BAM_Complexity:
    @classmethod
    def setup(self):
        self.bam_location = TEST_PATH + '/data/complexity.bam'
        self.reference_location =  TEST_PATH + '/data/complexity_reference.fasta'
        self.output_location_bam =  TEST_PATH + '/data/output_bam.csv'

    def test_complexity_bam(self):
        
        runner = CliRunner()
        result = runner.invoke(complexity.bam, [self.reference_location,\
                self.bam_location, "1",'--output_location', \
                self.output_location_bam])
        
        # If method ran successfully the exit code is 0.
        assert result.exit_code == 0
        # output checks for print messages at the end of method.
        # the bam method in complexity has no print message.
        assert result.output == ""

        # Check if output file is created
        assert os.path.exists(TEST_PATH + '/data/output_bam.csv') == True
        
        # Check to see if expected values are found in csv
        # file that we created. We will look at the last row.
        for line in open(TEST_PATH + '/data/output_bam.csv'):
            # split each row into a list
            csv_row = line.split(',')
        
        # Check if last row has these values for each column
        assert csv_row[0] == '199'
        assert csv_row[1] == '3'
        assert csv_row[2] == '10'
        assert csv_row[3] == '1'
        assert csv_row[4] == '2'

# Test to see if fasta subcommand runs.
class Test_FASTA_Complexity:
    @classmethod
    def setup(self):
        self.fasta_location = TEST_PATH + '/data/complexity.fasta'
        self.output_location_fasta =  TEST_PATH + '/data/output_fasta.csv'


    def test_complexity_fasta(self):
        
        runner = CliRunner()
        result = runner.invoke(complexity.fasta, [self.fasta_location, "--output_location", self.output_location_fasta], catch_exceptions=False)
        
        # If method ran successfully the exit code is 0.
        assert result.exit_code == 0
        # output checks for print messages at the end of method.
        # the bam method in complexity has no print message.
        assert result.output == ""

        # Check if output file is created
        assert os.path.exists(TEST_PATH + '/data/output_fasta.csv') == True
        
        # Check to see if expected values are found in csv
        # file that we created. We will look at the last row.
        for line in open(TEST_PATH + '/data/output_fasta.csv'):
            csv_row = line.split(',')

        # Check if last row has these values for each column
        assert csv_row[0] == '0'
        assert csv_row[1] == '2'
        assert csv_row[2] == '2'
        assert csv_row[3] == '7'
        assert csv_row[4] == '7'


# Test each measurment
class Test_Measurements():
    
    @classmethod
    def setup(self):
        
        self.consensus =  "AAGAGTTGTGAGGAGTACTCACCTCAGTAGACAAGGAGAGCTA"


        self.haplotypes_list = [
            haplotype.Haplotype("AAGAGTTGTGAGGAGTACTCACCTCAGTAGACAAGGAGAGCTA", 6),
            haplotype.Haplotype("AAGAGTTGTGAGGAGTACTCACCTCAGTAGACAAGGAGAGCTA", 1),
            haplotype.Haplotype("TAGAGTTGTGAGGAGTACTCACCTCAGTAGACAAGGAGAGCTA", 1),
            haplotype.Haplotype("TAGAGTTGTGAGGAGTACTCACCTCAGTAGACAAGGAGAGCTA", 1),
            haplotype.Haplotype("CAGAGTTGTGAGGAGTACTCACCTCAGTAGACAAGGAGAGCTA", 1)
            ]

        self.sorted_hap = haplotype.sort_haplotypes(self.haplotypes_list, self.consensus)
        
        self.pileup = haplotype.build_pileup_from_haplotypes(self.sorted_hap)

        self.frequencies =  haplotype.build_frequencies(self.sorted_hap)

        self.distance_matrix = haplotype.build_distiance_matrix(self.sorted_hap)
        
        self.counts = haplotype.build_counts(self.sorted_hap)

    def test_number_of_mutations(self):

       unique_mutations = self.pileup.count_unique_mutations()

       assert unique_mutations == 2

    
    def test_number_of_polymorphic_sites(self):

        number_of_polymorphic_sites = self.pileup.count_polymorphic_sites()

        assert number_of_polymorphic_sites == 1

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
        
        assert float("{0:.3f}".format(minimum_mutation_frequency)) == 0.005

    def test_mutation_frequency(self):

        mutation_frequency = complexity.get_mutation_frequency(self.distance_matrix)
        assert float("{0:.2f}".format(mutation_frequency)) == 0.01

    def test_functional_attribute_diversity(self):

        FAD = complexity.get_FAD(self.distance_matrix)
        assert float("{0:.2f}".format(FAD)) == 0.37

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

        self.measurements = [[1,2,3,4,5]]
    

    def test_measurements_to_csv(self):

            file_directory = TEST_PATH + '/data/output_built_csv.csv'

            with open(file_directory, 'w')  as \
            complexity_file:
                complexity.measurement_to_csv(self.measurements, complexity_file)

            assert os.path.exists(file_directory) == True

            # Check to see if expected values are found in csv
            # file that we created. We will look at the last row.
            for line in open(TEST_PATH + '/data/output_built_csv.csv'):
                csv_row = line.split(',')

            # Check if last row has these values for each column
            assert csv_row[0] == '0'
            assert csv_row[1] == '1'
            assert csv_row[2] == '2'
            assert csv_row[3] == '3'
            assert csv_row[4] == '4'
            assert csv_row[5] == '5\r\n'




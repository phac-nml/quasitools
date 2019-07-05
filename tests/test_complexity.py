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
                haplotype.Haplotype("AAAATAAA", 2),
                haplotype.Haplotype("AAAAAAAA", 3),
                haplotype.Haplotype("AAAAAAAA", 4),
                haplotype.Haplotype("GAAAAAAC", 6)               
                ]

        expected_consensus =  "AAAAAAAA"
        result = haplotype.build_consensus_from_haplotypes( \
                 haplotypes_list)

        assert expected_consensus == result

    def test_consensus_degenerate(self):
        haplotypes_list = [
                haplotype.Haplotype("AAAATAAC", 1),
                haplotype.Haplotype("AAAATAAC", 1),
                haplotype.Haplotype("AAAATAAC", 1),
                haplotype.Haplotype("TAAAGAAA", 1)               
                ]
   
        expected_consensus =  "AAAATAAC"
        result = haplotype.build_consensus_from_haplotypes( \
                 haplotypes_list)

        assert expected_consensus == result

    def test_consensus_tie(self):
        haplotypes_list = [
                haplotype.Haplotype("AAAATAAC", 6),
                haplotype.Haplotype("AAAATAAC", 3),
                haplotype.Haplotype("AAAATAAC", 4),
                haplotype.Haplotype("TAAAGAAA", 13)               
                ]
        
        expected_consensus =  "AAAAGAAA"        
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
        self.bam_location2 = TEST_PATH + '/data/complexity_data.bam'
        self.output_location_bam =  TEST_PATH + '/data/output_bam.csv'
        self.output_location_bam_filter =  TEST_PATH + '/data/output_bam2.csv'
        self.output_location_bam_filter2 =  TEST_PATH + '/data/output_bam3.csv'
        self.output_location_bam_filter3 =  TEST_PATH + '/data/output_bam4.csv'
        self.output_location_bam_filter4 =  TEST_PATH + '/data/output_bam5.csv'

    def test_complexity_bam(self):
        
        runner = CliRunner()
        result = runner.invoke(complexity.bam, [self.reference_location,\
                self.bam_location, "1", '--output_location', \
                self.output_location_bam])
        
        # Check if output file is created
        assert os.path.exists(TEST_PATH + '/data/output_bam.csv')
        
        # Check to see if expected values are found in csv
        # file that we created. We will look at the last row.
        for line in open(TEST_PATH + '/data/output_bam.csv'):
            # split each row into a list
            csv_row = line.split(',')
        
        # Check if last row has these values for each column
        assert csv_row[0].strip()  == '199'
        assert csv_row[1].strip()  == '3'
        assert csv_row[2].strip()  == '10'
        assert csv_row[3].strip()  == '1'
        assert csv_row[4].strip()  == '2'
    
    # Test BAM complexity when filter of 25  is applied
    def test_complexity_filter_25(self):

        runner = CliRunner()
        result = runner.invoke(complexity.bam, [self.reference_location,\
                self.bam_location, "1","--haplotype_filter", \
                25, '--output_location', \
                self.output_location_bam_filter])
        
        # Check if output file is created
        assert os.path.exists(TEST_PATH + '/data/output_bam2.csv')
        
        # Check to see if expected values are found in csv
        # file that we created. We will look at the last row.
        for line in open(TEST_PATH + '/data/output_bam2.csv'):
            # split each row into a list
            csv_row = line.split(',')

        assert csv_row[0].strip()  == '199'
        assert csv_row[1].strip()  == '2'
        assert csv_row[2].strip()  == '9'
        assert csv_row[3].strip()  == '1'
        assert csv_row[4].strip()  == '1'
    
    # Test Complexity filter with a filter of 0 (i.e no filtering, everything should pass)
    def test_compexity_filter_0(self):

        runner = CliRunner()
        result = runner.invoke(complexity.bam, [self.reference_location,\
                self.bam_location, "1","--haplotype_filter", \
                0, '--output_location', \
                self.output_location_bam_filter2])
        
        # Check if output file is created
        assert os.path.exists(TEST_PATH + '/data/output_bam3.csv')
        
        # Check to see if expected values are found in csv
        # file that we created. We will look at the last row.
        for line in open(TEST_PATH + '/data/output_bam3.csv'):
            # split each row into a list
            csv_row = line.split(',')

        assert csv_row[0].strip()  == '199'
        assert csv_row[1].strip()  == '3'
        assert csv_row[2].strip()  == '10'
        assert csv_row[3].strip()  == '1'
        assert csv_row[4].strip()  == '2'

   # Test the complexity filter with a value of 100.
   # Since the BAM file has multiple haplotypes at every position, nothing should pass the filter.
    def test_compexity_filter_100(self):

        runner = CliRunner()
        result = runner.invoke(complexity.bam, [self.reference_location,\
                self.bam_location, "1","--haplotype_filter", \
                100, '--output_location', \
                self.output_location_bam_filter3])
        
        # Check if output file is created
        assert os.path.exists(TEST_PATH + '/data/output_bam4.csv')
        
        # Check to see if expected values are found in csv
        # file that we created. We will look at the last row.
        for line in open(TEST_PATH + '/data/output_bam4.csv'):
            # split each row into a list
            csv_row = line.split(',')

        assert csv_row[0].strip()  == '199'
        assert csv_row[1].strip()  == ''
        assert csv_row[2].strip()  == ''
        assert csv_row[3].strip()  == ''
        assert csv_row[4].strip()  == ''

    # Test the complexity filter with a value of 100.
    # Since the BAM file has only one haplotype at every position, everything should pass filter.
    def test_compexity_filer_100_2(self):
        runner = CliRunner()
        result = runner.invoke(complexity.bam, [self.reference_location,\
                self.bam_location2, "1","--haplotype_filter", \
                100, '--output_location', \
                self.output_location_bam_filter4])
        
        # Check if output file is created:
        assert os.path.exists(TEST_PATH + '/data/output_bam5.csv')
        
        # Check to see if the expected values are found in the CSV 
        # file that we created. We will only look at the last row of this file.
        for line in open(TEST_PATH + '/data/output_bam5.csv'):
            # split each row into a list
            csv_row = line.split(',')
        
        # This is just the row number, nothing was processed.
        assert csv_row[0].strip()  == '199'
        
        # Row numbers
        assert csv_row[1].strip()  == '1'
        assert csv_row[2].strip()  == '1'
        assert csv_row[3].strip()  == '0'
        assert csv_row[4].strip()  == '0'



# Test to see if the fasta subcommand runs.
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
        # Check if output file is created
        assert os.path.exists(TEST_PATH + '/data/output_fasta.csv')
        
        # Check to see if expected values are found in csv
        # file that we created. We will look at the last row.
        for line in open(TEST_PATH + '/data/output_fasta.csv'):
            csv_row = line.split(',')

        # Check if last row has these values for each column
        assert csv_row[0].strip()  == '0'
        assert csv_row[1].strip()  == '2'
        assert csv_row[2].strip()  == '2'
        assert csv_row[3].strip()  == '7'
        assert csv_row[4].strip()  == '7'


# Test each measurment
class Test_Measurements():
    
    @classmethod
    def setup(self):
        
        self.consensus =  "AAGAGTTGTGAGGAGTACTCACCTCAGTAGACAAGGAGAGCTA"

        self.haplotypes_list = [
            haplotype.Haplotype("AAGAGTTGTGAGGAGTACTCACCTCAGTAGACAAGGAGAGCTA", 1),
            haplotype.Haplotype("AAGAGTTGTGAGGAGTACTCACCTCAGTAGACAAGGAGAGCTA", 1),
            haplotype.Haplotype("GAGAGTTGTGAGGAGTACTCACCTCAGTAGACAAGGAGAGCTA", 1),
            haplotype.Haplotype("TAGAGTTGTGAGGAGTACTCACCTCAGTAGACAAGGAGAGCTA", 6),
            haplotype.Haplotype("CAGAGTTGTGAGGAGTACTCACCTCAGTAGACAAGGAGAGCTA", 6)
            ]

        self.sorted_hap = haplotype.sort_haplotypes(self.haplotypes_list, self.consensus)
        
        self.pileup = haplotype.build_pileup_from_haplotypes(self.sorted_hap)

        self.frequencies =  haplotype.build_frequencies(self.sorted_hap)

        self.distance_matrix = haplotype.build_distiance_matrix(self.sorted_hap)
        
        self.counts = haplotype.build_counts(self.sorted_hap)

    def test_number_of_mutations(self):

       unique_mutations = self.pileup.count_unique_mutations()

       assert unique_mutations == 3

    
    def test_number_of_polymorphic_sites(self):

        number_of_polymorphic_sites = self.pileup.count_polymorphic_sites()

        assert number_of_polymorphic_sites == 1

    def test_shannon_entropy(self):
        
        shannon_entropy = calculate.shannon_entropy(self.frequencies)

        assert float("{0:.2f}".format(shannon_entropy))== 1.27

    def test_shannon_entropy_normalized_to_n(self):

        Hs = 1.23

        Hsn = complexity.get_shannon_entropy_normalized_to_n(self.sorted_hap, Hs)

        assert float("{0:.2f}".format(Hsn)) == 0.45

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
        assert float("{0:.2f}".format(FAD)) == 0.42

    def test_sample_nucleotide_diversity_entity(self):

        SND = complexity.get_sample_nucleotide_diversity(self.distance_matrix, self.frequencies, self.sorted_hap)

        assert float("{0:.2f}".format(SND)) == 0.02

    def test_maximum_mutation(self):

        maximum_mutation_frequency =  complexity.get_maximum_mutation_frequency(self.counts, self.distance_matrix, self.frequencies)
 
        assert float("{0:.2f}".format(maximum_mutation_frequency)) == 0.02

    def test_population_nucleotide_diversity(self):
        
        PND = complexity.get_population_nucleotide_diversity(self.distance_matrix, self.frequencies)

        assert float("{0:.2f}".format(PND)) == 0.02

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

            assert os.path.exists(file_directory)

            # Check to see if expected values are found in csv
            # file that we created. We will look at the last row.
            for line in open(TEST_PATH + '/data/output_built_csv.csv'):
                csv_row = line.split(',')

            # Check if last row has these values for each column
            assert csv_row[0].strip() == '0'
            assert csv_row[1].strip() == '1'
            assert csv_row[2].strip() == '2'
            assert csv_row[3].strip() == '3'
            assert csv_row[4].strip() == '4'
            assert csv_row[5].strip() == '5'


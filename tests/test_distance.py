"""
Copyright Government of Canada 2018

Written by: Matthew Fogel, Public Health Agency of Canada

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
from quasitools.parsers.mapped_read_parser import parse_mapped_reads_from_bam
from quasitools.parsers.reference_parser import parse_references_from_fasta
from quasitools.distance import Distance

TEST_PATH = os.path.dirname(os.path.abspath(__file__))

class TestDistance:
    @classmethod
    def setup_class(self):
        self.dist = Distance()
        #files for get_distance_matrix test function
        self.test_gdm_files = ('test1.bam', 'test2.bam', 'test3.bam', 'test4.bam',
                            'test5.bam', 'test6.bam', 'test7.bam', 'test8.bam')
        self.pileup = [[{'A': 1, 'T': 1, 'C': 1, 'G': 1}, {'A': 1, 'T': 1, 'C': 1},
                   {'T': 12}, {'C': 12}, {'G': 12}, {'A': 12}, {'T': 12}, {'C': 12}, {'G': 12}, {'A': 12}], #test1
                  [{'A': 1, 'T': 1, 'C': 1, 'G': 1}, {'A': 1, 'T': 1, 'C': 1000000},
                   {'T': 12}, {'C': 12}, {'G': 12}, {'A': 12}, {'T': 12}, {'C': 12}, {'G': 12}, {'A': 12}], #test2
                  [{'A': 1, 'T': 1, 'C': 1, 'G': 1}, {'A': 1, 'T': 1000000, 'C': 1},
                   {'T': 12}, {'C': 12}, {'G': 12}, {'A': 12}, {'T': 12}, {'C': 12}, {'G': 12}, {'A': 12}], #test3
                  [{'A': 1, 'T': 1, 'C': 1, 'G': 1}, {'A': 1000000, 'T': 1, 'C': 1},
                   {'T': 12}, {'C': 12}, {'G': 12}, {'A': 12}, {'T': 12}, {'C': 12}, {'G': 12}, {'A': 12}], #test4
                  [{'A': 1, 'T': 1, 'C': 1, 'G': 1}, {'A': 1000000, 'T': 1, 'C': 1000000},
                   {'T': 12}, {'C': 12}, {'G': 12}, {'A': 12}, {'T': 12}, {'C': 12}, {'G': 12}, {'A': 12}], #test5
                  [{'A': 1, 'T': 1, 'C': 1, 'G': 1}, {'A': 1000000, 'T': 1000000, 'C': 1},
                   {'T': 12}, {'C': 12}, {'G': 12}, {'A': 12}, {'T': 12}, {'C': 12}, {'G': 12}, {'A': 12}], #test6
                  [{'A': 1, 'T': 1, 'C': 1, 'G': 1}, {'A': 1, 'T': 1000000, 'C': 1000000},
                   {'T': 12}, {'C': 12}, {'G': 12}, {'A': 12}, {'T': 12}, {'C': 12}, {'G': 12}, {'A': 12}], #test7
                  [{'A': 1, 'T': 1, 'C': 1, 'G': 1}, {'A': 1000000, 'T': 1000000, 'C': 1000000},
                   {'T': 12}, {'C': 12}, {'G': 12}, {'A': 12}, {'T': 12}, {'C': 12}, {'G': 12}, {'A': 12}]] #test8

        self.pileup2 = [[{'A': 1, 'T': 1, 'C': 1, 'G': 1}, {'A': 1, 'T': 1, 'C': 1},
                   {'T': 12}, {'C': 12}, {'G': 12}, {'A': 12}, {'T': 12}, {'C': 12}, {'G': 12}, {'A': 12}], #test1
                  [{'A': 1, 'T': 1, 'C': 1, 'G': 1}, {'A': 1, 'T': 1, 'C': 1000000},
                   {'T': 12}, {'C': 12}, {'G': 12}, {'A': 12}, {'T': 12}, {'C': 12}, {'G': 12}, {'A': 12}], #test2
                  [{'A': 1, 'T': 1, 'C': 1, 'G': 1}, {'A': 1, 'T': 1000000, 'C': 1},
                   {'T': 12}, {'C': 12}, {'G': 12}, {'A': 12}, {'T': 12}, {'C': 12}, {'G': 12}, {'A': 12}], #test3
                  [{'A': 1, 'T': 1, 'C': 1, 'G': 1}, {'A': 1000000, 'T': 1, 'C': 1},
                   {'T': 12}, {'C': 12}, {'G': 12}, {'A': 12}, {'T': 12}, {'C': 12}, {'G': 12}, {'A': 12}], #test4
                  [{'A': 1, 'T': 1, 'C': 1, 'G': 1}, {'A': 1000000, 'T': 1, 'C': 1000000},
                   {'T': 12}, {'C': 12}, {'G': 12}, {'A': 12}, {'T': 12}, {'C': 12}, {'G': 12}, {'A': 12}], #test5
                  [{'A': 1, 'T': 1, 'C': 1, 'G': 1}, {'A': 1000000, 'T': 1000000, 'C': 1},
                   {'T': 12}, {'C': 12}, {'G': 12}, {'A': 12}, {'T': 12}, {'C': 12}, {'G': 12}, {'A': 12}], #test6
                  [{'A': 1, 'T': 1, 'C': 1, 'G': 1}, {'A': 1, 'T': 1000000, 'C': 1000000},
                   {'T': 12}, {'C': 12}, {'G': 12}, {'A': 12}, {'T': 12}, {'C': 12}, {'G': 12}, {'A': 12}], #test7
                  [{'A': 1, 'T': 1, 'C': 1, 'G': 1}, {'A': 1000000, 'T': 1000000, 'C': 1000000},
                   {'T': 12}, {'C': 12}, {'G': 12}, {'A': 12}, {'T': 12}, {'C': 12}, {'G': 12}, {'A': 12}]] #test8

        self.ones = [[{'A': 1, 'T': 1, 'C': 1}, {'T': 1}], #test 1
                      [{'A': 1, 'T': 1, 'C': 1}, {'T': 1}]] # test 2

        self.test_ones_files = ('test1.bam', 'test2.bam')

        #files for construct pileup test function
        self.test_cp_files = ((TEST_PATH+"/data/quasi1.bam"), (TEST_PATH+"/data/quasi2.bam"))
        self.test_cp_ref = TEST_PATH+"/data/hxb2_pol.fas"

    #end def

    def test_get_distance_matrix(self):
        normalized = self.dist.get_distance_matrix(self.pileup, 'normalize')
        assert(len(normalized) == 8)
        assert(len(normalized[0]) == 8)
        assert(normalized[0][7] == 1)
        assert(self.dist.convert_distance_to_csv(normalized, self.test_gdm_files) ==
        ("Quasispecies,test1.bam,test2.bam,test3.bam,test4.bam,test5.bam,test6.bam,test7.bam,test8.bam\n" +
"test1.bam,1.00000000,0.96329037,0.96329037,0.96329037,0.99043043,0.99043043,0.99043043,1.00000000\n" +
"test2.bam,0.96329037,1.00000000,0.89189249,0.89189249,0.97259768,0.91702091,0.97259768,0.96329037\n" +
"test3.bam,0.96329037,0.89189249,1.00000000,0.89189249,0.91702091,0.97259768,0.97259768,0.96329037\n" +
"test4.bam,0.96329037,0.89189249,0.89189249,1.00000000,0.97259768,0.97259768,0.91702091,0.96329037\n" +
"test5.bam,0.99043043,0.97259768,0.91702091,0.97259768,1.00000000,0.97142866,0.97142866,0.99043043\n" +
"test6.bam,0.99043043,0.91702091,0.97259768,0.97259768,0.97142866,1.00000000,0.97142866,0.99043043\n" +
"test7.bam,0.99043043,0.97259768,0.97259768,0.91702091,0.97142866,0.97142866,1.00000000,0.99043043\n" +
"test8.bam,1.00000000,0.96329037,0.96329037,0.96329037,0.99043043,0.99043043,0.99043043,1.00000000")
)

        dnormal = self.dist.get_distance_matrix(self.pileup2, None)
        assert(len(dnormal) == 8)
        assert(len(dnormal[0]) == 8)
        #assert(unnormalized[0][7] == 0.05089630)
        assert(self.dist.convert_distance_to_csv(dnormal, self.test_gdm_files) ==
        ("Quasispecies,test1.bam,test2.bam,test3.bam,test4.bam,test5.bam,test6.bam,test7.bam,test8.bam\n" +
"test1.bam,1.00000000,0.02940769,0.02940769,0.02940769,0.04156468,0.04156468,0.04156468,0.05089630\n" +
"test2.bam,0.02940769,1.00000000,0.00000200,0.00000200,0.70710749,0.00000212,0.70710749,0.57735142\n" +
"test3.bam,0.02940769,0.00000200,1.00000000,0.00000200,0.00000212,0.70710749,0.70710749,0.57735142\n" +
"test4.bam,0.02940769,0.00000200,0.00000200,1.00000000,0.70710749,0.70710749,0.00000212,0.57735142\n" +
"test5.bam,0.04156468,0.70710749,0.00000212,0.70710749,1.00000000,0.50000100,0.50000100,0.81649699\n" +
"test6.bam,0.04156468,0.00000212,0.70710749,0.70710749,0.50000100,1.00000000,0.50000100,0.81649699\n" +
"test7.bam,0.04156468,0.70710749,0.70710749,0.00000212,0.50000100,0.50000100,1.00000000,0.81649699\n" +
"test8.bam,0.05089630,0.57735142,0.57735142,0.57735142,0.81649699,0.81649699,0.81649699,1.00000000")
)

        matrixTwo = self.dist.get_distance_matrix(self.ones, 'normalize')
        assert(len(matrixTwo) == 2)
        assert(len(matrixTwo[0]) == 2)
        assert(matrixTwo[0][1] == 1)
        assert(matrixTwo[1][0] == 1)

    #end def

    def test_construct_pileup(self):
        bamPileup = self.dist.construct_pileup(self.test_cp_files, self.test_cp_ref)
        assert(len(bamPileup)==2)
        assert(len(bamPileup[0])==2844)
        assert(len(bamPileup[1])==2844)
        assert(bamPileup[0][0:10] == [{'C': 12}, {'C': 12}, {'T': 12}, {'C': 12}, {'G': 2, 'C': 3, 'T': 1, 'A': 6}, {'G': 12}, {'G': 12}, {'T': 12}, {'C': 12}, {'G': 2, 'C': 3, 'T': 1, 'A': 6}])
        assert(bamPileup[1][0:10] == [{'C': 12}, {'C': 12}, {'T': 12}, {'C': 12}, {'A': 6, 'C': 5, 'G': 1}, {'G': 12}, {'A': 4, 'G': 8}, {'T': 12}, {'C': 12}, {'A': 7, 'T': 1, 'C': 3, 'G': 1}])

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
        self.test_gdm_files = ['test1.bam', 'test2.bam', 'test3.bam', 'test4.bam',
                            'test5.bam', 'test6.bam', 'test7.bam', 'test8.bam']
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

        self.ones = [[{'A': 1, 'T': 1, 'C': 1}, {'T': 1}], #test 1
                      [{'A': 1, 'T': 1, 'C': 1}, {'T': 1}]] # test 2

        self.test_ones_files = ['test1.bam', 'test2.bam']

        self.pileup2 = [[{'A': 3, 'T': 3, 'C': 2, 'G': 4}, {'T': 6, 'C': 6},
                   {'T': 12}, {'C': 12}, {'G': 12}, {'A': 12}, {'T': 12}, {'C': 12}, {'G': 12}, {'A': 12}], #test1
                  [{'A': 3, 'T': 3, 'C': 3, 'G': 3}, {'T': 6, 'C': 6},
                   {'T': 12}, {'C': 12}, {'G': 12}, {'A': 12}, {'T': 12}, {'C': 12}, {'G': 12}, {'A': 12}], #test2
                  [{'A': 3, 'T': 3, 'C': 2, 'G': 4}, {'A': 1, 'T': 1, 'C': 10},
                   {'T': 12}, {'C': 12}, {'G': 12}, {'A': 12}, {'T': 12}, {'C': 12}, {'G': 12}, {'A': 12}], #test3
                  [{'A': 3, 'T': 3, 'C': 2, 'G': 4}, {'C': 12},
                   {'T': 12}, {'C': 12}, {'G': 12}, {'A': 12}, {'T': 12}, {'C': 12}, {'G': 12}, {'A': 12}], #test4
                  [{'A': 3, 'T': 3, 'C': 2, 'G': 4}, {'G': 6, 'C': 6},
                   {'T': 12}, {'C': 12}, {'G': 12}, {'A': 12}, {'T': 12}, {'C': 12}, {'G': 12}, {'A': 12}], #test5
                  [{'A': 3, 'T': 3, 'C': 2, 'G': 4}, {'T': 6, 'G': 6},
                   {'T': 12}, {'C': 12}, {'G': 12}, {'A': 12}, {'T': 12}, {'C': 12}, {'G': 12}, {'A': 12}], #test6
                  [{'A': 3, 'T': 3, 'C': 2, 'G': 4}, {'G': 6, 'A': 6},
                   {'T': 12}, {'C': 12}, {'G': 12}, {'A': 12}, {'T': 12}, {'C': 12}, {'G': 12}, {'A': 12}]] #test7

        #files for construct pileup test function
        self.test_cp_files = ((TEST_PATH+"/data/quasi1.bam"), (TEST_PATH+"/data/quasi2.bam"))
        self.test_cp_ref = TEST_PATH+"/data/hxb2_pol.fas"

    #end def

    def test_get_distance_matrix(self):
        normalized = self.dist.get_distance_matrix(self.pileup, self.test_gdm_files, 'normalize', None, None)
        assert(len(normalized) == 9)
        assert(len(normalized[0]) == 9)
        assert(normalized[1][8] == 1) #TODO: verify whether this works
        assert(normalized == [['Quasispecies', 'test1.bam', 'test2.bam', 'test3.bam', 'test4.bam', 'test5.bam', 'test6.bam', 'test7.bam', 'test8.bam'], ['test1.bam', 0.0, 0.963290365981529, 0.963290365981529, 0.963290365981529, 0.9904304301700023, 0.9904304301700023, 0.9904304301700023, 1.0], ['test2.bam', 0.963290365981529, 0.0, 0.891892493789242, 0.891892493789242, 0.9725976798731226, 0.9170209149268688, 0.9725976798731226, 0.963290365981529], ['test3.bam', 0.963290365981529, 0.891892493789242, 0.0, 0.891892493789242, 0.9170209149268688, 0.9725976798731226, 0.9725976798731226, 0.963290365981529], ['test4.bam', 0.963290365981529, 0.891892493789242, 0.891892493789242, 0.0, 0.9725976798731226, 0.9725976798731226, 0.9170209149268688, 0.963290365981529], ['test5.bam', 0.9904304301700023, 0.9725976798731226, 0.9170209149268688, 0.9725976798731226, 0.0, 0.9714286555101037, 0.9714286555101037, 0.9904304301700023], ['test6.bam', 0.9904304301700023, 0.9170209149268688, 0.9725976798731226, 0.9725976798731226, 0.9714286555101037, 0.0, 0.9714286555101037, 0.9904304301700023], ['test7.bam', 0.9904304301700023, 0.9725976798731226, 0.9725976798731226, 0.9170209149268688, 0.9714286555101037, 0.9714286555101037, 0.0, 0.9904304301700023], ['test8.bam', 1.0, 0.963290365981529, 0.963290365981529, 0.963290365981529, 0.9904304301700023, 0.9904304301700023, 0.9904304301700023, 0.0]])

        matrixTwo = self.dist.get_distance_matrix(self.ones, self.test_ones_files, 'normalize', None, None)
        assert(len(matrixTwo) == 3)
        assert(len(matrixTwo[0]) == 3)
        assert(matrixTwo[1][2] == 1)
        assert(matrixTwo[2][1] == 1)

    #end def

    def test_construct_pileup(self):
        bamPileup = self.dist.construct_pileup(self.test_cp_files, self.test_cp_ref)
        assert(len(bamPileup)==2)
        assert(len(bamPileup[0])==2844)
        assert(len(bamPileup[1])==2844)
        assert(bamPileup[0][0:10] == [{'C': 12}, {'C': 12}, {'T': 12}, {'C': 12}, {'G': 2, 'C': 3, 'T': 1, 'A': 6}, {'G': 12}, {'G': 12}, {'T': 12}, {'C': 12}, {'G': 2, 'C': 3, 'T': 1, 'A': 6}])
        assert(bamPileup[1][0:10] == [{'C': 12}, {'C': 12}, {'T': 12}, {'C': 12}, {'A': 6, 'C': 5, 'G': 1}, {'G': 12}, {'A': 4, 'G': 8}, {'T': 12}, {'C': 12}, {'A': 7, 'T': 1, 'C': 3, 'G': 1}])

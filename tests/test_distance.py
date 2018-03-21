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

    """
    CLASS VARIABLES
    """

    #files for testing pileup matrix
    pileup1 = [[{'A': 1, 'T': 1, 'C': 1, 'G': 1}, {'A': 1, 'T': 1, 'C': 1},
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

    pileup1_files = ('test1.bam', 'test2.bam', 'test3.bam', 'test4.bam',
                          'test5.bam', 'test6.bam', 'test7.bam', 'test8.bam')

    #expected output files for pileup1

    pileup1_normal_out = """Quasispecies,test1.bam,test2.bam,test3.bam,test4.bam,test5.bam,test6.bam,test7.bam,test8.bam\n
    test1.bam,1.00000000,0.96329037,0.96329037,0.96329037,0.99043043,0.99043043,0.99043043,1.00000000\n
    test2.bam,0.96329037,1.00000000,0.89189249,0.89189249,0.97259768,0.91702091,0.97259768,0.96329037\n
    test3.bam,0.96329037,0.89189249,1.00000000,0.89189249,0.91702091,0.97259768,0.97259768,0.96329037\n
    test4.bam,0.96329037,0.89189249,0.89189249,1.00000000,0.97259768,0.97259768,0.91702091,0.96329037\n
    test5.bam,0.99043043,0.97259768,0.91702091,0.97259768,1.00000000,0.97142866,0.97142866,0.99043043\n
    test6.bam,0.99043043,0.91702091,0.97259768,0.97259768,0.97142866,1.00000000,0.97142866,0.99043043\n
    test7.bam,0.99043043,0.97259768,0.97259768,0.91702091,0.97142866,0.97142866,1.00000000,0.99043043\n
    test8.bam,1.00000000,0.96329037,0.96329037,0.96329037,0.99043043,0.99043043,0.99043043,1.00000000"""

    pileup1_unnormal_out = """Quasispecies,test1.bam,test2.bam,test3.bam,test4.bam,test5.bam,test6.bam,test7.bam,test8.bam\n
    test1.bam,1.00000000,0.02940769,0.02940769,0.02940769,0.04156468,0.04156468,0.04156468,0.05089630\n
    test2.bam,0.02940769,1.00000000,0.00000200,0.00000200,0.70710749,0.00000212,0.70710749,0.57735142\n
    test3.bam,0.02940769,0.00000200,1.00000000,0.00000200,0.00000212,0.70710749,0.70710749,0.57735142\n
    test4.bam,0.02940769,0.00000200,0.00000200,1.00000000,0.70710749,0.70710749,0.00000212,0.57735142\n
    test5.bam,0.04156468,0.70710749,0.00000212,0.70710749,1.00000000,0.50000100,0.50000100,0.81649699\n
    test6.bam,0.04156468,0.00000212,0.70710749,0.70710749,0.50000100,1.00000000,0.50000100,0.81649699\n
    test7.bam,0.04156468,0.70710749,0.70710749,0.00000212,0.50000100,0.50000100,1.00000000,0.81649699\n
    test8.bam,0.05089630,0.57735142,0.57735142,0.57735142,0.81649699,0.81649699,0.81649699,1.00000000"""

    #files for testing pileup2 matrix of ones
    pileup2 = ([[{'A': 1, 'T': 1, 'C': 1}, {'T': 1}], #test 1
                [{'A': 1, 'T': 1, 'C': 1}, {'T': 1}]]) # test 2

    pileup2_files = ('test1.bam', 'test2.bam')

    pileup2_normal_out = """Quasispecies,test1.bam,test2.bam\n
    test1.bam,1.00000000,1.00000000\n
    test2.bam,1.00000000,1.00000000"""

    pileup2_unnormal_out = """Quasispecies,test1.bam,test2.bam\n
    test1.bam,1.00000000,1.00000000\n
    test2.bam,1.00000000,1.00000000"""

    tuple_list = [(True, pileup1, pileup1_files, pileup1_normal_out),
    (True, pileup2, pileup2_files, pileup2_normal_out),
    (False, pileup1, pileup1_files, pileup1_unnormal_out),
    (False, pileup2, pileup2_files, pileup2_unnormal_out)]

    #files for testing pileup matrix
    pileup3 = [[{'A': 1}, {'T': 2}, {'C': 3}, {'G': 4}, {'A': 5}, {'T': 6}, {'C': 7}, {'G': 8}], #test 1
    [{}, {'T': 2}, {'C': 3}, {'G': 4}, {'A': 5}, {'T': 6}, {'C': 7}, {'G': 8}], #test 2
    [{'A': 1}, {'T': 2}, {'C': 3}, {'G': 4}, {'A': 5}, {'T': 6}, {'C': 7}, {'G': 8}], #test 3
    [{'A': 1}, {'T': 2}, {'C': 3}, {'G': 4}, {'A': 5}, {'T': 6}, {'C': 7}, {'G': 8}], #test 4
    [{'A': 1}, {'T': 2}, {'C': 3}, {}, {'A': 5}, {'T': 6}, {'C': 7}, {'G': 8}], #test 5
    [{'A': 1}, {'T': 2}, {'C': 3}, {'G': 4}, {'A': 5}, {'T': 6}, {'C': 7}, {'G': 8}], #test 6
    [{'A': 1}, {'T': 2}, {'C': 3}, {'G': 4}, {'A': 5}, {'T': 6}, {}, {'G': 8}], #test 7
    [{'A': 1}, {'T': 2}, {'C': 3}, {'G': 4}, {'A': 5}, {'T': 6}, {'C': 7}, {'G': 8}]] #test 8

    pileup3_truncate_out = [[{'T': 2}, {'C': 3}, {'A': 5}, {'T': 6}, {'G': 8}],
                            [{'T': 2}, {'C': 3}, {'A': 5}, {'T': 6}, {'G': 8}],
                            [{'T': 2}, {'C': 3}, {'A': 5}, {'T': 6}, {'G': 8}],
                            [{'T': 2}, {'C': 3}, {'A': 5}, {'T': 6}, {'G': 8}],
                            [{'T': 2}, {'C': 3}, {'A': 5}, {'T': 6}, {'G': 8}],
                            [{'T': 2}, {'C': 3}, {'A': 5}, {'T': 6}, {'G': 8}],
                            [{'T': 2}, {'C': 3}, {'A': 5}, {'T': 6}, {'G': 8}],
                            [{'T': 2}, {'C': 3}, {'A': 5}, {'T': 6}, {'G': 8}]]

    pileup4 = [[{'A': 1, 'T': 0, 'C': 0, 'G': 0},  # test 1 a
                {'A': 0, 'T': 2, 'C': 0, 'G': 0},  # test 1 b
                {'A': 0, 'T': 0, 'C': 0, 'G': 0}],  # test 1 c
               [{'A': 0, 'T': 0, 'C': 0, 'G': 0},  # test 2 a
                {'A': 0, 'T': 2, 'C': 0, 'G': 0},  # test 2 b
                {'A': 0, 'T': 0, 'C': 0, 'G': 0}]]  # test 2 c

    pileup4_truncate_out = [[{'A': 0, 'T': 2, 'C': 0, 'G': 0}], #test 1
                            [{'A': 0, 'T': 2, 'C': 0, 'G': 0}]]

    truncate_list = [(pileup3, pileup3_truncate_out),
                     (pileup4, pileup4_truncate_out)]

    """
    TESTS
    """

    @pytest.fixture(scope="function", params=tuple_list)
    def matrix_fixture(self, request):
        """
        matrix_fixture - test fixture for the test_get_distance_matrix function

        INPUT:
            [LIST OF TUPLES]
            ---[BOOL] [normalize] # normalized or not
            ---[ARRAY] [pileup list]
            ---[ARRAY] [pileup_files] # file names corresponding to pileups
            ---[ARRAY] [pileup#_(normal_out/unnormal_out)] #csv formatted output

        RETURN:
            Tuple containing: (actual csv-formatted string output,
            expected csv-formatted string output)

        POST:
            [None]
        """
        dist = Distance()
        #save expected output (request.param[3])
        expected = request.param[3]

        #get similarity matrix based on pileup list (request.param[1])
        #and boolean normalize flag (request.param[0])
        matrix = dist.get_distance_matrix(request.param[1], request.param[0])

        #convert matrix  to csv, passing file list (request.param[2])
        csv_similarity = dist.convert_distance_to_csv(matrix, request.param[2])

        return (csv_similarity, request.param[3])

    #end def

    def test_get_distance_matrix(self, matrix_fixture):
        """
        test_get_distance_matrix - Checked that the actual output matches the
        expected output.

        INPUT:
            [FIXTURE] [matrix_fixture]

        RETURN:
            [None]

        POST:
            [None]
        """

        assert matrix_fixture[0] == matrix_fixture[1]
    #end def

    @pytest.fixture
    def pileup_fixture(self):
        """
        pileup_fixture - Passes file locations to construct_pileup function to
        calculate the bam pileup.

        INPUT:
            [NONE]

        RETURN:
            [Array] [bamPileup]

        POST:
            [None]
        """
        #files for construct pileup test function
        test_cp_files = ((TEST_PATH+"/data/quasi1.bam"), (TEST_PATH+"/data/quasi2.bam"))
        test_cp_ref = TEST_PATH+"/data/hxb2_pol.fas"
        dist = Distance()
        bamPileup = dist.construct_pileup(test_cp_files, test_cp_ref)
        return bamPileup
    #end def

    def test_construct_pileup(self, pileup_fixture):
        """
        test_construct_pileup - Checks that the pileup length and the first few
        indices of the pileup are correct.

        INPUT:
            [FIXTURE] [pileup_fixture]

        RETURN:
            [None]

        POST:
            [None]
        """

        assert len(pileup_fixture)==2
        assert len(pileup_fixture[0])==2844
        assert len(pileup_fixture[1])==2844
        assert pileup_fixture[0][0:10] == [{'C': 12}, {'C': 12}, {'T': 12}, {'C': 12}, {'G': 2, 'C': 3, 'T': 1, 'A': 6}, {'G': 12}, {'G': 12}, {'T': 12}, {'C': 12}, {'G': 2, 'C': 3, 'T': 1, 'A': 6}]
        assert pileup_fixture[1][0:10] == [{'C': 12}, {'C': 12}, {'T': 12}, {'C': 12}, {'A': 6, 'C': 5, 'G': 1}, {'G': 12}, {'A': 4, 'G': 8}, {'T': 12}, {'C': 12}, {'A': 7, 'T': 1, 'C': 3, 'G': 1}]
    #end def

    @pytest.fixture(params=truncate_list)
    def truncate_fixture(self, request):
        """
        truncate_fixture - Truncates output and passes result to test function to
        see if truncated output is as expected.

        INPUT:
            [NONE]

        RETURN:
            [Array] [truncated]

        POST:
            [None]
        """
        dist2 = Distance()
        truncated = dist2.truncate_output(request.param[0])
        return (truncated, request.param[1])
    #end def

    def test_truncate_output(self, truncate_fixture):
        """
        test_truncate_output - Checks that when some positions in the pileup are
        empty, the output is truncated correctly.

        INPUT:
            [ARRAY] [pileup_list]

        RETURN:
            [None]

        POST:
            [None]
        """

        assert truncate_fixture[0] == truncate_fixture[1]

    #end def

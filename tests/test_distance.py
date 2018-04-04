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
from quasitools.distance import Pileup_List
from quasitools.distance import DistanceMatrix

TEST_PATH = os.path.dirname(os.path.abspath(__file__))

class TestDistance:

    """
    CLASS VARIABLES
    """

    pileup0 = [[{'A': 1}, {'T': 2}, {'C': 3}, {'G': 4}, {'A': 5}], #test 1
    [{'A': 1}, {'T': 2}, {'C': 3}, {'G': 4}, {}], #test 2
    [{'A': 1}, {'T': 2}, {'C': 3}, {}, {'A': 5}], #test 3
    [{'A': 1}, {'T': 2}, {}, {'G': 4}, {}], #test 4
    [{'A': 1}, {}, {'C': 3}, {'G': 4}, {'A': 5}]] #test 5

    pileup0_files = ('test1.bam', 'test2.bam', 'test3.bam', 'test4.bam', 'test5.bam')

    pileup0_truncate_out = [[{'A': 1}], #test 1
    [{'A': 1}], #test 2
    [{'A': 1}], #test 3
    [{'A': 1}], #test 4
    [{'A': 1}]] #test 5

    pileup0_normal_out = """Quasispecies,test1.bam,test2.bam,test3.bam,test4.bam,test5.bam
test1.bam,1.00000000,1.00000000,1.00000000,1.00000000,1.00000000
test2.bam,1.00000000,1.00000000,1.00000000,1.00000000,1.00000000
test3.bam,1.00000000,1.00000000,1.00000000,1.00000000,1.00000000
test4.bam,1.00000000,1.00000000,1.00000000,1.00000000,1.00000000
test5.bam,1.00000000,1.00000000,1.00000000,1.00000000,1.00000000"""

    pileup0_unnormal_out = """Quasispecies,test1.bam,test2.bam,test3.bam,test4.bam,test5.bam
test1.bam,1.00000000,1.00000000,1.00000000,1.00000000,1.00000000
test2.bam,1.00000000,1.00000000,1.00000000,1.00000000,1.00000000
test3.bam,1.00000000,1.00000000,1.00000000,1.00000000,1.00000000
test4.bam,1.00000000,1.00000000,1.00000000,1.00000000,1.00000000
test5.bam,1.00000000,1.00000000,1.00000000,1.00000000,1.00000000"""

    pileup0_startpos = 0
    pileup0_endpos = 0

    #files for testing pileup matrix
    pileup1 = [[{'A': 1, 'T': 1, 'C': 1, 'G': 1}, {'A': 1, 'T': 1, 'C': 1}, {'T': 12}, {'C': 12}, {'G': 12}, {'A': 12}, {'T': 12}, {'C': 12}, {'G': 12}, {'A': 12}], #test1
    [{'A': 1, 'T': 1, 'C': 1, 'G': 1}, {'A': 1, 'T': 1, 'C': 1000000}, {'T': 12}, {'C': 12}, {'G': 12}, {'A': 12}, {'T': 12}, {'C': 12}, {'G': 12}, {'A': 12}], #test2
    [{'A': 1, 'T': 1, 'C': 1, 'G': 1}, {'A': 1, 'T': 1000000, 'C': 1}, {'T': 12}, {'C': 12}, {'G': 12}, {'A': 12}, {'T': 12}, {'C': 12}, {'G': 12}, {'A': 12}], #test3
    [{'A': 1, 'T': 1, 'C': 1, 'G': 1}, {'A': 1000000, 'T': 1, 'C': 1}, {'T': 12}, {'C': 12}, {'G': 12}, {'A': 12}, {'T': 12}, {'C': 12}, {'G': 12}, {'A': 12}], #test4
    [{'A': 1, 'T': 1, 'C': 1, 'G': 1}, {'A': 1000000, 'T': 1, 'C': 1000000}, {'T': 12}, {'C': 12}, {'G': 12}, {'A': 12}, {'T': 12}, {'C': 12}, {'G': 12}, {'A': 12}], #test5
    [{'A': 1, 'T': 1, 'C': 1, 'G': 1}, {'A': 1000000, 'T': 1000000, 'C': 1}, {'T': 12}, {'C': 12}, {'G': 12}, {'A': 12}, {'T': 12}, {'C': 12}, {'G': 12}, {'A': 12}], #test6
    [{'A': 1, 'T': 1, 'C': 1, 'G': 1}, {'A': 1, 'T': 1000000, 'C': 1000000}, {'T': 12}, {'C': 12}, {'G': 12}, {'A': 12}, {'T': 12}, {'C': 12}, {'G': 12}, {'A': 12}], #test7
    [{'A': 1, 'T': 1, 'C': 1, 'G': 1}, {'A': 1000000, 'T': 1000000, 'C': 1000000}, {'T': 12}, {'C': 12}, {'G': 12}, {'A': 12}, {'T': 12}, {'C': 12}, {'G': 12}, {'A': 12}]] #test8

    pileup1_files = ('test1.bam', 'test2.bam', 'test3.bam', 'test4.bam', 'test5.bam', 'test6.bam', 'test7.bam', 'test8.bam')

    #expected output files for pileup1

    pileup1_normal_out = """Quasispecies,test1.bam,test2.bam,test3.bam,test4.bam,test5.bam,test6.bam,test7.bam,test8.bam
test1.bam,1.00000000,0.96329037,0.96329037,0.96329037,0.99043043,0.99043043,0.99043043,1.00000000
test2.bam,0.96329037,1.00000000,0.89189249,0.89189249,0.97259768,0.91702091,0.97259768,0.96329037
test3.bam,0.96329037,0.89189249,1.00000000,0.89189249,0.91702091,0.97259768,0.97259768,0.96329037
test4.bam,0.96329037,0.89189249,0.89189249,1.00000000,0.97259768,0.97259768,0.91702091,0.96329037
test5.bam,0.99043043,0.97259768,0.91702091,0.97259768,1.00000000,0.97142866,0.97142866,0.99043043
test6.bam,0.99043043,0.91702091,0.97259768,0.97259768,0.97142866,1.00000000,0.97142866,0.99043043
test7.bam,0.99043043,0.97259768,0.97259768,0.91702091,0.97142866,0.97142866,1.00000000,0.99043043
test8.bam,1.00000000,0.96329037,0.96329037,0.96329037,0.99043043,0.99043043,0.99043043,1.00000000"""

    pileup1_unnormal_out = """Quasispecies,test1.bam,test2.bam,test3.bam,test4.bam,test5.bam,test6.bam,test7.bam,test8.bam
test1.bam,1.00000000,0.02940769,0.02940769,0.02940769,0.04156468,0.04156468,0.04156468,0.05089630
test2.bam,0.02940769,1.00000000,0.00000200,0.00000200,0.70710749,0.00000212,0.70710749,0.57735142
test3.bam,0.02940769,0.00000200,1.00000000,0.00000200,0.00000212,0.70710749,0.70710749,0.57735142
test4.bam,0.02940769,0.00000200,0.00000200,1.00000000,0.70710749,0.70710749,0.00000212,0.57735142
test5.bam,0.04156468,0.70710749,0.00000212,0.70710749,1.00000000,0.50000100,0.50000100,0.81649699
test6.bam,0.04156468,0.00000212,0.70710749,0.70710749,0.50000100,1.00000000,0.50000100,0.81649699
test7.bam,0.04156468,0.70710749,0.70710749,0.00000212,0.50000100,0.50000100,1.00000000,0.81649699
test8.bam,0.05089630,0.57735142,0.57735142,0.57735142,0.81649699,0.81649699,0.81649699,1.00000000"""

    pileup1_normal_angular_distance_out = """Quasispecies,test1.bam,test2.bam,test3.bam,test4.bam,test5.bam,test6.bam,test7.bam,test8.bam
test1.bam,0.00000000,0.17303053,0.17303053,0.17303053,0.08814309,0.08814309,0.08814309,0.00000000
test2.bam,0.17303053,0.00000000,0.29875524,0.29875524,0.14937762,0.26117362,0.14937762,0.17303053
test3.bam,0.17303053,0.29875524,0.00000000,0.29875524,0.26117362,0.14937762,0.14937762,0.17303053
test4.bam,0.17303053,0.29875524,0.29875524,0.00000000,0.14937762,0.14937762,0.26117362,0.17303053
test5.bam,0.08814309,0.14937762,0.26117362,0.14937762,0.00000000,0.15254569,0.15254569,0.08814309
test6.bam,0.08814309,0.26117362,0.14937762,0.14937762,0.15254569,0.00000000,0.15254569,0.08814309
test7.bam,0.08814309,0.14937762,0.14937762,0.26117362,0.15254569,0.15254569,0.00000000,0.08814309
test8.bam,0.00000000,0.17303053,0.17303053,0.17303053,0.08814309,0.08814309,0.08814309,0.00000000"""

    pileup1_unnormal_angular_distance_out = """Quasispecies,test1.bam,test2.bam,test3.bam,test4.bam,test5.bam,test6.bam,test7.bam,test8.bam
test1.bam,0.00000000,0.98127578,0.98127578,0.98127578,0.97353148,0.97353148,0.97353148,0.96758440
test2.bam,0.98127578,0.00000000,0.99999873,0.99999873,0.49999936,0.99999865,0.49999936,0.60817255
test3.bam,0.98127578,0.99999873,0.00000000,0.99999873,0.99999865,0.49999936,0.49999936,0.60817255
test4.bam,0.98127578,0.99999873,0.99999873,0.00000000,0.49999936,0.49999936,0.99999865,0.60817255
test5.bam,0.97353148,0.49999936,0.99999865,0.49999936,0.00000000,0.66666593,0.66666593,0.39182610
test6.bam,0.97353148,0.99999865,0.49999936,0.49999936,0.66666593,0.00000000,0.66666593,0.39182610
test7.bam,0.97353148,0.49999936,0.49999936,0.99999865,0.66666593,0.66666593,0.00000000,0.39182610
test8.bam,0.96758440,0.60817255,0.60817255,0.60817255,0.39182610,0.39182610,0.39182610,0.00000000"""

    #files for testing pileup2 matrix of ones
    pileup2 = ([[{'A': 1, 'T': 1, 'C': 1}, {'T': 1}], #test 1
                [{'A': 1, 'T': 1, 'C': 1}, {'T': 1}]]) # test 2

    pileup2_files = ('test1.bam', 'test2.bam')

    pileup2_normal_out = """Quasispecies,test1.bam,test2.bam
test1.bam,1.00000000,1.00000000
test2.bam,1.00000000,1.00000000"""

    pileup2_unnormal_out = """Quasispecies,test1.bam,test2.bam
test1.bam,1.00000000,1.00000000
test2.bam,1.00000000,1.00000000"""

    similarity_tuple_list = [(True, pileup0, pileup0_files, pileup0_normal_out, pileup0_startpos, pileup0_endpos),
    (True, pileup1, pileup1_files, pileup1_normal_out, None, None),
    (True, pileup2, pileup2_files, pileup2_normal_out, None, None),
    (False, pileup0, pileup0_files, pileup0_unnormal_out, pileup0_startpos, pileup0_endpos),
    (False, pileup1, pileup1_files, pileup1_unnormal_out, None, None),
    (False, pileup2, pileup2_files, pileup2_unnormal_out, None, None)]

    angular_distance_tuple_list = [(True, pileup1, pileup1_files, pileup1_normal_angular_distance_out, None, None),
    (False, pileup1, pileup1_files, pileup1_unnormal_angular_distance_out, None, None)]

    #files for testing pileup matrix
    pileup3 = [[{'A': 1}, {'T': 2}, {'C': 3}, {'G': 4}, {'A': 5}, {'T': 6}, {'C': 7}, {'G': 8}], #test 1
    [{}, {'T': 2}, {'C': 3}, {'G': 4}, {'A': 5}, {'T': 6}, {'C': 7}, {'G': 8}], #test 2
    [{'A': 1}, {'T': 2}, {'C': 3}, {'G': 4}, {'A': 5}, {'T': 6}, {'C': 7}, {'G': 8}], #test 3
    [{'A': 1}, {'T': 2}, {'C': 3}, {'G': 4}, {'A': 5}, {'T': 6}, {'C': 7}, {'G': 8}], #test 4
    [{'A': 1}, {'T': 2}, {'C': 3}, {}, {'A': 5}, {'T': 6}, {'C': 7}, {'G': 8}], #test 5
    [{'A': 1}, {'T': 2}, {'C': 3}, {'G': 4}, {'A': 5}, {'T': 6}, {'C': 7}, {'G': 8}], #test 6
    [{'A': 1}, {'T': 2}, {'C': 3}, {'G': 4}, {'A': 5}, {'T': 6}, {}, {'G':  8}], #test 7
    [{'A': 1}, {'T': 2}, {'C': 3}, {'G': 4}, {'A': 5}, {'T': 6}, {'C': 7}, {'G': 8}]] #test 8

    pileup3_trunc_all_out = [[{'T': 2}, {'C': 3}, {'A': 5}, {'T': 6}, {'G': 8}], #test 1
    [{'T': 2}, {'C': 3}, {'A': 5}, {'T': 6}, {'G': 8}], #test 2
    [{'T': 2}, {'C': 3}, {'A': 5}, {'T': 6}, {'G': 8}], #test 3
    [{'T': 2}, {'C': 3}, {'A': 5}, {'T': 6}, {'G': 8}], #test 4
    [{'T': 2}, {'C': 3}, {'A': 5}, {'T': 6}, {'G': 8}], #test 5
    [{'T': 2}, {'C': 3}, {'A': 5}, {'T': 6}, {'G': 8}], #test 6
    [{'T': 2}, {'C': 3}, {'A': 5}, {'T': 6}, {'G':  8}], #test 7
    [{'T': 2}, {'C': 3}, {'A': 5}, {'T': 6}, {'G': 8}]] #test 8

    pileup3_truncate_start = 1
    pileup3_truncate_end = 7

    pileup4 = [[{'A': 1, 'T': 0, 'C': 0, 'G': 0},  # test 1 a
                {'A': 0, 'T': 2, 'C': 0, 'G': 0},  # test 1 b
                {'A': 0, 'T': 0, '-': 4, 'G': 0}],  # test 1 c
               [{'-': 1, 'T': 0, 'C': 0, 'G': 0},  # test 2 a
                {'A': 0, 'T': 2, 'C': 0, 'G': 0},  # test 2 b
                {'A': 0, '-': 5, 'C': 0, 'G': 0}]]  # test 2 c

    pileup4_trunc_ends_out = [[{'A': 0, 'T': 2, 'C': 0, 'G': 0}], #test 1
                            [{'A': 0, 'T': 2, 'C': 0, 'G': 0}]]

    pileup4_truncate_start = 1
    pileup4_truncate_end = 1

    pileup5 = [[{'A': 1}, {'T': 2}, {'C': 3}, {'G': 4}, {'A': 5}], #test 1
    [{}, {'T': 2}, {'C': 3}, {'G': 4}, {'A': 5}], #test 2
    [{'A': 1}, {'T': 2}, {}, {'G': 4}, {'A': 5}], #test 3
    [{'A': 1}, {'T': 2}, {'C': 3}, {'G': 4}, {}], #test 4
    [{'A': 1}, {'T': 2}, {'C': 3}, {}, {'A': 5}]] #test 5

    pileup5_trunc_ends_out = [[{'T': 2}],
                            [{'T': 2}],
                            [{'T': 2}],
                            [{'T': 2}],
                            [{'T': 2}]]

    pileup5_truncate_start = 1
    pileup5_truncate_end = 1

    pileup6 = [[{'A': 1}, {'T': 2}, {'C': 3}, {'G': 4}, {'A': 5}], #test 1
    [{}, {'T': 2}, {'C': 3}, {'G': 4}, {'A': 5}], #test 2
    [{'A': 1}, {}, {'C': 3}, {'G': 4}, {'A': 5}], #test 3
    [{'A': 1}, {'T': 2}, {'C': 3}, {'G': 4}, {}], #test 4
    [{'A': 1}, {'T': 2}, {'C': 3}, {}, {'A': 5}]] #test 5

    pileup6_trunc_ends_out = [[{'C': 3}],
                            [{'C': 3}],
                            [{'C': 3}],
                            [{'C': 3}],
                            [{'C': 3}]]

    pileup6_truncate_start = 2
    pileup6_truncate_end = 2

    pileup7 = [[{'A': 1}, {'T': 2}, {'C': 3}, {'G': 4}, {'A': 5}], #test 1
    [{'A': 1}, {'T': 2}, {'C': 3}, {'G': 4}, {'A': 5}], #test 2
    [{'A': 1}, {'T': 2}, {'C': 3}, {'G': 4}, {'A': 5}], #test 3
    [{'A': 1}, {'T': 2}, {'C': 3}, {'G': 4}, {'A': 5}], #test 4
    [{'A': 1}, {'T': 2}, {'C': 3}, {}, {'A': 5}]] #test 5

    pileup7_trunc_ends_out = [[{'A': 1}, {'T': 2}, {'C': 3}, {'G': 4}, {'A': 5}], #test 1
    [{'A': 1}, {'T': 2}, {'C': 3}, {'G': 4}, {'A': 5}], #test 2
    [{'A': 1}, {'T': 2}, {'C': 3}, {'G': 4}, {'A': 5}], #test 3
    [{'A': 1}, {'T': 2}, {'C': 3}, {'G': 4}, {'A': 5}], #test 4
    [{'A': 1}, {'T': 2}, {'C': 3}, {}, {'A': 5}]] #test 5

    pileup7_truncate_start = 0
    pileup7_truncate_end = 4

    pileup8 = [[{'A': 1}, {'T': 2}, {'C': 3}, {'G': 4}, {'A': 5}], #test 1
    [{'A': 1}, {}, {'C': 3}, {'G': 4}, {'A': 5}], #test 2
    [{'A': 1}, {'T': 2}, {'C': 3}, {'G': 4}, {'A': 5}], #test 3
    [{'A': 1}, {'T': 2}, {'C': 3}, {'G': 4}, {'A': 5}], #test 4
    [{'A': 1}, {'T': 2}, {'C': 3}, {'G': 4}, {'A': 5}]] #test 5

    pileup8_trunc_ends_out = [[{'A': 1}, {'T': 2}, {'C': 3}, {'G': 4}, {'A': 5}], #test 1
    [{'A': 1}, {}, {'C': 3}, {'G': 4}, {'A': 5}], #test 2
    [{'A': 1}, {'T': 2}, {'C': 3}, {'G': 4}, {'A': 5}], #test 3
    [{'A': 1}, {'T': 2}, {'C': 3}, {'G': 4}, {'A': 5}], #test 4
    [{'A': 1}, {'T': 2}, {'C': 3}, {'G': 4}, {'A': 5}]] #test 5

    pileup8_truncate_start = 0
    pileup8_truncate_end = 4

    pileup9 = [[{}, {'T': 2}, {'C': 3}, {'G': 4}, {'A': 5}], #test 1
    [{'A': 1}, {'T': 2}, {'C': 3}, {'G': 4}, {}], #test 2
    [{}, {'T': 2}, {'C': 3}, {'G': 4}, {'A': 5}], #test 3
    [{}, {'T': 2}, {'C': 3}, {'G': 4}, {}], #test 4
    [{'A': 1}, {'T': 2}, {'C': 3}, {'G': 4}, {'A': 5}]] #test 5

    pileup9_trunc_ends_out = [[{'T': 2}, {'C': 3}, {'G': 4}], #test 1
    [{'T': 2}, {'C': 3}, {'G': 4}], #test 2
    [{'T': 2}, {'C': 3}, {'G': 4}], #test 3
    [{'T': 2}, {'C': 3}, {'G': 4}], #test 4
    [{'T': 2}, {'C': 3}, {'G': 4}]] #test 5

    pileup9_truncate_start = 1
    pileup9_truncate_end = 3

    pileup10 = [[{'A': 1}, {'T': 2}, {'C': 3}, {'G': 4}, {'A': 5}], #test 1
    [{'A': 1}, {'T': 2}, {'C': 3}, {'G': 4}, {}], #test 2
    [{'A': 1}, {'T': 2}, {'C': 3}, {}, {'A': 5}], #test 3
    [{'A': 1}, {'T': 2}, {}, {'G': 4}, {}], #test 4
    [{'A': 1}, {}, {'C': 3}, {'G': 4}, {'A': 5}]] #test 5

    pileup10_trunc_ends_out = [[{'A': 1}], #test 1
    [{'A': 1}], #test 2
    [{'A': 1}], #test 3
    [{'A': 1}], #test 4
    [{'A': 1}]] #test 5

    pileup10_truncate_start = 0
    pileup10_truncate_end = 0

    truncate_all_tuple = [(pileup3, pileup3_trunc_all_out)]

    truncate_ends_tuple = [(pileup4, pileup4_trunc_ends_out),
            (pileup5, pileup5_trunc_ends_out),
            (pileup6, pileup6_trunc_ends_out),
            (pileup7, pileup7_trunc_ends_out),
            (pileup8, pileup8_trunc_ends_out),
            (pileup9, pileup9_trunc_ends_out),
            (pileup10, pileup10_trunc_ends_out)]

    """
    TESTS
    """

    @pytest.fixture(scope="function", params=similarity_tuple_list)
    def similarity_matrix_fixture(self, request):
        """
        similarity_matrix_fixture - test fixture for the test_get_similarity_matrix function

        INPUT:
            [LIST OF TUPLES]
            ---[BOOL] [normalize] # normalized or not
            ---[ARRAY] [pileup list]
            ---[ARRAY] [pileup_files] # file names corresponding to pileups
            ---[ARRAY] [pileup#_(normal_out/unnormal_out)] #csv formatted output
            ---[INT or NONE] [startpos or default if NONE]
            ---[INT or NONE] [endpos or default if NONE]

        RETURN:
            Tuple containing: (actual csv-formatted string output,
            expected csv-formatted string output)

        POST:
            [None]
        """
        # if boolean normalize flag (request.param[0]) is true normalize
        pileup_util = Pileup_List(request.param[1])

        # if startpos is int and endpos is int (aka they are not None)
        if type(request.param[4]) is int and type(request.param[5]) is int:
            pileup_util.select_pileup_range(request.param[4],
                                            request.param[5])
        # end if

        if request.param[0] is True:
            pileup_util.normalize_pileup()

        # create matrix with pileup
        dist = DistanceMatrix(pileup_util.get_pileup())

        # get similarity matrix based on pileup list (request.param[1])
        matrix = dist.get_cosine_similarity_matrix()

        # convert matrix  to csv, passing file list (request.param[2])
        csv_similarity = dist.get_similarity_matrix_as_csv(request.param[2])

        dist = None

        return (csv_similarity, request.param[3])

    #end def

    def test_get_similarity_matrix(self, similarity_matrix_fixture):
        """
        test_get_similarity_matrix - Checked that the actual output matches the
        expected output.

        INPUT:
            [FIXTURE] [matrix_fixture]

        RETURN:
            [None]

        POST:
            [None]
        """

        assert similarity_matrix_fixture[0] == similarity_matrix_fixture[1]
    #end def

    @pytest.fixture(scope="function", params=angular_distance_tuple_list)
    def distance_matrix_fixture(self, request):
        """
        distance_matrix_fixture - test fixture for the test_get_distance_matrix function

        INPUT:
            [LIST OF TUPLES]
            ---[BOOL] [normalize] # normalized or not
            ---[ARRAY] [pileup list]
            ---[ARRAY] [pileup_files] # file names corresponding to pileups
            ---[ARRAY] [pileup#_(normal_out/unnormal_out)] #csv formatted output
            ---[INT or NONE] [startpos or default if NONE]
            ---[INT or NONE] [endpos or default if NONE]

        RETURN:
            Tuple containing: (actual csv-formatted string output,
            expected csv-formatted string output)

        POST:
            [None]
        """
        # if boolean normalize flag (request.param[0]) is true normalize
        pileup_util = Pileup_List(request.param[1])

        # if startpos is int and endpos is int (aka they are not None)
        if type(request.param[4]) is int and type(request.param[5]) is int:
            pileup_util.select_pileup_range(request.param[4],
                                            request.param[5])
        # end if

        if request.param[0] is True:
            pileup_util.normalize_pileup()

        # create matrix with pileup
        dist = DistanceMatrix(pileup_util.get_pileup())

        # get similarity matrix based on pileup list (request.param[1])

        matrix = dist.get_angular_cosine_distance_matrix()

        # convert matrix  to csv, passing file list (request.param[2])
        csv_similarity = dist.get_distance_matrix_as_csv(request.param[2])

        dist = None

        return (csv_similarity, request.param[3])

    #end def

    def test_get_distance_matrix(self, distance_matrix_fixture):
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

        assert distance_matrix_fixture[0] == distance_matrix_fixture[1]
    #end def

    @pytest.fixture
    def pileup_fixture(self):
        """
        pileup_fixture - Passes file locations to construct_array_of_pileups function to
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
        bamPileup = Pileup_List.construct_array_of_pileups(test_cp_files, test_cp_ref)
        return bamPileup.get_pileup()
    #end def

    def test_construct_array_of_pileup(self, pileup_fixture):
        """
        test_construct_array_of_pileups - Checks that the pileup length and the first few
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

    @pytest.fixture(scope="function",
                    params=truncate_ends_tuple)
    def test_truncate_output_fixture(self, request):
        """
        test_truncate_output_fixture - calls truncateoutput and returns
        a tuple that can be used for testing this function.

        INPUT:
            [TUPLE] [request.param] - a tuple containing:
            --pileup to be truncated
            --expected truncated pileup
            --expected truncated start position
            --expected truncated end position

        RETURN:
            [TUPLE]
            --tuple containing actual truncated pileup, start pos, end pos
            --expected truncated pileup

        POST:
            [None]
        """
        util = Pileup_List(request.param[0])
        util.truncate_output()
        truncated = util.get_pileup()

        return (truncated, request.param[1])

    def test_truncate_output(self, test_truncate_output_fixture):
        """
        test_truncate_output - Checks that the after truncating contiguous
        beginning and ending zero coverage positions from the pileup that the
        output is as expected.
        INPUT:
            [TUPLE]
            --tuple containing actual truncated pileup, start pos, end pos
            --expected truncated pileup

        RETURN:
            [None]

        POST:
            Checks that the expected outputs match the actual output
        """

        assert test_truncate_output_fixture[0] == test_truncate_output_fixture[1]

    @pytest.fixture(scope="function",
                    params=truncate_all_tuple)
    def test_truncate_all_output_fixture(self, request):
        """
        test_truncate_all_output_fixture - calls truncate_all_output and returns
        a tuple that can be used for testing this function.

        INPUT:
            [TUPLE] [request.param] - a tuple containing:
            --pileup to be truncated
            --expected truncated pileup

        RETURN:
            [TUPLE]
            --tuple containing actual truncated pileup
            --expected truncated pileup

        POST:
            [None]
        """
        util = Pileup_List(request.param[0])
        util.truncate_all_output()
        truncated = util.get_pileup()

        return (truncated, request.param[1])

    def test_truncate_all_output(self, test_truncate_all_output_fixture):
        """
        test_truncate_all_output - Checks that the after truncating all
        empty positions from the pileup that the output is as expected.

        INPUT:
            [TUPLE]
            --tuple containing actual truncated pileup, start pos, end pos
            --expected truncated pileup
            --expected truncated start position
            --expected truncated end position

        RETURN:
            [None]

        POST:
            Checks that the expected outputs match the actual output
        """

        pileup = [[{'A': 1}, {'T': 2}, {'C': 3}, {'G': 4}, {'A': 5}], #test 1
        [{'A': 1}, {'T': 2}, {'C': 3}, {'G': 4}, {'A': 5}], #test 2
        [{'A': 1}, {'T': 2}, {'C': 3}, {'G': 4}, {'A': 5}], #test 3
        [{'A': 1}, {'T': 2}, {'C': 3}, {'G': 4}, {'A': 5}], #test 4
        [{'A': 1}, {'T': 2}, {'C': 3}, {}, {}]] #test 5

        startpos = 2
        endpos = 4

        util = Pileup_List(pileup)
        util.select_pileup_range(startpos, endpos)
        util.truncate_all_output()
        truncated = util.get_pileup()

        assert truncated == [[{'C': 3}], #test 1
        [{'C': 3}], #test 2
        [{'C': 3}], #test 3
        [{'C': 3}], #test 4
        [{'C': 3}]]

        startpos = 3
        endpos = 3

        util = Pileup_List(pileup)
        util.select_pileup_range(startpos, endpos)
        util.truncate_all_output()
        truncated = util.get_pileup()

        assert truncated == [[], #test 1
        [], #test 2
        [], #test 3
        [], #test 4
        []]

        assert test_truncate_all_output_fixture[0] == test_truncate_all_output_fixture[1]

    @pytest.fixture(scope="function")
    def truncate_ends_select_range_fixture(self, request):
        """
        truncate_ends_and_select_range_fixture - Checks if the program works
        as expected when truncating contiguous start and end regions after
        first selecting a specified range.

        INPUT:
            [None]
        RETURN:
            TODO
        POST:
            TODO
        """
        pileup11 = [[{}, {'T': 2}, {'C': 3}, {'G': 4}, {'A': 5}], #test 1
        [{'A': 1}, {'T': 2}, {'C': 3}, {'G': 4}, {}], #test 2
        [{}, {'T': 2}, {'C': 3}, {'G': 4}, {'A': 5}], #test 3
        [{}, {'T': 2}, {}, {'G': 4}, {}], #test 4
        [{'A': 1}, {'T': 2}, {'C': 3}, {'G': 4}, {'A': 5}]] #test 5

        user_startpos = 2
        user_endpos = 4

        pileup11_selected_out = [[{'C': 3}, {'G': 4}, {'A': 5}], #test 1
        [{'C': 3}, {'G': 4}, {}], #test 2
        [{'C': 3}, {'G': 4}, {'A': 5}], #test 3
        [{}, {'G': 4}, {}], #test 4
        [{'C': 3}, {'G': 4}, {'A': 5}]] #test 5

        pileup11_trunc_ends_out = [[{'G': 4}], #test 1
        [{'G': 4}], #test 2
        [{'G': 4}], #test 3
        [{'G': 4}], #test 4
        [{'G': 4}]] #test 5

        return (pileup11, user_startpos, user_endpos, pileup11_selected_out, pileup11_trunc_ends_out)

    def test_truncate_ends_select_range(self, truncate_ends_select_range_fixture):
        """
        test_truncate_ends_select_range - tests whether works as expected
        when truncating ends and selecting custom pileup range

        INPUT:
            [truncate_ends_select_range_fixture]
        RETURN:
            [None]
        POST:
            [None]
        """

        pileup = truncate_ends_select_range_fixture[0]
        startpos = truncate_ends_select_range_fixture[1]
        endpos = truncate_ends_select_range_fixture[2]
        pileup_select_range_out = truncate_ends_select_range_fixture[3]
        pileup_truncate_ends_out = truncate_ends_select_range_fixture[4]

        util = Pileup_List(pileup)
        util.select_pileup_range(startpos, endpos)
        select_pileup = util.get_pileup()
        # assert that the pileup positions before startpos and after endpos
        # have been ignored
        assert select_pileup == pileup_select_range_out

        util.truncate_output()
        truncated_pileup = util.get_pileup()

        # assert that the pileup is truncated now as expected
        assert truncated_pileup == pileup_truncate_ends_out

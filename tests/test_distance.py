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
from quasitools.distance import DistanceMatrix
from quasitools.pileup import Pileup_List
from quasitools.pileup import Pileup

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

    pileup0_normal_angular_distance_out = """Quasispecies,test1.bam,test2.bam,test3.bam,test4.bam,test5.bam
test1.bam,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000
test2.bam,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000
test3.bam,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000
test4.bam,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000
test5.bam,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000"""

    pileup0_unnormal_angular_distance_out = """Quasispecies,test1.bam,test2.bam,test3.bam,test4.bam,test5.bam
test1.bam,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000
test2.bam,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000
test3.bam,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000
test4.bam,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000
test5.bam,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000"""

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

    pileup2_normal_angular_distance_out = """Quasispecies,test1.bam,test2.bam
test1.bam,0.00000000,0.00000000
test2.bam,0.00000000,0.00000000"""

    pileup2_unnormal_angular_distance_out = """Quasispecies,test1.bam,test2.bam
test1.bam,0.00000000,0.00000000
test2.bam,0.00000000,0.00000000"""

    """
    TESTS
    """

    @classmethod
    def setup_class(self):
        self.expected_csv_distance = ""
        self.expected_csv_similarity = ""

    @pytest.fixture(scope="function", params=[(True, pileup0, pileup0_files, pileup0_normal_out, pileup0_normal_angular_distance_out, pileup0_startpos, pileup0_endpos),
    (True, pileup1, pileup1_files, pileup1_normal_out, pileup1_normal_angular_distance_out, None, None),
    (True, pileup2, pileup2_files, pileup2_normal_out, pileup2_normal_angular_distance_out, None, None),
    (False, pileup0, pileup0_files, pileup0_unnormal_out, pileup0_unnormal_angular_distance_out, pileup0_startpos, pileup0_endpos),
    (False, pileup1, pileup1_files, pileup1_unnormal_out, pileup1_unnormal_angular_distance_out, None, None),
    (False, pileup2, pileup2_files, pileup2_unnormal_out, pileup2_unnormal_angular_distance_out, None, None)])
    def matrix(self, request):
        """
        matrix - test fixture for test_get_similarity_matrix function
                 and test_get_distance_matrix function

        INPUT:
            [LIST OF TUPLES]
            request.param[0]---[BOOL] [normalize] # normalized or not
            request.param[1]---[ARRAY] [pileup list]
            request.param[2]---[ARRAY] [pileup_files] # file names corresponding to pileups
            request.param[3]---[ARRAY] normalized or unnormalized similarity csv-format output
            request.param[4]---[ARRAY] normalized or unnormalized distance csv-format output
            request.param[5]---[INT or NONE] [startpos or default if NONE]
            request.param[6]---[INT or NONE] [endpos or default if NONE]

        RETURN:
            [DistanceMatrix] [matrix with the pileup to be used]

        POST:
            self.expected_csv_distance is now a csv representation of the
            expected distance that should be calculated from this matrix.

            self.expected_csv_similarity is now a csv representation of the
            expected similarity that should be calculated from this matrix.
        """
        pileups = Pileup_List([Pileup(bam) for bam in request.param[1]])

        # if startpos is int and endpos is int (aka they are not None)
        if type(request.param[5]) is int and type(request.param[6]) is int:
            pileups.select_pileup_range(request.param[5], request.param[6])

        # if boolean normalize flag (request.param[0]) is true normalize
        if request.param[0] is True:
            pileups.normalize_pileup()

        # create matrix with pileup
        dist = DistanceMatrix(pileups.get_pileups_as_numerical_array(), request.param[2])

        self.expected_csv_similarity = request.param[3]
        self.expected_csv_distance = request.param[4]

        return dist

    #end def

    def test_get_similarity_matrix(self, matrix):
        """
        test_get_similarity_matrix - Checked that the actual output matches the
        expected output.

        INPUT:
            [FIXTURE] [matrix] - returns DistanceMatrix object.

        RETURN:
            [None]

        POST:
            [None]
        """
        csv_similarity = matrix.get_similarity_matrix_as_csv()
        assert csv_similarity == self.expected_csv_similarity
    #end def

    def test_get_distance_matrix(self, matrix):
        """
        test_get_distance_matrix - Checked that the actual output matches the
        expected output.

        INPUT:
            [FIXTURE] [matrix] - returns DistanceMatrix object.

        RETURN:
            [None]

        POST:
            [None]
        """
        csv_distance = matrix.get_distance_matrix_as_csv()
        assert csv_distance == self.expected_csv_distance
    #end def

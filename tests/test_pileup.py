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

from quasitools.parsers.mapped_read_parser import parse_mapped_reads_from_bam, parse_pileup_list_from_bam
from quasitools.parsers.reference_parser import parse_references_from_fasta
from quasitools.pileup import Pileup_List
from quasitools.pileup import Pileup
from quasitools import cli
from quasitools.commands.cmd_distance import dist

from click.testing import CliRunner

TEST_PATH = os.path.dirname(os.path.abspath(__file__))

class TestPileups:

    """
    CLASS VARIABLES
    """

    pileup1 = [[{'A': 1}, {'T': 2}, {'C': 3}, {'G': 4}, {'A': 5}], #test 1
    [{'A': 1}, {'T': 2}, {'C': 3}, {'G': 4}, {'A': 5}], #test 2
    [{'A': 1}, {'T': 2}, {'C': 3}, {'G': 4}, {'A': 5}], #test 3
    [{'A': 1}, {'T': 2}, {'C': 3}, {'G': 4}, {'A': 5}], #test 4
    [{'A': 1}, {'T': 2}, {'C': 3}, {}, {}]] #test 5

    pileup1_startpos_1 = 2
    pileup1_endpos_1 = 4

    pileup1_startpos_2 = 3
    pileup1_endpos_2 = 3

    pileup1_truncated_1 = [[{'C': 3}], [{'C': 3}], [{'C': 3}], [{'C': 3}], [{'C': 3}]]

    pileup1_truncated_2 = [[], [], [], [], []]

    pileup2 = [[{}, {'T': 2}, {'C': 3}, {'G': 4}, {'A': 5}], #test 1
    [{'A': 1}, {'T': 2}, {'C': 3}, {'G': 4}, {}], #test 2
    [{}, {'T': 2}, {'C': 3}, {'G': 4}, {'A': 5}], #test 3
    [{}, {'T': 2}, {}, {'G': 4}, {}], #test 4
    [{'A': 1}, {'T': 2}, {'C': 3}, {'G': 4}, {'A': 5}]] #test 5

    pileup2_startpos = 2
    pileup2_endpos = 4

    pileup2_select_range_expected_out = [[{'C': 3}, {'G': 4}, {'A': 5}], #test 1
    [{'C': 3}, {'G': 4}, {}], #test 2
    [{'C': 3}, {'G': 4}, {'A': 5}], #test 3
    [{}, {'G': 4}, {}], #test 4
    [{'C': 3}, {'G': 4}, {'A': 5}]] #test 5

    pileup2_truncate_expected_out = [[{'G': 4}], [{'G': 4}], [{'G': 4}], [{'G': 4}], [{'G': 4}]]

    #files for testing pileup matrix
    pileup3 = [[{'A': 1}, {'T': 2}, {'C': 3}, {'G': 4}, {'A': 5}, {'T': 6}, {'C': 7}, {'G': 8}], #test 1
    [{}, {'T': 2}, {'C': 3}, {'G': 4}, {'A': 5}, {'T': 6}, {'C': 7}, {'G': 8}], #test 2
    [{'A': 1}, {'T': 2}, {'C': 3}, {'G': 4}, {'A': 5}, {'T': 6}, {'C': 7}, {'G': 8}], #test 3
    [{'A': 1}, {'T': 2}, {'C': 3}, {'G': 4}, {'A': 5}, {'T': 6}, {'C': 7}, {'G': 8}], #test 4
    [{'A': 1}, {'T': 2}, {'C': 3}, {}, {'A': 5}, {'T': 6}, {'C': 7}, {'G': 8}], #test 5
    [{'A': 1}, {'T': 2}, {'C': 3}, {'G': 4}, {'A': 5}, {'T': 6}, {'C': 7}, {'G': 8}], #test 6
    [{'A': 1}, {'T': 2}, {'C': 3}, {'G': 4}, {'A': 5}, {'T': 6}, {}, {'G':  8}], #test 7
    [{'A': 1}, {'T': 2}, {'C': 3}, {'G': 4}, {'A': 5}, {'T': 6}, {'C': 7}, {'G': 8}]] #test 8

    pileup3_remove_out = [[{'T': 2}, {'C': 3}, {'A': 5}, {'T': 6}, {'G': 8}], #test 1
    [{'T': 2}, {'C': 3}, {'A': 5}, {'T': 6}, {'G': 8}], #test 2
    [{'T': 2}, {'C': 3}, {'A': 5}, {'T': 6}, {'G': 8}], #test 3
    [{'T': 2}, {'C': 3}, {'A': 5}, {'T': 6}, {'G': 8}], #test 4
    [{'T': 2}, {'C': 3}, {'A': 5}, {'T': 6}, {'G': 8}], #test 5
    [{'T': 2}, {'C': 3}, {'A': 5}, {'T': 6}, {'G': 8}], #test 6
    [{'T': 2}, {'C': 3}, {'A': 5}, {'T': 6}, {'G':  8}], #test 7
    [{'T': 2}, {'C': 3}, {'A': 5}, {'T': 6}, {'G': 8}]] #test 8

    pileup3_num_start_truncated = 1
    pileup3_num_end_truncated = 0

    pileup4 = [[{'A': 1, 'T': 0, 'C': 0, 'G': 0},  # test 1 a
                {'A': 0, 'T': 2, 'C': 0, 'G': 0},  # test 1 b
                {'A': 0, 'T': 0, '-': 4, 'G': 0}],  # test 1 c
               [{'-': 1, 'T': 0, 'C': 0, 'G': 0},  # test 2 a
                {'A': 0, 'T': 2, 'C': 0, 'G': 0},  # test 2 b
                {'A': 0, '-': 5, 'C': 0, 'G': 0}]]  # test 2 c

    pileup4_truncated_out = [[{'A': 0, 'T': 2, 'C': 0, 'G': 0}], #test 1
                            [{'A': 0, 'T': 2, 'C': 0, 'G': 0}]]

    pileup4_num_start_truncated = 1
    pileup4_num_end_truncated = 1

    pileup5 = [[{'A': 1}, {'T': 2}, {'C': 3}, {'G': 4}, {'A': 5}], #test 1
    [{}, {'T': 2}, {'C': 3}, {'G': 4}, {'A': 5}], #test 2
    [{'A': 1}, {'T': 2}, {}, {'G': 4}, {'A': 5}], #test 3
    [{'A': 1}, {'T': 2}, {'C': 3}, {'G': 4}, {}], #test 4
    [{'A': 1}, {'T': 2}, {'C': 3}, {}, {'A': 5}]] #test 5

    pileup5_truncated_out = [[{'T': 2}],
                            [{'T': 2}],
                            [{'T': 2}],
                            [{'T': 2}],
                            [{'T': 2}]]

    pileup5_num_start_truncated = 1
    pileup5_num_end_truncated = 3

    pileup6 = [[{'A': 1}, {'T': 2}, {'C': 3}, {'G': 4}, {'A': 5}], #test 1
    [{}, {'T': 2}, {'C': 3}, {'G': 4}, {'A': 5}], #test 2
    [{'A': 1}, {}, {'C': 3}, {'G': 4}, {'A': 5}], #test 3
    [{'A': 1}, {'T': 2}, {'C': 3}, {'G': 4}, {}], #test 4
    [{'A': 1}, {'T': 2}, {'C': 3}, {}, {'A': 5}]] #test 5

    pileup6_truncated_out = [[{'C': 3}],
                            [{'C': 3}],
                            [{'C': 3}],
                            [{'C': 3}],
                            [{'C': 3}]]

    pileup6_num_start_truncated = 2
    pileup6_num_end_truncated = 2

    pileup7 = [[{'A': 1}, {'T': 2}, {'C': 3}, {'G': 4}, {'A': 5}], #test 1
    [{'A': 1}, {'T': 2}, {'C': 3}, {'G': 4}, {'A': 5}], #test 2
    [{'A': 1}, {'T': 2}, {'C': 3}, {'G': 4}, {'A': 5}], #test 3
    [{'A': 1}, {'T': 2}, {'C': 3}, {'G': 4}, {'A': 5}], #test 4
    [{'A': 1}, {'T': 2}, {'C': 3}, {}, {'A': 5}]] #test 5

    pileup7_truncated_out = [[{'A': 1}, {'T': 2}, {'C': 3}, {'G': 4}, {'A': 5}], #test 1
    [{'A': 1}, {'T': 2}, {'C': 3}, {'G': 4}, {'A': 5}], #test 2
    [{'A': 1}, {'T': 2}, {'C': 3}, {'G': 4}, {'A': 5}], #test 3
    [{'A': 1}, {'T': 2}, {'C': 3}, {'G': 4}, {'A': 5}], #test 4
    [{'A': 1}, {'T': 2}, {'C': 3}, {}, {'A': 5}]] #test 5

    pileup7_num_start_truncated = 0
    pileup7_num_end_truncated = 0

    pileup8 = [[{'A': 1}, {'T': 2}, {'C': 3}, {'G': 4}, {'A': 5}], #test 1
    [{'A': 1}, {}, {'C': 3}, {'G': 4}, {'A': 5}], #test 2
    [{'A': 1}, {'T': 2}, {'C': 3}, {'G': 4}, {'A': 5}], #test 3
    [{'A': 1}, {'T': 2}, {'C': 3}, {'G': 4}, {'A': 5}], #test 4
    [{'A': 1}, {'T': 2}, {'C': 3}, {'G': 4}, {'A': 5}]] #test 5

    pileup8_truncated_out = [[{'A': 1}, {'T': 2}, {'C': 3}, {'G': 4}, {'A': 5}], #test 1
    [{'A': 1}, {}, {'C': 3}, {'G': 4}, {'A': 5}], #test 2
    [{'A': 1}, {'T': 2}, {'C': 3}, {'G': 4}, {'A': 5}], #test 3
    [{'A': 1}, {'T': 2}, {'C': 3}, {'G': 4}, {'A': 5}], #test 4
    [{'A': 1}, {'T': 2}, {'C': 3}, {'G': 4}, {'A': 5}]] #test 5

    pileup8_num_start_truncated = 0
    pileup8_num_end_truncated = 0

    pileup9 = [[{}, {'T': 2}, {'C': 3}, {'G': 4}, {'A': 5}], #test 1
    [{'A': 1}, {'T': 2}, {'C': 3}, {'G': 4}, {}], #test 2
    [{}, {'T': 2}, {'C': 3}, {'G': 4}, {'A': 5}], #test 3
    [{}, {'T': 2}, {'C': 3}, {'G': 4}, {}], #test 4
    [{'A': 1}, {'T': 2}, {'C': 3}, {'G': 4}, {'A': 5}]] #test 5

    pileup9_truncated_out = [[{'T': 2}, {'C': 3}, {'G': 4}], #test 1
    [{'T': 2}, {'C': 3}, {'G': 4}], #test 2
    [{'T': 2}, {'C': 3}, {'G': 4}], #test 3
    [{'T': 2}, {'C': 3}, {'G': 4}], #test 4
    [{'T': 2}, {'C': 3}, {'G': 4}]] #test 5

    pileup9_num_start_truncated = 1
    pileup9_num_end_truncated = 1

    pileup10 = [[{'A': 1}, {'T': 2}, {'C': 3}, {'G': 4}, {'A': 5}], #test 1
    [{'A': 1}, {'T': 2}, {'C': 3}, {'G': 4}, {}], #test 2
    [{'A': 1}, {'T': 2}, {'C': 3}, {}, {'A': 5}], #test 3
    [{'A': 1}, {'T': 2}, {}, {'G': 4}, {}], #test 4
    [{'A': 1}, {}, {'C': 3}, {'G': 4}, {'A': 5}]] #test 5

    pileup10_truncated_out = [[{'A': 1}], #test 1
    [{'A': 1}], #test 2
    [{'A': 1}], #test 3
    [{'A': 1}], #test 4
    [{'A': 1}]] #test 5

    pileup10_num_start_truncated = 0
    pileup10_num_end_truncated = 4

    """
    TESTS
    """
    @classmethod
    def setup_class(self):
        self.bam1 = TEST_PATH + '/data/quasi1.bam'
        self.bam2 = TEST_PATH + '/data/quasi2.bam'
        self.test_cp_files = (self.bam1, self.bam2)
        self.test_cp_ref = TEST_PATH+'/data/hxb2_pol.fas'
        self.references = parse_references_from_fasta(self.test_cp_ref)

    def test_construct_pileup_list(self):
        """
        test_construct_pileup_list - Checks that the pileup length and the
        first few indices of the pileup are correct.

        INPUT:
            [None]

        RETURN:
            [None]

        POST:
            [None]
        """

        bamPileup = parse_pileup_list_from_bam(self.references, self.test_cp_files)
        pileup_as_array = bamPileup.get_pileups_as_array()
        pileup_as_numerical_array = bamPileup.get_pileups_as_numerical_array()

        assert len(pileup_as_array)==2
        assert len(pileup_as_array[0])==2844
        assert len(pileup_as_array[1])==2844
        assert len(pileup_as_numerical_array[0])==(2844 * 4)
        assert len(pileup_as_numerical_array[1])==(2844 * 4)
        assert pileup_as_array[0][0:10] == [{'C': 12}, {'C': 12}, {'T': 12}, {'C': 12}, {'G': 2, 'C': 3, 'T': 1, 'A': 6}, {'G': 12}, {'G': 12}, {'T': 12}, {'C': 12}, {'G': 2, 'C': 3, 'T': 1, 'A': 6}]
        assert pileup_as_array[1][0:10] == [{'C': 12}, {'C': 12}, {'T': 12}, {'C': 12}, {'A': 6, 'C': 5, 'G': 1}, {'G': 12}, {'A': 4, 'G': 8}, {'T': 12}, {'C': 12}, {'A': 7, 'T': 1, 'C': 3, 'G': 1}]
    #end def

    @pytest.mark.parametrize("pileup,expected_truncated_pileup,expected_left_pos_truncated,expected_right_pos_truncated",
                             [(pileup4, pileup4_truncated_out, pileup4_num_start_truncated, pileup4_num_end_truncated),
                              (pileup5, pileup5_truncated_out, pileup5_num_start_truncated, pileup5_num_end_truncated),
                              (pileup6, pileup6_truncated_out, pileup6_num_start_truncated, pileup6_num_end_truncated),
                              (pileup7, pileup7_truncated_out, pileup7_num_start_truncated, pileup7_num_end_truncated),
                              (pileup8, pileup8_truncated_out, pileup8_num_start_truncated, pileup8_num_end_truncated),
                              (pileup9, pileup9_truncated_out, pileup9_num_start_truncated, pileup9_num_end_truncated),
                              (pileup10, pileup10_truncated_out, pileup10_num_start_truncated, pileup10_num_end_truncated)])
    def test_truncate_output(self, pileup, expected_truncated_pileup, expected_left_pos_truncated, expected_right_pos_truncated):
        """
        test_truncate_output - Checks that the expected truncated outputs
        matches the actual output.

        INPUT:
            [2D ARRAY OF DICTIONARIES] [pileup] # to be truncated
            [2D ARRAY OF DICTIONARIES] [expected_truncated_pileup]
            [2D ARRAY OF DICTIONARIES] [expected_left_pos_truncated]
            # number of contiguous left positions that were truncated
            [2D ARRAY OF DICTIONARIES] [expected_right_pos_truncated]
            # number of contiguous right positions that were truncated

        RETURN:
            [None]

        POST:
            Checks that the expected outputs match the actual output
        """
        pileups = Pileup_List([Pileup(bam) for bam in pileup])
        pileups.truncate_output()
        truncated = pileups.get_pileups_as_array()

        assert truncated == expected_truncated_pileup
        assert pileups.get_num_left_positions_truncated() == expected_left_pos_truncated
        assert pileups.get_num_right_positions_truncated() == expected_right_pos_truncated

    @pytest.mark.parametrize("pileup,expected_truncated_pileup,expected_left_pos_truncated,expected_right_pos_truncated",
                             [(pileup3, pileup3_remove_out, pileup3_num_start_truncated, pileup3_num_end_truncated)])
    def test_remove_no_coverage(self, pileup, expected_truncated_pileup, expected_left_pos_truncated, expected_right_pos_truncated):
        """
        test_remove_no_coverage - Checks that the after truncating all
        empty positions from the pileup that the output is as expected

        INPUT:
            [2D ARRAY OF DICTIONARIES] [pileup] # to be truncated
            [2D ARRAY OF DICTIONARIES] [expected_remove_no_coverage_pileup]
            [2D ARRAY OF DICTIONARIES] [expected_left_pos_truncated]
            # number of contiguous left positions that were truncated
            [2D ARRAY OF DICTIONARIES] [expected_right_pos_truncated]
            # number of contiguous right positions that were truncated

        RETURN:
            [None]

        POST:
            Checks that the expected outputs match the actual output
        """
        pileups = Pileup_List([Pileup(bam) for bam in pileup])
        pileups.remove_no_coverage()
        truncated = pileups.get_pileups_as_array()

        assert truncated == expected_truncated_pileup
        assert pileups.get_num_left_positions_truncated() == expected_left_pos_truncated
        assert pileups.get_num_right_positions_truncated() == expected_right_pos_truncated

    @pytest.mark.parametrize("pileup,startpos,endpos,expected_remove_no_coverage",
                             [(pileup1, 2, 4, pileup1_truncated_1),
                              (pileup1, 3, 3, pileup1_truncated_2)])
    def test_select_pileup_range_and_remove_no_coverage(self, pileup, startpos, endpos, expected_remove_no_coverage):
        """
        select_pileup_range_and_truncate_output - Checks if the program works
        as expected when removing all no coverage regions after first selecting
        a specified range.

        INPUT:
            [2D ARRAY of DICTIONARIES] [pileup]
            [INT] [startpos]
            [INT] [endpos]
            [2D ARRAY OF DICTIONARIES] [expected_remove_no_coverage]

        RETURN:
            TODO
        POST:
            TODO
        """

        pileups = Pileup_List([Pileup(bam) for bam in pileup])
        pileups.select_pileup_range(startpos, endpos)
        pileups.remove_no_coverage()
        truncated = pileups.get_pileups_as_array()

        assert truncated == expected_remove_no_coverage

    @pytest.mark.parametrize("pileup,startpos,endpos,pileup_select_range_expected_out,pileup_truncate_expected_out",
                             [(pileup2, 2, 4, pileup2_select_range_expected_out, pileup2_truncate_expected_out)])
    def select_pileup_range_and_truncate_output(self, pileup, startpos, endpos, pileup_select_range_expected_out, pileup_truncate_expected_out):
        """
        select_pileup_range_and_truncate_output - Checks if the program works
        as expected when truncating contiguous start and end regions after
        first selecting a specified range.

        INPUT:
            [2D ARRAY of DICTIONARIES] [pileup]
            [INT] [startpos]
            [INT] [endpos]
            [2D ARRAY OF DICTIONARIES] [pileup_select_range_expected_out]
            [2D ARRAY OF DICTIONARIES] [pileup_truncate_expected_out]
        RETURN:
            TODO
        POST:
            TODO
        """

        pileups = Pileup_List([Pileup(bam) for bam in pileup])
        pileups.select_pileup_range(startpos, endpos)
        select_pileup = pileups.get_pileups_as_array()

        # assert that the pileup positions before startpos and after endpos
        # have been ignored
        assert select_pileup == pileup_select_range_expected_out

        pileups.truncate_output()
        truncated_pileup = pileups.get_pileups_as_array()

        # assert that the pileup is truncated now as expected
        assert truncated_pileup == pileup_truncate_expected_out

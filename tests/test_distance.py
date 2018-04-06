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
from quasitools.distance import Pileup
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

    pileup3_num_start_truncated = 1
    pileup3_num_end_truncated = 0

    pileup4 = [[{'A': 1, 'T': 0, 'C': 0, 'G': 0},  # test 1 a
                {'A': 0, 'T': 2, 'C': 0, 'G': 0},  # test 1 b
                {'A': 0, 'T': 0, '-': 4, 'G': 0}],  # test 1 c
               [{'-': 1, 'T': 0, 'C': 0, 'G': 0},  # test 2 a
                {'A': 0, 'T': 2, 'C': 0, 'G': 0},  # test 2 b
                {'A': 0, '-': 5, 'C': 0, 'G': 0}]]  # test 2 c

    pileup4_trunc_ends_out = [[{'A': 0, 'T': 2, 'C': 0, 'G': 0}], #test 1
                            [{'A': 0, 'T': 2, 'C': 0, 'G': 0}]]

    pileup4_num_start_truncated = 1
    pileup4_num_end_truncated = 1

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

    pileup5_num_start_truncated = 1
    pileup5_num_end_truncated = 3

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

    pileup6_num_start_truncated = 2
    pileup6_num_end_truncated = 2

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

    pileup7_num_start_truncated = 0
    pileup7_num_end_truncated = 0

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

    pileup8_num_start_truncated = 0
    pileup8_num_end_truncated = 0

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

    pileup9_num_start_truncated = 1
    pileup9_num_end_truncated = 1

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

    pileup10_num_start_truncated = 0
    pileup10_num_end_truncated = 4

    pileup11 = [[{'A': 1}, {'T': 2}, {'C': 3}, {'G': 4}, {'A': 5}], #test 1
    [{'A': 1}, {'T': 2}, {'C': 3}, {'G': 4}, {'A': 5}], #test 2
    [{'A': 1}, {'T': 2}, {'C': 3}, {'G': 4}, {'A': 5}], #test 3
    [{'A': 1}, {'T': 2}, {'C': 3}, {'G': 4}, {'A': 5}], #test 4
    [{'A': 1}, {'T': 2}, {'C': 3}, {}, {}]] #test 5

    pileup11_startpos_1 = 2
    pileup11_endpos_1 = 4

    pileup11_startpos_2 = 3
    pileup11_endpos_2 = 3

    pileup11_truncated_1 = [[{'C': 3}], [{'C': 3}], [{'C': 3}], [{'C': 3}], [{'C': 3}]]

    pileup11_truncated_2 = [[], [], [], [], []]

    pileup12 = [[{}, {'T': 2}, {'C': 3}, {'G': 4}, {'A': 5}], #test 1
    [{'A': 1}, {'T': 2}, {'C': 3}, {'G': 4}, {}], #test 2
    [{}, {'T': 2}, {'C': 3}, {'G': 4}, {'A': 5}], #test 3
    [{}, {'T': 2}, {}, {'G': 4}, {}], #test 4
    [{'A': 1}, {'T': 2}, {'C': 3}, {'G': 4}, {'A': 5}]] #test 5

    pileup12_startpos = 2
    pileup12_endpos = 4

    pileup12_select_range_expected_out = [[{'C': 3}, {'G': 4}, {'A': 5}], #test 1
    [{'C': 3}, {'G': 4}, {}], #test 2
    [{'C': 3}, {'G': 4}, {'A': 5}], #test 3
    [{}, {'G': 4}, {}], #test 4
    [{'C': 3}, {'G': 4}, {'A': 5}]] #test 5

    pileup12_truncate_expected_out = [[{'G': 4}], [{'G': 4}], [{'G': 4}], [{'G': 4}], [{'G': 4}]]

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
            ---[BOOL] [normalize] # normalized or not
            ---[ARRAY] [pileup list]
            ---[ARRAY] [pileup_files] # file names corresponding to pileups
            ---[ARRAY] normalized or unnormalized similarity csv-format output
            ---[ARRAY] normalized or unnormalized distance csv-format output
            ---[INT or NONE] [startpos or default if NONE]
            ---[INT or NONE] [endpos or default if NONE]

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

    def test_construct_array_of_pileups(self):
        """
        test_construct_array_of_pileups - Checks that the pileup length and the
        first few indices of the pileup are correct.

        INPUT:
            [None]

        RETURN:
            [None]

        POST:
            [None]
        """

        test_cp_files = ((TEST_PATH+"/data/quasi1.bam"), (TEST_PATH+"/data/quasi2.bam"))
        test_cp_ref = TEST_PATH+"/data/hxb2_pol.fas"
        bamPileup = Pileup_List.construct_array_of_pileups(test_cp_files, test_cp_ref)
        pileup_as_array = bamPileup.get_pileups_as_array()

        assert len(pileup_as_array)==2
        assert len(pileup_as_array[0])==2844
        assert len(pileup_as_array[1])==2844
        assert pileup_as_array[0][0:10] == [{'C': 12}, {'C': 12}, {'T': 12}, {'C': 12}, {'G': 2, 'C': 3, 'T': 1, 'A': 6}, {'G': 12}, {'G': 12}, {'T': 12}, {'C': 12}, {'G': 2, 'C': 3, 'T': 1, 'A': 6}]
        assert pileup_as_array[1][0:10] == [{'C': 12}, {'C': 12}, {'T': 12}, {'C': 12}, {'A': 6, 'C': 5, 'G': 1}, {'G': 12}, {'A': 4, 'G': 8}, {'T': 12}, {'C': 12}, {'A': 7, 'T': 1, 'C': 3, 'G': 1}]
    #end def

    @pytest.mark.parametrize("pileup,expected_truncated_pileup,expected_left_pos_truncated,expected_right_pos_truncated",
                             [(pileup4, pileup4_trunc_ends_out, pileup4_num_start_truncated, pileup4_num_end_truncated),
                              (pileup5, pileup5_trunc_ends_out, pileup5_num_start_truncated, pileup5_num_end_truncated),
                              (pileup6, pileup6_trunc_ends_out, pileup6_num_start_truncated, pileup6_num_end_truncated),
                              (pileup7, pileup7_trunc_ends_out, pileup7_num_start_truncated, pileup7_num_end_truncated),
                              (pileup8, pileup8_trunc_ends_out, pileup8_num_start_truncated, pileup8_num_end_truncated),
                              (pileup9, pileup9_trunc_ends_out, pileup9_num_start_truncated, pileup9_num_end_truncated),
                              (pileup10, pileup10_trunc_ends_out, pileup10_num_start_truncated, pileup10_num_end_truncated)])
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
                             [(pileup3, pileup3_trunc_all_out, pileup3_num_start_truncated, pileup3_num_end_truncated)])
    def test_remove_no_coverage(self, pileup, expected_truncated_pileup, expected_left_pos_truncated, expected_right_pos_truncated):
        """
        test_remove_no_coverage - Checks that the after truncating all
        empty positions from the pileup that the output is as expected

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
        pileups.remove_no_coverage()
        truncated = pileups.get_pileups_as_array()

        assert truncated == expected_truncated_pileup
        assert pileups.get_num_left_positions_truncated() == expected_left_pos_truncated
        assert pileups.get_num_right_positions_truncated() == expected_right_pos_truncated

    @pytest.mark.parametrize("pileup,startpos,endpos,expected_remove_no_coverage",
                             [(pileup11, 2, 4, pileup11_truncated_1),
                              (pileup11, 3, 3, pileup11_truncated_2)])
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
                             [(pileup12, 2, 4, pileup12_select_range_expected_out, pileup12_truncate_expected_out)])
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

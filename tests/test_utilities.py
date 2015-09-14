"""
Copyright Government of Canada 2015

Written by: Eric Enns, National Microbiology Laboratory, Public Health Agency of Canada

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
from quasitools.utilities import *
from quasitools.reference import Reference
from pysam import AlignedSegment

def test_sam_alignment_to_padded_alignment():
    alignment = AlignedSegment()
    alignment.reference_start = 0
    alignment.query_sequence = 'AGCTTAGCTAGCTACCTATATCTTGGTCTTGGCCG'
    alignment.cigartuples = ((0,10), (2,1), (0,25))
    ref = Reference('test', 'AGCTTAGCTAAGCTACCTATATCTTGGTCTTGGCCG')

    (pad_ref, pad_match, pad_query) = sam_alignment_to_padded_alignment(alignment, ref)

    assert pad_ref == 'AGCTTAGCTAAGCTACCTATATCTTGGTCTTGGCCG'
    assert pad_match == '|||||||||| |||||||||||||||||||||||||'
    assert pad_query == 'AGCTTAGCTA-GCTACCTATATCTTGGTCTTGGCCG'

def test_pairwise_alignment_to_differences():
    ref = Reference('test', 'AGCTTAGCTAAGCTACCTATATCTTGGTCTTGGCCG')
    pad_ref = 'AGCTTAGCTAAGCTACCTATATCTTGGTCTTGGCCG'
    pad_query = 'AGCTTAGCTA-GCTACCTATATCTTGGTCTTGGCCG'
    ref_start = 0

    differences = pairwise_alignment_to_differences(pad_ref, pad_query, ref_start)

    assert differences == {10: '-'}

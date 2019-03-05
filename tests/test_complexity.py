"""
Copyright Government of Canada 2019

Written by: Ahmed Kidwai, Public Health Agency of Canada

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
import quasitools.haplotype as haplotype


class TestComplexity():

    # Create a list of haplotypes.
    haplotype_list_01= [haplotype.Haplotype("AAAAAAA"), 
            haplotype.Haplotype("AAAAAAA"),  
            haplotype.Haplotype("AAAAAAA")]
    haplotype_list_02 = [haplotype.Haplotype(""), haplotype.Haplotype("")]
    haplotype_list_03 = [haplotype.Haplotype("TAG")]

    # What i expect the output to be based on how we define a consensus sequence.
    expected_consensus_01 = "AAAAAAA" 
    expected_consensus_02 = ""
    expected_consensus_03 = "TAG"
    
    @classmethod
    def setup_class(self):
        self.expected_consensus = ""
        self.haplotype_list = ""
    
    @pytest.fixture(scope="function", params=[
        (haplotype_list_01, expected_consensus_01), 
        (haplotype_list_02, expected_consensus_02),
        (haplotype_list_03, expected_consensus_03)])
    def consensus(self,request):
        
        self.haplotype_list = request.param[0]
        self.expected_consensus = request.param[1]

        return self.expected_consensus
    
    def test_build_consensus(self,consensus):
    
        result_list = haplotype.build_consensus_from_haplotypes(self.haplotype_list)
        result = ''.join(result_list)
        assert result == consensus




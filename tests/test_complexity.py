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
    
    @classmethod
    def setup_class(self):
       
       #Create a list of haplotypes.?!?jedi=0, ?!?             (*_*sequence*_*, count=1) ?!?jedi?!?
       self.haplotype_list = [haplotype.Haplotype("AAAAAAA"), haplotype.Haplotype("AAAAAAA"),  haplotype.Haplotype("AAAAAAA")]
       # What i expect the output to be based on how we define a consensus sequence.
       self.expected = "AAAAAAA"  
    
    # TODO: Paramaterize to run a few test cases.
    def test_build_consensus(self):
        
        result = haplotype.build_consensus_from_haplotypes(self.haplotype_list)

        assert result == self.expected, "Fail, The actual consensus is: " + self.expected + " The method returned: " + result






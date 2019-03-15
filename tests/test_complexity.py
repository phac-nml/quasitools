"""
# =============================================================================

Copyright Government of Canada 2019

Written by: Ahmed Kidwai, Public Health Agency of Canada,
    National Microbiology Laboratory

Funded by the National Micriobiology Laboratory and the Genome Canada / Alberta
    Innovates Bio Solutions project "Listeria Detection and Surveillance
    using Next Generation Genomics"

Licensed under the Apache License, Version 2.0 (the "License"); you may not use
this file except in compliance with the License. You may obtain a copy of the
License at:

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software distributed
under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR
CONDITIONS OF ANY KIND, either express or implied. See the License for the
specific language governing permissions and limitations under the License.

# =============================================================================
"""
import os
import pytest
import quasitools.commands.cmd_complexity as Complexity
from click.testing import CliRunner


class Test_BAM_Complexity:
    @classmethod
    def setup(self):
        self.bam_location = 'tests/data/complexity.bam'
        self.reference_location =  'tests/data/complexity_reference.fasta'
   
    def test_complexity_bam(self):

        runner = CliRunner()
        result = runner.invoke(Complexity.bam, [self.reference_location, self.bam_location, "50", "what_in_the_fuck.csv"])
        assert result.exit_code == 2
        #assert result.output == "test"
        

            

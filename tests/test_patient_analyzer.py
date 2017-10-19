"""
Copyright Government of Canada 2017

Written by: Camy Tran, National Microbiology Laboratory,
            Public Health Agency of Canada

Licensed under the Apache License, Version 2.0 (the "License"); you may not use
this work except in compliance with the License. You may obtain a copy of the
License at:

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software distributed
under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR
CONDITIONS OF ANY KIND, either express or implied. See the License for the
specific language governing permissions and limitations under the License.
"""

import pdb
import os
import pytest
from quasitools.patient_analyzer import PatientAnalyzer
from collections import defaultdict


TEST_PATH = os.path.dirname(os.path.abspath(__file__))

class TestPatientAnalyzer:
    def setup_class(self):
        options = {}
        options["quiet"] = False

        self.patient_analyzer = PatientAnalyzer(id="test", reads=TEST_PATH+"/data/reads_w_K103N.fastq",
                                 reference=TEST_PATH+"/data/hxb2_pol.fas", consensus_pct=50, 
                                 output_dir=TEST_PATH+"/testpa", genes_file=TEST_PATH+"/data/hxb2_pol.bed",
                                 mutation_db=TEST_PATH+"/data/mutation_db.tsv", options=options)

        assert(self.patient_analyzer.reference == TEST_PATH+"/data/hxb2_pol.fas")
    
    def test_filter_reads(self):
        filters = defaultdict(dict)
        filters["reads"]["length_cutoff"] = 100
        filters["reads"]["score_cutoff"] = 30
        filters["reads"]["ns"] = 1

        self.patient_analyzer.filter_reads(filters["reads"])
    
    def test_downsample_reads(self):
        self.patient_analyzer.downsample_reads(50)
        self.patient_analyzer.downsample_reads(1000)

    def test_analyze_reads(self):
        filters = defaultdict(dict)
        filters["variant_filtering"] = defaultdict(dict)
        filters["variant_filtering"]["error_rate"] = 0.0021
        filters["variant_filtering"]["min_qual"] = 30
        filters["variant_filtering"]["min_dp"] = 100
        filters["variant_filtering"]["min_ac"] = 5

        filters["mutations"]["min_freq"] = 0.01

        generate_consensus = True

        self.patient_analyzer.analyze_reads(filters, generate_consensus)

        

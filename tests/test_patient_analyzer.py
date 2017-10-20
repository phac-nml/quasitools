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
import shutil
import pytest
from quasitools.patient_analyzer import PatientAnalyzer
from collections import defaultdict


TEST_PATH = os.path.dirname(os.path.abspath(__file__))
READS = TEST_PATH + "/data/reads_w_K103N.fastq"
REFERENCE = TEST_PATH + "/data/hxb2_pol.fas"
GENES_FILE = TEST_PATH + "/data/hxb2_pol.bed"
MUTATION_DB = TEST_PATH + "/data/mutation_db.tsv"
OUTPUT_DIR = TEST_PATH + "/test_patient_analyzer_output"


class TestPatientAnalyzer:
    @classmethod
    def setup_class(self):
        # Test defaults
        self.patient_analyzer = PatientAnalyzer(id="test",reads=READS,
                                 reference=REFERENCE,
                                 output_dir=OUTPUT_DIR, 
                                 genes_file=GENES_FILE,
                                 mutation_db=MUTATION_DB,
                                 quiet=False, consensus_pct=20)

    def test_filter_reads(self):
        # test default values
        filters = defaultdict(dict)
        filters["length_cutoff"] = 100
        filters["score_cutoff"] = 30
        filters["ns"] = 1

        self.patient_analyzer.filter_reads(filters)

    def test_downsample_reads(self):
        # test default
        self.patient_analyzer.downsample_reads(10000)

        # test other
        self.patient_analyzer.downsample_reads(1000)

        # test disable
        self.patient_analyzer.downsample_reads(-1)

    def test_analyze_reads(self):
        # test defaults
        filters = defaultdict(dict)
        filters["error_rate"] = 0.0021
        filters["min_qual"] = 30
        filters["min_dp"] = 100
        filters["min_ac"] = 5
        filters["min_freq"] = 0.01

        generate_consensus = True
        reporting_threshold = 20

        self.patient_analyzer.analyze_reads(filters, reporting_threshold, generate_consensus)

        # Remove the output directory so that multiple tests (python 2.x, 3.x, etc.)
        # can run without erroring out with "File/directory exists"
        shutil.rmtree("%s" % OUTPUT_DIR)


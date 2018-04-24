"""
Copyright Government of Canada 2018

Written by: Matthew Fogel, National Microbiology Laboratory,
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

import os
import shutil
from collections import defaultdict
import pytest
from quasitools.patient_analyzer import PatientAnalyzer
import Bio.SeqIO

TEST_PATH = os.path.dirname(os.path.abspath(__file__))
TEST_PATH = os.path.dirname(os.path.abspath(__file__))
READS = TEST_PATH + "/data/reads_w_K103N.fastq"
FORWARD = TEST_PATH + "/data/forward.fastq"
REVERSE = TEST_PATH + "data/reverse.fastq"
REFERENCE = TEST_PATH + "/data/hxb2_pol.fas"
GENES_FILE = TEST_PATH + "/data/hxb2_pol.bed"
MUTATION_DB = TEST_PATH + "/data/mutation_db.tsv"
OUTPUT_DIR = TEST_PATH + "/test_patient_analyzer_output"


class TestQualityControl:
    @classmethod
    def setup_class(self):
        # Test defaults
        self.quality_control = QualityControl()

    def test_filter_reads(self):
        read_filters = defaultdict(dict)

        length_cutoff = 100
        score_cutoff = 30

        # read_filters["TRIMMING"] = "trimming"
        # read_filters["MASKING"] = "masking"

        read_filters["length_cutoff"] = length_cutoff
        read_filters["mean_cutoff"] = score_cutoff
        read_filters["ns"] = True
        read_filters["minimum_quality"] = min_qual

        assert self.quality_control.get_amount_filtered()["status"] == 0
        self.quality_control.filter_reads(filters)
        self.filtered_reads = "%s/filtered.fastq" % output_dir
        seq_rec_obj = Bio.SeqIO.parse(self.filtered_reads, "fastq")

        for seq in seq_rec_obj:
            avg_score = filters["score_cutoff"] + 1
            avg_score = (float(sum(seq.letter_annotations['phred_quality'])) /
                         float(len(seq.letter_annotations['phred_quality'])))

            # check that length and score are both over threshold
            assert self.quality_control.get_amount_filtered()["status"] == 1

            assert len(seq.seq) >= read_filters["length_cutoff"] and \
                avg_score >= read_filters["score_cutoff"]

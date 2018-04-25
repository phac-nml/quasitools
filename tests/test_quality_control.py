"""
Copyright Government of Canada 2018

Written by: Matthew Fogel and Camy Tran, National Microbiology Laboratory,
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
from quasitools.quality_control import QualityControl
import Bio.SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

TEST_PATH = os.path.dirname(os.path.abspath(__file__))
READS = TEST_PATH + "/data/reads_w_K103N.fastq"
OUTPUT_DIR = TEST_PATH + "/test_quality_control_output"
FILTERED_DIR = OUTPUT_DIR + "/filtered.fastq"

class TestQualityControl:
    @classmethod
    def setup_class(self):
        # Test defaults
        self.quality_control = QualityControl()
        if not os.path.isdir(OUTPUT_DIR):
            os.mkdir(OUTPUT_DIR)

    def test_passes_filters(self):
        failed_status = {0: "success", 1: "length", 2: "score", 3: "ns"}

        quality_filters = defaultdict(dict)

        length_cutoff = 2
        score_cutoff = 30
        min_qual = 30

        quality_filters["length_cutoff"] = length_cutoff
        quality_filters["mean_cutoff"] = score_cutoff
        quality_filters["ns"] = True
        quality_filters["minimum_quality"] = min_qual

        # sample Biopython read for testing
        seq = Seq("GATC")
        seq_record = SeqRecord(seq)
        seq_record.letter_annotations["phred_quality"] = [35, 3, 5, 7]
        key = self.quality_control.passes_filters(seq_record,
                                                  quality_filters)
        assert failed_status.get(key) == "score" # did not pass filters due to
                                                 # score

        seq = Seq("G")
        seq_record = SeqRecord(seq)
        seq_record.letter_annotations["phred_quality"] = [35]
        key = self.quality_control.passes_filters(seq_record,
                                                  quality_filters)
        assert failed_status.get(key) == "length" # did not pass filters due to
                                                  # score

        # test where multiple characters are masked due to low quality
        seq = Seq("GNNN")
        seq_record = SeqRecord(seq)
        seq_record.letter_annotations["phred_quality"] = [40, 40, 40, 40]
        key = self.quality_control.passes_filters(seq_record,
                                                  quality_filters)
        assert failed_status.get(key) == "ns" # did not pass filters due to
                                              # score

    def test_filter_reads(self):

        # test without iterative trimming or masking enabled

        read_filters = defaultdict(dict)

        length_cutoff = 100
        score_cutoff = 30
        min_qual = 30

        # read_filters["TRIMMING"] = "trimming"
        # read_filters["MASKING"] = "masking"

        read_filters["length_cutoff"] = length_cutoff
        read_filters["score_cutoff"] = score_cutoff
        read_filters["ns"] = True
        read_filters["minimum_quality"] = min_qual

        assert self.quality_control.get_amount_filtered()["status"] == 0
        self.quality_control.filter_reads(READS, FILTERED_DIR, read_filters)
        seq_rec_obj = Bio.SeqIO.parse(FILTERED_DIR, "fastq")

        for seq in seq_rec_obj:
            avg_score = read_filters["score_cutoff"] + 1
            avg_score = (float(sum(seq.letter_annotations['phred_quality'])) /
                         float(len(seq.letter_annotations['phred_quality'])))

            # check that length and score are both over threshold
            assert self.quality_control.get_amount_filtered()["status"] == 1

            assert len(seq.seq) >= read_filters["length_cutoff"] and \
                avg_score >= read_filters["score_cutoff"]

        # test with iterative trimming and masking enabled

        quality_filters = defaultdict(dict)

        quality_filters["TRIMMING"] = "trimming"
        quality_filters["MASKING"] = "masking"

        length_cutoff = 2
        score_cutoff = 30
        min_qual = 30

        quality_filters["length_cutoff"] = length_cutoff
        quality_filters["mean_cutoff"] = score_cutoff
        quality_filters["ns"] = True
        quality_filters["minimum_quality"] = min_qual

        # sample Biopython reads for testing
        seq_record_list = []
        seq = Seq("GATC")
        seq_record = SeqRecord(seq)
        seq_record.letter_annotations["phred_quality"] = [35, 3, 5, 7]
        seq_record_list.append(seq_record)

        seq = Seq("CGTA")
        seq_record = SeqRecord(seq)
        seq_record.letter_annotations["phred_quality"] = [4, 3, 5, 7]
        seq_record_list.append(seq_record)

        seq = Seq("CGTA")
        seq_record = SeqRecord(seq)
        seq_record.letter_annotations["phred_quality"] = [29, 29, 29, 29]
        seq_record_list.append(seq_record)

        INPUT_DIR = TEST_PATH+"/data/sample.fastq"

        Bio.SeqIO.write(seq_record_list, INPUT_DIR, "fastq")

        self.quality_control.filter_reads(INPUT_DIR,
                                          FILTERED_DIR,
                                          quality_filters)

        seq_rec_obj = Bio.SeqIO.parse(FILTERED_DIR, "fastq")
        assert(sum(1 for seq in seq_rec_obj) == 0)

        # more sample Biopython reads for testing
        seq_record_list = []
        seq = Seq("GATC")
        seq_record = SeqRecord(seq)
        seq_record.id = "first"
        seq_record.letter_annotations["phred_quality"] = [35, 30, 50, 70]
        seq_record_list.append(seq_record)

        seq = Seq("CGTA")
        seq_record = SeqRecord(seq)
        seq_record.id = "second"
        seq_record.letter_annotations["phred_quality"] = [3, 3, 5, 7]
        seq_record_list.append(seq_record)

        seq = Seq("CGTA")
        seq_record = SeqRecord(seq)
        seq_record.id = "third"
        seq_record.letter_annotations["phred_quality"] = [30, 33, 35, 30]
        seq_record_list.append(seq_record)

        INPUT_DIR = TEST_PATH+"/data/sample2.fastq"

        Bio.SeqIO.write(seq_record_list, INPUT_DIR, "fastq")

        self.quality_control.filter_reads(INPUT_DIR,
                                          FILTERED_DIR,
                                          quality_filters)

        seq_rec_obj = Bio.SeqIO.parse(FILTERED_DIR, "fastq")
        assert(sum(1 for seq in seq_rec_obj) == 2)

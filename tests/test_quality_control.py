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

    def test_get_median_score(self):
        seq = Seq("GATC")
        seq_record = SeqRecord(seq)
        seq_record.id = "first"
        seq_record.letter_annotations["phred_quality"] = [35, 30, 50, 70]

        assert self.quality_control.get_median_score(seq_record) == 40

        seq = Seq("CAT")
        seq_record = SeqRecord(seq)
        seq_record.id = "first"
        seq_record.letter_annotations["phred_quality"] = [35, 30, 50]

        assert self.quality_control.get_median_score(seq_record) == 30

    def test_get_mean_score(self):
        seq = Seq("GATC")
        seq_record = SeqRecord(seq)
        seq_record.id = "first"
        seq_record.letter_annotations["phred_quality"] = [20, 40, 60, 80]

        assert self.quality_control.get_mean_score(seq_record) == 50

    def test_trim_reads(self):
        seq = Seq("GATC")
        seq_record = SeqRecord(seq)
        seq_record.id = "first"
        seq_record.letter_annotations["phred_quality"] = [30, 30, 3, 1]

        advanced_filters = defaultdict(dict)

        advanced_filters["trimming"] = True

        advanced_filters["length_cutoff"] = 2
        advanced_filters["mean_cutoff"] = 30
        advanced_filters["ns"] = True
        advanced_filters["minimum_quality"] = 30

        trimmed_read = self.quality_control.trim_read(seq_record,
                                                      advanced_filters)

        assert len(trimmed_read.seq) == 2

    def test_mask_reads(self):
        seq = Seq("GATC")
        seq_record = SeqRecord(seq)
        seq_record.id = "first"
        seq_record.letter_annotations["phred_quality"] = [30, 30, 3, 1]

        advanced_filters = defaultdict(dict)

        advanced_filters["masking"] = True

        advanced_filters["length_cutoff"] = 2
        advanced_filters["mean_cutoff"] = 30
        advanced_filters["ns"] = True
        advanced_filters["minimum_quality"] = 30

        self.quality_control.mask_read(seq_record, advanced_filters)

        assert len(seq_record.seq) == 4
        assert seq_record.seq[3] == 'N'
        assert seq_record.seq[2] == 'N'
        assert seq_record.seq[1] == 'A'
        assert seq_record.seq[0] == 'G'

    def test_passes_filters(self):
        failed_status = {0: "success", 1: "length", 2: "score", 3: "ns"}

        advanced_filters = defaultdict(dict)

        length_cutoff = 2
        score_cutoff = 30
        min_qual = 30

        advanced_filters["length_cutoff"] = length_cutoff
        advanced_filters["mean_cutoff"] = score_cutoff
        advanced_filters["ns"] = True
        advanced_filters["minimum_quality"] = min_qual

        # sample Biopython read for testing
        seq = Seq("GATC")
        seq_record = SeqRecord(seq)
        seq_record.letter_annotations["phred_quality"] = [35, 3, 5, 7]
        key = self.quality_control.passes_filters(seq_record,
                                                  advanced_filters)
        assert failed_status.get(key) == "score" # did not pass filters due to
                                                 # score

        seq = Seq("G")
        seq_record = SeqRecord(seq)
        seq_record.letter_annotations["phred_quality"] = [35]
        key = self.quality_control.passes_filters(seq_record,
                                                  advanced_filters)
        assert failed_status.get(key) == "length" # did not pass filters due to
                                                  # score

        # test where multiple characters are masked due to low quality
        seq = Seq("GNNN")
        seq_record = SeqRecord(seq)
        seq_record.letter_annotations["phred_quality"] = [40, 40, 40, 40]
        key = self.quality_control.passes_filters(seq_record,
                                                  advanced_filters)
        assert failed_status.get(key) == "ns" # did not pass filters due to
                                              # score

    def test_filter_reads(self):

        # test with iterative trimming and masking enabled

        advanced_filters = defaultdict(dict)

        advanced_filters["trimming"] = True
        advanced_filters["masking"] = True

        length_cutoff = 2
        score_cutoff = 30
        min_qual = 30

        advanced_filters["length_cutoff"] = length_cutoff
        advanced_filters["mean_cutoff"] = score_cutoff
        advanced_filters["ns"] = True
        advanced_filters["minimum_quality"] = min_qual

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
                                          advanced_filters)

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
                                          advanced_filters)

        seq_rec_obj = Bio.SeqIO.parse(FILTERED_DIR, "fastq")
        assert(sum(1 for seq in seq_rec_obj) == 2)

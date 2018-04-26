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

TRIMMING = "trimming"
MASKING = "masking"
MASK_CHARACTER = "N"
MINIMUM_QUALITY = "minimum_quality"
LENGTH_CUTOFF = "length_cutoff"
MEDIAN_CUTOFF = "median_cutoff"
MEAN_CUTOFF = "mean_cutoff"

class TestQualityControl:
    @classmethod
    def setup_class(self):
        # Test defaults
        self.quality_control = QualityControl()
        if not os.path.isdir(OUTPUT_DIR):
            os.mkdir(OUTPUT_DIR)

    @pytest.fixture
    def filters(self):
        filtering = defaultdict(dict)

        filtering["trimming"] = True
        filtering["masking"] = True

        filtering["length_cutoff"] = 2
        filtering["mean_cutoff"] = 30
        filtering["ns"] = True
        filtering["minimum_quality"] = 30

        return filtering

    def test_get_median_score(self):
        """
        test_get_median_score - Checks that the function correctly calculates
        median score

        INPUT:
            [None]

        RETURN:
            [None]

        POST:
            [None]
        """
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
        """
        test_get_mean_score - Checks that the function correctly calculates
        mean score

        INPUT:
            [None]

        RETURN:
            [None]

        POST:
            [None]
        """
        seq = Seq("GATC")
        seq_record = SeqRecord(seq)
        seq_record.id = "first"
        seq_record.letter_annotations["phred_quality"] = [20, 40, 60, 80]

        assert self.quality_control.get_mean_score(seq_record) == 50

    def test_trim_reads(self, filters):
        """
        test_trim_reads - Checks that the reads have been trimmed until all
        filters pass or read falls below length cutoff

        INPUT:
            [dict] [filters] # provided with fixture

        RETURN:
            [None]

        POST:
            [None]
        """
        seq = Seq("GATC")
        seq_record = SeqRecord(seq)
        seq_record.id = "first"
        seq_record.letter_annotations["phred_quality"] = [30, 30, 3, 1]

        filters[MASKING] == False
        filters[TRIMMING] = True
        trimmed_read = self.quality_control.trim_read(seq_record, filters)

        assert len(trimmed_read.seq) == 2

    def test_mask_reads(self, filters):
        """
        test_mask_reads - Checks that all low coverage regions have been
        masked with an N.

        INPUT:
            [dict] [filters] # provided with fixture

        RETURN:
            [None]

        POST:
            [None]
        """
        seq = Seq("GATC")
        seq_record = SeqRecord(seq)
        seq_record.id = "first"
        seq_record.letter_annotations["phred_quality"] = [30, 30, 3, 1]

        filters[TRIMMING] = False
        filters[MASKING] == True
        self.quality_control.mask_read(seq_record, filters)

        assert len(seq_record.seq) == 4
        assert seq_record.seq[0:4] == 'GANN'

    def test_passes_filters(self, filters):
        """
        test_passes_filters - Checks that the correct status is returned if
        the read fails to pass at least one of the filters

        INPUT:
            [dict] [filters] # provided with fixture

        RETURN:
            [None]

        POST:
            [None]
        """
        failed_status = {0: "success", 1: "length", 2: "score", 3: "ns"}

        # sample Biopython read for testing
        seq = Seq("GATC")
        seq_record = SeqRecord(seq)
        seq_record.letter_annotations["phred_quality"] = [35, 3, 5, 7]
        key = self.quality_control.passes_filters(seq_record, filters)
        assert failed_status.get(key) == "score" # did not pass filters due to
                                                 # score

        seq = Seq("G")
        seq_record = SeqRecord(seq)
        seq_record.letter_annotations["phred_quality"] = [35]
        key = self.quality_control.passes_filters(seq_record, filters)
        assert failed_status.get(key) == "length" # did not pass filters due to
                                                  # score

        # test where multiple characters are masked due to low quality
        seq = Seq("GNNN")
        seq_record = SeqRecord(seq)
        seq_record.letter_annotations["phred_quality"] = [40, 40, 40, 40]
        key = self.quality_control.passes_filters(seq_record, filters)
        assert failed_status.get(key) == "ns" # did not pass filters due to
                                              # score

    def test_filter_reads_with_mean_score(self, filters):
        """
        test_filter_reads_with_mean_score - Checks that the filter_reads
        function performs as expected with iterative trimming and/or masking
        enabled. This includes testing with the mean used for score.

        INPUT:
            [dict] [filters] # provided with fixture

        RETURN:
            [None]

        POST:
            [None]
        """
        # test with iterative trimming and masking enabled

        # sample Biopython reads for testing
        seq_record_list = []
        seq = Seq("GATC")
        seq_record = SeqRecord(seq)
        seq_record.letter_annotations["phred_quality"] = [35, 3, 5, 7]
        seq_record_list.append(seq_record)

        seq_record = SeqRecord(seq)
        seq_record.letter_annotations["phred_quality"] = [4, 3, 5, 7]
        seq_record_list.append(seq_record)

        seq_record = SeqRecord(seq)
        seq_record.letter_annotations["phred_quality"] = [29, 29, 29, 29]
        seq_record_list.append(seq_record)

        seq_record = SeqRecord(seq)
        seq_record.letter_annotations["phred_quality"] = [35, 30, 50, 70]
        seq_record_list.append(seq_record)

        seq_record = SeqRecord(seq)
        seq_record.letter_annotations["phred_quality"] = [30, 33, 35, 30]
        seq_record_list.append(seq_record)

        seq_record = SeqRecord(seq)
        seq_record.letter_annotations["phred_quality"] = [30, 30, 30, 30]
        seq_record_list.append(seq_record)

        seq_record = SeqRecord(seq)
        seq_record.letter_annotations["phred_quality"] = [29, 29, 31, 31]
        seq_record_list.append(seq_record)

        INPUT_DIR = TEST_PATH+"/data/sample.fastq"

        Bio.SeqIO.write(seq_record_list, INPUT_DIR, "fastq")

        self.quality_control.filter_reads(INPUT_DIR, FILTERED_DIR, filters)

        seq_rec_obj = Bio.SeqIO.parse(FILTERED_DIR, "fastq")
        assert(sum(1 for seq in seq_rec_obj) == 4)

        self.quality_control.filter_reads(INPUT_DIR, FILTERED_DIR, filters)

        seq_rec_obj = Bio.SeqIO.parse(FILTERED_DIR, "fastq")

        for read in seq_rec_obj:
            mean = (float(sum(read.letter_annotations['phred_quality'])) /
                    float(len(read.letter_annotations['phred_quality'])))

            assert mean >= filters["mean_cutoff"]

            for base_pos in range(0, len(read.seq)):
                if (read.letter_annotations['phred_quality'][base_pos]
                    < filters['minimum_quality']):
                    assert read.seq[base_pos] == MASK_CHARACTER

    def test_filter_reads_with_median_score(self, filters):
        """
        test_filter_reads_with_median_score - Checks that the filter_reads
        function performs as expected with iterative trimming and/or masking
        enabled. This includes testing with the median used for score.

        INPUT:
            [dict] [filters] # provided with fixture

        RETURN:
            [None]

        POST:
            [None]
        """
        # test with iterative trimming and masking enabled
        # tests with median_cutoff used instead of mean_cutoff
        filters["median_cutoff"] = 30
        del filters["mean_cutoff"]

        INPUT_DIR = TEST_PATH+"/data/sample2.fastq"

        # more sample Biopython reads for testing
        seq_record_list = []
        seq = Seq("GATC")

        seq_record = SeqRecord(seq)
        seq_record.letter_annotations["phred_quality"] = [35, 30, 50, 10]
        seq_record_list.append(seq_record)

        seq_record = SeqRecord(seq)
        seq_record.letter_annotations["phred_quality"] = [29, 30, 29, 29]
        seq_record_list.append(seq_record)

        seq_record = SeqRecord(seq)
        seq_record.letter_annotations["phred_quality"] = [29, 29, 30, 29]
        seq_record_list.append(seq_record)

        seq_record = SeqRecord(seq)
        seq_record.letter_annotations["phred_quality"] = [29, 29, 29, 29]
        seq_record_list.append(seq_record)


        Bio.SeqIO.write(seq_record_list, INPUT_DIR, "fastq")
        self.quality_control.filter_reads(INPUT_DIR, FILTERED_DIR, filters)
        seq_rec_obj = Bio.SeqIO.parse(FILTERED_DIR, "fastq")

        for read in seq_rec_obj:
            length = len(read.seq)

            scores = list(read.letter_annotations['phred_quality'])
            if length % 2 == 0:
                median_score = ((scores[int((length - 1) // 2)] +
                                 scores[int((length - 1) // 2) + 1]) / 2)
            else:
                median_score = scores[int((length - 1) // 2)]

            assert self.quality_control.get_median_score(read) >= \
                   filters.get("median_cutoff")

            for base_pos in range(0, len(read.seq)):
                if (read.letter_annotations['phred_quality'][base_pos]
                    < filters['minimum_quality']):
                    assert read.seq[base_pos] == MASK_CHARACTER

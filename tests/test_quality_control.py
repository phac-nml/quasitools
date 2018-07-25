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

# globals

from quasitools.quality_control import TRIMMING

# used for masking
from quasitools.quality_control import MASKING
from quasitools.quality_control import MASK_CHARACTER
from quasitools.quality_control import MIN_READ_QUAL

# used in quality_control.passes_filters
from quasitools.quality_control import LENGTH_CUTOFF
from quasitools.quality_control import MEDIAN_CUTOFF
from quasitools.quality_control import MEAN_CUTOFF
from quasitools.quality_control import NS

# used in qaulity_control.passes_filters
from quasitools.quality_control import PASS
from quasitools.quality_control import FAIL_LENGTH
from quasitools.quality_control import FAIL_SCORE
from quasitools.quality_control import FAIL_NS

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

        filtering[TRIMMING] = True
        filtering[MASKING] = True

        filtering[LENGTH_CUTOFF] = 2
        filtering[MEAN_CUTOFF] = 30
        filtering[MIN_READ_QUAL] = 30
        #filtering[NS] = True

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

    def test_filter_and_trim_reads(self, filters):
        """
        test_filter_and_trim_reads - Checks that the reads have been trimmed
        until all filters pass or read falls below length cutoff and that this
        result is reflected in the filtered read saved to file by filter reads

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

        INPUT_DIR = TEST_PATH + "/data/sample_filter.fastq"

        Bio.SeqIO.write(seq_record, INPUT_DIR, "fastq")

        self.quality_control.filter_reads(INPUT_DIR, FILTERED_DIR, filters)

        # there should be only one object in the seq_rec_obj iterable
        seq_rec_obj = Bio.SeqIO.parse(FILTERED_DIR, "fastq")
        recs = list(seq_rec_obj)
        assert len(recs) == 1
        assert len(recs[0].seq) == 2

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

    def test_filter_and_mask_reads(self, filters):
        """
        test_filter_and_mask_reads - Checks that all low coverage regions have
        been masked with an N and that this result is
        reflected in the filtered read saved to file by filter reads

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
        filters[NS] = False
        # we are not testing any other functions that will be called when
        # score < MEAN_CUTOFF, only testing masking, so we set below filter val
        # to 0
        filters[MEAN_CUTOFF] = 0

        INPUT_DIR = TEST_PATH + "/data/sample_mask.fastq"

        Bio.SeqIO.write(seq_record, INPUT_DIR, "fastq")

        self.quality_control.filter_reads(INPUT_DIR, FILTERED_DIR, filters)

        seq_rec_obj = Bio.SeqIO.parse(FILTERED_DIR, "fastq")
        recs = list(seq_rec_obj)
        assert len(recs) == 1
        # there should be only one object in the seq_rec_obj iterable
        assert len(recs[0].seq) == 4
        assert recs[0].seq[0:4] == 'GANN'

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
        filters[NS] = False
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
        assert key == FAIL_SCORE # did not pass filters due
                                                 # to score

        seq = Seq("G")
        seq_record = SeqRecord(seq)
        seq_record.letter_annotations["phred_quality"] = [35]
        key = self.quality_control.passes_filters(seq_record, filters)
        assert key == FAIL_LENGTH # did not pass filters due
                                                  # to score

        # test where multiple characters are already masked due to low quality
        filters[NS] = True  # turn on filtering of ns in the sequence
        filters[MASKING] = False  # turn off masking of low coverage regions
        seq = Seq("GNNN")
        seq_record = SeqRecord(seq)
        seq_record.letter_annotations["phred_quality"] = [80, 20, 20, 20]
        key = self.quality_control.passes_filters(seq_record, filters)
        assert key == FAIL_NS # did not pass filters due to
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

        INPUT_DIR = TEST_PATH+"/data/sample_mean.fastq"

        Bio.SeqIO.write(seq_record_list, INPUT_DIR, "fastq")

        self.quality_control.filter_reads(INPUT_DIR, FILTERED_DIR, filters)

        seq_rec_obj = Bio.SeqIO.parse(FILTERED_DIR, "fastq")
        assert(sum(1 for seq in seq_rec_obj) == 4)

        self.quality_control.filter_reads(INPUT_DIR, FILTERED_DIR, filters)

        seq_rec_obj = Bio.SeqIO.parse(FILTERED_DIR, "fastq")

        for read in seq_rec_obj:
            mean = (float(sum(read.letter_annotations['phred_quality'])) /
                    float(len(read.letter_annotations['phred_quality'])))

            assert mean >= filters[MEAN_CUTOFF]

            for base_pos in range(0, len(read.seq)):
                if (read.letter_annotations['phred_quality'][base_pos]
                    < filters[MIN_READ_QUAL]):
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
        filters[MEDIAN_CUTOFF] = 30
        del filters[MEAN_CUTOFF]

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
                   filters.get(MEDIAN_CUTOFF)

            for base_pos in range(0, len(read.seq)):
                if (read.letter_annotations['phred_quality'][base_pos]
                    < filters[MIN_READ_QUAL]):
                    assert read.seq[base_pos] == MASK_CHARACTER

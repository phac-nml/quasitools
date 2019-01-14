"""
Copyright Government of Canada 2017 - 2018

Written by: Camy Tran and Matthew Fogel, National Microbiology Laboratory,
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

# globals
TEST_PATH = os.path.dirname(os.path.abspath(__file__))
READS = TEST_PATH + "/data/reads_w_K103N.fastq"
FORWARD = TEST_PATH + "/data/forward.fastq"
REVERSE = TEST_PATH + "/data/reverse.fastq"
REFERENCE = TEST_PATH + "/data/hxb2_pol.fas"
BED4_FILE = TEST_PATH + "/data/hxb2_pol.bed"
MUTATION_DB = TEST_PATH + "/data/mutation_db.tsv"
OUTPUT_DIR = TEST_PATH + "/test_patient_analyzer_output"
FILTERED_DIR = OUTPUT_DIR + "/filtered.fastq"

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

# used in patient_analyzer
from quasitools.patient_analyzer import ERROR_RATE
from quasitools.patient_analyzer import MIN_VARIANT_QUAL
from quasitools.patient_analyzer import MIN_AC
from quasitools.patient_analyzer import MIN_DP
from quasitools.patient_analyzer import MIN_FREQ

class TestPatientAnalyzer:
    @classmethod
    def setup_class(self):
        # Test defaults
        self.patient_analyzer = PatientAnalyzer(id="test",reads=READS,
                                 reference=REFERENCE,
                                 output_dir=OUTPUT_DIR,
                                 BED4_file=BED4_FILE,
                                 mutation_db=MUTATION_DB,
                                 quiet=False, consensus_pct=20)


    def test_combine_reads(self):
        # Combine the fwd and reverse into one file
        reads = "%s/combined_reads.fastq" % OUTPUT_DIR
        cat_cmd = "cat %s %s > %s" % (FORWARD, REVERSE, reads)
        os.system(cat_cmd)

        assert os.path.isfile("%s/combined_reads.fastq" % OUTPUT_DIR)

    def test_filter_reads(self):

        # tests for filtering of reads without iterative trimming or masking
        # of coverage regions enabled

        quality_filters = defaultdict(dict)

        quality_filters[LENGTH_CUTOFF] = 100
        quality_filters[MEAN_CUTOFF] = 30
        quality_filters[MASKING] = True
        quality_filters[MIN_READ_QUAL] = 30

        status = self.patient_analyzer.filter_reads(quality_filters)

        assert status # assert that status is true (filtering has occured)

        seq_rec_obj = Bio.SeqIO.parse(FILTERED_DIR, "fastq")

        for seq in seq_rec_obj:
            avg_score = quality_filters[MEAN_CUTOFF] + 1
            avg_score = (float(sum(seq.letter_annotations['phred_quality'])) /
                         float(len(seq.letter_annotations['phred_quality'])))

            # check that length and score are both over threshold
            assert len(seq.seq) >= quality_filters[LENGTH_CUTOFF] and \
                avg_score >= quality_filters[MEAN_CUTOFF]

        # patient_analyzer.filter_reads calls quality_control.filter_reads
        # more tests for filtering reads can be found in test_quality_control

    def test_generate_bam(self):
        assert not os.path.isfile("%s/align.bam" % OUTPUT_DIR)
        fasta_id = os.path.basename(self.patient_analyzer.reads).split('.')[0]
        self.patient_analyzer.generate_bam(fasta_id)

        assert os.path.isfile("%s/align.bam" % OUTPUT_DIR)

        if os.path.isfile("%s/align.bam" % OUTPUT_DIR):
            os.remove("%s/align.bam" % OUTPUT_DIR)

    def test_analyze_reads(self):
        quality_filters = defaultdict(dict)

        # test defaults
        variant_filters = defaultdict(dict)
        variant_filters[ERROR_RATE] = 0.0021
        variant_filters[MIN_VARIANT_QUAL] = 30
        variant_filters[MIN_DP] = 100
        variant_filters[MIN_AC] = 5
        variant_filters[MIN_FREQ] = 0.01

        generate_consensus = True
        reporting_threshold = 20

        fasta_id = os.path.basename(self.patient_analyzer.reads).split('.')[0]

        self.patient_analyzer.analyze_reads(fasta_id,
                                            variant_filters,
                                            reporting_threshold,
                                            generate_consensus)

        assert os.path.isfile("%s/consensus.fasta" % OUTPUT_DIR)
        assert os.path.isfile("%s/coverage_file.csv" % OUTPUT_DIR)
        assert os.path.isfile("%s/dr_report.csv" % OUTPUT_DIR)
        assert os.path.isfile("%s/hydra.vcf" % OUTPUT_DIR)
        assert os.path.isfile("%s/mutation_report.aavf" % OUTPUT_DIR)
        assert os.path.isfile("%s/stats.txt" % OUTPUT_DIR)

        assert not os.path.isfile("%s/tmp.bam" % OUTPUT_DIR)
        assert not os.path.isfile("%s/tmp.sam" % OUTPUT_DIR)

        # Remove the output directory so that multiple tests (python 2.x, 3.x, etc.)
        # can run without erroring out with "File/directory exists"
        shutil.rmtree("%s" % OUTPUT_DIR)

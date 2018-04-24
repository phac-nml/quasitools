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

import os
import shutil
from collections import defaultdict
import pytest
from quasitools.patient_analyzer import PatientAnalyzer
import Bio.SeqIO


TEST_PATH = os.path.dirname(os.path.abspath(__file__))
READS = TEST_PATH + "/data/reads_w_K103N.fastq"
FORWARD = TEST_PATH + "/data/forward.fastq"
REVERSE = TEST_PATH + "data/reverse.fastq"
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


    def test_combine_reads(self):
        # Combine the fwd and reverse into one file
        reads = "%s/combined_reads.fastq" % OUTPUT_DIR
        cat_cmd = "cat %s %s > %s" % (FORWARD, REVERSE, reads)
        os.system(cat_cmd)

        assert os.path.isfile("%s/combined_reads.fastq" % OUTPUT_DIR)

    def test_generate_bam(self):
        assert not os.path.isfile("%s/align.bam" % OUTPUT_DIR)
        fasta_id = os.path.basename(self.patient_analyzer.reads).split('.')[0]

        self.patient_analyzer.generate_bam(fasta_id)

        assert os.path.isfile("%s/align.bam" % OUTPUT_DIR)

    def test_analyze_reads(self):
        quality_filters = defaultdict(dict)

        length_cutoff = 100
        score_cutoff = 30
        min_qual = 30

        # read_filters["TRIMMING"] = "trimming"
        # read_filters["MASKING"] = "masking"

        quality_filters["length_cutoff"] = length_cutoff
        quality_filters["mean_cutoff"] = score_cutoff
        quality_filters["ns"] = True
        quality_filters["minimum_quality"] = min_qual

        # test defaults
        variant_filters = defaultdict(dict)
        variant_filters["error_rate"] = 0.0021
        variant_filters["min_qual"] = min_qual
        variant_filters["min_dp"] = 100
        variant_filters["min_ac"] = 5
        variant_filters["min_freq"] = 0.01

        generate_consensus = True
        reporting_threshold = 20

        fasta_id = os.path.basename(self.patient_analyzer.reads).split('.')[0]

        self.patient_analyzer.analyze_reads(fasta_id, quality_filters,
                                            variant_filters,
                                            reporting_threshold,
                                            generate_consensus)

        assert os.path.isfile("%s/consensus.fasta" % OUTPUT_DIR)
        assert os.path.isfile("%s/coverage_file.csv" % OUTPUT_DIR)
        assert os.path.isfile("%s/dr_report.csv" % OUTPUT_DIR)
        assert os.path.isfile("%s/hydra.vcf" % OUTPUT_DIR)
        assert os.path.isfile("%s/mutation_report.hmcf" % OUTPUT_DIR)
        assert os.path.isfile("%s/stats.txt" % OUTPUT_DIR)

        assert not os.path.isfile("%s/tmp.bam" % OUTPUT_DIR)
        assert not os.path.isfile("%s/tmp.sam" % OUTPUT_DIR)

        # Remove the output directory so that multiple tests (python 2.x, 3.x, etc.)
        # can run without erroring out with "File/directory exists"
        shutil.rmtree("%s" % OUTPUT_DIR)

"""
Copyright Government of Canada 2017 - 2026

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
import subprocess
import shutil
import tempfile
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
        self.tmp_dir = tempfile.mkdtemp()
        self.output_dir = os.path.join(self.tmp_dir, "test_patient_analyzer_output")
        self.filtered_dir = os.path.join(self.output_dir, "filtered.fastq")
        self.patient_analyzer = PatientAnalyzer(id="test",reads=READS,
                                 reference=REFERENCE,
                                 output_dir=self.output_dir,
                                 BED4_file=BED4_FILE,
                                 mutation_db=MUTATION_DB,
                                 quiet=False, consensus_pct=20)

    @classmethod
    def teardown_class(self):
        shutil.rmtree(self.tmp_dir, ignore_errors=True)

    def test_combine_reads(self):
        # Combine the fwd and reverse into one file
        reads = "%s/combined_reads.fastq" % self.output_dir
        cat_cmd = "cat %s %s > %s" % (FORWARD, REVERSE, reads)
        os.system(cat_cmd)

        assert os.path.isfile("%s/combined_reads.fastq" % self.output_dir)

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

        seq_rec_obj = Bio.SeqIO.parse(self.filtered_dir, "fastq")

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
        assert not os.path.isfile("%s/align.bam" % self.output_dir)
        fasta_id = os.path.basename(self.patient_analyzer.reads).split('.')[0]
        self.patient_analyzer.generate_bam(fasta_id)

        assert os.path.isfile("%s/align.bam" % self.output_dir)

        if os.path.isfile("%s/align.bam" % self.output_dir):
            os.remove("%s/align.bam" % self.output_dir)

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

        assert os.path.isfile("%s/consensus.fasta" % self.output_dir)
        assert os.path.isfile("%s/coverage_file.csv" % self.output_dir)
        assert os.path.isfile("%s/dr_report.csv" % self.output_dir)
        assert os.path.isfile("%s/hydra.vcf" % self.output_dir)
        assert os.path.isfile("%s/mutation_report.aavf" % self.output_dir)
        assert os.path.isfile("%s/stats.txt" % self.output_dir)

        assert not os.path.isfile("%s/tmp.bam" % self.output_dir)
        assert not os.path.isfile("%s/tmp.sam" % self.output_dir)

    def test_analyze_reads_reraises_bam_error(self, monkeypatch):
        def fail(_):
            raise RuntimeError("bam failed")

        monkeypatch.setattr(self.patient_analyzer, "generate_bam", fail)
        with pytest.raises(RuntimeError, match="bam failed"):
            self.patient_analyzer.analyze_reads(
                "sample", defaultdict(dict), 1, False)

    def test_generate_bam_reraises_tool_errors_and_closes_files(
            self, monkeypatch):
        class Process:
            returncode = 0

        monkeypatch.setattr(
            subprocess,
            "Popen",
            lambda *args, **kwargs: Process(),
        )
        analyzer = PatientAnalyzer.__new__(PatientAnalyzer)
        analyzer.output_dir = self.output_dir
        analyzer.reference = REFERENCE
        analyzer.filtered_reads_dir = self.filtered_dir
        for failing_tool in ("bowtie2", "samtools view", "samtools sort",
                             "samtools index"):
            def process_output(proc, log_fn, name, failing_tool=failing_tool):
                if name == failing_tool:
                    raise RuntimeError(name)

            monkeypatch.setattr(
                analyzer, "process_tool_output", process_output)
            with pytest.raises(RuntimeError, match=failing_tool):
                analyzer.generate_bam("sample")

        monkeypatch.setattr(analyzer, "process_tool_output", lambda *args: None)
        monkeypatch.setattr(os, "unlink", lambda path: None)
        assert analyzer.generate_bam("sample").endswith("align.bam")

    def test_process_tool_output_records_errors(self, tmp_path):
        class Process:
            returncode = 1

            def communicate(self):
                return b'', b'failure'

        log = tmp_path / 'tool.log'
        with pytest.raises(Exception, match='tool returned'):
            self.patient_analyzer.process_tool_output(
                Process(), str(log), 'tool')
        assert 'failure' in log.read_text()

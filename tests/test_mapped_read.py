"""
Copyright Government of Canada 2015-2026

Written by: Eric Enns, National Microbiology Laboratory, Public Health Agency of Canada

Licensed under the Apache License, Version 2.0 (the "License"); you may not use
this work except in compliance with the License. You may obtain a copy of the
License at:

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software distributed
under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR
CONDITIONS OF ANY KIND, either express or implied. See the License for the
specific language governing permissions and limitations under the License.
"""

import pytest
from quasitools.mapped_read import MappedRead, MappedReadCollection
from quasitools.reference import Reference
from quasitools.parsers.mapped_read_parser import (
    parse_haplotypes_from_bam_range,
    parse_haplotypes_from_fasta,
    get_index_in_list,
    is_consecutive_list,
    parse_mapped_reads_from_bam,
    parse_pileup_from_fasta,
    create_positional_array,
)

class TestMappedRead:
    @classmethod
    def setup_class(self):
        self.mr = MappedRead('read1', 10, 300, {'110': 'g', '300': 'tagt'}, 100, 390, '+')

    def test_query_length(self):
        assert self.mr.query_length() == 291

    def test_codon_start(self):
        assert self.mr.codon_start(0) == 102
        assert self.mr.codon_start(1) == 100
        assert self.mr.codon_start(2) == 101

    def test_codon_end(self):
        assert self.mr.codon_end(0) == 389
        assert self.mr.codon_end(1) == 390
        assert self.mr.codon_end(2) == 388

class TestMappedReadCollection:
    @classmethod
    def setup_class(self):
        self.reference = Reference('test', 'AGCTTAGCTAAGCTACCTATATCTTGGTCTTGGCCG')
        self.mapped_reads_obj = MappedReadCollection(self.reference)

        for i in range(0,5):
            self.mapped_reads_obj.mapped_reads["{0}_{1}".format('read%s' % i, 1)] = MappedRead('read%i' % i, 0, 34, {10: '-'}, 0, 35, '+' )

        for i in range(5,100):
            self.mapped_reads_obj.mapped_reads["{0}_{1}".format('read%s' % i, 1)] = MappedRead('read%i' % i, 0, 34, {}, 0, 35, '+' )

    def test_from_bam(self):
        reference = Reference('hxb2_pol', 'CCTCAGGTCACTCTTTGGCAACGACCCCTCGTCACAATAAAGATAGGGGGGCAACTAAAGGAAGCTCTATTAGATACAGGAGCAGATGATACAGTATTAGAAGAAATGAGTTTGCCAGGAAGATGGAAACCAAAAATGATAGGGGGAATTGGAGGTTTTATCAAAGTAAGACAGTATGATCAGATACTCATAGAAATCTGTGGACATAAAGCTATAGGTACAGTATTAGTAGGACCTACACCTGTCAACATAATTGGAAGAAATCTGTTGACTCAGATTGGTTGCACTTTAAATTTTCCCATTAGCCCTATTGAGACTGTACCAGTAAAATTAAAGCCAGGAATGGATGGCCCAAAAGTTAAACAATGGCCATTGACAGAAGAAAAAATAAAAGCATTAGTAGAAATTTGTACAGAGATGGAAAAGGAAGGGAAAATTTCAAAAATTGGGCCTGAAAATCCATACAATACTCCAGTATTTGCCATAAAGAAAAAAGACAGTACTAAATGGAGAAAATTAGTAGATTTCAGAGAACTTAATAAGAGAACTCAAGACTTCTGGGAAGTTCAATTAGGAATACCACATCCCGCAGGGTTAAAAAAGAAAAAATCAGTAACAGTACTGGATGTGGGTGATGCATATTTTTCAGTTCCCTTAGATGAAGACTTCAGGAAGTATACTGCATTTACCATACCTAGTATAAACAATGAGACACCAGGGATTAGATATCAGTACAATGTGCTTCCACAGGGATGGAAAGGATCACCAGCAATATTCCAAAGTAGCATGACAAAAATCTTAGAGCCTTTTAGAAAACAAAATCCAGACATAGTTATCTATCAATACATGGATGATTTGTATGTAGGATCTGACTTAGAAATAGGGCAGCATAGAACAAAAATAGAGGAGCTGAGACAACATCTGTTGAGGTGGGGACTTACCACACCAGACAAAAAACATCAGAAAGAACCTCCATTCCTTTGGATGGGTTATGAACTCCATCCTGATAAATGGACAGTACAGCCTATAGTGCTGCCAGAAAAAGACAGCTGGACTGTCAATGACATACAGAAGTTAGTGGGGAAATTGAATTGGGCAAGTCAGATTTACCCAGGGATTAAAGTAAGGCAATTATGTAAACTCCTTAGAGGAACCAAAGCACTAACAGAAGTAATACCACTAACAGAAGAAGCAGAGCTAGAACTGGCAGAAAACAGAGAGATTCTAAAAGAACCAGTACATGGAGTGTATTATGACCCATCAAAAGACTTAATAGCAGAAATACAGAAGCAGGGGCAAGGCCAATGGACATATCAAATTTATCAAGAGCCATTTAAAAATCTGAAAACAGGAAAATATGCAAGAATGAGGGGTGCCCACACTAATGATGTAAAACAATTAACAGAGGCAGTGCAAAAAATAACCACAGAAAGCATAGTAATATGGGGAAAGACTCCTAAATTTAAACTGCCCATACAAAAGGAAACATGGGAAACATGGTGGACAGAGTATTGGCAAGCCACCTGGATTCCTGAGTGGGAGTTTGTTAATACCCCTCCCTTAGTGAAATTATGGTACCAGTTAGAGAAAGAACCCATAGTAGGAGCAGAAACCTTCTATGTAGATGGGGCAGCTAACAGGGAGACTAAATTAGGAAAAGCAGGATATGTTACTAATAGAGGAAGACAAAAAGTTGTCACCCTAACTGACACAACAAATCAGAAGACTGAGTTACAAGCAATTTATCTAGCTTTGCAGGATTCGGGATTAGAAGTAAACATAGTAACAGACTCACAATATGCATTAGGAATCATTCAAGCACAACCAGATCAAAGTGAATCAGAGTTAGTCAATCAAATAATAGAGCAGTTAATAAAAAAGGAAAAGGTCTATCTGGCATGGGTACCAGCACACAAAGGAATTGGAGGAAATGAACAAGTAGATAAATTAGTCAGTGCTGGAATCAGGAAAGTACTATTTTTAGATGGAATAGATAAGGCCCAAGATGAACATGAGAAATATCACAGTAATTGGAGAGCAATGGCTAGTGATTTTAACCTGCCACCTGTAGTAGCAAAAGAAATAGTAGCCAGCTGTGATAAATGTCAGCTAAAAGGAGAAGCCATGCATGGACAAGTAGACTGTAGTCCAGGAATATGGCAACTAGATTGTACACATTTAGAAGGAAAAGTTATCCTGGTAGCAGTTCATGTAGCCAGTGGATATATAGAAGCAGAAGTTATTCCAGCAGAAACAGGGCAGGAAACAGCATATTTTCTTTTAAAATTAGCAGGAAGATGGCCAGTAAAAACAATACATACTGACAATGGCAGCAATTTCACCGGTGCTACGGTTAGGGCCGCCTGTTGGTGGGCGGGAATCAAGCAGGAATTTGGAATTCCCTACAATCCCCAAAGTCAAGGAGTAGTAGAATCTATGAATAAAGAATTAAAGAAAATTATAGGACAGGTAAGAGATCAGGCTGAACATCTTAAGACAGCAGTACAAATGGCAGTATTCATCCACAATTTTAAAAGAAAAGGGGGGATTGGGGGGTACAGTGCAGGGGAAAGAATAGTAGACATAATAGCAACAGACATACAAACTAAAGAATTACAAAAACAAATTACAAAAATTCAAAATTTTCGGGTTTATTACAGGGACAGCAGAAATCCACTTTGGAAAGGACCAGCAAAGCTCCTCTGGAAAGGTGAAGGGGCAGTAGTAATACAAGATAATAGTGACATAAAAGTAGTGCCAAGAAGAAAAGCAAAGATCATTAGGGATTATGGAAAACAGATGGCAGGTGATGATTGTGTGGCAAGTAGACAGGATGAGGATTAG')
        mrc = parse_mapped_reads_from_bam(reference, 'tests/data/test1.bam')

        assert len(mrc.mapped_reads) == 6308

    def test_pileup(self):
        assert self.mapped_reads_obj.pileup() == [{'A': 100}, {'G': 100}, {'C':100}, {'T':100}, {'T':100}, {'A':100}, {'G':100}, {'C':100}, {'T':100}, {'A':100}, {'A':95, '-':5}, {'G':100}, {'C':100}, {'T':100}, {'A':100}, {'C':100}, {'C':100}, {'T':100}, {'A':100}, {'T':100}, {'A':100}, {'T':100}, {'C':100}, {'T':100}, {'T':100}, {'G':100}, {'G':100}, {'T':100}, {'C':100}, {'T':100}, {'T':100}, {'G':100}, {'G':100}, {'C':100}, {'C':100}, {'G':100}]

    def test_to_consensus(self):
        assert self.mapped_reads_obj.to_consensus(100) == 'AGCTTAGCTAAGCTACCTATATCTTGGTCTTGGCCG'
        assert self.mapped_reads_obj.to_consensus(20) == 'AGCTTAGCTAAGCTACCTATATCTTGGTCTTGGCCG'

    def test_coverage_without_supplied_pileup(self):
        coverage = self.mapped_reads_obj.coverage()
        assert coverage == [100] * 10 + [95] + [100] * 25
        assert coverage[10] == 95

    def test_to_consensus_handles_indel_iupac_branches(self, monkeypatch):
        monkeypatch.setattr(
            self.mapped_reads_obj,
            'pileup',
            lambda indels=False: [{'AC': 1}, {'A': 1}, {'AC': 1}],
        )
        monkeypatch.setattr(
            self.mapped_reads_obj,
            'coverage',
            lambda pileup: [100, 100, 100],
        )

        assert self.mapped_reads_obj.to_consensus(100) == 'AAA'
        monkeypatch.setattr(
            self.mapped_reads_obj,
            'pileup',
            lambda indels=False: [{'A': 1}, {'A': 1}, {'ACGT': 1}],
        )
        assert self.mapped_reads_obj.to_consensus(100) == 'AAA'

    def test_mask_unconfident_differences(self):
        read = MappedRead('read', 0, 1, {0: 'a', 1: 'n'}, 0, 2, '+')
        self.mapped_reads_obj.mapped_reads = {'read': read}
        variants = type('Variants', (), {
            'variants': {'test': {
                1: {'a': type('Variant', (), {'filter': 'FAIL'})()},
                2: {'n': type('Variant', (), {'filter': 'PASS'})()},
            }}
        })()

        self.mapped_reads_obj.mask_unconfident_differences(variants)

        assert read.differences == {0: 'a', 1: 'n'}


def test_parse_pileup_and_haplotypes_from_fasta(tmp_path):
    fasta = tmp_path / 'reads.fasta'
    fasta.write_text('>one\nA-C\n>two\nA-C\n')

    assert parse_pileup_from_fasta(str(fasta)).pileup == [
        {'A': 2}, {}, {'C': 2}]
    haplotypes = parse_haplotypes_from_fasta(str(fasta))
    assert len(haplotypes) == 1
    assert haplotypes[0].count == 2


def test_parse_haplotypes_from_bam_range_skips_nonconsecutive_read(monkeypatch):
    class Read:
        reference_start = 0
        cigartuples = ()
        query_alignment_sequence = 'ABC'

        def get_overlap(self, start, end):
            return end - start

    class Samfile:
        def fetch(self, reference, start, end):
            return [Read()]

    monkeypatch.setattr(
        'quasitools.parsers.mapped_read_parser.create_positional_array',
        lambda _: [0, 2, 3],
    )

    assert parse_haplotypes_from_bam_range(
        Samfile(), Reference('test', 'ACG'), None, 0, 3) == []


def test_haplotype_range_skips_unmapped_and_deleted_reads(monkeypatch):
    class Read:
        reference_start = 0
        cigartuples = ()
        query_alignment_sequence = 'ABC'

        def get_overlap(self, start, end):
            return end - start

    class Samfile:
        def fetch(self, reference, start, end):
            return [Read()]

    monkeypatch.setattr(
        'quasitools.parsers.mapped_read_parser.create_positional_array',
        lambda _: [0, 0, 1],
    )
    assert parse_haplotypes_from_bam_range(
        Samfile(), Reference('test', 'ACG'), None, 0, 2) == []

    monkeypatch.setattr(
        'quasitools.parsers.mapped_read_parser.create_positional_array',
        lambda _: [0, 1],
    )
    assert parse_haplotypes_from_bam_range(
        Samfile(), Reference('test', 'ACG'), None, 0, 3) == []

    monkeypatch.setattr(
        'quasitools.parsers.mapped_read_parser.create_positional_array',
        lambda _: [1, 2],
    )
    assert parse_haplotypes_from_bam_range(
        Samfile(), Reference('test', 'ACG'), None, 0, 2) == []

    monkeypatch.setattr(
        'quasitools.parsers.mapped_read_parser.create_positional_array',
        lambda _: [0, 2, 3],
    )
    indices = iter([0, 2])
    monkeypatch.setattr(
        'quasitools.parsers.mapped_read_parser.get_index_in_list',
        lambda values, value: next(indices),
    )
    assert parse_haplotypes_from_bam_range(
        Samfile(), Reference('test', 'ACG'), None, 0, 3) == []


def test_haplotype_parser_helpers():
    assert get_index_in_list([1, 2], 3) == -1
    assert not is_consecutive_list([1, 3])


def test_create_positional_array_handles_insertions_and_deletions():
    assert create_positional_array([(0, 2), (1, 2), (2, 1), (0, 1)]) == [
        0, 1, 1, 1, 3]

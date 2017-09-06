"""
Copyright Government of Canada 2015

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

import decimal
import itertools
import operator
import pysam

from quasitools.utilities import *

class MappedRead(object):
    def __init__(self, seq_id, query_start, query_end, differences, ref_start, ref_end, overlap, identity, strand):
        self.seq_id      = seq_id
        self.query_start = query_start
        self.query_end   = query_end
        self.differences = differences
        self.ref_start   = ref_start
        self.ref_end     = ref_end
        self.overlap     = overlap
        self.identity    = identity
        self.strand      = strand

    def query_length(self):
        """Calculate and return the length of read mapped to the reference.

        >>> mr = MappedRead('read1', 10, 300, {'110': 'g', '300': 'tagt'}, 100, 390, 97, 99.3127, '+')
        >>> print(mr.query_length())
        291
        """
        return self.query_end + 1 - self.query_start

    def codon_start(self, frame):
        """Calculate and return the first position of the first codon for the given frame in reference base position.

        >>> mr = MappedRead('read1', 10, 300, {'110': 'g', '300': 'tagt'}, 100, 390, 96.6667, 99.3103, '+')
        >>> print(mr.codon_start(0))
        102
        >>> print(mr.codon_start(1))
        100
        >>> print(mr.codon_start(2))
        101
        """
        codon_start = self.ref_start

        while codon_start % 3 != frame:
            codon_start += 1

        return codon_start

    def codon_end(self, frame):
        """Calculate and return the last position of the last codon for the given frame in reference base position.

        >>> mr = MappedRead('read1', 10, 300, {'110': 'g', '300': 'tagt'}, 100, 390, 96.6667, 99.3103, '+')
        >>> print(mr.codon_end(0))
        389
        >>> print(mr.codon_end(1))
        390
        >>> print(mr.codon_end(2))
        388
        """
        codon_end = self.ref_end

        while (codon_end-2) % 3 != frame:
            codon_end -= 1

        return codon_end

class MappedReadCollection(object):
    def __init__(self, reference):
        self.mapped_reads = {}
        self.reference = reference

    @classmethod
    def from_bam(cls, reference, overlap_cutoff, identity_cutoff, bam):
        """Build the MappedReadCollection object from a bam."""
        obj = cls(reference)

        sam = pysam.AlignmentFile(bam, "rb")

        for alignment in sam.fetch(reference=reference.name):
            overlap = decimal.Decimal(alignment.query_alignment_length) / alignment.query_length * 100

            if overlap >= overlap_cutoff:
                padded_alignment = sam_alignment_to_padded_alignment(alignment, reference)

                matches = 0
                if alignment.has_tag('NM'):
                    matches = alignment.query_length - alignment.get_tag('NM')
                else:
                    matches = padded_alignment[1].count('|')

                #TODO not an accurate identity, as I am not distinguishing between an alignment match vs a sequence match
                identity = decimal.Decimal(matches) / len(padded_alignment[1]) * 100

                if identity >= identity_cutoff:
                    direct = '+'
                    if alignment.flag & 16:
                        direct = '-'
                    #TODO only calculate differences when identity < 100
                    differences = pairwise_alignment_to_differences(padded_alignment[0], padded_alignment[2], alignment.reference_start)

                    mapped_read = MappedRead(alignment.query_name, alignment.query_alignment_start, alignment.query_alignment_end-1,
                            differences, alignment.reference_start, alignment.reference_end-1, overlap, identity, direct)
                    obj.mapped_reads["{0}_{1}".format(alignment.query_name, '1' if direct == '+' else '2')] = mapped_read

        return obj

    def pileup(self):
        """Build and return a pileup from the object."""
        pileup = [{} for i in range(0,len(self.reference.seq))]

        for name, mapped_read in self.mapped_reads.items():
            for i in range(mapped_read.ref_start,mapped_read.ref_end+1):
                if i not in mapped_read.differences:
                    pileup[i][self.reference.sub_seq(i,i).upper()] = pileup[i].get(self.reference.sub_seq(i,i).upper(), 0) + 1
                elif len(mapped_read.differences[i]) == 1:
                    if mapped_read.differences[i] == '-':
                        pileup[i]['-'] = pileup[i].get('-', 0) + 1
                    else:
                        pileup[i][mapped_read.differences[i].upper()] = pileup[i].get(mapped_read.differences[i].upper(), 0) + 1
                else:
                    difference = mapped_read.differences[i]
                    difference = difference.replace('.', self.reference.sub_seq(i,i).upper())
                    pileup[i][difference.upper()] = pileup[i].get(difference.upper(), 0) + 1

        return pileup


    def to_consensus(self, percentage):
        """Generates and returns a consensus sequence

        >>> from quasitools.reference import Reference
        >>> ref = Reference('hxb2_pol', 'CCTCAGGTCACTCTTTGGCAACGACCCCTCGTCACAATAAAGATAGGGGGGCAACTAAAGGAAGCTCTATTAGATACAGGAGCAGATGATACAGTATTAGAAGAAATGAGTTTGCCAGGAAGATGGAAACCAAAAATGATAGGGGGAATTGGAGGTTTTATCAAAGTAAGACAGTATGATCAGATACTCATAGAAATCTGTGGACATAAAGCTATAGGTACAGTATTAGTAGGACCTACACCTGTCAACATAATTGGAAGAAATCTGTTGACTCAGATTGGTTGCACTTTAAATTTTCCCATTAGCCCTATTGAGACTGTACCAGTAAAATTAAAGCCAGGAATGGATGGCCCAAAAGTTAAACAATGGCCATTGACAGAAGAAAAAATAAAAGCATTAGTAGAAATTTGTACAGAGATGGAAAAGGAAGGGAAAATTTCAAAAATTGGGCCTGAAAATCCATACAATACTCCAGTATTTGCCATAAAGAAAAAAGACAGTACTAAATGGAGAAAATTAGTAGATTTCAGAGAACTTAATAAGAGAACTCAAGACTTCTGGGAAGTTCAATTAGGAATACCACATCCCGCAGGGTTAAAAAAGAAAAAATCAGTAACAGTACTGGATGTGGGTGATGCATATTTTTCAGTTCCCTTAGATGAAGACTTCAGGAAGTATACTGCATTTACCATACCTAGTATAAACAATGAGACACCAGGGATTAGATATCAGTACAATGTGCTTCCACAGGGATGGAAAGGATCACCAGCAATATTCCAAAGTAGCATGACAAAAATCTTAGAGCCTTTTAGAAAACAAAATCCAGACATAGTTATCTATCAATACATGGATGATTTGTATGTAGGATCTGACTTAGAAATAGGGCAGCATAGAACAAAAATAGAGGAGCTGAGACAACATCTGTTGAGGTGGGGACTTACCACACCAGACAAAAAACATCAGAAAGAACCTCCATTCCTTTGGATGGGTTATGAACTCCATCCTGATAAATGGACAGTACAGCCTATAGTGCTGCCAGAAAAAGACAGCTGGACTGTCAATGACATACAGAAGTTAGTGGGGAAATTGAATTGGGCAAGTCAGATTTACCCAGGGATTAAAGTAAGGCAATTATGTAAACTCCTTAGAGGAACCAAAGCACTAACAGAAGTAATACCACTAACAGAAGAAGCAGAGCTAGAACTGGCAGAAAACAGAGAGATTCTAAAAGAACCAGTACATGGAGTGTATTATGACCCATCAAAAGACTTAATAGCAGAAATACAGAAGCAGGGGCAAGGCCAATGGACATATCAAATTTATCAAGAGCCATTTAAAAATCTGAAAACAGGAAAATATGCAAGAATGAGGGGTGCCCACACTAATGATGTAAAACAATTAACAGAGGCAGTGCAAAAAATAACCACAGAAAGCATAGTAATATGGGGAAAGACTCCTAAATTTAAACTGCCCATACAAAAGGAAACATGGGAAACATGGTGGACAGAGTATTGGCAAGCCACCTGGATTCCTGAGTGGGAGTTTGTTAATACCCCTCCCTTAGTGAAATTATGGTACCAGTTAGAGAAAGAACCCATAGTAGGAGCAGAAACCTTCTATGTAGATGGGGCAGCTAACAGGGAGACTAAATTAGGAAAAGCAGGATATGTTACTAATAGAGGAAGACAAAAAGTTGTCACCCTAACTGACACAACAAATCAGAAGACTGAGTTACAAGCAATTTATCTAGCTTTGCAGGATTCGGGATTAGAAGTAAACATAGTAACAGACTCACAATATGCATTAGGAATCATTCAAGCACAACCAGATCAAAGTGAATCAGAGTTAGTCAATCAAATAATAGAGCAGTTAATAAAAAAGGAAAAGGTCTATCTGGCATGGGTACCAGCACACAAAGGAATTGGAGGAAATGAACAAGTAGATAAATTAGTCAGTGCTGGAATCAGGAAAGTACTATTTTTAGATGGAATAGATAAGGCCCAAGATGAACATGAGAAATATCACAGTAATTGGAGAGCAATGGCTAGTGATTTTAACCTGCCACCTGTAGTAGCAAAAGAAATAGTAGCCAGCTGTGATAAATGTCAGCTAAAAGGAGAAGCCATGCATGGACAAGTAGACTGTAGTCCAGGAATATGGCAACTAGATTGTACACATTTAGAAGGAAAAGTTATCCTGGTAGCAGTTCATGTAGCCAGTGGATATATAGAAGCAGAAGTTATTCCAGCAGAAACAGGGCAGGAAACAGCATATTTTCTTTTAAAATTAGCAGGAAGATGGCCAGTAAAAACAATACATACTGACAATGGCAGCAATTTCACCGGTGCTACGGTTAGGGCCGCCTGTTGGTGGGCGGGAATCAAGCAGGAATTTGGAATTCCCTACAATCCCCAAAGTCAAGGAGTAGTAGAATCTATGAATAAAGAATTAAAGAAAATTATAGGACAGGTAAGAGATCAGGCTGAACATCTTAAGACAGCAGTACAAATGGCAGTATTCATCCACAATTTTAAAAGAAAAGGGGGGATTGGGGGGTACAGTGCAGGGGAAAGAATAGTAGACATAATAGCAACAGACATACAAACTAAAGAATTACAAAAACAAATTACAAAAATTCAAAATTTTCGGGTTTATTACAGGGACAGCAGAAATCCACTTTGGAAAGGACCAGCAAAGCTCCTCTGGAAAGGTGAAGGGGCAGTAGTAATACAAGATAATAGTGACATAAAAGTAGTGCCAAGAAGAAAAGCAAAGATCATTAGGGATTATGGAAAACAGATGGCAGGTGATGATTGTGTGGCAAGTAGACAGGATGAGGATTAG')
        >>> mrs = MappedReadCollection.from_bam(ref, 65, 75, 'tests/data/test1.bam')
        >>> print(mrs.to_consensus(20))
        GAAGAAATGAGTTTGCCAGGAAGATGGAAACCAAAAATGATAGGGGGAATTGGAGGTTTTATCAAAGTAAGACAGTATGATCAGATACTCATAGAAATCTGTGGACATAAAGCTATAGGTACAGTATTAGTAGGACCTACACCTGTCAACATAATTGGAAGAAATCTGTTGACTCAGATTGGTTGCACTTTAAATTTTCCCATTAGCCCTATTGAGACTGTACCAGTAAAATTAAAGCCAGGAATGGATGGCCCAAAAGTTAAACAATGGCCATTGACAGAAGAAAAAATAAAAGCATTAGTAGAAATTTGTACAGAGATGGAAAAGGAAGGGAAAATTTCAAAAATTGGGCCTGAAAATCCATACAATACTCCAGTATTTGCCATAAAGAAAAAAGACAGTACTAAATGGAGAAAATTAGTAGATTTCAGAGAACTTAATAAGAGAACTCAAGACTTCTGGGAAGTTCAATTAGGAATACCACATCCCGCAGGGTTAAAAAAGAAAAAATCAGTAACAGTACTGGATGTGGGTGATGCATATTTTTCAGTTCCCTTAGATGAAGACTTCAGGAAGTATACTGCATTTACCATACCTAGTATAAACAATGAGACACCAGGGATTAGATATCAGTACAATGTGCTTCCACAGGGATGGAAAGGATCACCAGCAATATTCCAAAGTAGCATGACAAAAATCTTAGAGCCTTTTAGAAAACAAAATCCAGACATAGTTATCTATCAATACATGGATGATTTGTATGTAGGATCTGACTTAGAAATAGGGCAGCATAGAACAAAAATAGAGGAGCTGAGACAACATCTGTTGAGGTGGGGACTTACCACACCAGACAAAAAACATCAGAAAGAACCTCCATTCCTTTGGATGGGTTATGAACTCCATCCTGATAAATGGACAGTACAGCCTATAGTGCTGCCAGAAAAAGACAGCTGGACTGTCAATGACATACAGAAGTTAGTGGGGAAATTGAATTGGGCAAGTCAGATTTACCCAGGGATTAAAGTAAGGCAATTATGTAAACTCCTTAGAGGAACCAAAGCACTAACAGAAGTAATACCACTAACAGAAGAAGCAGAGCTAGAACTGGCAGAAAACAGAGAGATTCTAAAAGAACCAGTACATGGAGTGTATTATGACCCATCAAAAGACTTAATAGCAGAAATACAGAAGCAGGGGCAAGGCCAATGGACATATCAAATTTATCAAGAGCCATTTAAAAATCTGAAAACAGGAAAATATGCAAGAATGAGGGGTGCCCACACTAATGATGTAAAACAATTAACAGAGGCAGTGCAAAAAATAACCACAGAAAGCATAGTAATATGGGGAAAGACTCCTAAATTTAAACTGCCCATACAAAAGGAAACATGGGAAACATGGTGGACAGAGTATTGGCAAGCCACCTGGATTCCTGAGTGGGAGTTTGTTAATACCCCTCCCTTAGTGAAATTATGGTACCAGTTAGAGAAAGAACCCATAGTAGGAGCAGAAACCTTCTATGTAGATGGGGCAGCTAACAGGGAGACTAAATTAGGAAAAGCAGGATATGTTACTAATAGAGGAAGACAAAAAGTTGTCACCCTAACTGACACAACAAATCAGAAGACTGAGTTACAAGCAATTTATCTAGCTTTGCAGGATTCGGGATTAGAAGTAAACATAGTAACAGACTCACAATATGCATTAGGAATCATTCAAGCACAACCAGATCAAAGTGAATCAGAGTTAGTCAATCAAATAATAGAGCAGTTAATAAAAAAGGAAAAGGTCTATCTGGCATGGGTACCAGCACACAAAGGAATTGGAGGAAATGAACAAGTAGATAAATTAGTCAGTGCTGGAATCAGGAAAGTACTATTTTTAGATGGAATAGATAAGGCCCAAGATGAACATGAGAAATATCACAGTAATTGGAGAGCAATGGCTAGTGATTTTAACCTGCCACCTGTAGTAGCAAAAGAAATAGTAGCCAGCTGTGATAAATGTCAGCTAAAAGGAGAAGCCATGCATGGACAAGTAGACTGTAGTCCAGGAATATGGCAACTAGATTGTACACATTTAGAAGGAAAAGTTATCCTGGTAGCAGTTCATGTAGCCAGTGGATATATAGAAGCAGAAGTTATTCCAGCAGAAACAGGGCAGGAAACAGCATATTTTCTTTTAAAATTAGCAGGAAGATGGCCAGTAAAAACAATACATACTGACAATGGCAGCAATTTCACCGGTGCTACGGTTAGGGCCGCCTGTTGGTGGGCGGGAATCAAGCAGGAATTTGGAATTCCCTACAATCCCCAAAGTCAAGGAGTAGTAGAATCTATGAATAAAGAATTAAAGAAAATTATAGGACAGGTAAGAGATCAGGCTGAACATCTTAAGACAGCAGTACAAATGGCAGTATTCATCCACAATTTTAAAAGAAAAGGGGGGATTGGGGGGTACAGTGCAGGGGAAAGAATAGTAGACATAATAGCAACAGACATACAAACTAAAGAATTACAAAAACAAATTACAAAAATTCAAAATTTTCGGGTTTATTACAGGGACAGCAGAAATCCACTTTGGAAAGGACCAGCAAAGCTCCTCTGGAAAGGTGAAGGGGCAGTAGTAATACAAGATAATAGTGACATAAAAGTAGTGCCAAGAAGAAAAGCAAAGATCATTAGGGAT
        """
        iupac = {
                'A': 'A',
                'C': 'C',
                'G': 'G',
                'T': 'T',
                'AC': 'M',
                'AG': 'R',
                'AT': 'W',
                'CG': 'S',
                'CT': 'Y',
                'GT': 'K',
                'ACG': 'V',
                'ACT': 'H',
                'AGT': 'D',
                'CGT': 'B',
                'ACGT': 'N'
        }

        pileup = self.pileup()

        # do not include indels in coverage calculations
        coverage = [sum([v if not k.startswith('+') and not k.startswith('-') else 0 for k,v in pileup[pos].items()])for pos in range(0,len(pileup))]
        start, end = None, None
        consensus = ''
        for pos in range(0,len(pileup)):
            if coverage[pos] >= 100:
                if start is None:
                    start = pos
                end = pos

            con_bp = ''

            if coverage[pos] > 0:
                # if percentage is 100, we are generating a consensus sequence of the most common sequence
                if percentage == 100:
                    con_bp_with_indel = sorted([(v, k) for k,v in pileup[pos].items()], reverse=True)[0][1]
                    if len(con_bp_with_indel) > 1:
                        #TODO recalculate consensus base if we can't incorporate insertion
                        if pos % 3 != 2:
                            consensus += con_bp_with_indel[:1]
                        else:
                            if len(con_bp_with_indel) % 3 != 1:
                                consensus += con_bp_with_indel[:1]
                            else:
                                consensus += con_bp_with_indel
                    else:
                        consensus += con_bp_with_indel
                # else we are generating a sanger-like consensus sequence
                else:
                    pilep_wo_indels = [(k,sum([v for k,v in list(g)])) for k,g in itertools.groupby(sorted([(k,v) if len(k) == 0 else (k[:1],v) for k,v in pileup[pos].items()]), operator.itemgetter(0))]
                    for token, token_count in sorted(pileup[pos].items()):
                        if token != '-':
                            if decimal.Decimal(token_count) / coverage[pos] * 100 >= percentage:
                                con_bp += token

                    if len(con_bp.upper()) == 0:
                        consensus += 'N'
                    else:
                        consensus += iupac[con_bp.upper()]
            else:
                consensus += 'n'

        consensus_seq = ''
        if start is not None:
            while start % 3 != 0:
                start += 1
            while end % 3 != 2:
                end -= 1
            consensus_seq = consensus[start:end+1]

        return consensus_seq

if __name__ == '__main__':
    import doctest
    doctest.testmod()

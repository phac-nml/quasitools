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

from quasitools.mapped_read import MappedRead
from quasitools.utilities import *

class MappedReads(object):
    def __init__(self, reference):
        self.mapped_reads = {}
        self.reference = reference

    @classmethod
    def from_bam(cls, reference, overlap_cutoff, identity_cutoff, bam):
        """Build the MappedReads object from a bam."""
        obj = cls(reference)

        sam = pysam.AlignmentFile(bam, "rb")

        for alignment in sam.fetch(reference=reference.name):
            overlap = decimal.Decimal(alignment.query_alignment_length) / alignment.query_length * 100

            if overlap >= overlap_cutoff:
                padded_alignment = sam_alignment_to_padded_alignment(alignment, reference)

                #TODO not an accurate identity, as I am not distinguishing between an alignment match vs a sequence match
                identity = decimal.Decimal(padded_alignment[1].count('|')) / len(padded_alignment[1]) * 100

                if identity >= identity_cutoff:
                    direct = '+'
                    if alignment.flag & 16:
                        direct = '-'
                    #TODO only calculate differences when identity < 100
                    differences = pairwise_alignment_to_differences(padded_alignment[0], padded_alignment[2], alignment.reference_start)

                    mapped_read = MappedRead(alignment.query_name, alignment.query_alignment_start, alignment.query_alignment_end-1,
                            differences, alignment.reference_start, alignment.reference_end-1, overlap, identity, direct)
                    obj.mapped_reads["{0}_{1}".format(alignment.query_name, '1' if direct is '+' else '2')] = mapped_read

        return obj

    def pileup(self):
        """Build and return a pileup from the object."""
        pileup = [{} for i in range(0,len(self.reference.seq))]

        for name, mapped_read in self.mapped_reads.items():
            for i in range(mapped_read.ref_start,mapped_read.ref_end+1):
                if i not in mapped_read.differences:
                    pileup[i][self.reference.sub_seq(i,i).upper()] = pileup[i].get(self.reference.sub_seq(i,i).upper(), 0) + 1
                elif len(mapped_read.differences[i]) == 1:
                    if mapped_read.differences[i] is '-':
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
        >>> mrs = MappedReads.from_bam(ref, 65, 75, 'tests/data/test1.bam')
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
                        if token is not '-':
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

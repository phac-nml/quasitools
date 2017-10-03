"""
Copyright Government of Canada 2015-2017

Written by: Eric Enns, Eric Chubaty, National Microbiology Laboratory, Public Health Agency of Canada

Licensed under the Apache License, Version 2.0 (the "License"); you may not use
this work except in compliance with the License. You may obtain a copy of the
License at:

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software distributed
under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR
CONDITIONS OF ANY KIND, either express or implied. See the License for the
specific language governing permissions and limitations under the License.
"""

import pdb
import decimal
from collections import defaultdict


class MappedRead(object):

    def __init__(self, seq_id, query_start, query_end, differences, ref_start,
                 ref_end, strand):

        self.seq_id = seq_id
        self.query_start = query_start
        self.query_end = query_end
        self.differences = differences
        self.ref_start = ref_start
        self.ref_end = ref_end
        self.strand = strand

    def query_length(self):
        """Calculate and return the length of read mapped to the reference.

        >>> mr = MappedRead('read1', 10, 300, {'110': 'g', '300': 'tagt'}, 100,
        >>>        390, 97, 99.3127, '+')
        >>> print(mr.query_length())
        291
        """
        return self.query_end + 1 - self.query_start

    def codon_start(self, frame):
        """Calculate and return the first position of the first codon for the
        given frame in reference base position.

        >>> mr = MappedRead('read1', 10, 300, {'110': 'g', '300': 'tagt'}, 100,
        >>>        390, 96.6667, 99.3103, '+')
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
        """Calculate and return the last position of the last codon for the
        given frame in reference base position.

        >>> mr = MappedRead('read1', 10, 300, {'110': 'g', '300': 'tagt'}, 100,
        >>>        390, 96.6667, 99.3103, '+')
        >>> print(mr.codon_end(0))
        389
        >>> print(mr.codon_end(1))
        390
        >>> print(mr.codon_end(2))
        388
        """
        codon_end = self.ref_end

        while (codon_end - 2) % 3 != frame:
            codon_end -= 1

        return codon_end


class MappedReadCollection(object):

    def __init__(self, reference):
        self.mapped_reads = {}
        self.reference = reference

    def pileup(self, indels=True):
        """Build and return a pileup from the object."""
        pileup = [{} for i in range(0, len(self.reference.seq))]

        for name, mapped_read in self.mapped_reads.items():
            for i in range(mapped_read.ref_start, mapped_read.ref_end + 1):
                if i not in mapped_read.differences:
                    pileup[i][self.reference.sub_seq(i, i).upper()] = \
                        pileup[i].get(self.reference.sub_seq(i, i).upper(), 0)\
                        + 1
                elif len(mapped_read.differences[i]) == 1:
                    if mapped_read.differences[i] == '-':
                        pileup[i]['-'] = pileup[i].get('-', 0) + 1
                    else:
                        pileup[i][mapped_read.differences[i].upper()] = \
                            pileup[i].get(
                                mapped_read.differences[i].upper(), 0
                        ) + 1
                else:
                    difference = mapped_read.differences[i]
                    difference = difference.replace(
                        '.', self.reference.sub_seq(i, i).upper())
                    if not indels:
                        difference = difference[:1]
                    pileup[i][difference.upper()] = \
                        pileup[i].get(difference.upper(), 0) + 1

        return pileup

    def coverage(self, pileup=None):
        if pileup is None:
            pileup = self.pileup()

        coverage = [0 for pos in range(0, len(pileup))]
        for pos in range(0, len(pileup)):
            for k, v in pileup[pos].items():
                if not k.startswith('-'):
                    coverage[pos] += v

        return coverage

    def to_consensus(self, percentage):
        """Generates and returns a consensus sequence

        >>> ref = Reference('hxb2_pol',
        >>>     'GAAGAAATGAGTTTGCCAGGAAGATGGAAACCAAAAATGATAGGGGGAATTGGAGGTTTT')
        >>> mrs = MappedReadCollection.from_bam(ref, 65, 75,
        >>>     'tests/data/test1.bam')
        >>> print(mrs.to_consensus(20))
        GAAGAAATGAGTTTGCCAGGAAGATGGAAACCAAAAATGATAGGGGGAATTGGAGGTTTT
        """
        iupac = {'A': 'A',
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
                 'ACGT': 'N'}

        pileup = self.pileup(indels=(percentage == 100))

        # do not include indels in coverage calculations
        coverage = self.coverage(pileup)
        start, end = None, None
        consensus = ''
        for pos in range(0, len(pileup)):
            if coverage[pos] >= 100:
                if start is None:
                    start = pos
                end = pos

            con_bp = ''

            if coverage[pos] > 0:
                # if percentage is 100, we are generating a consensus sequence
                # of the most common sequence
                if percentage == 100:
                    con_bp_with_indel = sorted(
                        [(v, k) for k, v in pileup[pos].items()],
                        reverse=True)[0][1]
                    if len(con_bp_with_indel) > 1:
                        # TODO recalculate consensus base if we can't
                        #      incorporate insertion
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
                    for token, token_count in sorted(pileup[pos].items()):
                        if token != '-' and token.upper() != 'N':
                            if decimal.Decimal(token_count) / coverage[pos] * \
                                    100 >= percentage:
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
            consensus_seq = consensus[start:end + 1]

        return consensus_seq
    def mask_unconfident_differences_from_file(self, variant_file):
        """Prepare variants object from a vcf file"""

        variants = defaultdict(lambda: defaultdict(dict))

        # Read in and parse the variants file
        # The file uses 1 as the first position but 0 is the first position in
        # mapped reads.
        with open(variant_file, "r") as input:
            for line in input:
                if line[0] != "#":
                    chrom, pos, id, ref, alt, qual, filter, info = \
                        line.rstrip().split("\t")

                    variants[int(pos) - 1][alt]["filter"] = filter
        
        self.mask_unconfident_differences(variants)

    def mask_unconfident_differences_from_obj(self, variants_obj):
        """Process variants_obj"""

        variants = defaultdict(lambda: defaultdict(dict))

        for rid in variants_obj.variants:
            pdb.set_trace()
            for pos in variants_obj.variants[rid]:
                for alt_allele, variant in sorted(
                        variants_obj.variants[rid][pos].items()):
                    variants[int(pos) -1][alt_allele]["filter"] = variant.filter

        self.mask_unconfident_differences(variants)

    def mask_unconfident_differences(self, variants):
        """Mask unconfident differences by changing their case to lower"""

        for name, mapped_read in self.mapped_reads.items():
            for pos in mapped_read.differences:
                # mapped_read.differences[pos] will be a string of length 1.
                # or more.
                # If we have a substitution/deletion, it will be of length 1.
                # If we have an insertion, it will be of length >= 2 with the
                # first position being a substitution/deletion or a match (
                # indicated with a '.')
                substitution = mapped_read.differences[pos][:1]
                insertion = mapped_read.differences[pos][1:]

                if substitution != "." and substitution != "-":
                    if (substitution.lower() not in variants[pos] or
                            variants[pos][substitution.lower()]["filter"]
                            != "PASS"):

                        substitution = substitution.lower()

                mapped_read.differences[pos] = substitution + insertion


if __name__ == '__main__':
    import doctest
    doctest.testmod()

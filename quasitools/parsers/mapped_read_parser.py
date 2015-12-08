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
import pysam
from quasitools.utilities import sam_alignment_to_padded_alignment, \
    pairwise_alignment_to_differences
from quasitools.mapped_read import MappedRead, MappedReadCollection


def parse_mapped_reads_from_bam(reference, overlap_cutoff, identity_cutoff,
                                bam):
    """Parse MappedRead mrcects from a bam file and produce a MappedReadCollection.
    """
    mrc = MappedReadCollection(reference)
    sam = pysam.AlignmentFile(bam, "rb")

    for alignment in sam.fetch(reference=reference.name):
        overlap = decimal.Decimal(alignment.query_alignment_length) \
            / alignment.query_length * 100

        if overlap >= overlap_cutoff:
            padded_alignment = sam_alignment_to_padded_alignment(alignment,
                                                                 reference)

            matches = 0
            if alignment.has_tag('NM'):
                matches = alignment.query_length - alignment.get_tag('NM')
            else:
                matches = padded_alignment[1].count('|')

            # TODO not an accurate identity, as I am not distinguishing between
            #     an alignment match vs a sequence match
            identity = decimal.Decimal(matches) / len(padded_alignment[1])\
                * 100

            if identity >= identity_cutoff:
                direct = '+'
                if alignment.flag & 16:
                    direct = '-'
                # TODO only calculate differences when identity < 100
                differences = pairwise_alignment_to_differences(
                    padded_alignment[0], padded_alignment[2],
                    alignment.reference_start)

                mapped_read = MappedRead(alignment.query_name,
                                         alignment.query_alignment_start,
                                         alignment.query_alignment_end-1,
                                         differences,
                                         alignment.reference_start,
                                         alignment.reference_end-1,
                                         overlap, identity, direct)
                read_id = "{0}_{1}".format(alignment.query_name,
                                           '1' if direct == '+' else '2')
                mrc.mapped_reads[read_id] = mapped_read

    return mrc

if __name__ == '__main__':
    import doctest
    doctest.testmod()

"""
Copyright Government of Canada 2015-2017

Written by: Eric Enns, National Microbiology Laboratory,
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

import pysam
from quasitools.utilities import sam_alignment_to_padded_alignment, \
    pairwise_alignment_to_differences
from quasitools.mapped_read import MappedRead, MappedReadCollection
from quasitools.pileup import Pileup, Pileup_List

REVERSE_COMPLEMENTED = 16
FORWARD = '+'
REVERSE = '-'


def parse_mapped_reads_from_bam(reference, bam):
    """Parse MappedRead mrcects from a bam file and produce a MappedReadCollection.
    """
    mrc = MappedReadCollection(reference)
    sam = pysam.AlignmentFile(bam, "rb")

    for alignment in sam.fetch(reference=reference.name):
        padded_alignment = sam_alignment_to_padded_alignment(alignment,
                                                             reference)

        direct = FORWARD
        if alignment.flag & REVERSE_COMPLEMENTED:
            direct = REVERSE
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
                                 direct)
        # generate read_id, such that pair end data does not have the same key
        # in the hash, by adding a 1 for forward read and 2 for reverse read
        read_id = "{0}_{1}".format(alignment.query_name,
                                   '1' if direct == FORWARD else '2')
        mrc.mapped_reads[read_id] = mapped_read

    return mrc


def parse_pileup_from_bam(references, bam_location):

    pileup = []
    samfile = pysam.AlignmentFile(bam_location, "rb")

    for reference in references:

        for column in samfile.pileup(reference=reference.name):

            dictionary = {}

            for read in column.pileups:

                if not read.is_del and not read.is_refskip:
                    # query position is None if is_del or is_refskip is set.

                    base = read.alignment.query_sequence[read.query_position]

                    dictionary[base] = dictionary.get(base, 0) + 1

            pileup.append(dictionary)

    return Pileup(pileup)


def parse_pileup_list_from_bam(references, file_list):
    """
    Create a Pileup_List object.

    INPUT:
        [TUPLE] [references] - references tuple

        [FILE LOCATION TUPLE] [file_list] - files names which represent
                                            a pileup

    RETURN:
        [Pileup_List] - a new object containing a list of Pileup objects.

    POST:
        [None]

    """

    pileups = []

    for bam_location in file_list:

        pileup = parse_pileup_from_bam(references, bam_location)
        pileups.append(pileup)

    return Pileup_List(pileups)


if __name__ == '__main__':
    import doctest
    doctest.testmod()

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
import Bio
from Bio import SeqIO

from quasitools.utilities import sam_alignment_to_padded_alignment, \
    pairwise_alignment_to_differences
from quasitools.mapped_read import MappedRead, MappedReadCollection
from quasitools.pileup import Pileup, Pileup_List

REVERSE_COMPLEMENTED = 16
FORWARD = '+'
REVERSE = '-'
GAP = '-'

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
    """
    PARSE PILEUP FROM BAM


    PURPOSE
    -------

    Constructs a Pileup obect from reference objects and a BAM file.


    INPUT
    -----

    [LIST (REFERENCE)] [references]
        A list of quasitools Reference objects.


    [BAM FILE LOCATION)] [bam_location]
        The file location of the aligned BAM file from which to build the
        pileup object.


    RETURN
    ------

    [Pileup]
        A new pileup object constructed from the information in the Reference
        object(s) and the BAM file.

    """

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
    PARSE PILEUP LIST FROM BAM


    PURPOSE
    -------

    Constructs a Pileup_List object from Reference objects and multiple BAM
    files. The Pileup_List will contain multiple Pileup objects, one
    associated with each BAM file. The Reference objects must be correspond to
    every BAM file.


    INPUT
    -----

    [LIST (REFERENCE)] [references]
        A list of quasitools Reference objects.

    [LIST (BAM FILE LOCATIONS)] [file_list]
        A list of BAM file locations, each corresponding to one alignment
        pileup. All BAM files must each correspond to the same associated
        References objects.


    RETURN
    ------

    [Pileup_List]
        A new Pileup_List object representing a collection of Pileup objects.


    """

    pileups = []

    for bam_location in file_list:

        pileup = parse_pileup_from_bam(references, bam_location)
        pileups.append(pileup)

    return Pileup_List(pileups)


def parse_pileup_from_fasta(reads_location, gaps=False):
    """
    # ========================================================================

    PARSE PILEUP FROM FASTA


    PURPOSE
    -------

    Parses an aligned FASTA file and returns a Pileup file corresponding to
    the aligned FASTA file.


    INPUT
    -----

    [(FASTA) FILE LOCATION] [reads_location]
        The file location of the aligned FASTA file.

    [BOOLEAN] [gaps]
        Whether or not to include gaps in the pileup. This is default by
        false.


    RETURN
    ------

    [Pileup]
        A new pileup object constructed from the information in the aligned
        FASTA file.

    # ========================================================================
    """

    pileup = []
    reads = Bio.SeqIO.parse(reads_location, "fasta")

    read = next(reads)

    for i in range(len(read)):

        pileup.append({})

    while read:

        for i in range(len(read)):

            base = read[i]

            if pileup[i].get(base):
                pileup[i][base] += 1

            else:
                pileup[i][base] = 1

        read = next(reads, None)

    # Remove the gaps from the pileup.
    if not gaps:

        for position in pileup:

            position.pop(GAP, None)

    return Pileup(pileup)


if __name__ == '__main__':
    import doctest
    doctest.testmod()

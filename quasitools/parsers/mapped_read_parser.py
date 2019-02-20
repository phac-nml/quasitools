"""
Copyright Government of Canada 2015-2018

Written by: Eric Enns, National Microbiology Laboratory,
            Public Health Agency of Canada
Updated by: Ahmed Kidwai

Modified by: Matthew Fogel and Eric Marinier, National Microbiology Laboratory,
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
import Bio.SeqIO

from quasitools.utilities import sam_alignment_to_padded_alignment, \
    pairwise_alignment_to_differences
from quasitools.mapped_read import MappedRead, MappedReadCollection
from quasitools.pileup import Pileup, Pileup_List
from quasitools.haplotype import Haplotype, \
        build_consensus_from_haplotypes, sort_haplotypes

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
                                 alignment.query_alignment_end - 1,
                                 differences,
                                 alignment.reference_start,
                                 alignment.reference_end - 1,
                                 direct)
        # generate read_id, such that pair end data does not have the same key
        # in the hash, by adding a 1 for forward read and 2 for reverse read
        read_id = "{0}_{1}".format(alignment.query_name,
                                   '1' if direct == FORWARD else '2')
        mrc.mapped_reads[read_id] = mapped_read

    return mrc


def parse_haplotypes_from_bam(
        samfile,
        reference,
        bam_location,
        start,
        k,
        length):
    """""
    #========================================================================
    PARSE HAPLOTYPES FROM BAM

    PURPOSE
    -------

    Get the haplotypes

    INPUT
    -------
    [LIST (REFERENCE] [references]
        - a list of quasitool reference objects
    [BAM_FILE_LOCATION] [bam_location]
        - The aligned BAM FILE from which we'll
          retrieve our haplotypes.
    [INT] [start]
        - the starting region
    [INT] [k]
        - the number of nucleotides per haplotype.

    RETURN
    -------
    [STRING] [haplotyess]

    COMMENTS)
    -------
    - lets do it.

    #========================================================================
    """

    haplotypes = {}
    # TODO parse reference to get the sequenceID (should be a method already).
    reads = samfile.fetch("AF033819.3", start, start + k)

    for read in reads:

        read_sequence = read.get_forward_sequence()
        haplotype_start = start - read.reference_start
        haplotype_end = haplotype_start + k

        if read.get_overlap(start, start + k) == k:

            sequence = str(read_sequence[haplotype_start: haplotype_end])

            if sequence in haplotypes:
                haplotype = haplotypes.get(sequence)
                haplotype.count += 1
            else:

                haplotypes[sequence] = Haplotype(sequence)

    haplotypes_list = list(haplotypes.values())

    # TODO move the following two function calls  to complexity.
    haplotype_consensus = build_consensus_from_haplotypes(haplotypes_list)
    haplotypes_sorted = sort_haplotypes(haplotypes_list, haplotype_consensus)

    return haplotypes_sorted


def parse_haplotypes_called(
        references,
        reference,
        bam_location,
        k):
    """""
    #========================================================================

    PARSE HAPLOTYPES CALLED

    PURPOSE
    -------

    Caller method for the haplotype parser.

    INPUT
    -------
    [LIST (REFERENCE] [references]
        - a list of quasitool reference objects
    [REFERENCE_LOCATION][reference]
        - The location of the reference file
    [BAM_FILE_LOCATION] [bam_location]
        - The aligned BAM FILE from which we'll
          retrieve our haplotypes.
    [INT] [start]
        - the starting region
    [INT] [k]
        - the number of nucleotides per haplotype.

    RETURN
    -------
    [DICTIONARY of HAPLOTYPES] [haplotyess]

    COMMENTS
    -------
    from the reference in references we can find the number of nucleotides
    which we store as length.

    In the for loop I used 0 and 100 as our range to run the program quick.
    We could customize it to use start and length-k+1 for our range.

    we obtain the sequence by calling the driver method
    parse_haplotypes_from_bam.
    If the length of the sequence is equal to k
    we can create a haplotype object
    (this was to fix a bug where the first haplotype object was always empty)

    create haplotypes, store values as list, sort them and then return list.

    #========================================================================
    """

    haplotypes = []
    samfile = pysam.AlignmentFile(bam_location, "rb")
    for reference in references:
        length = len(reference.seq)
        #  placeholder to pass CI test

        for i in range(0, length - k + 1):

            haplotype_list = (
                parse_haplotypes_from_bam(
                    samfile,
                    reference,
                    bam_location,
                    i, k, length))

            haplotypes.append(haplotype_list)

    return haplotypes


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

    # PySam bases:
    A = 0
    C = 1
    G = 2
    T = 3

    pileup = []
    samfile = pysam.AlignmentFile(bam_location, "rb")

    for reference in references:

        coverage = samfile.count_coverage(
            contig=reference.name, start=0, stop=len(reference.seq),
            quality_threshold=0)

        for column in range(len(coverage[0])):

            dictionary = {}

            if coverage[A][column] > 0:
                dictionary["A"] = coverage[A][column]

            if coverage[C][column] > 0:
                dictionary["C"] = coverage[C][column]

            if coverage[G][column] > 0:
                dictionary["G"] = coverage[G][column]

            if coverage[T][column] > 0:
                dictionary["T"] = coverage[T][column]

            pileup.append(dictionary)

    return Pileup(pileup)


def parse_pileup_list_from_bam(references, file_list):
    """
    # ========================================================================

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

    # ========================================================================
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


def parse_haplotypes_from_fasta(reads_location, consensus):
    """
    # ========================================================================

    PARSE HAPLOTYPES FROM READS


    PURPOSE
    -------

    Builds a list of Haplotype objects from aligned FASTA reads.


    INPUT
    -----

    [FILE LOCATION] [reads_location]
        The location of the aligned FASTA reads.

    [STRING] [consensus]
        The consensus sequence of the pileup.


    RETURN
    ------

    [HAPLOTYPE LIST]
        A list of Haplotype objects, defined by the aligned FASTA reads.

    # ========================================================================
    """

    haplotypes = {}  # (sequence, Haplotype)

    reads = Bio.SeqIO.parse(reads_location, "fasta")

    for read in reads:

        sequence = str(read.seq)
        if sequence in haplotypes:

            haplotype = haplotypes.get(sequence)
            haplotype.count += 1

        else:

            haplotypes[sequence] = Haplotype(sequence, consensus)

    haplotypes_list = list(haplotypes.values())
    haplotypes_sorted = sort_haplotypes(haplotypes_list)

    return haplotypes_sorted


if __name__ == '__main__':
    import doctest
    doctest.testmod()

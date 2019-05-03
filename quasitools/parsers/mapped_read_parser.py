"""
Copyright Government of Canada 2015-2019

Written by: Eric Enns, National Microbiology Laboratory,
            Public Health Agency of Canada

Modified by: Ahmed Kidwai, Matthew Fogel and Eric Marinier,
            National Microbiology Laboratory, Public Health Agency of Canada

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
from quasitools.haplotype import Haplotype


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
        references,
        reference,
        bam_location,
        k):
    """""
    # ========================================================================

    PARSE HAPLOTYPES FROM BAM

    PURPOSE
    -------

    This is the function to call to generate a list of haplotypes of length k,
    across the reference, when the input is a BAM file.

    INPUT
    -------
    [LIST (REFERENCE] [references]
        - A list of quasitool Reference objects.
    [REFERENCE_LOCATION][reference]
        - The location of the reference file.
    [BAM_FILE_LOCATION] [bam_location]
        - The aligned BAM file from which we'll retrieve our haplotypes.
    [INT] [k]
        - The length we want our starting position take reads from.

    RETURN
    -------
    [LIST of HAPLOTYPES] [haplotypes]
        An 2D list containing a list of unsorted haplotypes, for each
        position in the reference sequenence until, reference length - k + 1.

    # ========================================================================
    """

    haplotypes = []
    samfile = pysam.AlignmentFile(bam_location, "rb")

    for reference in references:

        length = len(reference.seq)
        for i in range(0, length - k + 1):
            haplotype_list = (
                parse_haplotypes_from_bam_range(
                    samfile,
                    reference,
                    bam_location,
                    i, k))

            haplotypes.append(haplotype_list)

    return haplotypes


def parse_haplotypes_from_bam_range(
        samfile,
        reference,
        bam_location,
        start,
        k):
    """""
    # ========================================================================
    PARSE HAPLOTYPES FROM BAM

    PURPOSE
    -------

    Builds and returns an unsorted list of Haplotype objects from start to
    start + k.


    INPUT
    -------
    [LIST (REFERENCE] [references]
        - A list of Reference objects associated with the provided BAM
          file location.
    [BAM file location] [bam_location]
        - The aligned BAM FILE from which we'll retrieve our haplotypes.
    [INT] [start]
        - The starting 0-based reference postion.
    [INT] [k]
        - The length of the haplotypes to generate.
    RETURN
    -------
    [LIST] [haplotyess]
        - Unsorted list of Haplotype objects from start to start + k.

    # ========================================================================
    """

    haplotypes = {}

    reads = samfile.fetch(reference.name, start, start + k)

    for read in reads:

        positional_array = create_positional_array(read.cigartuples)

        # Shift the positional array to align with the reference:
        positional_array = [x + read.reference_start for x in positional_array]
        read_sequence = read.query_alignment_sequence

        haplotype_start = get_index_in_list(positional_array, start)
        haplotype_end = get_index_in_list(
                positional_array, start + k - 1) + 1

        # check if read maps to the reference.
        if haplotype_start < 0 or haplotype_end < 0:

            continue

        # Checks for deletions.
        if haplotype_end - haplotype_start != k:

            continue

        # checks for inserts
        if not is_consecutive_list(
                    positional_array[haplotype_start:haplotype_end]):
            continue

        # Checking the read covers the entire region:
        if read.get_overlap(start, start + k) == k:

            sequence = str(read_sequence[haplotype_start: haplotype_end])
            if sequence in haplotypes:
                haplotype = haplotypes.get(sequence)
                haplotype.count += 1
            else:
                haplotypes[sequence] = Haplotype(sequence)

    haplotypes_list = list(haplotypes.values())

    return haplotypes_list


def is_consecutive_list(list_of_integers):
    """
    # ========================================================================

    IS CONSECUTIVE LIST


    PURPOSE
    -------

    Reports if elments in a list increase in a consecutive order.


    INPUT
    -----

    [[List]] [list_of_integers]
        - A list of integers.

    Return
    ------
    [BOOLEAN]
        - Returns true is a list is consecutive or false if the same
          number appears consecutively.

    # ========================================================================
    """

    for i in range(1, len(list_of_integers)):
        if list_of_integers[i] - list_of_integers[i - 1] != 1:
            return False

    return True


def get_index_in_list(list_of_integers, value):
    """
    # ========================================================================

    GET INDEX IN LIST


    PURPOSE
    -------

    Checks if a list contains a specific value.

    INPUT
    -----

    [[List]] [list_of_integers]
       - A list of integers.
    [(INT) VALUE] [value]
       - The value we we are searching for in the list of integers.

    POST
    ------
    [[INDEX] [INT]]:
        - If the value is found in list of integers we return the index,
          otherwise return -1.


    # ========================================================================
    """

    try:
        return list_of_integers.index(value)

    except ValueError:
        return -1


def create_positional_array(cigar_tuples):
    """
    # ========================================================================

    CREATE POSITIONAL ARRAY


    PURPOSE
    -------

    Create a positional array that maps positions in a
    CIGAR tuple to a list.

    Ex. CIGAR tuple is: [(0, 4), (2, 1) (1, 2)]

        Positional Array is Initialized to Empty.
        position (an int) starts at 0.

        We look at each item in the CIGAR tuple where
        the first item is the operation (ex. match, delete, insert)
        and the second item is number of bases involved in the operation.

        The returned array maps positions the read (as a list indicies)
        to relative positions in the reference. This returned list of
        relative positions starts at 0.


        If we have a match we append the current reletive position
        of the reference to the positional array (which represents
        positions in the read) and then we will increase the relative
        position in the reference. This process is repeated for the
        length of the match.

        If the operation is a insertion we appending the positional array
        with the left anchored relative position of the insert in
        the reference. This proccess is repeated for the length of the insert.
        This means the same relative position is appended multiple times.

        If the operation is a deletion we will increase the relative position
        in the reference by the length of the operation.
        This means no value gets appended to the positional array.

        So for the CIGAR tuple list above we would get a positional
        array that looks as follows:

        1. Looking at first tuple in the list:
            The tuple's operation is 0 (i.e a match).
            positional_array = [0, 1, 2, 3]
            position: 4

        2. Looking at second tuple in the list:
            The tuple's operation is 2 (i.e a delete)
            positional_array: [0, 1, 2, 3] (didn't change)
            position: 5

        3. Looking at the third tuple in the list:
            The tuple's operation is 1 (i.e an insert)
            positional_array = [0, 1, 2, 3, 4,4]
            position: 5

    INPUT
    -----

    [[CIGAR] TUPLE] [cigar_tuples]
       - A list containing the CIGAR tuples. (operation, length).

    Return
    ------
    [[LIST] [INT]]
        - A positional array that maps CIGAR tuples to the read.

    # ========================================================================
    """

    positional_array = []
    OPERATION = 0
    LENGTH = 1
    position = 0  # 0-based

    MATCH = 0
    INSERT = 1
    DELETE = 2

    for tup in cigar_tuples:

        if tup[OPERATION] == MATCH:

            for i in range(tup[LENGTH]):

                positional_array.append(position)  # consume read
                position = position + 1  # consume reference

        if tup[OPERATION] == INSERT:

            for i in range(tup[LENGTH]):

                positional_array.append(position - 1)  # consume read

        if tup[OPERATION] == DELETE:

            position += tup[LENGTH]  # consume reference

    return positional_array


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


def parse_haplotypes_from_fasta(reads_location):
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

            haplotypes[sequence] = Haplotype(sequence)

    haplotypes_list = list(haplotypes.values())

    return haplotypes_list


if __name__ == '__main__':
    import doctest
    doctest.testmod()

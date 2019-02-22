"""
# =============================================================================

Copyright Government of Canada 2019

Written by: Eric Marinier, Public Health Agency of Canada,
    National Microbiology Laboratory

Funded by the National Micriobiology Laboratory and the Genome Canada / Alberta
    Innovates Bio Solutions project "Listeria Detection and Surveillance
    using Next Generation Genomics"

Licensed under the Apache License, Version 2.0 (the "License"); you may not use
this file except in compliance with the License. You may obtain a copy of the
License at:

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software distributed
under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR
CONDITIONS OF ANY KIND, either express or implied. See the License for the
specific language governing permissions and limitations under the License.

# =============================================================================
"""


import numpy

import quasitools.calculate as calculate
from quasitools.pileup import Pileup


class Haplotype:
    """
    # ========================================================================

    HAPLOTYPE

    # ========================================================================
    """

    def __init__(self, sequence, count=1):
        """
        # ====================================================================

        INIT

        # ====================================================================
        """

        self.sequence = sequence
        self.count = count

    def __eq__(self, other):
        """
        # ====================================================================

        EQUALS

        # ====================================================================
        """

        # Override the default Equals behavior:

        if isinstance(other, self.__class__):
            return self.sequence == other.sequence

        return False

    def __ne__(self, other):
        """
        # ====================================================================

        NOT EQUALS

        # ====================================================================
        """
        # Override the default Unequal behavior

        if isinstance(other, self.__class__):
            return self.sequence != other.sequence

        return False


def sort_haplotypes(haplotypes, consensus):
    """
    # ========================================================================

    SORT HAPLOTYPES


    PURPOSE
    -------

    Sorts a list of haplotypes according their hamming distance from the
    consensus sequence.


    INPUT
    -----

    [HAPLOTYPE LIST] [haplotypes]
        The list of haplotypes to sort.

    [String] consensus


    RETURN
    ------

    [HAPLOTYPE LIST]
        A list of sorted haplotypes, according to their number of mutations
        from the consensus sequence of all the sequences.

    # ========================================================================
    """

    sorted_haplotypes = []
    hap_list = []
    for haplotype in haplotypes:

        # creates a list of tuples that contains the haplotype and its hamming
        # distance
        hap_list.append(
            (haplotype,
             calculate.hamming_distance(
                 haplotype.sequence,
                 consensus)))

    # Sort list in ascending order based on tuple position 1.
    sorted_list = \
        sorted(hap_list, key=lambda items: items[1], reverse=False)

    # Pull first tuple item (the haplotype) from each item in list
    sorted_haplotypes = [x[0] for x in sorted_list]

    return sorted_haplotypes


def build_consensus_from_haplotypes(haplotypes):

    pileup = build_pileup_from_haplotypes(haplotypes)

    consensus = pileup.build_consensus()

    return consensus


def build_pileup_from_haplotypes(haplotypes):

    pileup_list = []

    if haplotypes:

        length = len(haplotypes[0].sequence)

        for i in range(0, length):
            pileup_list.append({})

        for haplotype in haplotypes:

            for i in range(0, length):

                base = haplotype.sequence[i]

                if pileup_list[i].get(base):
                    pileup_list[i][base] += 1
                else:
                    pileup_list[i][base] = 1

    # No checks for gaps because we shouldn't have any as the reads overlap.

    pileup = Pileup(pileup_list)

    return pileup


def build_distiance_matrix(haplotypes):
    """
    # ========================================================================

    BUILD DISTANCE MATRIX


    PURPOSE
    -------

    Builds a distance matrix of all the haplotypes.


    INPUT
    -----

    [HAPLOTYPE LIST] [haplotypes]
        The list of Haplotypes from which to build the distance matrix.


    RETURN
    ------

    [2D ARRAY]
        A distance matrix of the passed haplotypes.

    # ========================================================================
    """

    x = len(haplotypes)

    matrix = numpy.zeros(shape=(x, x))

    for i in range(0, len(matrix)):

        for j in range(0, len(matrix)):

            matrix[i][j] = calculate_distance(haplotypes[i], haplotypes[j])

    return matrix


def calculate_distance(haplotype1, haplotype2):
    """
    # ========================================================================

    CALCULATE DISTANCE


    PURPOSE
    -------

    Calculates the distance between two haplotypes.


    INPUT
    -----

    [HAPLOTYPE] [haplotype1]
        The first of two haplotypes to calculate the distance between.

    [HAPLOTYPE] [haplotype2]
        The second of two haplotypes to calculate the distance between.


    RETURN
    ------

    [FLOAT]
        The genetic distance between the two passed haplotypes.

    # ========================================================================
    """

    hamming_distance = \
        calculate.hamming_distance(haplotype1.sequence, haplotype2.sequence)
    genetic_distance = \
        float(hamming_distance) / float(len(haplotype1.sequence))

    return genetic_distance


def calculate_total_clones(haplotypes):
    """
    # ========================================================================

    CALCULATE TOTAL CLONES


    PURPOSE
    -------

    Calculates the total number of clones accross multiple haplotypes. Note
    that there may be multiple clones associated with each haplotype.


    INPUT
    -----

    [HAPLOTYPE LIST] [haplotypes]
        The list of Haplotypes from which to build the distance matrix.


    RETURN
    ------

    [INT]
        The total number of clones accross all the haplotypes in the list.

    # ========================================================================
    """

    total = 0

    for haplotype in haplotypes:

        total += haplotype.count

    return total


def build_counts(haplotypes):
    """
    # ========================================================================

    BUILD COUNTS


    PURPOSE
    -------

    Builds the a list of the counts of each haplotype in a passed haplotype
    list.

    Example:

    haplotype1.sequence = "AAA"
    haplotype1.count = 3

    haplotype1.sequence = "CGC"
    haplotype1.count = 5

    haplotype1.sequence = "TCC"
    haplotype1.count = 1

    haplotypes = [haplotype1, haplotype2, haplotype3]

    build_counts(haplotypes) -> [3, 5, 1]


    INPUT
    -----

    [HAPLOTYPE LIST] [haplotypes]
        The list of Haplotypes.


    RETURN
    ------

    [INT LIST]
        A list of the counts of each haplotype, in the same order as the
        original haplotype list.

    # ========================================================================
    """

    # haplotypes = sort_haplotypes(haplotypes)
    counts = []

    for haplotype in haplotypes:

        count = haplotype.count
        counts.append(count)

    return counts


def build_frequencies(haplotypes):
    """
    # ========================================================================

    BUILD FREQUENCIES


    PURPOSE
    -------

    Builds the a list of the frequencies of each haplotype in a passed
    haplotype list.

    Example:

    haplotype1.sequence = "AAA"
    haplotype1.count = 3

    haplotype1.sequence = "CGC"
    haplotype1.count = 5

    haplotype1.sequence = "TCC"
    haplotype1.count = 1

    haplotypes = [haplotype1, haplotype2, haplotype3]

    build_counts(haplotypes) -> [3/9, 5/9, 1/9]


    INPUT
    -----

    [HAPLOTYPE LIST] [haplotypes]
        The list of Haplotypes.


    RETURN
    ------

    [FLOAT LIST]
        A list of the frequencies of each haplotype, in the same order as the
        original haplotype list.

    # ========================================================================
    """

    counts = build_counts(haplotypes)
    total = calculate_total_clones(haplotypes)
    frequencies = []

    for count in counts:

        frequency = float(count) / float(total)
        frequencies.append(frequency)

    return frequencies

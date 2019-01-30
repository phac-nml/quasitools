"""
# =============================================================================

Copyright Government of Canada 2018

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

import math

GAP = '-'


def shannon_entropy(frequencies):
    """
    # ========================================================================

    SHANNON ENTROPY


    PURPOSE
    -------

    Calculates the Shannon entropy of a list of frequencies.


    INPUT
    -----

    [FLOAT LIST] [frequencies]
        The list of frequencies. These frequencies should sum to 1.


    RETURN
    ------

    [FLOAT]
        The Shannon entropy of the frequencies.

    # ========================================================================
    """

    entropy = 0

    for frequency in frequencies:

        entropy -= float(frequency) * math.log(float(frequency))

    return entropy


def minimum_mutation_frequency(M, N, a):
    """
    # ========================================================================

    MINIMUM MUTATION FREQUENCY


    PURPOSE
    -------

    Calculates the minimum mutation frequency.


    INPUT
    -----

    [INT] [M]
        The number of mutations.

    [INT] [N]
        The total number of clones (reads) sampled from the viral
        quasispecies.

    [INT] [a]
        The length of the amplicons.


    RETURN
    ------

    [FLOAT]
        The minimum mutation frequency.

    # ========================================================================
    """

    Mfmin = float(M) / (float(N) * float(a))

    return Mfmin


def mutation_frequency(H, D):
    """
    # ========================================================================

    MUTATION FREQUENCY


    PURPOSE
    -------

    Calculates the mutation frequency.


    INPUT
    -----

    [INT] [H]
        The number of haplotypes.

    [2D ARRAY] [D]
        A distance matrix of haplotypes pair-wise genetic distances
        (fraction of nt differences).


    RETURN
    ------

    [FLOAT]
        The mutation frequency.

    # ========================================================================
    """

    sumd = 0

    for i in range(0, H):

        sumd += D[0][i]

    Mfe = float(sumd) / float(H)

    return Mfe


def maximum_mutation_frequency(H, F, D):
    """
    # ========================================================================

    MAXIMUM MUTATION FREQUENCY


    PURPOSE
    -------

    Calculates the maximum mutation frequency.


    INPUT
    -----

    [INT] [H]
        The number of haplotypes.

    [FLOAT LIST] [F]
        A list of (relative) frequencies.

    [2D ARRAY] [D]
        A distance matrix of haplotypes pair-wise genetic distances
        (fraction of nt differences).


    RETURN
    ------

    [FLOAT]
        The maximum mutation frequency.

    # ========================================================================
    """

    Mfmax = 0

    for i in range(0, H):

        Mfmax += F[i] * D[0][i]

    return Mfmax


def sample_nucleotide_diversity_entity(H, D):
    """
    # ========================================================================

    SAMPLE NUCLEOTIDE DIVERSITY (ENTITY-LEVEL)


    PURPOSE
    -------

    Calculates the sample nucleotide diversity.


    INPUT
    -----

    [INT] [H]
        The number of haplotypes.

    [2D ARRAY] [D]
        A distance matrix of haplotypes pair-wise genetic distances
        (fraction of nt differences).


    RETURN
    ------

    [FLOAT]
        The entity-level sample nucleotide diversity.

    # ========================================================================
    """

    sum_substitutions = 0

    for i in range(0, H):

        for j in range(0, H):

            sum_substitutions += D[i][j]

    diversity = float(sum_substitutions) / (float(H) * float(H - 1))

    return diversity


def population_nucleotide_diversity(H, p, D):
    """
    # ========================================================================

    POPULATION NUCLEOTIDE DIVERSITY


    PURPOSE
    -------

    Calculates the population nucleotide diversity.


    INPUT
    -----

    [INT] [H]
        The number of haplotypes.

    [FLOAT] [p]
        A list of (relative) frequencies.

    [2D ARRAY] [D]
        A distance matrix of haplotypes pair-wise genetic distances
        (fraction of nt differences).


    RETURN
    ------

    [FLOAT]
        The population nucleotide diversity.

    # ========================================================================
    """

    diversity = 0

    for i in range(0, H):

        for j in range(0, H):

            diversity += p[i] * D[i][j] * p[j]

    return diversity


def sample_nucleotide_diversity(N, H, p, D):
    """
    # ========================================================================

    SAMPLE NUCLEOTIDE DIVERSITY


    PURPOSE
    -------

    Calculates the sample nucleotide diversity.


    INPUT
    -----

    [INT] [N]
        The total number of clones (reads) sampled from the viral
        quasispecies.

    [INT] [H]
        The number of haplotypes.

    [FLOAT] [p]
        A list of (relative) frequencies.

    [2D ARRAY] [D]
        A distance matrix of haplotypes pair-wise genetic distances
        (fraction of nt differences).


    RETURN
    ------

    [FLOAT]
        The sample nucleotide diversity.

    # ========================================================================
    """

    diversity = 0

    for i in range(0, H):

        for j in range(0, H):

            diversity += p[i] * D[i][j] * p[j]

    diversity *= (float(N) / float(N - 1))

    return diversity


def simpson_index(H, P):
    """
    # ========================================================================

    SIMPSON INDEX


    PURPOSE
    -------

    Calculates the Simpson index.


    INPUT
    -----

    [INT] [H]
        The number of haplotypes.

    [FLOAT] [P]
        A list of (relative) frequencies.


    RETURN
    ------

    [FLOAT]
        The Simpson index.

    # ========================================================================
    """

    index = 0

    for i in range(0, H):

        index += float(P[i]) * float(P[i])

    return index


def gini_simpson_index(H, P):
    """
    # ========================================================================

    GINI-SIMPSON INDEX


    PURPOSE
    -------

    Calculates the Gini-Simpson index.


    INPUT
    -----

    [INT] [H]
        The number of haplotypes.

    [FLOAT] [P]
        A list of (relative) frequencies.


    RETURN
    ------

    [FLOAT]
        The Gini-Simpson index.

    # ========================================================================
    """

    return (1 - simpson_index(H, P))


def hill_number(H, P, Q):
    """
    # ========================================================================

    HILL NUMBER


    PURPOSE
    -------

    Calculates the simpson index.


    INPUT
    -----

    [INT] [H]
        The number of haplotypes.

    [FLOAT] [P]
        A list of (relative) frequencies.

    [INT] [Q]
        The particular Hill number to calculate.


    RETURN
    ------

    [FLOAT]
        The Hill number for the passed Q.

    # ========================================================================
    """

    # Undefined at one, exponent of entropy:
    if Q == 1:

        entropy = shannon_entropy(P)
        number = math.exp(entropy)
        return number

    else:

        number = 0

        for i in range(0, H):

            number += math.pow(P[i], Q)

        number = math.pow(number, (1.0 / (1.0 - Q)))

        return number


def FAD(H, D):
    """
    # ========================================================================

    FUNCTIONAL ATTRIBUTE DIVERSITY


    PURPOSE
    -------

    Calculates the functional attribute diversity.


    INPUT
    -----

    [INT] [H]
        The number of haplotypes.

     [2D ARRAY] [D]
        A distance matrix of haplotypes pair-wise genetic distances
        (fraction of nt differences).


    RETURN
    ------

    [FLOAT]
        The functional attribute diversity.

    # ========================================================================
    """

    number = 0

    for i in range(0, H):

        for j in range(0, H):

            number += D[i][j]

    return number


def hamming_distance(sequence1, sequence2):
    """
    # ========================================================================

    HAMMING DISTANCE


    PURPOSE
    -------

    Calculates the Hamming distance between two sequences.


    INPUT
    -----

    [STRING] [sequence1]
        The first of two sequences.

    [STRING] [sequence2]
        The second of two sequences.


    RETURN
    ------

    [INT]
        The Hamming distance between the two passed sequences.

    # ========================================================================
    """

    if len(sequence1) != len(sequence2):
        raise ValueError(
            "Hamming Distance is undefined for sequences of unequal length.")

    distance = 0

    for i in range(0, len(sequence1)):

        if (sequence1[i] != sequence2[i]) \
                and (sequence1[i] != GAP) and (sequence2[i] != GAP):

            distance += 1

    return distance

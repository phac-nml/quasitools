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
import numpy

GAP = '-'

"""
# =============================================================================
# =============================================================================
"""
def shannon_entropy(frequencies):

    entropy = 0

    for frequency in frequencies:

        entropy -= float(frequency) * math.log(float(frequency))

    return entropy

"""
# =============================================================================
# =============================================================================
"""
def minimum_mutation_frequency(M, N, a):

    Mfmin = float(M) / (float(N) * float(a))

    return Mfmin

"""
# =============================================================================
# =============================================================================
"""
def mutation_frequency(H, D):

    sumd = 0

    for i in range(0, H):

        sumd += D[0][i]

    Mfe = float(sumd) / float(H)

    return Mfe

"""
# =============================================================================
# =============================================================================
"""
def maximum_mutation_frequency(H, F, D):

    Mfmax = 0

    for i in range(0, H):

        Mfmax += F[i] * D[0][i]

    return Mfmax

"""
# =============================================================================
# =============================================================================
"""
def sample_nucleotide_diversity_entity(H, D):

    sum_substitutions = 0

    for i in range(0, H):

        for j in range(0, H):

            sum_substitutions += D[i][j]

    diversity = float(sum_substitutions) / (float(H) * float(H - 1))

    return diversity

"""
# =============================================================================
# =============================================================================
"""
def population_nucleotide_diversity(H, p, D):

    diversity = 0

    for i in range(0, H):

        for j in range(0, H):

            diversity += p[i] * D[i][j] * p[j]

    return diversity


"""
# =============================================================================
# =============================================================================
"""
def sample_nucleotide_diversity(N, H, p, D):

    diversity = 0

    for i in range(0, H):

        for j in range(0, H):

            diversity += p[i] * D[i][j] * p[j]

    diversity *= (float(N) / float(N - 1))

    return diversity


"""
# =============================================================================
# =============================================================================
"""
def simpson_index(H, P):

    index = 0

    for i in range(0, H):

        index += float(P[i]) * float(P[i])

    return index

"""
# =============================================================================
# =============================================================================
"""
def gini_simpson_index(H, P):

    return (1 - simpson_index(H, P))

"""
# =============================================================================
# =============================================================================
"""
def hill_number(H, P, Q):

    # Undefined at one, exponent of entropy:
    if Q == 1:

        entropy = shannon_entropy(P)
        number = math.exp(entropy)
        return number

    else:

        number = 0

        for i in range(0, H):

            number += math.pow(P[i], Q)

        number = math.pow(number, (1 / (1 - Q) ))

        return number

"""
# =============================================================================
# =============================================================================
"""
def FAD(H, D):

    number = 0

    for i in range(0, H):

        for j in range(0, H):

            number += D[i][j]

    return number

"""
# =============================================================================
# =============================================================================
"""
def normalize(frequencies):

    normalized = []

    total = sum(frequencies)

    for frequency in frequencies:

        normalized.append(float(frequency) / float(total))

    return normalized

"""
# =============================================================================
# =============================================================================
"""
def hamming_distance(sequence1, sequence2):

    if len(sequence1) != len(sequence2):
        raise ValueError("Hamming Distance is undefined for sequences of unequal length.")

    distance = 0

    for i in range(0, len(sequence1)):

        if (sequence1[i] != sequence2[i]) and (sequence1[i] != GAP) and (sequence2[i] != GAP):

            distance += 1

    return distance




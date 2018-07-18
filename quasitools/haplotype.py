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

from operator import attrgetter

import Bio
from Bio import SeqIO
import numpy

import quasitools.calculate

class Haplotype:

    def __init__(self, sequence, consensus, count=1):

        self.sequence = sequence
        self.count = count
        self.mutations = calculate.hamming_distance(sequence, consensus)

    def __eq__(self, other):
        # Override the default Equals behavior:

        if isinstance(other, self.__class__):
            return self.sequence == other.sequence

        return False
 
    def __ne__(self, other):
        # Override the default Unequal behavior

        if isinstance(other, self.__class__):
            return self.sequence != other.sequence

        return False

"""
# =============================================================================
# =============================================================================
"""
def build_from_reads(reads_location, consensus):

    haplotypes = {} # (sequence, Haplotype)

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

"""
# =============================================================================
# =============================================================================
"""
def sort_haplotypes(haplotypes):

    sorted_haplotypes = sorted(haplotypes, key=attrgetter('mutations'), reverse=False)

    return sorted_haplotypes

"""
# =============================================================================
# =============================================================================
"""
def build_distiance_matrix(haplotypes):

    haplotypes = sort_haplotypes(haplotypes)
    x = len(haplotypes)

    matrix = numpy.zeros(shape=(x, x))

    for i in range(0, len(matrix)):

        for j in range(0, len(matrix)):

            matrix[i][j] = calculate_distance(haplotypes[i], haplotypes[j])

    return matrix
"""
# =============================================================================
# =============================================================================
"""
def calculate_distance(haplotype1, haplotype2):

    hamming_distance = calculate.hamming_distance(haplotype1.sequence, haplotype2.sequence)
    genetic_distance = hamming_distance / len(haplotype1.sequence)

    return genetic_distance

"""
# =============================================================================
# =============================================================================
"""
def calculate_total_clones(haplotypes):

    total = 0

    for haplotype in haplotypes:

        total += haplotype.count

    return total

"""
# =============================================================================
# =============================================================================
"""
def build_counts(haplotypes):

    haplotypes = sort_haplotypes(haplotypes)
    counts = []

    for haplotype in haplotypes:

        count = haplotype.count
        counts.append(count)

    return counts

"""
# =============================================================================
# =============================================================================
"""
def build_frequencies(haplotypes):

    counts = build_counts(haplotypes)
    total = calculate_total_clones(haplotypes)
    frequencies = []

    for count in counts:

        frequency = float(count) / float(total)
        frequencies.append(frequency)

    return frequencies



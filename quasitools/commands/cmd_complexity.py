"""
Copyright Government of Canada 2019

Written by: Eric Marinier and Ahmed Kidwai, Public Health Agency of Canada

Licensed under the Apache License, Version 2.0 (the "License"); you may not use
this work except in compliance with the License. You may obtain a copy of the
License at:

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software distributed
under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR
CONDITIONS OF ANY KIND, either express or implied. See the License for the
specific language governing permissions and limitations under the License.
"""

__version__ = '0.1.1'

import os
import click
import math
import copy
import csv
import quasitools.calculate as calculate
import quasitools.haplotype as haplotype
from quasitools.parsers.reference_parser import parse_references_from_fasta
from quasitools.parsers.mapped_read_parser import parse_haplotypes_called

UNDEFINED = ""

# Measurement Positions
NUMBER_OF_HAPLOTYPES = 0
NUMBER_OF_POLYMORPHIC_SITES = 1
NUMBER_OF_MUTATIONS = 2
SHANNON_ENTROPY_NUMBER = 3
SHANNON_ENTROPY_NUMBER_LOCALIZED_TO_N = 4
SHANNON_ENTROPY_NUMBER_LOCALIZED_TO_H = 5
SIMPSON_INDEX = 6
GINI_SIMPSON_INDEX = 7
HILL_NUMBER_0 = 8
HILL_NUMBER_1 = 9
HILL_NUMBER_2 = 10
HILL_NUMBER_3 = 11
MINIMUM_MUTATION_FREQUENCY = 12
MUTATION_FREQUENCY_FREQUENCY = 13
FUNCTIONAL_ATTRIBUTE_DIVERSITY = 14
SAMPLE_NUCLEOTIDE_DIVERSITY_Entity = 15
MAXIMUM_MUTATION_FREQUENCY = 16
POPULATION_NUCLEOTIDE_DIVERSITY = 17
SAMPLE_NUCLEOTIDE_DIVERSITY = 18


# Dictionary of Names
MEASUREMENTS_NAMES = {
    NUMBER_OF_HAPLOTYPES: "Number of Haplotypes",
    NUMBER_OF_POLYMORPHIC_SITES: "Number of Polymorphic Sites",
    NUMBER_OF_MUTATIONS: "Number of Mutations",
    SHANNON_ENTROPY_NUMBER: "Shannon Entropy",
    SHANNON_ENTROPY_NUMBER_LOCALIZED_TO_N: "Shannon Entropy Localized to N",
    SHANNON_ENTROPY_NUMBER_LOCALIZED_TO_H: "Shannon Entropy Localized to H",
    SIMPSON_INDEX: 'Simpson Index',
    GINI_SIMPSON_INDEX: "Gini Simpson Index",
    HILL_NUMBER_0: "Hill Number #0",
    HILL_NUMBER_1: "HIll Number #1",
    HILL_NUMBER_2: "Hill Number #2",
    HILL_NUMBER_3: "Hill Number #3",
    MINIMUM_MUTATION_FREQUENCY: "Minimum Mutation Frequency",
    MUTATION_FREQUENCY_FREQUENCY: "Mutation Frequency",
    FUNCTIONAL_ATTRIBUTE_DIVERSITY: "Functional Attribute Diversity",
    SAMPLE_NUCLEOTIDE_DIVERSITY_Entity: "Sample Nucleotide Diversity Entity",
    MAXIMUM_MUTATION_FREQUENCY: "Maximum Mutation Frequency",
    POPULATION_NUCLEOTIDE_DIVERSITY: "Population Nucleotide Diversity",
    SAMPLE_NUCLEOTIDE_DIVERSITY: "Sample Nucleotide Diversity",
}

HILL_NUMBER = 4


@click.command(
    'complexity', short_help='Calculates various quasispecies complexity \
    measures.')
@click.argument('reference', nargs=1, required=True,
                type=click.Path(exists=True, file_okay=True, dir_okay=False))
@click.argument('bam', nargs=1,
                type=click.Path(exists=True, file_okay=True, dir_okay=False))
@click.argument('k')
def cli(reference, bam, k):
    """

    Reports the complexity of a quasispecies sequenced through next
    generation sequencing using several measures
    outlined in the following work:

    Gregori, Josep, et al. "Viral quasispecies complexity measures."
    Virology 493 (2016): 227-237.

    """

    click.echo("Using file %s as reference" % reference)
    click.echo("Reading input from file(s)  %s" % bam)
    click.echo("Checking for overlaps at every position in reference" +
               " at %s intervals" % k)

    complexity(reference, bam, int(k))


def complexity(reference, bam, k):
    """
    # ========================================================================

    COMPLEXITY


    PURPOSE
    -------

    Computes various complexity measures for the quasispecies.


    INPUT
    -----

    [(BAM) FILE LOCATION] [bam]
        The file location of a bam file

    [(REFERENCE) FILE LOCATION] [refernece]
        the file location of the reference file
    [INT] k
        provides the sequence length for our reads from a given starting
        position


    RETURN
    ------

    [NONE]

    Comments
    ----

    The complexity computation and reporting will be completed once
    this method has run its course and be stored in CSV file.

    # ========================================================================
    """

    references = parse_references_from_fasta(reference)
    haplotype_list = parse_haplotypes_called(
        references, reference, bam, k)

    measurements_list = []

    for i in range((len(haplotype_list) - k + 1)):
        haplotypes = haplotype_list[i]

        measurements = [UNDEFINED for x in range(len(MEASUREMENTS_NAMES))]
        measurements_list.append(measurements)

        if not haplotypes:
            continue

        haplotype_consensus = haplotype.build_consensus_from_haplotypes(
            haplotypes)
        sorted_haplotypes = haplotype.sort_haplotypes(
            haplotypes, haplotype_consensus)

        pileup = haplotype.build_pileup_from_haplotypes(sorted_haplotypes)

        distance_matrix = haplotype.build_distiance_matrix(sorted_haplotypes)
        counts = haplotype.build_counts(sorted_haplotypes)
        frequencies = haplotype.build_frequencies(sorted_haplotypes)

        '''
        Set the Incidence - Entity Level
        '''
        measurements[NUMBER_OF_HAPLOTYPES] = \
            get_number_of_haplotypes(sorted_haplotypes)

        measurements[NUMBER_OF_POLYMORPHIC_SITES] = \
            get_number_of_polymorphic_sites(pileup)

        measurements[NUMBER_OF_MUTATIONS] = get_number_of_mutations(pileup)

        '''
        Set the Abundance - Molecular Level
        '''

        shannon_entropy = get_shannon_entropy(sorted_haplotypes, frequencies)

        measurements[SHANNON_ENTROPY_NUMBER] = shannon_entropy

        measurements[SHANNON_ENTROPY_NUMBER_LOCALIZED_TO_N] = \
            get_shannon_entropy_localized_to_n(
            sorted_haplotypes, shannon_entropy)

        measurements[SHANNON_ENTROPY_NUMBER_LOCALIZED_TO_H] = \
            get_shannon_entropy_localized_to_h(
            sorted_haplotypes, shannon_entropy)

        measurements[SIMPSON_INDEX] = \
            get_simpson_index(frequencies)

        measurements[GINI_SIMPSON_INDEX] = \
            get_gini_simpson_index(frequencies)

        hill_numbers_list = get_hill_numbers(frequencies)

        '''
        loops through all positions in the hills_numbers_list then adds
        the subsequent position of hill number to the index of the
        measurements list in charge of holding the specific hill number.

        '''
        for k in range(len(hill_numbers_list)):
            measurements[HILL_NUMBER_0 + k] = (hill_numbers_list[k])

        '''
        Functional,  Indidence - Entity Level
        '''
        measurements[MINIMUM_MUTATION_FREQUENCY] = \
            get_minimum_mutation_frequency(
            sorted_haplotypes, pileup)

        measurements[MINIMUM_MUTATION_FREQUENCY] = get_mutation_frequency(
            distance_matrix)

        measurements[FUNCTIONAL_ATTRIBUTE_DIVERSITY] = \
            get_FAD(distance_matrix)

        measurements[SAMPLE_NUCLEOTIDE_DIVERSITY_Entity] = \
            get_sample_nucleotide_diversity_entity(
            distance_matrix, frequencies)
        '''

        Functional, Abundance - Molecular Level
        '''
        measurements[MAXIMUM_MUTATION_FREQUENCY] = \
            get_maximum_mutation_frequency(
            counts, distance_matrix, frequencies)

        measurements[POPULATION_NUCLEOTIDE_DIVERSITY] = \
            get_population_nucleotide_diversity(
            distance_matrix, frequencies)
        '''

        Other
        '''
        measurements[SAMPLE_NUCLEOTIDE_DIVERSITY] = \
            get_sample_nucleotide_diversity(
            distance_matrix, frequencies, sorted_haplotypes)

    '''
    Measurement to CSV
    '''
    measurement_to_csv(measurements_list)


def get_sample_nucleotide_diversity(
        distance_matrix,
        frequencies,
        haplotypes):
    """
    # ========================================================================

    GETSAMPLE NUCLEOTIDE DIVERSITY


    PURPOSE
    -------

    Reports the sample nucleotide diversity.


    INPUT
    -----

    [HAPLOTYPE LIST] [haplotypes]
        A list of Haplotype objects.

    [FLOAT LIST] [frequencies]
        A list of (relative) frequencies of the Haplotypes.

    [2D ARRAY] [distance_matrix]
        A two dimensional array, representing the distance matrix of distances
        between the haplotypes.

        This is expected to be calculated in a similar manner as:
            haplotype.build_distiance_matrix(haplotypes)


    RETURN
    ------

    [INT] [snd]
        the sample nucleotide diversity.

    # ========================================================================
    """

    N = haplotype.calculate_total_clones(haplotypes)
    H = len(frequencies)
    P = frequencies
    D = distance_matrix

    if N > 1:
        snd = calculate.sample_nucleotide_diversity(N, H, P, D)
    else:
        snd = UNDEFINED

    return snd


def get_population_nucleotide_diversity(
        distance_matrix, frequencies):
    """
    # ========================================================================
    GET POPULATION NUCLEOTIDE DIVERSITY

    PURPOSE
    -------
    Returns the population nucleotide diversity.

    INPUT
    -----
    [2D ARRAY] [distance_matrix]
        A two dimensional array, representing the distance matrix of distances
        between the haplotypes.
    [FLOAT LIST] [frequencies]
        A list of (relative) frequencies of the Haplotypes.


    RETURN
    ------

    [INT] [pnd]
        The population nucleotide diversity

    # ========================================================================
    """

    H = len(frequencies)
    P = frequencies
    D = distance_matrix

    pnd = calculate.population_nucleotide_diversity(H, P, D)

    return pnd


def get_maximum_mutation_frequency(
        counts,
        distance_matrix,
        frequencies):
    """
    # ========================================================================

    GET MAXMIMUM MUTATION FREQUENCY


    PURPOSE
    -------

    Returns the maximum mutation frequency of the haplotypes.


    INPUT
    -----

    [INT LIST] [counts]
        A sorted list of haplotype counts, from the counts of the most
        abundant to the counts of the least abundant haplotype.

    [FLOAT LIST] [frequencies]
        A list of (relative) frequencies of the Haplotypes.

    [2D ARRAY] [distance_matrix]
        A two dimensional array, representing the distance matrix of distances
        between the haplotypes.

        This is expected to be calculated in a similar manner as:
            haplotype.build_distiance_matrix(haplotypes)

    RETURN
    ------

    [INT] [maximum_mutation_frequency]
        The maximum mutation frequency

    # ========================================================================
    """

    H = len(counts)
    F = frequencies
    D = distance_matrix

    maximum_mutation_frequency = calculate.maximum_mutation_frequency(H, F, D)

    return maximum_mutation_frequency


def get_sample_nucleotide_diversity_entity(
        distance_matrix, frequencies):
    """
    # ========================================================================

    GET SAMPLE NUCLEOTIDE DIVERSITY (ENTITY LEVEL)


    PURPOSE
    -------

    Returns the entity-level sample nucleotide diversity.


    INPUT
    -----

    [FLOAT LIST] [frequencies]
        A list of (relative) frequencies of the Haplotypes.

    [2D ARRAY] [distance_matrix]
        A two dimensional array, representing the distance matrix of distances
        between the haplotypes.

        This is expected to be calculated in a similar manner as:
            haplotype.build_distiance_matrix(haplotypes)


    RETURN
    ------

    [INT] snde
        sample nucleotide diversity entity

    # ========================================================================
    """

    H = len(frequencies)
    D = distance_matrix

    if H > 1:
        snde = calculate.sample_nucleotide_diversity_entity(H, D)
    else:
        snde = UNDEFINED

    return snde


def get_FAD(distance_matrix):
    """
    # ========================================================================

    REPORT FUNCTIONAL ATTRIBUTE DIVERSITY


    PURPOSE
    -------

    Reports the functional attribute diversity.


    INPUT
    -----

    [2D ARRAY] [distance_matrix]
        A two dimensional array, representing the distance matrix of distances
        between the haplotypes.

        This is expected to be calculated in a similar manner as:
            haplotype.build_distiance_matrix(haplotypes)


    RETURN
    ------

    [INT] [fad]
        the functional attribute diversity

    # ========================================================================
    """

    H = len(distance_matrix)
    D = distance_matrix

    fad = calculate.FAD(H, D)

    return fad


def get_mutation_frequency(distance_matrix):
    """
    # ========================================================================

    GET MUTATION FREQUENCY


    PURPOSE
    -------

    Reports the mutation frequency of the haplotypes.


    INPUT
    -----

    [2D ARRAY] [distance_matrix]
        A two dimensional array, representing the distance matrix of distances
        between the haplotypes.

        This is expected to be calculated in a similar manner as:
            haplotype.build_distiance_matrix(haplotypes)


    RETURN
    ------

    [INT] [mutation_frequency]
        the mutation frequency

    # ========================================================================
    """

    H = len(distance_matrix)
    D = distance_matrix
    mutation_frequency = calculate.mutation_frequency(H, D)

    return mutation_frequency


def get_minimum_mutation_frequency(haplotypes, pileup):
    """
    # ========================================================================

    Get MINIMUM MUTATION FREQUENCY


    PURPOSE
    -------

    Returns the minimum mutation frequency of the haplotypes.


    INPUT
    -----

    [PILEUP] [pileup]
        A Pileup object, which represents the pileup of aligned reads.

    [HAPLOTYPE LIST] [haplotypes]
        A list of Haplotype objects.


    RETURN
    ------

    [INT][minimum_mutation_frequency]
        the minimum mutation frequency

    # ========================================================================
    """

    M = pileup.count_unique_mutations()
    N = haplotype.calculate_total_clones(haplotypes)
    a = len(haplotypes[0].sequence)

    minimum_mutation_frequency = calculate.minimum_mutation_frequency(M, N, a)

    return minimum_mutation_frequency


def get_number_of_haplotypes(haplotypes):
    """
    # ========================================================================

    GET NUMBER OF HAPLOTYPES


    PURPOSE
    -------

    Returns the number of (unique) haplotypes.


    INPUT
    -----

    [HAPLOTYPE LIST] [haplotypes]
        A list of Haplotype objects.


    RETURN
    ------

    [INT] [len(haplotypes)]
        unique number of haplotypes

    # ========================================================================
    """

    return len(haplotypes)


def get_number_of_polymorphic_sites(pileup):
    """
    # ========================================================================

    GET NUMBER OF POLYMORPHIC SITES


    PURPOSE
    -------

    Returns the number of polymorphic sites.


    INPUT
    -----

    [PILEUP] [pileup]
        A Pileup object, which represents the pileup of aligned reads.


    RETURN
    ------

    [INT] (pileup.count_polymorphic_sites())
        number of polymorphic sites in pileup

    # ========================================================================
    """

    return pileup.count_polymorphic_sites()


def get_number_of_mutations(pileup):
    """
    # ========================================================================

    GET NUMBER OF MUTATIONS


    PURPOSE
    -------

    Returns the number of mutations.


    INPUT
    -----

    [PILEUP] [pileup]
        A Pileup object, which represents the pileup of aligned reads.


    RETURN
    ------

    [INT] pileup.count_unique_mutations()
        the number of mutations

    # ========================================================================
    """

    return pileup.count_unique_mutations()


def get_shannon_entropy(haplotypes, frequencies):
    """
    # ========================================================================

    GET SHANNON ENTROPY


    PURPOSE
    -------

    Returns the Shannon entropy of the haplotypes.


    INPUT
    -----

    [HAPLOTYPE LIST] [haplotypes]
        A list of Haplotype objects, for which to report the Shannon entropy.

    [FLOAT LIST] [frequencies]
        A list of (relative) frequencies of the Haplotypes.


    RETURN
    ------

    [INT] [hs]
        shannon entropy


    POST
    ----

    The Shannon entropy

    # ========================================================================
    """

    Hs = calculate.shannon_entropy(frequencies)

    return Hs


def get_shannon_entropy_localized_to_n(haplotypes, Hs):
    """
    # ========================================================================

    GET SHANNON ENTROPY LOCALIZED TO N


    PURPOSE
    -------

    Returns the Shannon entropy of the haplotypes localized to N.


    INPUT
    -----

    [HAPLOTYPE LIST] [haplotypes]
        A list of Haplotype objects, for which to report the Shannon entropy.

    [FLOAT LIST] [frequencies]
        A list of (relative) frequencies of the Haplotypes.


    RETURN
    ------

    [INT] [hsn]
        shannon entropy localized to n

    # ========================================================================
    """

    N = haplotype.calculate_total_clones(haplotypes)

    if float(math.log(N)) != 0:
        Hsn = float(Hs) / float(math.log(N))  # entropy localized to N
    else:
        Hsn = 0

    return Hsn


def get_shannon_entropy_localized_to_h(haplotypes, Hs):
    """
     # ========================================================================

     GET SHANNON ENTROPY LOCALIZED TO H


     PURPOSE
     -------

     Returns the Shannon entropy of the haplotypes localized to N.


     INPUT
     -----

     [HAPLOTYPE LIST] [haplotypes]
         A list of Haplotype objects, for which to report the Shannon entropy.

     [FLOAT LIST] [frequencies]
         A list of (relative) frequencies of the Haplotypes.


     RETURN
     ------

     [INT] [hsn]
         shannon entropy localized to h

     # ========================================================================
     """

    H = len(haplotypes)

    if float(math.log(H)) != 0:
        Hsh = float(Hs) / float(math.log(H))  # entropy localized to H
    else:
        Hsh = 0

    return Hsh


def get_simpson_index(frequencies):
    """
    # ========================================================================

    GET SIMPSON INDEX


    PURPOSE
    -------

    Returns the simpson index.


    INPUT
    -----

    [FLOAT LIST] [frequencies]
        A list of (relative) frequencies of the Haplotypes.


    RETURN
    ------

    [INT] [simpson_index]

    # ========================================================================
    """

    H = len(frequencies)
    P = frequencies

    simpson_index = calculate.simpson_index(H, P)

    return simpson_index


def get_gini_simpson_index(frequencies):
    """
    # ========================================================================

    GET GINI-SIMPSON INDEX


    PURPOSE
    -------

    Returns the Gini-Simpson index.


    INPUT
    -----

    [FLOAT LIST] [frequencies]
        A list of (relative) frequencies of the Haplotypes.


    RETURN
    ------

    [INT] [gini_simpson_index]

    # ========================================================================
    """

    H = len(frequencies)
    P = frequencies
    gini_simpson_index = calculate.gini_simpson_index(H, P)

    return gini_simpson_index


def get_hill_numbers(frequencies):
    """
    # ========================================================================

    GET HILL NUMBERS


    PURPOSE
    -------

    RETURN the 0, 1, 2, and 3 Hill numbers.


    INPUT
    -----

    [FLOAT LIST] [frequencies]
        A list of (relative) frequencies of the Haplotypes.


    RETURN
    ------

    [FLOAT LIST] [list_of_hill_numbers]

    # ========================================================================
    """

    P = frequencies
    H = len(frequencies)

    NUMBER_OF_HILL_NUMBERS = 4

    list_of_hill_numbers = []

    for hill_number_pos in range(NUMBER_OF_HILL_NUMBERS):
        if len(P) >= hill_number_pos + 1:
            list_of_hill_numbers.append(calculate.hill_number(
                H, P, hill_number_pos))

    return list_of_hill_numbers


def measurement_to_csv(measurements_list):

    measurements_col_titles = ["Position"]

    for i in range(len(MEASUREMENTS_NAMES)):
        measurements_col_titles.append(MEASUREMENTS_NAMES[i])

    file_name = "complexity_outputs.csv"

    for position in range(len(measurements_list)):

        measurements = copy.deepcopy(measurements_list[position])

        with open(file_name, 'a') as complexity_data:

            writer = csv.writer(complexity_data)
            # if file empty add the column titles (will be first row)
            if os.stat(file_name).st_size == 0:
                writer.writerow(measurements_col_titles)
                # will always add the measurements values
            writer.writerow([position] + measurements)
        complexity_data.close()

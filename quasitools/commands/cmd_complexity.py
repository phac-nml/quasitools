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
import quasitools.constants.con_complexity as constant
from quasitools.parsers.reference_parser import parse_references_from_fasta
from quasitools.parsers.mapped_read_parser import \
    parse_haplotypes_called,\
    parse_haplotypes_from_fasta


@click.group(invoke_without_command=False)
@click.pass_context
def cli(ctx):
    pass


# Multiple Aligned FASTA.
@cli.command(
    'fasta', short_help='Calculates various quasispecies complexity ' +
    'measures on a multiple aligned FASTA file.')
@click.argument('fasta_location', nargs=1,
                type=click.Path(exists=True, file_okay=True, dir_okay=False))
def fasta(fasta_location):

    haplotypes = parse_haplotypes_from_fasta(fasta_location)
    measurements_list = []

    measurements = conduct_measurements(haplotypes)

    measurements_list.append(measurements)

    measurement_to_csv(measurements_list)


# NGS Data from BAM and its corresponding reference file
@cli.command(
    'bam', short_help="Calculates various quasispecies complexity " +
    "measures on next generation sequenced data from a BAM file " +
    "and it's corresponding reference file.")
@click.argument('reference', nargs=1, required=True,
                type=click.Path(exists=True, file_okay=True, dir_okay=False))
@click.argument('bam', nargs=1,
                type=click.Path(exists=True, file_okay=True, dir_okay=False))
@click.argument('k')
def bam(reference, bam, k):
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
    k = int(k)
    references = parse_references_from_fasta(reference)
    haplotype_list = parse_haplotypes_called(
        references, reference, bam, k)

    measurements_list = []

    for i in range((len(haplotype_list) - k + 1)):
        haplotypes = haplotype_list[i]
        measurements = conduct_measurements(haplotypes)
        measurements_list.append(measurements)

    '''
    Measurement to CSV
    '''
    measurement_to_csv(measurements_list)


def conduct_measurements(haplotypes):

    measurements = [constant.UNDEFINED for x in range(
        len(constant.MEASUREMENTS_NAMES))]

    if not haplotypes:
        return measurements

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
    measurements[constant.NUMBER_OF_HAPLOTYPES] = \
        get_number_of_haplotypes(sorted_haplotypes)

    measurements[constant.NUMBER_OF_POLYMORPHIC_SITES] = \
        get_number_of_polymorphic_sites(pileup)

    measurements[constant.NUMBER_OF_MUTATIONS] = get_number_of_mutations(
        pileup)

    '''
    Set the Abundance - Molecular Level:
    '''

    shannon_entropy = get_shannon_entropy(sorted_haplotypes, frequencies)

    measurements[constant.SHANNON_ENTROPY_NUMBER] = shannon_entropy

    measurements[constant.SHANNON_ENTROPY_NUMBER_LOCALIZED_TO_N] = \
        get_shannon_entropy_localized_to_n(
        sorted_haplotypes, shannon_entropy)

    measurements[constant.SHANNON_ENTROPY_NUMBER_LOCALIZED_TO_H] = \
        get_shannon_entropy_localized_to_h(
        sorted_haplotypes, shannon_entropy)

    measurements[constant.SIMPSON_INDEX] = \
        get_simpson_index(frequencies)

    measurements[constant.GINI_SIMPSON_INDEX] = \
        get_gini_simpson_index(frequencies)

    hill_numbers_list = get_hill_numbers(
        frequencies, constant.HILL_NUMBER_LENGTH)

    '''
    loops through all positions in the hills_numbers_list then adds
    the subsequent position of hill number to the index of the
    measurements list in charge of holding the specific hill number.

    '''
    for k in range(len(hill_numbers_list)):
        measurements[constant.HILL_NUMBER_0 + k] = (hill_numbers_list[k])

    '''
    Functional,  Indidence - Entity Level
    '''
    measurements[constant.MINIMUM_MUTATION_FREQUENCY] = \
        get_minimum_mutation_frequency(
        sorted_haplotypes, pileup)

    measurements[constant.MUTATION_FREQUENCY] = get_mutation_frequency(
        distance_matrix)

    measurements[constant.FUNCTIONAL_ATTRIBUTE_DIVERSITY] = \
        get_FAD(distance_matrix)

    measurements[constant.SAMPLE_NUCLEOTIDE_DIVERSITY_Entity] = \
        get_sample_nucleotide_diversity_entity(
        distance_matrix, frequencies)
    '''

    Functional, Abundance - Molecular Level
    '''
    measurements[constant.MAXIMUM_MUTATION_FREQUENCY] = \
        get_maximum_mutation_frequency(
        counts, distance_matrix, frequencies)

    measurements[constant.POPULATION_NUCLEOTIDE_DIVERSITY] = \
        get_population_nucleotide_diversity(
        distance_matrix, frequencies)

    '''
    Other
    '''
    measurements[constant.SAMPLE_NUCLEOTIDE_DIVERSITY] = \
        get_sample_nucleotide_diversity(
        distance_matrix, frequencies, sorted_haplotypes)

    return measurements


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

    [FLOAT] [snd]
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
        snd = constant.UNDEFINED

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

    [FLOAT] [pnd]
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
        A haplotype counts, from the counts of the most
        abundant to the counts of the least abundant haplotype.

    [FLOAT LIST] [frequencies]
        A list of (relative) frequencies of the Haplotypes.

    [2D ARRAY] [distance_matrix]
        A two dimensional array, representing the distance matrix of distances
        between the sorted haplotypes.

        This is expected to be calculated in a similar manner as:
            haplotype.build_distiance_matrix(haplotypes)

    RETURN
    ------

    [FLOAT] [maximum_mutation_frequency]
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
        between the sorted haplotypes.

        This is expected to be calculated in a similar manner as:
            haplotype.build_distiance_matrix(haplotypes)


    RETURN
    ------

    [FLOAT] snde
        sample nucleotide diversity entity

    # ========================================================================
    """

    H = len(frequencies)
    D = distance_matrix

    if H > 1:
        snde = calculate.sample_nucleotide_diversity_entity(H, D)
    else:
        snde = constant.UNDEFINED

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
        between the sorted haplotypes.

        This is expected to be calculated in a similar manner as:
            haplotype.build_distiance_matrix(haplotypes)


    RETURN
    ------

    [FLOAT] [fad]
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
        between the sorted haplotypes.

        This is expected to be calculated in a similar manner as:
            haplotype.build_distiance_matrix(haplotypes)


    RETURN
    ------

    [FLOAT] [mutation_frequency]
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

    Returns the minimum mutation frequency of the sorted haplotypes.


    INPUT
    -----

    [PILEUP] [pileup]
        A Pileup object, which represents the pileup of aligned reads.

    [HAPLOTYPE LIST] [haplotypes]
        A list of sorted Haplotype objects.


    RETURN
    ------

    [FLOAT][minimum_mutation_frequency]
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
        A list of sorted Haplotype objects.


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
        A list of sorted Haplotype objects, for which to
        report the Shannon entropy.

    [FLOAT LIST] [frequencies]
        A list of (relative) frequencies of the Haplotypes.


    RETURN
    ------

    [FLOAT] [hs]
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
        A list of sorted Haplotype objects, for which to
        report the Shannon entropy.

    [FLOAT LIST] [frequencies]
        A list of (relative) frequencies of the Haplotypes.


    RETURN
    ------

    [FLOAT] [hsn]
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
         A list of sorted Haplotype objects, for which to
         report the Shannon entropy.

     [FLOAT LIST] [frequencies]
         A list of (relative) frequencies of the Haplotypes.


     RETURN
     ------

     [FLOAT] [hsn]
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

    [FLOAT] [simpson_index]

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

    [FLOAT] [gini_simpson_index]

    # ========================================================================
    """

    H = len(frequencies)
    P = frequencies
    gini_simpson_index = calculate.gini_simpson_index(H, P)

    return gini_simpson_index


def get_hill_numbers(frequencies, end, start=0):
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

    list_of_hill_numbers = []

    # Get an index pos between start to end(the range of hill numbers we want)
    for hill_number_pos in range(start, end):
        # Make sure a hill number is valid at that index
        # Hill num valid if len of frequency >= than the index pos+1.
        if H >= hill_number_pos + 1:
            list_of_hill_numbers.append(calculate.hill_number(
                H, P, hill_number_pos))

    return list_of_hill_numbers


def measurement_to_csv(measurements_list):

    measurements_col_titles = ["Position"]

    for i in range(len(constant.MEASUREMENTS_NAMES)):
        measurements_col_titles.append(constant.MEASUREMENTS_NAMES[i])

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

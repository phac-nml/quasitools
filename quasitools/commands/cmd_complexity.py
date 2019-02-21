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
import time
import copy
import csv
import quasitools.calculate as calculate
import quasitools.haplotype as haplotype
from quasitools.parsers.reference_parser import parse_references_from_fasta
from quasitools.parsers.mapped_read_parser import parse_haplotypes_called

UNDEFINED = "-"

# Measurement Positions
NUMBER_OF_HAPLOTYPES = 0
NUMBER_OF_POLYMORPHIC_SITES = 1
NUMBER_OF_MUTATIONS = 2
SHANNON_ENTROPY_NUMBER = 3
SHANNON_ENTROPY_NUMBER_LOCALIZED_TO_N = 4
SHANNON_ENTROPY_NUMBER_LOCALIZED_TO_H = 5
SIMPSON_INDEX = 6
GINI_SIMPSON_INDEX = 7
# HILL_NUMBER_0 = 8
# HILL_NUMBER_1 = 9
# HILL_NUMBER_2 = 10
# HILL_NUMBER_3 = 11
MINIMUM_MUTATION_FREQUENCY = 8  # 12
MUTATION_FREQUENCY_FREQUENCY = 9  # 13
FUNCTIONAL_ATTRIBUTE_DIVERSITY = 10  # 14
SAMPLE_NUCLEOTIDE_DIVERSITY_Entity = 11  # 15
MAXIMUM_MUTATION_FREQUENCY = 12  # 16
POPULATION_NUCLEOTIDE_DIVERSITY = 13  # 17
SAMPLE_NUCLEOTIDE_DIVERSITY = 14  # 18


# Dictionary of Names
MEASUREMENTS_NAMES = {
    'NUMBER_OF_HAPLOTYPES': "Number of Haplotypes",
    'NUMBER_OF_POLYMORPHIC_SITES': "Number of Polymorphic Sites",
    'NUMBER_OF_MUTATIONS': "Number of Mutations",
    'SHANNON_ENTROPY_NUMBER': "Shannon Entropy",
    'SHANNON_ENTROPY_NUMBER_LOCALIZED_TO_N': "Shannon Entropy Localized to N",
    'SHANNON_ENTROPY_NUMBER_LOCALIZED_TO_H': "Shannon Entropy Localized to H",
    'SIMPSON_INDEX': 'Simpson Index',
    'GINI_SIMPSON_INDEX': "Gini Simpson Index",
    # 'HILL_NUMBER_0': "Hill Number #0",
    # 'HILL_NUMBER_1': "HIll Number #1",
    # 'HILL_NUMBER_2': "Hill Number #2",
    # 'HILL_NUMBER_3': "Hill Number #3",
    'MINIMUM_MUTATION_FREQUENCY': "Minimum Mutation Frequency",
    'MUTATION_FREQUENCY_FREQUENCY': "Mutation Frequency",
    'FUNCTIONAL_ATTRIBUTE_DIVERSITY': "Functional Attribute Diversity",
    'SAMPLE_NUCLEOTIDE_DIVERSITY_Entity': "Sample Nucleotide Diversity Entity",
    'MAXIMUM_MUTATION_FREQUENCY': "Maximum Mutation Frequency",
    'POPULATION_NUCLEOTIDE_DIVERSITY': "Population Nucleotide Diversity",
    'SAMPLE_NUCLEOTIDE_DIVERSITY': "Sample Nucleotide Diversity",
}

# Creates an empty list.
EMPTY_LIST = ["-" for x in range(len(MEASUREMENTS_NAMES))]


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
    start_time = time.time()
    click.echo("\nStarting...")

    click.echo("Using file %s as reference" % reference)
    click.echo("Reading input from file(s)  %s" % bam)
    click.echo("Checking for overlaps at every position in reference" +
               " at %s intervals" % k)

    complexity(reference, bam, int(k))

    click.echo("\nComplete!")
    print("\nTime to run--- %s seconds ---" % (time.time() - start_time))


def complexity(reference, bam, k):
    """
    # ========================================================================

    COMPLEXITY


    PURPOSE
    -------

    Computes various complexity measures for the quasispecies.


    INPUT
    -----

    [(FASTA) FILE LOCATION] [fasta]
        The file location of an aligned FASTA file for which to calculate the
        complexity measures.


    RETURN
    ------

    [NONE]

    Comments
    ----

    The complexity computation and reporting will be completed once
    this method has run.

    # ========================================================================
    """

    references = parse_references_from_fasta(reference)
    haplotype_list = parse_haplotypes_called(
        references, reference, bam, k)

    # measurements_list = [ [] for x in range(len(haplotype_list)-k+1)]
    measurements_list = [[] for x in range(10)]
    measurements = [0 for x in range(len(MEASUREMENTS_NAMES))]

    # for i in range(len(haplotype_list)-k+1):
    # TODO The 10 is only for testing purpose replace with code above for
    # production.
    for i in range(10):
        haplotypes = haplotype_list[i]

        if not haplotypes:
            continue
        pileup = haplotype.build_pileup_from_haplotypes(haplotypes)

        distance_matrix = haplotype.build_distiance_matrix(haplotypes)
        counts = haplotype.build_counts(haplotypes)
        frequencies = haplotype.build_frequencies(haplotypes)

        '''
        Set the Incidence - Entity Level
        '''
        measurements[NUMBER_OF_HAPLOTYPES] = \
            get_number_of_haplotypes(haplotypes)

        measurements[NUMBER_OF_POLYMORPHIC_SITES] = \
            get_number_of_polymorphic_sites(pileup)

        measurements[NUMBER_OF_MUTATIONS] = get_number_of_mutations(pileup)

        '''
        Set the Abundance - Molecular Level
        '''

        shannon_entropy = get_shannon_entropy(haplotypes, frequencies)

        measurements[SHANNON_ENTROPY_NUMBER] = shannon_entropy

        measurements[SHANNON_ENTROPY_NUMBER_LOCALIZED_TO_N] = \
            get_shannon_entropy_localized_to_n(haplotypes, shannon_entropy)

        measurements[SHANNON_ENTROPY_NUMBER_LOCALIZED_TO_H] = \
            get_shannon_entropy_localized_to_h(
            haplotypes, shannon_entropy)

        measurements[SIMPSON_INDEX] = \
            get_simpson_index(frequencies)

        measurements[GINI_SIMPSON_INDEX] = \
            get_gini_simpson_index(frequencies)

        # returns a list of hill numbers
        # hill_numbers =  get_hill_numbers(frequencies)

        # for i in range (len(hill_numbers)):
        #    for hill_number_pos in range (HILL_NUMBER_0, HILL_NUMBER_3):
        #       measurements[hill_number_pos] = hill_numbers[i]

        '''
        Functional,  Indidence - Entity Level
        '''
        measurements[MINIMUM_MUTATION_FREQUENCY] = \
            get_minimum_mutation_frequency(
            haplotypes, pileup)

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
            distance_matrix, frequencies, haplotypes)

        measurements_list[i] = copy.deepcopy(measurements)

    for measurements in measurements_list:
        print(measurements)

    '''
    Measurement to CSV
    '''
    measurement_to_csv(measurements_list)


def get_sample_nucleotide_diversity(
        distance_matrix,
        frequencies,
        haplotypes):
    """"
    #========================================================================

    GET SAMPLE NUCLEOTIDE DIVERSITY


    PURPOSE
    -------

    Reports the diversity of nucleotides


    INPUT
    -------

    [2D ARRAY] [distance_matrix]
        A two dimensional array, representing the distance matrix of distances
        between the haplotypes.
    [FLOAT LIST] [frequencies]
        A list of (relative) frequencies of the Haplotypes.
    [HAPLOTYPE LIST] [haplotypes]
        A list of Haplotype objects.
    [DICTIONARY LIST] [measurements]
        Current Stateof measurements

    RETURN
    -------

    [DICTIONARY LIST] [measurements]

    # ========================================================================
    """

    N = haplotype.calculate_total_clones(haplotypes)
    H = len(frequencies)
    P = frequencies
    D = distance_matrix

    if N > 1:
        SND = calculate.sample_nucleotide_diversity(N, H, P, D)
    else:
        SND = UNDEFINED

    return SND


def get_population_nucleotide_diversity(
        distance_matrix, frequencies):
    """
    # ========================================================================
    REPORT POPULATION NUCLEOTIDE DIVERSITY

    PURPOSE
    -------
    Reports the population nucleotide diversity.

    INPUT
    -----
    [2D ARRAY] [distance_matrix]
        A two dimensional array, representing the distance matrix of distances
        between the haplotypes.
    [FLOAT LIST] [frequencies]
        A list of (relative) frequencies of the Haplotypes.


    RETURN
    ------

    [NONE]

    Comments
    ----

    The population sample nucleotide diversity will be reported to output.
    # ========================================================================
    """

    H = len(frequencies)
    P = frequencies
    D = distance_matrix

    PND = calculate.population_nucleotide_diversity(H, P, D)

    return PND


def get_maximum_mutation_frequency(
        counts,
        distance_matrix,
        frequencies):
    """
    # ========================================================================

    GET MAXMIMUM MUTATION FREQUENCY

    PURPOSE
    -------

    gets the maximum mutation frequency of the haplotypes.

    INPUT
    -----

    [INT LIST] [counts]
        A sorted list of haplotype counts, from the counts of the most
        abundant to the counts of the least abundant haplotype.
    [2D ARRAY] [distance_matrix]
        A two dimensional array, representing the distance matrix of distances
        between the haplotypes.
    [FLOAT LIST] [frequencies]
        A list of (relative) frequencies of the Haplotypes.
    [DICTIONARY LIST] [measurements]

    RETURN
    ------

    [DICTIONARY LIST] [measurements]


    POST
    ----
    The maximum mutation frequency will be reported to output.
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

    gets the entity-level sample nucleotide diversity.

    INPUT
    -----

    [2D ARRAY] [distance_matrix]
        A two dimensional array, representing the distance matrix of distances
        between the haplotypes.
        This is expected to be calculated in a similar manner as:
            haplotype.build_distiance_matrix(haplotypes)
    [FLOAT LIST] [frequencies]
        A list of (relative) frequencies of the Haplotypes.
    DICTIONARY LIST] [measurements]


    RETURN
    ------

    [DICTIONARY LIST] [measurements]


    # ========================================================================
    """

    H = len(frequencies)
    D = distance_matrix

    if H > 1:
        PND = calculate.sample_nucleotide_diversity_entity(H, D)
    else:
        PND = UNDEFINED

    return PND


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
    [DICTIONARY LIST] [measurements]

    RETURN
    ------

    [DICTIONARY LIST] [measurements]


    # ========================================================================
    """

    H = len(distance_matrix)
    D = distance_matrix

    FAD = calculate.FAD(H, D)

    return FAD


def get_mutation_frequency(distance_matrix):
    """
    # ========================================================================

    REPORT MUTATION FREQUENCY

    PURPOSE
    -------

    Reports the mutation frequency of the haplotypes.

    INPUT
    -----

    [2D ARRAY] [distance_matrix]
        A two dimensional array, representing the distance matrix of distances
        between the haplotypes.
    [DICTIONARY LIST] [measurements]

    RETURN
    ------

    [DICTIONARY LIST] [measurements]


    # ========================================================================
    """

    H = len(distance_matrix)
    D = distance_matrix
    mutation_frequency = calculate.mutation_frequency(H, D)

    return mutation_frequency


def get_minimum_mutation_frequency(haplotypes, pileup):
    """
    # ========================================================================

    REPORT MINIMUM MUTATION FREQUENCY

    PURPOSE
    -------
    Reports the minimum mutation frequency of the haplotypes.

    INPUT
    -----

    [PILEUP] [pileup]
        A Pileup object, which represents the pileup of aligned reads.
    [HAPLOTYPE LIST] [haplotypes]
        A list of Haplotype objects.
    [DICTIONARY LIST] [measurements]


    RETURN
    ------

    [DICTIONARY LIST] [measurements]


    # ========================================================================
    """

    M = pileup.count_unique_mutations()
    N = haplotype.calculate_total_clones(haplotypes)
    a = len(haplotypes[0].sequence)

    minimum_mutation_frequency = calculate.minimum_mutation_frequency(M, N, a)

    return minimum_mutation_frequency


def get_number_of_haplotypes(haplotypes):
    """""
    #========================================================================

    GET NUMBER OF HAPLOTYPES


    PURPOSE
    -------

    Reports the number of unique haplotypes


    INPUT
    -------

    [DICTIONARY LIST] [measurements]

    [NONE]


    RETURN
    -------

    [DICTIONARY LIST] [measurements]

    # ========================================================================
    """

    return len(haplotypes)


def get_number_of_polymorphic_sites(pileup):
    """""
    #========================================================================

    GET NUMBER OF POLYMORPHIC SITES


    PURPOSE
    -------

    Reports the number of polymorphic sites

    INPUT
    -------

    [DICTIONARY LIST] [measurements]
        Current measurements.

    [PILEUP] [pileup]
        A pileup object that represents the pileup of aligned reads

    RETURN
    -------

    [DICTIONARY LIST] [measurements]
        -Appended measurements after polymorphic sites added

    COMMENTS
    -------
    [N/A]

    #========================================================================
    """

    return pileup.count_polymorphic_sites()


def get_number_of_mutations(pileup):
    """""
    #========================================================================

    GET NUMBER OF MUTATIONS

    PURPOSE
    -------

    Gets the number of polymorphic sites

    INPUT
    -------

    [DICTIONARY LIST] [measurements]
        The current state of measurements.
    [PILEUP] [pileup]
        A pileup object that represents the pileup of aligned reads

    RETURN
    -------

    [DICTIONARY LIST] [measurements]
        The appended list of measurements now with number of mutations.

    COMMENTS
    -------
    [N/A]

    #========================================================================
    """

    return pileup.count_unique_mutations()


def get_shannon_entropy(haplotypes, frequencies):
    """""
    #========================================================================

    GET SHANNON ENTROPY

    PURPOSE
    -------

    Reports the shannon Entropy of the haplotpyes

    INPUT
    -------

    [FLOAT LIST] [frequencies]
        - A list of (relative) frequencies) of the haplotpyes.

    [HAPLOTYPES] [haplotypes]
        - A list of the haplotypes

    [DICTIONARY LIST] [measurements]
        - The current state of measurements.

    RETURN
    -------
    [DICTIONARY LIST] [measurements]
        - The appended list of measurements now with shannon entropy.

    COMMENTS
    -------
    [N/A]

    #========================================================================
    """

    Hs = calculate.shannon_entropy(frequencies)

    return Hs


def get_shannon_entropy_localized_to_n(haplotypes, Hs):
    """""
    #========================================================================

    GET SHANNON ENTROPY

    PURPOSE
    -------

    Reports the shannon Entropy of the haplotpyes

    INPUT
    -------

    [FLOAT LIST] [frequencies]
        - A list of (relative) frequencies) of the haplotpyes.

    [HAPLOTYPES] [haplotypes]
        - A list of the haplotypes

    [DICTIONARY LIST] [measurements]
        - The current state of measurements.

    RETURN
    -------
    [DICTIONARY LIST] [measurements]
        - The appended list of measurements now with shannon entropy.

    COMMENTS
    -------
    [N/A]

    #========================================================================
    """

    N = haplotype.calculate_total_clones(haplotypes)

    if float(math.log(N)) != 0:
        Hsn = float(Hs) / float(math.log(N))  # entropy localized to N
    else:
        Hsn = 0

    return Hsn


def get_shannon_entropy_localized_to_h(haplotypes, Hs):
    """""
    #========================================================================

    GET SHANNON ENTROPY

    PURPOSE
    -------

    Reports the shannon Entropy of the haplotpyes

    INPUT
    -------

    [FLOAT LIST] [frequencies]
        - A list of (relative) frequencies) of the haplotpyes.

    [HAPLOTYPES] [haplotypes]
        - A list of the haplotypes

    [DICTIONARY LIST] [measurements]
        - The current state of measurements.

    RETURN
    -------
    [DICTIONARY LIST] [measurements]
        - The appended list of measurements now with shannon entropy.

    COMMENTS
    -------
    [N/A]

    #========================================================================
    """

    H = len(haplotypes)

    if float(math.log(H)) != 0:
        Hsh = float(Hs) / float(math.log(H))  # entropy localized to H
    else:
        Hsh = 0

    return Hsh


def get_simpson_index(frequencies):
    """""
    #========================================================================

    REPORT SIMPSON INDEX

    PURPOSE
    -------

    Reports the Simpson Index

    INPUT
    -------
    [FLOAT LIST] frequencies
        A list of (relative) frequencies) of the haplotpyes.

    [DICTIONARY LIST] [measurements]
        The current state of measurements.

    RETURN
    -------

    [DICTIONARY LIST] [measurements]
        The appended list of measurements now with number of mutations.

    COMMENTS
    -------
    [N/A]

    #========================================================================
    """

    H = len(frequencies)
    P = frequencies

    simpson_index = calculate.simpson_index(H, P)

    return simpson_index


def get_gini_simpson_index(frequencies):
    """""
    #========================================================================

    GET GINI SIMPSON INDEX

    PURPOSE
    -------
    Reports the Gini-Simpson Index

    INPUT
    -------

    [FLOAT LIST] [frequencies]
        a list of relative frequencies of the haplotypes
    [DICTIONARY LIST] [measurements]

    RETURN
    -------

    [DICTIONARY LIST] [measurements]
        - The appended list of measurements now with the gini simpson index.

    COMMENTS
    -------

    #========================================================================
    """

    H = len(frequencies)
    P = frequencies
    gini_simpson_index = calculate.gini_simpson_index(H, P)

    return gini_simpson_index


def get_hill_numbers(frequencies):
    """""
    #========================================================================

    GET HILL NUMBERS

    PURPOSE
    -------

    Reports the 0, 1, 2, and 3 hill numbers.

    INPUT
    -------

    [FLOAT LIST] [frequencies]
        a list of relative frequencies of the haplotypes.
    [DICTIONARY LIST] [measurements]

    RETURN
    -------

    [DICTIONARY LIST] [measurements]
        - The appended list of measurements now with the hill numbers.

    COMMENTS
    -------

    #========================================================================
    """

    P = frequencies
    H = len(frequencies)

    NUMBER_OF_HILL_NUMBERS = 3

    # initialize list of hill numbers from 0 to 3
    list_of_hill_numbers = [0 for x in range(NUMBER_OF_HILL_NUMBERS)]

    for hill_number_pos in range(NUMBER_OF_HILL_NUMBERS):
        if len(P) >= hill_number_pos + 1:
            list_of_hill_numbers[hill_number_pos] = calculate.hill_number(
                H, P, hill_number_pos)

    # if len(P) >= 1:
    #     hill0 = calculate.hill_number(H, P, 0)
    # if len(P) >= 2:
    #     hill1 = calculate.hill_number(H, P, 1)
    # if len(P) >= 3:
    #     hill2 = calculate.hill_number(H, P, 2)
    # if len(P) >= 4:
    #     hill3 = calculate.hill_number(H, P, 3)

    return list_of_hill_numbers


def measurmentSummary(measurements):
    """""
    #========================================================================

    Measurement Summary

    PURPOSE
    -------

    Print all measurements.

    INPUT
    -------

    [DICTIONARY LIST] [measurements]

    RETURN
    -------

    [N/A]

    COMMENTS
    -------
    Checks keys, if key belongs to new measurement header we'll print
    the line break, and a the title of the header. Else just print the
    measurment.
    Accessing the measurement values in a for loop resulted in the each
    value being stored around '' within a {}.
    I stripped those out for a neater print.
    #========================================================================
    """

    for measurement in measurements:
        key = str(measurement.keys())
        if key == "dict_keys(['Number of haplotypes'])":
            click.echo("")
            click.echo("-----------------------------------------")
            click.echo("Incidence - Entity Level")
            click.echo("-----------------------------------------")
            click.echo(str(measurement).strip("{}").replace("'", ""))
        elif key == "dict_keys(['Shannon Entropy (Hs)'])":
            click.echo("")
            click.echo("-----------------------------------------")
            click.echo("Abundance - Molecular Level")
            click.echo("-----------------------------------------")
            click.echo(str(measurement).strip("{}").replace("'", ""))
        elif key == "dict_keys(['Minimum mutation frequency'])":
            click.echo("")
            click.echo("-----------------------------------------")
            click.echo("Functional, Indidence - Entity Level")
            click.echo("-----------------------------------------")
            click.echo(str(measurement).strip("{}").replace("'", ""))
        elif key == "dict_keys(['Maximum mutation frequency (Mf max)'])":
            click.echo("")
            click.echo("-----------------------------------------")
            click.echo("Functional, Abundance - Molecular Level")
            click.echo("-----------------------------------------")
            click.echo(str(measurement).strip("{}").replace("'", ""))
        elif key == "dict_keys(['Sample Nucleotide Diversity (^PI)'])":
            click.echo("")
            click.echo("-----------------------------------------")
            click.echo("Miscellaneous")
            click.echo("-----------------------------------------")
            click.echo(str(measurement).strip("{}").replace("'", ""))
        else:
            click.echo(str(measurement).strip("{}").replace("'", ""))


def measurement_to_csv(measurements_list):

    measurements_col_titles = ["Position"]

    # Itterates length of the dictionary 0 to 14
    for i in range(len(MEASUREMENTS_NAMES)):
        for key, value in MEASUREMENTS_NAMES.items():
            if eval(key) == i:
                measurements_col_titles.append(value)

    print(measurements_col_titles)

    file_name = "complexity_outputs.csv"

    for position in range(len(measurements_list)):

        # if measurement list at position i is not empty
        measurements = [UNDEFINED for x in range(1000)]
        if measurements_list[position]:
            measurements = copy.deepcopy(measurements_list[position])

        with open(file_name, 'a') as complexity_data:

            writer = csv.writer(complexity_data)
            # if file empty add the column titles (will be first row)
            if os.stat(file_name).st_size == 0:
                writer.writerow(measurements_col_titles)
                # will always add the measurements values
            writer.writerow([position] + measurements)
        complexity_data.close()

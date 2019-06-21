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

import click
import math
import csv
import sys

import quasitools.calculate as calculate
import quasitools.haplotype as haplotype
import quasitools.constants.complexity_constants as constant
from quasitools.parsers.reference_parser import parse_references_from_fasta
from quasitools.parsers.mapped_read_parser import \
    parse_haplotypes_from_bam, \
    parse_haplotypes_from_fasta

# Here we are setting up a nesting of Click commands. This allows us to run
# quasitools complexity only when grouped with a subcommand (BAM or FASTA)
@click.group(invoke_without_command=False)
@click.pass_context
def cli(ctx):
    '''
    Reports the per-amplicon (fasta) or k-mer complexity of the pileup,
    for each k-mer position in the reference complexity (bam and reference)
    of a quasispecies using several measures outlined in the following work:

    Gregori, Josep, et al. "Viral quasispecies complexity measures."
    Virology 493 (2016): 227-237.
    '''

    pass

# Subcommand for when a multi-alligned fasta file is provided
# When the fasta subcommand is called we will obtain complexity
# report of per-amplicon complexity.


@cli.command(
    'fasta', short_help='Calculates various quasispecies complexity ' +
    'measures on a multiple aligned FASTA file.')
@click.argument('fasta_location', nargs=1,
                type=click.Path(exists=True, file_okay=True, dir_okay=False))
@click.option('-o', '--output_location', type=click.Path(exists=False),
              help="Output the " +
              "quasispecies complexity in CSV format to the specified file.")
def fasta(fasta_location, output_location):
    '''
    Reports the per-amplicon (fasta)  or k-mer complexity of the pileup,
    for each k-mer position in the reference complexity (bam and reference)
    of a quasispecies using several measures outlined in the following work:

    Gregori, Josep, et al. "Viral quasispecies complexity measures."
    Virology 493 (2016): 227-237.
    '''

    """
    # ========================================================================

    FASTA COMPLEXITY


    PURPOSE
    -------

    Creates a report of k-mer complexity of the pileup, for each k-mer position
    in the reference.

    INPUT
    -----

    [(FASTA) FILE LOCATION] [fasta_location]
        The file location of an aligned FASTA file for which to calculate the
        complexity measures.
    [(OUTPUT) FILE LOCATION] [output_location]
        The location of the output file.

    RETURN
    ------

    [NONE]


    POST
    ----

    The complexity computation will be completed and the results will be
    stored in CSV file.

    # ========================================================================
    """

    haplotypes = parse_haplotypes_from_fasta(fasta_location)

    measurements = measure_complexity(haplotypes)

    # if the output_location is specificed open it as complexit_file, if not
    # specified, complexity_file is set as sys.stdout.
    with open(output_location, 'w') if output_location else sys.stdout as \
            complexity_file:
        measurement_to_csv([measurements], complexity_file)


# NGS Data from BAM and its corresponding reference file.
# When the bam subcommand is called we will produce a report of
# k-mer complexity of the pileup, for each k-mer position in the reference
@cli.command(
    'bam', short_help="Calculates various quasispecies complexity " +
    "measures on next generation sequenced data from a BAM file " +
    "and it's corresponding reference file.")
@click.argument('reference_location', required=True,
                type=click.Path(exists=True, file_okay=True, dir_okay=False))
@click.argument('bam_location', nargs=1,
                type=click.Path(exists=True, file_okay=True, dir_okay=False))
@click.argument('k')
@click.argument('haplotype_threshold')
@click.option('-o', '--output_location', type=click.Path(exists=False),
              help="Output the " +
              "quasispecies complexity in CSV format to the specified file.")
def bam(reference_location, bam_location, k,
        haplotype_threshold, output_location):
    '''
    Reports the per-amplicon (fasta)  or k-mer complexity of the pileup,
    for each k-mer position in the reference complexity (bam and reference)
    of a quasispecies using several measures outlined in the following work:

    Gregori, Josep, et al. "Viral quasispecies complexity measures."
    Virology 493 (2016): 227-237.
    '''

    """
    # ========================================================================

    BAM COMPLEXITY


    PURPOSE
    -------

    Create a report of  k-mer complexity of the pileup, for each k-mer position
    in the reference.


    INPUT
    -----

    [(BAM) FILE LOCATION] [bam_location]
        The file location of a bam file.
    [(REFERENCE) FILE LOCATION] [reference_location]
        The file location of the reference file.
    [INT] k
        Provides the sequence length for our reads from a given starting
        position.:
    [FLOAT] haplotype_threshold:
        The threshold percentage that a haplotype must appear for
              it to be included in each positional haplotype list.

    [(OUTPUT) FILE LOCATION] [output_location]
        The location of the output file.

    RETURN
    ------

    [NONE]


    POST
    ----

    The complexity computation will be completed and the results will be stored
    in CSV file or std.out.

    # ========================================================================
    """
    k = int(k)
    percent = 100

    haplotype_threshold = float(haplotype_threshold) / float(percent)

    references = parse_references_from_fasta(reference_location)
    # A list where each position contains a list of haplotypes of length k
    # starting at that position in the reference.
    haplotype_list = parse_haplotypes_from_bam(
        references, reference_location, bam_location, k)

    measurements_list = []

    for i in range(len(haplotype_list)):

        haplotypes = haplotype_list[i]
        print("Haplotypes at Position: " + str(i) + " Before filter")
        for hap_not_filtered in haplotypes:
            print("Haplotype: " + str(hap_not_filtered.sequence) +
                  " Count: " + str(hap_not_filtered.count))
        # Remove haplotypes below threshold.

        # Get total number of haplotypes for each position.
        total_haplotypes = haplotype.calculate_total_clones(haplotypes)
        # Add haplotypes within threshold to new haplotypes list
        haplotypes_in_threshold = [hap for hap in haplotypes if (
            float(hap.count) / float(total_haplotypes)) > haplotype_threshold]

        print("Haplotypes at Position: " + str(i) + " After filter")
        for haplotypes_filtered in haplotypes_in_threshold:
            print("Haplotype: " + str(haplotypes_filtered.sequence) +
                  " Count: " + str(haplotypes_filtered.count))
        print("\n")
        print("\n")
        print("\n")

        measurements = measure_complexity(haplotypes_in_threshold)
        measurements_list.append(measurements)

    # if the output_location is specificed open it as complexit_file, if not
    # specified, complexity_file is set as sys.stdout.
    with open(output_location, 'w') if output_location else sys.stdout as \
            complexity_file:
        measurement_to_csv(measurements_list, complexity_file)


def measure_complexity(haplotypes):
    """"
    #========================================================================

    MEASURE COMPLEXITY

    PURPOSE
    -------

    Calculate a number of complexity measurements for Haplotype objects.


    INPUT
    -------

    [HAPLOTYPE LIST] [haplotypes]
        - An unsorted  list of Haplotype objects.


    RETURN
    -------

    [LIST] [measurements]
        - A list of complexity measurements. List position is specified by
          constant values defined in constants/complexity_constanty.py

    #========================================================================
    """

    # Initialize measurements list to length of the measurements name
    # dictionary.

    measurements = [constant.UNDEFINED for x in range(
        len(constant.MEASUREMENTS_NAMES))]

    # If no haplotypes we will return measurements as its initialized state.
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

    # Set the Incidence - Entity Level #

    measurements[constant.NUMBER_OF_HAPLOTYPES] = \
        get_number_of_haplotypes(sorted_haplotypes)

    measurements[constant.HAPLOTYPE_POPULATION] = \
        get_haplotype_population(counts)

    measurements[constant.NUMBER_OF_POLYMORPHIC_SITES] = \
        get_number_of_polymorphic_sites(pileup)

    measurements[constant.NUMBER_OF_MUTATIONS] = get_number_of_mutations(
        pileup)

    # Set the Abundance - Molecular Level: #

    shannon_entropy = get_shannon_entropy(sorted_haplotypes, frequencies)

    measurements[constant.SHANNON_ENTROPY_NUMBER] = shannon_entropy

    measurements[constant.SHANNON_ENTROPY_NUMBER_NORMALIZED_TO_N] = \
        get_shannon_entropy_normalized_to_n(
        sorted_haplotypes, shannon_entropy)

    measurements[constant.SHANNON_ENTROPY_NUMBER_NORMALIZED_TO_H] = \
        get_shannon_entropy_normalized_to_h(
        sorted_haplotypes, shannon_entropy)

    measurements[constant.SIMPSON_INDEX] = \
        get_simpson_index(frequencies)

    measurements[constant.GINI_SIMPSON_INDEX] = \
        get_gini_simpson_index(frequencies)

    hill_numbers_list = get_hill_numbers(
        frequencies, constant.HILL_NUMBER_LENGTH)

    '''
     We iterate through the number of elements in the Hill number list,
     at each iteration we place the Hill number at element i into measurements
     at element HILL_NUMBER_0+i.
    '''
    for i in range(len(hill_numbers_list)):
        measurements[constant.HILL_NUMBER_0 + i] = (hill_numbers_list[i])

    # Functional Incidence - Entity Level #

    measurements[constant.MINIMUM_MUTATION_FREQUENCY] = \
        get_minimum_mutation_frequency(
        sorted_haplotypes, pileup)

    measurements[constant.MUTATION_FREQUENCY] = get_mutation_frequency(
        distance_matrix)

    measurements[constant.FUNCTIONAL_ATTRIBUTE_DIVERSITY] = \
        get_FAD(distance_matrix)

    measurements[constant.SAMPLE_NUCLEOTIDE_DIVERSITY_ENTITY] = \
        get_sample_nucleotide_diversity_entity(
        distance_matrix, frequencies)

    # Functional Abundance - Molecular Level #

    measurements[constant.MAXIMUM_MUTATION_FREQUENCY] = \
        get_maximum_mutation_frequency(
        counts, distance_matrix, frequencies)

    measurements[constant.POPULATION_NUCLEOTIDE_DIVERSITY] = \
        get_population_nucleotide_diversity(
        distance_matrix, frequencies)

    # Other #

    measurements[constant.SAMPLE_NUCLEOTIDE_DIVERSITY] = \
        get_sample_nucleotide_diversity(
        distance_matrix, frequencies, sorted_haplotypes)

    return measurements


def get_haplotype_population(counts):
    """
    # ========================================================================

    GET HAPLOTYPE POPULATION


    PURPOSE
    -------

    Returns the population of haplotypes.


    INPUT
    -----

    [INT LIST] [counts]
        A haplotype counts, from the counts of the most abundant to the counts
        of the least abundant haplotype.


    RETURN
    ------

    [INT]
        The haplotype population which is defined as the number of
        unique haplotypes.

    # ========================================================================
    """

    return sum(counts)


def get_sample_nucleotide_diversity(
        distance_matrix,
        frequencies,
        haplotypes):
    """
    # ========================================================================

    GET SAMPLE NUCLEOTIDE DIVERSITY


    PURPOSE
    -------

    Reports the sample nucleotide diversity.


    INPUT
    -----

    [HAPLOTYPE LIST] [haplotypes]
        A list of Haplotype objects.

    [FLOAT LIST] [frequencies]
        A list of (relative) frequencies of the haplotypes.

    [2D ARRAY] [distance_matrix]
        A two dimensional array, representing the distance matrix of distances
        between the haplotypes.

        This is expected to be calculated in a similar manner as:
            haplotype.build_distiance_matrix(haplotypes)

    RETURN
    ------

    [FLOAT] [snd]
        The sample nucleotide diversity.

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
        The population nucleotide diversity.

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
        A haplotype counts, from the counts of the most abundant to the counts
        of the least abundant haplotype.

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
        The maximum mutation frequency.

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
        A list of (relative) frequencies of the haplotypes.

    [2D ARRAY] [distance_matrix]
        A two dimensional array, representing the distance matrix of distances
        between the sorted haplotypes.

        This is expected to be calculated in a similar manner as:
            haplotype.build_distiance_matrix(haplotypes)


    RETURN
    ------

    [FLOAT] snde
        Sample nucleotide diversity entity.

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
        The functional attribute diversity.

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
        The mutation frequency.

    # ========================================================================
    """

    H = len(distance_matrix)
    D = distance_matrix
    mutation_frequency = calculate.mutation_frequency(H, D)

    return mutation_frequency


def get_minimum_mutation_frequency(haplotypes, pileup):
    """
    # ========================================================================

    GET  MINIMUM MUTATION FREQUENCY


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
        The minimum mutation frequency.

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

    Returns the number of (unique) Haplotype objects.


    INPUT
    -----

    [HAPLOTYPE LIST] [haplotypes]
        A list of sorted Haplotype objects.


    RETURN
    ------

    [INT] [len(haplotypes)]
        The unique number of Haplotype objects.

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
        The number of polymorphic sites in pileuip.

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
        The number of mutations.

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
        The Shannon Entropy of the haplotypes, determined from their
        frequencies.

    # ========================================================================
    """

    Hs = calculate.shannon_entropy(frequencies)

    return Hs


def get_shannon_entropy_normalized_to_n(haplotypes, Hs):
    """
    # ========================================================================

    GET SHANNON ENTROPY NORMALIZED TO N


    PURPOSE
    -------

    Returns the Shannon entropy of the haplotypes normalized to N.


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
        The Shannon entropy, normalized to the number of clones.
        This is the sum of all counts of each haplotype.

    # ========================================================================
    """

    N = haplotype.calculate_total_clones(haplotypes)

    if float(math.log(N)) != 0:
        Hsn = float(Hs) / float(math.log(N))  # entropy normalized to N.
    else:
        Hsn = 0

    return Hsn


def get_shannon_entropy_normalized_to_h(haplotypes, Hs):
    """
     # ========================================================================

     GET SHANNON ENTROPY NORMALIZED TO H


     PURPOSE
     -------

     Returns the Shannon entropy of the haplotypes normalized to N.


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
        The Shannon entropy, normalized to the number of haplotypes.
        This is the length of the haplotypes list.

     # ========================================================================
     """

    H = len(haplotypes)

    if float(math.log(H)) != 0:
        Hsh = float(Hs) / float(math.log(H))  # entropy normalized to H.
    else:
        Hsh = 0

    return Hsh


def get_simpson_index(frequencies):
    """
    # ========================================================================

    GET SIMPSON INDEX


    PURPOSE
    -------

    Returns the Simpson Index.


    INPUT
    -----

    [FLOAT LIST] [frequencies]
        A list of (relative) frequencies of the Haplotypes.


    RETURN
    ------

    [FLOAT] [simpson_index]
        The Simpson Index.

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
        The Gini-Simpson index.

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

    Return the 0, 1, 2, and 3 Hill numbers.


    INPUT
    -----

    [FLOAT LIST] [frequencies]
        A list of (relative) frequencies of the Haplotypes.

    [INT] [start]
        An integer value to designate the starting position of our
        Hill numbers list.
    [INT] [end]
        An integer value to designate the end position of our
        Hill numbers list.


    RETURN
    ------

    [FLOAT LIST] [list_of_hill_numbers]
        A list of Hill numbers.

    # ========================================================================
    """

    P = frequencies
    H = len(frequencies)

    list_of_hill_numbers = []

    # Iterate over the range of hill numbers.
    for hill_number_pos in range(start, end):
        # Make sure a hill number is valid at that index
        # A Hill number (Hn) can be calculate if we have at least Hn + 1
        # haplotypes. That is, the length the haplotype frequencies list (H)
        # must be at least (Hn + 1).
        if H >= hill_number_pos + 1:
            list_of_hill_numbers.append(calculate.hill_number(
                H, P, hill_number_pos))
        else:
            break

    return list_of_hill_numbers


def measurement_to_csv(measurements_list, complexity_file):
    """
    # ========================================================================

    MEASUREMENTS TO CSV


    PURPOSE
    -------

    Report a number of complexity measurements.


    INPUT
    -----

    [[List]] [measurements_list]
        A two dimensional list that contains a number of complexity
        measurements each position contains haplotypes of length position + k.
    [(OUTPUT) FILE] [complexity_file]
        The output file.

    POST
    ------
        A CSV file that reports a number of complexity measurements for a
        starting postion to k (BAM) or per amplicon (FASTA). Data is written
        into complexity_file which is either a user defined file or sys.stdout
        when the file is not specified.

    # ========================================================================
    """

    measurements_col_titles = ["Position"]

    for i in range(len(constant.MEASUREMENTS_NAMES)):
        measurements_col_titles.append(constant.MEASUREMENTS_NAMES[i])

    writer = csv.writer(complexity_file)
    writer.writerow(measurements_col_titles)

    for position in range(len(measurements_list)):
        measurements = measurements_list[position]
        writer.writerow([position] + measurements)

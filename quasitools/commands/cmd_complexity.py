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
import quasitools.calculate as calculate
import quasitools.haplotype as haplotype
from quasitools.parsers.reference_parser import parse_references_from_fasta

# Remove this
import time

import csv

from quasitools.parsers.mapped_read_parser import parse_haplotypes_called

UNDEFINED = ""


@click.command(
    'complexity', short_help='Calculates various quasispecies complexity \
    measures.')
@click.argument('reference', nargs=1, required=True,
                type=click.Path(exists=True, file_okay=True, dir_okay=False))
@click.argument('bam', nargs=1,
                type=click.Path(exists=True, file_okay=True, dir_okay=False))
# TODO is their a better name than k?
@click.argument('k')
@click.pass_context
def cli(ctx, reference, bam, k):
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
    click.echo("Checking for overlaps at every position in reference \
            at %s intervals" % k)

    complexity(ctx, reference, bam, int(k))

    click.echo("\nComplete!")
    print("\nTime to run--- %s seconds ---" % (time.time() - start_time))


def complexity(ctx, reference, bam, k):
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

        measurements = []

        '''
        Set the Incidence - Entity Level
        '''
        measurements = get_number_of_haplotypes(haplotypes, measurements)
        measurements = get_number_of_polymorphic_sites(measurements, pileup)
        measurements = get_number_of_mutations(measurements, pileup)

        '''
        Set the Abundance - Molecular Level
        '''
        measurements = measurements = get_shannon_entropy(
            haplotypes, frequencies, measurements)
        measurements = get_simpson_index(frequencies, measurements)
        measurements = get_gini_simpson_index(frequencies, measurements)

        measurements = get_hill_numbers(measurements, frequencies)

        '''
        Functional,  Indidence - Entity Level
        '''
        measurements = get_minimum_mutation_frequency(
            haplotypes, measurements, pileup)
        measurements = get_mutation_frequency(distance_matrix, measurements)
        measurements = get_FAD(distance_matrix, measurements)
        measurements = get_sample_nucleotide_diversity_entity(
            distance_matrix, frequencies, measurements)
        '''

        Functional, Abundance - Molecular Level
        '''
        measurements = get_maximum_mutation_frequency(
            counts, distance_matrix, frequencies, measurements)
        measurements = get_population_nucleotide_diversity(
            distance_matrix, frequencies, measurements)
        '''

        Other
        '''
        measurements = get_sample_nucleotide_diversity(
            distance_matrix, frequencies, haplotypes, measurements)

        '''
        Measurement Summary
        '''
        measurmentSummary(measurements)

        # creates two line break
        click.echo('\n' * 2)

        '''
        Measurement to CSV
        '''
        measurement_to_csv(measurements, i)


def get_sample_nucleotide_diversity(
        distance_matrix,
        frequencies,
        haplotypes,
        measurements):
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
    measurements.append({'Sample Nucleotide Diversity (^PI)': SND})

    return measurements


def get_population_nucleotide_diversity(
        distance_matrix, frequencies, measurements):
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
    measurements.append({'Population Nucleotide Diversity (PI)': PND})

    return measurements


def get_maximum_mutation_frequency(
        counts,
        distance_matrix,
        frequencies,
        measurements):
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
    measurements.append(
        {'Maximum mutation frequency (Mf max)': maximum_mutation_frequency})

    return measurements


def get_sample_nucleotide_diversity_entity(
        distance_matrix, frequencies, measurements):
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

    measurements.append({'Population Nucleotide Diversity': PND})

    return measurements


def get_FAD(distance_matrix, measurements):
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

    measurements.append({'Functional attribute diversity': FAD})

    return measurements


def get_mutation_frequency(distance_matrix, measurements):
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

    measurements.append({'Mutation frequency': mutation_frequency})
    return measurements


def get_minimum_mutation_frequency(haplotypes, measurements, pileup):
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

    measurements.append(
        {'Minimum mutation frequency': minimum_mutation_frequency})

    return measurements


def get_number_of_haplotypes(haplotypes, measurements):
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

    haplotype_num = len(haplotypes)

    measurements.append({'Number of haplotypes': haplotype_num})

    return measurements


def get_number_of_polymorphic_sites(measurements, pileup):
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
    number_of_polymorphic_sites = pileup.count_polymorphic_sites()
    measurements.append(
        {'Number of polymorphic sites': number_of_polymorphic_sites})

    return measurements


def get_number_of_mutations(measurements, pileup):
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

    number_of_mutations = pileup.count_unique_mutations()
    measurements.append({'Number of mutation': number_of_mutations})

    return measurements


def get_shannon_entropy(haplotypes, frequencies, measurements):
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
    H = len(haplotypes)
    N = haplotype.calculate_total_clones(haplotypes)

    if float(math.log(N)) != 0:
        Hsn = float(Hs) / float(math.log(N))  # entropy localized to N
        Hsh = float(Hs) / float(math.log(H))  # entropy localized to H
    else:
        Hsn = 0
        Hsh = 0

    measurements.append({'Shannon Entropy (Hs)': Hs})
    measurements.append({'Shannon Entropy Localized to N (Hsn)': Hsn})
    measurements.append({'Shannon Entropy Localized to H (Hsh)': Hsh})

    return measurements


def get_simpson_index(frequencies, measurements):
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
    measurements.append({'Simpson Index': simpson_index})

    return measurements


def get_gini_simpson_index(frequencies, measurements):
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
    measurements.append({'Gini Simpson Index': gini_simpson_index})

    return measurements


def get_hill_numbers(measurements, frequencies):
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

    hill0 = hill1 = hill2 = hill3 = UNDEFINED

    if len(P) >= 1:
        hill0 = calculate.hill_number(H, P, 0)
    if len(P) >= 2:
        hill1 = calculate.hill_number(H, P, 1)
    if len(P) >= 3:
        hill2 = calculate.hill_number(H, P, 2)
    if len(P) >= 4:
        hill3 = calculate.hill_number(H, P, 3)

    measurements.append({'Hill Number 0': hill0})
    measurements.append({'Hill Number 1': hill1})
    measurements.append({'Hill Number 2': hill2})
    measurements.append({'Hill Numbwe 3': hill3})

    return measurements


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


def measurement_to_csv(measurements, count):
    '''
    What's going on:

        Passed an list of dictionaries (measurements) and the
        haplotype position (count). My goal is to print the CSV file.

    Expected output:
        The first row contains column titles and subsequent rows
        contain values of each measurement.

    Approach:
        Storing my measurments in a dictionary paid off here.
        I can extract the keys and store them in one list
        I can extract the values and store them in a second list

        Each time measure to csv is called I can append the list
        with the values list.

    Problem:

        I only want to add column titles once. and then only append values.
        Because haplotypes could start at any position from 0 to length
        I have to be careful of where I start.

    Fix:
        Since column titles is the first row, the file should be empty
        adding it. If it's empty it shouldn't exist yet therefore have size of
        0. if a file is empty we can add the column titles, if it's not
        empty we have already added column titles and can skip to
        adding values.

    Problem:
        If file exists in directory already (i.e it was run before) then
        we will just append to the old file. This could confuse people.

    Possible Fixes:
        - Warn user to remove old file if already created
          (not very user friendly)
        - Give user option to name file, and warn not to use
          same file name... can't trust user though.
        - Force a re-write everytime the program runs
            - what if user is testing multiple data sets.
        - Give user option to name the file and overwrite old files(?)
   '''
    # stores the column titles (keys)
    measurements_col_titles = ["Position"]
    # stores the values in each row (values)
    measurements_values = []

    file_name = "complexity_outputs.csv"

    # Goes through every measurment dictionary in the list.
    for measurement in measurements:

        # Extracts the key and value from each item.
        ''' A key is part of the column titles and get appended to
            measurements_col_titles. '''
        # A value is part of the subsequent rows and get aooebded to
        # measurements_values
        for key, value in measurement.items():

            measurements_col_titles.append(key)
            measurements_values.append(value)

    # We have two lists, we will append both or one (just values if column
    # titles exist).
    with open(file_name, 'a') as complexity_data:

        '''
        TODO maybe move the insert value outside of the file  writer
        within outer forloop.
        '''
        measurements_values.insert(0, count)

        # call the instance method that returns the delimited string,
        writer = csv.writer(complexity_data)

        # if file empty add the column titles (will be first row)
        if os.stat(file_name).st_size == 0:
            writer.writerow(measurements_col_titles)
        # will always add the measurements values
        writer.writerow(measurements_values)

    complexity_data.close()
    measurements_values.pop(0)

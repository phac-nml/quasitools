"""
Copyright Government of Canada 2018

Written by: Matthew Fogel, Public Health Agency of Canada

Licensed under the Apache License, Version 2.0 (the "License"); you may not use
this work except in compliance with the License. You may obtain a copy of the
License at:

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software distributed
under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR
CONDITIONS OF ANY KIND, either express or implied. See the License for the
specific language governing permissions and limitations under the License.
"""

import click
from quasitools.pileup import Pileup_List
from quasitools.distance import DistanceMatrix


@click.command('distance', short_help='Calculate the evolutionary distance'
               'between viral quasispecies using angular cosine distances.')
@click.argument('reference', nargs=1, required=True,
                type=click.Path(exists=True, file_okay=True, dir_okay=False))
@click.argument('bam', nargs=-1,
                type=click.Path(exists=True, file_okay=True, dir_okay=False))
@click.option('-n/-dn', '--normalize/--dont_normalize', default=True,
              help="Normalize read count data so that the read counts per " +
              "4-tuple (A, C, T, G) sum to one.")
@click.option('-od/-os', '--output_distance/--output_similarity', default=True,
              help="Output an angular distance matrix (by default), or " +
              "output a cosine similarity matrix (cosine similarity is not a" +
              " metric)")
@click.option('-s', '--startpos', type=int, help="Set the start base " +
              "position of the reference to use in the distance or " +
              "similarity calculation. Start position is one-indexed. " +
              " I.e. it must be between one and the total length of the " +
              "reference sequence(s), inclusive.")
@click.option('-e', '--endpos', type=int, help="Set the end base position" +
              " of the reference to use in the distance or similarity "
              "calculation. End position is one-indexed. " +
              " I.e. it must be between one and the total length of the " +
              "reference sequence(s), inclusive.")
@click.option('-o', '--output', type=click.File('w'), help="Output the " +
              "quasispecies distance or similarity matrix in CSV format in a" +
              " file.")
@click.option('-t', '--truncate', 'no_coverage', flag_value='truncate',
              default=True, help="Ignore contiguous start and end pileup" +
              " regions with no coverage.")
@click.option('-r', '--remove_no_coverage', 'no_coverage',
              flag_value='remove_no_coverage',
              help="Ignore all pileup regions with no coverage.")
@click.option('-k', '--keep_no_coverage', 'no_coverage',
              flag_value='keep_no_coverage',
              help="Do not ignore pileup regions with no coverage.")
@click.pass_context
def cli(ctx, reference, bam, normalize, output_distance, startpos, endpos,
        output, no_coverage):
    """Quasitools distance produces a measure of evolutionary distance [0 - 1]
       between quasispecies, computed using the angular cosine distance
       function defined below.

       Cosine similarity = (u * v) / ( ||u|| * ||v|| )

       Angular Cosine Distance = 2 * ACOS(Cosine similarity) / PI

       The tool outputs by default an angular cosine distance matrix.
       Use the flag defined below to instead output a similarity matrix.

       By default the data is normalized and start and end regions of the
       pileup with no coverage are truncated.

       It is possible to remove all no coverage pileup regions, including inner
       regions, or keep all no coverage regions in the pileup.

       Normalization is done dividing base read counts (A, C, T, G) inside
       every 4-tuple by the sum of the read counts inside the same tuple.
       The normalized read counts inside each 4-tuple sum to one."""

    click.echo("Using file %s as reference" % (reference))
    for file in bam:
        click.echo("Reading input from file(s)  %s" % (file))

    click.echo(dist(ctx, reference, bam, normalize, output_distance,
               startpos, endpos, output, no_coverage))


def dist(ctx, reference, bam, normalize, output_distance, startpos, endpos,
         output, no_coverage):

    """
    dist - Performs the main part of the program

    INPUT:
        [CONTEXT] [ctx]
        [FASTA FILE LOCATION] [reference]
        [BAM FILE LOCATION] [bam]
        [BOOL] [normalize/dont_normalize]
        [BOOL] [output_distance/output_similarity]
        [INT] [startpos]
        [INT] [endpos]
        [STRING] [output]
            Output the CSV-formatted matrix output in a file
            instead of in the terminal.
        [STRING] [truncate/remove_no_coverage/keep_no_coverage]
            Options to truncate low-coverage regions on the ends of the pileup,
            ignore all low coverage regions, or keep all low coverage regions

    RETURN:
        [STRING] [error message(s) if applicable. Otherwise "Complete!"]

    POST:
        The calling function, cli, will click.echo the string that is returned.

    """

    if len(bam) < 2:
        return ("\nError: At least two bam file locations are" +
                " required to perform quasispecies distance comparison")
    # indicate if the start or end position is < 0 or a priori invalid
    if type(startpos) == int and int(startpos) < 1:
        return "\nError: Start position must be >= 1."
    if type(endpos) == int and int(endpos) < 1:
        return "\nError: End position must be >= 1."
    if (type(startpos) == int and type(endpos) == int and (startpos > endpos)):
        return ("\nError: Start position must be <= end position")
    pileups = Pileup_List.construct_pileup_list(bam, reference)

    if startpos is None:
        startpos = 1
    if endpos is None:
        endpos = pileups.get_pileup_length()

    if pileups.get_pileup_length() == 0:
        return ("Error: Empty pileup was produced from BAM files." +
                "Halting program")
    click.echo("The start position is %d." % startpos)
    click.echo("The end position is %d." % endpos)
    click.echo("Constructed pileup from reference.")
    # click.echo the number of positions in pileup
    click.echo("The pileup covers %d positions before modifications." %
               pileups.get_pileup_length())
    # indicate whether the user-specified start and end position is out
    # of bounds (comparing to actual number of positions in pileup)
    if startpos > pileups.get_pileup_length():
        return ("\nError: Start position must be less than or equal to " +
                " number of nucleotide base positions in pileup" +
                " (%s)." % pileups.get_pileup_length())
    if endpos > pileups.get_pileup_length():
        return ("\nError: End position must be less than or equal to length " +
                "of nucleotide base positions in pileup" +
                " (%s)." % pileups.get_pileup_length())

    # we convert the start and end positions from one-based indexing to
    # zero-based indexing which is expected by distance.py and pileup.py
    startpos -= 1
    endpos -= 1

    # if there is no errors so far, proceed with running program
    modified = modify_pileups(ctx, normalize, startpos, endpos, no_coverage,
                              pileups)

    if (no_coverage is not 'keep_no_coverage') and (len(modified) == 0):
        return ("Error: Entire pileup was truncated due to " +
                "lack of coverage. Halting program")
    dist = DistanceMatrix(modified, bam)
    if output_distance:
        click.echo("Outputting an angular cosine distance matrix.")
        if output:
            output.write(dist.get_distance_matrix_as_csv())
        else:
            click.echo(dist.get_distance_matrix_as_csv())
    else:
        click.echo("Outputting a cosine similarity matrix.")
        if output:
            output.write(dist.get_similarity_matrix_as_csv())
        else:
            click.echo(dist.get_similarity_matrix_as_csv())
    # end if
    return "Complete!"
# end def


def modify_pileups(ctx, normalize, startpos, endpos, no_coverage, pileups):

    """
    modify_pileups - Performs normalization, truncation and/or selecting the
                     range of the pileup, if these options are enabled.

    INPUT:
        [CONTEXT] [ctx]
        [BOOL] [normalize/dont_normalize]
        [INT] [startpos]
        [INT] [endpos]
        [STRING] [no_coverage] [truncate/remove_no_coverage/keep_no_coverage]
            Options to truncate low-coverage regions on the ends of the pileup,
            ignore all low coverage regions, or keep all low coverage regions
        [PILEUP_LIST] [pileups]
            The pileups to be modified

    RETURN:
        [2D ARRAY] [pileups]
            The modified pileups, returned as a two-dimensional array.

    POST:
        The calling function, dist, will use the 2D array that is returned.

    """

    startpos = int(startpos)
    endpos = int(endpos)
    pileups.select_pileup_range(startpos, endpos)
    # converting startpos and endpos back to one-based indexing for click.echo
    click.echo(("The pileup covers %s positions after selecting " +
               "range between original pileup positions %d and %d.")
               % (pileups.get_pileup_length(), startpos + 1, endpos + 1))
    if normalize:
        pileups.normalize_pileups()
    old_length = pileups.get_pileup_length()
    if no_coverage is not 'keep_no_coverage':
        if no_coverage is 'truncate':
            pileups.truncate_output()
            click.echo("Truncating positions with no coverage that " +
                       "are contiguous with the start or end " +
                       "position of the pileup only.")
        elif no_coverage is 'remove_no_coverage':
            pileups.remove_no_coverage()
            click.echo("Truncating all positions with no coverage.")
        # end if
        click.echo("%d positions were truncated on the left" %
                   pileups.get_num_left_positions_truncated())
        click.echo("%d positions were truncated on the right" %
                   pileups.get_num_right_positions_truncated())
        click.echo("%d positions were removed in total from the pileup" %
                   (old_length - pileups.get_pileup_length()))
    # end if
    return pileups.get_pileups_as_numerical_array()

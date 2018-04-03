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
from quasitools.distance import Pileup_List
from quasitools.distance import DistanceMatrix


@click.command('distance', short_help='Calculate the evolutionary distance'
               'between viral quasispecies using cosine similarity.')
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
              "similarity calculation. Start positions must be greater than " +
              "zero and less than the length of the base pileup (number of " +
              "positions to be compared with the reference.)")
@click.option('-e', '--endpos', type=int, help="Set the end base position" +
              " of the reference to use in the distance or similarity "
              "calculation. End positions must be greater than zero and less" +
              " than the length of the base pileup (number of positions to" +
              " be compared with the reference.)")
@click.option('-o', '--output', type=click.File('w'), help="Output the " +
              "quasispecies distance or similarity matrix in CSV format in a" +
              " file.")
@click.option('-te', '--truncate_ends', 'truncate', flag_value='truncate_ends',
              default=True, help="Ignore contiguous start and end pileup" +
              " regions with no coverage.")
@click.option('-ta', '--truncate_all', 'truncate', flag_value='truncate_all',
              help="Ignore all pileup regions with no coverage.")
@click.option('-dt', '--dont_truncate', 'truncate', flag_value='dont_truncate',
              help="Do not ignore pileup regions with no coverage.")
@click.pass_context
def cli(ctx, reference, bam, normalize, output_distance, startpos, endpos,
        output, truncate):
    """Quasitools distance produces a measure of evolutionary distance [0 - 1]
       between quasispecies, computed using the angular cosine distance
       function defined below.

       Cosine similarity = (u * v) / ( ||u|| * ||v|| )

       Angular Cosine Distance = 2 * ACOS(Cosine similarity) / PI

       The tool outputs by default an angular cosine distance matrix.
       Use the flag defined below to instead output a similarity matrix.

       By default the data is normalized and start and end regions of the
       pileup with no coverage are truncated.

       It is possible to truncate all pileup regions, including inner
       regions, with no coverage, or turn truncation off completely.

       Normalization is done dividing base read counts (A, C, T, G) inside
       every 4-tuple by the sum of the read counts inside the same tuple.
       The normalized read counts inside each 4-tuple sum to one."""
    message = ""
    click.echo("Using file %s as reference" % (reference))
    for file in bam:
        click.echo("Reading input from file(s)  %s" % (file))
    if len(bam) < 2:
        message += ("\nError: At least two bam file locations are" +
                    " required to perform quasispecies distance comparison")
    # indicate if the start or end position is < 0 or a priori invalid
    if type(startpos) == int and int(startpos) < 0:
        message += ("\nError: Start position must be 0 or greater.")
    if type(endpos) == int and int(endpos) < 0:
        message += ("\nError: End position must be 0 or greater.")
    if (type(startpos) == int and
            type(endpos) == int and
            int(startpos) > int(endpos)):
                message += ("\nError: Start position must be less than" +
                            " or equal to end position")
    if message == "":  # if no error messages have been created
        pileups = Pileup_List.construct_array_of_pileups(bam,
                                                              reference)
        pileups = Pileup_List(pileups)

        if startpos > endpos:
            message += ("Error: Empty pileup was produced from BAM files." +
                        "Halting program")
        if message == "":
            click.echo("The start position is %d." % startpos)
            click.echo("The end position is %d." % endpos)
            click.echo("Constructed pileup from reference.")
            # click.echo the number of positions in pileup
            if truncate is not 'dont_truncate':
                click.echo("The pileup covers %d positions before truncation."
                           % len(pileups.get_pileup_length()))
            else:
                click.echo("The pileup covers %d positions.")
            # indicate whether the user-specified start and end position is out
            # of bounds (comparing to actual number of positions in pileup)
            if startpos >= pileups.get_pileup_length():
                message += ("\nError: Start position must be less than " +
                            " number of nucleotide base positions in pileup" +
                            " (%s)." % len(pileups.get_pileup_length()))
            if endpos >= pileups.get_pileup_length():
                message += ("\nError: End position must be less than length " +
                            "of nucleotide base positions in pileup" +
                            " (%s)." % len(pileups.get_pileup_length()))
            # if there is no errors so far, proceed with running program
            if startpos is not None or endpos is not None:
                if startpos is None:
                    startpos = 0
                if endpos is None:
                    endpos = pileups.get_pileup_length() - 1
                pileups.select_pileup_range(startpos, endpos)
                click.echo("The pileup covers %d positions after selecting " +
                           "range between original pileup positions %d and %d."
                           % pileups.get_pileup_length(),
                           startpos,
                           endpos)
            if normalize:
                pileups.normalize_pileup()
            if truncate is not 'dont_truncate':
                if truncate is 'truncate_ends':
                    pileups.truncate_output()
                    click.echo("Truncating positions with no coverage that " +
                               "are contiguous with the start or end " +
                               "position of the pileup only.")
                elif truncate is 'truncate_all':
                    pileups.truncate_all_output()
                    click.echo("Truncating all positions with no coverage.")
                # end if
                click.echo("The pileup covers %d positions after truncation."
                           % pileups.get_pileup_length())
                if pileups.get_pileup_length() == 0:
                    message += ("Error: Entire pileup was truncated due to " +
                                "lack of coverage. Halting program")
            # end if
            if message == "":
                dist = DistanceMatrix(pileups.get_pileup())
                if output_distance:
                    click.echo("Outputting an angular cosine distance matrix.")
                    matrix = dist.get_angular_cosine_distance_matrix(startpos,
                                                                     endpos)
                else:
                    click.echo("Outputting a cosine similarity matrix.")
                    matrix = dist.get_cosine_similarity_matrix(startpos,
                                                               endpos)
                if output:
                    output.write(dist.get_matrix_as_csv(matrix, bam))
                else:
                    # click.echo(matrix)
                    click.echo(dist.get_matrix_as_csv(matrix, bam))
                # end if
                click.echo("Complete!")
            # end if
        # end if
    # end if
    click.echo(message)
# end def

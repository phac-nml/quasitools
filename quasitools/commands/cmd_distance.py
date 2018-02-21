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
from quasitools.parsers.mapped_read_parser import parse_mapped_reads_from_bam
from quasitools.parsers.reference_parser import parse_references_from_fasta
from quasitools.distance import Distance

@click.command('distance', short_help='Calculate the evolutionary distance'
               'between viral quasispecies using cosine similarity.')
@click.argument('reference', nargs=1, required=True,
                type=click.Path(exists=True, file_okay=True, dir_okay=False))
@click.argument('bam', nargs=-1,
                type=click.Path(exists=True, file_okay=True, dir_okay=False))
@click.option('--normalize/--dontnormalize', '-n/-d', default=True,
              help="Normalize read count data so that the read counts per " +
              "4-tuple (A, C, T, G) sum to one.")
@click.option('--startpos', '-s', type=int, help="Set the start base position" +
            " of the reference to use in the distance calculation.", default=-1)
@click.option('--endpos', '-e', type=int, help="Set the end base position of " +
              "the reference to use in the distance calculation.", default=-1)
@click.option('-o', '--output', type=click.File('w'), help="Output the " +
              "quasispecies distance matrix in CSV format in a file.")
@click.pass_context
def cli(ctx, reference, bam, normalize, startpos, endpos, output):
    """This script outputs the evolutionary distance [0 - 1] between
       quasispecies, computed using the cosine similarity function.
       It takes as input multiple bam files containing reads from viral
       quasispecies and a reference file. It outputs a pairwise distance
       matrix containing the distances between each viral quasispecies. This
       can later be saved as a CSV file.

       By default the data is normalized if not specified explicitly.
       This is done dividing base read counts (A, C, T, G) inside every 4-tuple
       by the sum of the read counts inside the same tuple. The normalized
       read counts inside each 4-tuple sum to one."""
    click.echo("Using file %s as reference" % (reference))
    for file in bam:
        print("Reading input from file(s)  %s" % (file))
    if len(bam)<2:
        click.echo("Error: At least two bam file locations are"
        + " required to perform quasispecies distance comparison")
    else:
        viralDist = Distance()
        pileup_list = viralDist.construct_pileup(bam, reference)
        matrix = viralDist.get_distance_as_csv(startpos, endpos, pileup_list, bam, normalize)
        if output:
            output.write(matrix)
        else:
            click.echo(matrix)
        #end if
        print("Complete!")
    #end if

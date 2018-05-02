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

import os
from collections import defaultdict
import click
from quasitools.quality_control import QualityControl


@click.command('quality', short_help='Perform quality control on FASTQ reads.')
@click.argument('forward', required=True,
                type=click.Path(exists=True, file_okay=True, dir_okay=False))
@click.argument('reverse', required=False,
                type=click.Path(exists=True, file_okay=True, dir_okay=False))
@click.option('-o', '--output_dir', required=True,
              type=click.Path(exists=False, file_okay=False, dir_okay=True,
                              writable=True))
@click.option('-tr', '--trim_reads', is_flag=True,
              help='Iteratively trim reads based on filter values if enabled. '
                   'Remove reads which do not meet filter values if disabled.')
@click.option('-mr', '--mask_reads', is_flag=True,
              help='Mask low coverage regions in reads based on filter values')
@click.option('-mq', '--min_qual', default=30, help='Minimum quality for '
              'read to be considered if masking is enabled.')
@click.option('-lc', '--length_cutoff', default=100,
              help='Reads which fall short of the specified length '
                   'will be filtered out.')
@click.option('-sc', '--score_cutoff', default=30,
              help='Reads that have an average quality score less than the '
                   'specified score will be filtered out.')
@click.option('-me/-mn', '--median_score/--mean_score', 'score_type',
              default=True,
              help='Use either median score (default) or mean score for '
              'score cutoff value')
@click.option('-n', '--ns', is_flag=True, help='Flag to enable the '
              'filtering of n\'s')
@click.pass_context
def cli(ctx, forward, reverse, output_dir, trim_reads, mask_reads, min_qual,
        length_cutoff, score_cutoff, score_type, ns):
    """Quasitools quality performs quality control on FASTQ reads and outputs
       the filtered FASTQ reads in the specified directory."""

    os.mkdir(output_dir)
    reads = forward

    # Combine the fwd and reverse reads into one fastq file
    if reverse:
        reads = "%s/combined_reads.fastq" % output_dir
        cat_cmd = "cat %s %s > %s" % (forward, reverse, reads)
        os.system(cat_cmd)

    quality_filters = defaultdict(dict)

    if trim_reads:
        quality_filters["trimming"] = True

    if mask_reads:
        quality_filters["masking"] = True

    quality_filters["length_cutoff"] = length_cutoff

    if score_type == "median_score":
        quality_filters["median_cutoff"] = score_cutoff
    elif score_type == "mean_score":
        quality_filters["mean_cutoff"] = score_cutoff
    # end if

    if ns:
        quality_filters["ns"] = True
    else:
        quality_filters["ns"] = False

    quality_filters["minimum_quality"] = min_qual

    quality_control = QualityControl()
    quality_control.filter_reads(quality_filters)

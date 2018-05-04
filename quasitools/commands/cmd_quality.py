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

# GLOBALS

from quasitools.quality_control import TRIMMING
from quasitools.quality_control import MASKING
from quasitools.quality_control import MIN_READ_QUAL
from quasitools.quality_control import LENGTH_CUTOFF
from quasitools.quality_control import MEDIAN_CUTOFF
from quasitools.quality_control import MEAN_CUTOFF
from quasitools.quality_control import NS


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
              help='Mask low coverage regions in reads based on filter'
              ' values.')
@click.option('-rq', '--min_read_qual', default=30, help='Minimum quality for '
              'positions in read if masking is enabled.')
@click.option('-lc', '--length_cutoff', default=100,
              help='Reads which fall short of the specified length '
                   'will be filtered out.')
@click.option('-sc', '--score_cutoff', default=30,
              help='Reads that have a median or mean quality score (depending'
                   ' on the score type specified) less than the score cutoff'
                   'value will be filtered out.')
@click.option('-me/-mn', '--median/--mean', 'score_type',
              default=True,
              help='Use either median score (default) or mean score for the '
              'score cutoff value.')
@click.option('-n', '--ns', is_flag=True, help='Flag to enable the '
              'filtering of n\'s.')
@click.pass_context
def cli(ctx, forward, reverse, output_dir, trim_reads, mask_reads,
        min_read_qual, length_cutoff, score_cutoff, score_type, ns):
    """Quasitools quality performs quality control on FASTQ reads and outputs
       the filtered FASTQ reads in the specified directory."""

    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)
    reads = forward

    # Combine the fwd and reverse reads into one fastq file
    if reverse:
        reads = "%s/combined_reads.fastq" % output_dir
        cat_cmd = "cat %s %s > %s" % (forward, reverse, reads)
        os.system(cat_cmd)

    quality_filters = defaultdict(dict)

    if trim_reads:
        quality_filters[TRIMMING] = True

    if mask_reads:
        quality_filters[MASKING] = True

    quality_filters[LENGTH_CUTOFF] = length_cutoff

    if score_type == "median":
        quality_filters[MEDIAN_CUTOFF] = score_cutoff
    elif score_type == "mean":
        quality_filters[MEAN_CUTOFF] = score_cutoff
    # end if

    if ns:
        quality_filters[NS] = True
    else:
        quality_filters[NS] = False

    quality_filters[MIN_READ_QUAL] = min_read_qual

    filter_dir = "%s/filtered.fastq" % output_dir
    quality_control = QualityControl()
    quality_control.filter_reads(reads, filter_dir, quality_filters)

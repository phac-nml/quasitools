"""
Copyright Government of Canada 2017 - 2018

Written by: Camy Tran and Matthew Fogel, National Microbiology Laboratory,
            Public Health Agency of Canada

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
from quasitools.patient_analyzer import PatientAnalyzer

# GLOBALS

from quasitools.quality_control import TRIMMING
from quasitools.quality_control import MASKING
from quasitools.quality_control import MIN_READ_QUAL
from quasitools.quality_control import LENGTH_CUTOFF
from quasitools.quality_control import MEDIAN_CUTOFF
from quasitools.quality_control import MEAN_CUTOFF
from quasitools.quality_control import NS
from quasitools.patient_analyzer import ERROR_RATE
from quasitools.patient_analyzer import MIN_VARIANT_QUAL
from quasitools.patient_analyzer import MIN_DP
from quasitools.patient_analyzer import MIN_AC
from quasitools.patient_analyzer import MIN_FREQ

BASE_PATH = os.path.abspath(os.path.join(os.path.abspath(__file__),
                                         os.pardir, os.pardir, "data"))
REFERENCE = os.path.join(BASE_PATH, "hxb2_pol.fas")
GENES_FILE = os.path.join(BASE_PATH, "hxb2_pol.bed")
MUTATION_DB = os.path.join(BASE_PATH, "mutation_db.tsv")


@click.command('hydra', short_help='Identify HIV Drug Resistance in a next '
               'generation sequencing dataset.')
@click.argument('forward', required=True,
                type=click.Path(exists=True, file_okay=True, dir_okay=False))
@click.argument('reverse', required=False,
                type=click.Path(exists=True, file_okay=True, dir_okay=False))
@click.option('-o', '--output_dir', required=True,
              type=click.Path(exists=False, file_okay=False, dir_okay=True,
                              writable=True))
@click.option('-m', '--mutation_db',
              type=click.Path(exists=True, file_okay=True, dir_okay=False),
              default=MUTATION_DB)
@click.option('-rt', '--reporting_threshold', default=1,
              type=click.IntRange(1, 100, clamp=True),
              help='Minimum mutation frequency percent to report.')
@click.option('-gc', '--generate_consensus',
              help='Generate a mixed base consensus sequence.', is_flag=True)
@click.option('-cp', '--consensus_pct', default=20,
              type=click.IntRange(1, 100, clamp=True),
              help='Minimum percentage a base needs to be incorporated '
              'into the consensus sequence.')
@click.option('-q', '--quiet', is_flag=True,
              help='Suppress all normal output.')
@click.option('-tr', '--trim_reads', is_flag=True,
              help='Iteratively trim reads based on filter values if enabled. '
                   'Remove reads which do not meet filter values if disabled.')
@click.option('-mr', '--mask_reads', is_flag=True,
              help='Mask low coverage regions in reads based on filter'
              ' values. When this option and --ns are both enabled, filtering'
              ' of n\'s will be performed before masking of low coverage'
              ' regions.')
@click.option('-rq', '--min_read_qual', default=30, help='Minimum quality for '
              'a position in a read to be masked.')
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
              'filtering of n\'s. When this option and --mask_reads are both'
              'enabled, filtering of n\'s will be performed before masking'
              ' of low coverage regions.')
@click.option('-e', '--error_rate', default=0.0021,
              help='Error rate for the sequencing platform.')
@click.option('-vq', '--min_variant_qual', default=30, help='Minimum quality '
              'for variant to be considered later on in the pipeline.')
@click.option('-md', '--min_dp', default=100,
              help='Minimum required read depth for variant to be considered'
              ' later on in the pipeline.')
@click.option('-ma', '--min_ac', default=5,
              help='The minimum required allele count for variant to be '
              'considered later on in the pipeline.')
@click.option('-mf', '--min_freq', default=0.01,
              help='The minimum required frequency for mutation to be '
              'considered in drug resistance report.')
@click.option('-i', '--id',
              help='Specify FASTA sequence identifier to be used in the '
              'consensus report.')
@click.pass_context
def cli(ctx, output_dir, forward, reverse, mutation_db, reporting_threshold,
        generate_consensus, consensus_pct, quiet, trim_reads, mask_reads,
        min_read_qual, length_cutoff, score_cutoff, score_type, ns, error_rate,
        min_variant_qual, min_dp, min_ac, min_freq, id):

    os.mkdir(output_dir)
    reads = forward

    # The fasta_id is used as the sequence id in the consensus report
    # and as the RG-ID in the bt2-generated bam file.
    # It defaults to the forward fasta file name.
    if id:
        fasta_id = id
    else:
        fasta_id = os.path.basename(forward).split('.')[0]

    # Combine the fwd and reverse reads into one fastq file
    if reverse:
        reads = "%s/combined_reads.fastq" % output_dir
        cat_cmd = "cat %s %s > %s" % (forward, reverse, reads)
        os.system(cat_cmd)

        # If user did not specify an id, append name of reverse fasta file
        if not id:
            fasta_id += ("_%s" % os.path.basename(reverse).split('.')[0])

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

    variant_filters = defaultdict(dict)
    variant_filters[ERROR_RATE] = error_rate
    variant_filters[MIN_VARIANT_QUAL] = min_variant_qual
    variant_filters[MIN_DP] = min_dp
    variant_filters[MIN_AC] = min_ac
    variant_filters[MIN_FREQ] = min_freq

    patient_analyzer = PatientAnalyzer(id=REFERENCE[REFERENCE.rfind('/')+1:],
                                       output_dir=output_dir,
                                       reads=reads, reference=REFERENCE,
                                       genes_file=GENES_FILE,
                                       mutation_db=mutation_db, quiet=quiet,
                                       consensus_pct=consensus_pct)

    patient_analyzer.filter_reads(quality_filters)
    patient_analyzer.analyze_reads(fasta_id, variant_filters,
                                   reporting_threshold, generate_consensus)

"""
Copyright Government of Canada 2017

Written by: Camy Tran, National Microbiology Laboratory,
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

BASE_PATH = os.path.dirname(os.path.realpath(__file__))
REFERENCE = BASE_PATH + "/data/hxb2_pol.fas"
GENES_FILE = BASE_PATH + "/data/hxb2_pol.bed"
MUTATION_DB = BASE_PATH + "/data/mutation_db.tsv"


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
              help='minimum mutation frequency percent to report.')
@click.option('-tc', '--target_coverage', default=-1,
              help='downsample input to target coverage.'
              'Defaults to -1 to disable downsampling.')
@click.option('-gc', '--generate_consensus',
              help='Generate a mixed base consensus sequence.', is_flag=True)
@click.option('-cp', '--consensus_pct', default=20,
              type=click.IntRange(1, 100, clamp=True),
              help='minimum percentage a base needs to be incorporated'
              'into the consensus sequence.')
@click.option('-q', '--quiet', is_flag=True)
@click.option('-lc', '--length_cutoff', default=100)
@click.option('-sc', '--score_cutoff', default=30)
@click.option('-n', '--ns', is_flag=True)
@click.option('-e', '--error_rate', default=0.0021)
@click.option('-mq', '--min_qual', default=30)
@click.option('-md', '--min_dp', default=100)
@click.option('-ma', '--min_ac', default=5)
@click.option('-mf', '--min_freq', default=0.01)
@click.pass_context
def cli(ctx, output_dir, forward, reverse, mutation_db, reporting_threshold,
        target_coverage, generate_consensus, consensus_pct,
        quiet, length_cutoff, score_cutoff, ns,
        error_rate, min_qual, min_dp, min_ac, min_freq):

    os.mkdir(output_dir)
    reads = forward

    # Combine the fwd and reverse reads into one fastq file
    if reverse:
        reads = "%s/combined_reads.fastq" % output_dir
        cat_cmd = "cat %s %s > %s" % (forward, reverse, reads)
        os.system(cat_cmd)

    patient_analyzer = PatientAnalyzer(id=REFERENCE[REFERENCE.rfind('/')+1:],
                                       output_dir=output_dir,
                                       reads=reads, reference=REFERENCE,
                                       genes_file=GENES_FILE,
                                       mutation_db=mutation_db, quiet=quiet,
                                       consensus_pct=consensus_pct)

    read_filters = defaultdict(dict)
    read_filters["length_cutoff"] = length_cutoff
    read_filters["score_cutoff"] = score_cutoff
    if ns:
        read_filters["ns"] = True
    else:
        read_filters["ns"] = False

    patient_analyzer.filter_reads(read_filters)

    if target_coverage != -1:
        patient_analyzer.downsample_reads(target_coverage)

    variant_filters = defaultdict(dict)
    variant_filters["error_rate"] = error_rate
    variant_filters["min_qual"] = min_qual
    variant_filters["min_dp"] = min_dp
    variant_filters["min_ac"] = min_ac
    variant_filters["min_freq"] = min_freq

    patient_analyzer.analyze_reads(variant_filters, reporting_threshold,
                                   generate_consensus)

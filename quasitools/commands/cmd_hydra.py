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

import pdb
import click
import os
from quasitools.patient_analyzer import PatientAnalyzer
from collections import defaultdict


@click.command('hydra', short_help='Provides a pipeline for identifying drug resistance.')
@click.argument('output_dir', required=True, type=click.Path(exists=False))
@click.argument('reads', required=True, type=click.Path(exists=True))
@click.argument('mutation_db', required=False, type=click.Path(exists=True), default="quasitools/data/mutation_db.tsv")
@click.option('-r', '--reporting_threshold', default=1, type=click.IntRange(1, 100, clamp=True),
              help='minimum mutation frequency percent to report.')
@click.option('-t', '--target_coverage', default=10000,
              help='downsample input to target coverage. Set to -1 to disable downsampling.')
@click.option('-g', '--generate_consensus', help='Generate a mixed base consensus sequence.', is_flag=True)
@click.option('-cp', '--consensus_pct', default=20, type=click.IntRange(1, 100, clamp=True),
              help='Minimum percentage a base needs to be incorporated into the consensus sequence.')
@click.option('-c', '--combine', is_flag=True)
@click.option('-f', '--filter', default='_R1_|_R2_')
@click.option('-q', '--quiet', is_flag=True)
@click.option('-p', '--processes', default=8)
@click.pass_context
def cli(ctx, output_dir, reads, mutation_db, reporting_threshold,
        target_coverage, generate_consensus, consensus_pct,
        combine, filter, quiet, processes):
    click.echo("Running hydra cmd")
    
    options = {}
    options["reporting_threshold"] = reporting_threshold
    options["target_coverage"] = target_coverage
    options["generate_consensus"] = generate_consensus
    options["consensus_pct"] = consensus_pct
    options["combine"] = combine
    options["filter"] = filter
    options["quiet"] = quiet
    options["processes"] = processes
    
    options["reporting_threshold"] = reporting_threshold

    patient_analyzer = PatientAnalyzer(id=1, output_dir=output_dir,
        reads=reads, reference="quasitools/data/hxb2_pol.fas", 
        genes_file="quasitools/data/hxb2_pol.bed",
        mutation_db=mutation_db, options=options)

    filters = defaultdict(dict)
    filters["reads"]["length_cutoff"] = 100
    filters["reads"]["score_cutoff"] = 30
    filters["reads"]["ns"] = 1

    patient_analyzer.filter_reads(filters["reads"])

    patient_analyzer.downsample_reads(target_coverage)

    filters = defaultdict(dict)
    filters["variant_filtering"] = defaultdict(dict)
    filters["variant_filtering"]["error_rate"] = 0.0021
    filters["variant_filtering"]["min_qual"] = 30
    filters["variant_filtering"]["min_dp"] = 100
    filters["variant_filtering"]["min_ac"] = 5

    filters["mutations"]["min_freq"] = 0.01

    patient_analyzer.analyze_reads(filters, generate_consensus)


    




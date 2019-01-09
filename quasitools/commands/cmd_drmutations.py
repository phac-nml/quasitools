"""
Copyright Government of Canada 2017

Written by: Cole Peters, Eric Chubaty, National Microbiology Laboratory,
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

import click
from quasitools.cli import pass_context
from quasitools.aa_census import AACensus
from quasitools.aa_variant import AAVariantCollection
from quasitools.mutations import MutationDB
from quasitools.parsers.reference_parser import parse_references_from_fasta
from quasitools.parsers.mapped_read_parser import parse_mapped_reads_from_bam
from quasitools.parsers.genes_file_parser import parse_genes_file
from quasitools.parsers.nt_variant_file_parser \
    import parse_nt_variants_from_vcf


@click.command('find_mutations',
               short_help='Identifies amino acid mutations.')
@click.argument('bam', required=True,
                type=click.Path(exists=True, file_okay=True, dir_okay=False))
@click.argument('reference', required=True,
                type=click.Path(exists=True, file_okay=True, dir_okay=False))
@click.argument('variants', required=True,
                type=click.Path(exists=True, file_okay=True, dir_okay=False))
@click.argument('genes_file', required=True,
                type=click.Path(exists=True, file_okay=True, dir_okay=False))
@click.argument('mutation_db', required=True,
                type=click.Path(exists=True, file_okay=True, dir_okay=False))
@click.option('-f', '--min_freq', default=0.01,
              help='the minimum required frequency.')
@click.option('-t', '--reporting_threshold', default=1,
              help='the minimum percentage required for an entry in the drug'
              'resistant report.')
@click.option('-o', '--output', type=click.File('w'))
@pass_context
def cli(ctx, bam, reference, variants, genes_file, min_freq, mutation_db,
        reporting_threshold, output):
    rs = parse_references_from_fasta(reference)

    mapped_read_collection_arr = []
    for r in rs:
        # Create a MappedReadCollection object
        mapped_read_collection_arr.append(parse_mapped_reads_from_bam(r, bam))

    variants_obj = parse_nt_variants_from_vcf(variants, rs)

    # Mask the unconfident differences
    for mrc in mapped_read_collection_arr:
        mrc.mask_unconfident_differences(variants_obj)

    # Parse the genes from the gene file
    genes = parse_genes_file(genes_file, rs[0].name)

    # Determine which frames our genes are in
    frames = set()

    for gene in genes:
        frames.add(genes[gene]['frame'])

    # Create an AACensus object
    aa_census = AACensus(reference, mapped_read_collection_arr, genes, frames)

    # Create AAVar collection and print the aavf file
    aa_vars = AAVariantCollection.from_aacensus(aa_census)

    # Filter for mutant frequency
    aa_vars.filter('mf' + str(min_freq), 'freq<' + str(min_freq), True)

    # Build the mutation database
    mutation_db = MutationDB(mutation_db, genes)

    # Generate the mutation report
    if output:
        output.write(aa_vars.report_dr_mutations(mutation_db,
                                                 reporting_threshold))
        output.close()
    else:
        click.echo(aa_vars.report_dr_mutations(mutation_db,
                                               reporting_threshold))

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

import click
from quasitools.parsers.mutation_report_parser \
    import parse_mutations_from_hmcf


@click.command('drreport', short_help='Identifies drug resistances.')
@click.argument('hmcf', required=True,
                type=click.Path(exists=True, file_okay=True, dir_okay=False))
@click.pass_context
def cli(ctx, hmcf):
    # simply generate a simplied AAmutation object
    # it only needs position-mutationAA
    # perhaps use a map

    mutation_list = parse_mutations_from_hmcf(hmcf)

    click.echo(mutation_list)

    # create a class that will handle this work:
    # parse the drug resistance object from the xml
    # 

    # 



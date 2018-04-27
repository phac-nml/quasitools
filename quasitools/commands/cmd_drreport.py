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
import os
from quasitools.parsers.mutation_report_parser import parse_mutations_from_hmcf
from quasitools.parsers.algorithm_parser import parse_drugs_from_xml
from quasitools.drug import Drug, DrugCollection
from quasitools.evaluated_drug import EvaluatedDrug, EvaluatedDrugCollection

BASE_PATH = os.path.abspath(os.path.join(os.path.abspath(__file__),
                                         os.pardir, os.pardir, "data"))
DEFAULT_XML = os.path.join(BASE_PATH, "HIVDB_8.4.xml")


@click.command('drreport', short_help='Identifies drug resistances.')
@click.argument('hmcf', required=True,
                type=click.Path(exists=True, file_okay=True, dir_okay=False))
@click.argument('xml', required=False,
                type=click.Path(exists=True, file_okay=True, dir_okay=False))
@click.option('-o', '--output', type=click.File('w'))
@click.pass_context
def cli(ctx, hmcf, xml, output):
    # generate a simplied AAmutation object
    # it only needs position-mutationAA
    # perhaps use a map

    if xml:
        drugs = parse_drugs_from_xml(xml)
    else:
        drugs = parse_drugs_from_xml(DEFAULT_XML)

    mutation_list = parse_mutations_from_hmcf(hmcf)

    # click.echo(mutation_list)

    # create a class that will handle this work:
    # parse the drug resistance object from the xml
    # 

    # For each drug in drugs, check is the mutation_list is resistant

    # for drug in drugs.drug_list:
    #     print(drug.to_csv_entry())

    drug_resistance_list = drugs.evaluate_drug_resistances(mutation_list)

    if output:
        output.write(drug_resistance_list.to_csv_string())
    else:
        click.echo(drug_resistance_list.to_csv_string())



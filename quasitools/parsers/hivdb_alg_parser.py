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

import xml.etree.ElementTree as ET
from collections import defaultdict
# from quasitools.grammar.parser import Parser

def parse_drugs_from_xml(xml_file):
    tree = ET.parse(xml_file)

    drugs = defaultdict(lambda: defaultdict(dict))

    root = tree.getroot()

    # parse each drug: name, rule
    # rule will be scores [ {95K: 15}, {97A: 10} , etc]
    # and max scores [ [

    # {'148': { {H: 60}, {K: 60}, {N: 10} }
    for drug in root.findall('DRUG'):
        name = drug.find('NAME').text
        rule = drug.find('RULE').find('CONDITION').text
        print ("%s %s" % (name, rule))
        

    return drugs

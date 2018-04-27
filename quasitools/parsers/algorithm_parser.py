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
from quasitools.drug import Drug, DrugCollection

def parse_drugs_from_xml(xml_file):
    """Parses an algorithm in xml format and generates a DrugCollection object"""
    
    tree = ET.parse(xml_file)
    root = tree.getroot()

    drug_list = []

    # parse each drug: name, rule
    # rule will be scores [ {95K: 15}, {97A: 10} , etc]
    # and max scores [ [

    # {'148': { {H: 60}, {K: 60}, {N: 10} }
    for drug in root.findall('DRUG'):
        name = drug.find('NAME').text
        rule_string = drug.find('RULE').find('CONDITION').text.strip()

        rules = parse_rules_from_string(rule_string)

        drug_list.append(Drug(name, rules))
        # print ("%s %s" % (name, rule))

    return DrugCollection(drug_list)

def parse_rules_from_string(rule_string):
    """A helper method that parses rules for a Drug object from a string
    containing rules"""

    rules = ""

    lines = rule_string.split('\n')

    for line in lines:
        line = line.strip()
        # print("line: %s" % line)
        tokens = line.split(" ")
        # print(tokens[0])

        # handle first line, which will start with "SCORE FROM"
        if tokens[0] == "SCORE":
            line = line.split("(", 1)[1]
            line = line.split(",", 1)[0]
            # print(line)

        # handle max scores
        elif tokens[0] == "MAX":
            line = line.split("(", 1)[1]
            line = line.split(")", 1)[0]
            # print(line)

        # handle conjunctions
        elif tokens[0][0] == '(':
            line = line.split("(", 1)[1]
            line = line.split(",", 1)[0]
            # print(line)

        # handle regular rule
        else:
            line = line.split(",", 1)[0]
            # print(line)



    return rules

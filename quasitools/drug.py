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

from collections import defaultdict
from quasitools.evaluated_drug import EvaluatedDrug, EvaluatedDrugCollection


class Drug(object):
    def __init__(self, name, rules):
        self.name = name
        self.rules = rules

    def calculate_score(self, mutations):
        # go through each mutation 
        # go through each rule in rules?

        # sum the score
        # Susceptible: Total score 0 to 9
        # Potential low-level resistance: Total score 10 to 14
        # Low-level resistance: Total score 15 to 29
        # Intermediate resistance: Total score 30 to 59
        # High-level resistance: Total score >= 60
        total_score = 0

        # for pos in mutations:
        #   if pos in self.rules:
        #         for mutant in mutations[pos]:
        #             if mutations[pos][mutant]
    
        for pos in self.rules:
            if pos in mutations:
                for mutant in self.rules[pos]:
                    if mutant in mutations[pos]:
                        total_score += self.rules[pos][mutant]

        return total_score

    # temporary, for testing
    def to_csv_entry(self):
            return "%s,%s\n" % (
            self.name, self.rules
        )

    @classmethod
    def determine_resistance_level(cls, score):
        if score < 0:
            resistance = "Not resistant"
        elif score >= 0 and score <= 9:
            resistance = "Susceptible"
        elif score >= 10 and score <= 14:
            resistance = "Potential low-level resistance"
        elif score >= 15 and score <= 29:
            resistance = "Low-level resistance"
        elif score >= 30 and score <= 59:
            resistance ="Intermediate resistance"
        else:
            resistance = "High-level resistance"
        
        return resistance


class DrugCollection(object):
    """Collection of drugs that are parsed from an xml-format algorithm"""
    
    def __init__(self, drug_list):
        self.drug_list = drug_list

    def evaluate_drug_resistances(self, mutations):
        evaluated_drug_list = []

        # for mutation in mutations:
        for drug in self.drug_list:
            total_score = drug.calculate_score(mutations)
            evaluated_drug_list.append(EvaluatedDrug(drug.name,
                                                     total_score,
                                                     Drug.determine_resistance_level(total_score)))

        evaluated_drugs = EvaluatedDrugCollection(evaluated_drug_list)

        return evaluated_drugs


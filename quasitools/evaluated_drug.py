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


class EvaluatedDrug(object):
    def __init__(self, name, total_score, resistance_level):
        self.name = name
        self.total_score = total_score
        self.resistance_level = resistance_level
    
    def to_csv_entry(self):
        return "%s,%s,%s\n" % (
            self.name, self.total_score, self.resistance_level
        )
        

"""
    Collection of drugs for which a mutant is resistant to.
"""
class EvaluatedDrugCollection(object):
    def __init__(self, drug_list):
        self.drug_list = drug_list

    def to_csv_string(self):
        """"Build a string representation of our EvaluatedDrugCollection
        objects (i.e. a csv file)."""

        report = ("#Drug Name,Total Score,Drug Resistance Level\n")

        for drug in self.drug_list:
            report += drug.to_csv_entry()

        return report
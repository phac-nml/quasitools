"""
Copyright Government of Canada 2017

Written by: Cole Peters, National Microbiology Laboratory, Public Health Agency of Canada

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


class Mutation(object):
    def __init__(self, gene, wildtype, gene_pos, frequency, filter):
        self.gene = gene
        self.wildtype = wildtype
        self.gene_pos = gene_pos
        self.frequency = frequency
        self.filter = filter


class MutationDBEntry(object):
    def __init__(self, gene, wildtype, gene_pos, category, surveillance,
                 comment):
        self.gene = gene
        self.wildtype = wildtype
        self.gene_pos = gene_pos
        self.category = category
        self.surveillance = surveillance
        self.comment = comment


class MutationDB(object):
    def __init__(self, mutation_db, genes):
        self.mutation_db = mutation_db
        self.genes = genes
        self.mutations = defaultdict(dict)

        self._build()

    def _build(self):
        """Builds the mutations dictionary for the mutation database object"""

        with open(self.mutation_db, "r") as input:
            for line in input:
                if line[0] != "#":
                    gene, wildtype, gene_pos, mutation, category, \
                        surveillance, comment = line.rstrip().split("\t")

                    ref_pos = \
                        int(gene_pos) - 1 + (self.genes[gene]["start"] // 3)

                    if surveillance == "Yes" or surveillance == "No":
                        db_entry = MutationDBEntry(gene, wildtype, gene_pos,
                                                   category, surveillance,
                                                   comment)

                        self.mutations[ref_pos][mutation] = db_entry
                    else:
                        raise ValueError("Mutation Database is incorrectly"
                                         "formatted.")

    def positions(self):
        """Returns a list of positions that contain a drug resistant mutation
        """

        return sorted(self.mutations.keys())

    def mutations_at(self, pos):
        """Returns a list of drug resistant mutations found at the given
        position.
        """

        return self.mutations[pos]

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

import pytest
import os
from quasitools.mutations import MutationDB
from quasitools.parsers.genes_file_parser import parse_BED4_file


TEST_PATH = os.path.dirname(os.path.abspath(__file__))


class TestMutationDB:
    @classmethod
    def setup_class(self):
        BED4_file = TEST_PATH + "/data/hxb2_pol.bed"
        mutation_db_file = TEST_PATH + "/data/mutation_db.tsv"

        genes = parse_BED4_file(BED4_file, "hxb2_pol")

        self.mutation_db = MutationDB(mutation_db_file, genes)

    def test_mutations_at(self):

        # Our pos is one less than the mutation db file since our first
        # position is 0 but the files is 1
        pos = 9
        mutations = ["F", "I", "V", "R", "Y"]

        mutations_at_pos = self.mutation_db.mutations_at(pos)

        assert len(mutations_at_pos) == 5

        for mutation in mutations:
            assert mutation in mutations_at_pos

    def test_positions(self):

        positions = self.mutation_db.positions()

        assert len(positions) == 96

        assert positions[0] == 9

        assert positions[30] == 142

        assert positions[-1] == 921

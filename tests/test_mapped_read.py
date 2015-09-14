"""
Copyright Government of Canada 2015

Written by: Eric Enns, National Microbiology Laboratory, Public Health Agency of Canada

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
from quasitools.mapped_read import MappedRead

class TestMappedRead:
    @classmethod
    def setup_class(self):
        self.mr = MappedRead('read1', 10, 300, {'110': 'g', '300': 'tagt'}, 100, 390, 97, 99.3127, '+')

    def test_query_length(self):
        assert self.mr.query_length() == 291

    def test_codon_start(self):
        assert self.mr.codon_start(0) == 102
        assert self.mr.codon_start(1) == 100
        assert self.mr.codon_start(2) == 101

    def test_codon_end(self):
        assert self.mr.codon_end(0) == 389
        assert self.mr.codon_end(1) == 390
        assert self.mr.codon_end(2) == 388

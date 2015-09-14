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
from quasitools.references import References

class TestReferences:
    def test_from_fasta(self):
        references = References.from_fasta('tests/data/ref1.fasta')

        assert list(references.references) == ['ref1']
        assert references.references['ref1'].seq == 'AGCATGTTAGATAAGATAGCTGTGCTAGTAGGCAGTCAGCGCCAT'

    def test_sub_seq(self):
        references = References.from_fasta('tests/data/ref1.fasta')

        assert references.sub_seq('ref1',1,5) == 'GCATG'

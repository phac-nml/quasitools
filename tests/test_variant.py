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
from quasitools.variant import Variant

class TestVariant:
    @classmethod
    def setup_class(self):
        self.variant = Variant('hxb2_pol',1,ref='g',alt='t',qual='30',filter='PASS',info={'DP':400,'AC':12,'AF':0.03})

    def test_str(self):
        assert str(self.variant) == 'hxb2_pol\t1\t.\tg\tt\t30\tPASS\tDP=400;AC=12;AF=0.0300'

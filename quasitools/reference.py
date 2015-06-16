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

class Reference(object):
    def __init__(self, name, seq):
        self.name = name
        self.seq = seq

    def sub_seq(self, start, end):
        """Returns a portion of the sequence

        >>> r = Reference('gattaca', 'gattaca')
        >>> print(r.sub_seq(1,5))
        attac
        >>> print(r.sub_seq(1,1))
        a
        """
        return self.seq[start:end+1]

if __name__ == '__main__':
    import doctest
    doctest.testmod()

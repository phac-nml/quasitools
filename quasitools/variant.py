"""
Copyright Government of Canada 2015-2017

Written by: Eric Enns, National Microbiology Laboratory,
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


class Variant(object):

    def __init__(self, chrom, pos, id='.', ref='', alt='', qual='.',
                 filter='.', info='.'):
        self.chrom = chrom
        self.pos = pos
        self.id = id
        self.ref = ref
        self.alt = alt
        self.qual = qual
        self.filter = filter
        self.info = info


class VariantCollection(object):

    def __init__(self, references):
        self.variants = defaultdict(lambda: defaultdict(dict))
        self.references = references
        self.filters = {}

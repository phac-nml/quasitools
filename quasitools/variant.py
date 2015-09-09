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

class Variant(object):
    def __init__(self, chrom, pos, id='.', ref='', alt='', qual='.', filter='.', info='.'):
        self.chrom  = chrom
        self.pos    = pos
        self.id     = id
        self.ref    = ref
        self.alt    = alt
        self.qual   = qual
        self.filter = filter
        self.info   = info

    def __info_to_str(self):
        """Convert info dict to info string for vcf entry."""
        return "DP=%i;AC=%i;AF=%0.4f" % (self.info['DP'],self.info['AC'],self.info['AF'])

    def __str__(self):
        """Build a string representation of our Variant object (i.e. an entry in a vcf file)."""
        return "%s\t%i\t%s\t%s\t%s\t%s\t%s\t%s" % (self.chrom, self.pos, self.id, self.ref, self.alt, self.qual, self.filter, self.__info_to_str())

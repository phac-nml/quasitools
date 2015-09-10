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

import re
from collections import defaultdict
from datetime import date
from scipy.stats import poisson
from numpy import log10
from quasitools.mapped_reads import MappedReads
from quasitools.variant import Variant

class Variants(object):
    def __init__(self, error_rate, references):
        self.variants = defaultdict(lambda: defaultdict(dict) )
        self.error_rate = error_rate
        self.references =  references
        self.filters = {}

    @classmethod
    def from_bam(cls, references, error_rate, overlap_cutoff, identity_cutoff, bam):
        """Build the Variants object from a bam."""

        mapped_reads_arr = []

        for rid, r in references.references.items():
            #create MappedReads object
            mapped_reads_arr.append(MappedReads.from_bam(r, overlap_cutoff, identity_cutoff, bam))

        obj = cls.from_mapped_reads(error_rate, references, *mapped_reads_arr)

        return obj

    @classmethod
    def from_mapped_reads(cls, error_rate, references, *args):
        """Build the Variants object from a MappedReads object"""
        obj = cls(error_rate, references)

        for mapped_reads in args:
            pileup = mapped_reads.pileup()
            rid = mapped_reads.reference.name

            #do not include indels in coverage calculations
            coverage = [sum([v if not k.startswith('+') and not k.startswith('-') else 0 for k,v in pileup[pos].items()])for pos in range(0,len(pileup))]

            for pos in range(0,len(pileup)):
                for event, event_count in pileup[pos].items():
                    alt_allele = event.lower()
                    if len(event) > 1:
                        alt_allele = event[:1].lower()

                    if alt_allele != '-' and alt_allele != mapped_reads.reference.sub_seq(pos,pos).lower():
                        if rid in obj.variants and pos+1 in obj.variants[rid] and alt_allele in obj.variants[rid][pos+1]:
                            obj.variants[rid][pos+1][alt_allele].info['AC'] += event_count
                            obj.variants[rid][pos+1][alt_allele].info['AF'] = obj.variants[rid][pos+1][alt_allele].info['AC'] / coverage[pos]
                        else:
                            variant = Variant(chrom=mapped_reads.reference.name, pos=pos+1, ref=mapped_reads.reference.sub_seq(pos,pos).lower(), alt=alt_allele, info={'DP':coverage[pos],'AC':event_count,'AF':event_count / coverage[pos]})
                            obj.variants[rid][pos+1][alt_allele] = variant

                for alt_allele, variant in obj.variants[rid][pos+1].items():
                    variant.qual = obj.__calculate_variant_qual(variant.info['AC']-1, variant.info['DP'])

        return obj

    def __str__(self):
        """Build a string representation of our Variants object (i.e. a vcf file)."""
        d = date.today()

        report = "##fileformat=VCFv4.2\n";
        report += "##fileDate=%s\n" % (d.strftime("%Y%m%d"));
        report += "##source=quasitools\n";

        #print contig info per reference
        for rid, reference in self.references.references.items():
            report += "##contig=<ID=%s,length=%i>\n" % (reference.name, len(reference.seq))

        report += "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">\n"
        report += "##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele Count\">\n"
        report += "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency\">\n"

        for id, filter in self.filters.items():
            report += "##FILTER=<ID=%s,Description=\"Set if %s; %s\">\n" % (id,filter['result'],filter['expression'])

        report += "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"

        for rid in self.variants:
            for pos in self.variants[rid]:
                for alt_allele, variant in sorted(self.variants[rid][pos].items()):
                    if variant.qual > 0:
                        report += "\n" + str(variant)

        return report

    def __calculate_variant_qual(self, variant_count, coverage):
        """Calculate variant qual using poisson distribution."""
        avg_errors = coverage * self.error_rate

        prob = poisson.cdf(variant_count, avg_errors)

        qual = 100
        if prob < 1:
            qual = int(min((-10) * log10(1-prob), 100))

        return qual

    def filter(self, id, expression, result):
        """Apply filter to variants given an id, expression and result."""
        self.filters[id] = {'expression':expression,'result':result}

        #only allow simple expressions for the time being i.e. DP>30
        (attribute, operator, value) = re.split('([><=!]+)', expression)

        for rid in self.variants:
            for pos in self.variants[rid]:
                for alt_allele, variant in self.variants[rid][pos].items():
                    attribute_value = None
                    if hasattr(variant, attribute.lower()):
                        attribute_value = eval("variant.%s" % attribute.lower())
                    else:
                        attribute_value = variant.info[attribute.upper()]

                    if eval("%s %s %s" % (attribute_value,operator,value)) != result:
                        if variant.filter == '.':
                            variant.filter = 'PASS'
                    else:
                        if variant.filter == '.' or variant.filter == 'PASS':
                            variant.filter = id
                        else:
                            variant.filter += ";%s" % id

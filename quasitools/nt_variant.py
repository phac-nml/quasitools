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
from datetime import date
from scipy.stats import poisson
from numpy import log10
from quasitools.mapped_read import MappedReadCollection
from quasitools.variant import Variant, VariantCollection

class NTVariant(Variant):
    def __info_to_str(self):
        """Convert info dict to info string for vcf entry."""
        return "DP=%i;AC=%i;AF=%0.4f" % (self.info['DP'],
                                         self.info['AC'],
                                         self.info['AF'])

    def to_vcf_entry(self):
        return "%s\t%i\t%s\t%s\t%s\t%s\t%s\t%s" % (self.chrom, self.pos,
                                                   self.id, self.ref, self.alt,
                                                   self.qual, self.filter,
                                                   self.__info_to_str())

class NTVariantCollection(VariantCollection):
    @classmethod
    def from_mapped_read_collections(cls, error_rate, references, *args):
        """Build the NTVariantCollection from any number of MappedReadCollection objects"""
        obj = cls(references)

        for mapped_read_collection in args:
            pileup = mapped_read_collection.pileup()
            rid = mapped_read_collection.reference.name

            #do not include indels in coverage calculations
            coverage = [sum([v if not k.startswith('+') and not k.startswith('-') else 0 for k,v in pileup[pos].items()]) for pos in range(0,len(pileup))]

            for pos in range(0,len(pileup)):
                for event, event_count in pileup[pos].items():
                    alt_allele = event.lower()
                    if len(event) > 1:
                        alt_allele = event[:1].lower()

                    if alt_allele != '-' and alt_allele != mapped_read_collection.reference.sub_seq(pos,pos).lower():
                        if rid in obj.variants and pos+1 in obj.variants[rid] and alt_allele in obj.variants[rid][pos+1]:
                            obj.variants[rid][pos+1][alt_allele].info['AC'] += event_count
                            obj.variants[rid][pos+1][alt_allele].info['AF'] = obj.variants[rid][pos+1][alt_allele].info['AC'] / coverage[pos]
                        else:
                            variant = NTVariant(chrom=mapped_read_collection.reference.name, pos=pos+1, ref=mapped_read_collection.reference.sub_seq(pos,pos).lower(), alt=alt_allele, info={'DP':coverage[pos],'AC':event_count,'AF':event_count / coverage[pos]})
                            obj.variants[rid][pos+1][alt_allele] = variant

                for alt_allele, variant in obj.variants[rid][pos+1].items():
                    variant.qual = obj.__calculate_variant_qual(error_rate, variant.info['AC'], variant.info['DP'])

        return obj

    def to_vcf_file(self):
        """Build a string representation of our Variants object (i.e. a vcf file)."""
        d = date.today()

        report = "##fileformat=VCFv4.2\n";
        report += "##fileDate=%s\n" % (d.strftime("%Y%m%d"));
        report += "##source=quasitools\n";

        #print contig info per reference
        for reference in self.references:
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
                        report += "\n" + variant.to_vcf_entry()

        return report

    def __calculate_variant_qual(self, error_rate, variant_count, coverage):
        """Calculate variant qual using poisson distribution."""
        avg_errors = coverage * error_rate

        prob = poisson.cdf(variant_count-1, avg_errors)

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

"""
Copyright Government of Canada 2015-2017

Written by: Camy Tran, Eric Enns, National Microbiology Laboratory,
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

from quasitools.nt_variant import NTVariant, NTVariantCollection


def parse_nt_variants_from_vcf(variant_file, references):
    """Build variants object from a vcf file"""

    obj = NTVariantCollection(references)

    # Read in and parse the variants file
    # The file uses 1 as the first position but 0 is the first position in
    # mapped reads.
    with open(variant_file, "r") as input:
        for line in input:
            if line[0] != "#":
                chrom, pos, var_id, ref, alt, qual, var_filter, info = \
                    line.rstrip().split("\t")

                dp, ac, af = info.rstrip().split(';')

                dp = dp.rstrip().split('=')[1]
                ac = ac.rstrip().split('=')[1]
                af = af.rstrip().split('=')[1]

                variant_obj = NTVariant(chrom=chrom,
                                        id=var_id,
                                        pos=int(pos),
                                        ref=ref,
                                        alt=alt,
                                        qual=qual,
                                        filter=var_filter,
                                        info={
                                            'DP': dp,
                                            'AC': ac,
                                            'AF': af
                                        })

                obj.variants[chrom][int(pos)][alt] = variant_obj

    return obj

"""
Copyright Government of Canada 2017

Written by: Camy Tran, National Microbiology Laboratory,
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

import pytest
import os
from quasitools.parsers.nt_variant_file_parser import parse_nt_variants_from_vcf
from quasitools.parsers.mapped_read_parser import parse_mapped_reads_from_bam
from quasitools.parsers.reference_parser import parse_references_from_fasta
from quasitools.nt_variant import NTVariant, NTVariantCollection


TEST_PATH = os.path.dirname(os.path.abspath(__file__))


class TestNtVariantParser:
    def test_valid_vcf_file(self):
        """Tests to ensure that valid vcf files are parsed properly."""

        reference = TEST_PATH + \
            "/data/hxb2_pol.fas"
        bam = TEST_PATH + "/data/align.bam"

        rs = parse_references_from_fasta(reference)

        mapped_read_collection_arr = []
        for r in rs:
            # Create a MappedReadCollection object
            mapped_read_collection_arr.append(
                parse_mapped_reads_from_bam(r, bam))

        variants_obj = NTVariantCollection(rs)

        for i in range(0, 20):
            variant = NTVariant(chrom="hxb2_pol",
                                pos=i,
                                id=".",
                                ref='a',
                                alt='t',
                                qual="50",
                                filter="PASS",
                                info={"DP": "300",
                                      "AC": "1",
                                      "AF": "0.0025"}
                                )

            variants_obj.variants["hxb2_pol"][i]['t'] = variant

        # Create a valid vcf file
        valid_vcf_file = TEST_PATH + "/data/valid_vcf_file.vcf"

        with open(valid_vcf_file, "w+") as f:
            f.write("##fileformat=VCFv4.2\n"
                    "##fileDate=20171005\n"
                    "##source=quasitools\n"
                    "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">\n"
                    "##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele Count\">\n"
                    "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency\">\n"
                    "##FILTER=<ID=q30,Description=\"Quality below 30\">\n"
                    "##FILTER=<ID=dp100,Description=\"Read depth below 100\">\n"
                    "##FILTER=<ID=ac5,Description=\"Allele count below 5\">\n"
                    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
                    )

            for rid in variants_obj.variants:
                for pos in variants_obj.variants[rid]:
                    for alt in variants_obj.variants[rid][pos]:
                        variant = variants_obj.variants[rid][pos][alt]
                        f.write("\n%s\t%i\t%s\t%s\t%s\t%s\t%s" % (
                            variant.chrom, int(variant.pos),
                            variant.id, variant.ref,
                            variant.alt, variant.qual,
                            variant.filter))
                        f.write("\tDP=%i;AC=%i;AF=%0.4f" % (int(variant.info["DP"]),
                                                            int(variant.info["AC"]),
                                                            float(variant.info["AF"])))

        parsed_nt_var = parse_nt_variants_from_vcf(valid_vcf_file, rs)

        # Check equality of parsed NTVariantCollection vs. the valid NTVariantCollection
        for rid in parsed_nt_var.variants:
            for pos in parsed_nt_var.variants[rid]:
                for alt in parsed_nt_var.variants[rid][pos]:
                    parsed_variant = parsed_nt_var.variants[rid][pos][alt]
                    variant = variants_obj.variants[rid][pos][alt]

                    assert parsed_variant.chrom == variant.chrom
                    assert parsed_variant.pos == variant.pos
                    assert parsed_variant.id == variant.id
                    assert parsed_variant.ref == variant.ref
                    assert parsed_variant.alt == variant.alt
                    assert parsed_variant.qual == variant.qual
                    assert parsed_variant.filter == variant.filter
                    assert parsed_variant.info["DP"] == variant.info["DP"]
                    assert parsed_variant.info["AC"] == variant.info["AC"]
                    assert parsed_variant.info["AF"] == variant.info["AF"]

        os.remove(valid_vcf_file)

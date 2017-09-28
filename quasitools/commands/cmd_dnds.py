"""
Copyright Government of Canada 2017

Written by: Camy Tran, National Microbiology Laboratory, Public Health Agency of Canada

Licensed under the Apache License, Version 2.0 (the "License"); you may not use
this work except in compliance with the License. You may obtain a copy of the
License at:

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software distributed
under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR
CONDITIONS OF ANY KIND, either express or implied. See the License for the
specific language governing permissions and limitations under the License.
"""

import pdb
import click
from collections import defaultdict
from Bio.Seq import Seq
from numpy import log
from quasitools.parsers.reference_parser import parse_references_from_fasta
from quasitools.parsers.genes_file_parser import parse_genes_file

@click.command('dnds', short_help='Calculate the dn/ds value for each region in a bed file.')
@click.argument('csv', required=True, type=click.Path(exists=True))
@click.argument('reference', required=True, type=click.Path(exists=True))
@click.argument('offset', required=True, type=int)
@click.argument('output', required=True)

@click.pass_context
def cli(ctx, csv, reference, offset, output):
    click.echo("Running dnds command...")

    rs = parse_references_from_fasta(reference)
    ref_seq = rs[0].seq
    # Parse the genes - maybe move this to a separate parser file?

    genes = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: defaultdict(int))))
    with open(csv, "r") as f:
        for line in f:
            if line[0] != "#":
                gene, gene_start_end, nt_start, nt_end, ref_codon, mutant_codon, \
                    ref_aa, mutant_aa, coverage, mutant_freq, mutant_type, ns_count, s_count  = \
                    line.rstrip().split(",")
                gene_start, gene_end = gene_start_end.split('-') 
                
                gene_start = int(gene_start)
                gene_end = int(gene_end)
                nt_start = int(nt_start)
                nt_end = int(nt_end)
                mutant_freq = float(mutant_freq)
                ns_count = float(ns_count)
                s_count = float(s_count)

                genes[gene]['start'] = gene_start
                genes[gene]['end'] = gene_end

                if ns_count > 0:
                    genes[gene][nt_start-gene_start]['NS'][ns_count] += mutant_freq/100.0
                if s_count > 0:
                    genes[gene][nt_start-gene_start]['S'][s_count] += mutant_freq/100.0
    
    f.close()

    report = "#gene,pn,ps,pn_sites,ps_sites,dn/ds\n"     

    for gene in genes:
       # pdb.set_trace()
        s_sites = 0
        ns_sites = 0
        gene_seq = ref_seq[ (genes[gene]['start'] - offset): (genes[gene]['end'] - offset + 1) ]

        pn = 0
        ps = 0

        ns_ncod = 0
        s_ncod = 0

        pn_ncod = 0
        ps_ncod = 0

        for i in range(0, len(gene_seq)-1, 3): # step through this and hydra to see if boundaries are right
            codon = gene_seq[i:i+3]
            aa = Seq(codon).translate()[0]
            stop = 0
            non_syn = 0

            # synonymous sites only occur at first and third position in a codon
            for j in range(0, 3):
                for nt in ('a', 'c', 'g', 't'):
                    if nt.lower() != codon[j:j+1].lower():
                        mod_codon = codon[:j] + nt + codon[j+1:]
                        mod_aa = Seq(mod_codon).translate()[0]
                        if mod_aa.upper() != aa.upper():
                            non_syn += 1/3.0
 
            ns_sites += non_syn
            s_sites += 3-non_syn

            if non_syn > 0:
                ns_ncod += 1
            if 3-non_syn > 0:
                s_ncod += 1

            pni = 0
            psi = 0
            
            if 'NS' in genes[gene][i]:  #hasattr(genes[gene][i], 'NS'):
                for count in genes[gene][i]['NS']:
                    pni += genes[gene][i]['NS'][int(count)] * count/non_syn
                pn += pni
                pn_ncod += 1

            if 'S' in genes[gene][i]: #hasattr(genes[gene][i], 'S'):
                for count in genes[gene][i]['S']:
                    psi += genes[gene][i]['S'][int(count)] * count/(3-non_syn)
                ps += psi
                ps_ncod += 1
        if pn_ncod > 0 and ps_ncod > 0:
            pn = pn/pn_ncod
            ps = ps/ps_ncod

            if (1-(4*pn/3.0) > 0) and (1-(4*ps/3.0) > 0):
                dn = -(3/4.0)*log(1-(4*pn/3.0))
                ds = -(3/4.0)*log(1-(4*ps/3.0))
                report += "%s,%0.4f,%0.4f,%i,%i,%0.4f\n" % \
                          (gene, pn, ps, pn_ncod, ps_ncod, dn/ds)
            else:
                report += "%s,%0.4f,%0.4f,%i,%i,N/A\n" % \
                          (gene, pn, ps, pn_ncod, ps_ncod)
    
    output_file = open(output, "w")
    output_file.write(report)
    output_file.close()






    








    

#! /usr/bin/env python3

# The GPL v3 License

#    Copyright (C) 2016 University of Geneva.
#    #
#    # Author: Andrew Brown <andrew.brown@unige.ch>

#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.

#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.

#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.


import argparse
import gzip
import pysam
from os import remove
from os.path import exists
from subprocess import call
import sys
import shutil

def intersection(lst1, lst2): 
    return [item for item in lst1 if item in lst2]

def get_snp(snp, samples, vcf_file):
    chrom, pos = snp.split('_')[:2]
    tabix_fetch = vcf_file.fetch(chrom, int(pos) - 1, int(pos))
    for snp_line in tabix_fetch:
        if snp_line.id == snp:
            if 'DS' in snp_line.format:
                return (snp + '\t' + '\t'.join(str(snp_line.samples[z]['DS'])
                                               for z in samples)
                        + '\n')
            elif 'GT' in snp_line.format:
                snp_line0 = [snp_line.samples[z]['GT'][0] for z in samples]
                snp_line1 = [snp_line.samples[z]['GT'][1] for z in samples]
                snp_line = [snp_line.samples[z]['GT'] for z in samples]
                mean_dosage = str(sum(z[0] + z[1] for z in snp_line
                                      if not z[0] is None) /
                                  sum(
                                      1 for z in samples
                                      if not z[0] is None))
                genos = [z[0] + z[1]
                                      if not z[0] is None
                                      else mean_dosage for z in snp_line]
                return(snp + '\t' + '\t'.join(map(str, genos)) + '\n')
#                return (snp + '\t' + '\t'.join(str(genos) + '\n')
#                return (snp + '\t' +
#                        '\t'.join(snp_line0+snp_line1)
#                                  if not z[0] is None
#                                  else mean_dosage for z in expr_ids)
#                        +'\n')
            else:
                print("Both GT and DS fields are missing.", file=sys.stderr)
                sys.exit(1)
    print("Can't find SNP %s." % snp, file=sys.stderr)
    sys.exit(1)


def call_fastqtl(args, tss, cov_temp, temp):
    if not args.exon:
        fast_cmd = [args.fastQTL,
              '--vcf', args.vcf,
              '--bed', args.bed,
              '--cov', cov_temp,
              '-W', str(args.W),
              '--region', '%s:%s-%s' % (tss[0], tss[1], tss[1]),
              '--permute', str(args.perm),
              '--seed', str(args.seed),
              '--include-samples', str(args.include_samples),
              '--exclude-phenotypes', str('/gpfs/data/gtex-group/mpavia/methylation/data/support_files/problematic.cpgs.txt'),
              '--maf-threshold', str(0.01),
              '--out', temp]
    else:
        fast_cmd = [args.fastQTL,
              '--vcf', args.vcf,
              '--bed', args.bed,
              '--cov', cov_temp,
              '--region', '%s:%s-%s' % (tss[0], tss[1], tss[1]),
              '--permute', str(args.perm),
              '--seed', str(args.seed),
              '--grp', args.exon,
              '--out', temp]
    if args.normal:
        fast_cmd.extend("--normal")
    call(fast_cmd)

def run_analysis(args):
    if not args.exon:
        snp_position = 6
        beta_position = 16
    else:
        snp_position = 7
        beta_position = 17

    if args.genes == 0:
        forward_file = args.tmpdir+"/"+args.forward
        backward_file = args.tmpdir+"/"+args.backward
    else:
        forward_file = "%s.%s" % (args.tmpdir+"/"+args.forward, args.jobnumber)
        backward_file = "%s.%s" % (args.tmpdir+"/"+args.backward, args.jobnumber)

    cov_temp = forward_file + ".cov_temp"
    temp = forward_file + ".temp"

    if not args.overwrite:
        for x in [forward_file, backward_file, cov_temp, temp]:
            if exists(x):
                print("Running analysis would overwrite file %s." % x,
                      file=sys.stderr)
                sys.exit(1)

    genes = {}
    genes_tss = {}

    vcf_file = pysam.VariantFile(args.vcf)

    snp_ids = list(vcf_file.header.samples)

    if not args.include_samples:
        print("Provide sample list file",
                      file=sys.stderr)
        sys.exit(1)
    else:
        with open(args.include_samples) as f:
            samples = f.read().splitlines()

    with gzip.open(args.bed, 'rt') as f:
        expr_ids = f.readline().strip().split()[4:]
        for line in (line.strip().split() for line in f):
            genes_tss[line[3]] = (line[0], line[2])
    samples = intersection(samples, intersection(snp_ids, expr_ids))

    if sum(1 for x in expr_ids if x not in snp_ids) != 0:
        print("Individuals with expression do not have genotype data.", file=sys.stderr)
        print("These individuals are: %s." %
              ", ".join(x for x in expr_ids if x not in snp_ids), file = sys.stderr)
        #sys.exit(1)

    if args.cov:
        with open(args.cov) as f:
            cov_ids = f.readline().strip().split()
            if sum(1 for x in expr_ids if x not in cov_ids[1: ]) != 0:
                print("Individuals with expression do not have covariate data.",
                      file=sys.stderr)
                print("These individuals are: %s." %
                      ", ".join(x for x in expr_ids if x not in cov_ids[1: ]),
                      file = sys.stderr)
            #    sys.exit(1)
            samples = intersection(samples, cov_ids[1: ])
            places = [0]
            places.extend([cov_ids.index(x) for x in samples])
            cov_iter = (line.strip().split() for line in f)
            cov_iter = ('\t'.join(x[y] for y in places) for x in cov_iter)
            cov_string = 'id\t' + '\t'.join(samples) + '\n'
            cov_string += '\n'.join(cov_iter) + '\n'

    if args.threshold == 1:
        with open(args.results) as f:
            threshold = max(float(x.strip().split()[beta_position]) for x in f)
            print(threshold)
    else:
        threshold = args.threshold

    with open(args.results) as eqtl_file:
        if args.genes == 0:
            results_iter = (x.strip().split() for x in eqtl_file)
        else:
            results_iter = (y.strip().split() for x, y in enumerate(eqtl_file)
                    if x >= (args.jobnumber - 1) * args.genes and
                    x < args.jobnumber * args.genes)
        with open(forward_file, 'w') as results_file:
            for line in results_iter:
                results_file.write('\t'.join(line) + "\t1\n")
                genes[line[0]] = {line[snp_position]:
                                  (get_snp(line[snp_position], samples, vcf_file), 1)}
                discovered = True
                iteration = 1
                with open(cov_temp, 'wt') as cov_file:
                    _ = cov_file.write(cov_string)
                    _ = cov_file.write(
                        genes[line[0]][line[snp_position]][0] + '\n')
                    cov_file.flush()
                    while discovered:
                        discovered = False
                        iteration += 1
                        call_fastqtl(args, genes_tss[line[0]], cov_temp, temp)
                        with open(temp) as tissue_file:
                            for tissue_line in (line.strip().split()
                                                for line in tissue_file):
                                if tissue_line[0] == line[0] and float(tissue_line[beta_position]) <= threshold:
                                    snp = get_snp(
                                        tissue_line[snp_position], samples, vcf_file)
                                    genes[line[0]][tissue_line[snp_position]] = (snp, iteration)
                                    _ = cov_file.write(snp + '\n')
                                    cov_file.flush()
                                    #destination=args.tmpdir+"/Covariates.fwd."+line[0]+"."+line[snp_position]+".txt"
                                    #destination=args.tmpdir+"/Covariates.fwd."+line[0]+".txt"
                                    destination=args.tmpdir+"/Covariates."+line[0]+"."+args.forward
                                    dest = shutil.copyfile(cov_temp, destination)
                                    _ = results_file.write('\t'.join(tissue_line) +
                                                '\t%d\n' % iteration)
                                    discovered = True


    print("Backward starting", file=sys.stderr)
    with open(backward_file, 'w') as beqtl_file:
        with open(forward_file) as feqtl_file:
            for line in (line.strip().split() for line in feqtl_file):
                gene = genes[line[0]]
                if len(gene) == gene[line[snp_position]][1]:
                    _ = beqtl_file.write('\t'.join(line) + '\n')
                else:
                    with open(cov_temp, 'wt') as cov_file:
                        _ = cov_file.write(cov_string)
                        _ = cov_file.write('\n'.join(gene[x][0]
                                                     for x in gene
                                                     if x != line[snp_position])
                                           + '\n')
                        cov_file.flush()
                        call_fastqtl(args, genes_tss[line[0]], cov_temp, temp)
                        with open(temp) as tissue_file:
                            for tissue_line in (line.strip().split()
                                                for line in tissue_file):
                                if tissue_line[0] == line[0] and float(tissue_line[beta_position]) <= threshold:
                                    _ = beqtl_file.write('\t'.join(tissue_line)
                                                           + '\t%d\n' % gene[line[snp_position]][1])


    if exists(cov_temp):
        remove(cov_temp)

    if exists(temp):
        remove(temp)

if __name__=="__main__":
    parser = argparse.ArgumentParser(
        description="Wrapper around fastQTL to do multiple eQTL mapping.")

    parser.add_argument("--vcf",
                        help="Specify vcf file.",
                        required=True)
    parser.add_argument("--bed",
                        help="Specify bed file.",
                        required=True)
    parser.add_argument("--results",
                        help="Specify list of all eQTLs.",
                        required=True)
    parser.add_argument("--forward",
                        help="File to write forward stage results to.",
                        default="Forward")
    parser.add_argument("--backward",
                        help="File to write backwards stage results to.",
                        default="Backward")
    parser.add_argument("--cov",
                        help="Covariates file.",
                        default="")
    parser.add_argument("--include-phenotypes",
                        help="File specifying phenotypes to include.",
                        required=False)
    parser.add_argument("--include-samples",
                        help="File specifying samples to include.",
                        required=False)
    parser.add_argument("--exon",
                        help="Specify exon groupings for analysis.")
    parser.add_argument("--perm",
                        help="Number of permutations.",
                        type=int,
                        default=10000)
    parser.add_argument("-W",
                        help="Cis window.",
                        type=int,
                        default=1000000)
    parser.add_argument("--seed",
                        help="Seed number.",
                        type=int,
                        default=1461167480)
    parser.add_argument("--normal",
                        help="Pass --normal flag to fastQTL.",
                        action='store_true',
                        default=False)
    parser.add_argument("--genes",
                        help="Number of genes analysed in each chunk.",
                        type=int,
                        default=0)
    parser.add_argument("--jobnumber",
                        help="Chunk to run.",
                        type=int,
                        default=1)
    parser.add_argument("--threshold",
                        help="Threshold on adjusted P values.",
                        type=float,
                        default=1)
    parser.add_argument("--fastQTL",
                        help="Path to fastQTL.",
                        default="fastQTL.static")
    parser.add_argument("--tmpdir",
                        help="Path to outputdir.",
                        default=".")
    parser.add_argument("--overwrite",
                        help="Allow overwriting of files.",
                        action='store_true',
                        default=False)

    args = parser.parse_args()

    run_analysis(args)

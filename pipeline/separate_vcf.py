#!/usr/bin/env python

"""
Separate VCF with multiple samples into individual VCFs with only one sample.
"""

import sys
import os
import re
from subprocess import check_call
from collections import defaultdict
from argparse import ArgumentParser

base_dir = os.path.dirname(__file__) or '.'
lib_dir = os.path.abspath(os.path.join(base_dir, 'lib_cftr'))
sys.path.insert(0, lib_dir)

from common.misc import open_file, have_file
from cftr import EXE, RESOURCE

gatkJar = EXE['gatk_jar']
gatkExe = ["java", "-Xmx4g", "-jar", gatkJar]
if not os.path.isfile(gatkJar):
    sys.exit("Could not find {}\n".format(gatkJar))

# Create separate VCFs from joint VCF
#gatk.sh -R ../../../resources/ref/CFTR_genomic_refseq.fa -T SelectVariants
#--excludeNonVariants --excludeFiltered --variant 48samplesL.gatk-merged.vcf -o
#Sample01_S1-CFTR_genomic_refseq.vcf -sn Sample01_S1


def separate_vcf(vcffile, ref, outdir, args):
    sys.stderr.write("\nReading {}\n".format(vcffile))
    fh = open_file(vcffile)
    fieldlist = [ l for l in fh.readlines() if l.startswith('#CHROM') ]
    fh.close()
    if fieldlist:
        fields = fieldlist[0].rstrip().split("\t")
        samples = fields[9:]
    else:
        sys.stderr.write("Could not find VCF header in {}".format(vcffile))
        sys.exit(1)
    sys.stderr.write("  Found {} samples.\n".format(len(samples)))
    outvcfs = []
    logfh = sys.stderr
    if args.logdir:
        logfile = os.path.basename(vcffile).replace('.vcf','') +\
                  ".separate_vcf.log"
        logfile = os.path.join(args.logdir, logfile)
        logfh = open(logfile, 'w')
    for sample in samples:
        outvcf = "{}.separated.vcf".format(sample)
        if outdir: outvcf = os.path.join(outdir, outvcf)
        run_select_variants(vcffile, outvcf, sample, ref, logfh, args)
    if args.logdir:
        logfh.close()
    return samples

def run_select_variants(vcffile, outvcf, sample, ref, logfh, args):
    logfh.write("\n-- SelectVariants --\n")
    cmd = gatkExe[:]
    cmd.extend(['-T', 'SelectVariants', '--excludeNonVariants']) 
    cmd.extend(['-R', ref])  
    cmd.extend(['--variant', vcffile])  
    cmd.extend(['-o', outvcf])  
    cmd.extend(['-sn', sample])  
    logfh.write("  Sample {}:\t{}\n".format(sample, outvcf))
    if have_file(outvcf, args.force):
        sys.stderr.write("  Already have {}.\n".format(outvcf))
    else: 
        logfh.write(" ".join(cmd)+"\n")
        logfh.flush()
        check_call(cmd, stderr=logfh)
        if os.path.isfile(outvcf):
            sys.stderr.write("  Created {}\n".format(outvcf))
        else:
            sys.stderr.write("  Failed to create {}\n".format(outvcf))
    return outvcf


if __name__ == '__main__':
    descr = "Separate VCF with multiple samples to one sample per VCF file."
    parser = ArgumentParser(description=descr)
    parser.add_argument("vcffiles", nargs="+", help="VCF file(s)")
    parser.add_argument("-o", "--outdir", help="Directory for output files.")
    parser.add_argument("-f", "--force", default=False, action="store_true",
                        help="Overwrite existing files.")
    parser.add_argument("--logdir", help="Directory for log files.")

    if len(sys.argv)<2:
        parser.print_help()
        sys.exit()
    args = parser.parse_args()

    if not os.path.isfile(RESOURCE['ref_fa']):
        sys.exit("Ref file {} not found.\n".format(RESOURCE['ref_fa']))
    for vcffile in args.vcffiles:
        separate_vcf(vcffile, RESOURCE['ref_fa'], args.outdir, args)


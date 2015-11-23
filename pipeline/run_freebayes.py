#!/usr/bin/env python

"""
Run freebayes
"""

import sys
import os
import re
import multiprocessing as mp
from subprocess import check_call
from collections import defaultdict
from argparse import ArgumentParser

base_dir = os.path.dirname(__file__) or '.'
lib_dir = os.path.abspath(os.path.join(base_dir, 'lib_cftr'))
sys.path.insert(0, lib_dir)

from common.misc import open_file, have_file, timestamp
from common.profile import profile
from cftr import EXE, RESOURCE

freebayesExe = EXE['freebayes']
if not os.path.isfile(freebayesExe):
    sys.exit("Could not find {}\n".format(freebayesExe))
vcffilterExe = os.path.join(base_dir, 'vcf_filter_freebayes.py')
if not os.path.isfile(vcffilterExe):
    sys.exit("Could not find {}\n".format(vcffilterExe))

def get_outlabel(bamfile, outdir=None):
    label = os.path.basename(bamfile).rstrip('bam').rstrip('.')
    if outdir:
        if os.path.isdir(outdir):
            label = os.path.join(outdir, label)
        else:
            msg = "ERROR: Outdir {} does not exist.\n".format(outdir)
            sys.stderr.write(msg)
            sys.exit(1)
    return label

@profile
def run_freebayes(bamfile, reffile, bedfile, logfh, args):
    label = get_outlabel(bamfile, args.outdir)
    outvcf = label + '.freebayes.vcf'
    cmd = [ freebayesExe, '-f', reffile, '-t', bedfile ]
    cmd.extend(['-b', bamfile, '-v', outvcf])  
#    cmd.extend(['--max-complex-gap', '5',])  
    if args.basequal:
        cmd.extend(['-q', str(args.basequal)])  
    if args.minaltcount: 
        cmd.extend(['-C', str(args.minaltcount)])
    logfh.write("\nCMD: "+" ".join(cmd)+"\n")
    if have_file(outvcf, args.force):
        sys.stderr.write("  Already have {}.\n".format(outvcf))
    else: 
        logfh.flush()
        check_call(cmd, stderr=logfh)
    if os.path.isfile(outvcf):
        sys.stderr.write("  Freebayes result in {}\n".format(outvcf))
    else:
        sys.stderr.write("  Failed to create {}\n".format(outvcf))
        sys.exit(1)
    return outvcf

def filter_freebayes_vcf(vcffile, outvcf, logfh, args):
#    cmd = [ vcffilterExe, '-f', "QUAL > 20", "-f", "DP > 10" ]
    cmd = [  vcffilterExe, '-s']
    if args.outdir:
        cmd.extend(['-o', args.outdir])
    if args.altbasequal:
        cmd.extend(['-a', str(args.altbasequal)])
    if args.dp:
        cmd.extend(['-d', str(args.dp)])
    if args.qual:
        cmd.extend(['-q', str(args.qual)])
    cmd.append(vcffile)
    logfh.write("\nCMD: "+" ".join(cmd)+"\n")
    if have_file(outvcf, args.force):
        sys.stderr.write("  Already have {}.\n".format(outvcf))
    else: 
        logfh.flush()
        check_call(cmd, stderr=logfh)
    if os.path.isfile(outvcf):
        sys.stderr.write("  Filtered FreeBayes result {}\n".format(outvcf))
    else:
        sys.stderr.write("  Failed to create {}\n".format(outvcf))
        sys.exit(1)
    return outvcf

def call_variants_freebayes(bamfile, reffile, bedfile, args):
    sys.stderr.write("Bam file: {}\n".format(bamfile))
    if args.logdir:
        label = get_outlabel(bamfile, args.logdir)
        outlog = label +'.freebayes.log'
        logfh = open(outlog, 'w')
    else:
        logfh = sys.stderr
    tag = "{}:\t".format(get_outlabel(bamfile))
    logfh.write("{}start time {}\n".format(tag, timestamp()))
    vcffile = run_freebayes(bamfile, reffile, bedfile, logfh, args)
    outvcf = vcffile.replace('.vcf', '')+".filtered.vcf"
    filteredvcf = filter_freebayes_vcf(vcffile, outvcf, logfh, args)
    vcfs.append(filteredvcf)
    logfh.write("{}finish time {}\n".format(tag, timestamp()))
    if args.logdir:
        logfh.close()
    return vcfs

if __name__ == '__main__':
    descr = "Run freebayes on bamfile(s)."
    parser = ArgumentParser(description=descr)
    parser.add_argument("bamfiles", nargs="+", help="BAM file(s)")
    parser.add_argument("-o", "--outdir", help="Directory for output files.")
    parser.add_argument("-b", "--basequal", type=int, 
                        help="Set --min-base-quality to this value.")
    parser.add_argument("-c", "--minaltcount", type=int, 
                        help="Set --min-alternate-count to this value.",)
    parser.add_argument("-a", "--altbasequal", type=int, default=20,
                        help="Min alt base quality to pass (default: 20).",)
    parser.add_argument("-d", "--dp", type=int, default=20,
                        help="Min DP to pass (default: 20).",)
    parser.add_argument("-q", "--qual", type=int, default=20,
                        help="Min QUAL to pass (default: 20).",)
    parser.add_argument("-f", "--force", default=False, action="store_true",
                        help="Overwrite existing files.")
    parser.add_argument("-p", "--processes", default=1, type=int,
                        help="Number of processes to run.")
    parser.add_argument("--logdir", help="Directory for log files.")

    if len(sys.argv)<2:
        parser.print_help()
        sys.exit()
    args = parser.parse_args()

    vcfs = []
    reffile = RESOURCE['ref_fa']
    if not os.path.isfile(reffile):
        sys.exit("Ref file {} not found.\n".format(reffile))
    bedfile = RESOURCE['covered_roi_bed']
    if not os.path.isfile(bedfile):
        sys.exit("BED file {} not found.\n".format(bedfile))
    sys.stderr.write("\nRunning freebayes on {} bam files using {}.\n".format(
                   len(args.bamfiles), "{} process{}".format(args.processes,
                   'es' if args.processes>1 else '')))
    pool = mp.Pool(processes=args.processes)
    results = [ pool.apply_async(call_variants_freebayes, args=[bamfile, 
                reffile, bedfile, args]) for bamfile in args.bamfiles ]
    vcfs = [p.get() for p in results ]


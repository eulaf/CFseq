#!/usr/bin/env python

"""
Run flagstat CFTR bam files and results compile into table.

Compiles stats will include:
1. QC passed reads
2. Mapped reads
3. Properly paired reads
"""

import sys
import os
import statistics
import pysam
from collections import defaultdict
from argparse import ArgumentParser

base_dir = os.path.dirname(__file__) or '.'
lib_dir = os.path.abspath(os.path.join(base_dir, '..', 'lib_cftr'))
sys.path.insert(0, lib_dir)

from common.misc import have_file

def get_sample_name(bamfile):
    sample = os.path.basename(bamfile).split('-')[0]
    return sample

def perc_from_flagstat_line(l):
    (num, x) = l.rstrip().split(' ', 1)
    perc = x.split('(',1)[1].split(':',1)[0]
    return perc

def run_flagstat(bamfile):
    stats = {}
    for l in pysam.flagstat(bamfile):
        if 'QC-passed' in l:
            stats['QC-passed reads'] = l.rstrip().split(' ', 1)[0]
        elif 'mapped' in l:
            stats['% Mapped'] = perc_from_flagstat_line(l)
        elif 'properly paired' in l:
            stats['% Properly paired'] = perc_from_flagstat_line(l)
            break
    sys.stderr.write("{}\n".format(stats))
    return stats

def print_summary(bamstats, outfile="summary.txt"):
    bamfiles = sorted(bamstats.keys())
    samples = map(get_sample_name, bamfiles)
    sys.stderr.write("\nWriting {}\n".format(outfile))
    with open(outfile, "w") as ofh:
        ofh.write("\t"+"\t".join(samples)+"\n")
        for f in ('QC-passed reads', '% Mapped', '% Properly paired'):
            row = [ bamstats[b][f] for b in bamfiles ]
            ofh.write(f + "\t" + "\t".join(row)+"\n") 
        ofh.write("\n")
#        ofh.write("\nBam files:\n"+"\n".join(bamfiles)+"\n")
        
if __name__ == '__main__':
    descr = "Create table of flagstat values for CFTR bam file(s)."
    parser = ArgumentParser(description=descr)
    parser.add_argument("bamfiles", nargs="+", help="BAM file(s)")
    parser.add_argument("-o", "--outfile", default="summary.txt",
                        help="Name for summary output file.")
    parser.add_argument("-f", "--force", default=False, action="store_true",
                        help="Overwrite existing files.")

    if len(sys.argv)<2:
        parser.print_help()
        sys.exit()
    args = parser.parse_args()

    if have_file(args.outfile, args.force):
        sys.stderr.write("  Already have {}\n".format(args.outfile))
        sys.exit()
    bamstats = {}
    for bamfile in args.bamfiles:
        sys.stderr.write("\nReading bam file: {}\n".format(bamfile))
        bamstats[bamfile] = run_flagstat(bamfile)
    print_summary(bamstats, args.outfile)

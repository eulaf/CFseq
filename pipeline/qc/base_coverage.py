#!/usr/bin/env python

"""
Count reads per base within the ROI for each bam file.
Calculate percent of bases covered by at least 20, 100 reads.
Calculate minimum coverage to cover all bases in ROI.
"""

import sys
import os
import pysam
import statistics
from collections import defaultdict
from argparse import ArgumentParser

base_dir = os.path.dirname(__file__) or '.'
lib_dir = os.path.abspath(os.path.join(base_dir, '..', 'lib_cftr'))
sys.path.insert(0, lib_dir)

from common.misc import have_file
import cftr

THRESHOLDS = (100, 200, 500, 1000)

def parse_roi(state=False):
    bedfile = cftr.RESOURCE['analysis_roi_bed'] if state else \
              cftr.RESOURCE['covered_roi_bed']
    bed = cftr.parse_bedfile(bedfile)
    regions = []
    for chrom in bed.keys():
        for start in sorted(bed[chrom].keys()):
            name = bed[chrom][start]['NAME']
            regions.append([chrom, start, bed[chrom][start]['END'], name])
    sys.stderr.write("ROI: {} regions\n".format(len(regions)))
    return regions

def get_coverage(sample, bamfile, roi):
#    sys.stderr.write("Bam file {}\n".format(bamfile))
    bampysam = pysam.AlignmentFile(bamfile, "rb")
    basecoverage = {}
    total = 0
    for region in roi:
        (ref, start, end, name) = region
        total += end - start + 1
#        sys.stderr.write("{}\t{}-{}\n".format(name, start, end))
        for pileupcol in bampysam.pileup(ref, start, end, stepper="samtools",
                                         max_depth=99999):
            pos = pileupcol.reference_pos+1 # convert to 1-based pos
            if start <= pos <= end:
                basecoverage[pos] = pileupcol.nsegments
    cov_vals = basecoverage.values()
#    total = len(cov_vals)
    mincov = min(cov_vals)
    coverage = {'lowest_coverage':mincov, 'total_roi_bases':total,}
    sys.stderr.write("{}\tbases={}\tmin={}".format(sample, total, mincov))
    for thresh in THRESHOLDS:
        coverage[thresh] = len([ n for n in cov_vals if n>thresh])*100.0/total
        sys.stderr.write("\t{}={:.1f}%".format(thresh, coverage[thresh]))
    sys.stderr.write("\n")
    return coverage

def compile_data(cov, outfile):
    fields = ['total_roi_bases', 'lowest_coverage', ]
    samples = sorted(cov.keys())
    sys.stderr.write("\nWriting {}\n".format(outfile))
    with open(outfile, 'w') as ofh:
        ofh.write("\t"+"\t".join(samples)+"\n")
        for f in fields:
            row = [ str(cov[s][f]) if f in cov[s] else '0' for s in samples ]
            ofh.write(f + "\t" + "\t".join(row)+"\n")
        for thresh in THRESHOLDS:
            row = [ "{:.1f}%".format(cov[s][thresh]) if thresh in cov[s] 
                    else 0 for s in samples ]
            ofh.write("coverage>{}".format(thresh)+"\t"+"\t".join(row)+"\n")
        ofh.write("\n")


if __name__ == '__main__':
    descr = "Count amplicon coverage and calculate coverage uniformity."
    parser = ArgumentParser(description=descr)
    parser.add_argument("bamfiles", nargs="+", help="Bam files.")
    parser.add_argument("-o", "--outfile", default="basecoverage.txt",
                        help="Name for output file.")
    parser.add_argument("-s", "--state", default=False, action="store_true",
                        help="Use state ROI.")
    parser.add_argument("-f", "--force", default=False, action="store_true",
                        help="Overwrite existing files.")
    parser.add_argument("-d", "--debug", default=False, action="store_true",
                        help="Keep intermediate files.")

    if len(sys.argv)<2:
        parser.print_help()
        sys.exit()
    args = parser.parse_args()

    if have_file(args.outfile, args.force):
        sys.stderr.write("  Already have {}\n".format(args.outfile))
        sys.exit()
    roi = parse_roi(args.state)
    cov = {}
    for bamfile in sorted(args.bamfiles):
        sample = os.path.basename(bamfile).split('-')[0].split('.')[0]
        cov[sample] = get_coverage(sample, bamfile, roi)
    compile_data(cov, args.outfile)

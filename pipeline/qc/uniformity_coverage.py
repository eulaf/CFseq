#!/usr/bin/env python

"""
Count reads per amplicon and then calculating uniformity of coverage.

Reads per amplicon is determined by finding expected reads in the
primer2reads.txt file and then finding if the read mapped to the expected
region in the bam file.

Uniformity of coverage is calculated as percent of amplicons with
coverage over 0.2x and 0.5x mean amplicon coverage.
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

from common.misc import open_file, have_file
import cftr

def parse_roi():
    bed = cftr.parse_bedfile(cftr.RESOURCE['amplicon_bed'])
    regions = []
    for chrom in bed.keys():
        for start in sorted(bed[chrom].keys()):
            name = bed[chrom][start]['NAME']
            regions.append([chrom, start, bed[chrom][start]['END'], name])
    sys.stderr.write("ROI: {} regions\n".format(len(regions)))
    return regions

def get_primer2reads_files(primer2readsdir):
    sys.stderr.write("\nReading dir {}\n".format(primer2readsdir))
    p2rfiles = [ f for f in os.listdir(primer2readsdir) if
              f.endswith('.primer2reads.txt') ]
    p2rdict = {}
    for p2rfile in p2rfiles:
        sample = p2rfile.split('.')[0]
        p2rdict[sample] = os.path.join(primer2readsdir, p2rfile)
    sys.stderr.write("  Primer2reads.txt files: {} files\n".format(
                     len(p2rdict.keys())))
    return p2rdict

def parse_primer2reads_file(p2rfile):
    sys.stderr.write("\nReading file {}\n".format(p2rfile))
    p2r_content = defaultdict(list)
    with open(p2rfile, 'r') as fh:
        head = fh.readline()
        for line in fh:
            (primer, num_reads, readstr) = line.rstrip('\n\r').split("\t")
            reads = [ r for r in readstr.split(", ") ]
            name = ''
            if primer.endswith('_F'):
                name = primer.rstrip('_F')
            elif primer.endswith('_R'):
                name = primer.rstrip('_R')
            else:
                sys.stderr.write("  Skipping {}\n".format(primer))
                continue
            p2r_content[name].extend(reads)
    sys.stderr.write("  Got {} amplicons\n".format(len(p2r_content)))
    return p2r_content

def get_amplicon_coverage(sample, bamfile, p2r, roi):
    sys.stderr.write("Bam file {}\n".format(bamfile))
    bampysam = pysam.AlignmentFile(bamfile, "rb")
    amp_readcounts = {}
    num_mapped = bampysam.mapped
    for region in roi:
        (ref, start, end, name) = region
        if name in p2r:
            expected_num = len(p2r[name])
            expected_reads = set(p2r[name])
            num = 0
            for alignedsegment in bampysam.fetch(ref, start, end):
                rdnum = '/R1' if alignedsegment.is_read1 else '/R2'
                if alignedsegment.query_name+rdnum in expected_reads:
                    num += 1
            amp_readcounts[name] = num
    coverages = amp_readcounts.values()
    covtot = sum(coverages)
    covmin = min(coverages)
    covmax = max(coverages)
    covmean = statistics.mean(coverages)
    amp_readcounts['total'] = str(covtot)
    amp_readcounts['max/min'] = "{:.1f}".format(covmax/covmin \
                                 if covmin>0 else 0)
    amp_readcounts['mean'] = "{:.0f}".format(covmean)
    for frac in (0.2, 0.5):
        frac_covmean = covmean*frac
        amp_readcounts[str(frac)+'X mean'] = "{:.0f}".format(frac_covmean)
        numgtmn = len([ n for n in coverages if n > frac_covmean])
        amp_readcounts['%>{}X mean'.format(frac)] = "{:.1f}%".format(
            numgtmn*100.0/len(coverages) if coverages else 0)
    sys.stderr.write("  {}:\tmapped={}\ttot={}\tmin={}\tmax={}\tmean={:.0f}\n".\
                  format(name, num_mapped, covtot, covmin, covmax, covmean))
    return amp_readcounts

def compile_data(ampcov, roi, outfile):
    fields = ['max/min', 'mean', '0.2X mean', '%>0.2X mean',
              '0.5X mean', '%>0.5X mean',]
    samples = sorted(ampcov.keys())
    sys.stderr.write("\nWriting {}\n".format(outfile))
    with open(outfile, 'w') as ofh:
        ofh.write("Amplicon\t"+"\t".join(samples)+"\n")
        for region in roi:
            (ref, start, end, name) = region
            row = [ str(ampcov[s][name]) if name in ampcov[s] else '0' \
                    for s in samples ]
            ofh.write(name + "\t" + "\t".join(row)+"\n")
        ofh.write("\n")
        for f in fields:
            row = [ str(ampcov[s][f]) if f in ampcov[s] else '0' \
                    for s in samples ]
            ofh.write(f + "\t" + "\t".join(row)+"\n")


if __name__ == '__main__':
    descr = "Count amplicon coverage and calculate coverage uniformity."
    parser = ArgumentParser(description=descr)
    parser.add_argument("primer2readsdir", 
                        help="Directory with primer2reads.txt files.")
    parser.add_argument("bamfiles", nargs="+", help="Bam files.")
    parser.add_argument("-o", "--outfile", default="uniformity.txt",
                        help="Name for output file.")
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
    roi = parse_roi()
    p2rdict = get_primer2reads_files(args.primer2readsdir)
    ampcov = {}
    for bamfile in sorted(args.bamfiles):
        sample = os.path.basename(bamfile).split('-')[0].split('.')[0]
        if not sample in p2rdict:
            sys.stderr.write("No primer2read file for {}\n".format(bamfile))
            continue
        p2r = parse_primer2reads_file(p2rdict[sample])
        ampcov[sample] = get_amplicon_coverage(sample, bamfile, p2r, roi)
    compile_data(ampcov, roi, args.outfile)

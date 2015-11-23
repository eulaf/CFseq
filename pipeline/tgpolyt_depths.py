#!/usr/bin/env python

"""
Count reads supporting each TG-polyT variant.
"""

import sys
import os
import pysam
from collections import defaultdict
from argparse import ArgumentParser

base_dir = os.path.dirname(__file__) or '.'
lib_dir = os.path.abspath(os.path.join(base_dir, 'lib_cftr'))
sys.path.insert(0, lib_dir)

from common.misc import open_file, have_file

# location of TG-polyT in CFTR reference
REGION = ['CFTR', 87824, 87852, 'TG-polyT']
HET_FREQ = 0.25
HOM_FREQ = 0.55

HOMCALLS = []
HETCALLS = []

def get_reads_covering_region(bamfile, region, mapq=1):
    """Get reads covering TG-polyT region."""
    bampysam = pysam.Samfile(bamfile, "rb")
    (refseq, s, e, region_name) = region
    # Get only mapped reads and reads with mapping quality > 0
    mapped_reads = [ r for r in bampysam.fetch(refseq, s, e) \
                     if not (r.is_unmapped or r.mapq < mapq) ]
    totreads = len(mapped_reads)
    sys.stderr.write("Total of {} reads cover region {}\n".format(totreads,
                     region))
    return mapped_reads

def get_tgpolyt_bases(read, s, e, flank=1):
    """Get the TG-polyT bases based on read alignment to reference."""
    refpos = read.get_reference_positions()
    if refpos[0] >= s or refpos[-1] <= e: 
        # TG-T not fully covered by read, so skip
        return ''
    i_range = [ i for i in refpos if s<=i<=e ]
    i_start = refpos.index(i_range[0])
    i_end = refpos.index(i_range[-1])
    seq = read.query_alignment_sequence
    tgpolyt = seq[i_start:i_end]
    while tgpolyt.startswith('T') or tgpolyt.startswith('G'):
        if i_start==0: break
        i_start -= 1
        tgpolyt = seq[i_start:i_end]
    while tgpolyt.endswith('T'):
        if i_end==len(refpos)-1: break
        i_end += 1
        tgpolyt = seq[i_start:i_end]
    return tgpolyt

def count_tg_t(subseq):
    """Count number of TG's and T's in subseq and return (TG)#-#T string"""
    try:
        xsubseq = subseq.replace('TG','X')
        start = xsubseq.index('X')
    except ValueError:
        return ""
    base = 'X'
    num = {'X':0, 'T':0}
    for b in xsubseq[start:]:
        if b==base: 
            num[base] += 1
        elif b=='T' and base=='X':
            base = 'T'
            num[base] += 1
        else:
            break
#    if 9 <= num['X'] <= 13 and num['T'] in (5,7,9):
    return "(TG){}-{}T".format(num['X'], num['T'])
#    else:
#        return ""

def count_tgpolyt(mapped_reads, region, debug=False):
    (refseq, s, e, region_name) = region
    tgpolyt = defaultdict(int)
    totreads = 0
    for read in mapped_reads:
        name = read.query_name
        subseq = get_tgpolyt_bases(read, s, e)
        if not subseq: continue
        tgt = count_tg_t(subseq)
        if not tgt: continue
        tgpolyt[tgt] += 1
        totreads += 1
    sys.stderr.write("{} TG-polyT reads counted\n".format(totreads))
    if debug:
        for tgt in sorted(tgpolyt.keys(), key=lambda x:(len(x),x)):
            vf = float(tgpolyt[tgt])/totreads
            sys.stderr.write("{}\t{:>5}\t{:.3f}\n".format(tgt, tgpolyt[tgt], 
                             vf))
    return (tgpolyt, totreads)

def hethom_call(tgtlist, vflist):
    hethom = []
    for (tgt, vf) in zip(tgtlist, vflist):
        if vf>HOM_FREQ: 
            HOMCALLS.append(vf)
            if tgt=='11-7T':
                hethom.append('hom')
            else:
                hethom.append('hom alt')
        elif vf>HET_FREQ: 
            HETCALLS.append(vf)
            hethom.append('het')
        else:
            hethom.append('--')
    return '/'.join(hethom)

def tgt_sortkey(x):
    (tg,t) = x.lstrip('(TG)').rstrip('T').split('-')
    return (int(t), int(tg))

def report_results(ofh, sfh, sample, tgpolyt, totreads, args):
    summtgt = []
    summvf = []
    summvarfreq = []
    summreads = []
    for tgt in sorted(tgpolyt.keys(), key=tgt_sortkey):
        numreads = tgpolyt[tgt]
        varfreq = float(numreads)/totreads
        vf = "{:.3f}".format(varfreq)
#        print "{}\t{}\t{}\t{}".format(tgt, vf, tgpolyt[tgt], totreads)
        row = (sample, tgt, vf, numreads, totreads)
        ofh.write("\t".join([str(r) for r in row])+"\n")
        if varfreq < args.lowfreq: continue
        if numreads < args.minreads: continue
        summtgt.append(tgt.lstrip('(TG)'))
        summvf.append(vf)
        summvarfreq.append(varfreq)
        summreads.append(str(numreads))
    sfh.write("{}\t{}\t{}\t{}\t{}\n".format(sample, 
        '/'.join(summtgt), hethom_call(summtgt, summvarfreq), 
        '/'.join(summvf), '/'.join(summreads)))

def get_samplename(bamfile):
    sample = os.path.basename(bamfile).split('-')[0].split('.')[0]
    return sample

if __name__ == '__main__':
    descr = "Count read depth for TG-polyT region of  CFTR bam file(s)."
    parser = ArgumentParser(description=descr)
    parser.add_argument("bamfiles", nargs="+", help="BAM file(s)")
    parser.add_argument("-l", "--lowfreq", type=float, default=HET_FREQ,
                        help="Low frequency cut-off.")
    parser.add_argument("-m", "--minreads", type=int, default=100,
                        help="Minimum reads.")
    parser.add_argument("-o", "--outlabel", default="tgpolyt_results",
                        help="Label for output files.")
    parser.add_argument("--debug", action="store_true", default=False,
                        help="Debug mode.")
    parser.add_argument("-f", "--force", default=False, action="store_true",
                        help="Overwrite existing files.")

    if len(sys.argv)<2:
        parser.print_help()
        sys.exit()
    args = parser.parse_args()

    outfile = args.outlabel + ".tgpolyt_counts.txt"
    summaryfile = args.outlabel + ".tgpolyt.txt"
    sys.stderr.write("Writing {}\n".format(outfile))
    sys.stderr.write("Writing {}\n".format(summaryfile))
    if have_file(outfile, args.force) and have_file(summaryfile, args.force):
        sys.stderr.write("Already have {} and {}\n".format(
                         outfile, summaryfile))
        sys.exit()
    outfields = ['sample', 'TG-polyT', 'frequency', 'num_reads', 
                 'tot_reads']
    summfields = ['sample', 'TG-polyT', 'hom_het', 'frequency', 
                  'num_reads']
    with open(outfile, 'w') as ofh, open(summaryfile, 'w') as sfh:
        ofh.write("\t".join(outfields)+"\n")
        sfh.write("\t".join(summfields)+"\n")
        for bamfile in args.bamfiles:
            sample = get_samplename(bamfile)
            sys.stderr.write("\nReading bam file: {}\n".format(bamfile))
            sys.stderr.write("  Sample: {}\n".format(sample))
            reads = get_reads_covering_region(bamfile, REGION)
            (tgpolyt, totreads) = count_tgpolyt(reads, REGION, args.debug)
            report_results(ofh, sfh, sample, tgpolyt, totreads, args)
        hetvf = sorted(HETCALLS)
        homvf = sorted(HOMCALLS)
        ofh.write("Hom call range: {:.3f}-{:.3f}\n".format(homvf[0], homvf[-1]))
        ofh.write("Het call range: {:.3f}-{:.3f}\n".format(hetvf[0], hetvf[-1]))
        sys.stderr.write("Hom call range: {:.3f}-{:.3f}\n".format(homvf[0], homvf[-1]))
        sys.stderr.write("Het call range: {:.3f}-{:.3f}\n".format(hetvf[0], hetvf[-1]))

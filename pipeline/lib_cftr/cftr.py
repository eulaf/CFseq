#!/usr/bin/env python

import os
import sys
from collections import defaultdict
from common.misc import open_file

base_dir = os.path.dirname(__file__) or '.'
bin_dir = os.path.abspath(os.path.join(base_dir,'..','..','bin'))
ref_dir = os.path.abspath(os.path.join(base_dir,'..','..','resources'))
EXE = {
    'annovar_convert': os.path.join(bin_dir, 'convert2annovar.pl'),
    'annovar_table': os.path.join(bin_dir, 'table_annovar.pl'),
    'annovar_variation': os.path.join(bin_dir, 'annotate_variation.pl'),
    'bwa': os.path.join(bin_dir, 'bwa'),
    'cm': os.path.join(bin_dir, 'cross_match.manyreads'),
    'freebayes': os.path.join(bin_dir, 'freebayes'),
    'gatk_jar': os.path.join(bin_dir, 'GenomeAnalysisTK.jar'),
    'picard_jar': os.path.join(bin_dir, 'picard.jar'),
    'samtools': os.path.join(bin_dir, 'samtools'),
}
RESOURCE = {
    'ref_fa':os.path.join(ref_dir, 'ref', 'CFTR_genomic_refseq.fasta'),
    'primer_fa': os.path.join(ref_dir, 'primer_info', 'CFTR_primers.fa'),
    'cftr_db': os.path.join(ref_dir, 'cftr_db', 'CFTR_db-variant.txt'),
    'dbsnp_vcf': os.path.join(ref_dir, 'variants',
                              'dbsnp_138.cftr.localref.vcf'),
    'dbsnp_edited': os.path.join(ref_dir, 'variants',
                              'dbsnp_138.cftr.edited.vcf'),
    '1000G_vcf': os.path.join(ref_dir, 'variants',
                              '1000G_phase1.indels.cftr.sites.localref.vcf'),
    'Mills_vcf': os.path.join(ref_dir, 'variants',
              'Mills_and_1000G_gold_standard.indels.cftr.sites.localref.vcf'),
    'covered_roi': os.path.join(ref_dir, 'roi', 'CFTR_covered.interval_list'),
    'covered_roi_bed': os.path.join(ref_dir, 'roi', 'CFTR_covered.bed'),
    'analysis_roi_bed': os.path.join(ref_dir, 'roi', 'CFTR_roi.bed'),
    'amplicon_bed': os.path.join(ref_dir, 'roi', 'CFTR_amplicons.bed'),
    'annovar_ref': os.path.join(ref_dir, 'annovar_dbs'),
}

def hg19_to_CFTR(chrom, pos):
    """Convert an hg19 location to position on CFTR.  All positions should be
    in CFTR region chr7:117100838-117358025."""
    if chrom=='chr7' or chrom=='7' or chrom==7:
        chrom = 'CFTR'
        pos_cftr = int(pos) - 117100838 + 1
    else:
        pos_cftr = pos
    return (chrom, pos_cftr)

def CFTR_to_hg19(chrom, pos):
    """Convert a CFTR location to position on hg19.  All positions should be
    in CFTR region chr7:117100838-117358025."""
    if chrom=='CFTR':
        chrom = 'chr7'
        pos_hg19 = int(pos) - 1 + 117100838
    else:
        pos_hg19 = pos
    return (chrom, pos_hg19)

def parse_bedfile(bedfile):
    """Parse ROI from bed file."""
    sys.stderr.write("\nReading {}\n".format(bedfile))
    fh = open_file(bedfile)
    lines = fh.readlines()
    fh.close()

    roi = defaultdict(dict)
    fields = ['CHROM', 'START', 'END', 'NAME']
    numlines = 0
    for line in lines:
        vals = line.rstrip().split("\t")
        d = dict(zip(fields, vals))
        d['START'] = int(d['START']) + 1 # convert to 1-based
        d['END'] = int(d['END'])
        roi[d['CHROM']][d['START']] = d
        numlines += 1
    sys.stderr.write("  Got {} lines in {} chrom\n".format(numlines,
                     len(roi.keys())))
    return roi

def in_roi(chrom, loc, roi):
    pos = int(loc)
    roi_flag = chrom
    if chrom in roi:
        roi_flag = 'N'
        for s in sorted(roi[chrom].keys()):
            e = roi[chrom][s]['END']
            if pos >= s:
                if pos <=e:
                    roi_flag = "Y = " + roi[chrom][s]['NAME']
            elif pos > e:
                break
    return roi_flag
   
def rev_complement(seq):
    def rev_complement(seq):
        intab =  'AGCTURYSWKMBVDHagcturyswkmbvdh'
        outtab = 'TCGAAYRSWMKVBHDtcgaayrswmkvbhd';
        trantab = string.maketrans(intab, outtab)
        return seq.translate(trantab)[::-1]

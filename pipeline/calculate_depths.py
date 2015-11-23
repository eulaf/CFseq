#!/usr/bin/env python

"""
Calculate depths from bam file for VCF variant positions.
"""

import sys
import os
import re
import pysam
import multiprocessing as mp
from pipes import quote
from collections import defaultdict
from argparse import ArgumentParser

base_dir = os.path.dirname(__file__) or '.'
lib_dir = os.path.abspath(os.path.join(base_dir, 'lib_cftr'))
sys.path.insert(0, lib_dir)

from common.misc import open_file, have_file, timestamp
from cftr import RESOURCE

PILECHAR_FWD = 'ACTGN.'#"".join([ chr(ord('A')+i) for i in xrange(26) ]) + "."
PILECHAR_REV = 'actgn,'#"".join([ chr(ord('a')+i) for i in xrange(26) ]) + ","

def organize_input_files(inputfilelist, delim = r'[-\.]'):
    """Pairs bam and vcf files.  For a bam file named samplename.bam, 
    the possible vcf filenames to pair with are: samplename.vcf, 
    samplename.gatk.vcf or samplename.gatk-cohort.vcf."""
    infiles = defaultdict(dict)
    for infile in inputfilelist:
        (filename, ext) = os.path.splitext(os.path.basename(infile))
        (label, other) = re.split(delim, filename, 1)
        ext = ext.lstrip('.')
        if ext not in ('bam', 'vcf'):
            sys.stderr.write("Unrecognized file type ({}): {}\n".format(
                ext, infile))
            continue
        infiles[label][ext] = infile
    num = 0
    for label in infiles.keys():
        if 'bam' in infiles[label] and 'vcf' in infiles[label]:
            num += 1
        else:
            files = ', '.join(infiles[label].values())
            sys.stderr.write("Not properly paired: {}\n".format(files))
            infiles.pop(label, None)
    sys.stderr.write("Total of {} paired bam and vcf files.\n".format(num))
    return infiles

def parse_vcf(vcffile, logfh):
    """Parse sample VCF file."""
    logfh.write("Parsing vcffile: {} -".format(vcffile))
    fh = open_file(vcffile)
    lines = fh.readlines()
    fh.close()
    while lines: # skip vcf header
        if not lines[0].startswith('##'):
            break
        lines.pop(0)
    fields = lines.pop(0).lstrip('#').rstrip().split("\t")
    variants = {}
    for line in lines:
        vals = line.rstrip().split("\t")
        d = dict(zip(fields, vals))
        if d['QUAL'][0].isdigit():
            pos = int(d['POS'])
            variants[pos] = d
    logfh.write(" {} variants\n".format(len(variants)))
    # special case for counting (TG)9-9T
#    if 87842 in variants and 87823 not in variants:
#        variants[87823] = {'CHROM':'CFTR', 'POS':'87823', 
#                           'REF':'ATG', 'ALT':'A' }
#    if 87844 in variants and 87845 not in variants:
#        variants[87845] = {'CHROM':'CFTR', 'POS':'87845', 
#                           'REF':'G', 'ALT':'T' }
    # special case for counting (TG)13-5T
#    if 87846 in variants and 87823 not in variants:
#        variants[87823] = {'CHROM':'CFTR', 'POS':'87823', 
#                           'REF':'A', 'ALT':'ATG' }
#    if 87847 in variants and variants[87847]['REF']=='T' and\
#        variants[87847]['ALT']=='TG':
#        variants[87847]['ALT'] += ',G'
    return (fields, variants)

def get_outlabel(filename, ext, outdir):
    label = os.path.basename(filename).rstrip(ext).rstrip('.')
    if outdir:
        label = os.path.join(outdir, label)
    return label

def get_alt_patt(ref, alt, pos, logfh, debug):
    patt = alt
    reflen = len(ref)
    altlen = len(alt)
    if reflen==altlen: # substitution
        if altlen==1:
            return patt
        else:
            patt = ''
            for i,b in enumerate(ref):
                if not b==alt[i]: patt += alt[i]
            if len(patt)>1:
                logfh.write("What do I do?  pos{}, ref{}, alt{}\n".\
                                 format(pos, ref, alt))
                return 'xx'
    elif reflen>altlen: # deletion
        num = reflen - altlen
        if ref.startswith(alt):
            delseq = ref[altlen:]
        else:
            delseq = ''
            offset = 0
            try:
                for i,b in enumerate(ref):
                    if not b==alt[i-offset]: 
                        delseq += b
                        offset += 1
            except:
                logfh.write("Deletion problem at {}: ref={}, alt={}\n".\
                            format(pos, ref, alt))
                sys.stderr.write("Deletion problem at {}: ref={}, alt={}\n".\
                            format(pos, ref, alt))
                return 'xxdel'
        if not len(delseq)==num:
            logfh.write("Deletion problem at {}: ref={}, alt={}\n".\
                            format(pos, ref, alt))
            return 'xxdelnum'
        patt = "-{}{}".format(num, delseq)
        if debug:
            logfh.write("{} DEL ref={} alt={} patt={}\n".format(pos,
                             ref, alt, patt))
    else: # reflen < altlen; insertion
        num = altlen - reflen
        if alt.startswith(ref):
            insseq = alt[reflen:]
        else:
            insseq = ''
            offset = 0
            try:
                for i,b in enumerate(alt):
                    if not b==ref[i-offset]: 
                        insseq += b
                        offset += 1
            except:
                logfh.write("Insertion problem at {}: ref={}, alt={}\n".\
                            format(pos, ref, alt))
                sys.stderr.write("Insertion problem at {}: ref={}, alt={}\n".\
                            format(pos, ref, alt))
                return 'xxins'
        if not len(insseq)==num:
            logfh.write("Insertion problem at {}: ref={}, alt={}\n".\
                            format(pos, ref, alt))
            return 'xxinsnum'
        patt = "+{}{}".format(num, insseq)
        if debug:
            logfh.write("{} INS ref={} alt={} patt={}\n".format(pos, 
                             ref, alt, patt))
    return patt

def split_pile_bases(bases_str):
    read_bases = []
    bases = [ r for r in re.split(r'([+-]|\d+|.)', bases_str) 
                    if len(r)>0 ]
    while bases:
        b = bases.pop(0)
        if b in '-+':
            for i in range(int(bases[0])+1):
                b += bases.pop(0)
        elif b=='^':
            bases.pop(0)
            continue
        elif b == '$':
            continue
        read_bases.append(b)
    return read_bases

def get_depths_from_pileline(variant, pileinfo, pos, logfh, debug):
    basealigns = split_pile_bases(pileinfo.bases)
    bases = [ b for b in basealigns if not b[0] in '-+']
    depths = { 'read_depth':len(bases), }
    altstr = ''
    alts = variant['ALT'].upper().split(',')
    altcounts = defaultdict(dict)
    for alt in alts:
        indelflag = False
        if alt=='<NON_REF>': 
            altpattF = '.'
            altpattR = ','
        else:
            altpatt = get_alt_patt(variant['REF'], alt, pos, logfh, debug)
            altpattF = altpatt.upper()
            altpattR = altpatt.lower()
            if len(altpatt)>1: indelflag = True
        altcounts['fwd'][alt] = basealigns.count(altpattF)
        altcounts['rev'][alt] = basealigns.count(altpattR)
        altcounts['tot'][alt] = altcounts['fwd'][alt] + altcounts['rev'][alt]
        alt_freq = float(altcounts['tot'][alt])/depths['read_depth']\
                                 if depths['read_depth'] > 0 else 0
        altcounts['freq'][alt] = "{:.3f}".format(alt_freq)
        altstr += ' {}={}/{}:{}/{} '.format(alt, altpattF, altpattR,
                            altcounts['fwd'][alt], altcounts['rev'][alt])
        # If indel alt count is low, then probably wrong so don't report
        # (this is because gatk will do local realignment)
        if indelflag and alt_freq<0.20: 
            altcounts['freq'][alt] = '-'
            altcounts['tot'][alt] = '-'
            altcounts['fwd'][alt] = '-'
            altcounts['rev'][alt] = '-'
            if debug:
                logfh.write("indelflag {}\t{}/{} = {}\n".format(indelflag, 
                                 variant['REF'], alt, alt_freq))
    depths['alt_depth'] = '/'.join([ str(altcounts['tot'][alt]) \
                                       for alt in alts ])
    depths['alt_depth_fwd'] = '/'.join([ str(altcounts['fwd'][alt]) \
                                       for alt in alts ])
    depths['alt_depth_rev'] = '/'.join([ str(altcounts['rev'][alt]) \
                                       for alt in alts ])
    depths['alt_AF'] = '/'.join([ altcounts['freq'][alt] for alt in alts ])
    if debug:
        logfh.write("{} tot {}\tref {}\t{}\n".format(pos, 
                    depths['read_depth'], variant['REF'], altstr))
    return depths

class PileLine:
    def __init__(self, pileline):
        d = pileline.split("\t")
        self.pos = int(d[1])
        self.ref = d[2]
        self.num_reads = int(d[3])
        self.bases = d[4] if len(d)>4 else ''

def add_fwd_rev_depths(depths, variants, bamfile, args):
    if not variants: return depths
    bampysam = pysam.AlignmentFile(bamfile, "rb")
    varpos = sorted(variants.keys())
    stepper_val = "all" if args.allreads else "samtools"
    bampileup = bampysam.pileup('CFTR', varpos[0], varpos[-1],
                                stepper=stepper_val, max_depth=99999)
    for pos in varpos:
        for pileupcol in bampileup:
            pilepos = pileupcol.reference_pos+1 # convert to 1-based pos
            if pilepos == pos:
                read_depth = pileupcol.nsegments
                numrev = len([ r for r in pileupcol.pileups \
                               if r.alignment.is_reverse ])
                depths[pos].update({ 'read_depth':read_depth,
                                'read_depth_rev':numrev,
                                'read_depth_fwd':read_depth-numrev, })
                break
    return depths

def find_variant_depths(variants, bamfile, logfh, args):
    depths = {}
    if not variants: return depths
    logfh.write("Creating pileup from bamfile: {}\n".format(bamfile))
    pysam_pileup = pysam.mpileup("-B", "-d99999", '-Q0', 
                                 "-f", RESOURCE['ref_fa'], bamfile)
    varpos = sorted(variants.keys())
    for pos in varpos:
        variant = variants[pos]
        for pileline in pysam_pileup:
            pileinfo = PileLine(pileline)
            if pileinfo.pos == pos:
                depths[pos] = get_depths_from_pileline(variant, pileinfo,
                                                  pos, logfh, args.debug)
                break
    return depths

def print_depths(outfile, variants, depths, vcffields, logfh):
    fields = vcffields[:5] 
    fields2 = ['read_depth', 'read_depth_fwd', 'read_depth_rev', 'alt_AF', 
               'alt_depth']
    rows = []
    for pos in sorted(variants.keys()):
        v = variants[pos]
        row = [ v[f] if f in v else '.' for f in fields ]
        pos = int(v['POS'])
        num_reads = 0
        if pos in depths:
            row.extend([str(depths[pos][f]) if f in depths[pos] else '' \
                       for f in fields2 ])
        else:
            row.extend(['']*len(fields2))
        rows.append(row)
    logfh.write("Writing {}\n".format(outfile))
    with open(outfile, 'w') as ofh:
        ofh.write("\t".join(fields + fields2) + "\n")
        for row in rows:
            ofh.write("\t".join(row) + "\n")

def calculate_depths(label, infiles, args):
    bamfile = infiles[label]['bam']
    vcffile = infiles[label]['vcf']
    outlabel = get_outlabel(bamfile, 'bam', args.outdir)
    loglabel = get_outlabel(bamfile, 'bam', args.logdir) if args.logdir \
               else outlabel
    logfile = loglabel + ".calculate_depths.log"
    outfile = outlabel + '.depths{}.txt'.format(args.tag)
    sys.stderr.write("  Calculating depths for {}\n".format(outlabel))
    sys.stderr.write("    Log file: {}\n".format(logfile))
    logfh = open(logfile, 'w')
    logfh.write("Bam file: {}\n".format(bamfile))
    logfh.write("VCF file: {}\n\n".format(vcffile))
    try:
        if have_file(outfile, args.force, stderr=logfh):
            logfh.write("  Already have {}\n".format(outfile))
        else:
            logfh.write("Start time: {}\n".format(timestamp()))
            (vcffields, variants) = parse_vcf(vcffile, logfh)
            depths = find_variant_depths(variants, bamfile, logfh, args)
            add_fwd_rev_depths(depths, variants, bamfile, args)
            print_depths(outfile, variants, depths, vcffields, logfh)
            logfh.write("End time: {}\n".format(timestamp()))
    except Exception, e:
        e.args += (vcffile, )
        raise
    finally:
        logfh.close()
    if os.path.isfile(outfile):
        sys.stderr.write("    Depth file {}\n".format(outfile))
    else:
        sys.stderr.write("    Failed writing depth file {}\n".format(outfile))
    return outfile


if __name__ == '__main__':
    descr = "Count read depths at variant positions."
    parser = ArgumentParser(description=descr)
    parser.add_argument("inputfiles", nargs="+", 
                        help="Pairs of bam and vcf files")
    parser.add_argument("-t", "--tag", default='', 
                        help="Tag to add to output file names.")
    parser.add_argument("-o", "--outdir", help="Directory for output files.")
    parser.add_argument("-a", "--allreads", default=False, 
                        action="store_true", help="Count anomalous reads.")
    parser.add_argument("-f", "--force", default=False, action="store_true",
                        help="Overwrite existing files.")
    parser.add_argument("-p", "--processes", default=1, type=int,
                        help="Number of processes to run.")
    parser.add_argument("--debug", default=False, action="store_true",
                        help="Debug mode.")
    parser.add_argument("--logdir", help="Directory for log files.")

    if len(sys.argv)<2:
        parser.print_help()
        sys.exit()
    args = parser.parse_args()

    infiles = organize_input_files(args.inputfiles)
    sys.stderr.write("\nCalculating depths for {} samples using {}.\n".format(
           len(infiles), "{} process{}".format(args.processes, 'es' if 
           args.processes>1 else '')))
    pool = mp.Pool(processes=args.processes)
    results = [ pool.apply_async(calculate_depths, args=[label, infiles, 
                args]) for label in sorted(infiles.keys()) ]
    files = [p.get() for p in results ]


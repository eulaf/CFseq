#!/usr/bin/env python

"""
Convert sample VCFs to spreadsheet format with useful information in 
new columns.
"""

import sys
import os
import re
from collections import defaultdict
from argparse import ArgumentParser

base_dir = os.path.dirname(__file__) or '.'
lib_dir = os.path.abspath(os.path.join(base_dir, 'lib_cftr'))
sys.path.insert(0, lib_dir)

from common.misc import open_file, have_file
import cftr

def parse_vcf(vcffile):
    """Parse variants in joint vcf file."""
    sys.stderr.write("VCF file {} -".format(vcffile))
    fh = open_file(vcffile)
    lines = fh.readlines()
    fh.close()
    while lines: # skip vcf header
        if not lines[0].startswith('##'):
            break
        lines.pop(0)
    fields = lines.pop(0).lstrip('#').rstrip().split("\t")
    variants = []
    for line in lines:
        vals = line.rstrip().split("\t")
        d = dict(zip(fields, vals))
        variants.append(d)
    sys.stderr.write(" {} variants\n".format(len(variants)))
    return (fields, variants)

def create_pos_key(chrom, pos):
    poskey = "{}:{:>9s}".format(chrom, str(pos))
    return poskey

def calculate_AF(idata):
    dp =  idata['DP'] if 'DP' in idata else ''
    saflist =  idata['SAF'].split(',')
    sarlist =  idata['SAR'] .split(',')
    aflist = []
    for i, saf in enumerate(saflist):
        ad = int(saf) + int(sarlist[i])
        af = "{:.3f}".format(ad/float(dp)) if dp else '-'
        aflist.append(af)
    return ('/'.join(aflist), dp)

def add_avg_base_quality(sample_row, idata):
    quallist = []
    if 'AO' in idata and 'QA' in idata:
        aolist = idata['AO'].split(',')
        qalist = idata['QA'].split(',')
        qaaolist = []
        for (ao, qa) in zip(aolist, qalist):
            if float(ao)>0:
                qaao = float(qa)/float(ao)
                quallist.append(qaao)
                qaaolist.append("{:.3f}".format(qaao))
        sample_row['QA/AO'] = '/'.join(qaaolist)
    else: sample_row['QA/AO'] = ""
    if quallist:
        sample_row['MaxAvgAltQual'] = "{:.3f}".format(max(quallist))

def flatten_vcf_data(fields, vcfdata, flatdata=defaultdict(dict)):
    infofields = fields[:9]
    samples = fields[9:]
    ifields = [] #'DP', 'SRF', 'SRR', 'SAF', 'SAR','ABP', 'QA', 'AO']
    for d in vcfdata:
        if d['CHROM']=='chr7':
            (d['CHROM'], pos) = cftr.hg19_to_CFTR(d['CHROM'], d['POS'])
            d['POS'] = str(pos)
        poskey = create_pos_key(d['CHROM'], d['POS'])
        ids = d['ID'].split(';')
        dbsnp = [ rs for rs in ids if rs.startswith('rs') ]
        d['hg19_coordinates'] = "{}:{}".format(*cftr.CFTR_to_hg19(
                                d['CHROM'], d['POS']))
        d['ID'] = ';'.join(dbsnp) if dbsnp else '-'
        d['hg19_ID'] = '-' 
        for sample in samples:
            sfields = d['FORMAT'].split(':')
            sdata = dict(zip(sfields, d[sample].split(':')))
            idata = dict([ v.split('=') for v in d['INFO'].split(';') \
                           if '=' in v ])
            if sdata['GT']=='0/0' or sdata['GT']=='./.':
                continue
            sample_row = d.copy()
            sample_row['Sample ID'] = sample
            gtvals = [int(i) for i in sdata['GT'].split('/')]
            sample_row['hom/het'] = 'hom' if gtvals[0]==gtvals[1] else 'het'
            (sample_row['alt_AF'], sample_row['DP']) = calculate_AF(idata)
#            for f in ifields: sample_row[f] = idata[f]
            add_avg_base_quality(sample_row, idata)
            sample_row['format_values'] = ':'.join([ sdata[f] for f in sfields ])
            flatdata[sample][poskey] = sample_row
    ifields.extend(['QA/AO',])
    newfields = ['Sample ID',]+infofields+['format_values', 'hom/het', 
                 'ROI (Y/N)', 'c./p.(AnnoVar)', 'DP', 'alt_AF',]+\
                 ifields + ['hg19_coordinates', ]
    return (flatdata, newfields)

def get_annovar_data(annovar_files):
    if not annovar_files: return {}
    sys.stderr.write("\nParsing annovar data {}\n".format(annovar_files))
    afiles = []
    for ext in ('.variant_function', '.exonic_variant_function'):
        files = [ f for f in annovar_files if f.endswith(ext) ]
        if files:
            afiles.append(files[0])
            if not os.path.isfile(files[0]):
                sys.stderr.write("  Not a valid file: {}\n".format(files[0]))
                sys.exit(1)
        else:
            sys.stderr.write("  Could not find annovar {} file\n".format(ext))
            sys.exit(1)
    sys.stderr.write("  Reading {}\n".format(afiles[0]))
    with open(afiles[0], 'r') as fh: #variant_function file
        data_rows = [ [l.rstrip().split("\t"),] for l in fh.readlines() \
                      if len(l)>1 ]
    sys.stderr.write("    Read {} lines\n".format(len(data_rows)))
    sys.stderr.write("  Reading {}\n".format(afiles[1]))
    with open(afiles[1], 'r') as fh: #exonic_variant_function file
        numdat = 0
        for l in fh.readlines():
            d = l.rstrip().split("\t")
            if not d: continue
            linenum = int(d.pop(0).lstrip('line'))
            data_rows[linenum-1].append(d)
            numdat += 1
        sys.stderr.write("    Read {} lines\n".format(numdat))
    annovar_data = defaultdict(dict)
    for row in data_rows:
        (chrom_hg19, pos_hg19) = row[0][2:4]
        (chrom, pos) = cftr.hg19_to_CFTR(chrom_hg19, pos_hg19)
        poskey = create_pos_key(chrom, pos)
        if len(row)>1: #has exonic_variant_function data
            annovar_data[poskey][row[1][1].rstrip(',')] = row[1][0]
        elif ':' in row[0][1]:
            annovar_data[poskey][row[0][1]] = row[0][0]
    return annovar_data

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
                sys.stderr.write("pos{} not in {}-{}\n".format(pos, s, e))
                break
    return roi_flag
   
def create_spreadsheet(fields, flatdata, annovar_data, roi, outfile,
                       rejectfile, args):
    sys.stderr.write("\nWriting {}\n".format(outfile))
    numlines = 0
    filtered = defaultdict(list)
    afpatt = re.compile('[,0-9\.]+$')
    with open(outfile, 'w') as ofh:
        ofh.write("#"+"\t".join(fields)+"\n")
        for sample in sorted(flatdata.keys()):
            for poskey in sorted(flatdata[sample].keys()):
                d = flatdata[sample][poskey]
                if poskey in annovar_data:
                    d['c./p.(AnnoVar)'] = ';'.join(annovar_data[poskey].keys())
                d['ROI (Y/N)'] = cftr.in_roi(d['CHROM'], d['POS'], roi)
                row = [ d[f] if f in d else '-' for f in fields ]
                ofh.write("\t".join(row) + "\n")
                numlines += 1
    sys.stderr.write("    Wrote {} data lines\n".format(numlines))
    if filtered:
        sys.stderr.write("Putting rejects in {}\n".format(rejectfile))
        numreject = 0
        with open(rejectfile, 'w') as ofh:
            ofh.write("\t".join(fields)+"\n")
            for fkey in filtered.keys():
                sys.stderr.write("Filtered {} variants for {}.\n".format(
                len(filtered[fkey]), fkey))
                for d in filtered[fkey]:
                    row = [ d[f] if f in d else '-' for f in fields ]
                    ofh.write("\t".join(row) + "\n")
                    numreject += 1
        sys.stderr.write("    Wrote {} rejected variants\n".format(numreject))
    else:
        sys.stderr.write("No filtered variants.\n")


def get_outlabel(vcffiles, args):
    label = args.label if args.label else os.path.basename(vcffiles[0])\
            .rstrip('vcf').rstrip('.')
    if args.outdir:
        if os.path.isdir(args.outdir):
            label = os.path.join(args.outdir, label)
        else:
            msg = "ERROR: Outdir {} does not exist.\n".format(args.outdir)
            sys.stderr.write(msg)
            sys.exit(1)
    return label

if __name__ == '__main__':
    descr = """Convert freebayes' sample VCFs into one friendly 
            spreadsheet format."""
    parser = ArgumentParser(description=descr)
    parser.add_argument("vcffiles", nargs="+", 
                        help="VCF files to convert")
    parser.add_argument("-a", "--annovar", action="append",
                        help="Add annovar data from joint VCF.  Specify" +\
                        " .variant_function and" +\
                        " .exonic_variant_function files.", )
    parser.add_argument("-o", "--outdir", help="Directory for output file.",)
    parser.add_argument("-l", "--label", help="Label for output file.",)
    parser.add_argument("-f", "--force", default=False, action='store_true',
                        help="Overwrite existing files.")

    if len(sys.argv)<2:
        parser.print_help()
        sys.exit()
    args = parser.parse_args()

    outlabel = get_outlabel(args.vcffiles, args)
    outfile = outlabel + ".results.txt"
    rejectfile = outlabel + ".rejects.txt"
    if have_file(outfile, args.force):
        sys.stderr.write("  Already have {}.\n".format(outfile))
    else:
        flatdata = defaultdict(dict)
        for vcffile in args.vcffiles:
            (fields, vcfdata) = parse_vcf(vcffile)
            (flatdata, newfields) = flatten_vcf_data(fields, vcfdata, flatdata)
        annovar_data = get_annovar_data(args.annovar)
        bedfile = cftr.RESOURCE['analysis_roi_bed']
        if not os.path.isfile(bedfile): 
            sys.exit("BED file {} not found\n".format(bedfile))
        roi = cftr.parse_bedfile(bedfile)
        create_spreadsheet(newfields, flatdata, annovar_data, roi, outfile,
                           rejectfile, args)


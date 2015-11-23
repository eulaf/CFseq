#!/usr/bin/env python

"""
Filter freebayes vcf file(s).
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

def create_pos_key(chrom, pos):
    poskey = "{}:{:>9s}".format(chrom, pos)
    return poskey

def create_ref_alt_key(ref, alt):
    return "{}:{}".format(ref, alt)

def parse_dbsnp():
    snpfile = cftr.RESOURCE['dbsnp_edited']
    dbsnp = defaultdict(dict)
    with open(snpfile, 'r') as fh:
        for line in fh.readlines():
            if line.startswith('#'): continue
            else:
                vals = line.split("\t")
                poskey = create_pos_key(vals[0], vals[1])
                for alt in vals[4].split(','):
                    ref_alt = create_ref_alt_key(vals[3], alt)
                    dbsnp[poskey][ref_alt] = vals[2]
    sys.stderr.write("dbSNP: {} positions\n".format(len(dbsnp.keys())))
    return dbsnp

def parse_vcf(vcffile):
    fh = open_file(vcffile)
    lines = fh.readlines()
    fh.close()
    header = []
    while lines: # vcf header
        if lines[0].startswith('#'):
            header.append(lines.pop(0))
        else:
            break
    fields = header[-1].lstrip('#').rstrip().split("\t")
    vcfinfo = {}
    for line in lines:
        vals = line.rstrip().split("\t")
        d = dict(zip(fields, vals))
        poskey = create_pos_key(d['CHROM'], d['POS'])
        if poskey in vcfinfo:
            sys.stderr.write("  Variant pos {} duplicated\n".format(poskey))
        vcfinfo[poskey] = d
    return (header, fields, vcfinfo)

def add_values_to_vcf_cell(d, field, newvalstr, prepend=False):
    if d[field] and len(d[field]) > 1:
        # Already have values.  Check to make sure new values are not already
        # in list before including.
        new_val_list = newvalstr.split(';')
        if prepend: new_val_list.reverse()
        for newval in new_val_list:
            if newval not in d[field]:
                if prepend:
                    d[field] = newval + ";" + d[field]
                else:
                    d[field] += ";" + newval
    else:
        d[field] = newvalstr

def add_dbsnp(vcfinfo, dbsnp):
    numadded = 0
    for poskey in vcfinfo.keys():
        if poskey in dbsnp:
            v = vcfinfo[poskey]
            ids = []
            for alt in v['ALT'].split(','):
                ref_alt = create_ref_alt_key(v['REF'], alt)
                if ref_alt in dbsnp[poskey]:
                    ids.append(dbsnp[poskey][ref_alt])
                elif len(alt)>1 and alt[-1]==v['REF'][-1]:
                    ref_alt = create_ref_alt_key(v['REF'][:-1], alt[:-1])
                    if ref_alt in dbsnp[poskey]:
                        ids.append(dbsnp[poskey][ref_alt])
                    else:
                        ids.append('.')
                else:
                    ids.append('.')
            if [ rs for rs in ids if rs.startswith('rs') ]:
                add_values_to_vcf_cell(v, 'ID', ",".join(ids))
                numadded += 1
    sys.stderr.write("  Added {} snp IDs\n".format(numadded))

def max_avg_alt_base_quality(idata):
    quallist = [0]
    if 'AO' in idata and 'QA' in idata:
        aolist = idata['AO'].split(',')
        qalist = idata['QA'].split(',')
        qaaolist = []
        for (ao, qa) in zip(aolist, qalist):
            if float(ao)>0:
                qaao = float(qa)/float(ao)
                quallist.append(qaao)
    return max(quallist)

def is_float(val):
    try:
        float(val)
    except ValueError:
        return False
    return True

def add_filter_info_line(args, header):
    filterinfo = "##" + os.path.basename(sys.argv[0])
    filteropts = []
    if args.qual:
        filteropts.append("QUAL < {}".format(args.qual))
    if args.dp:
        filteropts.append("DP < {}".format(args.dp))
    if args.altbasequal:
        filteropts.append("QA/AO < {}".format(args.altbasequal))
    if filteropts:
        filterinfo += ": remove "+" or ".join(filteropts)
    if args.dbsnp:
        filterinfo += "; annotating dbsnp"
    for i, line in enumerate(header):
        if line.startswith('##INFO'):
            header.insert(i, filterinfo + "\n")
            break

def filter_vcf(header, fields, vcfinfo, outvcf, dbsnp, args):
    numlines = 0
    add_filter_info_line(args, header)
    pos_removed = []
    with open(outvcf, 'w') as ofh:
        ofh.write(''.join(header))
        for poskey in sorted(vcfinfo.keys()):
            d = vcfinfo[poskey]
            idata = dict([ v.split('=') for v in d['INFO'].split(';') \
                           if '=' in v ])
            if args.qual and 'QUAL' in d and is_float(d['QUAL']):
                if float(d['QUAL']) < args.qual: 
                    pos_removed.append(d['POS'])
                    continue
            if args.dp and 'DP' in idata and idata['DP'].isdigit():
                if int(idata['DP']) < args.dp: 
                    pos_removed.append(d['POS'])
                    continue
            if args.altbasequal:
                abq = max_avg_alt_base_quality(idata)
                if abq < args.altbasequal: 
                    pos_removed.append(d['POS'])
                    continue
            row = [ d[f] if f in d else '.' for f in fields ]
            ofh.write("\t".join(row) + "\n")
            numlines += 1
    sys.stderr.write("  Wrote {} variant lines. Skipped {}.\n".format(
                     numlines, len(pos_removed)))
    if pos_removed: 
        sys.stderr.write("  Skipped: {}\n".format(", ".join(pos_removed)))
    return pos_removed

if __name__ == '__main__':
    descr = """Filter freebayes VCF files."""
    parser = ArgumentParser(description=descr)
    parser.add_argument("vcffiles", nargs="+", 
                        help="VCF files to convert")
    parser.add_argument("-o", "--outdir", help="Directory for output files.",)
    parser.add_argument("-a", "--altbasequal", type=int, 
                        help="Min alt base quality to pass.",)
    parser.add_argument("-d", "--dp", type=int, 
                        help="Min DP to pass.",)
    parser.add_argument("-q", "--qual", type=int,
                        help="Min QUAL to pass.",)
    parser.add_argument("-s", "--dbsnp", default=False, action='store_true',
                        help="Add dbsnp to ID field.",)
    parser.add_argument("-l", "--label", default='filtered',
                        help="Label for output files.",)
    parser.add_argument("-f", "--force", default=False, action='store_true',
                        help="Overwrite existing files.")

    if len(sys.argv)<2:
        parser.print_help()
        sys.exit()
    args = parser.parse_args()

    dbsnp = parse_dbsnp() if args.dbsnp else {}
    pos_removed = []
    for vcffile in args.vcffiles:
        sys.stderr.write("Processing {}.\n".format(vcffile))
        outlabel = os.path.basename(vcffile).rstrip('vcf').rstrip('.')
        outvcf = outlabel + '.{}.vcf'.format(args.label)
        if args.outdir:
            outvcf = os.path.join(args.outdir, outvcf)
        if have_file(outvcf, args.force):
            sys.stderr.write("  Already have {}.\n".format(outvcf))
        else:
            sys.stderr.write("Writing {}\n".format(outvcf))
            (header, fields, vcfinfo) = parse_vcf(vcffile)
            if args.dbsnp: add_dbsnp(vcfinfo, dbsnp)
            pos_rm = filter_vcf(header, fields, vcfinfo, outvcf, dbsnp, args)
            pos_removed.extend(pos_rm)
    sys.stderr.write("\nPositions removed: {}\n".format(", ".join(
        [str(j) for j in sorted([int(i) for i in set(pos_removed)])])))

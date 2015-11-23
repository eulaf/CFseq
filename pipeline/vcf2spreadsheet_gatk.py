#!/usr/bin/env python

"""
Convert joint sample VCF to spreadsheet format with useful information in 
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

def DPcalculate_AF(sdata, infofield):
    if not 'DP' in sdata:
        dpinfo = [ f for f in infofield.split(';') if 'DP=' in f ]
        if dpinfo:
            sdata['DP'] = dpinfo[0].split('=')[1]
            print "DP="+sdata['DP']
        else:
            sdata['DP'] = ''
    dp = int(sdata['DP']) if sdata['DP'].isdigit() else ''
    aflist = []
    if 'FREQ' in sdata: # varscan
        aflist.append(str(float(sdata['FREQ'].rstrip('%'))/100))
    else:
        adlist = sdata['AD'].split(',')[1:]
        for ad in adlist:
            if dp and ad.isdigit():
                af = "{:.3f}".format(int(ad)*1.0/dp)
                aflist.append(af)
            else:
                aflist.append('-')
    dp = '-' if not type(dp)==int else str(dp)
    return (dp, ','.join(aflist))


def flatten_vcf_data(fields, vcfdata, flatdata=defaultdict(dict)):
    infofields = fields[:9]
    samples = fields[9:]
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
            altchoices = [d['REF'],] + d['ALT'].split(',')
            altvals = set([ altchoices[i] for i in gtvals if i>0 ])
            sample_row['ALT'] = ','.join(altvals)
            sample_row['hom/het'] = 'hom' if gtvals[0]==gtvals[1] else 'het'
            (sample_row['DP'], sample_row['AF']) = DPcalculate_AF(sdata,
                                                   d['INFO'])
            sample_row['QD'] = idata['QD'] if 'QD' in idata else ''
            sample_row['FS'] = idata['FS'] if 'FS' in idata else ''
            newgtvals = []
            if gtvals[0]==0:
                sdata['GT'] = '0/1'
            elif gtvals[0]==gtvals[1]:
                sdata['GT'] = '1/1'
            else:
                sdata['GT'] = '1/2'
            sample_row['format_values'] = ':'.join([ sdata[f] for f in sfields ])
            flatdata[sample][poskey] = sample_row
    newfields = ['Sample ID',]+infofields+['format_values', 'hom/het', 
                 'ROI (Y/N)', 'c./p.(AnnoVar)', 'DP', 'AF', 'QD', 'FS', 
                 'hg19_coordinates', ]
    return (flatdata, newfields)

def add_depth_info(depthdir, flatdata, fieldlist, args):
    if not depthdir: return []
    depthfiles = [ f for f in os.listdir(depthdir) if \
                   f.endswith('.depths.txt') ]
    newfields = ['read_depth', 'read_depth_fwd', 'read_depth_rev', 'alt_AF',
                 'alt_depth']
    altfields = [f for f in newfields if f.startswith('alt')]
    i = fieldlist.index('hg19_coordinates')
    fieldlist[i:i] = newfields
    for sample in sorted(flatdata.keys()):
        dfile = get_depth_file_for_sample(sample, depthfiles)
        if not dfile: continue
        with open(os.path.join(depthdir, dfile), 'r') as fh:
            fields = fh.readline().lstrip('#').rstrip().split("\t")
            sampledepths = {}
            for line in fh.readlines():
                vals = line.rstrip().split("\t")
                d = dict(zip(fields, vals))
                poskey = create_pos_key(d['CHROM'], d['POS'])
                sampledepths[poskey] = d
            for poskey in sorted(flatdata[sample].keys()):
                if poskey in sampledepths:
                    for f in newfields:
                        flatdata[sample][poskey][f] = sampledepths[poskey][f]
                else:
                    for f in newfields:
                        flatdata[sample][poskey][f] = '-'

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
    with open(afiles[0], 'r') as fh:
        data_rows = [ [l.rstrip().split("\t"),] for l in fh.readlines() \
                      if len(l)>1 ]
    sys.stderr.write("    Read {} lines\n".format(len(data_rows)))
    sys.stderr.write("  Reading {}\n".format(afiles[1]))
    with open(afiles[1], 'r') as fh:
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
        if len(row[0])>10:
            (chrom_hg19, pos_hg19) = row[0][10:12]
        else:
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
   
def get_depth_file_for_sample(sample, depthfiles, debug=False):
    dfile = [ f for f in depthfiles if f.startswith(sample) ]
    if debug: dfile = [ f for f in dfile if 'debug' in f ]
    else: dfile = [ f for f in dfile if not 'debug' in f ]
    if dfile:
        if len(dfile)==1:
            return dfile[0]
        else:
            sys.stderr.write("Too many matches to {}: {}\n".format(sample,
                             dfile))
    else:
        sys.stderr.write("No depth file for {}\n".format(sample))
    return None

def is_float(val):
    try:
        float(val)
    except ValueError:
        return False
    return True

def af_less_than_minfreq(af_str, minfreq):
    aflist = af_str.split(',')
    is_less = True
    for af in aflist:
        if is_float(af) and float(af) >= minfreq:
            is_less = False
    return is_less

def filter_info(args):
    filterinfo = []
    if args.maxfs: filterinfo.append('FS>{}'.format(args.maxfs))
    if args.minqd: filterinfo.append('QD<{}'.format(args.minqd))
    if args.mindp: filterinfo.append('DP<{}'.format(args.mindp))
    if args.depths and args.mindepth: 
        filterinfo.append('read_depth<{}'.format(args.mindepth))
    if args.minfreq: filterinfo.append('alt_AF<{}'.format(args.minfreq))
    return ' || '.join(filterinfo)

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

def create_spreadsheet(fields, flatdata, annovar_data, roi, outfile,
                       rejectfile, args):
    sys.stderr.write("\nWriting {}\n".format(outfile))
    numlines = 0
    mindepth = args.mindepth if (args.depths and args.mindepth) else 0
    filtered = defaultdict(list)
    afpatt = re.compile('[,0-9\.]+$')
    with open(outfile, 'w') as ofh:
        ofh.write("#Filtering: {}\n".format(filter_info(args)))
        ofh.write("#"+"\t".join(fields)+"\n")
        for sample in sorted(flatdata.keys()):
            for poskey in sorted(flatdata[sample].keys()):
                d = flatdata[sample][poskey]
                if poskey in annovar_data:
                    d['c./p.(AnnoVar)'] = ';'.join(annovar_data[poskey].keys())
                d['ROI (Y/N)'] = cftr.in_roi(d['CHROM'], d['POS'], roi)
                if not args.nofilter and d['FILTER'] not in ('PASS','.'):
                    filtered['FILTER'].append(d)
                    continue
                if args.maxfs and ('FS' in d) and is_float(d['FS']):
                    if (float(d['FS']) > args.maxfs):
                        add_values_to_vcf_cell(d, 'FILTER', 'highFS')
                        filtered['highFS'].append(d)
                        continue
                if args.minqd and 'QD' in d and is_float(d['QD']):
                    if float(d['QD']) < args.minqd:
                        add_values_to_vcf_cell(d, 'FILTER', 'lowQD')
                        filtered['lowQD'].append(d)
                        continue
                if args.mindp and 'DP' in d and d['DP'].isdigit():
                    if int(d['DP']) < args.mindp: 
                        add_values_to_vcf_cell(d, 'FILTER', 'lowDP')
                        filtered['lowDP'].append(d)
                        continue
                if 'alt_AF' in d and afpatt.match(d['alt_AF']):
                    if af_less_than_minfreq(d['alt_AF'], args.minfreq):
                        add_values_to_vcf_cell(d, 'FILTER', 'low_alt_AF')
                        filtered['alt_AF'].append(d)
                        continue
                elif af_less_than_minfreq(d['AF'], args.minfreq):
                    add_values_to_vcf_cell(d, 'FILTER', 'lowAF')
                    filtered['AF'].append(d)
                    continue
                if 'read_depth' in d and int(d['read_depth']) < mindepth:
                    add_values_to_vcf_cell(d, 'FILTER', 'low_read_depth')
                    filtered['read_depth'].append(d)
                    continue
                row = [ d[f] if f in d else '-' for f in fields ]
                ofh.write("\t".join(row) + "\n")
                numlines += 1
    sys.stderr.write("    Wrote {} data lines\n".format(numlines))
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
    descr = """Convert joint sample VCF and/or multiple VCFs into one friendly 
            spreadsheet format."""
    parser = ArgumentParser(description=descr)
    parser.add_argument("vcffiles", nargs="+", 
                        help="VCF files to convert")
    parser.add_argument("-a", "--annovar", action="append",
                        help="Add annovar data from joint VCF.  Specify" +\
                        " .variant_function and" +\
                        " .exonic_variant_function files.", )
    parser.add_argument("-d", "--depths", 
                        help="Directory to find depths.txt files",)
    parser.add_argument("-l", "--label", help="Label for output file.",)
    parser.add_argument("-m", "--mindepth", type=int, default=20,
                        help="Minimum coverage to report.",)
    parser.add_argument("-o", "--outdir", help="Directory for output file.",)
    parser.add_argument("-x", "--nofilter", default=False, action='store_true',
                        help="Include variants that do not pass filter.")
    parser.add_argument("-q", "--minfreq", type=int, default=0.2,
                        help="Minimum AF to report.",)
    parser.add_argument("--mindp", type=int, default=5,
                        help="Minimum DP to report (default 5)",)
    parser.add_argument("--minqd", type=float, default=2.0,
                        help="Minimum QD to report (default 2.0).",)
    parser.add_argument("--maxfs", type=float, default=30.0,
                        help="Maximum FS to report (default 30.0).",)
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
        add_depth_info(args.depths, flatdata, newfields, args)
        annovar_data = get_annovar_data(args.annovar)
        bedfile = cftr.RESOURCE['analysis_roi_bed']
        if not os.path.isfile(bedfile): 
            sys.exit("BED file {} not found\n".format(bedfile))
        roi = cftr.parse_bedfile(bedfile)
        create_spreadsheet(newfields, flatdata, annovar_data, roi, outfile,
                           rejectfile, args)


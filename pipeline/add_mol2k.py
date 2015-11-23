#!/usr/bin/env python

"""
Find mol2k variants in spreadsheet.  Writes a new spreadsheet with
column of mol2k variants found and a new mol2k spreadsheet with
'Seen' column updated.
"""

import sys
import os
from collections import defaultdict
from argparse import ArgumentParser

base_dir = os.path.dirname(__file__) or '.'
lib_dir = os.path.abspath(os.path.join(base_dir, 'lib_cftr'))
sys.path.insert(0, lib_dir)

from common.misc import have_file
import cftr

def create_pos_key(chrom, pos):
    if not chrom or not pos: return None
    poskey = "{}:{:>9s}".format(chrom, str(pos))
    return poskey

def create_ref_alt_key(ref, alt):
    return "{}:{}".format(ref, alt)

def parse_mol2k(mol2k_file, reset_seen):
    mol2k = defaultdict(dict)
    num = 0
    sys.stderr.write("\nParsing {}\n".format(mol2k_file))
    with open(mol2k_file, 'r') as fh:
        lines = fh.readlines()
    if len(lines)==1: lines = lines[0].split("\r")
    fields = lines.pop(0).lstrip('#').rstrip().split("\t")
    for line in lines:
        vals = line.rstrip().split("\t")
        d = dict(zip(fields, vals))
        poskey = create_pos_key(d['CHROM'], d['POS']) 
        refaltkey = create_ref_alt_key(d['REF'], d['ALT'])
        if not poskey: continue
        if not 'Seen' in d or reset_seen:
            d['Seen'] = 'N'
        hg19loc = cftr.CFTR_to_hg19(d['CHROM'], d['POS'])
        d['hg19_coordinates'] = "{}:{}".format(*hg19loc)
        mol2k[poskey][refaltkey] = d
        num += 1
    sys.stderr.write("  Got {} variants\n".format(num))
    if not 'Seen' in fields: fields.insert(fields.index('Unique_ID'), 'Seen')
    if not 'hg19_coordinates' in fields: fields.insert(
                               fields.index('Unique_ID'), 'hg19_coordinates')
    if not 'rsID' in fields: fields.insert(fields.index('Legacy_name')+1, 
                                           'rsID')
    return (mol2k, fields)

def find_mol2k_variants(spreadsheet, mol2k, outfile):
    num = 0
    numseen = 0
    sys.stderr.write("\nParsing {}\n".format(spreadsheet))
    sys.stderr.write("Writing {}\n".format(outfile))
    molpath_field = 'MolPath db (Y/N)'
    hg19id_field = 'hg19_ID'
    cdna_field = 'cDNA_name'
    legacy_field = 'Legacy_name'
    with open(spreadsheet, 'r') as fh, open(outfile, 'w') as ofh:
        lines = fh.readlines()
        while lines[0].startswith('#'):
            fields = lines.pop(0).lstrip('#').rstrip().split("\t")
        if not molpath_field in fields: fields.append(molpath_field)
        i = fields.index(molpath_field)
        if not hg19id_field in fields: fields.insert(i, hg19id_field)
        if not legacy_field in fields: fields.insert(i, legacy_field)
        if not cdna_field in fields: fields.insert(i, cdna_field)
        ofh.write("#"+"\t".join(fields)+"\n")
        for line in lines:
            vals = line.rstrip().split("\t")
            d = dict(zip(fields, vals))
            poskey = create_pos_key(d['CHROM'], d['POS']) 
            if not poskey: continue
            d[molpath_field] = 'N' 
            if poskey in mol2k:
                alts = d['ALT'].split(',')
                hg19ids = []
                cdna_names = []
                legacy_names = []
                molpath_yn = []
                for alt in alts:
                    ref = d['REF']
                    refaltkey = create_ref_alt_key(ref, alt)
                    while refaltkey not in mol2k[poskey] and len(ref)>1 \
                          and ref[-1]==alt[-1]:
                        ref = ref[:-1]
                        alt=alt[:-1]
                        refaltkey = create_ref_alt_key(ref, alt)
                    if refaltkey in mol2k[poskey]:
                        mol2k[poskey][refaltkey]['Seen'] = 'Y' 
                        if d['ID'] and d['ID'].startswith('rs'):
                            rsID = mol2k[poskey][refaltkey]['rsID']
                            if rsID.startswith('rs') and d['ID'] not in rsID:
                                mol2k[poskey][refaltkey]['rsID']+= ";"+d['ID']
                            else:
                                mol2k[poskey][refaltkey]['rsID'] = d['ID']
                        molpath_yn.append('Y')
                        cdna = mol2k[poskey][refaltkey]['cDNA_name']
                        cdna_names.append(cdna if cdna else '.')
                        lname = mol2k[poskey][refaltkey]['Legacy_name']
                        legacy_names.append(lname if lname else '.')
                        hg19id = mol2k[poskey][refaltkey]['Unique_ID']
                        hg19ids.append(hg19id if hg19id else '.')
                        comment = mol2k[poskey][refaltkey]['Comment-CS']
                        if len(comment)>1: 
                            molpath_yn[-1] += ' = '+ comment
                    else:
                        molpath_yn.append('N')
                        cdna_names.append('-')
                        legacy_names.append('-')
                        hg19ids.append('-')
                d[hg19id_field] = ';'.join(hg19ids)
                d[cdna_field] = ';'.join(cdna_names)
                d[legacy_field] = ';'.join(legacy_names)
                d[molpath_field] = ';'.join(molpath_yn)
                numseen += 1
            row = [ d[f] if f in d else '' for f in fields ]
            ofh.write("\t".join(row)+"\n")
            num += 1
    sys.stderr.write("  Checked {} variants\n".format(num))
    sys.stderr.write("  Found {} mol2k variants\n".format(numseen))

def print_mol2k_seen(mol2k, fields, outfile):
    sys.stderr.write("Writing {}\n".format(outfile))
    num = 0
    with open(outfile, 'w') as ofh:
        ofh.write("#"+"\t".join(fields)+"\n")
        for poskey in sorted(mol2k.keys()):
            for refaltkey in sorted(mol2k[poskey].keys()):
                d = mol2k[poskey][refaltkey]
                row = [ str(d[f]) if f in d else '.' for f in fields ]
                ofh.write("\t".join(row)+"\n")
                num += 1
    sys.stderr.write("  Wrote {} data rows\n".format(num))

if __name__ == '__main__':
    descr = """Find Mol2k variants in spreadsheet."""
    parser = ArgumentParser(description=descr)
#    parser.add_argument("mol2k", help="Mol2k file of variants.")
    parser.add_argument("spreadsheets", nargs="+",
                        help="Spreadsheet to check")
    parser.add_argument("-o", "--outdir", help="Directory for output.")
    parser.add_argument("-s", "--seen", default=False, action='store_true',
                        help="Reset Seen counter.")
    parser.add_argument("-f", "--force", default=False, action='store_true',
                        help="Overwrite existing files.")

    if len(sys.argv)<2:
        parser.print_help()
        sys.exit()
    args = parser.parse_args()

    mol2kfile = cftr.RESOURCE['cftr_db']
    if not os.path.isfile(mol2kfile):
        sys.exit("Could not find mol2k resource file: {}\n".format(mol2kfile))
    (mol2k, mol2k_fields) = parse_mol2k(mol2kfile, args.seen)
    seenfile = os.path.basename(mol2kfile).rstrip('txt').\
                       rstrip('csv').rstrip('.') + ".seen.txt"
    if args.outdir:
        seenfile = os.path.join(args.outdir, seenfile)
    for spreadsheet in args.spreadsheets:
        outfile = spreadsheet.rstrip('txt').rstrip('.') + ".mol2k.txt"
        if args.outdir: 
            outfile = os.path.join(args.outdir, os.path.basename(outfile))
        if have_file(outfile, args.force):
            sys.stderr.write("  Already have {}.\n".format(outfile))
        else:
            find_mol2k_variants(spreadsheet, mol2k, outfile)
    if have_file(seenfile, args.force):
        sys.stderr.write("  Already have {}.\n".format(seenfile))
    else:
        print_mol2k_seen(mol2k, mol2k_fields, seenfile)

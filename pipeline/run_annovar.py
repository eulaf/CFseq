#!/usr/bin/env python

"""
Run annovar gene-based annotation on variants in vcf file.
"""

import sys
import os
import re
import subprocess
from collections import defaultdict
from argparse import ArgumentParser

base_dir = os.path.dirname(__file__) or '.'
ref_dir = os.path.abspath(os.path.join(base_dir, '..', 'resources'))
lib_dir = os.path.abspath(os.path.join(base_dir, 'lib_cftr'))
sys.path.insert(0, lib_dir)

from common.misc import open_file, have_file
import cftr

annovarInputExe = cftr.EXE['annovar_convert']
annovarTableExe = cftr.EXE['annovar_table']
annovarExe = cftr.EXE['annovar_variation']
annovarRef = cftr.RESOURCE['annovar_ref']
annovarDBs = {
    'refGene':          'g',
    'clinvar_20140929': 'f',
    'cosmic70':         'f',
    'esp6500siv2_all':  'f',
    'snp138':           'f',
    'ljb26_all':        'f',
}

def annovar_protocol(dbdict, refdir):
    available_dbs = sorted([ db for db in dbdict.keys() if \
        os.path.isfile(os.path.join(refdir, 'hg19_'+db+'.txt')) ])
    glist = [ db for db in available_dbs if dbdict[db]=='g' ]
    flist = [ db for db in available_dbs if dbdict[db]=='f' ]
    put_at_end = [ db for db in flist if db.startswith('ljb') ]
    for db in put_at_end:
        flist.remove(db)
        flist.append(db)
    protocol_list = ','.join(glist + flist)
    operation_list = ','.join([dbdict[k] for k in glist + flist])
    return (protocol_list, operation_list)

def create_annovar_input_file(vcffiles, outlabel, args):
    outfile = outlabel + '.avinput';
    cmd = [annovarInputExe, '-format', 'vcf4', '-allsample', '-withfreq', ]
    cmd.extend([ '-includeinfo', ])  
    sys.stderr.write("\nCreating annovar input file: "+" ".join(cmd)+"\n")
    if have_file(outfile, args.force):
        sys.stderr.write("  Already have {}.\n".format(outfile))
    else: 
        lines = []
        for vcffile in vcffiles:
            sys.stderr.write("    Running: {} {}\n".format(" ".join(cmd),
                             vcffile))
            output = subprocess.check_output(cmd + [vcffile,])
            if not output:
                sys.stderr.write("  No output.\n")
            else:
                for line in output.split("\n"):
                    v = line.split("\t")[0:5]
                    if len(v)==5: 
                        (chrom, v[1]) = cftr.CFTR_to_hg19(v[0], v[1])
                        (v[0], v[2]) = cftr.CFTR_to_hg19(v[0], v[2])
                        lines.append("\t".join(map(str, v)))
        uniqlines = sorted(set(lines))
        with open(outfile, 'w') as ofh:
            ofh.write("\n".join(uniqlines))
        if not os.path.isfile(outfile):
            sys.stderr.write("  Failed to create {}\n".format(outfile))
            sys.exit(1)
    return outfile

def run_annovar(annov_input, outlabel, refdir, args):
    outfile1 = outlabel + '-hgvs.variant_function'
    outfile2 = outlabel + '-hgvs.exonic_variant_function'
    cmd = [annovarExe, '-build', 'hg19', '-hgvs', '-out', outlabel+'-hgvs', ]
    cmd.extend([ annov_input, refdir ])  
    sys.stderr.write("\nRunning annovar: "+" ".join(cmd)+"\n")
    if have_file(outfile1, args.force) and have_file(outfile2, args.force):
        sys.stderr.write("  Already have {} and {}.\n".format(outfile1,
                         outfile2))
    else: 
        subprocess.check_call(cmd)
        if not os.path.isfile(outfile1):
            sys.stderr.write("  Failed to create {}\n".format(outfile1))
            sys.exit(1)
        if not os.path.isfile(outfile2):
            sys.stderr.write("  Failed to create {}\n".format(outfile2))
            sys.exit(1)
    return (outfile1, outfile2)

def table_annovar(annov_input, outlabel, refdir, args):
    outfile = outlabel + '.hg19_multianno.txt'
    outfile1 = outlabel + '.variant_function'
    outfile2 = outlabel + '.exonic_variant_function'
    (protocol_list, operation_list) = annovar_protocol(annovarDBs, refdir)
    cmd = [ annovarTableExe, annov_input, refdir, '-out', outlabel, ]
    cmd.extend(['-buildver', 'hg19', '-out', outlabel, '-nastring', '.'])
    cmd.extend(['-protocol', protocol_list, '-operation', operation_list])
    sys.stderr.write("\nRunning table_annovar: "+" ".join(cmd)+"\n")
    if have_file(outfile, args.force):
        sys.stderr.write("  Already have {}.\n".format(outfile))
    else: 
        subprocess.check_call(cmd)
        if not os.path.isfile(outfile):
            sys.stderr.write("  Failed to create {}\n".format(outfile))
            sys.exit(1)
    return (outfile)


if __name__ == '__main__':
    descr = "Run annovar gene-based annotation on variants in vcf file."
    parser = ArgumentParser(description=descr)
    parser.add_argument("vcffiles", nargs="+", help="VCF file(s)")
    parser.add_argument("-l", "--label", default="annovar", 
                        help="Label for output files.")
    parser.add_argument("-o", "--outdir", help="Directory for output files.")
    parser.add_argument("-f", "--force", default=False, action="store_true",
                        help="Overwrite existing files.")

    if len(sys.argv)<2:
        parser.print_help()
        sys.exit()
    args = parser.parse_args()

    outlabel = os.path.join(args.outdir, args.label) if args.outdir \
               else args.label
    annov_input = create_annovar_input_file(args.vcffiles, outlabel, args)
#    table_annovar(annov_input, outlabel, annovarRef, args)
    run_annovar(annov_input, outlabel, annovarRef, args)

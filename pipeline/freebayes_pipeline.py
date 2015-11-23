#!/usr/bin/env python

"""
Run CFTR pipeline.
"""

import sys
import os
import datetime
import subprocess
from argparse import ArgumentParser

base_dir = os.path.dirname(__file__) or '.'
lib_dir = os.path.abspath(os.path.join(base_dir, 'lib_cftr'))
scripts_dir = os.path.abspath(base_dir)
sys.path.insert(0, lib_dir)

from common.misc import open_file, have_file, enumdate, timestamp
from common.profile import profile

SCRIPTS = {
    'run_freebayes': os.path.join(scripts_dir, 'run_freebayes.py'),
    'run_annovar': os.path.join(scripts_dir, 'run_annovar.py'),
    'create_spreadsheet': os.path.join(scripts_dir, 
                                       'vcf2spreadsheet_freebayes.py'),
}

def find_outfiles(outdir, infiles, ext, delim='.', num=1, debug=False,
                  quiet=False):
    files = [ f for f in os.listdir(outdir) if f.endswith(ext) ]
    outfiles = []
    for infile in infiles:
        basename = os.path.basename(infile).split(delim)
#        if debug: print basename
        matches = [ f for f in files if f.startswith(basename[0]) ]
        nummatches = len(matches)
        if nummatches==num:
            outfiles.extend(matches)
        elif nummatches>0 and not quiet:
            sys.stderr.write("**Bad number of matches for {}\n".format(
                             infile))
            sys.stderr.write("**Matches: {}\n".format(matches))
            sys.stderr.write("**Expected {} got {}\n".format(num, nummatches))
        elif not quiet:
            sys.stderr.write("No matching files for {}\n".format(infile))
    if not quiet:
        sys.stderr.write("Found {} {} file(s) in {}\n".format(len(outfiles),
                         ext, outdir))
    return [ os.path.join(outdir,f) for f in outfiles ]

@profile
def run_freebayes(bamfiles, fbdir, args):
    cmd = [SCRIPTS['run_freebayes'], '-o', fbdir,]
    cmd += ['-p', str(args.processes)]
    if args.logdir: cmd += ['--logdir', args.logdir]
    if args.force: cmd += ['-f',]
    cmd.extend(bamfiles)
    sys.stderr.write("\n==================================================\n")
    sys.stderr.write("Current time: {}\n".format(timestamp()))
    sys.stderr.write("Running freebayes: "+" ".join(cmd)+"\n")
    try:
        subprocess.check_call(cmd)
    except subprocess.CalledProcessError as e:
        sys.stderr.write("Error running freebayes pipeline.")
        raise
    outfiles = find_outfiles(fbdir, bamfiles, 'freebayes.filtered.vcf',
                             debug=args.debug)
    if not outfiles:
        sys.stderr.write("No freebayes .filtered.vcf files found.\n")
        sys.exit(1)
    return (outfiles)

@profile
def run_annovar(vcffiles, outdir, args):
    cmd = [SCRIPTS['run_annovar'], '-o', outdir]
    if args.label: cmd += ['-l', args.label + '.annovar']
    if args.force: cmd += ['-f',]
    cmd.extend(vcffiles)
    sys.stderr.write("\n==================================================\n")
    sys.stderr.write("Current time: {}\n".format(timestamp()))
    sys.stderr.write("Running annovar: "+" ".join(cmd)+"\n")
    try:
        subprocess.check_call(cmd)
    except subprocess.CalledProcessError as e:
        print e.output
        raise
    vffile = find_outfiles(outdir, [args.label,], 'annovar-hgvs.variant_function',
                           debug=args.debug)
    evffile = find_outfiles(outdir, [args.label,], 
                           'annovar-hgvs.exonic_variant_function',
                           debug=args.debug)
    if not vffile or not evffile:
        sys.exit(1)
    return (vffile[0], evffile[0])

@profile
def create_spreadsheet(label, samplevcfs, annovar_out, outdir):
    cmd = [SCRIPTS['create_spreadsheet'], '-l', label,'-o', outdir]
    for afile in annovar_out:
        cmd += ['-a', afile]
    if args.force: cmd += ['-f',]
    cmd.extend(samplevcfs)
    sys.stderr.write("\n==================================================\n")
    sys.stderr.write("Current time: {}\n".format(timestamp()))
    sys.stderr.write("Creating spreadsheet: "+" ".join(cmd)+"\n")
    try:
        subprocess.check_call(cmd)
    except subprocess.CalledProcessError as e:
        print e.output
        raise
    outfile = find_outfiles(outdir, [label,], 'results.txt', debug=args.debug)
    if not outfile:
        sys.exit(1)
    return outfile[0]


if __name__ == '__main__':
    descr = "Run GATK pipeline, including calculating actual depths and" +\
            " creating a spreadsheet of results annotated using annovar."
    parser = ArgumentParser(description=descr)
    parser.add_argument("bamfiles", nargs="+", help="BAM file(s).")
    parser.add_argument("-l", "--label", default=enumdate(d_only=True),
                        help="Label for output files.")
    parser.add_argument("-o", "--outdir", help="Directory for output files.")
#    parser.add_argument("--logfile", default=None,
#                        help="Save stderr to this log file (full path).")
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
    sys.stderr.write("{}\n".format(datetime.datetime.now()))
    vcffiles = run_freebayes(args.bamfiles, args.outdir, args)
    annovar_out = run_annovar(vcffiles, args.outdir, args)
    create_spreadsheet(args.label, vcffiles, annovar_out, args.outdir)
    sys.stderr.write("{}\n".format(datetime.datetime.now()))


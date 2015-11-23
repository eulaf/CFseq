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
    'run_gatk': os.path.join(scripts_dir, 'run_gatk.py'),
    'calculate_depths': os.path.join(scripts_dir, 'calculate_depths.py'),
    'separate_vcf': os.path.join(scripts_dir, 'separate_vcf.py'),
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
def run_gatk_all_bam(label, bamfiles, gatkdir, args):
    cmd = [SCRIPTS['run_gatk'], '-c', label, '-o', gatkdir,]
    cmd += ['-p', str(args.processes)]
    cmd += ['--logdir', args.logdir]
    if args.force: cmd += ['-f',]
    if args.debug: cmd += ['--debug',]
    if args.maxreads: cmd += ['-m', str(args.maxreads)]
    cmd.extend(bamfiles)
    sys.stderr.write("\n==================================================\n")
    sys.stderr.write("Current time: {}\n".format(timestamp()))
    sys.stderr.write("Running gatk: "+" ".join(cmd)+"\n")
    try:
        subprocess.check_call(cmd)
    except subprocess.CalledProcessError as e:
        print e.output
        raise
    outfiles = find_outfiles(gatkdir, bamfiles, 'gatk-cohort.vcf',
                             debug=args.debug)
    gvcf = find_outfiles(gatkdir, [label,], 'gatk-merged.vcf',
                         debug=args.debug)
    if not outfiles:
        sys.stderr.write("No gatk-cohort.vcf files found.\n")
        sys.exit(1)
    return (outfiles, gvcf[0])

@profile
def separate_vcf(vcf, gatkdir, bamfiles, args):
    cmd = [SCRIPTS['separate_vcf'], '-o', gatkdir, ]
    cmd += ['--logdir', args.logdir]
    if args.force: cmd += ['-f',]
    cmd.append(vcf)
    sys.stderr.write("\n==================================================\n")
    sys.stderr.write("Current time: {}\n".format(timestamp()))
    sys.stderr.write("Separating gVCF: "+" ".join(cmd)+"\n")
    try:
        subprocess.check_call(cmd)
    except subprocess.CalledProcessError as e:
        print e.output
        raise
    samplevcfs = find_outfiles(gatkdir, bamfiles, 'separated.vcf', delim='-',
                           debug=args.debug)
    if not samplevcfs:
        sys.exit(1)
    return samplevcfs

@profile
def calculate_depths(vcffiles, bamfiles, gatkdir, args):
    cmd = [SCRIPTS['calculate_depths'], '-o', gatkdir, ]
    if args.logdir: cmd += ['--logdir', args.logdir]
    cmd += ['-p', str(args.processes)]
    if args.force: cmd += ['-f',]
    cmd.extend(vcffiles)
    cmd.extend(bamfiles)
    sys.stderr.write("\n==================================================\n")
    sys.stderr.write("Current time: {}\n".format(timestamp()))
    sys.stderr.write("Calculating depths: "+" ".join(cmd)+"\n")
    try:
        subprocess.check_call(cmd)
    except subprocess.CalledProcessError as e:
        print e.output
        raise
    outfiles = find_outfiles(gatkdir, bamfiles, 'depths.txt', 
                             debug=args.debug)
    if not outfiles:
        sys.stderr.write("No depths.txt files found.\n")
        sys.exit(1)
    return outfiles


if __name__ == '__main__':
    descr = "Run GATK pipeline, including calculating actual depths."
    parser = ArgumentParser(description=descr)
    parser.add_argument("bamfiles", nargs="+", help="BAM file(s).")
    parser.add_argument("-l", "--label", default=enumdate(d_only=True),
                        help="Label for output files.")
    parser.add_argument("-o", "--outdir", help="Directory for gatk output.")
    parser.add_argument("-m", "--maxreads", type=int, default=15000,
                      help="Set gatk maxReadsInRegionPerSample to this value.")
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
    (vcffiles, gvcf) = run_gatk_all_bam(args.label, args.bamfiles, 
                                        args.outdir, args)
    samplevcfs = separate_vcf(gvcf, args.outdir, args.bamfiles, args)
    depthfiles = calculate_depths(samplevcfs, args.bamfiles, args.outdir, args)
    sys.stderr.write("{}\n".format(datetime.datetime.now()))


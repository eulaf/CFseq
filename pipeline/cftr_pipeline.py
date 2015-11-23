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

from common.profile import profile
from common.misc import open_file, have_file, check_directory, enumdate, timestamp

SCRIPTS = {
    'primer_trim': os.path.join(scripts_dir, 'primer_trim.py'),
    'run_bwa': os.path.join(scripts_dir, 'run_bwa_mem.py'),
    'preprocess_bams': os.path.join(scripts_dir, 'preprocess_bams.py'),
    'run_freebayes': os.path.join(scripts_dir, 'run_freebayes.py'),
    'gatk_pipeline': os.path.join(scripts_dir, 'gatk_pipeline.py'),
    'run_annovar': os.path.join(scripts_dir, 'run_annovar.py'),
    'fb_create_spreadsheet': os.path.join(scripts_dir, 
                                       'vcf2spreadsheet_freebayes.py'),
    'gatk_create_spreadsheet': os.path.join(scripts_dir, 
                                       'vcf2spreadsheet_gatk.py'),
    'add_mol2k': os.path.join(scripts_dir, 'add_mol2k.py'),
}

def find_outfiles(outdir, infiles, ext, delim='.', num=1, quiet=False):
    files = [ f for f in os.listdir(outdir) if f.endswith(ext) ]
    outfiles = []
    for infile in infiles:
        basename = os.path.basename(infile).split(delim)
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
def trim_all_fastq(fqfiles, aligndir, logdir, args):
    cmd = [SCRIPTS['primer_trim'], '-o', aligndir,]
    cmd += ['-p', str(args.processes)]
    if logdir: cmd += ['--logdir', logdir]
    if args.force: cmd += ['-f',]
    if args.label: cmd += ['-s', args.label + '.summary.txt' ]
    cmd.extend(fqfiles)
    sys.stderr.write("\n==================================================\n")
    sys.stderr.write("Current time: {}\n".format(timestamp()))
    sys.stderr.write("Trimming fastq files: "+" ".join(cmd)+"\n")
    outfiles = find_outfiles(aligndir, fqfiles, 'trimmed.fastq', quiet=True)
    if len(outfiles)==len(fqfiles):
        sys.stderr.write("Already have trimmed fq files.\n")
        return outfiles
    try:
        subprocess.check_call(cmd)
    except Exception, e:
        raise
    outfiles = find_outfiles(aligndir, fqfiles, 'trimmed.fastq')
    if not outfiles:
        sys.stderr.write("No trimmed.fastq files found.\n")
        sys.exit(1)
    return outfiles

@profile
def run_bwa_all_fastq(fqfiles, aligndir, logdir, args):
    cmd = [SCRIPTS['run_bwa'], '-o', aligndir, '-p', str(args.processes)]
    if logdir: cmd += ['--logdir', logdir]
    if args.force: cmd += ['-f',]
    if args.sam: cmd += ['-s',]
    cmd.extend(fqfiles)
    sys.stderr.write("\n==================================================\n")
    sys.stderr.write("Current time: {}\n".format(timestamp()))
    sys.stderr.write("Running bwa: "+" ".join(cmd)+"\n")
    try:
        subprocess.check_call(cmd)
    except subprocess.CalledProcessError as e:
        sys.stderr.write("Error running bwa.\n")
        raise
    delim = '_L001' if '_L001' in fqfiles[0] else '_R'
    ext = 'genomic_refseq.bam'
    outfiles = find_outfiles(aligndir, [ f for f in fqfiles if 'R1' in f ], 
                             ext, delim=delim)
    if not outfiles:
        sys.stderr.write("No bam files found.\n")
        sys.exit(1)
    return outfiles

@profile
def gatk_pipeline(label, bamfiles, gatkdir, logdir, args):
    cmd = [SCRIPTS['gatk_pipeline'], '-l', label, '-o', gatkdir,]
    cmd += ['-p', str(args.processes)]
    if logdir: cmd += ['--logdir', logdir]
    if args.maxreads: cmd += ['-m', str(args.maxreads)]
    if args.force: cmd += ['-f',]
    cmd.extend(bamfiles)
    sys.stderr.write("\n==================================================\n")
    sys.stderr.write("Current time: {}\n".format(timestamp()))
    sys.stderr.write("Running gatk pipeline: "+" ".join(cmd)+"\n")
    try:
        subprocess.check_call(cmd)
    except Exception, e:
        raise
    gvcf = find_outfiles(gatkdir, [label,], 'gatk-merged.vcf')
    vcffiles = find_outfiles(gatkdir, bamfiles, 'separated.vcf', delim='-')
    if not vcffiles:
        sys.exit(1)
    return (vcffiles, gvcf[0])

@profile
def run_freebayes(bamfiles, fbdir, logdir, args):
    cmd = [SCRIPTS['run_freebayes'], '-o', fbdir,]
    cmd += ['-p', str(args.processes)]
    if args.force: cmd += ['-f',]
    cmd.extend(bamfiles)
    sys.stderr.write("\n==================================================\n")
    sys.stderr.write("Current time: {}\n".format(timestamp()))
    sys.stderr.write("Running freebayes: "+" ".join(cmd)+"\n")
    logfile = os.path.join(logdir, args.label + '.freebayes.log')
    try:
        with open(logfile, 'w') as logfh:
            subprocess.check_call(cmd, stderr=logfh)
    except subprocess.CalledProcessError as e:
        sys.stderr.write("Error running freebayes."+\
                         " Please check logfile {}\n".format(logfile))
        raise
    vcffiles = find_outfiles(fbdir, bamfiles, 'freebayes.filtered.vcf')
    if not vcffiles:
        sys.stderr.write("No freebayes .filtered.vcf files found.\n")
        sys.exit(1)
    return (vcffiles)

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
    except Exception, e:
        raise
    vffile = find_outfiles(outdir, [args.label,], 'annovar-hgvs.variant_function')
    evffile = find_outfiles(outdir, [args.label,], 
                           'annovar-hgvs.exonic_variant_function')
    if not vffile or not evffile:
        sys.exit(1)
    return (vffile[0], evffile[0])

@profile
def fb_create_spreadsheet(label, samplevcfs, annovar_out, outdir):
    label += '.freebayes'
    cmd = [SCRIPTS['fb_create_spreadsheet'], '-l', label,'-o', outdir]
    for afile in annovar_out:
        cmd += ['-a', afile]
    if args.force: cmd += ['-f',]
    cmd.extend(samplevcfs)
    sys.stderr.write("\n==================================================\n")
    sys.stderr.write("Current time: {}\n".format(timestamp()))
    sys.stderr.write("Creating FreeBayes spreadsheet: "+" ".join(cmd)+"\n")
    try:
        subprocess.check_call(cmd)
    except subprocess.CalledProcessError as e:
        print e.output
        raise
    outfile = find_outfiles(outdir, [label,], 'freebayes.results.txt')
    if not outfile:
        sys.exit(1)
    return outfile[0]

@profile
def gatk_create_spreadsheet(label, samplevcfs, annovar_out, gatkdir, outdir):
    label += '.gatk'
    cmd = [SCRIPTS['gatk_create_spreadsheet'], '-l', label,'-o', outdir]
    cmd += ['-d', gatkdir]
    for afile in annovar_out:
        cmd += ['-a', afile]
    if args.force: cmd += ['-f',]
    cmd.extend(samplevcfs)
    sys.stderr.write("\n==================================================\n")
    sys.stderr.write("Current time: {}\n".format(timestamp()))
    sys.stderr.write("Creating GATK spreadsheet: "+" ".join(cmd)+"\n")
    try:
        subprocess.check_call(cmd)
    except subprocess.CalledProcessError as e:
        print e.output
        raise
    outfile = find_outfiles(outdir, [label,], 'gatk.results.txt')
    if not outfile:
        sys.exit(1)
    return outfile[0]

@profile
def annotate_spreadsheets(spreadsheets, outdir, args):
    cmd = [SCRIPTS['add_mol2k'], '-o', outdir]
    if args.force: cmd += ['-f',]
    cmd.extend(spreadsheets)
    sys.stderr.write("\n==================================================\n")
    sys.stderr.write("Current time: {}\n".format(timestamp()))
    sys.stderr.write("Annotating spreadsheets: "+" ".join(cmd)+"\n")
    try:
        subprocess.check_call(cmd)
    except subprocess.CalledProcessError as e:
        print e.output
        raise
    outfiles = find_outfiles(outdir, spreadsheets, 'mol2k.txt', delim='.results')
    if not outfiles:
        sys.exit(1)
    return outfiles


if __name__ == '__main__':
    descr = "Run CFTR pipeline."
    parser = ArgumentParser(description=descr)
    parser.add_argument("fqfiles", nargs="+", 
                        help="Fastq file(s). May be gzipped.")
    parser.add_argument("-l", "--label", default=enumdate(d_only=True),
                        help="Label for output files.")
    parser.add_argument("-m", "--maxreads", type=int, default=15000,
                        help="Set gatk maxReadsInRegionPerSample to this "+\
                        " value (default: 15000).")
    parser.add_argument("-o", "--outdir", help="Directory for output.")
    parser.add_argument("-p", "--processes", default=1, type=int,
                        help="Number of processes to run.")
    parser.add_argument("-r", "--resultsdir", default="results",
                        help="Output subdirectory for result files.")
    parser.add_argument("-s", "--sam", default=False, 
                        action="store_true", help="Keep sam files.")
    parser.add_argument("--gatkdir", default="gatk_output",
                        help="Output subdirectory for GATK files "+\
                        " (default: gatk_output).")
    parser.add_argument("--fbdir", default="freebayes_output",
                        help="Output subdirectory for freebayes files "+\
                        " (default: freebayes_output).")
    parser.add_argument("--logdir", default="logs",
                        help="Output subdirectory for log files.")
    parser.add_argument("--aligndir", default="alignfiles",
                        help="Output subdirectory for alignment files.")
    parser.add_argument("-f", "--force", default=False, action="store_true",
                        help="Overwrite existing files.")

    if len(sys.argv)<2:
        parser.print_help()
        sys.exit()
    args = parser.parse_args()

    sys.stderr.write("{}\n".format(datetime.datetime.now()))
    aligndir = check_directory(args.aligndir, args.outdir)
    logdir = check_directory(args.logdir, args.outdir) if args.logdir else None
    resultsdir = check_directory(args.resultsdir, args.outdir) if \
             args.resultsdir else None

    trimmedfq = trim_all_fastq(args.fqfiles, aligndir, logdir, args)
    bamfiles = run_bwa_all_fastq(trimmedfq, aligndir, logdir, args)

    fbdir = check_directory(args.fbdir, args.outdir)
    fb_vcfs = run_freebayes(bamfiles, fbdir, logdir, args)
    fb_annovar = run_annovar(fb_vcfs, fbdir, args)
    fb_spreadsheet = fb_create_spreadsheet(args.label, fb_vcfs, fb_annovar, 
                                           resultsdir)

    gatkdir = check_directory(args.gatkdir, args.outdir)
    (gatk_vcfs, gvcf) = gatk_pipeline(args.label, bamfiles, gatkdir, logdir, 
                                      args)
    gatk_annovar = run_annovar([gvcf,], gatkdir, args)
    ga_spreadsheet = gatk_create_spreadsheet(args.label, gatk_vcfs, 
                                         gatk_annovar, gatkdir, resultsdir)
    
    spreadsheets = annotate_spreadsheets([fb_spreadsheet, ga_spreadsheet], 
                                         resultsdir, args)
    sys.stderr.write("{}\n".format(datetime.datetime.now()))


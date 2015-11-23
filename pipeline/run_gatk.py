#!/usr/bin/env python

"""
Run gatk variant filter.
"""

import sys
import os
import re
import multiprocessing as mp
from subprocess import check_call
from collections import defaultdict
from argparse import ArgumentParser

base_dir = os.path.dirname(__file__) or '.'
lib_dir = os.path.abspath(os.path.join(base_dir, 'lib_cftr'))
sys.path.insert(0, lib_dir)

from common.profile import profile
from common.misc import open_file, have_file, timestamp
from cftr import EXE, RESOURCE

picardJar = EXE['picard_jar']
picardExe = ["java", "-Xmx2g", "-jar", picardJar]
if not os.path.isfile(picardJar):
    sys.exit("Could not find {}\n".format(picardJar))
gatkJar = EXE['gatk_jar']
gatkExe = ["java", "-Xmx4g", "-jar", gatkJar]
if not os.path.isfile(gatkJar):
    sys.exit("Could not find {}\n".format(gatkJar))

@profile
def index_bam(bamfile, args, logfh=sys.stderr):
    label = bamfile.rstrip('bam').rstrip('.')
    outidx = label+'.bai'
    cmd = picardExe[:] + ['BuildBamIndex', 'I='+bamfile]
    if not have_file(outidx, args.force, quiet=True, stderr=logfh):
        logfh.write("\nIndex bam: "+" ".join(cmd)+"\n")
        logfh.flush()
        check_call(cmd, stderr=logfh)
    if not os.path.isfile(outidx):
        logfh.write("  Failed to create {}\n".format(outidx))
        sys.stderr.write("  Failed to create {}\n".format(outidx))
        sys.exit(1)

def get_outlabel(bamfile, outdir):
    label = os.path.basename(bamfile).rstrip('bam').rstrip('.')
    if outdir:
        if os.path.isdir(outdir):
            label = os.path.join(outdir, label)
        else:
            msg = "ERROR: Outdir {} does not exist.\n".format(outdir)
            sys.stderr.write(msg)
            sys.exit(1)
    return label

@profile
def run_gatk(bamfile, label, ref, args, logfh=sys.stderr):
    caller = 'HaplotypeCaller'
    outvcf = label + '.gatk.vcf'
    if args.cohort:
        outvcf = outvcf.replace('.vcf', '-cohort.vcf')
    cmd = gatkExe[:] + ['-T', caller, '--genotyping_mode', 'DISCOVERY']
    cmd.extend(['-R', ref, '-I', bamfile, '-o', outvcf])  
    if args.intervals: cmd.extend(['-L', args.intervals])
    if args.dbsnp: cmd.extend(["--dbsnp", args.dbsnp])
    if args.maxreads:
        cmd.extend(['--maxReadsInRegionPerSample', str(args.maxreads)])
    if args.debug:
        bamout = label + '.gatk-debug.bam'
        cmd.extend(['-bamout', bamout])
    if args.cohort:
        cmd.extend(['-ERC', 'GVCF', '--variant_index_type', 'LINEAR'])
        cmd.extend(['--variant_index_parameter', '128000'])
#    if args.basequal:
#        cmd.extend(['-mbq', str(args.basequal)])
    logfh.write("\nGATK: "+" ".join(cmd)+"\n")
    if have_file(outvcf, args.force, stderr=logfh):
        logfh.write("  Already have {}.\n".format(outvcf))
    else: 
        logfh.flush()
        check_call(cmd, stderr=logfh)
    if not os.path.isfile(outvcf):
        logfh.write("  Failed to create {}\n".format(outvcf))
        sys.stderr.write("  Failed to create {}\n".format(outvcf))
        sys.exit(1)
    return outvcf

def call_variants_gatk(bamfile, ref, args):
    outlabel = get_outlabel(bamfile, args.outdir)
    loglabel = get_outlabel(bamfile, args.logdir if args.logdir else
               args.outdir)
    logfile = loglabel + ".gatk.log"
    logfh = open(logfile, 'w')
    sys.stderr.write("  Running GATK on {}\n".format(bamfile))
    sys.stderr.write("    Log file {}\n".format(logfile))
    try:
        logfh.write("Start time: {}\n".format(timestamp()))
        index_bam(bamfile, args, logfh)
        vcffile = run_gatk(bamfile, outlabel, ref, args, logfh)
        if os.path.isfile(vcffile):
            sys.stderr.write("    VCF file: {}\n".format(vcffile))
        logfh.write("Finish time: {}\n".format(timestamp()))
    except Exception, e:
        e.args += (bamfile, )
        sys.stderr.write("Error running gatk." +\
                         "  Check log file {}\n".format(logfile))
        raise
    finally:
        logfh.close()
    return vcffile

def cohort_merge_gvcfs(vcfs, ref, args):
    sys.stderr.write("Genotyping gVCFs: {} vcfs\n".format(len(vcfs)))
    label = args.cohort
    outlabel = os.path.join(args.outdir, label) if args.outdir else label
    outvcf = outlabel + '.gatk-merged.vcf'
    cmd = gatkExe[:] + ['-T', 'GenotypeGVCFs']
    cmd.extend(['-R', ref, '-o', outvcf])  
    if args.intervals: cmd.extend(['-L', args.intervals])
    if args.dbsnp: cmd.extend(["--dbsnp", args.dbsnp])
    variants = ("--variant "+" --variant ".join(vcfs)).split(' ')
    cmd += variants
    if have_file(outvcf, args.force):
        sys.stderr.write("  Already have {}.\n".format(outvcf))
    else: 
        loglabel = os.path.join(args.logdir, label) if args.logdir else label
        logfile = loglabel + '.gatk-merged.log'
        sys.stderr.write("    Log file {}\n".format(logfile))
        with open(logfile, 'w') as logfh:
            logfh.write("CMD: {}\n".format(cmd))
            logfh.flush()
            check_call(cmd, stderr=logfh)
    if os.path.isfile(outvcf):
        sys.stderr.write("  Merged gVCF: {}\n".format(outvcf))
    else:
        sys.stderr.write("  Failed to create {}\n".format(outvcf))
        sys.exit(1)
    return outvcf

if __name__ == '__main__':
    descr = "Run gatk on bam file(s)."
    parser = ArgumentParser(description=descr)
    parser.add_argument("bamfiles", nargs="+", help="BAM file(s)")
    parser.add_argument("-c", "--cohort", 
                        help="Run optimized workflow for cohorts."+\
                        "  Final merged vcf is named using given label.")
    parser.add_argument("-o", "--outdir", help="Directory for output files.")
    parser.add_argument("-m", "--maxreads", type=int, 
                        help="Set maxReadsInRegionPerSample to this value.")
#    parser.add_argument("-b", "--basequal", type=int,
#                        help="Set min_base_quality_score to this value.")
    parser.add_argument("-f", "--force", default=False, action="store_true",
                        help="Overwrite existing files.")
    parser.add_argument("-p", "--processes", default=1, type=int,
                        help="Number of processes to run.")
    parser.add_argument("--debug", default=False, action="store_true",
                        help="Create gatk debug bam files.")
    parser.add_argument("--logdir", help="Directory for log files.")
    parser.add_argument("--ref",help="Use this reference instead of default.")

    if len(sys.argv)<2:
        parser.print_help()
        sys.exit()
    args = parser.parse_args()

    if args.ref:
        ref = os.path.abspath(args.ref)
    else:
        ref = RESOURCE['ref_fa']
        args.intervals = RESOURCE['covered_roi']
        args.dbsnp = RESOURCE['dbsnp_vcf']
        if not os.path.isfile(ref):
            sys.exit("Ref file {} not found.\n".format(ref))
        if not os.path.isfile(args.intervals):
            sys.exit("Interval file {} not found.\n".format(args.intervals))
        if not os.path.isfile(args.dbsnp):
            sys.exit("dbSNP file {} not found.\n".format(args.dbsnp))
    sys.stderr.write("\nRunning GATK on {} bam files using {}.\n".format(
                   len(args.bamfiles), "{} process{}".format(args.processes,
                   'es' if args.processes>1 else '')))
    pool = mp.Pool(processes=args.processes)
    results = [ pool.apply_async(call_variants_gatk, args=[bamfile,ref,args])
                for bamfile in args.bamfiles ]
    vcfs = [p.get() for p in results ]
    if args.cohort:
        cohort_merge_gvcfs(vcfs, ref, args)


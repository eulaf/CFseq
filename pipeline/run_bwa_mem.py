#!/usr/bin/env python

"""
Run bwa on fastq files.
"""

import sys
import os
import subprocess
import multiprocessing as mp
from collections import defaultdict
from argparse import ArgumentParser

base_dir = os.path.dirname(__file__) or '.'
lib_dir = os.path.abspath(os.path.join(base_dir, 'lib_cftr'))
sys.path.insert(0, lib_dir)

from common.misc import open_file, have_file, timestamp, remove_file
from common.profile import profile
from cftr import EXE, RESOURCE

bwaExe = EXE['bwa']
if not os.path.isfile(bwaExe):
    sys.exit("Could not find {}\n".format(bwaExe))
picardJar = EXE['picard_jar']
picardExe = ["java", "-Xmx2g", "-jar", picardJar]
if not os.path.isfile(picardJar):
    sys.exit("Could not find {}\n".format(picardJar))
samtoolsExe = EXE['samtools']
if not os.path.isfile(samtoolsExe):
    sys.exit("Could not find {}\n".format(samtoolsExe))


def get_samplename_from_fqname(fqfile):
    label = os.path.basename(fqfile)
    if '_L001_' in label:
        sample = label.split('_L001_')[0]
    elif '_R1_' in label:
        sample = label.split('_R1_')[0]
    elif '_R2_' in label:
        sample = label.split('_R2_')[0]
    else:
        sys.stderr.write("Bad fqfile name: {}\n".format(fqfile))
        sys.exit(1)
    return sample

def pair_fq_files(fqfiles):
    """Pair fastq files assuming Illumina naming convention:
       sample_L001_R[12]_001.fastq"""
    fqpairs = defaultdict(list)
    for fqfile in sorted(fqfiles):
        samplename = get_samplename_from_fqname(fqfile)
        fqpairs[samplename].append(fqfile)
    return fqpairs

def get_bwa_version():
    cmd = [bwaExe, ]
    sp = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = sp.communicate()
    version = [ l for l in err.split("\n") if l.startswith('Version') ]
    return "BWA "+ version[0]

def run_bwa(sample, outlabel, ref, fqfiles, logfh, force):
    samfile = outlabel + ".sam"
    bamfile = outlabel + ".bam"
    logfh.write("Output sam: {}\n".format(samfile))
    readgroup = "\\t".join(['@RG', "ID:"+sample, "SM:"+sample, "PL:Illumina",
                           "LB:"+sample, "PU:unit1"]);
    cmd = [bwaExe, 'mem', '-M', '-R', readgroup, ref,] + fqfiles
    logfh.write(" ".join(cmd)+"\n")
    if have_file(samfile, force, stderr=logfh):
        logfh.write("Already have sam file: {}\n".format(samfile))
    elif have_file(bamfile, force, stderr=logfh):
        logfh.write("Already have bam file: {}\n".format(bamfile))
    else:
        logfh.flush()
        output = subprocess.check_output(cmd, stderr=logfh)
        with open(samfile, 'w') as ofh:
            ofh.write(output)
    return samfile

def create_sorted_bam(samfile, outlabel, logfh, force, keepsam):
    bamfile = outlabel + '.bam'
    logfh.write("\nSorting sam: {}\n".format(samfile))
    logfh.write("Creating bam: {}\n".format(bamfile))
    cmd = picardExe[:] + [ 'SortSam', 'I='+samfile, 'O='+bamfile,
           'SORT_ORDER=coordinate']
    logfh.write(" ".join(cmd)+"\n")
    if have_file(bamfile, force, stderr=logfh):
        logfh.write("Already have bam file: {}\n".format(bamfile))
    else:
        logfh.flush()
        subprocess.call(cmd, stderr=logfh)
    if not keepsam and have_file(bamfile, stderr=logfh):
        remove_file(samfile, stderr=logfh)
    return bamfile

def align_create_bam(sample, fqfiles, ref, reflabel, args):
    outlabel = "{}-{}".format(sample, reflabel)
    loglabel = os.path.join(args.logdir, outlabel) if args.logdir \
               else outlabel
    logfile = loglabel + ".bwa.log"
    sys.stderr.write("  Mapping {}\n".format(sample))
    sys.stderr.write("    Log file {}\n".format(logfile))
    if args.outdir: 
        outlabel = os.path.abspath(os.path.join(args.outdir, outlabel))
    logfh = open(logfile, 'w')
    try:
        logfh.write("Start time {}\n".format(timestamp()))
        samfile = run_bwa(sample, outlabel, ref, fqfiles, logfh, args.force)
        bamfile = create_sorted_bam(samfile, outlabel, logfh, args.force, 
                                    args.sam)
        logfh.write("Finish time {}\n".format(timestamp()))
    except Exception, e:
        e.args += (sample, )
        sys.stderr.write("Error running bwa." +\
                         "  Check log file {}\n".format(logfile))
        raise
    finally:
        logfh.close()
    if os.path.isfile(bamfile):
        sys.stderr.write("    Bam file {}\n".format(bamfile))
    return (sample, bamfile)

@profile
def create_bams(fqpairs, ref, reflabel, args):
    sys.stderr.write("\nMapping {} samples using {}.\n".format(len(fqpairs),
           "{} process{}".format(args.processes, 'es' if args.processes>1 
           else '')))
    pool = mp.Pool(processes=args.processes)
    results = [ pool.apply_async(align_create_bam, args=[sample,
                    fqpairs[sample], ref, reflabel, args]) for sample in 
                    sorted(fqpairs.keys()) ]
    bamfiles = dict([p.get() for p in results ])
    return bamfiles

if __name__ == '__main__':
    descr = "Run bwa on fastq files and create sorted bam files."
    descr += " Bam files are sorted by coordinates using Picard."
    parser = ArgumentParser(description=descr)
    parser.add_argument("fqfiles", nargs="+", 
                        help="Fastq file(s).  May be gzipped.")
    parser.add_argument("-s", "--sam", default=False, action="store_true",
                        help="Keep sam file.")
    parser.add_argument("-o", "--outdir", help="Directory for output files.")
    parser.add_argument("-f", "--force", default=False, action="store_true",
                        help="Overwrite existing files.")
    parser.add_argument("-p", "--processes", default=1, type=int,
                        help="Number of processes to run.")
    parser.add_argument("--logdir", help="Directory for log files.")
    parser.add_argument("--ref",help="Use this reference instead of default.")

    if len(sys.argv)<2:
        parser.print_help()
        sys.exit()
    args = parser.parse_args()

    ref = os.path.abspath(args.ref) if args.ref else RESOURCE['ref_fa']
    if not os.path.isfile(ref):
        sys.exit("Ref file {} not found.\n".format(ref))
    reflabel = os.path.basename(ref).rstrip('fast').rstrip('.')
    sys.stderr.write("Ref label: {}\n".format(reflabel))
    fqpairs = pair_fq_files(args.fqfiles)
    sys.stderr.write(get_bwa_version()+"\n")
    create_bams(fqpairs, ref, reflabel, args)


#!/usr/bin/env python

"""
Trim primers from fastq files.
Pairs fastq files according to Illumina naming convention:
samplename_L001_R[12]_001.fastq.

1. Identify primers in overlapping amplicon regions.  These are the primers 
we want to trim.
2. Determine length of longest primer (maxlength)
3. Create fasta of first maxlength+2 bases of each sequence in fastq
4. Align primer sequences against fasta file
5. For each fastq read, if matching primer is in overlapping region, then trim
primer.  If no matching primer found, trim maxlength+2 bases.  Otherwise, do
not trim.
"""

import sys
import os
import re
import multiprocessing as mp
import subprocess
from collections import defaultdict
from argparse import ArgumentParser

base_dir = os.path.dirname(__file__) or '.'
lib_dir = os.path.abspath(os.path.join(base_dir, 'lib_cftr'))
sys.path.insert(0, lib_dir)

from common.misc import open_file, have_file, have_files, timestamp
from common.math import mean, median
from common.profile import profile
from ParseFastA import FastAParser
from ParseFastQ import FastQParser
from cftr import EXE, RESOURCE

MAX_READS = 800000
CM_EXE = EXE['cm']
if not os.path.isfile(CM_EXE):
    sys.exit("Could not find {}\n".format(CM_EXE))

def primer_info(primerfa):
    """Reads a fasta file and returns (primer_info_dict, max_primer_length)
       where max_primer_length is the length of the longest primer and 
       primer_info_dict is keyed by primer name and values a dict containing 
       length of primer and location info: start, end, strand."""
    sys.stderr.write("Parsing primer file: {}\n".format(primerfa))
    inseq = FastAParser(primerfa)
    primerinfo = {}
    max_primer_len = 0
    for seqrec in inseq:
        primerlen = len(seqrec.seq)
        if max_primer_len < primerlen:
            max_primer_len = primerlen
        primerinfo[seqrec.id] = { 'len':primerlen }
        des = seqrec.description
        if des and ':' in des:
            (ref, loc, strand) = des.split(':')
            (start, end) = loc.split('-')
            primerinfo[seqrec.id]['start'] = int(start)
            primerinfo[seqrec.id]['end'] = int(end)
            primerinfo[seqrec.id]['strand'] = int(strand)
    numprimers = len(primerinfo.keys())
    sys.stderr.write("  Num primers: \t{} seqs\n".format(numprimers))
    sys.stderr.write("  Max primer length:\t{} bp\n".format(max_primer_len))
    return (primerinfo, max_primer_len)

def find_overlaps(pinfo, debug):
    amplicons = [ name.replace('_F','') for name in pinfo.keys() \
                  if name.endswith('_F') ]
    amplicons = sorted(amplicons, key=lambda x: 
            pinfo[x+'_F']['start'] if pinfo[x+'_F']['start'] < \
            pinfo[x+'_R']['start'] else pinfo[x+'_R']['start'])
    lastamp = []
    for ampname in amplicons:
        thisamp = []
        if pinfo[ampname+'_F']['start'] < pinfo[ampname+'_R']['start']:
            thisamp = [ampname+'_F', ampname+'_R']
        else:
            thisamp = [ampname+'_R', ampname+'_F']
        pinfo[thisamp[0]]['overlap'] = False
        pinfo[thisamp[1]]['overlap'] = False
        if lastamp:
            if pinfo[lastamp[1]]['end'] > pinfo[thisamp[0]]['start']:
                pinfo[thisamp[0]]['overlap'] = True
                pinfo[lastamp[1]]['overlap'] = True
        lastamp = thisamp
    if debug:
        for ampname in amplicons:
            ampF = ampname+'_F'
            ampR = ampname+'_R'
            sys.stderr.write("{}\t{:,d} - {:,d}\t{:,d} - {:,d}\t{} {}\n".\
                format(ampname, pinfo[ampF]['start'], pinfo[ampF]['end'], 
                pinfo[ampR]['start'], pinfo[ampR]['end'], 
                'overlap' if pinfo[ampF]['overlap'] else '-', 
                'overlap' if pinfo[ampR]['overlap'] else '-'))

def fastq_file_label(fqfile, outdir):
    if outdir and not os.path.isdir(outdir):
        sys.stderr.write("ERROR: Outdir {} does not exist\n".format(outdir))
        sys.exit(1)
    label = os.path.basename(fqfile)
    while '.fastq' in label or '.fq' in label:
        label = os.path.splitext(label)[0]
    if outdir:
        label = os.path.join(outdir, label)
    return label

def create_fasta_of_primer_region(fqfile, max_trim_len, outlabel, logfh, 
                                  force):
    """Creates fasta files of the first max_trim_len bases of
    each sequence in given fqfile.  Also, returns a list with names of
    all sequences in fqfile."""
    outfile = outlabel + ".primer_region{}.fa".format(max_trim_len)
    logfh.write("    Creating 5' fasta: {}\n".format(outfile))
    outfiles = [outfile,]
    if have_file(outfile, force, stderr=logfh):
        logfh.write("      Already have {}\n".format(outfile))
    else:
        logfh.write("      Writing {}\n".format(outfile))
        numseqs = 0
        totseqs = 0
        outfa = open(outfile, 'w')
        inseq = FastQParser(fqfile)
        for seqrec in inseq:
            numseqs += 1
            subseq = seqrec.seq[0:max_trim_len]
            outfa.write(">{}\n{}\n".format(seqrec.id, subseq))
            if numseqs==MAX_READS:
                outfa.close()
                logfh.write("      Wrote {} seqs\n".format(numseqs))
                numfiles = len(outfiles) + 1
                outfile = outlabel + ".primer_region{}-{}.fa".format(
                          max_trim_len, numfiles)
                logfh.write("      Writing {}\n".format(outfile))
                outfiles.append(outfile)
                outfa = open(outfile, 'w')
                totseqs += numseqs
                numseqs = 0
        logfh.write("      Wrote {} seqs\n".format(numseqs))
        if len(outfiles)>1:
            totseqs += numseqs
            logfh.write("  Wrote {} total seqs\n".format(totseqs))
    return outfiles

def run_aligner(queryfiles, primerfa, outlabel, logfh, force):
    outfile = outlabel + "-primers.cm.out"
    logfh.write("\nAlignment output in {}\n".format(outfile))
    if (have_file(outfile, force, stderr=logfh)):
        logfh.write("      Already have {}\n".format(outfile))
        return (outfile, None)
    try:
        with open(outfile, 'w') as ofh:
            for queryfile in queryfiles:
                logfh.write("    Aligning {}\n".format(queryfile))
                cmd = [CM_EXE, queryfile, primerfa, "-minscore", "12", 
                       "-minmatch", "8", "-tags", "-alignments"]
                logfh.write("      Running {}\n".format(" ".join(cmd)))
                logfh.flush()
                subprocess.check_call(cmd, stdout=ofh, stderr=logfh)
    except subprocess.CalledProcessError as e:
        sys.stderr.write("Error running cross_match for {}\n".format(outfile))
        raise
    if not os.path.isfile(outfile):
        sys.stderr.write("Output file {} not found.\n".format(outfile))
    return outfile

def parse_alignout(align_output):
    results = {}
    with open(align_output, 'r') as fh:
        aligns = [ l for l in fh.readlines() if l.startswith('ALIGN') ]
        for align in aligns:
            vals = align.split()
            if len(vals)==13:  # only analyze matches not complemented
                fragname = vals[5]
                start = int(vals[10])
                end = int(vals[11])
                left = int(vals[12].lstrip('(').rstrip(')'))
                perc_match = (end - start + 1)*100/(end + left)
                if perc_match < 80:
                    continue
                results[fragname] = {
                    'score': vals[1],
                    'start': start,
                    'end': end,
                    'left': left,
                    'primer': vals[9],
                    'perc_match': perc_match,
                }
    return results

def trim_primers(fqfile, alignout, max_trim_len, primerinfo, outlabel, logfh, 
                 args):
    """Returns trimmed fastq file and file with list of sequence names"""
    trimmedfq = outlabel + ".trimmed.fastq"
    seqfile = outlabel + ".seqlist.txt"
    logfh.write("    Trimming fq: {}\n".format(trimmedfq))
    if have_files([trimmedfq, seqfile], args.force, stderr=logfh):
        logfh.write("      Already have {}\n".format(trimmedfq))
        return (trimmedfq, seqfile)
    aligns = parse_alignout(alignout)
    seqlist = []
    with open(trimmedfq, 'w') as outfq:
        inseq = FastQParser(fqfile)
        for seqrec in inseq:
            seqlist.append(seqrec.id)
            if seqrec.id in aligns:
                primer = aligns[seqrec.id]['primer']
                if primerinfo[primer]['overlap']:
                    primerend = aligns[seqrec.id]['end'] +\
                                aligns[seqrec.id]['left']
                    subrec = seqrec[primerend:]
                    if args.debug:
                        logfh.write("{}\tTrimming\t{}\n".format(
                                         primer, seqrec.id))
                else:
                    if args.debug:
                        logfh.write("{}\tNot trimming\t{}\n".format(
                                         primer, seqrec.id))
                    subrec = seqrec
            else: #trim default max_primer_len+2
                subrec = seqrec[max_trim_len:]
            outfq.write("{}\n".format(subrec.fastq()))
    logfh.write("    Seq list: {}\n".format(seqfile))
    if have_file(seqfile, True, stderr=logfh):
        logfh.write("      Still have {}\n".format(seqfile))
        sys.exit()
    with open_file(seqfile, 'w') as ifh:
        ifh.write("\n".join(seqlist)+"\n")
    return (trimmedfq, seqfile)

def align_and_trim(fqfile, primerfa, primerinfo, max_trim_len, args):
    outlabel = fastq_file_label(fqfile, args.outdir)
    loglabel = fastq_file_label(fqfile, args.logdir) if args.logdir \
               else outlabel
    logfile = loglabel + ".primer_trim.log"
    sys.stderr.write("  Trimming {}\n".format(os.path.basename(fqfile)))
    sys.stderr.write("    Log file {}\n".format(logfile))
    logfh = open(logfile, 'w')
    try:
        logfh.write("Start time: {}\n".format(timestamp()))
        fafiles = create_fasta_of_primer_region(fqfile, max_trim_len, 
                                                outlabel, logfh, args.force)
        alignout = run_aligner(fafiles, primerfa, outlabel, logfh, args.force)
        (trimfq, seqnamefile) = trim_primers(fqfile, alignout, max_trim_len, 
                                            primerinfo, outlabel, logfh, args)
        if os.path.isfile(trimfq):
            sys.stderr.write("    Created {}\n".format(trimfq))
        logfh.write("Finish time: {}\n".format(timestamp()))
    except Exception, e:
        e.args += (fqfile, )
        sys.stderr.write("Error trimming primers." +\
                         "  Check log file {}\n".format(logfile))
        raise
    finally:
        logfh.close()
    files = { 'fafiles': fafiles, 
              'alignfile': alignout, 
              'seqnamefile': seqnamefile, 
              'trimmed': trimfq }
    return (fqfile, files)

@profile
def align_trim_all_fqfiles(fqfiles, primerfa, primerinfo, 
                                max_primer_len, args):
    numfiles = len(fqfiles)
    max_trim_len = max_primer_len+2
    sys.stderr.write("\nTrimming {} fastq files using {}.\n".format(numfiles,
           "{} process{}".format(args.processes, 'es' if args.processes>1 
           else '')))
    files = {}
    pool = mp.Pool(processes=args.processes)
    results = [ pool.apply_async(align_and_trim, args=[fqfile, primerfa,
                    primerinfo, max_trim_len, args]) for fqfile in fqfiles ]
    files = dict([p.get() for p in results ])
    return files

def pair_fqfiles(fqfiles, test=False):
    """Pair a list of fastq files based on name.  Expects Illumina naming 
       convention: samplename_L001_R[12]_001.fastq."""
    fqpairs = {}
    for fqfile in sorted(fqfiles):
        base = os.path.basename(fqfile)
        sample = base.split("_L00", 1)[0].split("_R1", 1)[0]\
                       .split("_R2", 1)[0]
        if not sample in fqpairs:
            fqpairs[sample] = []
        fqpairs[sample].append(fqfile)
    for sample in sorted(fqpairs.keys()):
        if len(fqpairs[sample]) > 2:
            sys.stderr.write("Bad number of files for {}\n".format(sample))
            sys.stderr.write("  "+"\n  ".join(fqpairs[sample])+"\n")
            sys.exit()
        if test:
            sys.stderr.write("%s\t%s\n" % (sample, fqpairs[sample]))
    return fqpairs

def parse_seqnames(seqnamefile):
    with open(seqnamefile, 'r') as fh:
        names = [ l.rstrip() for l in fh ]
    return names

def estimate_amplicon_size(primers, primerinfo, label=None):
    if primers[0].endswith('_F') and primers[1].endswith('_R'):
        (pf, pr) = primers
    elif primers[1].endswith('_F') and primers[0].endswith('_R'):
        (pr, pf) = primers
    else:
        return 0
    if primerinfo[pf]['strand'] < 0 and primerinfo[pr]['strand'] < 0:
        amplicon_size = primerinfo[pf]['end']-primerinfo[pr]['start']+1
    else:
        amplicon_size = primerinfo[pr]['end']-primerinfo[pf]['start']+1
    if label and amplicon_size < 360:
        sys.stderr.write("\nAmpsize: {} bp\t{}\n".format(amplicon_size, label))
        sys.stderr.write("{}\t{}\n".format(pf, primerinfo[pf]))
        sys.stderr.write("{}\t{}\n".format(pr, primerinfo[pr]))
    return amplicon_size

def pair_status(primers, primerinfo):
    if len(primers) == 1:
        return ('singleton', '')
    p1 = primers[0]
    p2 = primers[1]
    if p1 == '' and p2 == '':
        return ('unidentified', '')
    elif p1 == '' or p2 == '':
        return ('singleton', '')
    elif ((p1.endswith('_F') and p2.endswith('_F')) or
       (p1.endswith('_R') and p2.endswith('_R'))): 
           return ('misprime', '')
    ampsize = estimate_amplicon_size(primers, primerinfo)
    if 50 < ampsize < 1000:
        pname1 = p1.rstrip('_F').rstrip('_R') 
        pname2 = p2.rstrip('_F').rstrip('_R') 
        if pname1 == pname2:
            return ('paired-good', ampsize)
        else:
            return ('paired-other', ampsize)
    else:
        return ('misprime', '')

def arrange_by_fragments(fqpair, outfiles, primerinfo, logfh):
    """Create file listing each fragment and the primers each read of the
    fragment mapped to."""
    frag2primers = defaultdict(lambda: {'primers':{}, 'status':'', 
                               'ampsize':''})
    seqnames = {}
    for fqfile in fqpair:
        alignout = outfiles[fqfile]['alignfile']
        seqnamefile = outfiles[fqfile]['seqnamefile']
        alignresults = parse_alignout(alignout)
        seqnames[fqfile] = parse_seqnames(seqnamefile)
        numseqs = len(seqnames[fqfile])
        logfh.write("  {}: {} seqs\n".format(fqfile, numseqs))
        for fragname in seqnames[fqfile]:
            if fragname in alignresults:
                primer = alignresults[fragname]['primer']
                frag2primers[fragname]['primers'][fqfile] = primer
            else:
                frag2primers[fragname]['primers'][fqfile] = ''
    for fragname in frag2primers:
        primers = frag2primers[fragname]['primers'].values()
        (status, ampsize) = pair_status(primers, primerinfo)
        frag2primers[fragname]['status'] = status
        frag2primers[fragname]['ampsize'] = ampsize
    # check that all fqfiles have same number of seqs
    numseqs = 0
    for fqfile in fqpair:
        if not numseqs:
            numseqs = len(seqnames[fqfile])
        else:
            if not len(seqnames[fqfile]) == numseqs:
                err = "Different number of seqs.\n" +\
                      "   {} seqs: {}\n".format(numseqs, fqpair[0]) +\
                      "   {} seqs: {}\n".format(len(seqnames[fqfile]), fqfile)
                logfh.write(err)
                sys.stderr.write(err)
                sys.exit()
    return frag2primers

def create_fragment_report(fqpair, frag2primers, read_primer_file, force, 
                           logfh):
    logfh.write("\nCreating fragment report\n".format(read_primer_file))
    if have_file(read_primer_file, force, stderr=logfh):
        logfh.write("  Already have {}\n".format(read_primer_file))
        return
    readnum_patt = re.compile('.*_(R[12]).*')
    readnums = [ readnum_patt.sub('\\1_primer', fqfile) for fqfile in fqpair ]
    logfh.write("  Writing {}\n".format(read_primer_file))
    fragcounts = { 'tot_fragments':0, 'unidentified': 0, 'singleton': 0,
                   'paired-good': 0, 'paired-other': 0, 'misprime': 0, }
    with open(read_primer_file, 'w') as ofh:
        ofh.write("Fragment\t" + "\t".join(readnums) +\
                  "\tStatus\tEstimated_fragment_size\n")
        for fragname in sorted(frag2primers.keys()):
            row = [fragname,]
            for fqfile in fqpair:
                if fqfile in frag2primers[fragname]['primers']:
                    row.append(frag2primers[fragname]['primers'][fqfile])
                else:
                    row.append('')
            status = frag2primers[fragname]['status']
            fragcounts[status] += 1
            row.extend([status, frag2primers[fragname]['ampsize']])
            ofh.write("\t".join([str(r) for r in row])+"\n")
    fragcounts['tot_fragments'] = len(frag2primers.keys())
    logfh.write("{:<12}\t{:<6} fragments\t{:<5}%\n".format("Fragment status", 
                "Number", "Percent"))
    for k in fragcounts:
        perc = fragcounts[k]*100.0/fragcounts['tot_fragments']
        logfh.write("{:<12}\t{:>6} fragments\t{:5.1f}%\n".format(k, 
                    fragcounts[k], perc))
    logfh.write("{:<12}\t{:>6} fragments\n".format("Total", 
                fragcounts['tot_fragments']))
    fragcounts['paired'] = fragcounts['paired-good'] + \
                           fragcounts['paired-other']
    if os.path.isfile(read_primer_file):
        sys.stderr.write("    Fragment report {}\n".format(read_primer_file))
    return fragcounts

def bin_primer_reads(frag2primers, primerlist, logfh):
    logfh.write("\nCounting read primers\n")
    primerreads = { 'unidentified': [], 'misprime': [] }
    for primer in primerlist:
        primerreads[primer] = []
    readnum_patt = re.compile('.*_(R[12]).*')
    for fragname in frag2primers.keys():
        for fqfile in frag2primers[fragname]['primers']:
            readnum = readnum_patt.sub('\\1', fqfile)
            read = fragname + "/" + readnum
            primer = frag2primers[fragname]['primers'][fqfile]
            status = frag2primers[fragname]['status']
            if status=='unidentified' or status=='misprime':
                primerreads[status].append(read)
            elif primer == '':
                primerreads['unidentified'].append(read)
            else:
                primerreads[primer].append(read)
    return primerreads

def create_primer_report(primerreads, primerlist, primer_read_file, logfh,
                         force=False, debug=False):
    """Create file listing each primer and the reads that match it.
    Include counts of unidentified and mismatched reads.  Tally percent
    of total reads amplified by each primer pair."""
    logfh.write("\nCreating primer report\n".format(primer_read_file))
    if have_file(primer_read_file, force, stderr=logfh):
        logfh.write("  Already have {}\n".format(primer_read_file))
        return
    logfh.write("  Writing {}\n".format(primer_read_file))
    tot_reads = 0
    keylist = primerlist + ['unidentified', 'misprime']
    primercounts = {}
    primerkeys = []
    with open(primer_read_file, 'w') as ofh:
        ofh.write("Primer\tNum_reads\tReads\n")
        for primer in keylist:
            numreads = len(primerreads[primer])
            tot_reads += numreads
            readlist = ", ".join(primerreads[primer])
            primerpair = primer.rstrip('_F').rstrip('_R')
            if primerpair in primercounts:
                primercounts[primerpair] += numreads
            else:
                primerkeys.append(primerpair)
                primercounts[primerpair] = numreads
            ofh.write("{}\t{}\t{}\n".format(primer, numreads, readlist))
    if debug:
        logfh.write("{:<12}\t{:<6} reads\t{:<5}%\n".format(
            "Primer", "Number", "Percent"))
        for k in keylist:
            numreads = len(primerreads[k])
            perc = numreads*100.0/tot_reads
            logfh.write("{:<12}\t{:>6} reads\t{:5.2f}%\n".format(k, 
                numreads, perc))
        logfh.write("{:<12}\t{:<6} reads\n".format("Total", tot_reads))
    primercounts['tot_reads'] = tot_reads
    primerkeys.insert(0, 'tot_reads')
    if os.path.isfile(primer_read_file):
        sys.stderr.write("    Primer read report {}\n".format(primer_read_file))
    return (primercounts, primerkeys)

def assess_sample_primers(sample, fqpair, outfiles, primerlist, primerinfo, 
                          args):
    outlabel = fastq_file_label(sample, args.outdir)
    loglabel = fastq_file_label(sample, args.logdir) if args.logdir \
               else outlabel
    frag_primer_file = outlabel + ".read2primers.txt"
    primer_read_file = outlabel + ".primer2reads.txt"
    logfile = loglabel + ".primer_trim.log"
    logfh = open(logfile, 'w')
    try:
        sys.stderr.write("  Creating primer reports for {}\n".format(sample))
        logfh.write("Sample {}\n".format(sample))
        if have_files([frag_primer_file, primer_read_file], args.force,
                      stderr=logfh):
            logfh.write("  Already have {} and {}\n".format(frag_primer_file,
                        primer_read_file))
        frag2primers = arrange_by_fragments(fqpair, outfiles, primerinfo, 
                                            logfh)
        fragcounts = create_fragment_report(fqpair, frag2primers, 
                                            frag_primer_file, True, logfh)
        primerreads = bin_primer_reads(frag2primers, primerlist, logfh)
        (primercounts, primerkeys) = create_primer_report(primerreads, 
                      primerlist, primer_read_file, logfh, True, args.debug)
    except Exception, e:
        e.args += (sample, )
        raise
    finally:
        logfh.close()
    samplecounts = { 'fragcounts': fragcounts,
                     'primercounts': primercounts,
                     'primerkeys': primerkeys }
    return (sample, samplecounts)

@profile
def assess_all_primers(fqfiles, outfiles, primerinfo, args):
    fqpairs = pair_fqfiles(fqfiles)
    numpairs= len(fqpairs)
    sys.stderr.write("\nAssessing primers in {} samples\n".format(numpairs))
    primerlist = sorted(primerinfo.keys(), 
                        key=lambda x: primerinfo[x]['start'])
    pool = mp.Pool(processes=args.processes)
    results = [ pool.apply_async(assess_sample_primers, args=[sample, 
                fqpairs[sample], outfiles, primerlist, primerinfo, args]) 
                for sample in sorted(fqpairs.keys()) ]
    samplecounts = dict([p.get() for p in results ])
    return samplecounts

def print_summary(samplecounts, summary):
    samples = sorted(samplecounts)
    summfields = ['tot_fragments', 'paired', 'singleton', 'unidentified', 
                  'misprime']
    summfields2 = samplecounts[samples[0]]['primerkeys']
    sys.stderr.write("Writing {}\n".format(summary))
    with open(summary, 'w') as ofh:
        ofh.write("\t"+"\t".join(samples)+"\n")
        for field in summfields:
            row = [ str(samplecounts[s]['fragcounts'][field]) \
                    for s in samples ]
            ofh.write("{}\t{}\n".format(
                      field.title(), "\t".join(row)))
        for field in summfields2:
            row = [ str(samplecounts[s]['primercounts'][field]) \
                    for s in samples ]
            ofh.write("{}\t{}\n".format(field, "\t".join(row)))


if __name__ == '__main__':
    descr = "Trim primers from fastq files and provide statistics."
    parser = ArgumentParser(description=descr)
    parser.add_argument("fqfiles", nargs="+", help="Fastq file(s)")
    parser.add_argument("-s", "--summary", default="summary.txt",
                        help="Name for summary file.")
    parser.add_argument("-o", "--outdir", help="Directory for output files.")
    parser.add_argument("-f", "--force", default=False, action='store_true',
                        help="Overwrite existing files.")
    parser.add_argument("-p", "--processes", default=1, type=int,
                        help="Number of processes to run.")
    parser.add_argument("--debug", default=False, action='store_true',
                        help="Debug mode.")
    parser.add_argument("--logdir", help="Directory for log files.")

    if len(sys.argv)<2:
        parser.print_help()
        sys.exit()
    args = parser.parse_args()

    primerfa = RESOURCE['primer_fa']
    if not os.path.isfile(primerfa):
        sys.exit("Could not find resource {}\n".format(primerfa))
    (primerinfo, max_primer_len) = primer_info(primerfa)
    find_overlaps(primerinfo, args.debug)
    outfiles = align_trim_all_fqfiles(args.fqfiles, primerfa, 
                                      primerinfo, max_primer_len, args)
    if args.outdir: 
        args.summary = os.path.join(args.outdir, args.summary)
    if have_file(args.summary, args.force):
        sys.stderr.write("Already have {}.\n".format(args.summary))
    else:
        samplecounts = assess_all_primers(args.fqfiles, outfiles, primerinfo, 
                                          args)
        print_summary(samplecounts, args.summary)

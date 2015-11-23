#!/usr/bin/env python

import os
import sys
import gzip
import datetime
from contextlib import contextmanager

def open_file(filename, mode='r'):
    """Returns a filehandle to <filename>.  Opens both
       normal and gzipped files.
    """
    if filename.endswith('gz'):
        ifh = gzip.open(filename, mode+'b')
    else:
        ifh = open(filename, mode)
    return ifh

def have_file(filename, rmfile=False, quiet=False, stderr=sys.stderr):
    """Returns True if filename exists and rmfile=False.
       Returns False and removes any existing file if rmfile=True.
       Returns False if filename does not exist.
    """
    if os.path.isfile(filename) and rmfile:
        stderr.write("Removing {}\n".format(filename))
        try:
            os.remove(filename)
        except OSError as e:
            if not quiet:
                stderr.write("Failed to remove {}: {}\n".format(filename, 
                              e.strerror))
    havefile = True if os.path.isfile(filename) else False
    return havefile

def have_files(filenames, rmfile=False, quiet=False, stderr=sys.stderr):
    """filenames: list of files

       Returns True if all files listed in filenames exist and rmfile=False.
       Returns False and removes any existing files if rmfile=True.
       Returns False and removes any existing files if not all files exists and rmfile=False.
    """
    file_exists = []
    file_dne = []
    havefile = True
    for filename in filenames:
        if os.path.isfile(filename):
            file_exists.append(filename)
        else:
            havefile = False
            file_dne.append(filename)
    if file_dne or rmfile:
        for filename in file_exists:
            stderr.write("Removing {}\n".format(filename))
            try:
                os.remove(filename)
                havefile = False
            except OSError as e:
                if not quiet:
                    stderr.write("Failed to remove {}: {}\n".format(filename, 
                                  e.strerror))
    return havefile

def remove_file(filename, quiet=False, stderr=sys.stderr):
    have_file(filename, True, quiet=quiet, stderr=stderr)

def check_directory(subdir, outdir, stderr=sys.stderr):
    if outdir:
        subdir = os.path.join(outdir, subdir)
    if not os.path.isdir(subdir):
        os.mkdir(subdir)
    if not os.path.isdir(subdir):
        stderr.write("Subdirectory {} not found.\n".format(subdir))
        sys.exit(1)
    return subdir


def timestamp():
    return str(datetime.datetime.now()).split('.')[0]

# Returns current datetime as YYMoDD-HHMiSS.
# full=True: return YYYYMoDD-HHMiSS 
# sep='-': set the separator between date and time
# d_only=False: only return the date
# t_only=False: only return the time
def enumdate(dt=None, d_only=False, t_only=False, full=False, sep='-'):
    if not dt: dt = datetime.datetime.now()
    d = None
    t = None
    if d_only or not t_only:
        yr = dt.year
        if not full: yr = yr % 100
        d = "%d%02d%02d" % (yr, dt.month, dt.day)
    if t_only or not d_only:
        t = "%02d%02d%02d" % (dt.hour, dt.minute, dt.second)
    return sep.join([ p for p in [d, t] if p ])



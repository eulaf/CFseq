# Based on Python FastQ parser downloaded from https://scipher.wordpress.com/2010/05/06/simple-python-fastq-parser/

import gzip
import itertools

class FastARead():
   """A record corresponding to a single read and its associated metadata in a FastA file"""
   def __init__(self,seqHeader,seqStr):
       self.header = seqHeader
       self.seq = seqStr
       (seq_id, seq_desc) = seqHeader.split(None, 1)
       self.id = seq_id.lstrip('>')
       self.description = seq_desc
   def toString(self):
       return "%s\n%s" %(self.header,self.seq)

class FastAParser(object):
    """Returns a read-by-read fastA parser analogous to file.readline()"""
    def __init__(self,filePath):
        """Returns a read-by-read fastA parser analogous to file.readline().
        Exmpl: parser.next()
        -OR-
        Its an iterator so you can do:
        for read in parser:
            ... do something with rec ...

        read is a FastARead object with the fields:
        header,seq
        """
        if filePath.endswith('.gz') or filePath.endswith('.bgz'):
            self._file = gzip.open(filePath)
        else:
            self._file = open(filePath, 'rU')
        self._recordNumber = 0
        self._groups = itertools.groupby(self._file, 
                                         lambda l: l.startswith('>'))
        
    def __iter__(self):
        return self
    
    def next(self):
        """Reads in next element, parses, and does minimal verification.
        Returns: FastARead object with fields: (header,seq)"""
        # ++++ Get Next Four Lines ++++
        header = None
        seq = None
        for i in range(2):
            is_header, lines = self._groups.next()
            self._recordNumber += 1
            if is_header:
                header = ''.join([l.rstrip() for l in lines])
            else: 
                seq = ''.join([l.rstrip() for l in lines])
        # -- Make sure we got 4 full lines of data --
        assert header,\
               "** ERROR: It looks like I encountered a premature EOF or empty line.\n\
               Please check FastA file near record number %s and try again**" % (
               self._recordNumber)
        # -- Make sure we are in the correct "register" --
        assert header.startswith('>'),\
               "** ERROR: Record {} does not start with '>'. **".format(
               self._recordNumber)
        assert seq,\
               "** ERROR: No sequence for record number {}.**".format(
               self._recordNumber) 
        
        # ++++ Return fastA data as tuple ++++
        return FastARead(header,seq)

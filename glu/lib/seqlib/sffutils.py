import sys

import xml.etree.cElementTree as etree

from   collections             import defaultdict

import numpy as np

from   Bio               import SeqIO
from   Bio.SeqIO.SffIO   import _sff_read_roche_index_xml as sff_manifest


class SFFIndex(object):
  def __init__(self, sfffiles):
    self.sffindex = defaultdict(list)
    self.headers  = []

    for sfffile in sfffiles:
      self.add_sff(sfffile)

  def add_sff(self, sfffile):
    print >> sys.stderr,'Loading SFF index for',sfffile
    sffindex   = self.sffindex
    headers    = self.headers
    manifest   = sff_manifest(open(sfffile,'rb'))
    reads      = SeqIO.index(sfffile, 'sff')

    for run in etree.fromstring(manifest).findall('run'):
      run_name         = run.find('run_name').text
      run_path         = run.find('path').text
      run_type         = run.find('run_type').text
      accession_prefix = run.find('accession_prefix').text
      analysis_name    = run.find('analysis_name').text

      # Record the read group only the first time we see it
      if accession_prefix not in sffindex:
        read_group = ('@RG','ID:%s' % accession_prefix, 'PL:454', 'SM:%s' % run_name)
        headers.append(read_group)

      # Add this SFF file as a source of reads for this accession prefix
      # Most of the time this will be a unique mapping, but it is possible
      # to split a run into multiple SFF files.  If one was to align based
      # on those files, reads with the same accession will be present in the
      # same input alignment and multiple SFF files.  Fortunately, this
      # should not be typical.
      sffindex[accession_prefix].append(reads)

  def get_read(self, qname):
    prefix = qname[:9]

    for sff in self.sffindex.get(prefix,[]):
      rec = sff.get(qname)
      if rec:
        return rec
    return None


def get_qual(read):
  phred = read.letter_annotations['phred_quality']
  qual  = np.array(phred,dtype=np.uint8)
  qual += 33
  return qual.tostring()


def hard_trim(read,seq,qual,trim):
  seqlen = len(seq)
  hard_trim_left = hard_trim_right = 0

  if 'flowkey' in trim:
    hard_trim_left = 4
  if 'adapter' in trim:
    hard_trim_left = max(hard_trim_left,read.annotations['clip_qual_left'])
  if 'lowquality' in trim:
    hard_trim_right = seqlen-read.annotations['clip_qual_right']

  seq  =  seq[hard_trim_left:seqlen-hard_trim_right]
  qual = qual[hard_trim_left:seqlen-hard_trim_right]

  return seq,qual,hard_trim_left,hard_trim_right



# -*- coding: utf-8 -*-

from __future__ import division

__gluindex__  = True
__abstract__  = 'Convert Newbler 454PairAlign.txt and the corresponding SFF files into SAM/BAM format'
__copyright__ = 'Copyright (c) 2010, BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.'
__license__   = 'See GLU license for terms by running: glu license'
__revision__  = '$Id$'

import re
import sys

from   operator          import getitem, itemgetter
from   itertools         import izip, imap, groupby, repeat
from   subprocess        import Popen, PIPE

import numpy as np

from   glu.lib.fileutils import autofile,hyphen,table_writer,table_reader


CIGAR_map = { ('-','-'):'P' }
for a in 'NACGTacgt':
  CIGAR_map[a,'-'] = 'I'
  CIGAR_map['-',a] = 'D'
  for b in 'NACGTacgt':
    CIGAR_map[a,b] = 'M'

NDIFF_map = { ('-','-'): ('P',None) }
for a in 'NACGTacgt':
  NDIFF_map[a,'-'] = ('I',a)
  NDIFF_map['-',a] = ('D',a)
  # N,N is considered a mismatch(?)
  for b in 'NACGTacgt':
    NDIFF_map[a,b] = ('=' if a==b and a!='N' else 'X',b)


def make_cigar_py(query,ref,qstart,qstop,qlen):
  assert len(query)==len(ref)
  igar  = imap(getitem, repeat(CIGAR_map), izip(query,ref))
  cigar = ''.join('%d%s' % (len(list(run)),code) for code,run in groupby(igar))

  if qstart<qstop:
    clip1 = '%dH' % (qstart-1)    if qstart!=1    else ''
    clip2 = '%dH' % (qlen-qstop)  if qstop !=qlen else ''
    cigar = clip1+cigar+clip2
  else:
    clip1 = '%dH' % (qstop-1)     if qstop !=1    else ''
    clip2 = '%dH' % (qlen-qstart) if qstart!=qlen else ''
    cigar = clip2+cigar+clip1

  return cigar


def make_ndiff_py(query,ref):
  assert len(query)==len(ref)

  nd   = groupby(imap(getitem, repeat(NDIFF_map), izip(query,ref)),itemgetter(0))
  nm   = 0
  md   = ''
  eq   = 0

  for code,items in nd:
    if code=='P':
      continue
    elif code=='=':
      eq += len(list(items))
    elif code=='I':
      nm += len(list(items))
    else:
      md += '%d' % eq
      eq = 0

      bases = ''.join(b for c,b in items)
      nm   += len(bases)

      if code=='D':
        md   += '^'
      else:
        # For some silly reason mismatch bases must be 0 separated
        bases = '0'.join(bases)

      md += bases

  md += '%d' % eq

  return nm,md


# Try to import the optimized Cython version
# The Python version is pretty fast, but I wanted to play with Cython.
try:
  from samhelpers import make_cigar, make_ndiff
except ImportError:
  make_cigar = make_cigar_py
  make_ndiff = make_ndiff_py


class SFFIndex(object):
  def __init__(self, sfffiles):
    self.sffindex = {}
    self.headers  = []

    for sfffile in sfffiles:
      self.add_sff(sfffile)

  def add_sff(self, sfffile):
    from   Bio             import SeqIO
    from   Bio.SeqIO.SffIO import _sff_read_roche_index_xml as sff_manifest
    import xml.etree.cElementTree as etree

    print >> sys.stderr,'Loading SFF index for',sfffile
    sffindex   = self.sffindex
    headers    = self.headers
    manifest   = sff_manifest(open(sfffile,'rb'))
    reads      = SeqIO.index(sfffile, 'sff-trim')

    for run in etree.fromstring(manifest).findall('run'):
      run_name         = run.find('run_name').text
      run_path         = run.find('path').text
      run_type         = run.find('run_type').text
      accession_prefix = run.find('accession_prefix').text
      analysis_name    = run.find('analysis_name').text
      if accession_prefix in sffindex:
        raise ValueError('Duplicate SFF accession prefix found.  Reads may be not be unique.')
      read_group = ('@RG','ID:%s' % accession_prefix, 'PL:454', 'SM:%s' % run_name)
      headers.append(read_group)
      sffindex[accession_prefix] = reads

  def get_quality(self, qname, query, qstart, qstop):
    prefix = qname[:9]
    sff = self.sffindex.get(prefix)

    if not sff:
      return '*'

    rec      = sff[qname]
    phred    = rec.letter_annotations['phred_quality']
    sffqual  = np.array(phred,dtype=np.uint8)
    sffqual += 33
    sffqual  = sffqual.tostring()

    # Align the query to the original read to find the matching quality
    # score information.  This is complicated by the extra trimming done by
    # gsMapper.  We could obtain this information by parsing the
    # 454TrimStatus.txt, but it is easier to search for the sub-sequence in
    # the reference.  Ones hopes the read maps uniquely, but this is not
    # checked.

    # CASE 1: Forward read alignment
    if qstart<qstop:
      # Try using specified cut-points
      read = str(rec.seq)
      seq  = read[qstart-1:qstop]

      # If it matches, then compute quality
      if seq==query:
        qual  = sffqual[qstart-1:qstop]
      else:
        # otherwise gsMapper applied extra trimming, so we have to manually find the offset
        start = read.index(query)
        seq   = read[start:start+len(query)]
        if seq==query:
          #print >> sys.stderr,'MATCHED TYPE F2: name=%s, qstart=%d(%d), qstop=%d, qlen=%d, len.query=%d' % (qname,start+1,qstart,qstop,qlen,len(query))
          qual  = sffqual[start:start+len(query)]

    # CASE 2: Backward read alignment
    else:
      # Try using specified cut-points
      read = str(rec.seq.complement())
      seq  = read[qstop-1:qstart][::-1]
      read = read[::-1]

      # If it matches, then compute quality
      if seq==query:
        qual  = sffqual[qstop-1:qstart][::-1]
      else:
        # otherwise gsMapper applied extra trimming, so we have to manually find the offset
        start = read.index(query)
        seq   = read[start:start+len(query)]

        if seq==query:
          #print >> sys.stderr,'MATCHED TYPE R2: name=%s, qstart=%d, qstop=%d(%d), qlen=%d, len.query=%d' % (qname,qstart,start+1,qstop,qlen,len(query))
          qual = sffqual[::-1][start:start+len(query)]

    assert seq==query
    assert len(qual) == len(query)
    return qual


def pair_align_records(records, sffindex):
  split   = re.compile('[\t ,.]+').split
  mrnm    = '*'
  mpos    = 0
  isize   = 0
  mapq    = 60
  flag    = 0   # query is normalized by gsMapper to forward reference strand

  for line in records:
    assert line.startswith('>')

    fields = split(line)
    qname  = fields[0][1:]
    qstart = int(fields[1])
    qstop  = int(fields[2])
    qlen   = int(fields[4])
    rname  = fields[6]
    rstart = fields[7]       # Keep as string
    #rstop  = fields[8]      # Unused
    #rlen  = int(fields[10]) # Unused
    query  = split(records.next())[2]
    qq     = query.replace('-','')
    ref    = split(records.next())[2]
    cigar  = make_cigar(query, ref, qstart, qstop, qlen)
    qual   = sffindex.get_quality(qname, qq, qstart, qstop)
    nm,md  = make_ndiff(query,ref)
    opt    = ['NM:i:%d' % nm,'MD:Z:%s' % md]

    if qual!='*':
      opt.append('RG:Z:%s' % qname[:9])

    flag   = 0
    if qstop<qstart:
      flag |= 0x10

    yield [qname, flag, rname, rstart, mapq, cigar, mrnm, mpos, isize, qq, qual]+opt


def fixup_mapq(alignment):
  for qname,qalign in groupby(alignment,itemgetter(0)):
    qalign = list(qalign)

    if len(qalign)>1:
      # Set MAPQ to 1 for multiply aligned reads
      for row in qalign:
        row[4] = 1
        yield row
    else:
      yield qalign[0]


def load_contig_remap(filename):
  remap = {}
  data  = table_reader(filename, columns=['SRC_CONTIG','DST_CONTIG','OFFSET'],want_header=False)

  for src,dst,offset in data:
    offset = int(offset or 0)
    if src!=dst or offset:
      remap[src] = dst,offset

  return remap


def remap_contigs(alignment, remap):
  for row in alignment:
    contig = row[2]
    if contig in remap:
      dst,offset = remap[contig]
      row[2] = dst
      if offset:
        row[3] = int(row[3])+offset
    yield row


def option_parser():
  import optparse

  usage = 'usage: %prog [options] 454PairAlign.txt[.gz] [SFFfiles.sff..]'
  parser = optparse.OptionParser(usage=usage)

  parser.add_option('-r', '--reflist', dest='reflist', metavar='FILE',
                    help='Reference genome contig list')
  parser.add_option('-m', '--remapcontig', dest='remapcontig', metavar='FILE',
                    help='Contig remapping')
  parser.add_option('-o', '--output', dest='output', metavar='FILE', default='-',
                    help='Output SAM file')
  return parser


def main():
  parser = option_parser()
  options,args = parser.parse_args()

  if not args:
    parser.print_help(sys.stderr)
    sys.exit(2)

  alignfile = autofile(hyphen(args[0],sys.stdin))
  sffindex  = SFFIndex(args[1:])

  write_bam = options.output.endswith('.bam')

  if write_bam:
    if not options.reflist:
      raise ValueError('Conversion to BAM format requires a reference genome contig list (-r/--reflist)')

    # Creating the following two-stage pipeline deadlocks due to problems with subprocess
    # -- use the shell method below instead
    #sammer = Popen(['samtools','import',options.reflist,'-','-'],stdin=PIPE,stdout=PIPE,bufsize=-1)
    #bammer = Popen(['samtools','sort','-', options.output[:-4]], stdin=sammer.stdout,bufsize=-1)

    cmd    = 'samtools import "%s" - - | samtools sort - "%s"' % (options.reflist,options.output[:-4])
    bammer = Popen(cmd,stdin=PIPE,shell=True,bufsize=-1)
    out    = table_writer(bammer.stdin)
  else:
    out = table_writer(options.output,hyphen=sys.stdout)

  out.writerow(['@HD', 'VN:1.0'])

  if options.reflist:
    reflist = table_reader(options.reflist)
    for row in reflist:
      if len(row)<2:
        continue

      contig_name = row[0]
      contig_len  = int(row[1])

      out.writerow(['@SQ', 'SN:%s' % contig_name, 'LN:%d' % contig_len])

  rg_seen = set()
  for row in sffindex.headers:
    if row[1] not in rg_seen:
      rg_seen.add(row[1])
      out.writerow(row)

  print >> sys.stderr, 'Generating alignment from %s to %s' % (args[0],options.output)
  alignment = fixup_mapq(pair_align_records(alignfile, sffindex))

  if options.remapcontig:
    alignment = remap_contigs(alignment, load_contig_remap(options.remapcontig))

  out.writerows(alignment)

  if write_bam:
    print >> sys.stderr,'Finishing BAM encoding...'
    bammer.communicate()


if __name__=='__main__':
  if 1:
    main()
  else:
    try:
      import cProfile as profile
    except ImportError:
      import profile
    import pstats

    prof = profile.Profile()
    try:
      prof.runcall(main)
    finally:
      stats = pstats.Stats(prof)
      stats.strip_dirs()
      stats.sort_stats('time', 'calls')
      stats.print_stats(25)

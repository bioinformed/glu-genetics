# -*- coding: utf-8 -*-

from __future__ import division

__gluindex__  = True
__abstract__  = 'Convert Newbler 454PairAlign.txt and the corresponding SFF files into SAM/BAM format'
__copyright__ = 'Copyright (c) 2010, BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.'
__license__   = 'See GLU license for terms by running: glu license'
__revision__  = '$Id$'

import re
import sys
import random

import xml.etree.cElementTree as etree

from   operator          import getitem, itemgetter
from   itertools         import izip, imap, groupby, repeat
from   collections       import defaultdict
from   subprocess        import Popen, PIPE

import numpy as np

from   glu.lib.fileutils import autofile,hyphen,table_writer,table_reader

from   Bio               import SeqIO
from   Bio.SeqIO.SffIO   import _sff_read_roche_index_xml as sff_manifest


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


def make_cigar_py(query,ref):
  assert len(query)==len(ref)
  igar  = imap(getitem, repeat(CIGAR_map), izip(query,ref))
  cigar = ''.join('%d%s' % (len(list(run)),code) for code,run in groupby(igar))
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


def cigar_add_trim(cigar, trim_char, left, right):
  if left and right:
    cigar = '%d%s%s%d%s' % (left,trim_char,cigar,right,trim_char)
  elif left:
    cigar = '%s%s%s' % (left,trim_char,cigar)
  elif right:
    cigar = '%s%d%s' % (cigar,right,trim_char)
  return cigar


def pair_align_records(records, trim, sffindex):
  split   = re.compile('[\t ,.]+').split
  mrnm    = '*'
  mpos    = 0
  isize   = 0
  mapq    = 60
  flag    = 0   # query is normalized by gsMapper to forward reference strand

  trim_unaligned = 'unaligned' in trim

  # trimming unaligned bases superceeds all other trimming options
  if trim_unaligned:
    trim = False

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
    qseq   = query.replace('-','')
    ref    = split(records.next())[2]

    flag   = 0 if qstart<qstop else 0x10
    cigar  = make_cigar(query, ref)
    nm,md  = make_ndiff(query,ref)
    opt    = ['NM:i:%d' % nm,'MD:Z:%s' % md]

    read   = sffindex.get_read(qname)

    # Read not found (or SFF files not provided)
    if not read:
      seq  = qseq
      qual = '*'

      # No soft trimming, since all we have is the aligned portion of the read
      soft_trim_left = soft_trim_right = 0

      # Add hard trimming information, since we are told if any of the read
      # was trimmed by the aligner.  We do not know how much was trimmed for
      # the flowkey, adapter, and due to low base quality.
      if qstart<=qstop:
        hard_trim_left,hard_trim_right = qstart-1,qlen-qstop
      else:
        hard_trim_left,hard_trim_right = qlen-qstart,qstop-1

    # Read was found in one of the supplied SFF files
    else:
      seq    = read.seq
      qual   = get_qual(read)

      # Since we know the read name is valid, we can annotate that it came
      # from a read group for the region
      opt.append('RG:Z:%s' % qname[:9])

      # Add any requested hard trimming
      if trim:
        seq,qual,hard_trim_left,hard_trim_right = hard_trim(read,seq,qual,trim)
      else:
        hard_trim_left = hard_trim_right = 0

      # If aligned in reverse, flip everything to the forward strand
      if qstart>qstop:
        hard_trim_left,hard_trim_right = hard_trim_right,hard_trim_left
        seq  = seq.reverse_complement()
        qual = qual[::-1]

      # Align the query to the original read to find the matching quality
      # score information.  This is complicated by the extra trimming done by
      # gsMapper.
      seq    = str(seq)
      start  = seq.index(qseq)

      # Hard trim everything but the aligned portion, if requested
      if trim_unaligned:
        end  = start+len(qseq)
        hard_trim_left  += start
        hard_trim_right += len(seq)-end
        soft_trim_left,soft_trim_right = 0,0
        seq  =  seq[start:end]
        qual = qual[start:end]

      # Compute soft-trimming
      else:
        soft_trim_left,soft_trim_right = start,len(seq)-len(qseq)-start

    # Add soft and hard trimming codes to the CIGAR
    cigar = cigar_add_trim(cigar, 'S', soft_trim_left, soft_trim_right)
    cigar = cigar_add_trim(cigar, 'H', hard_trim_left, hard_trim_right)

    # Add custom tags to allow trimming downstream by indicating the aligned
    # read start and end position
    if soft_trim_left or soft_trim_right:
      opt += ['ZS:i:%d' % (soft_trim_left+1), 'ZE:i:%d' % (len(seq)-soft_trim_right-soft_trim_left) ]

    yield [qname, flag, rname, rstart, mapq, cigar, mrnm, mpos, isize, seq, qual]+opt


def handle_maligned(alignment, malign_action, pick_method, reads_seen):
  malign_action = malign_action.lower()
  pick_method   = pick_method.lower()

  if pick_method=='random':
    def _pick(aligns):
      i = random.randint(0,len(aligns)-1)
      aligns[0],aligns[i] = aligns[i],aligns[0]
      return aligns

  elif pick_method=='best':
    # Use sum of qualities?
    def _pick(aligns):
      aligns.sort(key=lambda a: (-len(a[9]),int(a[11].split(':')[-1])))
      return aligns

  else:
    raise ValueError('Invalid primary alignment action: %s' % pick_method)

  alignment = groupby(alignment,itemgetter(0))

  if malign_action=='keep-primary':
    for name,aligns in alignment:
      aligns = list(aligns)

      if len(aligns)>1:
        # Find longest read with fewest mismatches, set MAPQ to 1
        # and yield
        align    = _pick(aligns)[0]
        align[4] = 1
        yield align
      else:
        yield aligns[0]

  elif malign_action=='keep-all':
    for name,aligns in alignment:
      aligns = list(aligns)

      if len(aligns)>1:
        # Find longest read with fewest mismatches, set MAPQ to 1,
        # set all subsequent reads as non-primary, and yield
        aligns = _pick(aligns)

        for i,align in enumerate(aligns):
          align[4] = 1
          if i:
            align[1] |= 0x100

          yield align
      else:
        yield aligns[0]

  elif malign_action=='unalign':
    keep_opts = set(['RG','ZS','ZE'])

    for name,aligns in alignment:
      aligns = list(aligns)

      if len(aligns)>1:
        # Skip manually unaligning, since it introduces several
        # complications and can be done more better in handle_unaligned
        # FIXME: What happens when no SFF file is provided?
        if 0:
          align = _pick(aligns)[0]

          align[1]  = 0x04
          align[2]  = '*'
          align[3]  = 0
          align[4]  = '*'

          yield align[:11] + [ a for a in align[11:] if a[:2] in keep_opts ]

      else:
        yield aligns[0]

  elif malign_action=='drop':
    for name,aligns in alignment:
      aligns = list(aligns)

      if len(aligns)>1:
        # When non-unique, mark each alignment as seen, but do not yield.
        # This ensures that the read is excluded as from being re-added as
        # unaligned.
        for align in aligns:
          reads_seen.add(align[0])
      else:
        yield aligns[0]

  else:
    raise ValueError('Invalid multiply aligned read action: %s' % malign_action)


def handle_unaligned(alignment, action, trim, reads_seen, sfffiles):
  # Pass-through all reads and mark them as seen
  action = action.lower()

  if action=='keep':
    for align in alignment:
      reads_seen.add(align[0])
      yield align

  elif action=='drop':
    for align in alignment:
      if not align[1]&0x04:
        yield align

  else:
    raise ValueError('Invalid unaligned read action: %s' % action)

  if not sfffiles or action=='drop':
    return

  flag   = 0x04
  rname  = '*'
  pos    = 0
  mapq   = 0
  cigar  = '*'
  mrnm   = '*'
  mpos   = 0
  isize  = 0

  aligned_count   = len(reads_seen)
  unaligned_count = 0

  for filename in sfffiles:
    for read in SeqIO.parse(open(filename,'rb'), 'sff'):
      qname = read.id

      if qname in reads_seen:
        continue

      reads_seen.add(qname)

      seq    = str(read.seq)
      qual   = get_qual(read)
      rg     = 'RG:Z:%s' % qname[:9]
      start  = read.annotations['clip_qual_left']
      end    = read.annotations['clip_qual_right']

      if trim:
        seq,qual,hard_trim_left,hard_trim_right = hard_trim(read,seq,qual,trim)
        start -= hard_trim_left
        end   -= hard_trim_left

      unaligned_count += 1

      yield [qname, flag, rname, pos, mapq, cigar, mrnm, mpos, isize, seq, qual,
                    'ZS:i:%d' % (start+1), 'ZE:i:%d' % end, rg]

  print >> sys.stderr,'Processed %d aligned reads and %d unaligned reads.' % (aligned_count,unaligned_count)


def load_contig_remap(filename):
  remap = {}
  data  = table_reader(filename, columns=['SRC_CONTIG','DST_CONTIG','OFFSET'],want_header=False)

  for src,dst,offset in data:
    offset = int(offset or 0)
    if src!=dst or offset:
      remap[src] = dst,offset

  return remap


def remap_contigs(alignment, remap):
  for align in alignment:
    contig = align[2]
    if contig in remap:
      dst,offset = remap[contig]
      align[2] = dst
      if offset:
        align[3] = int(align[3])+offset
    yield align


def validate_trimming(options):
  options.trim = set(f.strip().lower() for f in options.trim.split(','))
  options.trim.discard('none')

  all_trim     = set(['unaligned', 'flowkey', 'adapter', 'lowquality'])
  allowed_trim = all_trim|set(['all'])
  extra_trim   = options.trim - allowed_trim
  if extra_trim:
    raise ValueError('Invalid read trim options specified: %s' % ', '.join(sorted(extra_trim)))

  if options.trim&all_trim == len(all_trim):
    options.trim.add('all')
  elif 'all' in options.trim:
    options.trim = allowed_trim


def option_parser():
  import optparse

  usage = 'usage: %prog [options] 454PairAlign.txt[.gz|.bz2] [SFFfiles.sff..]'
  parser = optparse.OptionParser(usage=usage)

  parser.add_option('--reflist', dest='reflist', metavar='FILE',
                    help='Reference genome contig list')
  parser.add_option('--remapcontig', dest='remapcontig', metavar='FILE',
                    help='Contig remapping')
  parser.add_option('--trim', dest='trim', metavar='ACTION', default='all',
                    help='Trim feature(s) of reads.  Comma separated list of: flowkey, adapter, lowquality, unaligned, all.  Default=all')
  parser.add_option('--maligned', dest='maligned', metavar='ACTION', default='keep-all',
                    help='Action to perform for multiply aligned reads: keep-primary, keep-all, unalign, drop.  Default=keep-all')
  parser.add_option('--mpick', dest='mpick', metavar='METHOD', default='best',
                    help='Method of selecting primary alignment when keeping multiply aligned reads: best, random.  Default=best')
  parser.add_option('--unaligned', dest='unaligned', metavar='ACTION', default='keep',
                    help='Action to perform for unaligned reads: keep, drop.  Default=keep')
  parser.add_option('-o', '--output', dest='output', metavar='FILE', default='-',
                    help='Output SAM file')
  return parser


def main():
  parser = option_parser()
  options,args = parser.parse_args()

  if not args:
    parser.print_help(sys.stderr)
    sys.exit(2)

  validate_trimming(options)

  alignfile = autofile(hyphen(args[0],sys.stdin))

  sfffiles  = args[1:]
  sffindex  = SFFIndex(sfffiles)

  write_bam = options.output.endswith('.bam')

  if write_bam:
    if not options.reflist:
      raise ValueError('Conversion to BAM format requires a reference genome contig list (--reflist)')

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
  reads_seen = set()
  alignment  = pair_align_records(alignfile, options.trim, sffindex)
  alignment  = handle_maligned(alignment,  options.maligned,  options.mpick, reads_seen)
  alignment  = handle_unaligned(alignment, options.unaligned, options.trim,  reads_seen, sfffiles)

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

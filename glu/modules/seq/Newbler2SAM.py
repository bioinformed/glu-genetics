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

from   operator                import itemgetter
from   itertools               import groupby
from   subprocess              import Popen, PIPE

from   Bio                     import SeqIO

from   glu.lib.utils           import namedtuple
from   glu.lib.fileutils       import autofile,hyphen,table_writer,table_reader

from   glu.lib.seqlib.sffutils import SFFIndex, get_qual, hard_trim
from   glu.lib.seqlib.cigar    import make_cigar, make_ndiff, cigar_add_trim


PAIR_ALIGN_TAB_HEADER = ['QueryAccno', 'QueryStart', 'QueryEnd', 'QueryLength',
                         'SubjAccno',  'SubjStart',  'SubjEnd',  'SubjLength',
                         'NumIdent', 'AlignLength', 'QueryAlign', 'SubjAlign']

AlignRecord = namedtuple('AlignRecord', ' '.join(PAIR_ALIGN_TAB_HEADER))


def read_pair_align(alignfile):
  line   = alignfile.next()
  fields = line.strip().split('\t')

  # If in tab-delimited format
  if fields == PAIR_ALIGN_TAB_HEADER:
    for line in alignfile:
      fields     = line.split('\t')
      fields[1]  = int(fields[1])
      fields[2]  = int(fields[2])
      fields[3]  = int(fields[3])
      fields[5]  = int(fields[5])
      fields[6]  = int(fields[6])
      fields[7]  = int(fields[7])
      fields[8]  = int(fields[8])
      fields[9]  = int(fields[9])
      fields[11] = fields[11].strip()

      yield AlignRecord._make(fields)

  # human-readable format
  elif line.startswith('>'):
    split = re.compile('[\t ,.()/]+').split

    def _record(line):
      fields = split(line)
      qname  = fields[0][1:]
      qstart = int(fields[1])
      qstop  = int(fields[2])
      qlen   = int(fields[4])
      rname  = fields[6]
      rstart = int(fields[7])
      rend   = int(fields[8])
      rlen   = int(fields[10])
      ident  = int(fields[11])
      alen   = int(fields[12])
      query  = split(alignfile.next())[2]
      ref    = split(alignfile.next())[2]

      return AlignRecord(qname,qstart,qstop,qlen,rname,rstart,rend,rlen,ident,alen,query,ref)

    yield _record(line)
    for line in alignfile:
      yield _record(line)

  else:
    raise ValueError('Unknown pairwise alignment format')


def pair_align_records(alignments, trim, sffindex):
  mrnm    = '*'
  mpos    = 0
  isize   = 0
  mapq    = 60
  flag    = 0   # query is normalized by gsMapper to forward reference strand

  trim_unaligned = 'unaligned' in trim

  # trimming unaligned bases superceeds all other trimming options
  if trim_unaligned:
    trim = False

  for align in alignments:
    qname  = align.QueryAccno
    qstart = align.QueryStart
    qstop  = align.QueryEnd
    qlen   = align.QueryLength
    rname  = align.SubjAccno
    rstart = align.SubjStart
    query  = align.QueryAlign
    qseq   = query.replace('-','')
    ref    = align.SubjAlign

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
      # FIXME: SeqIO records already have quality scores!!!
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
  from glu.lib.glu_argparse import GLUArgumentParser

  parser = GLUArgumentParser(description=__abstract__)

  parser.add_argument('pairalign',             help='454PairAlign.txt[.gz|.bz2] file')
  parser.add_argument('sfffile',   nargs='*',  help='Simple Flowgram Format (SFF) files')

  parser.add_argument('--reflist', metavar='FILE',
                    help='Reference genome contig list')
  parser.add_argument('--remapcontig', metavar='FILE',
                    help='Contig remapping')
  parser.add_argument('--trim', metavar='ACTION', default='all',
                    help='Trim feature(s) of reads.  Comma separated list of: flowkey, adapter, lowquality, unaligned, all.  Default=all')
  parser.add_argument('--maligned', metavar='ACTION', default='keep-all',
                    help='Action to perform for multiply aligned reads: keep-primary, keep-all, unalign, drop.  Default=keep-all')
  parser.add_argument('--mpick', metavar='METHOD', default='best',
                    help='Method of selecting primary alignment when keeping multiply aligned reads: best, random.  Default=best')
  parser.add_argument('--unaligned', metavar='ACTION', default='keep',
                    help='Action to perform for unaligned reads: keep, drop.  Default=keep')
  parser.add_argument('-o', '--output', metavar='FILE', default='-',
                    help='Output SAM file')
  return parser


def main():
  parser    = option_parser()
  options   = parser.parse_args()

  validate_trimming(options)

  alignfile = autofile(hyphen(options.pairalign,sys.stdin))
  sffindex  = SFFIndex(options.sfffile)

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

  print >> sys.stderr, 'Generating alignment from %s to %s' % (options.pairalign,options.output)
  reads_seen = set()

  alignment = read_pair_align(alignfile)
  alignment = pair_align_records(alignment, options.trim, sffindex)
  alignment = handle_maligned(alignment,  options.maligned,  options.mpick, reads_seen)
  alignment = handle_unaligned(alignment, options.unaligned, options.trim,  reads_seen, options.sfffile)

  if options.remapcontig:
    alignment = remap_contigs(alignment, load_contig_remap(options.remapcontig))

  out.writerows(alignment)

  if write_bam:
    print >> sys.stderr,'Finishing BAM encoding...'
    bammer.communicate()


if __name__=='__main__':
  main()

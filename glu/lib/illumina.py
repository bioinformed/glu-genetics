# -*- coding: utf-8 -*-

__abstract__  = 'utility functions and objects to read and parse Illumina data files'
__copyright__ = 'Copyright (c) 2009, BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.'
__license__   = 'See GLU license for terms by running: glu license'
__revision__  = '$Id$'

import os
import sys
import csv
import struct

from   glu.lib.utils         import chunk, namedtuple
from   glu.lib.fileutils     import autofile,parse_augmented_filename,guess_format,get_arg,namefile
from   glu.lib.sections      import read_sections

from   glu.lib.genolib.locus import STRAND_UNKNOWN

from   glu.lib.seqlib.strand import norm_snp_seq, complement_base


# Credit for a large proportion of the reverse-engineering in this module
# goes to the R/BioConductor project, who in turn were helped by Keith
# Baggerly.


ManifestRow = namedtuple('ManifestRow',
                         'ilmnid name design_strand alleles assay_type_id norm_id '
                         'addressA_id alleleA_probe_sequence addressB_id alleleB_probe_sequence '
                         'genome_version chromosome mapinfo ploidy species '
                         'source source_version source_strand source_sequence top_genomic_sequence '
                         'customer_strand genomic_strand')

assay_type_map = { 0 : 'Infinium II',
                   1 : 'Infinium I red channel',
                   2 : 'Infinium I green channel' }


IDATData = namedtuple('IDATData', 'filename version snp_count illumina_ids sds means bead_counts '
                                  'midblock red_green manifest barcode format label opa sampleid '
                                  'descr plate well runinfo')


IDAT_FIELD_CODES = { 1000 : 'snp_count',
                      102 : 'illumina_ids',
                      103 : 'sds',
                      104 : 'means',
                      107 : 'bead_counts',
                      200 : 'midblock',
                      300 : 'runinfo',
                      400 : 'red_green',
                      401 : 'manifest',
                      402 : 'barcode',
                      403 : 'format',
                      404 : 'label',
                      405 : 'opa',
                      406 : 'sampleid',
                      407 : 'descr',
                      408 : 'plate',
                      409 : 'well',
                      510 : 'unknown' }


def printable(s):
  return (''.join(map(myord,s or ''))).replace('\n','(\\n)')


def myord(x):
  import string
  if x in string.printable[:-5]:
    return x
  else:
    return '(%02d)' % ord(x)


def int2bin(n, count=32):
  return ''.join([str((n >> y) & 1) for y in range(count-1, -1, -1)])


try:
  from glu.lib._illumina import readstr

except ImportError:
  sys.stderr.write('WARNING: Using slow binary string reader.\n')

  def readstr(afile):
    '''
    String data are encoded as a sequence of one or more length bytes followed
    by the specified number of data bytes.

    The lower 7 bits of each length byte encodes the bits that comprise the
    length of the following byte string.  When the most significant bit it
    set, then an additional length byte follows with 7 additional high bits to
    be added to the current length.  The following string lengths are
    accommodated by increasing sequences of length bytes:

    length  maximum
    bytes   length
    ------  --------
      1       127 B
      2        16 KB
      3         2 MB
      4       256 MB
      5        32 GB

    While this seems like a sensible progression, there is some uncertainty
    about this interpretation, since the longest of string observed in the
    wild has been of length 6,264 with two length bytes.
    '''
    read = afile.read
    n = m = ord(read(1))

    if not n:
      return ''

    if m&0x80:
      shift = 7
      n     = m&0x7F

      while m&0x80:
        m      = ord(read(1))
        n     += (m&0x7F)<<shift
        shift += 7

    return read(n)


class IlluminaManifest(object):
  def __init__(self, filename, extra_args=None, **kwargs):
    if extra_args is None:
      args = kwargs
    else:
      args = extra_args
      args.update(kwargs)

    name   = parse_augmented_filename(filename,args)
    format = get_arg(args, ['format'], guess_format(name, ['csv','bpm'])) or 'csv'

    if extra_args is None and args:
      raise ValueError('Unexpected filename arguments: %s' % ','.join(sorted(args)))

    format = format.lower()
    if format == 'csv':
      self._load_csv(name)
    elif format == 'bpm':
      self._load_bpm(name)
    else:
      raise NotImplementedError("File format '%s' is not a supported Illumina manifest format" % format)

  def _load_csv(self,filename):
    self.filename = filename
    self.version  = 'csv'

    manifest = csv.reader(autofile(filename),dialect='csv')
    sections = read_sections(manifest)

    attrs = dict()
    while 1:
      heading,contents = sections.next()
      attrs.update(c[:2] for c in contents if len(c) > 1 and c[0])
      if 'Assay Format' in attrs:
        break

    def _getattr(names,default):
      for name in names:
        if name in attrs:
          return attrs[name]
      return default

    # Controls may be at the end, but we're going to ignore them for now
    self.controls    = None
    self.snp_count   = int(_getattr(['Loci Count','SNP Count'],0)) or None
    self.snp_entries = None
    self.snp_names   = None
    self.norm_ids    = None

    format = attrs['Assay Format']

    # Old style OPA manifest
    if format == 'Golden Gate':
      pass
    # Infinium or new Golden Gate
    # Known formats: Infinium,Infinium II,Infinium 2,Infinium HD Super,GoldenGate
    elif format.startswith('Infinium') or format == 'GoldenGate':
      heading,contents = sections.next()
      assert heading == 'Assay'
    else:
      raise ValueError('Unknown manifest format: %s' % format)

    def _manifest_rows():
      header = contents.next()
      fields = len(header)

      yield header

      for row in contents:
        n = fields-len(row)
        if n>0:
          row += ['']*n
        yield row

    self.manifest_data = _manifest_rows()

  def _load_bpm(self,filename):
    import numpy as np

    self.filename = filename
    stats         = os.stat(filename)
    filesize      = stats.st_size
    bpm           = autofile(filename,'rb')
    read          = bpm.read

    signature     = read(3)

    if signature != 'BPM':
      raise ValueError('Invalid BPM file signature: %s' % printable(signature))

    self.version = struct.unpack('<BL', read(5))

    if self.version != (1,4):
      raise ValueError('Invalid BPM version number (%d.%d)' % self.version)

    self.manifest_name = readstr(bpm)
    self.controls      = [ c.split(',') for c in readstr(bpm).split('\n') if c ]
    self.snp_count,    = struct.unpack('<L', read(4))
    self.snp_entries   = np.fromfile(bpm, dtype='<u4', count=self.snp_count)
    self.snp_names     = [ readstr(bpm) for i in xrange(self.snp_count) ]
    self.norm_ids      = np.fromfile(bpm, dtype='u1', count=self.snp_count)

    def _manifest_rows():
      header = ['IlmnID','Name','IlmnStrand','SNP', 'AssayTypeID', 'NormID',
                'AddressA_ID','AlleleA_ProbeSeq','AddressB_ID','AlleleB_ProbeSeq',
                'GenomeBuild','Chr','MapInfo','Ploidy','Species',
                'Source','SourceVersion','SourceStrand','SourceSeq',
                'TopGenomicSeq', 'CustomerStrand', 'GenomicStrand']

      yield header

      snp_count  = self.snp_count
      norm_ids   = self.norm_ids
      unpackL    = struct.Struct('<L').unpack
      unpackLL   = struct.Struct('<LL').unpack
      unpack3BLB = struct.Struct('<3BLB').unpack
      unpack4B4L = struct.Struct('<4B4L').unpack

      for i in xrange(snp_count):
        record_version,         = unpackL(read(4))

        if record_version not in (4,7,8):
          raise ValueError('Unsupported Illumina BPM record version: %d' % record_version)

        ilmnid                  = readstr(bpm)
        name                    = readstr(bpm)
        something1,something2,  \
        something3,snpnum,      \
        something4              = unpack3BLB(read(8))
        assert 1<=snpnum<=snp_count
        assert snpnum==snp_count-i
        design_strand           = intern(readstr(bpm))
        alleles                 = intern(readstr(bpm))
        chromosome              = intern(readstr(bpm))
        ploidy                  = intern(readstr(bpm))
        species                 = intern(readstr(bpm))
        mapinfo                 = readstr(bpm)
        top_genomic_sequence    = readstr(bpm)
        customer_strand         = intern(readstr(bpm))
        addressA_id,addressB_id = unpackLL(read(8))
        alleleA_probe_sequence  = readstr(bpm)
        alleleB_probe_sequence  = readstr(bpm)
        genome_version          = intern(readstr(bpm))
        source                  = intern(readstr(bpm))
        source_version          = intern(readstr(bpm))
        source_strand           = intern(readstr(bpm))
        source_sequence         = readstr(bpm)
        norm_id                 = norm_ids[snpnum-1]

        if record_version in (7,8):
          something5,something6,something7,assay_type_id, \
          something8,something9,something10,something11     = unpack4B4L(read(20))
        else:
          assay_type_id = 0

        genomic_strand = intern(readstr(bpm)) if record_version==8 else None

        norm_ids[snpnum-1] = norm_id = norm_id + 100*assay_type_id + 1

        if chromosome=='0':
          chromosome = ''

        yield ManifestRow(ilmnid,name,design_strand,alleles,assay_type_id,norm_id,
                          addressA_id,alleleA_probe_sequence,addressB_id,alleleB_probe_sequence,
                          genome_version,chromosome,mapinfo,ploidy,species,
                          source,source_version,source_strand,source_sequence,top_genomic_sequence,
                          customer_strand,genomic_strand)

        if 0:
          print 'SNP %d, SNPs left %d' % (i+1,snp_count-i-1)
          print 'record_version:',record_version
          print 'ilmnid:',ilmnid
          print 'name:',name
          print 'something1:',something1
          print 'something2:',something2
          print 'something3:',something3
          print 'snpnum:',snpnum
          print 'something4:',something4
          print 'norm_id:',self.norm_ids[snpnum-1]
          print 'snp_name:',self.snp_names[snpnum-1]
          print 'snp_entry:',self.snp_entries[snpnum-1]
          print 'design_strand:',design_strand
          print 'alleles:',alleles
          print 'chromosome:',chromosome
          print 'ploidy:',ploidy
          print 'species:',species
          print 'mapinfo:',mapinfo
          print 'top_genomic_sequence:',top_genomic_sequence
          print 'customer_strand:',customer_strand
          print 'addressA_id:',addressA_id
          print 'addressB_id:',addressB_id
          print 'alleleA_ProbeSeq:',alleleA_probe_sequence
          print 'alleleB_ProbeSeq:',alleleB_probe_sequence
          print 'genome_version:',genome_version
          print 'source:',source
          print 'source_version:',source_version
          print 'source_strand:',source_strand
          print 'source_sequence:',source_sequence
          print 'assay_type_id:',assay_type_id

          if record_version in (7,8):
            print 'something5:',something5
            print 'something6:',something6
            print 'something7:',something7
            print 'something8:',something8
            print 'something9:',something9
            print 'something10:',something10
            print 'something11:',something11

          if record_version==8:
            print 'genomic_strand:',printable(genomic_strand)

          print

      if 0: # Ignore stuff at the end of the file for now
        if filesize-bpm.tell():
          stuff, = unpackL(read(4))
          for i in range(stuff):
            line = readstr(bpm)
            print line

    if 0:
      print 'filename:',filename
      print 'version',self.version
      print 'manifest_name:',self.manifest_name
      print 'Controls:'
      for i,control in enumerate(self.controls):
        print '  %02d' % (i+1),printable(','.join(control))
      print 'snp_count:',self.snp_count
      print 'snp_entries:',self.snp_entries[:10],self.snp_entries[-10:]
      print 'snp_names:',self.snp_names[:10]+self.snp_names[-10:]
      print 'norm_ids:',self.norm_ids[:10],self.norm_ids[-10:]

    self.manifest_data = _manifest_rows()

  def __iter__(self):
    return self.manifest_data


def find_index(header,headings,optional=False):
  for h in headings:
    try:
      return header.index(h)
    except ValueError:
      pass

  if not optional:
    raise ValueError('Cannot find heading index')


def orient_manifest(manifest,targetstrand='customer',errorhandler=None):
  '''
  Parse an Illumina manifest to obtain mappings from A/B probes to alleles
  relative to a specified genomic strand.
  '''
  if errorhandler is None:
    def errorhandler(msg):
      raise ValueError(msg)

  gstrandmap   = {'+':'F','-':'R'}
  indelmap     = {'D':'-','I':'+','A':'A','C':'C','G':'G','T':'T'}
  strandmap    = {STRAND_UNKNOWN:STRAND_UNKNOWN,'+':'-','-':'+'}
  NA           = ('N','A')
  AB           = ('A','B')
  cnv          = ('-','+')
  targetstrand = targetstrand.lower()
  validstrand  = set(('top','bot','p','plus','m','minus'))
  plusminus    = set(('p','plus','m','minus'))

  assert targetstrand in ('ab','top','bottom','forward','reverse',
                          'real_forward','real_reverse',
                          'customer','anticustomer','design','antidesign')

  manifest    = iter(manifest)
  header      = manifest.next()
  assayid_idx = find_index(header,['IlmnID','Ilmn ID'])
  name_idx    = find_index(header,['Name'])
  alleles_idx = find_index(header,['SNP'])
  chrom_idx   = find_index(header,['Chr','CHR'])
  loc_idx     = find_index(header,['MapInfo'])
  cstrand_idx = find_index(header,['CustomerStrand','SourceStrand'])
  dstrand_idx = find_index(header,['IlmnStrand','Ilmn Strand'])
  gstrand_idx = find_index(header,['GenomicStrand'],optional=True)
  topseq_idx  = find_index(header,['TopGenomicSeq'])
  probea_idx  = find_index(header,['AlleleA_ProbeSeq'],optional=True)
  probeb_idx  = find_index(header,['AlleleB_ProbeSeq'],optional=True)

  for assay in manifest:
    lname   = assay[name_idx]
    alleles = assay[alleles_idx]
    chrom   = assay[chrom_idx] or None
    loc     = assay[loc_idx]
    cstrand = assay[cstrand_idx].lower()
    dstrand = assay[dstrand_idx].lower()
    topseq  = assay[topseq_idx]
    gstrand = assay[assayid_idx].split('_')[-2]

    if targetstrand in ('real_forward','real_reverse'):
      gstrand = gstrandmap.get(assay[gstrand_idx]) or 'U'

    if chrom and chrom.startswith('chr'):
      chrom = chrom[3:]
    if chrom=='Mt':
      chrom='M'

    if cstrand not in validstrand:
      errorhandler('Invalid customer strand %s for %s' % (cstrand,lname))
      yield lname,chrom,loc,STRAND_UNKNOWN,None
      continue

    if dstrand not in validstrand:
      errorhandler('Invalid design strand %s for %s' % (dstrand,lname))
      yield lname,chrom,loc,STRAND_UNKNOWN,None
      continue

    if gstrand not in '+-FRU':
      errorhandler('Unknown gstrand %s for %s' % (gstrand,lname))
      yield lname,chrom,loc,STRAND_UNKNOWN,None
      continue

    if len(alleles) != 5 or (alleles[0],alleles[2],alleles[4]) != ('[','/',']'):
      errorhandler('Invalid SNP alleles %s for %s' % (alleles,lname))
      yield lname,chrom,loc,STRAND_UNKNOWN,None
      continue

    a,b = alleles[1],alleles[3]

    if loc:
      loc = int(loc)

    if targetstrand=='ab' or a==b:
      yield lname,chrom,loc,STRAND_UNKNOWN,AB
      continue

    # CNV probes can simply be skipped
    if (a,b) == NA:
      yield lname,chrom,loc,STRAND_UNKNOWN,cnv
      continue

    # Indels are strand neutral
    elif cstrand in plusminus or dstrand in plusminus:
      a = indelmap[a]
      b = indelmap[b]
      yield lname,chrom,loc,None,(a,b)
      continue

    # SNP with A and B alleles on the design strand
    elif None not in (probea_idx,probeb_idx):
      probea = assay[probea_idx]
      probeb = assay[probeb_idx]
      if probea and probeb and (a,b) != (probea[-1],probeb[-1]):
        errorhandler('Design alleles do not match probes (%s,%s) != (%s,%s) for %s'
                         % (a,b,probea[-1],probeb[-1],lname))
        yield lname,chrom,loc,STRAND_UNKNOWN,None
        continue

    try:
      tstrand,aa,bb = norm_snp_seq(topseq)
      tstrand       = tstrand.lower()
      if tstrand != 'top':
        errorhandler('Supplied sequence is not correctly normalize to top strand for %s' % lname)
        yield lname,chrom,loc,STRAND_UNKNOWN,None
        continue

    except ValueError:
      tstrand,aa,bb = dstrand,a,b

    if dstrand!=tstrand:
      a,b = complement_base(a),complement_base(b)

    # a,b and aa,bb should both be normalized to tstrand and must be equal
    if (a,b) != (aa,bb):
      errorhandler('Assay alleles do not match sequence alleles (%s/%s != %s/%s) for %s' % (a,b,aa,bb,lname))
      yield lname,chrom,loc,STRAND_UNKNOWN,None
      continue

    if gstrand == '+':
      forward = True
      strand  = '+'
    elif gstrand == '-':
      forward = False
      strand  = '-'
    elif gstrand != 'U':
      # Get the strand orientation of the design sequence
      # Alleles are forward strand if the tstrand matches the design strand
      # and the design is on the forward strand or the converse of both
      # conditions is true.
      forward = (tstrand != dstrand) ^ (gstrand == 'F')
      strand  = '+' if forward else '-'
    else:
      if targetstrand in ('real_forward','real_reverse') and gstrand == 'U':
        raise ValueError("Unknown strand for assay '%s'" % lname)
      if targetstrand in ('forward','reverse') and gstrand == 'U':
        raise ValueError("Unknown strand for assay '%s'" % lname)
      strand = STRAND_UNKNOWN

    flip =    ((targetstrand == 'customer'     and tstrand != cstrand)
           or  (targetstrand == 'anticustomer' and tstrand == cstrand)
           or  (targetstrand == 'design'       and tstrand != dstrand)
           or  (targetstrand == 'antidesign'   and tstrand == dstrand)
           or  (targetstrand == 'top'          and tstrand != 'top'  )
           or  (targetstrand == 'bottom'       and tstrand != 'bot'  )
           or  (targetstrand == 'real_forward' and not forward       )
           or  (targetstrand == 'forward'      and not forward       )
           or  (targetstrand == 'real_reverse' and     forward       )
           or  (targetstrand == 'reverse'      and     forward       ))

    if flip:
      a,b = complement_base(a),complement_base(b)
      strand = strandmap[strand]

    yield lname,chrom,loc,strand,(a,b)


def create_abmap(manifest,genome):
  '''
  Create mapping from A/B probes to final alleles from the oriented manifest
  entries and updating the genome metadata for each locus.
  '''
  abmap = {}
  for lname,chrom,loc,strand,alleles in manifest:
    genome.merge_locus(lname, chromosome=chrom, location=loc, strand=strand)
    if alleles:
      abmap[lname] = alleles
  return abmap


def create_Illumina_abmap(filename,genome,targetstrand='customer',errorhandler=None):
  manifest = IlluminaManifest(filename)
  manifest = orient_manifest(manifest,targetstrand=targetstrand,errorhandler=errorhandler)
  return create_abmap(manifest,genome)


def manifest_snps(filename):
  '''
  Parse an Illumina manifest to obtain SNPs and locations.
  '''
  manifest    = IlluminaManifest(filename)
  manifest    = iter(manifest)
  header      = manifest.next()
  name_idx    = find_index(header,['Name'])
  chrom_idx   = find_index(header,['Chr'])
  loc_idx     = find_index(header,['MapInfo'])

  for assay in manifest:
    lname   = assay[name_idx]
    chrom   = assay[chrom_idx] or None
    loc     = assay[loc_idx]

    if chrom and chrom.startswith('chr'):
      chrom = chrom[3:]
    if chrom=='Mt':
      chrom='M'

    if loc:
      loc = int(loc)

    yield lname,chrom,loc


def read_Illumina_IDAT(filename):
  import numpy as np

  idat = autofile(filename,'rb')

  sig = idat.read(4)

  if sig != 'IDAT':
    raise ValueError('Invalid IDAT file signature')

  version, = struct.unpack('<L', idat.read(4))

  if version != 3:
    return dict(filename=filename,version=version)

  unknown0,field_count = struct.unpack('<LL', idat.read(8))

  fields = {}
  for i in xrange(field_count):
    field_code,field_offset = struct.unpack('<HQ', idat.read(10))
    field_name = IDAT_FIELD_CODES.get(field_code,'Unknown')
    if field_name in fields:
      raise ValueError('Invalid duplicated field %s in IDAT file' % field_name)
    fields[field_name] = field_offset

  if 0: # DEBUG TOC
    from operator import itemgetter

    s = os.stat(filename)
    filesize = s.st_size

    sfields = sorted(fields.items(),key=itemgetter(1,0))
    for i,(field_name,field_offset) in enumerate(sfields):
      if i+1 < len(sfields):
        field_size = sfields[i+1][1]-field_offset
      else:
        field_size = filesize - field_offset
      print '  ',field_name,field_offset,field_size

    print '  ','file size',filesize

  snp_count = illumina_ids = sds = means = bead_counts = midblock \
            = red_green = manifest = barcode = format = label     \
            = opa = sampleid = descr = plate = well = unknown = None
  runinfo   = []

  if 'snp_count' in fields:
    idat.seek(fields['snp_count'])
    snp_count, = struct.unpack('<L',idat.read(4))

  if 'illumina_ids' in fields:
    idat.seek(fields['illumina_ids'])
    illumina_ids = np.fromfile(idat, dtype='<u4', count=snp_count)

  if 'sds' in fields:
    idat.seek(fields['sds'])
    sds = np.fromfile(idat, dtype='<u2', count=snp_count)

  if 'means' in fields:
    idat.seek(fields['means'])
    means = np.fromfile(idat, dtype='<u2', count=snp_count)

  if 'bead_counts' in fields:
    idat.seek(fields['bead_counts'])
    bead_counts = np.fromfile(idat, dtype='b', count=snp_count)

  if 'midblock' in fields:
    idat.seek(fields['midblock'])
    midblock_entry_count, = struct.unpack('<L', idat.read(4))
    midblock = np.fromfile(idat, dtype='u4', count=midblock_entry_count)

  if 'red_green' in fields:
    idat.seek(fields['red_green'])
    red_green = struct.unpack('4B', idat.read(4))

  if 'manifest' in idat:
    idat.seek(fields['manifest'])
    manifest = readstr(idat)

  if 'barcode' in fields:
    idat.seek(fields['barcode'])
    barcode = readstr(idat)

  if 'format' in fields:
    idat.seek(fields['format'])
    format = readstr(idat)

  if 'label' in fields:
    idat.seek(fields['label'])
    label = readstr(idat)

  if 'opa' in fields:
    idat.seek(fields['opa'])
    opa = readstr(idat)

  if 'sampleid' in fields:
    idat.seek(fields['sampleid'])
    sampleid = readstr(idat)

  if 'descr' in fields:
    idat.seek(fields['descr'])
    descr = readstr(idat)

  if 'plate' in fields:
    idat.seek(fields['plate'])
    plate = readstr(idat)

  if 'well' in fields:
    idat.seek(fields['well'])
    well = readstr(idat)

  if 'unknown' in fields:
    idat.seek(fields['unknown'])
    unknown = readstr(idat)

  if 'runinfo' in fields:
    idat.seek(fields['runinfo'])

    runinfo_entry_count, = struct.unpack('<L', idat.read(4))

    for i in xrange(runinfo_entry_count):
     timestamp    = readstr(idat)
     entry_type   = readstr(idat)
     parameters   = readstr(idat)
     codeblock    = readstr(idat)
     code_version = readstr(idat)

     runinfo.append( (timestamp,entry_type,parameters,codeblock,code_version) )

  if 0: # DEBUG internal decoding
    print 'filename:',filename
    print 'idat_version:',version
    print 'unknown0:',unknown0
    print 'snp_count:',snp_count
    print 'Illumina IDs:',illumina_ids[:5]
    print 'SDs:',sds[:5]
    print 'Means:',means[:5]
    print 'Bead counts:',bead_counts[:5]
    print 'midblock:',midblock[:5]
    print 'RedGreen:',red_green
    print 'manifest:',manifest
    print 'barcode:',barcode
    print 'format:',format
    print 'label:',label
    print 'opa:',printable(opa)
    print 'sampleid:',sampleid
    print 'descr:',printable(descr)
    print 'plate:',plate
    print 'well:',well
    print 'unknown:',printable(unknown)

    for i,(timestamp,entry_type,parameters,codeblock,code_version) in enumerate(runinfo):
      print 'RunInfo %d:' % (i+1)
      print '  timestamp:',timestamp
      print '  entry_type:',entry_type
      print '  parameters:',parameters
      print '  codeblock:',codeblock
      print '  code_version:',code_version

    print

  return IDATData(filename,version,snp_count,illumina_ids,sds,means,bead_counts,
                  midblock,red_green,manifest,barcode,format,label,opa,sampleid,
                  descr,plate,well,runinfo)


def read_Illumina_LBD(filename,options):
  def sample_generator(sopts):
    import numpy as np

    include = sopts.include
    exclude = sopts.exclude

    for genos,scores in chunk(data,2):
      if genos[:6]!=scores[:6] or genos[6]!='calls' or scores[6]!='Score_Call':
        raise ValueError('Invalid Locus-by-DNA report found %s' % namefile(filename))

      sampleid = genos[0]

      if include is not None and sampleid not in include:
        continue

      if exclude is not None and sampleid     in exclude:
        continue

      genos = np.array(genos[8:],dtype='S1')
      try:
        scores = np.array(scores[8:],dtype=float)
      except ValueError:
        scores = np.array( [ float(s) if s!='NaN' else np.nan for s in scores[8:] ] )

      yield sampleid,genos,scores

  data = csv.reader(autofile(filename), dialect='csv')

  row = next(data)

  if row != ['OPA','LinkedGentrainFilePath']:
    raise ValueError('Invalid Locus-by-DNA report found %s' % namefile(filename))

  LBD_HEADER1 = ['oligoPoolId','recordType', 'data']
  LBD_HEADER2 = ['oligoPoolId','GTS LocusId','data']

  for row in data:
    if row == LBD_HEADER2:
      raise ValueError('Invalid Locus-by-DNA report found %s' % namefile(filename))
    elif row == LBD_HEADER1:
      break

  for row in data:
    if row == LBD_HEADER2:
      break
    elif row and row[1] == 'Gentrain Scores':
      gentrain = map(float,row[3:])
    elif row and row[1] == 'locusNames':
      loci = row[3:]

  for row in data:
    if row[:1] == ['instituteLabel']:
      break

  samples = sample_generator(options.samples)

  assert len(gentrain) == len(loci)

  return loci,gentrain,samples


def main():
  if 0:
    read_Illumina_IDAT('test/4666416106_R01C01_Grn.idat')

  if 0:
    read_Illumina_IDAT('test/4583987055_R02C01_Grn.idat')
    read_Illumina_IDAT('test/4583987055_R02C01_Red.idat')
    read_Illumina_IDAT('test/4196813065_A_Grn.idat')
    read_Illumina_IDAT('test/1495210050_A_Red.idat')

  if 0:
    import glob
    for filename in glob.iglob('/home/jacobske/projects/CGEMS/Scans/TGS/Old_Raw/GWAS_Raw_Data/DATA/AdvProstate/repeat/MEC/*.idat'):
      read_Illumina_IDAT(filename)

  if 0:
    import traceback
    for filename in open('idats'):
      try:
        read_Illumina_IDAT(filename.strip())
      except KeyboardInterrupt:
        break
      except:
        print 'ERROR: Cannot read %s' % filename
        print
        print 'Traceback:  %s\n' % (traceback.format_exc().replace('\n','\n  '))
        print

  if 0:
    manifest1 = IlluminaManifest('test/HumanHap300_(v1.0.0).bpm.gz')
    print 'filename:',manifest1.filename
    print 'SNPs:',manifest1.snp_count
    print 'rows:',len(list(manifest1))

  if 0:
    manifest2 = IlluminaManifest('test/HumanHap300_(v1.0.0).csv')
    print 'filename:',manifest2.filename
    print 'SNPs:',manifest2.snp_count
    print 'rows:',len(list(manifest2))

  if 0:
    list(IlluminaManifest('test/Human1M-Duov3_B.bpm'))
    list(IlluminaManifest('test/HumanHap300_(v1.0.0).bpm'))
    list(IlluminaManifest('test/HumanLinkage-12_E.bpm'))
    list(IlluminaManifest('test/HumanOmni1-Quad_v1-0_B_062509.bpm'))
    list(IlluminaManifest('test/HumanCNV370-Quadv3_C.bpm'))
    list(IlluminaManifest('test/CGEMS_P_F2_272225_A.bpm'))
    list(IlluminaManifest('test/Human610-Quadv1_B.bpm'))
    list(IlluminaManifest('test/Rare_Cancer_272049_A.bpm'))
    list(IlluminaManifest('test/240S_243991_Loci.bpm'))
    list(IlluminaManifest('test/Breast_Wide_Track_271628_A.bpm'))
    list(IlluminaManifest('test/CI1v410407_270414_A.bpm'))
    list(IlluminaManifest('test/NHL_B2B_271329_A.bpm'))
    list(IlluminaManifest('test/nsSNP_1x12_A_13371.bpm'))
    list(IlluminaManifest('test/Human660W-Quad_v1_A.bpm'))
    list(IlluminaManifest('test/HumanCNV370-Duo_v1-0_C.bpm'))


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

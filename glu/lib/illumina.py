# -*- coding: utf-8 -*-

__abstract__  = 'utility functions and objects to read and parse Illumina data files'
__copyright__ = 'Copyright (c) 2009, BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.'
__license__   = 'See GLU license for terms by running: glu license'
__revision__  = '$Id$'

import os
import csv
import struct
import traceback

import numpy as np

from   glu.lib.utils     import chunk
from   glu.lib.fileutils import autofile,parse_augmented_filename,guess_format,get_arg
from   glu.lib.sections  import read_sections


IDAT_FIELD_CODES = { 1000 : 'nSNPsRead',
                      102 : 'IlluminaID',
                      103 : 'SD',
                      104 : 'Mean',
                      107 : 'NBeads',
                      200 : 'MidBlock',
                      300 : 'RunInfo',
                      400 : 'RedGreen',
                      401 : 'Manifest',
                      402 : 'Barcode',
                      403 : 'ChipType',
                      404 : 'Stripe',
                      405 : 'Unknown1',
                      406 : 'SampleID',
                      407 : 'Unknown2',
                      408 : 'Plate',
                      409 : 'Well',
                      510 : 'Unknown3' }


def printable(s):
  return (''.join(map(myord,s or ''))).replace('\n','(\\n)')


def myord(x):
  import string
  if x in string.printable[:-5]:
    return x
  else:
    return '(%02d)' % ord(x)


def readstr(afile):
  '''
  String data are encoded as a sequence of one or more length bytes followed
  by the specified number of data bytes.

  If the high-bit of the first length byte is set, then a second length byte
  follows with the number of additional 128 character blocks.  This
  acommodates strings up to length 16,384 (128**2) without ambiguity.  It is
  unknown of this scheme scales to additional length bytes, since no strings
  longer than 6,264 bytes have been observed in the wild.
  '''
  read = afile.read
  n = ord(read(1))

  if not n:
    return ''

  if n&0x80:
    m = ord(read(1))
    if m:
      n += (m-1)<<7

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
    self.snp_types   = None

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

    self.manifest_data = contents


  def _load_bpm(self,filename):
    self.filename = filename
    stats         = os.stat(filename)
    filesize      = stats.st_size
    bpm           = autofile(filename,'rb')
    read          = bpm.read
    readline      = bpm.readline

    signature     = read(3)

    if signature != 'BPM':
      raise ValueError('Invalid BPM file signature: %s' % printable(sig))

    self.version = struct.unpack('<BL', read(5))

    if self.version != (1,4):
      raise ValueError('Invalid BPM version number (%d.%d)' % version)

    self.manifest_name = readstr(bpm)
    self.controls      = [ c.split(',') for c in readstr(bpm).split('\n') if c ]
    self.snp_count,    = struct.unpack('<L', read(4))
    self.snp_entries   = np.fromfile(bpm, dtype='<u4', count=self.snp_count)
    self.snp_names     = [ readstr(bpm) for i in xrange(self.snp_count) ]
    self.snp_types     = np.fromfile(bpm, dtype='u1', count=self.snp_count)

    def _manifest_rows():
      header = ['IlmnID','Name','IlmnStrand','SNP',
                'AddressA_ID','AlleleA_ProbeSeq','AddressB_ID','AlleleB_ProbeSeq',
                'GenomeBuild','Chr','MapInfo','Ploidy','Species',
                'Source','SourceVersion','SourceStrand','SourceSeq',
                'TopGenomicSeq']

      yield header

      snp_count  = self.snp_count
      unpackL    = struct.Struct('<L').unpack
      unpackLL   = struct.Struct('<LL').unpack
      unpack3BLB = struct.Struct('<3BLB').unpack

      for i in xrange(snp_count):
        record_version_maybe,   = unpackL(read(4))
        assert record_version_maybe in (4,7)
        ilmnid                  = readstr(bpm)
        name                    = readstr(bpm)
        something1,something2,  \
        something3,snpnum,      \
        something4              = unpack3BLB(read(8))
        assert 1<=snpnum<=snp_count
        assert snpnum==snp_count-i
        design_strand           = readstr(bpm)
        alleles                 = readstr(bpm)
        chromosome              = readstr(bpm)
        ploidy                  = readstr(bpm)
        species                 = readstr(bpm)
        mapinfo                 = readstr(bpm)
        top_genomic_sequence    = readstr(bpm)
        customer_strand         = readstr(bpm)
        addressA_id,addressB_id = unpackLL(read(8))
        alleleA_probe_sequence  = readstr(bpm)
        alleleB_probe_sequence  = readstr(bpm)
        genome_version          = readstr(bpm)
        source                  = readstr(bpm)
        source_version          = readstr(bpm)
        source_strand           = readstr(bpm)
        source_sequence         = readstr(bpm)

        if record_version_maybe == 7:
          endstuff = read(20)

        yield [ilmnid,name,design_strand,alleles,
               addressA_id,alleleA_probe_sequence,addressB_id,alleleB_probe_sequence,
               genome_version,chromosome,mapinfo,ploidy,species,
               source,source_version,source_strand,source_sequence,top_genomic_sequence]

        if 0:
          print 'SNP %d, SNPs left %d' % (i+1,snp_count-i-1)
          print 'record_version_maybe:',record_version_maybe
          print 'ilmnid:',ilmnid
          print 'name:',name
          print 'something1:',something1
          print 'something2:',something1
          print 'something3:',something1
          print 'snpnum:',snpnum
          print 'something4:',something1
          print 'snp_type:',self.snp_types[snpnum-1]
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
          print 'AlleleA_ProbeSeq',alleleA_probe_sequence
          print 'AlleleB_ProbeSeq',alleleB_probe_sequence
          print 'genome_version:',genome_version
          print 'source:',source
          print 'source_version:',source_version
          print 'source_strand:',source_strand
          print 'source_sequence:',source_sequence

          if record_version_maybe == 7:
            print 'version7_endstuff',printable(endstuff)

          print

      if 0: # Ignore stuff at the end of the file for now
        if filesize-bpm.tell():
          stuff, = unpackL(read(4))
          for i in range(stuff):
            line = readstr(bpm)
            print line

    if 0:
      print 'filename:',filename
      print 'version',version
      print 'manifest_name:',manifest_name
      print 'Controls:'
      for i,control in enumerate(controls):
        print '  %02d' % (i+1),printable(','.join(control))
      print 'snp_count:',snp_count
      print 'snp_entries:',snp_entries[:10],snp_entries[-10:]
      print 'snp_names:',snp_names[:10]+snp_names[-10:]
      print 'snp_types:',snp_types[:10],snp_types[-10:]

    self.manifest_data = _manifest_rows()

  def __iter__(self):
    return self.manifest_data


def read_Illumina_IDAT(filename):
  s = os.stat(filename)

  filesize = s.st_size

  idat = autofile(filename,'rb')

  sig = idat.read(4)

  if sig != 'IDAT':
    raise ValueError('Invalid IDAT file signature')

  version, = struct.unpack('<L', idat.read(4))

  if version != 3:
    raise ValueError('Invalid IDAT version number: %d' % version)

  unknown0,field_count = struct.unpack('<LL', idat.read(8))

  fields = {}
  for i in xrange(field_count):
    field_code,field_offset = struct.unpack('<HQ', idat.read(10))
    field_name = IDAT_FIELD_CODES.get(field_code,'Unknown')
    if field_name in fields:
      raise ValueError('Invalid duplicated field %s in IDAT file' % field_name)
    fields[field_name] = field_offset

  if 0:
    from operator import itemgetter
    sfields = sorted(fields.items(),key=itemgetter(1,0))
    for i,(field_name,field_offset) in enumerate(sfields):
      if i+1 < len(sfields):
        field_size = sfields[i+1][1]-field_offset
      else:
        field_size = filesize - field_offset
      print '  ',field_name,field_offset,field_size

    print '  ','file size',filesize

  snp_count = illumina_ids = sds = means = bead_counts = midblock \
            = red_green = manifest = barcode = chip_type = stripe \
            = unknown1 = sampleid = unknown2 = plate = well       \
            = unknown3 = None
  runinfo   = []

  if 'nSNPsRead' in fields:
    idat.seek(fields['nSNPsRead'])
    snp_count, = struct.unpack('<L',idat.read(4))

  if 'IlluminaID' in fields:
    idat.seek(fields['IlluminaID'])
    illumina_ids = np.fromfile(idat, dtype='<u4', count=snp_count)

  if 'SD' in fields:
    idat.seek(fields['SD'])
    sds = np.fromfile(idat, dtype='<u2', count=snp_count)

  if 'Mean' in fields:
    idat.seek(fields['Mean'])
    means = np.fromfile(idat, dtype='<u2', count=snp_count)

  if 'NBeads' in fields:
    idat.seek(fields['NBeads'])
    bead_counts = np.fromfile(idat, dtype='b', count=snp_count)

  if 'MidBlock' in fields:
    idat.seek(fields['MidBlock'])
    midblock_entry_count, = struct.unpack('<L', idat.read(4))
    midblock = np.fromfile(idat, dtype='u4', count=midblock_entry_count)

  if 'RedGreen' in fields:
    idat.seek(fields['RedGreen'])
    red_green = struct.unpack('4B', idat.read(4))

  if 'Manifest' in idat:
    idat.seek(fields['Manifest'])
    manifest = readstr(idat)

  if 'Barcode' in fields:
    idat.seek(fields['Barcode'])
    barcode = readstr(idat)

  if 'ChipType' in fields:
    idat.seek(fields['ChipType'])
    chip_type = readstr(idat)

  if 'Stripe' in fields:
    idat.seek(fields['Stripe'])
    stripe = readstr(idat)

  if 'Unknown1' in fields:
    idat.seek(fields['Unknown1'])
    unknown1 = readstr(idat)

  if 'SampleID' in fields:
    idat.seek(fields['SampleID'])
    sampleid = readstr(idat)

  if 'Unknown2' in fields:
    idat.seek(fields['Unknown2'])
    unknown2 = readstr(idat)

  if 'Plate' in fields:
    idat.seek(fields['Plate'])
    plate = readstr(idat)

  if 'Well' in fields:
    idat.seek(fields['Well'])
    well = readstr(idat)

  if 'Unknown3' in fields:
    idat.seek(fields['Unknown3'])
    unknown3 = readstr(idat)

  # Who cares about RunInfo anyway?
  if 'RunInfo' in fields:
    idat.seek(fields['RunInfo'])

    runinfo_entry_count, = struct.unpack('<L', idat.read(4))

    for i in xrange(runinfo_entry_count):
     timestamp    = readstr(idat)
     entry_type   = readstr(idat)
     parameters   = readstr(idat)
     codeblock    = readstr(idat)
     code_version = readstr(idat)

     runinfo.append( (timestamp,entry_type,parameters,codeblock,code_version) )

  if 0:
    print 'filename:',filename
    print 'version:',version
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
    print 'chip_type:',chip_type
    print 'stripe:',stripe
    print 'unknown1:',printable(unknown1)
    print 'sampleid:',sampleid
    print 'unknown2:',printable(unknown2)
    print 'plate:',plate
    print 'well:',well
    print 'unknown3:',printable(unknown3)

    for i,(timestamp,entry_type,parameters,codeblock,code_version) in enumerate(runinfo):
      print 'RunInfo %d:' % (i+1)
      print '  timestamp:',timestamp
      print '  entry_type:',entry_type
      print '  parameters:',parameters
      print '  codeblock:',codeblock
      print '  code_version:',code_version

    print


def read_Illumina_LBD(filename,options):
  def parse_gentrain():
    row = data.next()
    assert row[1] == 'Gentrain Scores'
    assert row[2] == ''
    return map(float,row[3:])

  def parse_loci():
    row = data.next()
    assert row[1] == 'locusNames'
    assert row[2] == ''
    return row[3:]

  def skiprows(n):
    for i in xrange(n):
      data.next()

  def sample_generator(sopts):
    for genos,scores in chunk(data,2):
      assert genos[:6] == scores[:6]
      assert  genos[6] == 'calls'
      assert scores[6] == 'Score_Call'

      sampleid = genos[0]

      if sopts.include is not None and sampleid not in sopts.include:
        continue

      if sopts.exclude is not None and sampleid     in sopts.exclude:
        continue

      genos = np.array(genos[8:],dtype='S1')
      try:
        scores = np.array(scores[8:],dtype=float)
      except ValueError:
        scores = np.array( [ float(s) if s!='NaN' else np.nan for s in scores[8:] ] )

      yield sampleid,genos,scores

  data = csv.reader(autofile(filename), dialect='csv')

  skiprows(11)
  gentrain = parse_gentrain()
  skiprows(3)
  loci    = parse_loci()
  skiprows(5)
  samples = sample_generator(options.samples)

  assert len(gentrain) == len(loci)

  return loci,gentrain,samples


def main():
  if 0:
    read_Illumina_IDAT('test/4583987055_R02C01_Grn.idat')
    read_Illumina_IDAT('test/4583987055_R02C01_Red.idat')
    read_Illumina_IDAT('test/4196813065_A_Grn.idat')
    #read_Illumina_IDAT('test/1495210050_A_Red.idat')

    import glob
    for filename in glob.iglob('/mnt/nfs/gigantor/vol1/GWAS/Scans/Prostate/3/builds/3/delivery/ATBC/genotypes/rawdata/*/*.idat'):
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

  if 1:
    manifest2 = IlluminaManifest('test/HumanHap300_(v1.0.0).csv')
    print 'filename:',manifest2.filename
    print 'SNPs:',manifest2.snp_count
    print 'rows:',len(list(manifest2))

  if 0:
    IlluminaManifest('test/Human1M-Duov3_B.bpm')
    IlluminaManifest('test/HumanHap300_(v1.0.0).bpm')
    IlluminaManifest('test/HumanLinkage-12_E.bpm')
    IlluminaManifest('test/HumanOmni1-Quad_v1-0_B_062509.bpm')
    IlluminaManifest('test/HumanCNV370-Quadv3_C.bpm')
    IlluminaManifest('test/CGEMS_P_F2_272225_A.bpm')
    IlluminaManifest('test/Human610-Quadv1_B.bpm')
    IlluminaManifest('test/Rare_Cancer_272049_A.bpm')
    IlluminaManifest('test/240S_243991_Loci.bpm')
    IlluminaManifest('test/Breast_Wide_Track_271628_A.bpm')
    IlluminaManifest('test/CI1v410407_270414_A.bpm')
    IlluminaManifest('test/NHL_B2B_271329_A.bpm')
    IlluminaManifest('test/nsSNP_1x12_A_13371.bpm')
    IlluminaManifest('test/Human660W-Quad_v1_A.bpm')
    IlluminaManifest('test/HumanCNV370-Duo_v1-0_C.bpm')


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

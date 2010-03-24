# -*- coding: utf-8 -*-

__abstract__  = 'utility functions and objects to read and parse Illumina data files'
__copyright__ = 'Copyright (c) 2009, BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.'
__license__   = 'See GLU license for terms by running: glu license'
__revision__  = '$Id$'

import os
import csv
import struct

from   glu.lib.utils     import chunk, namedtuple
from   glu.lib.fileutils import autofile,parse_augmented_filename,guess_format,get_arg
from   glu.lib.sections  import read_sections


# Credit for a large proportion of the reverse-engineering in this module
# goes to the R/BioConductor project, who in turn were helped by Keith
# Baggerly.


ManifestRow = namedtuple('ManifestRow',
                         'ilmnid name design_strand alleles assay_type_id norm_id '
                         'addressA_id alleleA_probe_sequence addressB_id alleleB_probe_sequence '
                         'genome_version chromosome mapinfo ploidy species '
                         'source source_version source_strand source_sequence top_genomic_sequence customer_strand')

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


def readstr(afile):
  '''
  String data are encoded as a sequence of one or more length bytes followed
  by the specified number of data bytes.

  If the high-bit of the first length byte is set, then a second length byte
  follows with the number of additional 128 character blocks.  This
  accommodates strings up to length 16,384 (128**2) without ambiguity.  It is
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

    self.manifest_data = contents

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
                'TopGenomicSeq', 'CustomerStrand']

      yield header

      snp_count  = self.snp_count
      unpack4B   = struct.Struct('<4B').unpack
      unpackL    = struct.Struct('<L').unpack
      unpackLL   = struct.Struct('<LL').unpack
      unpack4L   = struct.Struct('<4L').unpack
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
        norm_id                 = self.norm_ids[snpnum-1]

        if record_version_maybe == 7:
          something5,something6,something7,assay_type_id = unpack4B(read(4))
          something8,something9,something10,something11  = unpack4L(read(16))
        else:
          assay_type_id = 0

        self.norm_ids[snpnum-1] = norm_id = norm_id + 100*assay_type_id + 1

        yield ManifestRow(ilmnid,name,design_strand,alleles,assay_type_id,norm_id,
                          addressA_id,alleleA_probe_sequence,addressB_id,alleleB_probe_sequence,
                          genome_version,chromosome,mapinfo,ploidy,species,
                          source,source_version,source_strand,source_sequence,top_genomic_sequence,
                          customer_strand)

        if 0:
          print 'SNP %d, SNPs left %d' % (i+1,snp_count-i-1)
          print 'record_version_maybe:',record_version_maybe
          print 'ilmnid:',ilmnid
          print 'name:',name
          print 'something1:',something1
          print 'something2:',something2
          print 'something3:',something3
          print 'snpnum:',snpnum
          print 'something4:',something4
          print 'snp_type:',self.norm_ids[snpnum-1]
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

          if record_version_maybe == 7:
            print 'something5:',something5
            print 'something6:',something6
            print 'something7:',something7
            print 'something8:',something8
            print 'something9:',something9
            print 'something10:',something10
            print 'something11:',something11

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


def read_Illumina_IDAT(filename):
  import numpy as np

  s = os.stat(filename)

  filesize = s.st_size

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

  if 0:
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
    import numpy as np

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

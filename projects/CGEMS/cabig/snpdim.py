import csv
import sys

from biozilla.fileutils import autofile
from itertools          import count


DATA_HEADER_ILLUMINA = ['IlmnID','Name','IlmnStrand','SNP','AddressA_ID','AlleleA_ProbeSeq','AddressB_ID',
                        'AlleleB_ProbeSeq','Chr','MapInfo','Ploidy','Species','CustomerStrand','IlmnStrand',
                        'IllumicodeSeq','TopGenomicSeq']

DATA_HEADER_AFFY     = ['Probe Set ID','Affy SNP ID','dbSNP RS ID','Chromosome','Genome Version','DB SNP Version',
                        'Physical Position','Strand','ChrX pseudo-autosomal region','Cytoband','Flank','Allele A',
                        'Allele B','Associated Gene','Genetic Map','Microsatellite','Fragment Length Start Stop',
                        'Freq Asian','Freq AfAm','Freq Cauc','Het Asian','Het AfAm','Het Cauc','Num chrm Asian',
                        'Num chrm AfAm','Num chrm Cauc']

DIMOUT_HEADER        = ['CHROMOSOME','DBSNP_BUILD','DBSNPID','GENOME_BUILD',
                        'PHYSICAL_LOCATION','REFERENCE_SEQUENCE','REFERENCE_STRAND','SNPANNO_ID']

COMPLEMENT_ALLELE    = {'A':'T', 'T':'A', 'G':'C', 'C':'G', 'Y':'R', 'R':'Y', 'S':'S', 'W':'W','K':'M','M':'K',
                        'B':'V', 'V':'B', 'D':'H', 'H':'D','N':'N', '[':']', ']':'[', '/':'/'}

def reverse_complement(seq):
  rcseq = ''
  for a in seq:
    rcseq = COMPLEMENT_ALLELE[a] + rcseq
  return rcseq


def snp_dim_ilm(manifestfile,dimid,snpmap):
  manifestfile = csv.reader(autofile(manifestfile))
  for row in manifestfile:
    if 'IlmnID' in row:
      header = row
      break
  assert header == DATA_HEADER_ILLUMINA

  for row in manifestfile:
    if '[Controls]' in row:
      raise StopIteration

    snpid          = row[1]
    dbsnp_build    = '124'
    genome_build   = '35'
    snp_unique_str = snpid + '_' + dbsnp_build + '_' + genome_build

    if snp_unique_str in snpmap:
      snpannoid = snpmap[snp_unique_str]
    else:
      chr = row[8].upper()
      if 'chr' not in chr and 'M' not in chr:
        chr = 'chr' + chr

      physical_loc = row[9]
      ref_strand   = row[12]
      if ref_strand == 'Bot':
        ref_seq = reverse_complement(row[15].upper())
      elif ref_strand == 'Top':
        ref_seq = row[15].upper()
      snpannoid = dimid.next()
      snpmap[snp_unique_str] = snpannoid

      yield [chr,dbsnp_build,snpid,genome_build,physical_loc,ref_seq,ref_strand,snpannoid]


def snp_dim_aff(manifestfile,dimid,snpmap):

  def _check_id(ids):
    '''
    Create a mapping file between dbsnp ids and snpannoids;
    if unknown such as '---' used by Affy, use other unique ids for identification purpose only;
    If others, raise error.
    '''
    snpid = ids[1]
    if snpid[:2].upper() in ('RS','MI'):
      rowid        = snpid
    elif snpid == '---':
      rowid        = ids[0]
    else:
      raise ValueError, "Invalid dbsnp id in '%s'" % snpid
    return rowid


  manifestfile = csv.reader(autofile(manifestfile))
  for row in manifestfile:
    if 'Probe Set ID' in row:
      header = row
      break
  assert header == DATA_HEADER_AFFY

  for row in manifestfile:
    snpid = row[2]
    rowid = _check_id(row[1:3])

    dbsnp_build    = row[5].split(',')[0]
    genome_build   = '35'
    snp_unique_str = rowid + '_' + dbsnp_build + '_' + genome_build

    if snp_unique_str in snpmap:
      snpannoid = snpmap[snp_unique_str]
    else:
      chr = row[3].upper()
      if 'chr' not in chr and 'M' not in chr:
        chr = 'chr' + chr

      physical_loc = row[6]
      ref_seq      = row[10].upper()
      ref_strand   = row[7]
      snpannoid    = dimid.next()
      snpmap[snp_unique_str] = snpannoid

      yield [chr,dbsnp_build,snpid,genome_build,physical_loc,ref_seq,ref_strand,snpannoid]


def option_parser():
  import optparse

  usage = 'usage: %prog [options] file'
  parser = optparse.OptionParser(usage=usage)

  parser.add_option('-o', '--outfile',   dest='outfile',   metavar='FILE',   default= '-',
                    help='Output file for table snp_dim')
  parser.add_option('-m', '--snpmap',    dest='snpmap',    metavar='FILE',   default= '-',
                    help='Output file mapping dbsnp ids')
  parser.add_option('-c', '--inicount',  dest='inicount',  type = 'int')

  return parser



def main():
  parser = option_parser()
  options,args = parser.parse_args()

  dimoutfile = csv.writer(autofile(options.outfile,'w'))
  dimoutfile.writerow(DIMOUT_HEADER)
  snpmapoutfile  = csv.writer(autofile(options.snpmap,'w'))

  dimid = count(options.inicount)
  snpmap = {}

  datatype = 1
  for arg in args:
    if arg == 'ilm':
      datatype = 1
    elif arg == 'aff':
      datatype = 0
    else:
      if datatype:
        for dim in snp_dim_ilm(arg,dimid,snpmap):
          if dim:
            dimoutfile.writerow(dim)
      else:
        for dim in snp_dim_aff(arg,dimid,snpmap):
          if dim:
            dimoutfile.writerow(dim)

  for key,val in snpmap.iteritems():
    snpmapoutfile.writerow([key.split('_')[0],val,key])



if __name__=='__main__':
  main()


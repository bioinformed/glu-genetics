# -*- coding: utf-8 -*-
'''
Input genodata file format(command line argument)

  sdat
        sdat	l1	l2
        s1	g1	g2

  ldat
        ldat	s1	s2
        l1	g1	g2

  genotriples
        s1	l1	g1


Optional input file to indicate the regions

  regions (-e/--regions)
  [region]
  region1
  [samples]
  s1
  s2
  [loci]
  l1
  l2
  [region]
  region2
  [samples]
  s1
  s2
  s3
  [loci]
  l2
  l3
  l4

Optionsl input file to map the sample/locus to grouping variables

  samplegroup (-g/--samplegroup)
        Samples	Group
        s1	name1

  locusgroup (-G/--locusgroup)
        Samples	Group
        l1	name1
'''

__gluindex__  = False
__abstract__  = 'Compute assay and sample completion'
__copyright__ = 'Copyright (c) 2007-2009, BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.'
__license__   = 'See GLU license for terms by running: glu license'
__revision__  = '$Id$'

import sys
import unittest

from   textwrap             import fill
from   collections          import defaultdict

from   glu.lib.utils        import percent
from   glu.lib.fileutils    import map_reader,autofile,hyphen
from   glu.lib.genolib      import load_genostream
from   glu.lib.regionparser import load_regions
from   glu.lib.sections     import save_section, SectionWriter, save_metadata_section


def output_summary(out,samcomp,samempty,droppedsam,droppedloc,loccomp,locempty,nonmissing,genos_inf,genos_all):
  '''
  Output completion summary report
  '''

  totalsamples = set(samcomp) | droppedsam
  totalloci    = set(loccomp) | droppedloc
  infsamples   = set(samcomp) - samempty
  infloci      = set(loccomp) - locempty
  out.write('\nSamples     Total: %7d, Empty: %7d, Dropped: %7d, Informative: %7d\n' %
           (len(totalsamples),len(samempty),len(droppedsam),len(infsamples)))
  out.write('Loci        Total: %7d, Empty: %7d, Dropped: %7d, Informative: %7d\n' %
           (len(totalloci),len(locempty),len(droppedloc),len(infloci)))

  out.write('\nGLOBAL GENOTYPE COMPLETION RATE FOR ALL DATA:         %10d / %10d = %5.3f%%\n' % \
           (nonmissing,genos_all,percent(nonmissing,genos_all)))
  out.write('GLOBAL GENOTYPE COMPLETION RATE FOR INFORMATIVE DATA: %10d / %10d = %5.3f%%\n' % \
           (nonmissing,genos_inf,percent(nonmissing,genos_inf)))
  e = ' '*5
  empty = e+'(empty)'

  out.write('\nSamples with no data:\n')
  out.write(fill(', '.join(samempty), initial_indent=e,subsequent_indent=e) or empty)
  out.write('\n')
  out.write('\nLoci with no data:\n')
  out.write(fill(', '.join(locempty), initial_indent=e,subsequent_indent=e) or empty)
  out.write('\n')

  out.write('\nDropped samples:\n')
  out.write(fill(', '.join(droppedsam), initial_indent=e,subsequent_indent=e) or empty)
  out.write('\n')
  out.write('\nDropped loci:\n')
  out.write(fill(', '.join(droppedloc), initial_indent=e,subsequent_indent=e) or empty)
  out.write('\n')


def output_detail(out,name,comp):
  '''
  Output completion detail report by sample/locus or by group
  '''
  out.write('\n\nMISSING GENOTYPES BY %s\n' % name.upper())
  out.write('                                                         Informative                   All\n')
  out.write('  Rank     %-25s  Empty   Dropped     N / Total      %%           N / Total      %%  \n' % name)
  out.write('  -------  -------------------------  ------- ------- -----------------  ------  -----------------  ------\n')

  data = sorted( (vals[0],k,vals[2],vals[3],vals[4],vals[5]) for k,vals in comp.iteritems() )

  for rank,(nonmiss,id,empty,dropped,inf,total) in enumerate(data):
    d = (rank+1,id,empty,dropped,nonmiss,inf,percent(nonmiss,inf),nonmiss,total,percent(nonmiss,total))
    out.write('  %7d  %-25s  %7d %7d %7d / %7d  %6.2f  %7d / %7d  %6.2f\n' % d)


def output_group(out,name,comp):
  '''
  Output completion detail report by sample/locus or by group
  '''
  out.write('\n\nMISSING GENOTYPES BY %s\n' % name.upper())
  out.write('                                                         Informative                   All\n')
  out.write('  Rank     %-17s   Members Empty   Dropped    N / Total      %%           N / Total      %%  \n' % name)
  out.write('  -----  ------------------  ------- ------- ------- -------------------  ------  -------------------  ------\n')

  data = sorted( (vals[0],k,vals[2],vals[3],vals[4],vals[5],vals[6]) for k,vals in comp.iteritems() )

  for rank,(nonmiss,gname,empty,dropped,inf,total,mtotal) in enumerate(data):
    d = (rank+1,gname,mtotal,empty,dropped,nonmiss,inf,percent(nonmiss,inf),nonmiss,total,percent(nonmiss,total))
    out.write('  %5d  %-17s  %7d %7d %7d %8d /%9d  %6.2f  %8d /%9d  %6.2f\n' % d)


def save_elements(sw, elements, name, type):
  '''
  Writes a section for empty or dropped samples/loci to a file
  '''
  data = [['type', type]] + [[e] for e in elements]
  save_section(sw, name, data)


def save_summary(sw,comp,emptyset,dropped,name):
  '''
  Writes a section for completion summary to a file
  '''
  inf   = set(comp) - emptyset
  total = set(comp) | dropped

  data = [['type',        name],
          ['total',       len(total)],
          ['empty',       len(emptyset)],
          ['dropped',     len(dropped)],
          ['complete',    len(inf)]]

  save_section(sw, 'summary', data)


def save_group(sw, comp, name):
  '''
  Writes a section for completion statistics by specified groups to a file
  '''
  data = sorted( (vals[0],k,vals[2],vals[3],vals[4],vals[5],vals[6]) for k,vals in comp.iteritems() )
  data = [(gname,mtotal,empty,dropped,nonmiss,total) for nonmiss,gname,empty,dropped,inf,total,mtotal in data]
  data = [['data',name],['id','members','empty','dropped','completed','total']] + data
  save_section(sw, 'group', data)


def save_detail(sw, comp, name):
  '''
  Writes a section for completion statistics by sample/locus to a file
  '''
  data = sorted( (vals[0],k,vals[2],vals[3],vals[4],vals[5]) for k,vals in comp.iteritems() )
  data = [(id,1,empty,dropped,nonmiss,total) for nonmiss,id,empty,dropped,inf,total in data]
  data = [['type', name],['id','members', 'empty','dropped','completed','total']] + data
  save_section(sw, 'data', data)


def count_missing(genotriples):
  '''
  Count missing genotype by locus/sample

  @param genotriples: genotype triplets
  @type  genotriples: tuple

  >>> genotriples = ('1420_11', 'rs2070', 17),('1420_12', 'rs2070', 19),('1420_2', 'rs2070', 0)
  >>> map(dict,count_missing(genotriples))
  [{'1420_2': [0, 1], '1420_11': [1, 0], '1420_12': [1, 0]}, {'rs2070': [2, 1]}]
  '''
  samcomp = defaultdict(lambda: [0,0])
  loccomp = defaultdict(lambda: [0,0])

  for sample,locus,geno in genotriples:
    missing = not geno
    samcomp[sample][missing] += 1
    loccomp[locus][missing]  += 1
  return samcomp,loccomp


def process_expected(regions,loccomp,samcomp,samempty,locempty):
  '''
  Populate the statistics with additional information
  '''
  for locus,counts in loccomp.iteritems():
    if regions:
      expected = expected_samples(locus,regions)
      stats    = calculate(samempty,expected,samcomp,counts)
    else:
      inf      = len(samcomp)-len(samempty) if locus not in locempty else 0
      stats    = [len(samempty),0,inf,len(samcomp),len(samcomp)-sum(counts)]
    loccomp[locus].extend(stats)

  for sample,counts in samcomp.iteritems():
    if regions:
      expected = expected_loci(sample,regions)
      stats    = calculate(locempty,expected,loccomp,counts)
    else:
      inf      = len(loccomp)-len(locempty) if sample not in samempty else 0
      stats    = [len(locempty),0,inf,len(loccomp),len(loccomp)-sum(counts)]
    samcomp[sample].extend(stats)


def calculate(emptyset,expected,comp,counts):
  '''
  Compute the additional statistics
  '''
  empty   = len(emptyset & expected)
  dropped = len(expected - set(comp))
  inf     = sum(counts)  - empty
  total   = len(expected)
  return  [empty,dropped,inf,total]


def expected_samples(locus, regions):
  '''
  Compute the expected sample set per locus

  @param   locus: locus id
  @type    locus: str
  @param regions: list of genotyping regions
  @type  regions: list
  @return       : set of sample ids
  @rtype        : set
  '''
  samples = set()
  for name,rsample,rlocus in regions:
    if locus in rlocus:
      samples.update(rsample)
  return samples


def expected_loci(sample, regions):
  '''
  Compute the expected locus set per sample

  @param  sample: sample id
  @type   sample: str
  @param regions: list of genotyping regions
  @type  regions: list
  @return       : a set of locus ids
  @rtype        : set
  '''
  loci = set()
  for name,rsample,rlocus in regions:
    if sample in rsample:
      loci.update(rlocus)
  return loci


def dropped_all(samcomp,loccomp,regions):
  '''
  Compute dropped samples/loci
  '''
  if regions:
    expectedsam_all = set()
    expectedloc_all = set()
    for name,rsample,rlocus in regions:
      expectedsam_all.update(rsample)
      expectedloc_all.update(rlocus)

    droppedsam = expectedsam_all-set(samcomp)
    droppedloc = expectedloc_all-set(loccomp)
    return droppedsam,droppedloc
  else:
    return set(),set()


def build_emptyset(counts):
  '''
  Build a set of samples/loci without any genotypes
  '''
  return set(k for k,vals in counts.iteritems() if vals[0]==0)


def completion(genotriples,regions):
  '''
  Compute the completion statistics

  @param  genotriples: genotype triplets
  @type   genotriples: tuple
  @param      regions: list of genotyping regions
  @type       regions: list
  '''
  print >> sys.stderr, '\nComputing completion rates...'

  samcomp,loccomp = count_missing(genotriples)
  samempty  = build_emptyset(samcomp)
  locempty  = build_emptyset(loccomp)

  process_expected(regions,loccomp,samcomp,samempty,locempty)

  # FIXME: Not thrilled by all of the opaque positional references
  nonmissing = sum(vals[0] for vals in samcomp.itervalues())
  genos_inf  = sum(vals[4] for vals in samcomp.itervalues())
  genos_all  = sum(vals[5] for vals in samcomp.itervalues())

  return samcomp,samempty,loccomp,locempty,nonmissing,genos_inf,genos_all


def completion_group(comp,groupmap):
  '''
  Compute the completion statistics to report completion rates by specified groups

  @param      comp: the counts of missing and non-missing genotypes by sample/locus
  @type       comp: dict
  @param  groupmap: map sample/locus ids to grouping variables
  @type   groupmap: dict
  '''
  groupcomp = defaultdict(lambda: [0,0,0,0,0,0,0])

  for i,g in comp.iteritems():
    group = groupmap.get(i) or 'default'
    # FIXME: Not thrilled by all of the opaque positional references
    comp  = groupcomp[group]
    for i in range(6):
      comp[i] += g[i]
    comp[6] += 1

  return groupcomp


def option_parser():
  import optparse

  usage = 'usage: %prog [options] file'
  parser = optparse.OptionParser(usage=usage)

  parser.add_option('-f', '--format',          dest='format',          metavar='NAME',
                    help='Format of the input data. Values=sdat,ldat,hapmap,genotriple')
  parser.add_option('-g', '--genorepr',        dest='genorepr',        metavar='REP',
                    help='Input genotype representation')
  parser.add_option('-l', '--loci', dest='loci', metavar='FILE',
                    help='Locus description file and options')
  parser.add_option('-o', '--output',          dest='output',          metavar='FILE', default='-',
                    help='Output of completion report')
  parser.add_option('-e', '--regions', dest='regions', metavar='FILE',
                    help='Regions of genotypes expected to be genotyped. Used to compute overall completion.')
  parser.add_option(     '--samplegroup',     dest='samplegroup',     metavar='FILE',
                    help='Map the sample ids to the grouping variable')
  parser.add_option(     '--locusgroup',      dest='locusgroup',      metavar='FILE',
                    help='Map the locus ids to the grouping variable')
  parser.add_option(      '--tabularoutput',   dest='tabularoutput',   metavar='FILE',
                    help='Generate machine readable tabular output of results')

  return parser


class TestCompletion(unittest.TestCase):
  def setUp(self):
    self.regions = [('region1', set(['1420_6','1510_2','1510_3','1420_12']), set(['rs2070','rs619'])),
                    ('region2', set(['1420_6','1510_2','1510_3','1420_8','1420_9','1420_12']), set(['rs6048','rs7100']))]
    self.loci = ['rs619','rs2070','rs6048']
    self.samples = ['1420_12','1420_6','1420_8','1420_9']
    self.samcounts = {'1420_12':[5,1],'1420_6':[5,1],'1420_8':[0,4],'1420_9':[3,4]}

  def test_expected_samples(self):
    results = [set(['1420_6', '1420_12', '1510_2', '1510_3']),
               set(['1420_6', '1420_12', '1510_2', '1510_3']),
               set(['1420_6', '1510_2', '1510_3', '1420_8', '1420_9', '1420_12'])]
    expected = []
    for i,locus in enumerate(self.loci):
      expectedsample = expected_samples(locus,self.regions)
      self.assertEquals(expectedsample, results[i])
      expected.append(expectedsample)
    return expected


  def test_expected_loci(self):
    results = [set(['rs2070', 'rs6048', 'rs619', 'rs7100']),
               set(['rs2070', 'rs6048', 'rs619', 'rs7100']),
               set(['rs6048', 'rs7100']),
               set(['rs6048', 'rs7100'])]
    expected = []
    for i,sample in enumerate(self.samples):
      expectedlocus = expected_loci(sample,self.regions)
      self.assertEquals(expectedlocus, results[i])
      expected.append(expectedlocus)
    return expected



def main():
  parser = option_parser()
  options,args = parser.parse_args()

  if len(args) != 1:
    parser.print_help(sys.stderr)
    sys.exit(2)

  regions = None
  if options.regions:
    regions = list(load_regions(options.regions))

  outfile     = autofile(hyphen(options.output, sys.stdout),'w')
  genostream  = load_genostream(args[0],format=options.format,genorepr=options.genorepr,
                                       genome=options.loci,hyphen=sys.stdin)
  genotriples = genostream.as_genotriples()

  samcomp,samempty,loccomp,locempty,nonmissing,genos_inf,genos_all = completion(genotriples,regions)

  print >> sys.stderr, '\nWriting completion output...'

  droppedsam,droppedloc = dropped_all(samcomp,loccomp,regions)
  genos_all += len(droppedsam)*(len(loccomp)+len(droppedloc))
  output_summary(outfile,samcomp,samempty,droppedsam,droppedloc,loccomp,locempty,nonmissing,genos_inf,genos_all)

  output_detail(outfile,'Samples', samcomp)
  output_detail(outfile,'Loci', loccomp)

  if options.samplegroup:
    samplegroup = map_reader(options.samplegroup)
    groupcomp   = completion_group(samcomp,samplegroup)
    output_group(outfile,'Sample Group', groupcomp)

  if options.locusgroup:
    locusgroup = map_reader(options.locusgroup)
    groupcomp  = completion_group(loccomp,locusgroup)
    output_group(outfile,'Locus Group', groupcomp)

  if options.tabularoutput:
    sw = SectionWriter(options.tabularoutput)
    save_metadata_section(sw, analysis='completion', analysis_version='0.1', format_version='0.1')
    globalcomp    = set(samcomp)    | set(loccomp)
    globalempty   = set(samempty)   | set(locempty)
    globaldropped = set(droppedsam) | set(droppedloc)
    save_summary(sw,globalcomp,globalempty,globaldropped,'global')
    save_summary(sw,samcomp,samempty,droppedsam,'samples')
    save_summary(sw,loccomp,locempty,droppedloc,'loci')

    if options.samplegroup:
      samplegroup = map_reader(options.samplegroup)
      groupcomp = completion_group(samcomp,samplegroup)
      save_group(sw, groupcomp, 'samples')

    if options.locusgroup:
      locusgroup = map_reader(options.locusgroup)
      groupcomp = completion_group(loccomp,locusgroup)
      save_group(sw, groupcomp, 'loci')

    save_detail(sw, samcomp, 'samples')
    save_detail(sw, loccomp, 'loci')

    save_elements(sw, locempty,  'empty', 'loci')
    save_elements(sw, samempty,  'empty', 'samples')
    save_elements(sw, droppedsam, 'dropped', 'samples')
    save_elements(sw, droppedloc, 'dropped', 'loci')

  print >> sys.stderr, 'Done.\n'


def _test():
  import doctest
  return doctest.testmod()


if __name__ == '__main__':
  #unittest.main(argv=[sys.argv[0],'-q'])
  _test()
  main()

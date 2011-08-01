# -*- coding: utf-8 -*-

__gluindex__  = True
__program__   = 'HagZilla'
__authors__   = ['Kevin Jacobs (jacobs@bioinformed.com)']
__abstract__  = 'Simplified interface for tag-based assay design'
__copyright__ = 'Copyright (c) 2007-2009, BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.'
__license__   = 'See GLU license for terms by running: glu license'
__revision__  = '$Id$'

import os
import sys
import time
import select
import fcntl
import subprocess

from   itertools               import groupby
from   operator                import itemgetter

from   glu.lib.fileutils       import list_reader,table_reader
from   glu.modules.ld.tagzilla import locus_result_sequence
from   glu.lib.genedb          import open_genedb
from   glu.lib.genedb.queries  import query_snps_by_location, query_gene_by_name


#FIXME: Make even moe configurable
HAPPATH='/usr/local/share/hapmap/'

DATA = { 'hapmap28' : dict(   genome='genedb_ncbi36.3_dbsnp130',
                           genotypes=HAPPATH+'build28/glu/hapmap_%(POP)s_%(CHR)s_r28_fwd_nr_b36.lbat'),
         'hapmap27' : dict(   genome='genedb_ncbi36.3_dbsnp130',
                           genotypes=HAPPATH+'build27/glu/hapmap_%(POP)s_%(CHR)s_r27_fwd_nr_b36.lbat'),
         'hapmap26' : dict(   genome='genedb_ncbi36.3_dbsnp130',
                           genotypes=HAPPATH+'build26/forward/non-redundant/hapmap_%(POP)s_%(CHR)s_r26_fwd_nr_b36.lbat'),
         'hapmap23' : dict(   genome='genedb_ncbi36.3_dbsnp130',
                            pedigree=HAPPATH+'hapmap_pedigree.tsv',
                           genotypes=HAPPATH+'build23/rs_strand/non-redundant/genotypes_chr%(CHR)s_%(POP)s_r23a_nr.b36.txt.gz',
                              format='hapmap'),
         'hapmap22' : dict(   genome='genedb_ncbi36.3_dbsnp130',
                            pedigree=HAPPATH+'hapmap_pedigree.tsv',
                           genotypes=HAPPATH+'build22/rs_strand/non-redundant/genotypes_chr%(CHR)s_%(POP)s_r22_nr.b36.txt.gz',
                              format='hapmap'),
         'hapmap21' : dict(   genome='genome35-1.db',
                            pedigree=HAPPATH+'hapmap_pedigree.tsv',
                           genotypes=HAPPATH+'build21a/rs_strand/non-redundant/genotypes_chr%(CHR)s_%(POP)s_r21a_nr.txt.gz',
                              format='hapmap'),
         'hapmap19' : dict(   genome='genome35-1.db',
                            pedigree=HAPPATH+'hapmap_pedigree.tsv',
                           genotypes=HAPPATH+'build19/non-redundant/genotypes_chr%(CHR)s_%(POP)s.txt.gz',
                              format='hapmap') }


def project_chdir(*dirs):
  cwd = os.getcwd()

  for dir in dirs:
    try:
      os.mkdir(dir)
    except (OSError,IOError):
      pass
    #print 'DIR',dir
    os.chdir(dir)

  pwd = os.getcwd()
  return pwd


def project_path(*dirs):
  cwd = os.getcwd()
  try:
    return project_chdir(*dirs)
  finally:
    os.chdir(cwd)


def project_file(*dirs_and_file):
  dirs  = dirs_and_file[:-1]
  pfile = dirs_and_file[-1]
  return os.path.join(project_path(*dirs),pfile)


def get_genotypes(config,populations,chromosome):
  command = []

  populations = populations.split(',')
  command.append(' -M %s' % escape_spaces(','.join(populations)))

  for population in populations:
    if 'format' in config:
      command.append('-f %s' % escape_spaces(config['format']))
    if 'pedigree' in config:
      command.append('-p %s' % escape_spaces(config['pedigree']))
    if 'genotypes' in config:
      command.append( escape_spaces(config['genotypes'] % dict(POP=population,CHR=chromosome)) )

  return ' '.join(command)


def escape_spaces(s):
  return s.replace(' ','\\ ')


def run_tagzilla(config,outdir,project,gene,dprime,r2,populations,chromosome,snps,maf,include,exclude,
                         designfiles,designdefault,includefile,cmdprefix):

  prefix = '%s_%s_%s' % (project,gene,populations)

  if designfiles:
    designfiles =  [ os.path.abspath(d) for d in designfiles ]

  if includefile:
    includefile = os.path.abspath(includefile)

  cwd = os.getcwd()

  try:
    project_chdir(outdir, project, gene)

    command = []

    if cmdprefix:
      command.append(cmdprefix)

    command.append('glu ld.tagzilla')

    command.append('--filternonfounders')

    command.append('-s subset')
    file('subset','w').write('\n'.join(snps))

    if dprime:
      command.append('-d %s' % dprime)

    if r2:
      command.append('-r %s' % r2)

    if maf:
      command.append('-a %s' % maf)

    if include and includefile:
      command.append('-i include')
      incs = (set(list_reader(includefile)) | set(include)) & set(snps)
      file('include','w').write('\n'.join(sorted(incs)))
    elif includefile:
      command.append('-i %s' % escape_spaces(includefile))
    elif include:
      command.append('-i include')
      file('include','w').write('\n'.join(include))

    if exclude:
      command.append('-e exclude')
      file('exclude','w').write('\n'.join(exclude))

    for d in designfiles or []:
      command.append('-D %s' % escape_spaces(d))

    command.append('--designdefault=%f' % designdefault)

    command.append(get_genotypes(config,populations,chromosome))

    command += ['-C maxsnp','-O loci.out','-b sum.txt','-B bins.txt','2>&1']
    command  = ' '.join(command)
    #print command

    return subprocess.Popen(command, cwd=os.getcwd(), shell=True, stdout=subprocess.PIPE, bufsize=1)
  finally:
    os.chdir(cwd)


def get_tags(outdir,project,gene,include,exclude):
  tags     = set()
  snps     = set()
  obligate = set()
  exclset  = set()

  bins = list(locus_result_sequence(project_file(outdir,project,gene,'loci.out'),{},exclude))

  for pop,bin in bins:
    if bin.disposition == 'maximal-bin':
      snps.update(bin)
      tags.update(bin.recommended_tags)
    elif bin.disposition in ('obligate-untyped','obligate-typed'):
      obligate.update(bin.recommended_tags)
    elif bin.disposition == 'obligate-exclude':
      exclset.update(bin.recommended_tags)

  if include:
    obligate |= set(include)

  if exclude:
    obligate -= set(exclude)

  return snps,tags,obligate,exclset


def extend(s,n):
  m = len(s)
  if m < n:
    s = list(s)
    s.extend(['']*(n-m))
  return s


def strip(r):
  return [ f.strip() for f in r ]


def gene_margin(start, stop, strand, upstream=20000, downstream=10000):
  if None in (start,stop,strand):
    return None,None

  if start>stop:
    start,stop = stop,start

  if strand == '+':
    return start-upstream,stop+downstream
  elif strand == '-':
    return start-downstream,stop+upstream
  else:
    raise ValueError('Unknown gene orientation')


def clean_tasks(con,tasks,options):
  for task in tasks:
    project,population,gene,chromosome,start,stop,strand,dprime,r2,maf,include,exclude = task

    if not project:
      task[0] = 'default'

    if chromosome and start and stop:
      task[4] = int(start)
      task[5] = int(stop)
    else:
      try:
        alias,qgene,qchromosome,qstart,qstop,qstrand,qtype = query_gene_by_name(con, gene)
      except KeyError,e:
        sys.stderr.write('Invalid gene: %s\n' % e.args[0])
        continue

      if chromosome and chromosome != qchromosome:
        sys.stderr.write('Warning: Chromosome does not match for %s (%s != %s)\n' %
                                     (gene,chromosome,qchromosome))

      if None in (qstart,qstop):
        sys.stderr.write('Gene not mapped: %s\n' % gene)
        continue

      start,stop = gene_margin(qstart,qstop,qstrand,upstream=options.upstream,
                                                  downstream=options.downstream)

      task[3] = qchromosome
      task[4] = start
      task[5] = stop
      task[6] = qstrand

    task[10] = set(i.strip() for i in include.split(',') if i.strip())
    task[11] = set(e.strip() for e in exclude.split(',') if e.strip())

    snps = set(s[0] for s in query_snps_by_location(con, task[3], task[4], task[5]))
    task.append(snps)

    yield task


def setblocking(fd):
  '''
  set blocking mode
  '''
  flags = fcntl.fcntl(fd, fcntl.F_GETFL)|os.O_NONBLOCK
  fcntl.fcntl(fd, fcntl.F_SETFL, flags)


class SubprocessManager(object):
  '''
  Simple subprocess manager

  This object manages a set of concurrently running sub-processes, as
  returned by the subprocess module, and discards all output.  This allows
  all subprocesses to run until completion without blocking on output.

  The current implementation probably works on Linux only.
  '''
  def __init__(self):
    self.poll    = select.poll()
    self.taskmap = {}

  def add(self,proc):
    fd = proc.stdout.fileno()
    setblocking(fd)
    assert fd not in self.taskmap
    self.taskmap[fd] = proc
    self.poll.register(fd,select.POLLIN|select.POLLHUP)

  def run(self,timeout=None):
    done = []
    for fd,event in self.poll.poll(timeout):
      proc = self.taskmap[fd]
      if event&select.POLLIN:
        proc.stdout.read(8196)
      if event&select.POLLHUP:
        proc.wait()
        done.append(proc)
        del self.taskmap[fd]
        self.poll.unregister(fd)
    return done

  def __len__(self):
    return len(self.taskmap)


def option_parser():
  from glu.lib.glu_argparse import GLUArgumentParser

  parser = GLUArgumentParser(description=__abstract__)

  parser.add_argument('dataset', help='Name of reference data set')
  parser.add_argument('outdir',  help='Output directory name')
  parser.add_argument('jobfile', help='Tabular or delimited job definition file')

  parser.add_argument('-D', '--designscores', metavar='FILE', action='append',
                    help='Design scores (optional)')
  parser.add_argument('--designdefault', metavar='N', type=float, default=0,
                    help='Default design score for any locus not found in a design file')
  parser.add_argument('-i', '--include',   metavar='FILE',
                      help='Obligate tags (optional)')
  parser.add_argument('--upstream',   default=20000, type=int,  metavar='N',
                    help='upstream margin in bases')
  parser.add_argument('--downstream', default=10000, type=int,  metavar='N',
                    help='downstream margin in bases')
  parser.add_argument('--dryrun', action='store_true',
                          help='Perform a dryrun. Perform all processing except actually tagging')
  parser.add_argument('--maxqueue', default=1, type=int,
                          help='Maximum number of jobs to queue in parallel')
  parser.add_argument('--cmdprefix', metavar='VALUE',
                    help='Command prefix for running tagzilla (useful for scheduling jobs)')

  return parser


def main():
  parser  = option_parser()
  options = parser.parse_args()

  if options.dataset not in DATA:
    sys.stderr.write('[ERROR] Invalid dataset specified: %s\n' % options.dataset)
    sys.stderr.write('[ERROR] Choices are: %s\n' % ', '.join(sorted(DATA)))
    sys.exit(2)

  config = DATA[options.dataset]

  con = open_genedb(config.get('genome'))

  if not os.path.isdir(options.outdir):
    sys.stderr.write('[ERROR] Output directory is not a directory: %s\n' % options.outdir)
    sys.exit(1)

  taskfile   = table_reader(options.jobfile)
  header     = taskfile.next()
  tasks      = sorted( extend(strip(t),12) for t in taskfile )
  tasks      = sorted(clean_tasks(con,tasks,options))

  genecounts = {}
  genefail   = {}

  print '\t'.join(['Project','Population','Region','Chr','Start','Stop', 'Strand',
                   'SNPs', 'Qualifying SNPs', 'Tags', 'Ob. Include', 'Ob. Exclude'])

  manager=SubprocessManager()

  startdelay = 0.25

  for task in tasks:
    # Run queue until we have space for the next job
    while len(manager) >= options.maxqueue:
      manager.run()

    project,population,gene,chromosome,start,stop,strand,dprime,r2,maf,include,exclude,snps = task

    if options.dryrun:
      print '\t'.join(map(str,[project,population,gene,chromosome,start,stop,strand]))
      continue

    if manager:
      time.sleep(startdelay)

    proc=run_tagzilla(config,options.outdir,project,gene,dprime,r2,population,chromosome,snps,
                      maf,include,exclude,options.designscores,options.designdefault,
                      options.include,options.cmdprefix)
    manager.add(proc)


  if options.dryrun:
    return

  # Drain queue
  while manager:
    manager.run()

  # Reap results
  for project,ptasks in groupby(tasks,itemgetter(0)):
    ptags = set()
    for task in ptasks:
      project,population,gene,chromosome,start,stop,strand,dprime,r2,maf,include,exclude,snps = task
      qsnps,tags,obligate,excl = get_tags(options.outdir,project,gene,include,exclude)

      design = tags|obligate
      ptags.update(tags|obligate)

      designfile = file(project_file(options.outdir,project,gene,'design.lst'),'w')
      for rs in sorted(design):
        designfile.write('%s\n' % rs)

      genecounts[gene] = len(design)

      d = map(str,[project,population,gene,chromosome,start,stop,strand,
                   len(snps),len(qsnps),len(tags),len(obligate),len(excl)])
      print '\t'.join(d)

    if options.dryrun:
      continue

    designfile = file(project_file(options.outdir,project,'design.lst'),'w')
    for rs in ptags:
      designfile.write('%s\n' % rs)

    pdir = project_path(options.outdir,project)
    command = 'glu ld.binsum "%s"/*/loci.out > "%s"/sum.out 2>/dev/null' % (pdir,pdir)
    binsum = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, bufsize=1)
    sys.stderr.write(binsum.stdout.read())

    sumfile = file(project_file(options.outdir,project,'genesum.out'),'w')
    for gene,n in sorted(genecounts.iteritems()):
      sumfile.write('%s\t%d\n' % (gene,n))


if __name__ == '__main__':
  main()

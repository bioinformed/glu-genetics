# -*- coding: utf-8 -*-
'''
File:          hagzilla.py

Authors:       Kevin Jacobs (jacobs@bioinformed.com)

Created:       Wed May 17 08:51:57 EDT 2006

Abstract:      Front-end for running TagZilla for SNPPlex and Illumina assay design

Compatibility: Python 2.5 and above

Requires:

Revision:      $Id$
'''

__program__   = 'HagZilla'
__authors__   = ['Kevin Jacobs (jacobs@bioinformed.com)']
__copyright__ = 'Copyright (c) 2008, BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.'
__license__   = 'See GLU license for terms by running: glu license'

import os
import sys
import time
import select
import fcntl
import subprocess
import sqlite3

from   itertools                      import groupby,chain
from   operator                       import itemgetter

from   glu.lib.fileutils              import load_list,load_table
from   glu.modules.tagzilla.tagzilla  import launcher,locus_result_sequence


POPS = {'CEU'    : 'hapmap',
        'YRI'    : 'hapmap',
        'JPT'    : 'hapmap',
        'CHB'    : 'hapmap',
        'JPT+CHB': 'hapmap',
        'NHS'    : 'NHS',
       }

GENOMEDB='/usr/local/share/genedb/genome36-1.db'

#FIXME: Make configurable
if 1:
  HAP_PEDIGREES = '/usr/local/share/hapmap/peds'
  HAP_GENOTYPES = '/usr/local/share/hapmap/build22/rs_strand/non-redundant/genotypes_chr%s_%s_r22_nr.b36.txt.gz'
  HAP_FORMAT    = 'hapmap'
elif 0:
  HAP_PEDIGREES = '/usr/local/share/hapmap/peds'
  HAP_GENOTYPES = '/usr/local/share/hapmap/build21a/rs_strand/non-redundant/genotypes_chr%s_%s_r21a_nr.txt.gz'
  HAP_FORMAT    = 'hapmap'
else:
  HAP_PEDIGREES = '/usr/local/share/hapmap/peds'
  HAP_GENOTYPES = '/usr/local/share/hapmap/build19/non-redundant/genotypes_chr%s_%s.txt.gz'
  HAP_FORMAT    = 'hapmap'

NHS_PEDIGREES = None
NHS_GENOTYPES = '/u3/projects/CGEMS/Scans/Breast/1/current/genotypes/STUDY/subjects_STUDY_%s.ldat.gz'
NHS_MAP       = '/u3/projects/CGEMS/Scans/Breast/1/raw/snp.map'
NHS_FORMAT    = 'ldat'

PLCO_PEDIGREES = None
PLCO_GENOTYPES = '/u3/projects/CGEMS/Scans/Prostate/1/build2/genotypes/STUDY/subjects_STUDY_%s.ldat.gz'
PLCO_MAP       = '/u3/projects/CGEMS/Scans/Prostate/1/raw/snp.map'
PLCO_FORMAT    = 'ldat'


def query_genes_by_name(con,gene):
  sql = '''
  SELECT   a.Alias,s.featureName,s.chromosome,s.chrStart,s.chrEnd,s.orientation
  FROM     alias a, gene s
  WHERE    s.geneID = a.geneID
    AND    a.Alias %s
    AND    s.featureType='GENE'
    AND    s.submitGroup='reference'
  ORDER BY s.chromosome,MIN(s.chrStart,s.chrEnd);
  '''
  if '%' in gene:
    sql = sql % 'LIKE ?'
  else:
    sql = sql % '= ?'

  cur = con.cursor()
  cur.execute(sql, (gene,))
  return cur.fetchall()


def query_gene_by_name(con,gene):
  genes = query_genes_by_name(con,gene)
  if not genes:
    raise KeyError('Cannot find gene "%s"' % gene)
  elif len(genes) > 1:
    raise KeyError('Gene not unique "%s"' % gene)
  return genes[0]


def query_genes_by_location(con,chr,loc):
  sql = '''
  SELECT   s.featureName,s.featureName,s.chromosome,s.chrStart,s.chrEnd,s.orientation
  FROM     gene s
  WHERE    s.chromosome = %s
    AND    %d BETWEEN s.chrStart AND s.chrEnd
    AND    s.featureType='GENE'
    AND    s.submitGroup='reference'
  ORDER BY s.chromosome,MIN(s.chrStart,s.chrEnd);
  '''
  cur = con.cursor()
  cur.execute(sql, chr, loc)
  return cur.fetchall()


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


def get_snps(con, population, chromosome, start, stop):
  sql = '''
  SELECT   lname, chromosome, location
  FROM     snp
  WHERE    chromosome=?
    AND    location BETWEEN ? AND ?;
  '''
  cur = con.cursor()
  cur.execute(sql, (chromosome,start,stop) )
  return cur.fetchall()


def escape(s):
  return "'%s'" % s.replace("'","''")


def get_sequences(con,snps):
  cur = con.cursor()

  try:
    cur.execute('CREATE TABLE snpsequence (lname TEXT PRIMARY KEY, sequence TEXT);')
  except:
    pass

  sql = '''
  SELECT  lname,sequence
  FROM    snpsequence
  WHERE   lname IN (%s);
  '''
  sqlsnps = ','.join(escape(snp) for snp in snps if snp)
  cur.execute(sql % sqlsnps)
  results  = cur.fetchall()

  missing  = snps - set(lname for lname,seq in results)
  results2 = [ s for s in get_sequences_from_genewindow(missing) if len(s) == 2 ]

  missing -= set(lname for lname,seq in results2)
  results3 = ((lname,'') for lname in missing)

  sql = 'INSERT INTO snpsequence VALUES (?,?);'
  cur.executemany(sql, chain(results2,results3))

  con.commit()

  return ( (lname,seq) for lname,seq in chain(results,results2) if seq)


def get_sequences_from_genewindow(snps):
  snps = list(snps)
  if not snps:
    return []

  command = 'java -Xmx1000M -classpath sequence/gw_tools.jar:sequence/classes12.jar:sequence nci/cgf/annotator/tools/export/BatchExport'
  snpname = 'tmpFoo%d' % time.time()
  sequencename = 'tmpBar%d' % time.time()
  try:
    snpfile = file(snpname,'w')
    for snp in snps:
      snpfile.write('%s\n' % snp)
    snpfile.close()
    args = '%s %s sequence/sublist.txt sequence/oblig.txt %s 2>/dev/null' % (command,snpname,sequencename)
    subprocess.Popen(args, shell=True).communicate()
    return load_table(sequencename)
  finally:
    os.unlink(snpname)
    try:
      os.unlink(sequencename)
    except:
      pass


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


def get_genotypes(populations,chromosome):
  command = []

  populations = populations.split(',')
  command.append(' -M %s' % escape(','.join(populations)))

  if len(populations) > 1:
    command.append('--multimethod global')

  for population in populations:
    data = POPS[population]

    if data == 'hapmap':
      command.append('-f hapmap')
      command.append('-p %s' % escape(HAP_PEDIGREES))
      command.append(escape(HAP_GENOTYPES % (chromosome,population)))
    elif data == 'NHS':
      command.append('-f ldat')
      command.append('-l %s' % escape(NHS_MAP))
      command.append(NHS_GENOTYPES % chromosome)
    else:
      raise ValueError('Unknown population')

  return ' '.join(command)


def escape(s):
  return s.replace(' ','\\ ')


def run_tagzilla(outdir,project,gene,dprime,r2,populations,chromosome,snps,maf,include,exclude,
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

    command.append('python /usr/local/bin/glu tagzilla')

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
      incs = (set(load_list(includefile)) | set(include)) & set(snps)
      file('include','w').write('\n'.join(sorted(incs)))
    elif includefile:
      command.append('-i %s' % escape(includefile))
    elif include:
      command.append('-i include')
      file('include','w').write('\n'.join(include))

    if exclude:
      command.append('-e exclude')
      file('exclude','w').write('\n'.join(exclude))

    for d in designfiles or []:
      command.append('-D %s' % escape(d))

    command.append('--designdefault=%f' % designdefault)

    command.append(get_genotypes(populations,chromosome))

    command += ['-C maxsnp','-O loci.out','-b sum.txt','-B bins.txt','2>&1']
    command  = ' '.join(command)
    #print command

    return subprocess.Popen(command, cwd=os.getcwd(), shell=True, stdout=subprocess.PIPE)
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
        alias,qgene,qchromosome,qstart,qstop,qstrand = query_gene_by_name(con, gene)
      except KeyError,e:
        print >> sys.stderr, 'Invalid gene: %s' % e.args[0]
        continue

      if chromosome and chromosome != qchromosome:
        print >> sys.stderr, 'Warning: Chromosome does not match for %s (%s != %s)' % \
                                     (gene,chromosome,qchromosome)

      if None in (qstart,qstop):
        print >> sys.stderr, 'Gene not mapped: %s' % gene
        continue

      start,stop = gene_margin(qstart,qstop,qstrand,upstream=options.upstream,
                                                  downstream=options.downstream)

      task[3] = qchromosome
      task[4] = start
      task[5] = stop
      task[6] = qstrand

    task[10] = set(i.strip() for i in include.split(',') if i.strip())
    task[11] = set(e.strip() for e in exclude.split(',') if e.strip())

    snps = set(snp[0] for snp in get_snps(con, population, task[3], task[4], task[5]))
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
  import optparse

  usage = 'usage: %prog [options] outdir jobfile'
  parser = optparse.OptionParser(usage=usage)

  parser.add_option('-D', '--designscores', dest='designscores', metavar='FILE', action='append',
                    help='Design scores (optional)')
  parser.add_option('--designdefault', dest='designdefault', metavar='N', type='float', default=0,
                    help='Default design score for any locus not found in a design file')
  parser.add_option('-g', '--genomedb',   dest='genomedb', metavar='FILE', default=GENOMEDB,
                      help='Genome database file')
  parser.add_option('-i', '--include',   dest='include', metavar='FILE',
                      help='Obligate tags (optional)')
  parser.add_option('--upstream',   dest='upstream',   default=20000, type='int',  metavar='N',
                    help='upstream margin in bases')
  parser.add_option('--downstream', dest='downstream', default=10000, type='int',  metavar='N',
                    help='downstream margin in bases')
  parser.add_option('--getsequence', dest='getsequence', action='store_true',
                          help='Output design sequence from GeneWindow')
  parser.add_option('--dryrun', dest='dryrun', action='store_true',
                          help='Perform a dryrun. Perform all processing except actually tagging')
  parser.add_option('--maxqueue', dest='maxqueue', default=1, type='int',
                          help='Maximum number of jobs to queue in parallel')
  parser.add_option('--cmdprefix', dest='cmdprefix', metavar='VALUE',
                    help='Command prefix for running tagzilla (useful for scheduling jobs)')

  return parser


def main():
  parser = option_parser()
  options,args = parser.parse_args()

  if len(args) != 2:
    parser.print_help()
    return

  con = sqlite3.connect(options.genomedb)
  outdir = args[0]

  if not os.path.isdir(outdir):
    sys.stderr.write('[ERROR] Output directory is not a directory: %s\n' % outdir)

  taskfile   = load_table(args[1])
  header     = taskfile.next()
  tasks      = sorted( extend(strip(t),12) for t in taskfile )
  tasks      = sorted(clean_tasks(con,tasks,options))

  genecounts = {}
  genestuff  = {}
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

    proc=run_tagzilla(outdir,project,gene,dprime,r2,population,chromosome,snps,
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
      qsnps,tags,obligate,excl = get_tags(outdir,project,gene,include,exclude)

      design = tags|obligate
      ptags.update(tags|obligate)

      designfile = file(project_file(outdir,project,gene,'design.lst'),'w')
      for rs in sorted(design):
        designfile.write('%s\n' % rs)

      seqs = []
      if options.getsequence:
        designfile = file(project_file(outdir,project,gene,'design.fasta'),'w')
        for rs,seq in get_sequences(con, design):
          designfile.write('>%s\n%s\n' % (rs,seq))

      genecounts[gene] = len(design)
      genestuff[gene]  = seqs

      d = map(str,[project,population,gene,chromosome,start,stop,strand,
                   len(snps),len(qsnps),len(tags),len(obligate),len(excl)])
      print '\t'.join(d)

    if options.dryrun:
      continue

    designfile = file(project_file(outdir,project,'design.lst'),'w')
    for rs in ptags:
      designfile.write('%s\n' % rs)

    if options.getsequence:
      designfile = file(project_file(outdir,project,'design.fasta'),'w')
      for rs,seq in get_sequences(con, ptags):
        designfile.write('>%s\n%s\n' % (rs,seq))

    pdir = project_path(outdir,project)
    command = 'glu tagzilla.binsum "%s"/*/loci.out > "%s"/sum.out 2>/dev/null' % (pdir,pdir)
    print >> sys.stderr, subprocess.Popen(command, shell=True, stdout=subprocess.PIPE).stdout.read()

    sumfile = file(project_file(outdir,project,'genesum.out'),'w')
    for gene,n in sorted(genecounts.iteritems()):
      sumfile.write('%s\t%d\n' % (gene,n))

    sumfile = file(project_file(outdir,project,'genestuff.out'),'w')
    for gene,seqs in sorted(genestuff.iteritems()):
      for rs,seq in seqs:
        sumfile.write('%s\t%s\t%s\n' % (gene,rs,seq))


if __name__ == '__main__':
  main()

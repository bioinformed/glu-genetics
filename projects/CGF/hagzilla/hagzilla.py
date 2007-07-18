# -*- coding: utf-8 -*-
'''
File:          hagzilla.py

Authors:       Kevin Jacobs (jacobs@bioinformed.com)

Created:       Wed May 17 08:51:57 EDT 2006

Abstract:      Front-end for running TagZilla for SNPPlex and Illumina assay design

Compatibility: Python 2.4 and above

Requires:      No external dependencies, yet...

Revision:      $Id: $
'''

__program__   = 'HagZilla'
__authors__   = ['Kevin Jacobs (jacobs@bioinformed.com)']
__version__   = '1.1'
__copyright__ = 'Copyright 2006 Science Applications International Corporation ("SAIC")'

__license__ = '''
The software subject to this notice and license includes both human readable
source code form and machine readable, binary, object code form ("the
TagZilla Software"). The TagZilla Software was developed in conjunction with
the National Cancer Institute ("NCI") by NCI employees and employees or
contractors of SAIC. To the extent government employees are authors, any
rights in such works shall be subject to Title 17 of the United States Code,
section 105.

This TagZilla Software License (the "License") is between NCI and You. "You
(or "Your") shall mean a person or an entity, and all other entities that
control, are controlled by, or are under common control with the entity.
"Control" for purposes of this definition means (i) the direct or indirect
power to cause the direction or management of such entity, whether by
contract or otherwise, or (ii) ownership of fifty percent (50%) or more of
the outstanding shares, or (iii) beneficial ownership of such entity.

This License is granted provided that You agree to the conditions described
below. NCI grants You a non-exclusive, worldwide, perpetual, fully-paid-up,
no-charge, irrevocable, transferable and royalty-free right and license in
its rights in the TagZilla Software to (i) use, install, access, operate,
execute, copy, modify, translate, market, publicly display, publicly
perform, and prepare derivative works of the TagZilla Software; (ii)
distribute and have distributed to and by third parties the TagZilla
Software and any modifications and derivative works thereof; and (iii)
sublicense the foregoing rights set out in (i) and (ii) to third parties,
including the right to license such rights to further third parties. For
sake of clarity, and not by way of limitation, NCI shall have no right of
accounting or right of payment from You or Your sublicensees for the rights
granted under this License. This License is granted at no charge to You.

1. Your redistributions of the source code for the Software must retain the
   above copyright notice, this list of conditions and the disclaimer and
   limitation of liability of Article 6, below. Your redistributions in
   object code form must reproduce the above copyright notice, this list of
   conditions and the disclaimer of Article 6 in the documentation and/or
   other materials provided with the distribution, if any.

2. Your end-user documentation included with the redistribution, if any,
   must include the following acknowledgment: "This product includes
   software developed by SAIC and the National Cancer Institute." If You do
   not include such end-user documentation, You shall include this
   acknowledgment in the Software itself, wherever such third-party
   acknowledgments normally appear.

3. You may not use the names "The National Cancer Institute", "NCI" "Science
   Applications International Corporation" and "SAIC" to endorse or promote
   products derived from this Software. This License does not authorize You
   to use any trademarks, service marks, trade names, logos or product names
   of either NCI or SAIC, except as required to comply with the terms of
   this License.

4. For sake of clarity, and not by way of limitation, You may incorporate
   this Software into Your proprietary programs and into any third party
   proprietary programs. However, if You incorporate the Software into third
   party proprietary programs, You agree that You are solely responsible for
   obtaining any permission from such third parties required to incorporate
   the Software into such third party proprietary programs and for informing
   Your sublicensees, including without limitation Your end-users, of their
   obligation to secure any required permissions from such third parties
   before incorporating the Software into such third party proprietary
   software programs. In the event that You fail to obtain such permissions,
   You agree to indemnify NCI for any claims against NCI by such third
   parties, except to the extent prohibited by law, resulting from Your
   failure to obtain such permissions.

5. For sake of clarity, and not by way of limitation, You may add Your own
   copyright statement to Your modifications and to the derivative works,
   and You may provide additional or different license terms and conditions
   in Your sublicenses of modifications of the Software, or any derivative
   works of the Software as a whole, provided Your use, reproduction, and
   distribution of the Work otherwise complies with the conditions stated in
   this License.

6. THIS SOFTWARE IS PROVIDED "AS IS," AND ANY EXPRESSED OR IMPLIED
   WARRANTIES, (INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
   MERCHANTABILITY, NON-INFRINGEMENT AND FITNESS FOR A PARTICULAR PURPOSE)
   ARE DISCLAIMED. IN NO EVENT SHALL THE NATIONAL CANCER INSTITUTE, SAIC, OR
   THEIR AFFILIATES BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
   EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
   PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
   PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
   LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
   NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
   SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
'''

import os
import sys
import csv
import time
import sqlite3

from   itertools           import groupby,chain
from   operator            import itemgetter

from   biozilla.fileutils  import autofile,load_list
from   tagzilla            import launcher,locus_result_sequence


POPS = {'CEU'    : 'hapmap',
        'YRI'    : 'hapmap',
        'JPT'    : 'hapmap',
        'CHB'    : 'hapmap',
        'JPT+CHB': 'hapmap',
        'NHS'    : 'NHS',
       }

if 1:
  HAP_PEDIGREES = '/usr/local/share/hapmap/peds'
  HAP_GENOTYPES = '/usr/local/share/hapmap/build22/rs_strand/non-redundant/genotypes_chr%s_%s_r22_nr.b36.txt.gz'
  HAP_FORMAT    = 'hapmap'
else:
  HAP_PEDIGREES = '/usr/local/share/hapmap/peds'
  HAP_GENOTYPES = '/usr/local/share/hapmap/build21a/rs_strand/non-redundant/genotypes_chr%s_%s_r21a_nr.txt.gz'
  HAP_FORMAT    = 'hapmap'

NHS_PEDIGREES = None
NHS_GENOTYPES = '/u3/projects/CGEMS/Scans/Breast/1/current/genotypes/STUDY/subjects_STUDY_%s.ldat.gz'
NHS_MAP       = '/u3/projects/CGEMS/Scans/Breast/1/raw/snp.map'
NHS_FORMAT    = 'ldat'


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
    raise KeyError,'Cannot find gene "%s"' % gene
  elif len(genes) > 1:
    raise KeyError,'Gene not unique "%s"' % gene
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
    raise ValueError, 'Unknown gene orientation'


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
    args = '%s sequence/sublist.txt sequence/oblig.txt %s' % (snpname,sequencename)
    os.popen('%s %s 2>/dev/null' % (command,args)).read()
    return csv.reader(open(sequencename),dialect='excel-tab')
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
  if len(populations) > 1:
    command.append(' -M %s --multimethod global' % ','.join(populations))

  for population in populations:
    data = POPS[population]

    if data == 'hapmap':
      command.append('-f hapmap')
      command.append('-p %s' % HAP_PEDIGREES)
      command.append(HAP_GENOTYPES % (chromosome,population))
    elif data == 'NHS':
      command.append('-f ldat')
      command.append('-l %s' % NHS_MAP)
      command.append(NHS_GENOTYPES % chromosome)
    else:
      raise ValueError('Unknown population')

  return ' '.join(command)


def get_tags(outdir, project, gene, dprime, r2, populations, chromosome, snps, maf, include, exclude,
                     designfile, includefile):

  prefix = '%s_%s_%s' % (project,gene,populations)

  if designfile:
    designfile = os.path.abspath(designfile)

  if includefile:
    includefile = os.path.abspath(includefile)

  cwd = os.getcwd()

  try:
    project_chdir(outdir, project, gene)

    command = 'tagzilla'

    command += ' -s subset'
    file('subset','w').write('\n'.join(snps))

    if dprime:
      command += ' -d %s' % dprime

    if r2:
      command += ' -r %s' % r2

    if maf:
      command += ' -a %s' % maf

    if include and includefile:
      command += ' -i include'
      incs = (set(load_list(includefile)) | set(include)) & set(snps)
      file('include','w').write('\n'.join(sorted(incs)))
    elif includefile:
      command += ' -i %s' % includefile
    elif include:
      command += ' -i include'
      file('include','w').write('\n'.join(include))

    if exclude:
      command += ' -e exclude'
      file('exclude','w').write('\n'.join(exclude))

    if designfile:
      command += ' -D %s' % designfile

    command += ' ' + get_genotypes(populations,chromosome)

    command += ' -C maxsnp -O loci.out -b sum.txt -B bins.txt'
    command += ' 2>&1'

    #print command

    out = os.popen(command)
    err = out.read()
    bins = list(locus_result_sequence('loci.out',{},exclude))
  finally:
    os.chdir(cwd)

  tags = set()
  snps = set()
  obligate = set()
  exclset  = set()

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


def option_parser():
  import optparse

  usage = 'usage: %prog [options] genofile...'
  parser = optparse.OptionParser(usage=usage)

  parser.add_option('-d', '--design',   dest='design', metavar='FILE',
                      help='Design scores (optional)')
  parser.add_option('-g', '--genomedb',   dest='genomedb', metavar='FILE', default='genome.db',
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
  return parser


def hagzilla():
  parser = option_parser()
  options,args = parser.parse_args()

  if len(args) != 2:
    parser.print_help()
    return

  con = sqlite3.connect(options.genomedb)
  outdir = args[0]

  taskfile   = csv.reader(autofile(args[1]),dialect='excel-tab')
  header     = taskfile.next()
  tasks      = sorted( extend(strip(t),12) for t in taskfile )
  genecounts = {}
  genestuff  = {}
  genefail   = {}

  print '\t'.join(['Project','Population','Region','Chr','Start','Stop',
                   'SNPs', 'Qualifying SNPs', 'Tags', 'Ob. Include', 'Ob. Exclude'])

  for project,ptasks in groupby(tasks,itemgetter(0)):
    if not project:
      continue

    ptasks = list(ptasks)
    ptags = set()
    for project,population,gene,chromosome,start,stop,strand,dprime,r2,maf,include,exclude in ptasks:
      if chromosome or start or stop:
        start = int(start)
        stop  = int(stop)
      else:
        try:
          row = query_gene_by_name(con, gene)
        except KeyError,e:
          print >> sys.stderr, 'Invalid gene: %s' % e.args[0]
          continue

        if chromosome.strip() and chromosome.strip() != row[2]:
          print >> sys.stderr, 'Warning: Chromosome does not match for %s (%s != %s)' % \
                                       (gene,chromosome,row[2])

        chromosome = row[2]
        start      = row[3]
        stop       = row[4]
        strand     = row[5]

        start,stop = gene_margin(start,stop,strand,upstream=options.upstream,
                                                   downstream=options.downstream)

      include = set(i for i in include.split(',') if i)
      exclude = set(e for e in exclude.split(',') if e)

      if options.dryrun:
        print '\t'.join(map(str,[project,population,gene,chromosome,start,stop,strand]))
        continue

      if None in (start,stop):
        snps = qsnps = tags = obligate = undesignable = set()
      else:
        snps = set(snp[0] for snp in get_snps(con, population, chromosome, start, stop))

        qsnps,tags,obligate,excl = get_tags(outdir,project,gene,dprime,r2,population,
                                            chromosome,snps,maf,include,exclude,options.design,options.include)

      #print 'Tags',len(tags)
      #print 'Obligates',len(obligate)

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
      genestuff[gene] = seqs

      d = map(str,[project,population,gene,chromosome,start,stop,
                   len(snps),len(qsnps),len(tags),len(obligate),len(excl)])
      print '\t'.join(d)

    designfile = file(project_file(outdir,project,'design.lst'),'w')
    for rs in ptags:
      designfile.write('%s\n' % rs)

    if options.getsequence:
      designfile = file(project_file(outdir,project,'design.fasta'),'w')
      for rs,seq in get_sequences(con, ptags):
        designfile.write('>%s\n%s\n' % (rs,seq))

    pdir = project_path(outdir,project)
    print >> sys.stderr, os.popen('binsum "%s"/*/loci.out > "%s"/sum.out 2>/dev/null' % (pdir,pdir)).read()

    sumfile = file(project_file(outdir,project,'genesum.out'),'w')
    for gene,n in sorted(genecounts.iteritems()):
      sumfile.write('%s\t%d\n' % (gene,n))

    sumfile = file(project_file(outdir,project,'genestuff.out'),'w')
    for gene,seqs in sorted(genestuff.iteritems()):
      for rs,seq in seqs:
        sumfile.write('%s\t%s\t%s\n' % (gene,rs,seq))


if __name__ == '__main__':
  hagzilla()

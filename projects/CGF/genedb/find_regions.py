# -*- coding: utf-8 -*-
'''
File:          find_regions.py

Authors:       Kevin Jacobs (jacobs@bioinformed.com)

Created:       Wed May 17 08:51:57 EDT 2006

Abstract:      Front-end for running TagZilla for SNPPlex and Illumina assay design

Compatibility: Python 2.4 and above

Requires:      No external dependencies, yet...

Revision:      $Id: find_regions.py 315 2006-07-31 14:43:47Z jacobske $
'''

__program__   = 'find_regions'
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

import sqlite3
import time
import csv
import os

import seq_comp

from tagzilla import *

HAPMAP_PATH= '/home/jacobske/projects/CGEMS/hapmap'
GENOTYPES  = os.path.join(HAPMAP_PATH,'build20','non-redundant')
PEDIGREES  = os.path.join(HAPMAP_PATH,'pedigrees','peds')
POPS       = ['CEU','YRI','JPT','CHB','JPT+CHB']

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
  cur.execute(sql, [gene])
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
  WHERE    s.chromosome = ?
    AND    ? BETWEEN s.chrStart AND s.chrEnd
    AND    s.featureType='GENE'
    AND    s.submitGroup='reference'
  ORDER BY s.chromosome,MIN(s.chrStart,s.chrEnd);
  '''
  cur = con.cursor()
  cur.execute(sql, chr, loc)
  return cur.fetchall()


def gene_margin(gene, upstream=20000, downstream=10000):
  if None in tuple(gene[2:5]):
    return None,None
  if gene[5] == '+':
    return gene[-3]-upstream,gene[-2]+downstream
  elif gene[5] == '-':
    return gene[-3]-downstream,gene[-2]+upstream
  else:
    raise ValueError, 'Unknown gene orientation for %s=%s' % (gene[0],gene[5])


def get_snps(con, chromosome, start, stop):
  sql = '''
  SELECT   lname, chromosome, location
  FROM     snp
  WHERE    chromosome=?
    AND    location BETWEEN ? AND ?;
  '''
  cur = con.cursor()
  cur.execute(sql,[chromosome,start,stop])
  return cur.fetchall()


def escape(s):
  return "'%s'" % s.replace("'","''")


def extend(s,n):
  m = len(s)
  if m < n:
    s = list(s)
    s.extend(['']*(n-m))
  return s


def option_parser():
  import optparse

  usage = 'usage: %prog [options] genofile...'
  parser = optparse.OptionParser(usage=usage, add_help_option=False)

  parser.add_option('-h', '--help', dest='help', action='store_true',
                          help='show this help message and exit')
  parser.add_option('--license', dest='license', action='store_true',
                          help="show program's copyright and license terms and exit")
  return parser


def find_regions(options,args):
  con = sqlite3.connect('genome35_v3.db')

  out = csv.writer(sys.stdout,dialect='excel-tab')
  out.writerow( ['GENE','CHROMOSOME','GENE START','GENE END','LOCUS','LOCATION'] )

  for arg in args:
    genefile = csv.reader(autofile(arg),dialect='excel-tab')
    header= genefile.next()
    genes = [ line[0] for line in genefile if line and line[0] ]

    for gene in genes:
      try:
        row = query_gene_by_name(con, gene)
      except KeyError,e:
        print >> sys.stderr, 'Invalid gene: %s' % e.args[0]
        continue

      if gene != row[0]:
        gene = '%s (%s)' % (gene,row[0])

      chromosome = row[2]
      start,stop = gene_margin(row,20000,10000)
      snps = get_snps(con, chromosome, start, stop)

      for lname,chromosome,location in snps:
        out.writerow( [gene,chromosome,row[3],row[4],lname,location] )


def main():
  launcher(find_regions, option_parser, **globals())


if __name__ == '__main__':
  main()

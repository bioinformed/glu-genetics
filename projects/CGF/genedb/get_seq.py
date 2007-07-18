# -*- coding: utf-8 -*-
'''
File:          get_seq.py

Authors:       Kevin Jacobs (jacobs@bioinformed.com)

Created:       Wed May 17 08:51:57 EDT 2006

Abstract:

Compatibility: Python 2.4 and above

Requires:      No external dependencies, yet...

Revision:      $Id: get_seq.py 315 2006-07-31 14:43:47Z jacobske $
'''

__program__   = 'get_seq'
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
GENOTYPES  = os.path.join(HAPMAP_PATH,'build20/non-redundant')
PEDIGREES  = os.path.join(HAPMAP_PATH,'pedigrees', 'peds')
POPS       = ['CEU','YRI','JPT','CHB','JPT+CHB']

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


def get_seq(options,args):
  con = sqlite3.connect('genome35.db')

  snps = set()
  for arg in args:
    snps.update(row[0] for row in csv.reader(autofile(arg),dialect='excel-tab') if row)

  degenerate = set('wymkwsbdhvnx')
  dcount = 0
  n = 16
  for rs,seq in get_sequences(con, snps):
    i = seq.index('[')
    left = seq[i-n:i].lower()
    i = seq.index(']')
    right = seq[i+1:i+n+1].lower()
    bases = set(left+right)
    if bases & degenerate:
      dcount += 1
    print '>%s\n%s' % (rs,seq)

  print 'Degenerate:',dcount

def main():
  launcher(get_seq, option_parser, **globals())


if __name__ == '__main__':
  main()

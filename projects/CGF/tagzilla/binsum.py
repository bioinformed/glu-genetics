# -*- coding: utf-8 -*-
'''
File:          binsum.py

Authors:       Kevin Jacobs (jacobs@bioinformed.com)

Created:       November 8, 2005

Abstract:      Script to generate bin summaries from TagZilla output files.

Compatibility: Python 2.4 and above

Requires:      Matching version of tagzilla.py

Revision:      $Id: binsum.py 474 2006-12-26 20:01:57Z jacobske $
'''

__program__   = 'TagZilla binsum'
__authors__   = ['Kevin Jacobs (jacobs@bioinformed.com)']
__version__   = '1.2alpha'
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


import heapq
from   tagzilla import *


def imerge(*its):
  pqueue = []
  for i in map(iter, its):
    try:
      pqueue.append((i.next(), i.next))
    except StopIteration:
      pass

  heapq.heapify(pqueue)

  while pqueue:
    val, it = pqueue[0]
    yield val
    try:
      heapq.heapreplace(pqueue, (it(), it))
    except StopIteration:
      heapq.heappop(pqueue)


def option_parser():
  import optparse

  usage = 'usage: %prog [options] genofile...'
  parser = optparse.OptionParser(usage=usage, add_help_option=False)

  parser.add_option('-h', '--help', dest='help', action='store_true',
                          help='show this help message and exit')
  parser.add_option('--license', dest='license', action='store_true',
                          help="show program's copyright and license terms and exit")
  parser.add_option('-b', '--summary', dest='sumfile', metavar='FILE', default='-',
                          help="Output summary tables FILE (default='-' for standard out)")
  parser.add_option('-B', '--bininfo', dest='bininfo', metavar='FILE',
                          help='Output summary information about each bin to FILE')
  parser.add_option('-H', '--histomax', dest='histomax', metavar='N', type='int', default=10,
                          help='Largest bin size output in summary histogram output (default=10)')
  parser.add_option('-s', '--subset', dest='subset', metavar='FILE', default='',
                          help='File containing loci to filter bins.  Only bins that contain one or more of these loci will be included in the output.')
  parser.add_option('-t', '--targetbins', dest='targetbins', metavar='N', type='int', default=0,
                          help='Stop when N bins have been selected (default=0 for unlimited)')
  parser.add_option('-T', '--targetloci', dest='targetloci', metavar='N', type='int', default=0,
                          help='Stop when N loci have been tagged (default=0 for unlimited)')

  return parser


def binsum(options,args):
  infofile = None
  if options.bininfo:
    infofile = autofile(options.bininfo, 'w', hyphen=sys.stdout)

  if options.bininfo or options.sumfile:
    bininfo = BinInfo(infofile, options.histomax+1)
  else:
    bininfo = NullBinInfo()

  sumfile = autofile(options.sumfile, 'w', hyphen=sys.stdout)

  if [infofile,sumfile].count(sys.stdout) > 1:
    print >> sys.stderr, 'ERROR: More than one output file directed to standad out.'
    return

  subset  = set()
  exclude = set()

  if options.subset:
    read_snp_list(options.subset, subset)

  locusmap = {}
  results = [ locus_result_sequence(filename, locusmap, exclude) for filename in args ]
  results = imerge(*results)

  popdtags = {}
  populations = set()
  binned_loci = 0
  binnum = None
  for population,bin in results:
    if subset and not subset.intersection(bin):
      continue

    if binnum != bin.binnum:
      binnum = bin.binnum
      popdtags[bin.disposition] = popdtags.get(bin.disposition,0) + 1

    populations.add(population)

    bin_qualifier(bin, binned_loci, options)

    binned_loci += len(bin)
    bininfo.emit_bin(bin, locusmap, exclude, population)

  # Emit useful bin summary table
  for population in sorted(populations):
    bininfo.emit_summary(sumfile, population)

  if len(populations) > 1:
    bininfo.emit_multipop_summary(sumfile,popdtags)


def main():
  launcher(binsum, option_parser, **globals())


if __name__ == '__main__':
  main()

# -*- coding: utf-8 -*-
'''
File:          coverage.py

Authors:       Kevin Jacobs (jacobs@bioinformed.com)

Created:       November 8, 2005

Abstract:      Evaluate maximum coverage of a set of tags

Compatibility: Python 2.4 and above

Requires:      No external dependencies, yet...

Revision:      $Id: coverage.py 474 2006-12-26 20:01:57Z jacobske $
'''

__program__   = 'TagZilla coverage'
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

__accelerators__ = ['tagzillac']

from tagzilla import *


def option_parser():
  usage = 'usage: %prog [options] tagfile genofile...'
  parser = TagZillaOptionParser(usage=usage, add_help_option=False)

  parser.add_option('-h', '--help', dest='help', action='store_true',
                        help='show this help message and exit')
  parser.add_option('--license', dest='license', action='store_true',
                          help="show program's copyright and license terms and exit")
  parser.add_option('--profile', dest='profile', metavar='P', help=optparse.SUPPRESS_HELP)

  inputgroup = optparse.OptionGroup(parser, 'Input options')

  inputgroup.add_option('-f', '--format', dest='format', metavar='NAME', default='',
                          help='Format for genotype/pedigree or ld input data.  Values: hapmap (default), linkage, festa, prettybase, raw.')
  inputgroup.add_option('-l', '--loci', dest='loci', metavar='FILE',
                          help='Locus description file for input in Linkage format')
  inputgroup.add_option('-p', '--pedfile', dest='pedfile', metavar='FILE', action='append',
                          help='Pedigree file for HapMap or PrettyBase data files (optional)')
  inputgroup.add_option('-e', '--excludetag', dest='exclude', metavar='FILE', default='',
                          help='File containing loci that are excluded from being a tag')
  inputgroup.add_option('-s', '--subset', dest='subset', metavar='FILE', default='',
                          help='File containing loci that define the subset to be analyzed of the loci that are read')
  inputgroup.add_option('-R', '--range', dest='range', metavar='S-E,...', default='',
                          help='Ranges of genomic locations to analyze, specified as a comma seperated list of start and '
                               'end coordinates "S-E".  If either S or E is not specified, then the ranges are assumed '
                               'to be open.  The end coordinate is exclusive and not included in the range.')
  inputgroup.add_option('-D', '--designscores', dest='designscores', metavar='FILE', type='str', action='append',
                          help='Read in design scores or other weights to use as criteria to choose the optimal tag for each bin')
  inputgroup.add_option('-L', '--limit', dest='limit', metavar='N', type='int', default=0,
                          help='Limit the number of loci considered to N for testing purposes (default=0 for unlimited)')

  outputgroup = optparse.OptionGroup(parser, 'Output options')

  outputgroup.add_option('-o', '--output', dest='outfile', metavar='FILE', default='-',
                          help="Output tabular LD information for bins to FILE ('-' for standard out)")

  genoldgroup = optparse.OptionGroup(parser, 'Genotype and LD estimation options')

  genoldgroup.add_option('-a', '--minmaf', dest='maf', metavar='FREQ', type='float', default=0.05,
                          help='Minimum minor allele frequency (MAF) (default=0.05)')
  genoldgroup.add_option('-A', '--minobmaf', dest='obmaf', metavar='FREQ', type='float', default=None,
                          help='Minimum minor allele frequency (MAF) for obligate tags (defaults to -a/--minmaf)')
  genoldgroup.add_option('-c', '--mincompletion', dest='mincompletion', metavar='N', default=0, type='int',
                          help='Drop loci with less than N valid genotypes. Default=0')
  genoldgroup.add_option(      '--mincompletionrate', dest='mincompletionrate', metavar='N%', default=0, type='float',
                          help='Drop loci with completion rate less than N% (0-100). Default=0')
  genoldgroup.add_option('-m', '--maxdist', dest='maxdist', metavar='D', type='int', default=200,
                          help='Maximum inter-marker distance in kb for LD comparison (default=200)')
  genoldgroup.add_option('-P', '--hwp', dest='hwp', metavar='p', default=None, type='float',
                          help='Filter out loci that fail to meet a minimum signficance level (pvalue) for a '
                               'test Hardy-Weignberg proportion (no default)')

  bingroup = optparse.OptionGroup(parser, 'Binning options')

  bingroup.add_option('-d', '--dthreshold', dest='d', metavar='DPRIME', type='float', default=0.,
                          help='Minimum d-prime threshold to output (default=0)')
  bingroup.add_option('-r', '--rthreshold', dest='r', metavar='N', type='float', default=0,
                          help='Minimum r-squared threshold to output (default=0)')

  parser.add_option_group(inputgroup)
  parser.add_option_group(outputgroup)
  parser.add_option_group(genoldgroup)
  parser.add_option_group(bingroup)

  return parser


def coverage(options,args):
  tags    = set()
  subset  = set()
  exclude = set()

  read_snp_list(args[0][1], tags)

  args = args[1:]

  if options.subset:
    read_snp_list(options.subset, subset)

  if options.exclude:
    read_snp_list(options.exclude, exclude)

  designscores = build_design_score(options.designscores)

  locusmap = {}
  options.multipopulation = None
  ldpairs = generate_ldpairs(args, locusmap, set(), subset, tags, options)

  besttag = dict( (tag,(tag,1)) for tag in tags )
  missing = '',-1
  for pairs in ldpairs:
    for lname1,lname2,r2,dprime in pairs:
      if designscores:
        if designscores.get(lname1,0) < epsilon:
          exclude.add(lname1)
        if designscores.get(lname2,0) < epsilon:
          exclude.add(lname2)

      for l1,l2 in [(lname1,lname2),(lname2,lname1)]:
        if l2 in tags:
          best_locus,best_r2 = besttag.get(l1,missing)
          if l2 not in exclude and r2 > best_r2:
            besttag[l1] = l2,r2

  outfile = autofile(options.outfile, 'w', hyphen=sys.stdout)
  outfile.write('LNAME\tTAG\tBEST RSQUARED\n')
  for lname in locusmap:
    tag,r2 = '',''
    if lname in besttag:
      tag,r2 = besttag[lname]
      r2 = sfloat(r2)
    outfile.write('%s\t%s\t%s\n' % (lname,tag,r2))


def main():
  launcher(coverage, option_parser, **globals())


if __name__ == '__main__':
  main()

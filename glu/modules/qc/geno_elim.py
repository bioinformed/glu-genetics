# -*- coding: utf-8 -*-

__gluindex__  = False
__abstract__  = 'Perform genotype-elimination on a set of genotypes based on their relationships'
__copyright__ = 'Copyright (c) 2007-2009, BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.'
__license__   = 'See GLU license for terms by running: glu license'
__revision__  = '$Id$'

import sys

from   itertools         import izip,chain
from   collections       import defaultdict

from   glu.lib.fileutils import table_reader,table_writer
from   glu.lib.genolib   import load_genostream,geno_options


errbylochead1 = ['','LEVEL_1_ERRORS','','LEVEL_2_ERRORS','','LEVEL_3_ERRORS']
errbylochead2 = ['LOCUS','INDIVIDUAL','FAMILY','INDIVIDUAL','FAMILY','FAMILY']
errbypedhead  = ['FAMILY','INDIVIDUAL','LEVEL_1_ERRORS','LEVEL_2_ERRORS','LEVEL_3_ERRORS']


def option_parser():
  from glu.lib.glu_argparse import GLUArgumentParser

  parser = GLUArgumentParser(description=__abstract__)

  parser.add_argument('genotypes', help='Input genotype file')

  geno_options(parser,input=True)

  parser.add_argument('-o', '--errdetails',  default='-', metavar='FILE',
                    help="The output file containing genotype matrix after elimination, '-' for standard out")
  parser.add_argument('--locsum', default='-', metavar='FILE',
                    help="The output file containing summary of errors by locus, '-' for standard out")
  parser.add_argument('--locdet', default='-', metavar='FILE',
                    help="The output file containing details of errors by locus, '-' for standard out")
  parser.add_argument('--pedsum', default='-', metavar='FILE',
                    help="The output file containing summary of errors by pedigree, '-' for standard out")
  parser.add_argument('--peddet', default='-', metavar='FILE',
                    help="The output file containing details of errors by pedigree, '-' for standard out")

  return parser


def mkgeno(a,b):
  return min(a,b)+max(a,b)


def build_children_genos(p1geno,p2geno):
  '''
  Build a set of possible children genotypes from a pair of parents genotypes

  @param p1geno: genotype of one of the parents
  @type  p1geno: str
  @param p2geno: genotype of the other parent
  @type  p2geno: str
  @return      : a set of possible children genotypes
  @rtype       : set

  >>> p1geno,p2geno = 'AA','AG'
  >>> build_children_genos(p1geno,p2geno)
  set(['AA', 'AG'])
  '''
  return set(mkgeno(a,b) for a in p1geno for b in p2geno)


def build_uninformative_genoset(alleles):
  '''
  Build a set of possible unphased genotypes given a set of alleles
  @param alleles: a set of alleles
  @type  alleles: set
  @return       : a set of possible genotypes
  @rtype        : set

  >>> alleles = set(['A','G'])
  >>> build_uninformative_genoset(alleles)
  set(['AA', 'GG', 'AG'])
  '''
  return set(mkgeno(a,b) for a in alleles for b in alleles)


def build_nuclear_families(pedfile,individuals):
  '''
  Build nuclear families based on a pedigree file
  @param     pedfile: pedigree file
  @type      pedfile: str
  @param individuals: genotyped individuals
  @type  individuals: list
  @return           : nuclear families and the total number of built families
  @rtype            : tuple

  >>> import StringIO
  >>> contents = "fam101\\tind01\\t0\\t0\\n\
  fam101\\tind02\\t0\\t0\\n\
  fam101\\tind11\\tind01\\tind02\\n\
  fam101\\tind12\\tind01\\tind02\\n"
  >>> f = StringIO.StringIO(contents)
  >>> ind = ['fam101_ind01','fam101_ind02','fam101_ind11','fam101_ind12']
  >>> fams,n=build_nuclear_families(f,ind)
  >>> dict(fams)
  {('fam101_ind01', 'fam101_ind02'): ['fam101_ind11', 'fam101_ind12']}
  >>> n
  1
  '''

  def concat_id(fam,id):
    return '%s_%s' % (fam.strip(),id.strip())

  def get_field(row,n,default=None):
    try:
      return row[n]
    except IndexError:
      return default

  typed = {}
  for individual in individuals:
    typed[individual] = 1

  peds = table_reader(pedfile)

  nfam = defaultdict(list)
  famset = set()
  for ped in peds:
    if len(ped) < 2:
      continue

    fam,ind = ped[:2]

    parent1 = get_field(ped,2,'0')
    parent2 = get_field(ped,3,'0')

    if '0' in (parent1,parent2):
      continue

    if fam:
      ind=concat_id(fam,ind)
      parent1=concat_id(fam,parent1)
      parent2=concat_id(fam,parent2)
      if parent1 not in typed or parent2 not in typed or ind not in typed:
        continue

    parents = min(parent1,parent2),max(parent1,parent2)
    nfam[parents].append(ind)
    famset.add(fam)

  return nfam,len(famset)


def build_genodict(genos,individuals):
  '''
  Make a possible genotype list for each individual
  @param       genos: genotypes
  @type        genos: list
  @param individuals: individuals
  @type  individuals: list
  @return           : genotype dictionary for each individual
  @rtype            : dict

  >>> genos = ['AA','AG','GG','AA','','']
  >>> inds  = ['P1','P2','P3','P4','P5','P6']
  >>> genodict = build_genodict(genos,inds)
  >>> for ind,geno in genodict.iteritems():
  ...   print ind,geno
  P2 set(['AG'])
  P3 set(['GG'])
  P1 set(['AA'])
  P6 set(['AA', 'GG', 'AG'])
  P4 set(['AA'])
  P5 set(['AA', 'GG', 'AG'])
  '''

  alleles = set()
  for geno in genos:
    if geno:
      alleles.update(geno)

  missing = build_uninformative_genoset(alleles)

  genodict = {}
  for ind,geno in izip(individuals,genos):
    if geno:
      genodict[ind] = set([mkgeno(*geno)])
    else:
      genodict[ind] = missing.copy()

  return genodict


def genotype_elimination(nfams,genos,individuals,out):
  '''
  Implement the sequential Genotype-Elimination Algorithm

  @param       nfams: nuclear families
  @type        nfams: dict
  @param       genos: genotypes
  @type        genos: list
  @param individuals: individuals
  @type  individuals: list
  @param         out: output file
  @type          out: file
  '''

  def norm_geno(g):
    if not g:
      return ''
    else:
      return intern(''.join(g))

  errbyloc1 = {}
  errbyloc2 = {}
  errbyloc  = {}

  for ct,(lname,locusgenos) in enumerate(genos):
    print '------------------------line: %d' % (ct+1)
    locusgenos = map(norm_geno,locusgenos)
    genos = build_genodict(locusgenos,individuals)
    prev_genos = {}
    while genos!=prev_genos:
      prev_genos = genos.copy()
      for parents,children in nfams.iteritems():
        eliminate_nuclear_family(parents,children,lname,genos,errbyloc1,errbyloc2)

    build_errbyloc(individuals,lname,genos,errbyloc)
    if out:
      emit_gea_details(out, individuals, lname, genos)

  return errbyloc1,errbyloc2,errbyloc


def build_errbyloc(individuals, lname, genos, errbyloc):
  '''
  Output genotype errors by locus after elimination
  '''
  # FIXME: intern/split?
  for ind in individuals:
    if not genos[ind]:
      errbyloc[lname].add(intern(ind.split('_')[0]))


def emit_gea_details(outdetails, individuals, lname, genos):
  '''
  Output genotypes after elimination
  '''
  locusgenos=[lname]
  for ind in individuals:
    locusgenos.append(sorted(genos[ind]))
  outdetails.writerow(locusgenos)


def eliminate_nuclear_family(parents,children,lname,genos,errbyloc1,errbyloc2):
  '''
  Genotype elimination within each nuclear family

  @param   parents: ids for both parents
  @type    parents: tuple
  @param  children: ids for all the children
  @type   children: list
  @param     lname: locus name
  @type      lname: str
  @param     genos: current possible genotypes list for each person
  @type      genos: dict
  @param errbyloc1: ids for individuals with level 1 errors by locus
  @type  errbyloc1: dict
  @param errbyloc2: ids for individuals with level 2 errors by locus
  @type  errbyloc2: dict

  >>> parents = '1140_01','1140_02'
  >>> children = ['1140_11','1140_12']
  >>> lname = 'rs1167772'
  >>> genos = {'1140_01':set(['AA']),'1140_02':set(['AA','AG','GG']),'1140_11':set(['AG']),'1140_12':set(['AA','AG','GG'])}
  >>> errbyloc1=defaultdict(set)
  >>> errbyloc2=defaultdict(set)
  >>> eliminate_nuclear_family(parents,children,lname,genos,errbyloc1,errbyloc2)
  >>> sorted( (i,sorted(g)) for i,g in genos.iteritems() )
  [('1140_01', ['AA']), ('1140_02', ['AG', 'GG']), ('1140_11', ['AG']), ('1140_12', ['AA', 'AG'])]
  '''
  p1,p2 = parents

  do_elim = False
  for ind in chain(parents,children):
    if not genos[ind] or len(genos[ind]) > 1:
      do_elim = True
      break

  if do_elim:
    eliminate_single_genotypes(parents,children,genos)

    # Set saved genos to empty for all family members
    saved_genos = dict( (i,set()) for i in chain(parents,children) )

    # Eliminate p1 against p2 and children
    eliminate_parents(p1,genos[p1],p2,genos[p2],children,genos,saved_genos)

    # Eliminate unsaved genotypes of p2 against p1 and children
    p2genos   = genos[p2]-saved_genos[p2]
    eliminate_parents(p2,p2genos,p1,genos[p1],children,genos,saved_genos)

    # Eliminate children against parents
    eliminate_children(parents,children,genos,saved_genos)

  elif not level1_check(parents,children,lname,genos,errbyloc1):
    level2_check(parents,children,lname,genos,errbyloc2)


def level1_check(parents,children,lname,genos,errbyloc1):
  '''
  Find genotype errors at parent-offspring level
  '''

  def pair_check(p,c):
    return ( set(getitem(p)) & set(getitem(c)) )

  err = False
  for p in parents:
    for c in children:
      pgeno,cgeno=genos[p],genos[c]
      if not pair_check(pgeno,cgeno):
        errbyloc1[lname].add( (p,c) )
        err = True
  return err


def level2_check(parents,children,lname,genos,errbyloc2):
  '''
  Find genotype errors at parent-parent-offspring level
  '''
  p1geno,p2geno=getitem(genos[parents[0]]),getitem(genos[parents[1]])
  for child in children:
    if not compatible_genos(p1geno,p2geno,[child],genos):
      errbyloc2[lname].add( (parents[0],parents[1],child) )


def eliminate_parents(p1,p1genos,p2,p2genos,children,genos,saved_genos):
  '''
  Eliminate superfluous genotypes of parent1 with respect to parent2 and
  their children.  Parent2 must be processed separately, since compatibility
  is only thoroughly verified for parent1's genotypes.  Thus, this function
  may be called with parent1 and parent2 reversed and all saved genotypes of
  parent2 spared from consideration.

  @param           p1: first parent
  @type            p1: str
  @param      p1genos: genotypes of the first parent
  @type       p1genos: set
  @param           p2: second parent
  @type            p2: str
  @param      p2genos: genotypes of the second parent
  @type       p2genos: set
  @param     children: children of the nuclear family
  @type      children: list
  @param        genos: current possible genotypes list for each person
  @type         genos: dict
  @param  saved_genos: current saved genotypes set for each person
  @type   saved_genos: dict
  '''
  for p1geno in p1genos:
    for p2geno in p2genos:
      matchedgenos = compatible_genos(p1geno,p2geno,children,genos)

      if matchedgenos:
        saved_genos[p1].add(p1geno)
        saved_genos[p2].add(p2geno)
        for child,matched in matchedgenos.iteritems():
          saved_genos[child].update(matched)
        break

  # p1 is fully eliminated, so all unsaved genotypes may be discarded
  genos[p1] = saved_genos[p1]


def eliminate_children(parents,children,genos,saved_genos):
  '''
  Eliminate superfluous genotypes of the children

  @param       parents: ids for both parents
  @type        parents: tuple
  @param      children: ids for all the children
  @type       children: list
  @param         genos: current possible genotypes set for each person
  @type          genos: dictionary
  @param   saved_genos: current saved genotypes set for each person
  @type    saved_genos: dictionary
  '''
  for child in children:
    for cgeno in genos[child] - saved_genos[child]:
      calleles = set(cgeno)
      pgenos = ( (p1,p2) for p1 in saved_genos[parents[0]] if set(p1)&calleles
                         for p2 in saved_genos[parents[1]] if set(p2)&calleles)

      for p1geno,p2geno in pgenos:
        matchedgenos = compatible_genos(p1geno,p2geno,children,genos)
        if matchedgenos:
          for child,matched in matchedgenos.iteritems():
            saved_genos[child].update(matched)
        else:
          genos[child].discard(cgeno)
          saved_genos[child].discard(cgeno)

  # children are fully eliminated, so all unsaved genotypes may be discarded
  for child in children:
    genos[child] = saved_genos[child]


def compatible_genos(p1geno,p2geno,children,genos):
  '''
  Check to see if the parents and children genotypes are matching

  @param   p1geno: first parent genotype
  @type    p1geno: str
  @param   p2geno: second parent genotype
  @type    p2geno: str
  @param children: list of children
  @type  children: list
  @param    genos: individual to genotypes dictionary
  @type     genos: dict
  @return        : compatible or not and matched children genotypes
  @rtype         : tuple
  '''
  cgenoset = build_children_genos(p1geno,p2geno)
  matchedcgenos = {}
  for child in children:
    matchedcgenos[child] = cgenoset&genos[child]
    if not matchedcgenos[child]:
      return None
  return matchedcgenos


def eliminate_inconsistency(geno, genoset):
  '''
  Eliminate incompatible genotypes given a fixed genotype against a set of
  either parent or offspring genotypes such at least one allele from the
  fixed genotype must appear in each genotype of the genoset.  Any
  impossible genotypes are removed from the genotype set.

  @param    geno: Fixed genotype of either a parent or child in a nuclear family
  @type     geno: string
  @param genoset: Set of either parents or children of the fixed genotype
                  that are to be eliminated
  @type  genoset: set of string

  >>> g = set(['AA','AB','BB'])
  >>> eliminate_inconsistency('AA',g)
  >>> sorted(g)
  ['AA', 'AB']
  >>> g = set(['AA','AB','BB'])
  >>> eliminate_inconsistency('AB',g)
  >>> sorted(g)
  ['AA', 'AB', 'BB']
  >>> g = set(['AA','AB','BB'])
  >>> eliminate_inconsistency('CC',g)
  >>> sorted(g)
  []
  '''
  alleles = set(geno)
  for po_geno in list(genoset):
    po_alleles = set(po_geno)
    if not alleles & po_alleles:
      genoset.remove(po_geno)


def getitem(seq):
  assert len(seq)==1
  return iter(seq).next()


def eliminate_single_genotypes(parents,children,genos):
  '''
  Eliminate incompatible genotypes for ambiguous individuals based on the
  unambiguous family member genotypes with a single genotype

  @param  parents: ids for both parents
  @type   parents: tuple
  @param children: ids for all the children
  @type  children: list
  @param    genos: genotype list for each person
  @type     genos: dict

  >>> parents  = 'p1','p2'
  >>> children = ['c1','c2']
  >>> genos    = {'p1':set(['AA','AG','GG']),'p2':set(['AG']),'c1':set(['AA','AG','GG']),'c2':set(['AA'])}
  >>> err = eliminate_single_genotypes(parents,children,genos)
  >>> expected_genos = {'p2':set(['AG']),'c2':set(['AA']),'c1':set(['AA','AG','GG']),'p1':set(['AA','AG'])}
  >>> if expected_genos == genos:
  ...   print True
  True
  '''
  for p in parents:
    pgenos = genos[p]

    for child in children:
      cgenos = genos[child]

      # if parent is fixed, check the genotype against each child
      if len(pgenos) == 1:
        pgeno = getitem(pgenos)
        eliminate_inconsistency(pgeno,cgenos)

      # if child is fixed, check its genotypes against each parent
      elif len(cgenos) == 1:
        cgeno = getitem(cgenos)
        eliminate_inconsistency(cgeno,pgenos)


def invert_dict(errbyloc,simple=True):
  '''
  Invert a dict
  '''
  errbyped = defaultdict(set)
  for locus,fams in errbyloc.iteritems():
    if fams:
      for inds in fams:
        if simple:
          fam=getitem(inds)
          errbyped[inds].add(locus)
        else:
          for ind in inds:
            errbyped[ind].add(locus)
  return errbyped


def count_by_family(errdict):
  '''
  Convert error count by individual to by family
  '''
  newdict = defaultdict(set)

  for key,valset in errdict.iteritems():
    for vals in valset:
      for val in vals:
        newdict[key].add(val.split('_')[0])

  return newdict


def concat_element(dataset):
  '''
  Format a set of tuples into a list
  '''
  return ';'.join(' '.join(d) for d in dataset)


def make_errbyloc3(errbyloc,errbylocf1,errbylocf2):
  '''
  Retrieve level 3 errors by subtracting level 1 and 2 errors from the total errors
  '''
  # FIXME: This used to not copy.  Who knows what was intended.
  errbyloc3=errbyloc.copy()
  for key,vals in errbyloc3.iteritems():
    # FIXME: Is this guaranteed to be non-negative?
    errbyloc3[key] -= errbylocf2.get(key,0) + errbylocf1.get(key,0)
  return errbyloc3


def emit_errbyloc(locus,errbyloc1,errbylocf1,errbyloc2,errbylocf2,errbyloc3):
  '''
  Output genotype errors by locus
  '''
  dind1,dfam1,dind2,dfam2,dfam3='','','','',''
  sind1,sfam1,sind2,sfam2,sfam3=0,0,0,0,0

  if locus in errbyloc1:
    dind1=concat_element(errbyloc1[locus])
    dfam1=';'.join(errbylocf1[locus])
    sind1=len(errbyloc1[locus])
    sfam1=len(errbylocf1[locus])

  if locus in errbyloc2:
    dind2=concat_element(errbyloc2[locus])
    dfam2=';'.join(errbylocf2[locus])
    sind2=len(errbyloc2[locus])
    sfam2=len(errbyloc2[locus])

  if locus in errbyloc3:
    dfam3=';'.join(errbyloc3[locus])
    sfam3=len(errbyloc3[locus])

  rowd=[locus,dind1,dfam1,dind2,dfam2,dfam3]
  rows=[locus,sind1,sfam1,sind2,sfam2,sfam3]

  return rowd,rows


def emit_errbyfam(fam,errbypedf1,errbypedf2,errbyped3):
  '''
  Output genotype errors by family
  '''
  derrf1,derrf2,derrf3='','',''
  serrf1,serrf2,serrf3=0,0,0
  if fam in errbypedf1:
    derrf1=';'.join(errbypedf1[fam])
    serrf1=len(errbypedf1[fam])
  if fam in errbypedf2:
    derrf2=';'.join(errbypedf2[fam])
    serrf2=len(errbypedf2[fam])
  if fam in errbyped3:
    derrf3=';'.join(errbyped3[fam])
    serrf3=len(errbyped3[fam])
  rowd=[fam,'',derrf1,derrf2,derrf3]
  rows=[fam,'',serrf1,serrf2,serrf3]
  return rowd,rows


def emit_errbyind(ind,errbyped1,errbyped2):
  '''
  Output genotype errors by individual
  '''
  derr1,derr2='',''
  serr1,serr2=0,0
  if ind in errbyped1:
    derr1=';'.join(errbyped1[ind])
    serr1=len(errbyped1[ind])
  if ind in errbyped2:
    derr2=';'.join(errbyped2[ind])
    serr2=len(errbyped2[ind])
  rowd=['',ind,derr1,derr2,'']
  rows=['',ind,serr1,serr2,'']
  return rowd,rows


def output_file(filename,*headers):
  out = table_writer(filename,hyphen=sys.stdout)
  out.writerows(headers)
  return out


def main():
  parser  = option_parser()
  options = parser.parse_args()

  genos = load_genostream(options.genotypes,format=options.informat,
                          genorepr=options.ingenorepr,
                          genome=options.loci,
                          phenome=options.pedigree).as_ldat()

  if options.errdetails:
    outdetails = output_file(options.errdetails,['']+genos.samples)

  nfams,numoffams = build_nuclear_families(genos)

  errbyloc1,errbyloc2,errbyloc=genotype_elimination(nfams,genos,genos.samples,outdetails)
  errbylocf1=count_by_family(errbyloc1)
  errbylocf2=count_by_family(errbyloc2)

  errbyloc3=make_errbyloc3(errbyloc,errbylocf1,errbylocf2)

  # Output genotype errors by locus
  if options.locdet and options.locsum:
    outd = output_file(options.locdet,errbylochead1,errbylochead2)
    outs = output_file(options.locsum,errbylochead1,errbylochead2)
    for locus in sorted(set(chain(errbyloc1,errbyloc2,errbyloc3))):
      rowd,rows=emit_errbyloc(locus,errbyloc1,errbylocf1,errbyloc2,errbylocf2,errbyloc3)
      outd.writerow(rowd)
      outs.writerow(rows)

  # Invert dictionary from by locus to by family
  errbyped1  = invert_dict(errbyloc1,False)
  errbyped2  = invert_dict(errbyloc2,False)
  errbyped3  = invert_dict(errbyloc3)
  errbypedf1 = invert_dict(errbylocf1)
  errbypedf2 = invert_dict(errbylocf2)

  # Output genotype errors by pedigree
  if options.peddet and options.pedsum:
    outd = output_file(options.peddet,errbypedhead)
    outs = output_file(options.pedsum,errbypedhead)

    famdict={}
    for ind in sorted(set(chain(errbyped1,errbyped2))):
      fam=ind.split('_')[0]
      if fam not in famdict:
        famdict[fam]=1
        rowd,rows=emit_errbyfam(fam,errbypedf1,errbypedf2,errbyped3)
        outd.writerow(rowd)
        outs.writerow(rows)

      rowd,rows=emit_errbyind(ind,errbyped1,errbyped2)

      outd.writerow(rowd)
      outs.writerow(rows)


def _test():
  import doctest,geno_elim
  return doctest.testmod(geno_elim)


if __name__ == '__main__':
  _test()
  main()

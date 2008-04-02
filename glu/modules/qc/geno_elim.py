# -*- coding: utf-8 -*-
'''
File:          geno_elim.py

Authors:       Xiang Deng(dengx@mail.nih.gov)

Created:       Thr Aug  11 14:45:03 EDT 2006

Abstract:      implement genotype elimination algorithm

Compatibility: Python 2.5 and above

Requires:      glu

Revision:      $Id$
'''

__copyright__ = 'Copyright (c) 2008, BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.'
__license__   = 'See GLU license for terms by running: glu license'

import csv
import sys

from   itertools         import izip,chain
from   operator          import itemgetter

from   glu.lib.fileutils import autofile,hyphen
from   glu.lib.genolib   import load_genostream,snp


errbylochead1 = ['','LEVEL_1_ERRORS','','LEVEL_2_ERRORS','','LEVEL_3_ERRORS']
errbylochead2 = ['LOCUS','INDIVIDUAL','FAMILY','INDIVIDUAL','FAMILY','FAMILY']
errbypedhead  = ['FAMILY','INDIVIDUAL','LEVEL_1_ERRORS','LEVEL_2_ERRORS','LEVEL_3_ERRORS']


def option_parser():
  import optparse

  usage = 'usage: %prog [options] genofile'
  parser = optparse.OptionParser(usage=usage, add_help_option=False)

  parser.add_option('-g', '--genofile',    dest='genofile',                 metavar='FILE',
                    help="Genotype file")

  parser.add_option('-p', '--pedfile',     dest='pedfile',                  metavar='FILE',
                    help="Pedigree file")

  parser.add_option('-o', '--errdetails',  dest='errdetails',  default='-', metavar='FILE',
                    help="The output file containing genotype matrix after elimination, '-' for standard out")

  parser.add_option('-L', '--errbylocsum', dest='errbylocsum', default='-', metavar='FILE',
                    help="The output file containing summary of errors by locus, '-' for standard out")

  parser.add_option('-l', '--errbylocdet', dest='errbylocdet', default='-', metavar='FILE',
                    help="The output file containing details of errors by locus, '-' for standard out")

  parser.add_option('-F', '--errbypedsum', dest='errbypedsum', default='-', metavar='FILE',
                    help="The output file containing summary of errors by pedigree, '-' for standard out")

  parser.add_option('-f', '--errbypeddet', dest='errbypeddet', default='-', metavar='FILE',
                    help="The output file containing details of errors by pedigree, '-' for standard out")

  return parser


def mkgeno(a,b):
  return min(a,b)+max(a,b)


def build_children_genos(p1geno,p2geno):
  '''
  Build a set of possible children genotypes from a pair of parents genotypes

  @param p1geno: the genotype of one of the parents
  @type  p1geno: str
  @param p2geno: the genotype of the other parent
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
  Build nuclear familes based on a pedigree file
  @param     pedfile: the pedigree file
  @type      pedfile: str
  @param individuals: genotyped individuals
  @type  individuals: list
  @return           : the nuclear families and the total number of built families
  @rtype            : tuple

  >>> import StringIO
  >>> contents = "fam101\\tind01\\t0\\t0\\n\
  fam101\\tind02\\t0\\t0\\n\
  fam101\\tind11\\tind01\\tind02\\n\
  fam101\\tind12\\tind01\\tind02\\n"
  >>> f = StringIO.StringIO(contents)
  >>> ind = ['fam101_ind01','fam101_ind02','fam101_ind11','fam101_ind12']
  >>> build_nuclear_families(f,ind)
  ({('fam101_ind01', 'fam101_ind02'): ['fam101_ind11', 'fam101_ind12']}, 1)
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

  peds = csv.reader(autofile(pedfile),dialect='excel-tab')

  nfam = {}
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
    nfam.setdefault(parents, []).append(ind)
    famset.add(fam)

  return nfam,len(famset)


def build_genodict(genos,individuals):
  '''
  Make a possible genotype list for each individual
  @param       genos: genotypes
  @type        genos: list
  @param individuals: individuals
  @type  individuals: list
  @return           : the genotype dictionary for each individual
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

  @param       nfams: the nuclear families
  @type        nfams: dict
  @param       genos: genotypes
  @type        genos: list
  @param individuals: individuals
  @type  individuals: list
  @param         out: the output file
  @type          out: file
  '''

  def norm_geno(g):
    if not g:
      return ''
    else:
      return intern(''.join(g))

  errbyloc1 = {}
  errbyloc2 = {}
  errbyloc = {}
  ct=0
  for lname,locusgenos in genos:
    ct+=1
    print '------------------------line: %s' % ct
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

  return errbyloc1,errbyloc2,errbyloc,ct


def build_errbyloc(individuals, lname, genos, errbyloc):
  '''
  Output genotype errors by locus after elimination
  '''
  for ind in individuals:
    if not genos[ind]:
      errbyloc.setdefault(lname, set()).add(intern(ind.split('_')[0]))


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

  @param   parents: the ids for both parents
  @type    parents: tuple
  @param  children: the ids for all the children
  @type   children: list
  @param     lname: locus name
  @type      lname: str
  @param     genos: the current possible genotypes list for each person
  @type      genos: dict
  @param errbyloc1: the ids for individuals with level 1 errors by locus
  @type  errbyloc1: dict
  @param errbyloc2: the ids for individuals with level 2 errors by locus
  @type  errbyloc2: dict

  >>> parents = '1140_01','1140_02'
  >>> children = ['1140_11','1140_12']
  >>> lname = 'rs1167772'
  >>> genos = {'1140_01':set(['AA']),'1140_02':set(['AA','AG','GG']),'1140_11':set(['AG']),'1140_12':set(['AA','AG','GG'])}
  >>> errbyloc1={}
  >>> errbyloc2={}
  >>> eliminate_nuclear_family(parents,children,lname,genos,errbyloc1,errbyloc2)
  >>> genos == {'1140_02':set(['AG','GG']),'1140_01':set(['AA']),'1140_11':set(['AG']),'1140_12':set(['AA','AG'])}
  True
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
    return ( set(p) & set(c) )

  err = False
  for p in parents:
    for c in children:
      pgeno,cgeno=genos[p],genos[c]
      if not pair_check(first(pgeno),first(cgeno)):
        errbyloc1.setdefault(lname,set()).add( (p,c) )
        err = True
  return err


def level2_check(parents,children,lname,genos,errbyloc2):
  '''
  Find genotype errors at parent-parent-offspring level
  '''
  p1geno,p2geno=first(genos[parents[0]]),first(genos[parents[1]])
  for child in children:
    if not compatible_genos(p1geno,p2geno,[child],genos):
      errbyloc2.setdefault(lname,set()).add( (parents[0],parents[1],child) )


def eliminate_parents(p1,p1genos,p2,p2genos,children,genos,saved_genos):
  '''
  Eliminate superfluous genotypes of parent1 with respect to parent2 and
  their children.  Parent2 must be processed seperately, since compatibility
  is only throughly verified for parent1's genotypes.  Thus, this function
  may be called with parent1 and parent2 reversed and all saved genotypes of
  parent2 spared from consideration.

  @param           p1: the first parent
  @type            p1: str
  @param      p1genos: genotypes of the first parent
  @type       p1genos: set
  @param           p2: the second parent
  @type            p2: str
  @param      p2genos: genotypes of the second parent
  @type       p2genos: set
  @param     children: children of the nuclear family
  @type      children: list
  @param        genos: the current possible genotypes list for each person
  @type         genos: dict
  @param  saved_genos: the current saved genotypes set for each person
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

  @param       parents: the ids for both parents
  @type        parents: tuple
  @param      children: the ids for all the children
  @type       children: list
  @param         genos: the current possible genotypes set for each person
  @type          genos: dictionary
  @param   saved_genos: the current saved genotypes set for each person
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

  @param   p1geno: the first parent genotype
  @type    p1geno: str
  @param   p2geno: the second parent genotype
  @type    p2geno: str
  @param children: the list of children
  @type  children: list
  @param    genos: the individual to genotypes dictionary
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
  @parem genoset: Set of either parents or children of the fixed genotype
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


def first(seq):
  return iter(seq).next()


def eliminate_single_genotypes(parents,children,genos):
  '''
  Eliminate incompatible genotypes for ambiguous individuals based on the
  unambiguous family member genotypes with a single genotype

  @param  parents: the ids for both parents
  @type   parents: tuple
  @param children: the ids for all the children
  @type  children: list
  @param    genos: the genotype list for each person
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
        pgeno = first(pgenos)
        eliminate_inconsistency(pgeno,cgenos)

      # if child is fixed, check its genotypes against each parent
      elif len(cgenos) == 1:
        cgeno = first(cgenos)
        eliminate_inconsistency(cgeno,pgenos)


def invert_dict(errbyloc,simple=1):
  '''
  Invert a dict
  '''
  errbyped = {}
  for locus,fams in errbyloc.iteritems():
    if fams:
      for inds in fams:
        if simple:
          fam=first(inds)
          errbyped.setdefault(inds,set()).add(intern(locus))
        else:
          for ind in inds:
            errbyped.setdefault(ind,set()).add(intern(locus))
  return errbyped


def count_by_family(errdict):
  '''
  Convert error count by individual to by family
  '''
  newdict={}
  for key,valset in errdict.iteritems():
    for vals in valset:
      for val in vals:
        newdict.setdefault(key,set()).add(val.split('_')[0])
  return newdict


def concat_element(dataset):
  '''
  Format a set of tuples into a list
  '''
  datalist=[]
  for datatup in dataset:
    datalist.append(' '.join(datatup))
  return ';'.join(datalist)


def make_errbyloc3(errbyloc,errbylocf1,errbylocf2):
  '''
  Retrieve level 3 errors by subtracting level 1 and 2 errors from the total errors
  '''
  errbyloc3=errbyloc
  for key,vals in errbyloc3.iteritems():
    if key in errbylocf1:
      errbyloc3[key]=errbyloc3[key]-errbylocf1[key]
    if key in errbylocf2:
      errbyloc3[key]=errbyloc3[key]-errbylocf2[key]
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
  if errbypedf1.has_key(fam):
    derrf1=';'.join(errbypedf1[fam])
    serrf1=len(errbypedf1[fam])
  if errbypedf2.has_key(fam):
    derrf2=';'.join(errbypedf2[fam])
    serrf2=len(errbypedf2[fam])
  if errbyped3.has_key(fam):
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
  if errbyped1.has_key(ind):
    derr1=';'.join(errbyped1[ind])
    serr1=len(errbyped1[ind])
  if errbyped2.has_key(ind):
    derr2=';'.join(errbyped2[ind])
    serr2=len(errbyped2[ind])
  rowd=['',ind,derr1,derr2,'']
  rows=['',ind,serr1,serr2,'']
  return rowd,rows


def open_filehandle(file,*headers):
  out = autofile(hyphen(file,sys.stdout),'w')
  out = csv.writer(out,dialect='excel-tab')
  for header in headers:
    out.writerow(header)
  return out


def main():
  parser=option_parser()
  options,args=parser.parse_args()

  if not options.pedfile or not options.genofile:
    parser.print_help()
    return

  genos = iter(load_genostream(options.genofile,genorepr=snp).as_ldat())
  individuals = genos.next()

  if options.errdetails:
    outdetails = open_filehandle(options.errdetails,['']+individuals)

  nfams,numoffams = build_nuclear_families(options.pedfile,individuals)

  errbyloc1,errbyloc2,errbyloc,ctlocus=genotype_elimination(nfams,genos,individuals,outdetails)
  errbylocf1=count_by_family(errbyloc1)
  errbylocf2=count_by_family(errbyloc2)

  errbyloc3=make_errbyloc3(errbyloc,errbylocf1,errbylocf2)

  # Output genotype errorw by locus
  if options.errbylocdet and options.errbylocsum:
    outd = open_filehandle(options.errbylocdet,errbylochead1,errbylochead2)
    outs = open_filehandle(options.errbylocsum,errbylochead1,errbylochead2)
    for locus in sorted(set(chain(errbyloc1,errbyloc2,errbyloc3))):
      rowd,rows=emit_errbyloc(locus,errbyloc1,errbylocf1,errbyloc2,errbylocf2,errbyloc3)
      outd.writerow(rowd)
      outs.writerow(rows)

  # Invert dictionary from by locus to by family
  errbyped1=invert_dict(errbyloc1,0)
  errbyped2=invert_dict(errbyloc2,0)
  errbyped3=invert_dict(errbyloc3)
  errbypedf1=invert_dict(errbylocf1)
  errbypedf2=invert_dict(errbylocf2)

  # Output genotype errors by pedigree
  if options.errbypeddet and options.errbypedsum:
    outd = open_filehandle(options.errbypeddet,errbypedhead)
    outs = open_filehandle(options.errbypedsum,errbypedhead)

    famdict={}
    for ind in sorted(set(chain(errbyped1,errbyped2))):
      fam=ind.split('_')[0]
      if not famdict.has_key(fam):
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

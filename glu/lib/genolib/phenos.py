# -*- coding: utf-8 -*-
'''
File:          phenos.py

Authors:       Kevin Jacobs (jacobske@bioinformed.com)

Created:       2007-12-26

Abstract:      GLU phenotype model

Requires:      Python 2.5

Revision:      $Id$
'''

from collections               import defaultdict

from glu.lib.fileutils         import namefile,load_table,get_arg,tryint,parse_augmented_filename
from glu.lib.genolib.genoarray import model_from_alleles

UNKNOWN,MALE,FEMALE = range(3)

SEXMAP = {'':UNKNOWN, '?':UNKNOWN,
          '1':MALE,   'M':MALE,   'MALE':MALE,
          '2':FEMALE, 'F':FEMALE, 'FEMALE':FEMALE}

class Nothing(object): pass


# FIXME: Try various cases
def _get_header(filename,hmap,header,required=False):
  index = hmap.get(header,[])

  if not index:
    if required:
      raise ValueError('Cannot find %s header in locus file %s' % (header,namefile(filename)))
    return None

  elif len(index) != 1:
    raise ValueError('Ambiguous %s header in locus file %s' % (header,namefile(filename)))

  return index[0]


class Phenotypes(object):
  '''
  Locus metadata
  '''
  __slots__ = ('name','family','individual','parent1','parent2','phenoclass','sex')

  def __init__(self, name, family=None, individual=None, parent1=None, parent2=None,
                           phenoclass=None, sex=UNKNOWN):
    self.name       = name
    self.family     = family
    self.individual = individual
    self.parent1    = parent1
    self.parent2    = parent2
    self.phenoclass = phenoclass
    self.sex        = sex


class Phenome(object):
  '''
  Phenotype metadata storage

  FIXME: docstring
  '''

  def __init__(self, phenos=None):
    '''
    FIXME: docstring
    '''
    self.phenos = dict( (pheno.name,pheno) for pheno in phenos or [])

  def set_phenos(self, name, family=Nothing, individual=Nothing, parent1=Nothing, parent2=Nothing,
                             phenoclass=Nothing, sex=Nothing):
    '''
    FIXME: docstring
    '''
    pheno = self.phenos.get(name)
    if pheno is None:
      pheno = self.phenos[name] = Phenotypes(name)

    if family is not Nothing:
      pheno.family = family
    if individual is not Nothing:
      pheno.individual = individual
    if parent1 is not Nothing:
      pheno.parent1 = parent1
    if parent2 is not Nothing:
      pheno.parent2 = parent2
    if phenoclass is not Nothing:
      pheno.phenoclass = phenoclass
    if sex is not Nothing:
      pheno.sex = sex

  def merge_phenos(self, name, family=None, individual=None, parent1=None, parent2=None,
                              phenoclass=None, sex=UNKNOWN):
    '''
    FIXME: docstring
    '''
    pheno = self.phenos.get(name)

    # Fastpath: No merge needed
    if pheno is None:
      self.phenos[name] = Phenotypes(name,family,individual,parent1,parent2,phenoclass,sex)
      return

    # Slowpath: Must check for incompatibilities
    if family is not None:
      if pheno.family is None:
        pheno.family = family
      elif pheno.family is not family:
        raise ValueError('Incompatible family')

    if individual is not None:
      if pheno.individual is None:
        pheno.individual = individual
      elif pheno.individual is not individual:
        raise ValueError('Incompatible individual')

    if parent1 is not None:
      if pheno.parent1 is None:
        pheno.parent1 = parent1
      elif pheno.parent1 is not parent1:
        raise ValueError('Incompatible parent1')

    if parent2 is not None:
      if pheno.parent2 is None:
        pheno.parent2 = parent2
      elif pheno.parent2 is not parent2:
        raise ValueError('Incompatible parent2')

    if phenoclass is not None:
      if pheno.phenoclass is None:
        pheno.phenoclass = phenoclass
      elif pheno.phenoclass is not phenoclass:
        raise ValueError('Incompatible phenoclass')

    if sex is not UNKNOWN:
      if pheno.sex is UNKNOWN:
        pheno.sex = sex
      elif pheno.sex is not sex:
        raise ValueError('Incompatible sex')

  def get_phenos(self, name):
    '''
    FIXME: docstring
    '''
    pheno = self.phenos.get(name)

    if pheno is None:
      pheno = self.phenos[name] = Phenotypes(name)

    return pheno


def load_phenome_records(filename,extra_args=None,**kwargs):
  '''
  Load a pheno file

  @param     filename: file name or file object
  @type      filename: str or file object
  @type    extra_args: dict
  @return            : Phenome instance
  '''
  if extra_args is None:
    args = kwargs
  else:
    args = extra_args
    args.update(kwargs)

  rows   = load_table(filename,want_header=True,extra_args=args)
  header = rows.next()

  if extra_args is None and args:
    raise ValueError('Unexpected filename arguments: %s' % ','.join(sorted(args)))

  hmap = defaultdict(list)
  for i,h in enumerate(header):
    hmap[h].append(i)

  # FIXME: Add aliases
  ind_index     = _get_header(filename,hmap,'NAME',required=True)
  family_index  = _get_header(filename,hmap,'FAMILY')
  parent1_index = _get_header(filename,hmap,'PARENT1')
  parent2_index = _get_header(filename,hmap,'PARENT2')
  pheno_index   = _get_header(filename,hmap,'PHENO')
  sex_index     = _get_header(filename,hmap,'SEX')

  def _phenos():
    for i,row in enumerate(rows):
      if not row:
        continue

      n = len(row)

      if ind_index >= n:
        ind = ''
      else:
        ind = intern(row[ind_index].strip())

      if not ind:
        raise ValueError('Invalid phenotype record %d in %s (blank individual name)' % (i+1,namefile(filename)))

      family = individual = parent1 = parent2 = phenoclass = None
      sex = UNKNOWN

      if family_index is not None and family_index<n:
        family = intern(row[family].strip()) or None

      if parent1_index is not None and parent1_index<n:
        parent1 = intern(row[parent1].strip()) or None

      if parent2_index is not None and parent2_index<n:
        parent2 = intern(row[parent2].strip()) or None

      if pheno_index is not None and pheno_index<n:
        pheno = intern(row[pheno].strip()) or None

      if sex_index is not None and sex_index<n:
        sex = SEXMAP[row[sex].strip()]

      if family is not None:
        name = '%s:%s' % (ind,family)
      else:
        name = ind

      yield name,family,individual,parent1,parent2,phenoclass,sex

  return _pheno()


def load_phenome(filename,**kwargs):
  '''
  Return the default model and a sequence of Locus objects from an augmented
  locus description file

  FIXME: Add docstring and doctests
  '''
  phenos = load_locus_records(filename,**kwargs)

  phenome = Phenome()

  for p in phenos:
    phenome.merge_phenos(*p)

  return phenome


def _test():
  import doctest
  return doctest.testmod()


if __name__ == '__main__':
  _test()

# -*- coding: utf-8 -*-
'''
File:          phenos.py

Authors:       Kevin Jacobs (jacobske@bioinformed.com)

Created:       2007-12-26

Abstract:      GLU phenotype model

Requires:      Python 2.5

Revision:      $Id$
'''

__copyright__ = 'Copyright (c) 2008, BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.'
__license__   = 'See GLU license for terms by running: glu license'

from collections               import defaultdict

from glu.lib.fileutils         import namefile,load_table


SEX_UNKNOWN,SEX_MALE,SEX_FEMALE = None,0,1
PHENO_UNKNOWN,PHENO_UNAFFECTED,PHENO_AFFECTED = None,0,1

SEX_MAP = { '':SEX_UNKNOWN, '?':SEX_UNKNOWN,
           '1':SEX_MALE,    'M':SEX_MALE,     'MALE':SEX_MALE,
           '2':SEX_FEMALE,  'F':SEX_FEMALE, 'FEMALE':SEX_FEMALE}

PHENO_MAP = {  '':PHENO_UNKNOWN,    '?':PHENO_UNKNOWN,    'UNKNOWN':PHENO_UNKNOWN,
              '0':PHENO_UNAFFECTED, '-':PHENO_UNAFFECTED, 'UNAFFECTED':PHENO_UNAFFECTED,
              '1':PHENO_AFFECTED,   '+':PHENO_AFFECTED,   'AFFECTED':PHENO_AFFECTED}


class Nothing(object): pass


# FIXME: Try various cases, accept aliases
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
  __slots__ = ('name','family','individual','parent1','parent2','sex','phenoclass')

  def __init__(self, name, family=None, individual=None, parent1=None, parent2=None,
                           sex=SEX_UNKNOWN, phenoclass=None):
    self.name       = name
    self.family     = family
    self.individual = individual
    self.parent1    = parent1
    self.parent2    = parent2
    self.sex        = sex
    self.phenoclass = phenoclass


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
                             sex=Nothing, phenoclass=Nothing):
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
    if sex is not Nothing:
      pheno.sex = sex
    if phenoclass is not Nothing:
      pheno.phenoclass = phenoclass

  def merge_phenos(self, name, family=None, individual=None, parent1=None, parent2=None,
                               sex=SEX_UNKNOWN, phenoclass=None):
    '''
    FIXME: docstring
    '''
    pheno = self.phenos.get(name)

    # Fastpath: No merge needed
    if pheno is None:
      self.phenos[name] = Phenotypes(name,family,individual,parent1,parent2,sex,phenoclass)
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

    if sex is not SEX_UNKNOWN:
      if pheno.sex is SEX_UNKNOWN:
        pheno.sex = sex
      elif pheno.sex is not sex:
        raise ValueError('Incompatible sex')

    if phenoclass is not PHENO_UNKNOWN:
      if pheno.phenoclass is PHENO_UNKNOWN:
        pheno.phenoclass = phenoclass
      elif pheno.phenoclass is not phenoclass:
        raise ValueError('Incompatible phenoclass')

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

  # FIXME: Add aliases (MOTHER AND FATHER!!!)
  ind_index     = _get_header(filename,hmap,'NAME',required=True)
  family_index  = _get_header(filename,hmap,'FAMILY')
  parent1_index = _get_header(filename,hmap,'PARENT1')
  parent2_index = _get_header(filename,hmap,'PARENT2')
  sex_index     = _get_header(filename,hmap,'SEX')
  pheno_index   = _get_header(filename,hmap,'PHENO')

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

      family = individual = parent1 = parent2 = None
      sex    = SEX_UNKNOWN
      pheno  = PHENO_UNKNOWN

      if family_index is not None and family_index<n:
        family = intern(row[family].strip()) or None

      if parent1_index is not None and parent1_index<n:
        parent1 = intern(row[parent1_index].strip()) or None

      if parent2_index is not None and parent2_index<n:
        parent2 = intern(row[parent2_index].strip()) or None

      if sex_index is not None and sex_index<n:
        sex = SEX_MAP[row[sex_index].strip().upper()]

      if pheno_index is not None and pheno_index<n:
        pheno = PHENO_MAP.get(row[pheno_index].upper(),PHENO_UNKNOWN)

      if family is not None:
        name    = '%s:%s' % (family,ind)
        parent1 = '%s:%s' % (family,parent1) if parent1 else None
        parent2 = '%s:%s' % (family,parent2) if parent2 else None
      else:
        name = ind

      yield name,family,individual,parent1,parent2,sex,pheno

  return _phenos()


def load_phenome(filename,phenome=None,**kwargs):
  '''
  Return the default model and a sequence of Locus objects from an augmented
  locus description file

  FIXME: Add docstring and doctests
  '''
  phenos = load_phenome_records(filename,**kwargs)

  phenome = phenome or Phenome()

  for name,family,individual,parent1,parent2,sex,phenoclass in phenos:
    if parent1:
      phenome.merge_phenos(parent1)
    if parent2:
      phenome.merge_phenos(parent2)
    phenome.merge_phenos(name,family,individual,parent1,parent2,sex,phenoclass)

  return phenome


def merge_phenome_list(phenomes):
  phenome = phenomes[0]
  for ph in phenomes[1:]:
    for sample,pheno in ph.phenos.iteritems():
      if sample in phenome.phenos:
        phenome.merge_phenos(sample, file_phenos.family,  file_phenos.individual,
                                     file_phenos.parent1, file_phenos.parent2,
                                     file_phenos.sex,     file_phenos.phenoclass)
      else:
        phenome.phenos[sample] = pheno

  return phenome


def _test():
  import doctest
  return doctest.testmod()


if __name__ == '__main__':
  _test()

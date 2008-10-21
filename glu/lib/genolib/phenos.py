# -*- coding: utf-8 -*-

__abstract__  = 'GLU phenotype model'
__copyright__ = 'Copyright (c) 2008, BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.'
__license__   = 'See GLU license for terms by running: glu license'
__revision__  = '$Id$'


import sys


from collections               import defaultdict

from glu.lib.fileutils         import namefile,table_reader


SEX_UNKNOWN,SEX_MALE,SEX_FEMALE               = None,0,1
PHENO_UNKNOWN,PHENO_UNAFFECTED,PHENO_AFFECTED = None,0,1

SEX_MAP = { '':SEX_UNKNOWN, '?':SEX_UNKNOWN,
           '1':SEX_MALE,    'M':SEX_MALE,     'MALE':SEX_MALE,
           '2':SEX_FEMALE,  'F':SEX_FEMALE, 'FEMALE':SEX_FEMALE}

SEX_RMAP = {SEX_UNKNOWN:'?',SEX_MALE:'M',SEX_FEMALE:'F'}

PHENO_MAP = {  '':PHENO_UNKNOWN,    '?':PHENO_UNKNOWN,    'UNKNOWN':PHENO_UNKNOWN,
              '0':PHENO_UNAFFECTED, '-':PHENO_UNAFFECTED, 'UNAFFECTED':PHENO_UNAFFECTED,
              '1':PHENO_AFFECTED,   '+':PHENO_AFFECTED,   'AFFECTED':PHENO_AFFECTED}

PHENO_RMAP = {PHENO_UNKNOWN:'UNKNOWN',PHENO_UNAFFECTED:'UNAFFECTED',PHENO_AFFECTED:'AFFECTED'}


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

  def founder(self):
    return self.parent1 is None and self.parent2 is None

  def nonfounder(self):
    return self.parent1 is not None or self.parent2 is not None


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
                               sex=SEX_UNKNOWN, phenoclass=None, warn=False):
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
      elif pheno.family != family:
        msg = 'Phenotype record %s incompatible family (%s != %s)' % (name,pheno.family,family)
        if warn:
          sys.stderr.write('[WARNING] %s\n' % msg)
        else:
          raise ValueError(msg)

    if individual is not None:
      if pheno.individual is None:
        pheno.individual = individual
      elif pheno.individual != individual:
        msg = 'Phenotype record %s incompatible individual (%s != %s)' % (name,pheno.individual,individual)
        if warn:
          sys.stderr.write('[WARNING] %s\n' % msg)
        else:
          raise ValueError(msg)

    if parent1 is not None:
      if pheno.parent1 is None:
        pheno.parent1 = parent1
      elif pheno.parent1 != parent1:
        msg = 'Phenotype record %s incompatible parent1 (%s != %s)' % (name,pheno.parent1,parent1)
        if warn:
          sys.stderr.write('[WARNING] %s\n' % msg)
        else:
          raise ValueError(msg)

    if parent2 is not None:
      if pheno.parent2 is None:
        pheno.parent2 = parent2
      elif pheno.parent2 != parent2:
        msg = 'Phenotype record %s incompatible parent2 (%s != %s)' % (name,pheno.parent2,parent2)
        if warn:
          sys.stderr.write('[WARNING] %s\n' % msg)
        else:
          raise ValueError(msg)

    if sex is not SEX_UNKNOWN:
      if pheno.sex is SEX_UNKNOWN:
        pheno.sex = sex
      elif pheno.sex is not sex:
        msg = 'Phenotype record %s incompatible sex (%s != %s)' % (name,SEX_RMAP[pheno.sex],SEX_RMAP[sex])
        if warn:
          sys.stderr.write('[WARNING] %s\n' % msg)
        else:
          raise ValueError(msg)

    if phenoclass is not PHENO_UNKNOWN:
      if pheno.phenoclass is PHENO_UNKNOWN:
        pheno.phenoclass = phenoclass
      elif pheno.phenoclass is not phenoclass:
        msg = 'Phenotype record %s incompatible phenoclass (%s != %s)' \
                  % (name,PHENO_RMAP[pheno.phenoclass],PHENO_RMAP[phenoclass])
        if warn:
          sys.stderr.write('[WARNING] %s\n' % msg)
        else:
          raise ValueError(msg)

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

  rows   = table_reader(filename,want_header=True,extra_args=args)
  header = rows.next()

  if extra_args is None and args:
    raise ValueError('Unexpected filename arguments: %s' % ','.join(sorted(args)))

  hmap = defaultdict(list)
  for i,h in enumerate(header):
    hmap[h].append(i)

  name_index    = _get_header(filename,hmap,'NAME')
  ind_index     = _get_header(filename,hmap,'INDIVIDUAL')
  family_index  = _get_header(filename,hmap,'FAMILY')
  parent1_index = _get_header(filename,hmap,'PARENT1')
  parent2_index = _get_header(filename,hmap,'PARENT2')
  sex_index     = _get_header(filename,hmap,'SEX')
  pheno_index   = _get_header(filename,hmap,'PHENO')

  if name_index is None and ind_index is None:
    raise ValueError('Cannot find a NAME or INDIVIDUAL header in locus file %s' % namefile(filename))

  def _phenos():
    for i,row in enumerate(rows):
      if not row:
        continue

      n = len(row)

      name = intern(row[name_index].strip()) if name_index < n else ''
      ind  = intern(row[ind_index].strip())  if ind_index  < n else name

      if not name and not ind:
        raise ValueError('Invalid phenotype record %d in %s (blank individual name)' % (i+1,namefile(filename)))

      family = parent1 = parent2 = None
      sex    = SEX_UNKNOWN
      pheno  = PHENO_UNKNOWN

      if family_index is not None and family_index<n:
        family = intern(row[family_index].strip()) or None

      if parent1_index is not None and parent1_index<n:
        parent1 = intern(row[parent1_index].strip()) or None

      if parent2_index is not None and parent2_index<n:
        parent2 = intern(row[parent2_index].strip()) or None

      if sex_index is not None and sex_index<n:
        sex = SEX_MAP[row[sex_index].strip().upper()]

      if pheno_index is not None and pheno_index<n:
        pheno = PHENO_MAP.get(row[pheno_index].upper(),PHENO_UNKNOWN)

      if family is not None:
        ind     = '%s:%s' % (family,ind)
        parent1 = '%s:%s' % (family,parent1) if parent1 else None
        parent2 = '%s:%s' % (family,parent2) if parent2 else None

      if not name:
        name = ind

      yield name,family,ind,parent1,parent2,sex,pheno

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
        phenome.merge_phenos(sample, pheno.family,  pheno.individual,
                                     pheno.parent1, pheno.parent2,
                                     pheno.sex,     pheno.phenoclass)
      else:
        phenome.phenos[sample] = pheno

  return phenome


def _test():
  import doctest
  return doctest.testmod()


if __name__ == '__main__':
  _test()

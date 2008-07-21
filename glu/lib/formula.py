# -*- coding: utf-8 -*-

__abstract__  = 'Formula objects and parser for association testing'
__copyright__ = 'Copyright (c) 2008, BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.'
__license__   = 'See GLU license for terms by running: glu license'
__revision__  = '$Id$'

import re

from   ply    import lex,yacc

from   numpy  import exp,sqrt
from   scipy  import stats


try:
  import glu.lib.formula_parsetab as parsetab
except ImportError:
  parsetab = None


class TERM(object):
  def __init__(self, name):
    self.name  = name
    self.index = None

  def loci(self):
    if self.name is None:
      return []
    return [self.name]

  def terms(self):
    return [self]

  def expand_terms(self):
    return [self]

  def __mul__(self, other):
    t1 =  self.terms() if isinstance(self,  INTERACTION) else [self]
    t2 = other.terms() if isinstance(other, INTERACTION) else [other]
    return INTERACTION(t1 + t2)

  def __add__(self, other):
    t1 =  self.terms() if isinstance(self,  COMBINATION) else [self]
    t2 = other.terms() if isinstance(other, COMBINATION) else [other]
    return COMBINATION(t1 + t2)

  def formula(self, loci=None):
    return ' + '.join(self.term_names(loci))

  def effect_names(self, loci=None):
    return self.term_names(loci)

  def signature(self):
    return (type(self),self.name)

  def __hash__(self):
    return hash(self.signature())

  def __eq__(self,other):
    return self.signature() == other.signature()

  def __contains__(self, other):
    selfitems  = set(t.signature() for t in  self.terms())
    otheritems = set(t.signature() for t in other.terms())
    return otheritems.issubset(selfitems)

  def find(self, other):
    selfmap = dict( (t.signature(),t) for t in self.terms() )
    results = [ selfmap[t.signature()] for t in other.terms() ]

    # FIXME: wait until COMBINATION's signature is more clever
    #if not results:
    #  return NULL()
    #elif len(results) == 1:
    #  return results[0]
    #else:
    return COMBINATION(results)

  def odds_ratios(self,p):
    return exp(self.estimates(p))

  def se(self,c):
    return sqrt(self.var(c))

  def estimate_ci(self,p,c,alpha=0.95):
    a = stats.distributions.norm.ppf( (1+alpha)/2 )
    return [ (p-a*e,p+a*e) for p,e in
                        zip(self.estimates(p),self.se(c)) ]

  def odds_ratio_ci(self,p,c,alpha=0.95):
    return exp(self.estimate_ci(p,c,alpha=alpha))


ident_re = re.compile('[a-zA-Z_][a-zA-Z_:.0-9]*')

def quote_ident(name):
  if not ident_re.match(name):
    return "'%s'" % name
  else:
    return name


class PHENOTERM(TERM):
  def effects(self, loci, phenos, i):
    return [phenos[self.pindex]]

  def term_names(self, loci=None):
    return [quote_ident(self.name)]

  def __len__(self):
    return 1

  def loci(self):
    return []

  def indices(self):
    return [self.index]

  def estimates(self,p):
    return [p[self.index,0]]

  def var(self,c):
    i = self.index
    return [c[i,i]]


class INTERCEPT(TERM):
  def __init__(self):
    self.name = None

  def effects(self, loci, phenos, i):
    return [1]

  def term_names(self, loci=None):
    return ['_intercept']

  def __len__(self):
    return 1

  def indices(self):
    return [self.index]

  def estimates(self,p):
    return [p[self.index,0]]

  def var(self,c):
    i = self.index
    return [c[i,i]]


class NO_INTERCEPT(TERM):
  def __init__(self):
    self.name = None

  def effects(self, loci, phenos, i):
    return []

  def term_names(self, loci=None):
    return []

  def __len__(self):
    return 0


class GENOTERM(TERM):
  def effect_names(self, loci=None):
    if loci is None:
      return [quote_ident('%s:het' % self.name), quote_ident('%s:hom' % self.name)]

    lmodel = loci[self.name]
    if lmodel.genocount < 2:
      return []
    return [ quote_ident('%s:%s%s' % (self.name,lmodel.alleles[0],lmodel.alleles[1])),
             quote_ident('%s:%s%s' % (self.name,lmodel.alleles[1],lmodel.alleles[1])) ]


class NULL(GENOTERM):
  def __init__(self, name=None):
    self.name = name

  def effects(self, loci, phenos, i):
    if self.name is None:
      return []

    lmodel = loci[self.name]
    if lmodel.genos[i] not in lmodel.genomap:
      return [None]
    return []

  def term_names(self, loci=None):
    return []

  def __len__(self):
    return 0

  def indices(self):
    return []


class GENO(GENOTERM):
  def effects(self, loci, phenos, i):
    lmodel  = loci[self.name]
    genomap = lmodel.genomap
    geno    = lmodel.genos[i]

    if geno not in genomap:
      return [None]

    genos = [0,0]
    j = genomap[geno]
    if j:
      genos[j-1] = 1

    return genos

  def term_names(self, loci=None):
    if loci is None:
      return [quote_ident('%s:het' % self.name), quote_ident('%s:hom' % self.name)]

    lmodel = loci[self.name]
    if lmodel.genocount < 3:
      return []
    return [ quote_ident('%s:%s' % (self.name,''.join(g))) for g in lmodel.tests[1:] ]

  def indices(self):
    return [self.index,self.index+1]

  def estimates(self,p):
    return p[self.index:self.index+2,0].A.flatten()

  def var(self,c):
    i = self.index
    return [c[i,i],c[i+1,i+1]]

  def __len__(self):
    return 2


class ADOM(GENOTERM):
  def effects(self, loci, phenos, i):
    lmodel  = loci[self.name]
    genomap = loci[self.name].genomap
    geno    = lmodel.genos[i]
    if geno not in lmodel.genomap:
      return [None]
    return ([0,0],[1,1],[2,0])[lmodel.genomap[geno]]

  def term_names(self, loci=None):
    if loci is None:
      return [quote_ident('%s:trend' % self.name), quote_ident('%s:domdev' % self.name)]

    lmodel = loci[self.name]
    if lmodel.genocount < 3:
      return []
    return [ quote_ident('%s:trend:%s'     % (self.name,lmodel.alleles[1])),
             quote_ident('%s:domdev:%s%s'  % (self.name,lmodel.alleles[0],lmodel.alleles[1])) ]

  def indices(self):
    return [self.index,self.index+1]

  def estimates(self,p):
    return p[self.index:self.index+2,0].A.flatten()

  def var(self,c):
    i = self.index
    return [c[i,i],c[i+1,i+1]]

  def __len__(self):
    return 2


class TREND(GENOTERM):
  def effects(self, loci, phenos, i):
    lmodel = loci[self.name]
    return [ lmodel.genomap.get(lmodel.genos[i]) ]

  def term_names(self, loci=None):
    if loci is None:
      return [quote_ident('%s:trend' % self.name)]
    lmodel = loci[self.name]
    if lmodel.genocount < 2:
      return []
    return [quote_ident('%s:trend:%s' % (self.name,lmodel.alleles[1]))]

  def indices(self):
    return [self.index]

  def estimates(self,p):
    return [p[self.index,0],p[self.index,0]*2]

  def var(self,c):
    i = self.index
    return [c[i,i],4*c[i,i]]

  def __len__(self):
    return 1


class DOM(GENOTERM):
  def effects(self, loci, phenos, i):
    lmodel  = loci[self.name]
    genomap = lmodel.genomap
    geno    = lmodel.genos[i]
    if geno not in genomap:
      return [None]
    return [ int(genomap.get(geno)>0) ]

  def term_names(self, loci=None):
    if loci is None:
      return [quote_ident('%s:dom' % self.name)]
    lmodel = loci[self.name]
    if lmodel.genocount < 2:
      return []
    return [quote_ident('%s:dom:%s' % (self.name,lmodel.alleles[1])) ]

  def indices(self):
    return [self.index]

  def estimates(self,p):
    return [p[self.index,0],p[self.index,0]]

  def var(self,c):
    i = self.index
    return [c[i,i],c[i,i]]

  def __len__(self):
    return 1


class REC(GENOTERM):
  def effects(self, loci, phenos, i):
    lmodel  = loci[self.name]
    genomap = lmodel.genomap
    geno    = lmodel.genos[i]
    if geno not in genomap:
      return [None]
    return [ int(genomap.get(geno)>1) ]

  def term_names(self, loci=None):
    if loci is None:
      return [quote_ident('%s:rec' % self.name)]
    lmodel = loci[self.name]
    if lmodel.genocount < 2:
      return []
    return [quote_ident('%s:rec:%s%s' % (self.name,lmodel.alleles[1],lmodel.alleles[1])) ]

  def indices(self):
    return [self.index]

  def estimates(self,p):
    return [0,p[self.index,0]]

  def var(self,c):
    i = self.index
    return [0,c[i,i]]

  def __len__(self):
    return 1


class MISSING(GENOTERM):
  def effects(self, loci, phenos, i):
    if loci[self.name].genos[i]:
      return [1]
    else:
      return [0]

  def term_names(self, loci=None):
    if loci is not None and not loci[self.name].genocount:
      return []
    return [quote_ident('%s:missing' % self.name)]

  def indices(self):
    return [self.index]

  def estimates(self,p):
    return [p[self.index,0]]

  def var(self,c):
    i = self.index
    return [c[i,i]]

  def __len__(self):
    return 1


class NOT_MISSING(GENOTERM):
  def effects(self, loci, phenos, i):
    if loci[self.name].genos[i]:
      return [0]
    else:
      return [1]

  def term_names(self, loci=None):
    if self.loci is not None and not loci[self.name].genocount:
      return []
    return [quote_ident('%s:not_missing' % self.name)]

  def indices(self):
    return [self.index]

  def estimates(self,p):
    return [p[self.index,0]]

  def var(self,c):
    i = self.index
    return [c[i,i]]

  def __len__(self):
    return 1


class COMPOUNDTERM(TERM):
  def __init__(self, terms=None):
    terms = terms or []
    if not isinstance(terms,list):
      terms = list(terms)

    self.subterms = terms

  def loci(self):
    results = set()
    for term in self.subterms:
      results.update(term.loci())
    return sorted(results)

  def expand_terms(self):
    results = []
    for term in self.subterms:
      results.extend(term.expand_terms())
    return results

  def signature(self):
    sigs = tuple(sorted(t.signature() for t in self.subterms))
    return (type(self),) + sigs

  def __eq__(self,other):
    return self.signature() == other.signature()


class INTERACTION(COMPOUNDTERM):
  def indices(self):
    return range(self.index,self.index+len(self))

  def estimates(self,p):
    return p[self.index:self.index+len(self),0].A.flatten()

  def var(self,c):
    return [c[i,i] for i in self.indices()]

  def effects(self, loci, phenos, i):
    results = []
    for term in self.subterms:
      effects = term.effects(loci,phenos,i)
      if None in effects:
        return [None]
      if not results:
        results = list(effects)
      else:
        results = [ t1*t2 for t1 in results for t2 in effects ]
    return results

  def term_names(self, loci=None):
    results = []
    for term in self.subterms:
      names = term.term_names(loci)

      if not results:
        results.extend(names)
      else:
        results = [ '%s*%s' % (t1,t2) for t1 in results for t2 in names ]

    return results

  def effect_names(self, loci=None):
    results = []
    for term in self.subterms:
      names = term.effect_names(loci)

      if not results:
        results.extend(names)
      else:
        results = [ '%s*%s' % (t1,t2) for t1 in results for t2 in names ]

    return results

  def __len__(self):
    if not self.subterms:
      return 0

    l = 1
    for term in self.subterms:
      l *= len(term)
    return l


class COMBINATION(COMPOUNDTERM):
  def terms(self):
    return self.subterms

  def indices(self):
    results = []
    for term in self.subterms:
      results.extend(term.indices())
    return results

  def effects(self, loci, phenos, i):
    results = []
    for term in self.subterms:
      effects = term.effects(loci,phenos,i)
      if None in effects:
        return [None]
      results.extend(effects)
    return results

  def term_names(self, loci=None):
    results = []
    for term in self.subterms:
      results.extend(term.term_names(loci))
    return results

  def effect_names(self, loci=None):
    results = []
    for term in self.subterms:
      results.extend(term.effect_names(loci))
    return results

  def estimates(self,p):
    results = []
    for term in self.subterms:
      results.extend(term.estimates(p))
    return results

  def var(self,c):
    results = []
    for term in self.subterms:
      results.extend(term.var(c))
    return results

  def __len__(self):
    return sum(len(term) for term in self.subterms)


termmap = { 'GENO'           : GENO,
            'GENOTYPE'       : GENO,
            'ADDDOM'         : ADOM,
            'ADOM'           : ADOM,
            'TREND'          : TREND,
            'MULTIPLICATIVE' : TREND,
            'MULT'           : TREND,
            'DOMINANT'       : DOM,
            'DOM'            : DOM,
            'RECESSIVE'      : REC,
            'REC'            : REC,
            'MISSING'        : MISSING,
            'MISS'           : MISSING,
            'NOT_MISSING'    : NOT_MISSING,
            'NOT_MISS'       : NOT_MISSING,
            'NULL'           : NULL,
          }

def get_term(name):
  try:
    return termmap[name.upper()]
  except KeyError:
    raise KeyError('Unknown genetic model term: %s' % name)


#################################################################


class FormulaParser(object):
  def __init__(self, parsetab=parsetab):
    self.lexer  = lex.lex(module=self)
    self.parser = yacc.yacc(module=self,tabmodule=parsetab,start='formula')


  def parse(self, s):
    return self.parser.parse(s,lexer=self.lexer)


  tokens = ('ONE', 'ZERO',
            'PLUS','TIMES','EQUALS',
            'LPAREN','RPAREN', 'TERM', 'IDENT', 'QUOTED')

  # Tokens
  t_ONE     = r'1'
  t_ZERO    = r'0'
  t_PLUS    = r'\+'
  t_TIMES   = r'\*'
  t_EQUALS  = r'='
  t_LPAREN  = r'\('
  t_RPAREN  = r'\)'

  def t_IDENT(self, t):
    r'[a-zA-Z_][a-zA-Z_:.0-9]*'
    if t.value in termmap:
      t.type = 'TERM'
    return t

  def t_QUOTED(self, t):
    r'\'([^\\\n]|(\\.))*?\'|"([^\\\n]|(\\.))*?"'
    t.value = t.value[1:-1]
    if t.value in termmap:
      t.type = 'TERM'
    return t

  t_ignore = ' \t\r\n'

  def t_error(self, t):
    raise ValueError("Illegal character '%s'" % t.value[0])

  # Parsing rules

  precedence = [ ('left','PLUS' ),
                 ('left','TIMES') ]

  def p_formula(self, t):
    '''
    formula : name EQUALS expression
    '''
    t[0] = t[1],t[3]

  def p_formula2(self, t):
    '''
    formula : expression
    '''
    t[0] = None,t[1]

  def p_name(self, t):
    '''
    name : IDENT
         | QUOTED
    '''
    t[0] = t[1]

  def p_expression_binop(self, t):
    '''
    expression : expression PLUS  expression
               | expression TIMES expression
    '''
    if   t[2] == '+': t[0] = t[1] + t[3]
    elif t[2] == '*': t[0] = t[1] * t[3]

  def p_expression_one(self, t):
    '''
    expression : ONE
    '''
    t[0] = INTERCEPT()

  def p_expression_zero(self, t):
    '''
    expression : ZERO
    '''
    t[0] = NO_INTERCEPT()

  def p_expression_paren(self, t):
    '''
    expression : LPAREN expression RPAREN
    '''
    t[0] = t[2]

  def p_expression_var(self, t):
    '''
    expression : name
    '''
    t[0] = PHENOTERM(t[1])

  def p_expression_func(self, t):
    '''
    expression : TERM LPAREN name RPAREN
    '''
    t[0] = termmap[t[1]](t[3])

  def p_error(self, t):
    raise ValueError("Illegal token '%s'" % t)


def main():
  parser = FormulaParser(parsetab='formula_parsetab')

  while 1:
    try:
      s = raw_input('formula> ')
    except EOFError:
      break

    p,t = parser.parse(s)
    print t.formula()


if __name__=='__main__':
  main()

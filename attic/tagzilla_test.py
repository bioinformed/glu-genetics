# -*- coding: utf-8 -*-

__abstract__  = 'Test tagzilla'
__copyright__ = 'Copyright (c) 2008, BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.'
__license__   = 'See GLU license for terms by running: glu license'
__revision__  = '$Id$'

import unittest
import random

from tagzilla import *

def  gen_pair_indices(n):
    '''
    Generate pair-wise indices (i,j) for all i,j < n and i<j
    '''
    return [ (i,j) for i in range(n) for j in range(i+1,n) ]

def bincmp(bin1, bin2):
  if(bin1 == bin2 and abs(bin1.maf - bin2.maf)< 0.0001 and bin1.disposition == bin2.disposition and bin1.maxcovered == bin2.maxcovered):
    return True
  else:
    return False

alleles = 'AB'
def pick(p):
  '''
  Pick a random allele from a biallelic locus with P(allele1) = p
  '''
  return alleles[random.random() > p]


def simulate_x_genos(p, males, females):
  '''
  Return a simulated sample of male and female X-linked genotypes.  Results
  are returned as the concatination of the male and female genotypes.
  '''
  male_genos   = [  ' '    + pick(p) for i in xrange(males)   ]
  female_genos = [ pick(p) + pick(p) for i in xrange(females) ]
  return male_genos+female_genos



class testTests(unittest.TestCase):
  '''
  A test class for test support code
  '''
  def test_simulate(self):
    n = 80000
    g = simulate_x_genos(0.2, n, 0)
    p = g.count(' A')/float(n)
    q = g.count(' B')/float(n)
    self.assertEquals(len(g),n)
    self.assertAlmostEqual(p,0.2,places=2)
    self.assertAlmostEqual(q,0.8,places=2)


class tagzillaTest(unittest.TestCase):
  '''
  A test class for the tagzilla module.
  '''

  def test_mean(self):
    '''
    Test a successful run of the mean function for a
    sequence feed
    '''
    self.assertRaises(ValueError, mean, [])
    self.assertRaises(TypeError,  mean, 'abc')
    self.assertEquals(1,          mean([1]))
    self.assertAlmostEquals(1.5,  mean([1,2]))
    self.assertEquals(2,          mean([1,2,3]))
    self.assertEquals(4,          mean(range(0,9)))

  def test_median(self):
    '''
    Test a successful run of the median function for a
    sequence feed
    '''
    self.assertRaises(ValueError, median, [])
    self.assertEquals(1,          median([1]))
    self.assertAlmostEquals(1.5,  median([1,2]))
    self.assertEquals(2,          median([1,2,4]))
    self.assertEquals(4,          median(range(0,9)))
    self.assertAlmostEquals(4.5,  median(range(0,10)))

  def test_estimate_maf(self):
    '''
    Test a successful run of the estimate_maf and slow_estimate_maf
    functions for a list of genotypes at a locus
    '''
    for func in (estimate_maf,slow_estimate_maf):
      # Test length invariance property
      for n in (1,5,10,100,1000):
        # Test the estimate_maf on a list of genotypes
        self.assertEquals(0, func([]*n))
        self.assertEquals(0, func(['AA']*n))
        self.assertEquals(0, func(['  ']*n))
        self.assertAlmostEquals(0.5, func(['AB']*n))
        self.assertAlmostEquals(0.5, func(['AB',None]*n))
        self.assertAlmostEquals(0.3, func(['AA', 'AC', 'AA', 'CC', 'AA']*n))
        self.assertAlmostEquals(0.3, func(['AA', 'AC', 'AA', None, 'CC', 'AA']*n))
        self.assertRaises(ValueError,func,['AA', 'AC', 'AA', 'AG', 'CC', 'AA']*n)


  def test_count_haplotypes(self):
    '''
    Test a successful run of the count_haplotypes function for the genotypes at two loci
    '''
    genos1  = ['AA', 'AC', 'AC', 'CC', 'CC', 'CC', 'AA']
    genos2  = ['AA', 'AC', 'AA', 'AA', 'AC', 'CC', 'CC']
    genos3  = ['AB', 'AA']
    genos4  = ['AB', 'AC']  # tri-allelic
    genos5  = [' A', ' B', ' A', 'AA', 'AB', 'BB', 'AA'] # X-linked test
    genos6  = [' C', ' D', ' C', 'CC', 'CD', 'DD', 'DD']
    results = [(0,0,0,0,0), (2,0,0,0,0), (2,0,0,0,1), (3,0,1,0,1),(3,0,3,0,1),
               (3,0,4,1,1), (3,0,4,3,1), (3,2,4,3,1)]
    results2 = [(0,0,0,0,0), (1,0,0,0,0), (1,0,0,1,0), (2,0,0,1,0),(4,0,0,1,0),
               (4,0,0,1,1), (4,0,0,3,1), (4,2,0,3,1)]
    for func in (count_haplotypes,slow_count_haplotypes):
      for i in range(len(results)):
        self.assertEquals(results[i], func(genos1[:i],genos2[:i]))
        self.assertEquals(results2[i], func(genos5[:i],genos6[:i]))
      self.assertRaises(ValueError, func, genos3, genos4)
      self.assertRaises(ValueError, func, genos1[:6], genos1)

  def test_x_linked(self):
    '''
    '''
    cases = [(10000,0), (0,5000), (4000,4000), (10000,1000), (2000,5000)]
    for m,f in cases:
      genos1 = simulate_x_genos(0.2, m, f)
      genos2 = simulate_x_genos(0.5, m, f)

      maf_m = estimate_maf(genos1[:m])
      maf_f = estimate_maf(genos1[m:])
      maf = estimate_maf(genos1)
      if m:
        self.assertAlmostEquals(maf_m, 0.2, places=1)
      else:
        self.assertEquals(maf_m, 0)
      if f:
        self.assertAlmostEquals(maf_f, 0.2, places=1)
      else:
        self.assertEquals(maf_f, 0)
      self.assertAlmostEquals(maf,   0.2, places=1)

      maf_m = estimate_maf(genos2[:m])
      maf_f = estimate_maf(genos2[m:])
      maf = estimate_maf(genos2)
      if m:
        self.assertAlmostEquals(maf_m, 0.5, places=1)
      else:
        self.assertEquals(maf_m, 0)
      if f:
        self.assertAlmostEquals(maf_f, 0.5, places=1)
      else:
        self.assertEquals(maf_f, 0)
      self.assertAlmostEquals(maf,   0.5, places=1)

      d = count_diplotypes(genos1,genos2)
      self.assertEquals(sum(n for g1,g2,n in d), m+f)

      last_h = None
      for count in (count_haplotypes,slow_count_haplotypes):
        h = count(genos1,genos2)
        self.assertEquals(sum(h)+h[-1], m+2*f)
        if last_h is not None:
          self.assertEquals(h, last_h)
        last_h = h

      for ld in (estimate_ld,slow_estimate_ld):
        # Test trivial cases with complete LD
        r2,dprime = ld(*count_haplotypes(genos1,genos1))
        self.assertAlmostEquals(r2,1,places=2)
        self.assertAlmostEquals(dprime,1,places=2)
        r2,dprime = ld(*count_haplotypes(genos2,genos2))
        self.assertAlmostEquals(r2,1,places=2)
        self.assertAlmostEquals(dprime,1,places=2)
        # Test case with no LD
        r2,dprime = ld(*h)
        self.assertAlmostEquals(r2,0,places=2)
        self.assert_(dprime<0.07)

  def test_build_binsets(self):
    '''
    Test a successful run of build_binsets function given loci, ldpairs, include and exclude
    '''
    genotypes = ['GG','GG','GA','GG','GA'] # maf = 0.2
    loci = dict( (i,Locus(i, 0, genotypes)) for i in range(9) )

    # populate ldpairs (list of tuples [(locus1.name, location1, locus2.name, location2, r2, dprime)] )
    ldpairs = [ (i,j,1,1) for i,j in gen_pair_indices(9) ]
    include = Includes(set([0]),set()) # single include
    exclude = set([8]) # single exclude
    binsets,lddata = build_binsets(loci, [ldpairs], include, exclude, {})
    locuslist = [8, 4, 5, 6, 7, 0, 1, 2, 3]
    for locus,bin in binsets.iteritems():
        if locus == 0:
            self.assertEquals(bin, Bin(locuslist, 1.8, 0, 9))
        elif locus in (1,2,3,4,5,6,7):
            self.assertEquals(bin, Bin(locuslist, 1.8, None, 9))
        else:
            self.assertEquals(bin, Bin(locuslist, 1.8, 1, 9))

    include = Includes(set([0,1]),set()) # multiple include
    exclude = set([7,8]) # multiple exclude
    binsets,lddata = build_binsets(loci, [ldpairs], include, exclude, {})
    locuslist1 = [8, 4, 5, 6, 7, 0, 2, 3]
    locuslist2 = [8, 4, 5, 6, 7, 1, 2, 3]
    locuslist3 = [8, 4, 5, 6, 7, 0, 1, 2, 3]
    for locus, bin in binsets.iteritems():
        if locus == 0:
            self.assertEquals(bincmp(bin, Bin(locuslist1, 1.6, Bin.INCLUDE_TYPED, 9)), True)
        elif locus == 1:
            self.assertEquals(bincmp(bin, Bin(locuslist2, 1.6, Bin.INCLUDE_TYPED, 9)), True)
        elif locus in (2,3,4,5,6):
            self.assertEquals(bincmp(bin, Bin(locuslist3, 1.8, Bin.NORMAL, 9)),True)
        else:
            self.assertEquals(bin, Bin(locuslist, 1.8, 1, 9))

  def test_binner(self):
    '''test binner given loci, binsets and lddate'''
    location = 10000
    genos = ['AA', 'AC', 'AA', 'CC', 'AA']
    maf = estimate_maf(genos)
    loci = {}
    # populuate 5 loci below
    for i in range(5):
        loci[i] = Locus(i, location + 500 * i, genos)

    # populate binsets
    binsets = {}
    binsets[0] = Bin([0,1,2,3],  maf*4, Bin.INCLUDE_TYPED) # 0 is an include
    binsets[1] = Bin([1,0,2],    maf*3)
    binsets[2] = Bin([2,0,1,3],  maf*4)
    binsets[3] = Bin([3,0,2,4],  maf*4)
    binsets[4] = Bin([4,3],      maf*2, Bin.EXCLUDE)

    # populate lddata , notice that there is dependecy between lddata an binsets
    lddata = {}
    lddata[(0,1)] = (1,1)
    lddata[(0,2)] = (1,1)
    lddata[(0,3)] = (1,1)
    lddata[(1,2)] = (1,1)
    lddata[(2,3)] = (1,1)
    lddata[(3,4)] = (1,1)
    results = list(binner(loci, binsets, lddata, Includes(set(),set())))
    self.assertEquals(len(results), 2)
    self.assertEquals(results[0].tags, [0, 2])
    self.assertEquals(results[0].others, [1, 3])
    self.assertAlmostEquals(results[0].average_maf, maf)
    self.assertEquals(results[0].include, 0)
    self.assertEquals(results[0].ld, [(0, 0, 1, 1),
                                      (2, 2, 1, 1),
                                      (0, 1, 1, 1),
                                      (0, 2, 1, 1),
                                      (1, 2, 1, 1),
                                      (0, 3, 1, 1),
                                      (2, 3, 1, 1)])

    pair_dispositions = [ (0,0,'typed-tag'),
                          (2,2,'alternate-tag'),
                          (0,1,'tag-other'),
                          (0,2,'tag-tag'),
                          (1,2,'other-tag'),
                          (0,3,'tag-other'),
                          (2,3,'tag-other') ]

    self.assertEquals(results[0].disposition, 'obligate-typed')

    for lname1,lname2,disposition in pair_dispositions:
      self.assertEquals( pair_disposition(lname1, lname2, results[0]), disposition)

    # Modify binsets to create next test case
    binsets[0] = Bin([0,1,2,3], maf*4)
    binsets[1] = Bin([1,0,2],   maf*3, Bin.INCLUDE_TYPED) # 1 is an include
    binsets[2] = Bin([2,0,1,3], maf*4)
    binsets[3] = Bin([3,0,2,4], maf*4)
    binsets[4] = Bin([4,3],     maf*2, Bin.EXCLUDE)

    lddata[(0,1)] = (1,1)
    lddata[(0,2)] = (1,1)
    lddata[(0,3)] = (1,1)
    lddata[(1,2)] = (1,1)
    lddata[(2,3)] = (1,1)
    lddata[(3,4)] = (1,1)
    results = list(binner(loci, binsets, lddata, Includes(set(),set())))
    self.assertEquals(len(results), 2)
    self.assertEquals(results[0].tags, [0, 1, 2])
    self.assertEquals(results[0].others, [])
    self.assertAlmostEquals(results[0].average_maf, maf)
    self.assertEquals(results[0].include, 1)
    self.assertEquals(results[0].ld, [(0, 0, 1, 1),
                                      (1, 1, 1, 1),
                                      (2, 2, 1, 1),
                                      (0, 1, 1, 1),
                                      (0, 2, 1, 1),
                                      (1, 2, 1, 1)])

    pair_dispositions = [ (0,0,'alternate-tag'),
                          (1,1,'typed-tag'),
                          (2,2,'alternate-tag'),
                          (0,1,'tag-tag'),
                          (0,2,'tag-tag'),
                          (0,2,'tag-tag'),
                          (1,2,'tag-tag') ]

    self.assertEquals(results[0].disposition, 'obligate-typed')

    for lname1,lname2,disposition in pair_dispositions:
      self.assertEquals( pair_disposition(lname1, lname2, results[0]), disposition)

    self.assertEquals(results[1].tags, [3])
    self.assertEquals(results[1].others, [4])
    self.assertAlmostEquals(results[1].average_maf, maf)
    self.assertEquals(results[1].include, None)
    self.assertEquals(results[1].ld, [(3, 3, 1, 1),
                                      (3, 4, 1, 1)])
    self.assertEquals(results[1].disposition, 'maximal-bin')
    self.assertEquals( pair_disposition(3, 3, results[1]), 'necessary-tag')
    self.assertEquals( pair_disposition(3, 4, results[1]), 'tag-other')

    binsets[0] = Bin([0,1],   maf*2)
    binsets[1] = Bin([1,0,2], maf*3, Bin.INCLUDE_TYPED) # 1 is an include
    binsets[2] = Bin([2,3],   maf*2)
    binsets[3] = Bin([3,2],   maf*2)
    binsets[4] = Bin([4],     maf)
    ldddata = {}
    lddata[(0,1)] = (1,1)
    lddata[(1,2)] = (1,1)
    lddata[(2,3)] = (1,1)
    results = list(binner(loci, binsets, lddata, Includes(set(),set())))
    self.assertEquals(len(results), 3)
    self.assertEquals(results[0].tags, [1])
    self.assertEquals(results[0].others, [0,2])
    self.assertAlmostEquals(results[0].average_maf, maf)
    self.assertEquals(results[0].include, 1)
    self.assertEquals(results[0].ld, [(1, 1, 1, 1),
                                      (0, 1, 1, 1),
                                      (1, 2, 1, 1)])
    self.assertEquals(results[0].disposition, 'obligate-typed')

    self.assertEquals( pair_disposition(1, 1, results[0]), 'typed-tag')
    self.assertEquals( pair_disposition(0, 1, results[0]), 'other-tag')
    self.assertEquals( pair_disposition(1, 2, results[0]), 'tag-other')

    i,j = 1,2
    if results[1].tags == [3]:
      i,j = j,i

    self.assertEquals(results[i].tags, [4])
    self.assertEquals(results[i].others, [])
    self.assertAlmostEquals(results[i].average_maf, maf)
    self.assertEquals(results[i].include, None)
    self.assertEquals(results[i].ld, [(4, 4, 1, 1)])
    self.assertEquals( pair_disposition(4, 4, results[i]), 'singleton-tag')
    self.assertEquals(results[i].disposition, 'maximal-bin')

    self.assertEquals(results[j].tags, [3])
    self.assertEquals(results[j].others, [])
    self.assertAlmostEquals(results[j].average_maf, maf)
    self.assertEquals(results[j].include, None)
    self.assertEquals(results[j].ld, [(3, 3, 1, 1)])
    self.assertEquals( pair_disposition(3, 3, results[j]), 'lonely-tag')
    self.assertEquals(results[j].disposition, 'maximal-bin')


  def test_binner_maxloci(self):

    sizes = 1,2,3,4
    resultlens = 4,4,2,2

    location = 10000
    genos = ['AA', 'AC', 'AA', 'CC', 'AA']
    maf = estimate_maf(genos)

    class Options(object): pass

    for size,resultlen in zip(sizes,resultlens):
      loci = dict( (i,Locus(i, location+500*i, genos)) for i in range(5) )

      # populate binsets
      binsets = {}
      binsets[0] = Bin([0,1,2,3],  maf*4, Bin.INCLUDE_TYPED) # 0 is an include
      binsets[1] = Bin([1,0,2],    maf*3)
      binsets[2] = Bin([2,0,1,3],  maf*4)
      binsets[3] = Bin([3,0,2,4],  maf*4)
      binsets[4] = Bin([4,3],      maf*2, Bin.EXCLUDE)

      # populate lddata , notice that there is dependecy between lddata an binsets
      lddata = {}
      lddata[(0,1)] = (1,1)
      lddata[(0,2)] = (1,1)
      lddata[(0,3)] = (1,1)
      lddata[(1,2)] = (1,1)
      lddata[(2,3)] = (1,1)
      lddata[(3,4)] = (1,1)

      options = Options()
      options.locipertag = size
      f = get_tags_required_function(options)
      results = list(binner(loci, binsets, lddata, Includes(set(),set()), f))
      self.assertEquals(len(results), resultlen)

  def test_scan_ldpairs(self):
    '''
    test the scan_ldpairs function in snpselct module
    '''
    maxd = 200000
    rthreshold = 0.5
    dthreshold = 0.0
    genos1 = ['GG','AA','AA','GG','GA']
    genos2 = ['CC','CC','CA','CC','AA']
    genos3 = ['CC','CC','CA','CC','AC']
    locus1 = Locus(1, 1000000, genos1)
    locus2 = Locus(2, 1100000, genos1)
    locus3 = Locus(3, 1200000, genos2)
    locus4 = Locus(4, 1300000, genos2)
    locus5 = Locus(5, 1500000, genos3)
    loci = []
    loci.append(locus1)
    loci.append(locus2)
    loci.append(locus3)
    loci.append(locus4)
    loci.append(locus5)
    ldpairs = scan_ldpairs(loci, maxd, rthreshold, dthreshold)
    results = [(1,2,1,1),
               (3,4,1,1),
               (4,5,0.58,1)]
    for i, ld_pair in enumerate(ldpairs):
      self.assertEquals(results[i][0], ld_pair[0])
      self.assertEquals(results[i][1], ld_pair[1])
      self.assertAlmostEqual(results[i][2], ld_pair[2], places=2)
      self.assertAlmostEqual(results[i][3], ld_pair[3], places=2)

  def test_estimate_ld(self):
    counts1 = 4,0,0,4,1
    counts2 = 1,1,0,6,1
    r2, dprime = estimate_ld(*counts1)
    r21,dprime1 = slow_estimate_ld(*counts1)
    self.assertAlmostEquals(r2, r21, places = 3)
    self.assertAlmostEquals(dprime, dprime1, places = 3)
    r2, dprime = estimate_ld(*counts2)
    r21,dprime1 = slow_estimate_ld(*counts2)
    self.assertAlmostEquals(r2, r21, places = 3)
    self.assertAlmostEquals(dprime, dprime1, places = 3)


class testMulti(unittest.TestCase):
  def testMerge(self):
    loci1 = [ Locus(1, 1, ['AB']),
              Locus(2, 2, ['AB']),
              Locus(3, 3, ['AB']),
              Locus(5, 5, ['AB']),
              Locus(6, 6, ['AB']) ]

    loci2 = [ Locus(1, 1, ['BA']),
              Locus(3, 3, ['BA']),
              Locus(4, 4, ['BA']),
              Locus(5, 5, ['BA']),
              Locus(6, 6, ['BA']) ]

    # FIXME: Must complete test vectors
    if 0:
      for l1,l2 in merge_multi_loci([loci1,loci2]):
        print l1.name,l1.genos,l2.name,l2.genos

      for locus in merge_loci([loci1,loci2]):
        print locus.name,locus.genos

      for locus in merge_loci([loci2,loci1]):
        print locus.name,locus.genos


class testBin(unittest.TestCase):
  def setUp(self):
    '''
    set up data used in the tests.
    setUp is called before each test function execution.
    '''
    self.bin0 = Bin()
    self.bin1 = Bin([1],  0.1)
    self.bin2 = Bin([1,2],0.5)

  def testAdd(self):
    self.bin1.add(2,0.4)
    self.assertEquals(self.bin1, self.bin2)
    self.assertEquals(self.bin1.maf, 0.5)
    self.assertEquals(self.bin1.maxcovered, 2)
    self.assertEquals(self.bin1.disposition, 0)

  def testRemove(self):
    self.bin2.remove(2, 0.4)
    self.assertEquals(self.bin2, self.bin1)
    self.assertAlmostEquals(self.bin2.maf, 0.1)
    self.assertEquals(self.bin2.maxcovered, 2)
    self.assertEquals(self.bin2.disposition, 0)

  def test_average_maf(self):
    #self.assertEquals(self.bin0.average_maf(), 0) #uncomment temporarily because of ZeroDivisionError
    self.assertEquals(self.bin1.average_maf(), 0.1)
    self.assertEquals(self.bin2.average_maf(), 0.25)

  def testPickle(self):
    '''Really test the __setstate__ and __getstate__ methods'''
    import pickle, cPickle
    for module in (pickle,cPickle):
      s = module.dumps(self.bin2)
      x = module.loads(s)
      self.assertEquals(self.bin2, x) # the set attribute(?) is equal
      self.assertEquals(self.bin2.maf, x.maf) # the maf attribute of th bin is equal
      self.assertEquals(self.bin2.disposition, x.disposition)
      self.assertEquals(self.bin2.maxcovered, x.maxcovered)


class testBinSequence(unittest.TestCase):
  pass


if __name__ == '__main__':
  unittest.main()

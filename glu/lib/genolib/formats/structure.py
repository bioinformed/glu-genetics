# -*- coding: utf-8 -*-

from __future__ import with_statement

__abstract__  = 'GLU Structure genotype format output'
__copyright__ = 'Copyright (c) 2007-2009, BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.'
__license__   = 'See GLU license for terms by running: glu license'
__revision__  = '$Id$'

__all__       = ['StructureWriter', 'save_structure']

__genoformats__ = [
#  LOADER    SAVER              WRITER        PFORMAT    ALIAS     EXTS
  (None, 'save_structure', 'StructureWriter', 'sdat', 'structure', None) ]


from   itertools                 import chain

from   glu.lib.fileutils         import autofile,related_file,map_reader,   \
                                        parse_augmented_filename,get_arg

from   glu.lib.genolib.streams   import GenomatrixStream


mainparam_template = '''\
KEY PARAMETERS FOR THE PROGRAM structure.  YOU WILL NEED TO SET THESE
IN ORDER TO RUN THE PROGRAM.  VARIOUS OPTIONS CAN BE ADJUSTED IN THE
FILE extraparams.

"(int)" means that this takes an integer value.
"(B)"   means that this variable is Boolean
        (ie insert 1 for True, and 0 for False)
"(str)" means that this is a string (but not enclosed in quotes!)

Data File

#define INFILE %(INFILE)s   // name of input data file
#define OUTFILE %(OUTFILE)s       // name of output data file

#define NUMINDS %(NUMINDS)s        // number of diploid individuals in data file
#define NUMLOCI %(NUMLOCI)s        // number of loci in data file
#define LABEL %(LABEL)s           // Input file contains individual labels
#define POPDATA %(POPDATA)s           // Input file contains a population identifier
#define POPFLAG %(POPFLAG)s           // Input file contains a flag which says
                             //     whether to use popinfo when USEPOPINFO==1
#define PHENOTYPE %(PHENOTYPE)s           // Input file contains phenotype information
#define EXTRACOLS %(EXTRACOLS)s           // Number of additional columns of data
                             //       before the genotype data start.
#define PHASEINFO %(PHASEINFO)s           // the data for each individual contains a line
                                    indicating phase
#define MARKOVPHASE %(MARKOVPHASE)s         // the phase info follows a Markov model.

#define MISSING %(MISSING)s       // value given to missing genotype data
#define PLOIDY %(PLOIDY)s        // ploidy of data

#define ONEROWPERIND %(ONEROWPERIND)s        // store data for individuals in a single line
#define MARKERNAMES %(MARKERNAMES)s        // data file contains row of marker names
#define MAPDISTANCES %(MAPDISTANCES)s        // data file contains row of map distances
                             //     between loci

Program Parameters

#define MAXPOPS %(MAXPOPS)s           // number of populations assumed
#define BURNIN %(BURNIN)s       // length of burnin period
#define NUMREPS %(NUMREPS)s       // number of MCMC reps after burnin
'''

mainparam_required = ['MAXPOPS','OUTFILE']
mainparam_defaults = {
  #'INFILE'
  #'OUTFILE'
  #'NUMINDS'
  #'NUMLOCI'
  #'POPDATA'
  #'POPFLAG'
  #'MAXPOPS'
  'LABEL'       :  1,
  'PHENOTYPE'   :  0,
  'EXTRACOLS'   :  0,
  'PHASEINFO'   :  0,
  'MARKOVPHASE' :  0,
  'MISSING'     : -9,
  'PLOIDY'      :  2,
  'ONEROWPERIND':  0,
  'MARKERNAMES' :  1,
  'MAPDISTANCES':  0,
  'BURNIN'      :  5000,
  'NUMREPS'     :  10000,
}

extraparam_template = '''\
EXTRA PARAMS FOR THE PROGRAM structure.  THE MOST IMPORTANT PARAMS ARE
IN THE FILE mainparams.

"(int)" means that this takes an integer value.
"(d)"   means that this is a double (ie, a Real number such as 3.14).
"(B)"   means that this variable is Boolean
        (ie insert 1 for True, and 0 for False).

PROGRAM OPTIONS

#define FREQSCORR %(FREQSCORR)s        // allele frequencies are correlated among pops
#define ONEFST %(ONEFST)s        // assume same value of Fst for all subpopulations.

#define INFERALPHA %(INFERALPHA)s        // Infer ALPHA (the admixture parameter)
#define POPALPHAS %(POPALPHAS)s        // Individual alpha for each population

#define INFERLAMBDA %(INFERLAMBDA)s        // Infer LAMBDA (the allele frequencies parameter)
#define POPSPECIFICLAMBDA %(POPSPECIFICLAMBDA)s  // infer a separate lambda for each pop
                            //     (only if INFERLAMBDA=1).

#define NOADMIX %(NOADMIX)s        // Use no admixture model
#define LINKAGE %(LINKAGE)s        // Use the linkage model model
#define PHASED %(PHASED)s        // Data are in correct phase (required unless data are diploid)
                            //     If (LINKAGE=1, PHASED=0), then PHASEINFO can be used--this is an
                            //     extra line in the input file that gives phase probabilities.
                            //     When PHASEINFO =0 each value is set to 0.5, implying no phase information.

#define LOG10RMIN %(LOG10RMIN)s   // Log10 of minimum allowed value of r under linkage model
#define LOG10RMAX %(LOG10RMAX)s   // Log10 of maximum allowed value of r
#define LOG10RPROPSD %(LOG10RPROPSD)s   // standard deviation of log r in update
#define LOG10RSTART %(LOG10RSTART)s   // initial value of log10 r

#define COMPUTEPROB %(COMPUTEPROB)s        // Estimate the probability of the Data under
                            //     the model.  This is used when choosing the
                            //     best number of subpopulations.

#define ADMBURNIN %(ADMBURNIN)s      // initial period of burnin with admixture model.
                            // The admixture model normally has the best mixing properties,
                            // and therefore should be used as the first phase of the burnin.
                            // A significant burnin is highly recommended under the linkage
                            // model which otherwise can perform very badly.

USING PRIOR POPULATION INFO

#define USEPOPINFO %(USEPOPINFO)s        // Use prior population information to assist
                            //     clustering.  Need POPDATA=1.
#define GENSBACK %(GENSBACK)s        // For use when inferring whether an indiv-
                            //      idual is an immigrant, or has an immigrant an-
                            //      cestor in the past GENSBACK generations.  eg, if
                            //      GENSBACK==2, it tests for immigrant ancestry
                            //      back to grandparents.
#define MIGRPRIOR %(MIGRPRIOR)s      // prior prob that an individual is a migrant
                            //     (used only when USEPOPINFO==1).  This should
                            //     be small, eg 0.01 or 0.1.
#define PFROMPOPFLAGONLY %(PFROMPOPFLAGONLY)s   // only use individuals with POPFLAG=1 to update P.
                            //     This is to enable use of a reference set of
                            //     individuals for clustering additional "test"
                            //     individuals.

OUTPUT OPTIONS

#define PRINTKLD %(PRINTKLD)s       // Print estimated Kullback-Leibler diver-
                            //     gence to screen during the run
#define PRINTLAMBDA %(PRINTLAMBDA)s       // Print current value(s) of lambda to screen
#define PRINTQSUM %(PRINTQSUM)s       // Print summary of current population membership to screen



#define SITEBYSITE %(SITEBYSITE)s       // whether or not to print site by site results. Large!

#define PRINTQHAT %(PRINTQHAT)s       // Q-hat printed to a separate file.  Turn this
                            //     on before using STRAT.
#define UPDATEFREQ %(UPDATEFREQ)s       // frequency of printing update on the screen.
                            //       Set automatically if this is 0.
#define PRINTLIKES %(PRINTLIKES)s       // print current likelihood to screen every rep
#define INTERMEDSAVE %(INTERMEDSAVE)s       // number of saves to file during run

#define ECHODATA %(ECHODATA)s       // Print some of data file to screen to check
                            //     that the data entry is correct.

(NEXT 3 ARE FOR COLLECTING DISTRIBUTION OF Q:)

#define ANCESTDIST %(ANCESTDIST)s       // collect data about the distribution of an-
                                   cestry coefficients (Q) for each individual
#define NUMBOXES %(NUMBOXES)s      // the distribution of Q values is stored as
                                   a histogram with this number of boxes.
#define ANCESTPINT %(ANCESTPINT)s      // the size of the displayed probability
                                   interval on Q (values between 0.0--1.0)
PRIORS

#define ALPHA %(ALPHA)s       // Dirichlet parameter for degree of admixture
                            //     (this is the initial value if INFERALPHA==1).
#define FPRIORMEAN %(FPRIORMEAN)s      // Prior mean and SD of Fst for pops.
#define FPRIORSD %(FPRIORSD)s      // The prior is a Gamma distribution with these parameters

#define LAMBDA %(LAMBDA)s      // Dirichlet parameter for allele frequencies
#define UNIFPRIORALPHA %(UNIFPRIORALPHA)s     // use a uniform prior for alpha;
                            //     otherwise gamma prior
#define ALPHAMAX %(ALPHAMAX)s    // max value of alpha if uniform prior
#define ALPHAPRIORA %(ALPHAPRIORA)s    // (only if UNIFPRIORALPHA==0): alpha has a gamma
                            //     prior with mean A*B, and
#define ALPHAPRIORB %(ALPHAPRIORB)s   // variance A*B^2.  Suggest A=0.1, B=1.0

MISCELLANEOUS

#define ALPHAPROPSD %(ALPHAPROPSD)s    // SD of proposal for updating alpha
#define STARTATPOPINFO %(STARTATPOPINFO)s     // Use given populations as the initial condition
                            //     for population origins.  (Need POPDATA==1).  It
                            //     is assumed that the PopData in the input file
                            //     are between 1 and k where k<=MAXPOPS.
#define RANDOMIZE %(RANDOMIZE)s       // use new random seed for each run
#define METROFREQ %(METROFREQ)s      // Frequency of using Metropolis step to update
                            //      Q under admixture model (ie use the metr. move every
                            //      i steps).  If this is set to 0, it is never used.
                            //      (Proposal for each q^(i) sampled from prior.  The
                            //      goal is to improve mixing for small alpha.)
#define REPORTHITRATE %(REPORTHITRATE)s      // report hit rate if using METROFREQ
'''

extraparam_defaults = {
  'FREQSCORR' : 1,
  'ONEFST' : 0,
  'INFERALPHA' : 1,
  'POPALPHAS' : 0,
  'INFERLAMBDA' : 0,
  'POPSPECIFICLAMBDA' : 0,
  'NOADMIX' : 0,
  'LINKAGE' : 0,
  'PHASED' : 0,
  'LOG10RMIN' : -4.0,
  'LOG10RMAX' : 1.0,
  'LOG10RPROPSD' : 0.1,
  'LOG10RSTART' : -2.0,
  'COMPUTEPROB' : 1,
  'ADMBURNIN' : 500,
  'USEPOPINFO' : 1,
  'GENSBACK' : 1,
  'MIGRPRIOR' : 0.001,
  'PFROMPOPFLAGONLY' : 1,
  'PRINTKLD' : 1,
  'PRINTLAMBDA' : 1,
  'PRINTQSUM' : 1,
  'SITEBYSITE' : 0,
  'PRINTQHAT' : 0,
  'UPDATEFREQ' : 5,
  'PRINTLIKES' : 0,
  'INTERMEDSAVE' : 0,
  'ECHODATA' : 1,
  'ANCESTDIST' : 0,
  'NUMBOXES' : 1000,
  'ANCESTPINT' : 0.90,
  'ALPHA' : 1.0,
  'FPRIORMEAN' : 0.01,
  'FPRIORSD' : 0.05,
  'LAMBDA' : 1.0,
  'UNIFPRIORALPHA' : 1,
  'ALPHAMAX' : 20.0,
  'ALPHAPRIORA' : 0.05,
  'ALPHAPRIORB' : 0.001,
  'ALPHAPROPSD' : 0.025,
  'STARTATPOPINFO' : 0,
  'RANDOMIZE' : 1,
  'METROFREQ' : 10,
  'REPORTHITRATE' : 0,
}


class StructureWriter(object):
  '''
  Object to write Structure data

  >>> loci =           ('l1',      'l2',      'l3')
  >>> rows = [('s1', [('A','A'),(None,None),('C','T')]),
  ...         ('s2', [('A','G'), ('C','G'), ('C','C')]),
  ...         ('s3', [('G','G'),(None,None),('C','T')]) ]
  >>> genos = GenomatrixStream.from_tuples(rows,'sdat',loci=loci)
  >>> from cStringIO import StringIO
  >>> o = StringIO()
  >>> with StructureWriter(o,'structure',genos.loci,genos.genome,genos.phenome) as w:
  ...   genos=iter(genos)
  ...   w.writerow(*genos.next())
  ...   w.writerow(*genos.next())
  ...   w.writerows(genos)
  >>> print o.getvalue() # doctest: +NORMALIZE_WHITESPACE
  l1 l2 l3
  s1 1 -9 1
  s1 1 -9 2
  s2 1 1 1
  s2 2 2 1
  s3 2 -9 1
  s3 2 -9 2
  '''
  def __init__(self,filename,format,loci,genome,phenome,extra_args=None,**kwargs):
    '''
    @param     filename: file name or file object
    @type      filename: str or file object
    @param       format: data format string
    @type        format: str
    @param       header: column headings
    @type        header: list or str
    @param     genorepr: object representing the input/output encoding and
                         internal representation of genotypes
    @type      genorepr: UnphasedMarkerRepresentation or similar object
    '''
    if extra_args is None:
      args = kwargs
    else:
      args = extra_args
      args.update(kwargs)

    filename = parse_augmented_filename(filename,args)

    mainfile  = get_arg(args, ['mainparam'])
    extrafile = get_arg(args, ['extraparam'])
    popfile   = get_arg(args, ['populations','pops'])
    defpop    = get_arg(args, ['default_population', 'def_pop', 'defpop'])

    main = mainparam_defaults.copy()
    for param in chain(mainparam_required,mainparam_defaults):
      value = get_arg(args,[param])
      if value is not None:
        main[param] = value

    extra = extraparam_defaults.copy()
    for param in extraparam_defaults:
      value = get_arg(args,[param])
      if value is not None:
        extra[param] = value

    if extra_args is None and args:
      raise ValueError('Unexpected filename arguments: %s' % ','.join(sorted(args)))

    # Careful: file=<blank> is intended to suppress output
    if mainfile is None:
      mainfile = related_file(filename,'mainparam')
    if extrafile is None:
      extrafile = related_file(filename,'extraparam')

    main['INFILE'] = filename

    if main.get('OUTFILE') is None:
      main['OUTFILE'] = related_file(filename,'out') or 'structure.out'

    populations = None
    if popfile is not None:
      populations = map_reader(popfile)

    if defpop is None:
      defpop = '1'

    self.out         = autofile(filename,'w')
    self.loci        = loci
    self.genome      = genome
    self.main        = main
    self.extra       = extra
    self.mainfile    = mainfile
    self.extrafile   = extrafile
    self.populations = populations
    self.defpop      = defpop
    self.samples     = 0

    self.out.write(' '.join(self.loci))
    self.out.write('\n')

  def writerow(self, sample, genos):
    '''
    Write a row of genotypes given the row key and list of genotypes

    @param rowkey: row identifier
    @type  rowkey: str
    @param  genos: sequence of genotypes in an internal representation
    @type   genos: sequence
    '''
    out = self.out
    if out is None:
      raise IOError('Cannot write to closed writer object')

    if len(genos) != len(self.loci):
      raise ValueError('[ERROR] Internal error: Genotypes do not match header')

    missing = self.main['MISSING']

    prefix = [sample]
    if self.populations is not None:
      pop = self.populations.get(sample)
      if pop is None:
        prefix += [self.defpop,'0']
      else:
        prefix += [pop,'1']

    row1 = prefix+[ str(g.allele1_index or missing) for g in genos ]
    row2 = prefix+[ str(g.allele2_index or missing) for g in genos ]

    out.write(' '.join(row1))
    out.write('\n')
    out.write(' '.join(row2))
    out.write('\n')

    self.samples += 1

  def writerows(self, rows):
    '''
    Write rows of genotypes given pairs of row key and list of genotypes

    @param rows: sequence of pairs of row key and sequence of genotypes in
                 an internal representation
    @type  rows: sequence of (str,sequence)
    '''
    out = self.out
    if out is None:
      raise IOError('Cannot write to closed writer object')

    n = len(self.loci)
    defpop  = self.defpop
    missing = self.main['MISSING'] or '-9'

    for sample,genos in rows:
      if len(genos) != n:
        raise ValueError('[ERROR] Internal error: Genotypes do not match header')

      prefix = [sample]
      if self.populations is not None:
        pop = self.populations.get(sample)
        if pop is None:
          prefix += [defpop,'0']
        else:
          prefix += [pop,'1']

      row1 = prefix+[ str(g.allele1_index or missing) for g in genos ]
      row2 = prefix+[ str(g.allele2_index or missing) for g in genos ]

      out.write(' '.join(row1))
      out.write('\n')
      out.write(' '.join(row2))
      out.write('\n')

      self.samples += 1

  def close(self):
    '''
    Close the writer

    A closed writer cannot be used for further I/O operations and will
    result in an error if called more than once.
    '''
    if self.out is None:
      raise IOError('Writer object already closed')

    # FIXME: Closing out causes problems with StringIO objects used for
    #        testing
    #self.out.close()
    self.out = None

    self.main['NUMINDS'] = self.samples
    self.main['NUMLOCI'] = len(self.loci)

    if self.populations is not None:
      self.main.setdefault('POPDATA', 1)
      self.main.setdefault('POPFLAG', 1)
      self.main.setdefault('MAXPOPS', len(set(self.populations.itervalues())))
    else:
      self.main['POPDATA'] = 0
      self.main['POPFLAG'] = 0
      self.main.setdefault('MAXPOPS',1)

    if self.mainfile is not None:
      autofile(self.mainfile,'w').write(mainparam_template % self.main)

    if self.extrafile is not None:
      autofile(self.extrafile,'w').write(extraparam_template % self.extra)

  def __enter__(self):
    '''
    Context enter function
    '''
    return self

  def __exit__(self, *exc_info):
    '''
    Context exit function that closes the writer upon exit
    '''
    self.close()


def save_structure(filename,genos,format,extra_args=None,**kwargs):
  '''
  Write the genotype matrix data to file.

  @param     filename: file name or file object
  @type      filename: str or file object
  @param        genos: genomatrix stream
  @type         genos: sequence

  >>> from cStringIO import StringIO
  >>> o = StringIO()
  >>> loci =              ('l1',     'l2',    'l3')
  >>> rows = [('s1', [('A','A'),(None,None),('C','T')]),
  ...           ('s2', [('A','G'), ('C','G'), ('C','C')]),
  ...           ('s3', [('G','G'),(None,None),('C','T')]) ]
  >>> genos = GenomatrixStream.from_tuples(rows,'sdat',loci=loci)
  >>> save_structure(o,genos,'structure')
  >>> print o.getvalue() # doctest: +NORMALIZE_WHITESPACE
  l1 l2 l3
  s1 1 -9 1
  s1 1 -9 2
  s2 1 1 1
  s2 2 2 1
  s3 2 -9 1
  s3 2 -9 2
  '''
  if extra_args is None:
    args = kwargs
  else:
    args = extra_args
    args.update(kwargs)

  filename  = parse_augmented_filename(filename,args)

  mergefunc = get_arg(args, ['mergefunc'])

  genos = genos.as_sdat(mergefunc)

  with StructureWriter(filename, format, genos.loci, genos.genome, genos.phenome, extra_args=args) as writer:

    if extra_args is None and args:
      raise ValueError('Unexpected filename arguments: %s' % ','.join(sorted(args)))

    writer.writerows(genos)


def test():
  import doctest
  return doctest.testmod()


if __name__ == '__main__':
  test()

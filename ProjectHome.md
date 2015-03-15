# GLU: GENOTYPE LIBRARY AND UTILITIES #

Whole-genome association studies are generating unprecedented amounts of genotype data, frequently billions of genotypes per study, and require new and scalable computational approaches to address the storage, management, quality control, and genetic analysis. GLU is a framework and a software package that was designed around a set of novel conceptual approaches. GLU addresses the need for general and powerful tools that can scale to effectively handle trillions of genotypes.

# Key innovations #
  * Compressed binary genotype storage
  * Data sets are usually not loaded into main memory unless absolutely necessary.  Instead many of the analytical functions read the data in small chunks, perform the requested analysis, then move on to the next.
  * Integration with a high-level scripting language for easy customization and extension
  * Support for parallel processing and distributed computing (more on this in later versions)

# Data management features #
  * The ability to import, export, merge, and split genotype data among several common formats and standards
  * Filter based on powerful criteria for inclusion, exclusion
  * Rename and adjust sample and locus metadata

# Genotype quality assurance #
  * Estimation of assay completion and other summary statistics
  * Reproducibility and concordance
  * Verification of known and detection of unknown duplicate samples
  * Empirical sex determination
  * Testing for deviations from Hardy-Weinberg proportions, Mendelian inheritance patterns, non-random patterns of missing data.

# Analytic tools #
  * Fitting generalized linear models to test for phenotype/genotype association
  * Fast linkage disequilibrium (LD) estimation
  * Fast tag SNP estimation with the ability to augment SNPs from set panels with a optimal tag SNPs using flexible criteria including design scores from major genotyping vendors.

# Extensibility #
Viewed as a library or framework, GLU is designed to be highly extensible, so that it may be easily augmented, customized, and serve as a foundation for the rapid development of new applications.

# License Information #
GLU is released under a [modified BSD license](http://cgf.nci.nih.gov/glu/docs/LICENSE.html)

# Documentation (still a work in progress) #

[GLU 1.0b1 Documentation](http://cgf.nci.nih.gov/glu/docs/1.0b1)

[GLU 1.0b2 Development Documentation](http://cgf.nci.nih.gov/glu/docs/1.0b2)

# Genomic Annotation Files (GENEdb) #

  * [GENEdb for NCBI genome build 36.3 (hg18) and dbSNP 130](http://cgf.nci.nih.gov/glu/files/genedb_hg18_snp130_rtree.db)

  * [GENEdb for GRCh37 genome build (hg19) and dbSNP 131](http://cgf.nci.nih.gov/glu/files/genedb_hg19_snp131_rtree.db)

# Community #
  * [GLU users group](http://groups.google.com/group/glu-users) for discussion and help with GLU
  * [GLU developer group](http://groups.google.com/group/glu-dev) for those involved in improving and extending GLU

# Acknowledgments #
GLU is primarily developed and maintained by Kevin Jacobs <jacobs at
bioinformed dot com> to support the Cancer Genetic Markers of Susceptibility
([CGEMS](http://cgems.cancer.gov/)) project, an initiative by the National Cancer Institute ([NCI](http://cancer.gov/)) Division of Cancer Epidemiology and Genetics ([DCEG](http://dceg.cancer.gov/)) and run by the NCI Core Genotyping Facility ([CGF](http://cgf.nci.nih.gov/)) to identify genetic alterations that affect susceptible to prostate and breast cancer.

CGEMS is funded by NCI under Contract N01-CO-12400 by [SAIC-Frederick](http://saic.ncifcrf.gov/), a subsidiary of Science Applications International Corporation ([SAIC](http://www.saic.com/)).

[Contributors](http://glu-genetics.googlecode.com/svn/trunk/docs/contributors.rst)

## Widgets ##


&lt;wiki:gadget url="http://www.ohloh.net/p/19871/widgets/project\_languages.xml" border="1" height="220" width="340"/&gt;&lt;wiki:gadget url="http://www.ohloh.net/p/19871/widgets/project\_factoids.xml" height="175" width="340" border="1"/&gt;&lt;wiki:gadget url="http://www.ohloh.net/p/19871/widgets/project\_basic\_stats.xml" height="260" width="340" border="1"/&gt;
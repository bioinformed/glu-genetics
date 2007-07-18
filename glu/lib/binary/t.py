import tables

gfile = tables.openFile('test.h5',mode='w')
gfile.root._v_attrs.format = 'genotriple'

class TripleDesc(tables.IsDescription):
  sample = tables.UInt32Col(pos=1)
  locus  = tables.UInt32Col(pos=2)
  geno   = tables.UInt8Col(pos=3)

filters = tables.Filters(complevel=1, complib='lzo',shuffle=True)
genos = gfile.createTable(gfile.root, 'genotypes', TripleDesc,
                          'Sequence of encoded sample, locus, genotype triples',
                          filters=filters, expectedrows=5000000)

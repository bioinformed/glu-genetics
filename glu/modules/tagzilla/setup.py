from distutils.core import setup, Extension

module1 = Extension('tagzillac', sources = ['tagzillac.c'])
module2 = Extension('pqueue',    sources = ['pqueue.c'])

if __name__=='__main__':
  setup(name = 'tagzilla',
        version = '1.2alpha',
        author = 'Kevin Jacobs',
        author_email = 'jacobs@bioinformed.com',
        description = 'A robust and fast SNP binning and tagging program',
        py_modules = ['tagzilla','binsum','coverage','ldmatrix','surrogates'],
        scripts = ['tagzilla','binsum','coverage','ldmatrix','surrogates'],
        ext_modules = [module1, module2])

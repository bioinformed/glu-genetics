from distutils.core import setup, Extension

module1 = Extension('dupcheckc', sources = ['dupcheckc.c'])

setup (name = 'dupcheckc',
       version = '1.0',
       description = 'This is a dupcheckc package',
       ext_modules = [module1])
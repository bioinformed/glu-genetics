# -*- coding: utf-8 -*-
'''
File:          dbexceptions.py

Author:        Kevin Jacobs (jacobs@bioinformed.com)

Created:       January 21, 2001

Purpose:       Abstract database exceptions

Compatibility: Python 2.2 and above

Requires:

Revision:      $Id: dbexceptions.py 202 2006-06-13 18:17:19Z jacobske $

Copyright (c) 2005 BioInformed Consulting Services.
'''

dblib_exceptions = [ 'TooFewRowsError', 'TooManyRowsError' ]

dbapi_exceptions = [ 'Warning',
                     'Error',
                     'InterfaceError',
                     'DatabaseError',
                     'DataError',
                     'OperationalError',
                     'IntegrityError',
                     'InternalError',
                     'ProgrammingError',
                     'NotSupportedError' ]

#######################################################################################

class Warning                  (StandardError)       : pass
class Error                    (StandardError)       : pass
class   InterfaceError           (Error)             : pass
class   DatabaseError            (Error)             : pass
class     DataError                (DatabaseError)   : pass
class       TooFewRowsError          (DataError)     : pass
class       TooManyRowsError         (DataError)     : pass
class     OperationalError         (DatabaseError)   : pass
class     IntegrityError           (DatabaseError)   : pass
class     InternalError            (DatabaseError)   : pass
class     ProgrammingError         (DatabaseError)   : pass
class     NotSupportedError        (DatabaseError)   : pass


#######################################################################################

def __import_exceptions(module):
  for e in dbapi_exceptions:
    sub_exception    = getattr(module, e)
    global_exception = globals()[e]
    sub_exception.__bases__ += (global_exception,)

#######################################################################################

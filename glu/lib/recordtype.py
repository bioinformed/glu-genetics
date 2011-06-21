# -*- coding: utf-8 -*-

__all__ = ['recordtype']


import sys as _sys

from   keyword   import iskeyword as _iskeyword
from   itertools import izip      as _izip


def recordtype(typename, field_names, verbose=False, **default_kwds):
    '''Returns a new class with named fields.

    Based on http://code.activestate.com/recipes/576555-records/ by George Sakkis.

    @keyword defaults: A mapping from (a subset of) field names to default
        values.
    @keyword default: If provided, the default value for all fields without an
        explicit default in `defaults`.

    >>> Point = recordtype('Point', 'x y', default=0)
    >>> Point.__doc__           # docstring for the new class
    'Point(x, y)'
    >>> Point()                 # instantiate with defaults
    Point(x=0, y=0)
    >>> p = Point(11, y=22)     # instantiate with positional args or keywords
    >>> p[0] + p.y              # accessible by name and index
    33
    >>> p.x = 100
    >>> p[1] =200               # modifiable by name and index
    >>> p
    Point(x=100, y=200)
    >>> 100 in p                # Supports containment checks
    True
    >>> 200 in p
    True
    >>> 300 in p
    False
    >>> x, y = p                # unpack
    >>> x, y
    (100, 200)
    >>> p[:1],p[1:]             # slicing
    ([100], [200])
    >>> p[:] = 200,100
    >>> p
    Point(x=200, y=100)
    >>> d = p._asdict()         # convert to a dictionary
    >>> d['x']
    200
    >>> Point(**d) == p         # convert from a dictionary
    True
    '''
    # Parse and validate the field names.  Validation serves two purposes,
    # generating informative error messages and preventing template injection attacks.
    if isinstance(field_names, basestring):
        # names separated by whitespace and/or commas
        field_names = field_names.replace(',', ' ').split()
    field_names = tuple(map(str, field_names))

    for name in (typename,) + field_names:
        if not min(c.isalnum() or c=='_' for c in name):
            raise ValueError('Type names and field names can only contain alphanumeric characters and underscores: %r' % name)
        if _iskeyword(name):
            raise ValueError('Type names and field names cannot be a keyword: %r' % name)
        if name[0].isdigit():
            raise ValueError('Type names and field names cannot start with a number: %r' % name)

    seen_names = set()
    for name in field_names:
        if name.startswith('_'):
            raise ValueError('Field names cannot start with an underscore: %r' % name)
        if name in seen_names:
            raise ValueError('Encountered duplicate field name: %r' % name)
        seen_names.add(name)

    # determine the func_defaults of __init__
    defaults = default_kwds.pop('defaults', {})

    if 'default' in default_kwds:
        default = default_kwds.pop('default')
        init_defaults = tuple(defaults.get(f,default) for f in field_names)
    elif not defaults:
        init_defaults = None
    else:
        default_fields = field_names[-len(defaults):]
        if set(default_fields) != set(defaults):
            raise ValueError('Missing default parameter values')
        init_defaults = tuple(defaults[f] for f in default_fields)

    if default_kwds:
        raise ValueError('Invalid keyword arguments: %s' % default_kwds)

    # Create and fill-in the class template
    numfields = len(field_names)
    argtxt   = ', '.join(field_names)
    reprtxt  = ', '.join('%s=%%r' % f for f in field_names)
    dicttxt  = ', '.join('%r: self.%s' % (f,f) for f in field_names)
    tupletxt = repr(tuple('self.%s' % f for f in field_names)).replace("'",'')
    inittxt  = '; '.join('self.%s=%s' % (f,f) for f in field_names) or 'pass'
    itertxt  = '; '.join('yield self.%s' % f for f in field_names) or 'return iter([])'
    eqtxt    = ' and '.join('self.%s==other.%s' % (f,f) for f in field_names) or 'True'
    conttxt  = ' or '.join('item==self.%s' % f for f in field_names) or 'False'
    setstxt  = '%s = state' % tupletxt if numfields else 'pass'

    template = '''
class %(typename)s(object):
    '%(typename)s(%(argtxt)s)'

    __slots__  = %(field_names)r
    __hash__   = None

    @classmethod
    def _make(cls, iterable):
        'Make a new %(typename)s object from a sequence or iterable'
        return cls(*iterable)

    def __init__(self, %(argtxt)s):
        %(inittxt)s

    def __len__(self):
        return %(numfields)d

    def __iter__(self):
        %(itertxt)s

    def __getitem__(self, index):
        if isinstance(index, slice):
            return [ getattr(self, s) for s in self.__slots__[index] ]
        else:
            return getattr(self, self.__slots__[index])

    def __setitem__(self, index, value):
        if isinstance(index, slice):
            slots = self.__slots__[index]
            try:
                value_len = len(value)
            except TypeError:
                value = list(value)
                value_len = len(value)
            if value_len != len(slots):
                raise IndexError('Invalid length slice assignment')
            for s,v in _izip(slots,value):
                setattr(self, s, v)
        else:
            setattr(self, self.__slots__[index], value)

    def _asdict(self):
        'Return a new dict which maps field names to their values'
        return {%(dicttxt)s}

    def __repr__(self):
        return '%(typename)s(%(reprtxt)s)' %% %(tupletxt)s

    def __eq__(self, other):
        return isinstance(other, self.__class__) and %(eqtxt)s

    def __ne__(self, other):
        return not self==other

    def __contains__(self, item):
        return %(conttxt)s

    def __getstate__(self):
        return %(tupletxt)s

    def __setstate__(self, state):
        %(setstxt)s
''' % locals()

    if verbose:
        print template

    # Execute the template string in a temporary namespace
    namespace = dict( __name__='recordtype_%s' % typename,_izip=_izip)
    try:
        exec template in namespace
    except SyntaxError, e:
        raise SyntaxError(e.message + ':\n' + template)

    cls = namespace[typename]
    cls.__init__.im_func.func_defaults = init_defaults

    # For pickling to work, the __module__ variable needs to be set to the frame
    # where the named tuple is created.  Bypass this step in environments where
    # _sys._getframe is not defined (Jython for example).
    try:
        cls.__module__ = _sys._getframe(1).f_globals.get('__name__','__main__')
    except (AttributeError, ValueError):
        pass

    return cls


if __name__ == '__main__':
    import doctest
    TestResults = recordtype('TestResults', 'failed, attempted')
    print TestResults(*doctest.testmod())

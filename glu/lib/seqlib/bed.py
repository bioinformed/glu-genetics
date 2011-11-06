# -*- coding: utf-8 -*-

from __future__ import division

__abstract__  = 'Read and process features from BED files'
__copyright__ = 'Copyright (c) 2010, BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.'
__license__   = 'See GLU license for terms by running: glu license'
__revision__  = '$Id$'


import csv

from   collections import defaultdict
from   itertools   import groupby, count

import pysam

from   glu.lib.fileutils import autofile


GOOD,UNALIGNED,TOOSHORT,LOWOVERLAP = range(4)


def merge_features(features):
  features.sort()

  current_start,current_end,current_name = None,None,set()

  for feature_start,feature_end,feature_name in features:
    if current_start is None:
      current_start,current_end,current_name = feature_start,feature_end,set([feature_name])
    elif current_end<feature_start:
      yield current_start,current_end,','.join(sorted(current_name))
      current_start,current_end,current_name = feature_start,feature_end,set([feature_name])
    else:
      current_name.add(feature_name)
      current_end = max(current_end,feature_end)

  if current_start is not None:
    yield current_start,current_end,','.join(sorted(current_name))



def read_features(filename,merge=True):
  features = defaultdict(list)

  if not filename:
    return features

  bed = csv.reader(autofile(filename),dialect='excel-tab')

  generic_names = ('feature_%06d' % i for i in count(1))

  for row in bed:
    n = len(row)
    if n<3 or row[0].startswith('track ') or row[0].startswith('#'):
      continue

    contig,start,end = row[:3]

    if n>3:
      name = row[3]
    else:
      name = next(generic_names)

    features[contig].append( (int(start),int(end),name) )

  for contig in features:
    if merge:
      features[contig] = list(merge_features(features[contig]))
    else:
      features[contig] = sorted(features[contig])

  return features

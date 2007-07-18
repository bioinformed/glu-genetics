from utils import all


class union_find(object):
  def __init__(self, labels=None):
    if labels is not None:
      self.labels = dict( (l,l) for l in labels )
    else:
      self.labels = {}

  def find(self, i, create=True):
    labels = self.labels

    # Assign new label
    if i not in labels:
      if create:
        labels[i] = i
      return i

    # Find exemplar
    j = i
    while labels[j] != j:
      j = labels[j]

    # Apply path compression
    while labels[i] != i:
      labels[i],i = j,labels[i]

    return j

  def has_key(self, i):
    return i in self.labels

  def union(self, i1, i2):
    self.labels[self.find(i1)] = self.find(i2)

  def sets(self):
    return self.setmap().values()

  def setmap(self):
    sets = {}
    for i in self.labels:
      j = self.find(i)
      d = sets.get(j)
      if not d:
        d = sets[j] = set([j])
      d.add(i)
    return sets

  def __contains__(self,seq):
    '''Returns true if all items in seq belong to the same set'''
    seq = iter(seq)
    j = self.find(seq.next(),False)
    return all(j == self.find(i,False) for i in seq)


def build_dupsets(dups):
  uf = union_find()
  for i1,i2 in dups:
    uf.union(i1,i2)
  return uf.sets()


def main():
  dups = [(1,1),(2,2),(3,3),(4,4),(1,2),(1,3)]
  #dups = [(1,2),(1,3)]

  print build_dupsets(dups)


if __name__ == '__main__':
  main()

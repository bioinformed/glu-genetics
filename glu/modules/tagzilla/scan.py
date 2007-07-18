data = range(11)
monitor = 2,6
maxd = 2

def scan(monitor,data,maxd):
  n = len(data)
  monitor = sorted(monitor)

  pos = 0
  for m in monitor:
    while pos < n and data[pos] < m:
      pos += 1
    while pos < n and data[pos]-m <= maxd:
      yield data[pos]
      pos += 1

  pos = n-1
  for m in reversed(monitor):
    while pos >= 0 and data[pos] > m:
      pos -= 1
    while pos >= 0 and m-data[pos] <= maxd:
      yield data[pos]
      pos -= 1


def main():
  print sorted(set(scan( [2,9], range(11), 2 )))


if __name__ == '__main__':
  main()

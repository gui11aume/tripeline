#!/usr/bin/env python
# -*- coding:utf-8 -*-

import collections
import sys

from gzopen import gzopen

# Constants.
no_hits = object()
many_hits = object()
position_unknown = set([no_hits, many_hits])

# Collection of counting dictionaries.
brcd_counter = collections.defaultdict(collections.Counter)
pos_counter = collections.defaultdict(collections.Counter)

def main(mapfile):
   with gzopen(mapfile) as f:
      for line in f:
         items = line.split()
         brcd = items[0].split(':')[0]
         # The character '-' at the end of the line indicates
         # that there is no hit for the sequence.
         if items[-1] == '-':
            brcd_counter[brcd][no_hits] += 1
            continue
         # In case that there are several hits, they will be
         # separated by ",".
         try:
            (loc,) = items[-1].split(',')
         except ValueError:
            brcd_counter[brcd][many_hits] += 1
            continue
         (chrom,strand,pos,ignore) = loc.split(':')
         brcd_counter[brcd][(chrom,pos)] += 1
         pos_counter[(chrom,pos)][brcd] += 1

   # Find the Charlies.
   Charlies = set()
   for (pos,counts) in pos_counter.items():
      try:
         ((brcd1, alpha),(brcd2, beta)) = counts.most_common(2)
      except ValueError:
         continue
      if alpha < 8 * beta:
         Charlies.add(pos)
         print pos, counts

   print '----'

   # Find the Bobs.
   Bobs = set()
   for (brcd,counts) in brcd_counter.items():
      try:
         ((pos1, alpha),(pos2, beta)) = counts.most_common(2)
      except ValueError:
         ((pos1, alpha),) = counts.most_common(1)
         beta = 0
      if alpha < 8 * beta: Bobs.add(brcd)
      elif pos1 not in position_unknown: print brcd, pos1, alpha


if __name__ == '__main__':
   main(sys.argv[1])

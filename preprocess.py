#!/usr/bin/env python
# -*- coding:utf-8 -*-

import sys

from itertools import izip
from gzopen import gzopen

# Params.
min_brcd = 15
max_brcd = 25
min_genome = 15
transposon = 'TTCAACTGTA'
shift = len(transposon)

# FASTQ params
Q0 = 33 # '!' represents the lowest quality

def substr(string, start, length):
   return string[start:start+length]

def score_from_quality(qstring):
   '''Output the 5-th smallest quality value.'''
   return sorted([ord(c) for c in qstring])[4] - Q0

with gzopen(sys.argv[1]) as f, gzopen(sys.argv[2]) as g:
   # Aggregate iterator of f,g iterators -> izip(f,g)
   for lineno,(line1,line2) in enumerate(izip(f,g)):
      # Take only sequence and quality on lines 1 and 3 (mod 4).
      modulo = lineno % 4
      if modulo == 1:
         valid = False
         # Split on "CATG" and take the first fragment.
         # In case there is no "CATG", the barcode will be rejected
         # for being too long.
         brcd = line1.rstrip().split('CATG')[0]
         if not min_brcd < len(brcd) < max_brcd: continue
         try:
            gpos = line2.index(transposon) + shift
         except ValueError:
            continue 
         # Select the region from the end of the transposon to
         # the first "CATG", if any.
         genome = line2[gpos:].split('CATG')[0].rstrip()
         if len(genome) < min_genome: continue
         valid = True
      elif modulo == 3 and valid:
         qbrcd = score_from_quality(substr(line1, 0, len(brcd)))
         qgen = score_from_quality(substr(line2, gpos, len(genome)))
         sys.stdout.write('>%s:%d,%d\n%s\n' % (brcd,qbrcd,qgen,genome))

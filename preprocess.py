#!/usr/bin/env python
# -*- coding:utf-8 -*-

import sys

from itertools import izip
from gzopen import gzopen
from math import sqrt

from automata import PatternMatcher

# Params.
min_brcd = 15
max_brcd = 25
min_genome = 15

# The known parts of the sequences are matched with a Levenshtein
# automaton. On the reverse read, the end of the transposon
# corresponds to a 34 bp sequence ending as shown below. We allow
# up to 3 mismatches/indels. On the forward read, the only known
# sequence is the CATG after the barcode, which is matched exactly.
transposon = PatternMatcher('TGTATGTAAACTTCCGACTTCAACTGTA', 3)

# FASTQ params.
Q0 = 33 # '!' represents the lowest quality.
QM = 74 # 'J' maximum quality.
norm = 1000/float(QM-Q0)

def substr(string, start, length):
   return string[start:start+length]

def score_from_quality(qstring):
   # Geometric mean of the 4 smallest values.
   q = 1
   for num in sorted(qstring)[:4]:
      q *= ord(num) - Q0
   return int(sqrt(sqrt(q))*norm)


def main(fastq1, fastq2):
   with gzopen(fastq1) as f, gzopen(fastq2) as g:
      # Aggregate iterator of f,g iterators -> izip(f,g).
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
            # Use a Levenshtein automaton to find the transpsoson
            # sequence. Genomic position starts next position (equal
            # to 0 in case there is no match).
            gpos = transposon.end(line2) + 1
            if not gpos: continue
            # Select the region from the end of the transposon to
            # the first "CATG", if any.
            genome = line2[gpos:].split('CATG')[0].rstrip()
            if len(genome) < min_genome: continue
            valid = True
         elif modulo == 3 and valid:
            qbrcd = score_from_quality(substr(line1, 0, len(brcd)))
            qgen = score_from_quality(substr(line2, gpos, len(genome)))
            sys.stdout.write('>%s:%d,%d\n%s\n' % (brcd,qbrcd,qgen,genome))


if __name__ == '__main__':
   main(sys.argv[1], sys.argv[2])

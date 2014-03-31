#!/usr/bin/env python
# -*- coding:utf-8 -*-

import unittest
import preprocess

from automata import PatternMatcher

class TestPartialMatcher(unittest.TestCase):
   def test_(self):
      transposon = PatternMatcher('TGTATGTAAACTTCCGACTTCAACTGTA', 3)
      seq1 = 'TGTATGTAAACTTCCGACTTCAACTGTA'
      self.assertEqual(transposon.start(seq1), 0)
      self.assertEqual(transposon.end(seq1), 27)
      seq2 = 'TGTATGTAAACTTCCGACTTCAACTGTAx'
      self.assertEqual(transposon.start(seq2), 0)
      self.assertEqual(transposon.end(seq2), 27)
      seq3 = 'xTGTATGTAAACTTCCGACTTCAACTGTA'
      self.assertEqual(transposon.start(seq3), 1)
      self.assertEqual(transposon.end(seq3), 28)
      seq4 = 'TGTATGTAAACTTCCGxACTTCAACTGTA'
      self.assertEqual(transposon.start(seq4), 0)
      self.assertEqual(transposon.end(seq4), 28)
      seq5 = 'TAAGGTGTATGTAAACTTCCGACTTCAACTGTAGATCCCACAGCGATATCGACGTGA'
      self.assertEqual(transposon.start(seq5), 5)
      self.assertEqual(transposon.end(seq5), 32)

class TestPreprocess(unittest.TestCase):

   def test_score(self):
      qual = "C@@FFFFFDFFBFH@GHIGCGHICFEHBEHGDDHGIGIGIIGCDDGHFHII" \
             "FAEAAG>@EEHG9AAA;);;7;(.6;;=@B6;B-9@BB8AA<:A@CCAA"
      self.assertEqual(preprocess.score_from_quality(qual), 235)

   def test_susbtr(self):
      string = '123456789'
      self.assertEqual(preprocess.substr(string, 1, 3), '234')
      # Counter-intuitive properties of 'substr' (used here only
      # for documentation purposes, should not be used like this).
      self.assertEqual(preprocess.substr(string, 20, 10), '')
      self.assertEqual(preprocess.substr(string, -5, 10), '5')

if __name__ == '__main__':
   unittest.main()

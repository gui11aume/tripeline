#!/usr/bin/env python
# -*- coding:utf-8 -*-

import unittest
import preprocess

class TestPreprocess(unittest.TestCase):

   def test_score(self):
      qual = "C@@FFFFFDFFBFH@GHIGCGHICFEHBEHGDDHGIGIGIIGCDDGHFHII" \
             "FAEAAG>@EEHG9AAA;);;7;(.6;;=@B6;B-9@BB8AA<:A@CCAA"
      self.assertEqual(preprocess.score_from_quality(qual), 21)

   def test_susbtr(self):
      string = '123456789'
      self.assertEqual(preprocess.substr(string, 1, 3), '234')
      # Counter-intuitive properties of 'substr' (used here only
      # for documentation purposes, should not be used like this).
      self.assertEqual(preprocess.substr(string, 20, 10), '')
      self.assertEqual(preprocess.substr(string, -5, 10), '5')

if __name__ == '__main__':
   unittest.main()

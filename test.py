#!/usr/bin/env python

import os
import unittest
import tripeline
import warnings


class TestCollectIntegrations(unittest.TestCase):
   def test_high_level(self):
      if os.path.exists('testcase_insertions.txt'):
         os.unlink('testcase_insertions.txt')
      tripeline.collect_integrations('testcase_starcode.txt','testcase.map',('testcase_gDNA_starcode.txt','testcase_gDNA_spikes_starcode.txt'),('testcase_cDNA_starcode.txt','testcase_cDNA_spikes_starcode.txt'))

                                     
      self.assertTrue(os.path.exists('testcase_insertions.txt'))
      with open('testcase_insertions.txt') as f:
         lines = [line for line in f if line[0] != '#']
      with open('testcase_expected.txt') as f:
         expected = f.readlines()
      self.assertEqual(len(lines), len(expected))
      for (line1, line2) in zip(lines, expected):
         self.assertEqual(line1, line2)

   def test_gDNA_handling(self):
      if os.path.exists('testcase_insertions.txt'):
         os.unlink('testcase_insertions.txt')
      tripeline.collect_integrations('testcase_starcode.txt',
      'testcase.map',('testcase_gDNA_starcode.txt','testcase_gDNA_spikes_starcode.txt'),('testcase_cDNA_starcode.txt','testcase_cDNA_spikes_starcode.txt'))
      
      self.assertTrue(os.path.exists('testcase_insertions.txt'))
      with open('testcase_insertions.txt') as f:
         lines = [line for line in f if line[0] != '#']
      with open('testcase_expected.txt') as f:
         expected = f.readlines()
      self.assertEqual(len(lines), len(expected))
      for (line1, line2) in zip(lines, expected):
         self.assertEqual(line1, line2)

   def test_file_format(self):
      if os.path.exists('testcase_insertions.txt'):
         os.unlink('testcase_insertions.txt')
      with self.assertRaises(tripeline.FormatException) as cm:
         tripeline.collect_integrations('testcase_starcode.txt','testcase.map',('testcase_gDNA_wformat_starcode.txt','testcase_gDNA_spikes_starcode.txt'),('testcase_cDNA_starcode.txt','testcase_cDNA_spikes_starcode.txt'))
      self.assertEqual(str(cm.exception),"Input file with wrong format")
      
      with self.assertRaises(tripeline.FormatException) as cm:
         tripeline.collect_integrations('testcase_starcode.txt','testcase.map',('testcase_gDNA_wformat2_starcode.txt','testcase_gDNA_spikes_starcode.txt'),('testcase_cDNA_starcode.txt','testcase_cDNA_spikes_starcode.txt'))
      self.assertEqual(str(cm.exception),"Input file with wrong format")

      
if __name__ == '__main__':
    unittest.main()


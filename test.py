#!/usr/bin/env python

import os
import unittest
import tripeline
import warnings

# class TestFirstFunction(unittest.TestCase):
   
#    def test_imput_fastq(self):
#       '''Tests if the imput provided is fastq. Should write a Warning if the
#       file is in another format'''
#       #self.assertEqual()
#       return

#    def test_warn_fasta_present(self):
#       """Test that the user is warned when the fasta file will not
#       be generated."""
#       expected_warning = "Fasta already present."
#       if not os.path.exists('p8_RA_7077_ACTTGA_read1_20.fastq.gz'):
#          raise Exception('File not present (cannot perform test).')
#       if not os.path.exists('p8_RA_7077_ACTTGA_read2_20.fastq.gz'):
#          raise Exception('File not present (cannot perform test).')
#       with warnings.catch_warnings(record=True) as wlist:
#          fastq_file_1 = 'p8_RA_7077_ACTTGA_read1_20.fastq.gz'
#          fastq_file_2 = 'p8_RA_7077_ACTTGA_read2_20.fastq.gz'
#          tripeline.extract_reads_from_PE_fastq(fastq_file_1, fastq_file_2)
#       self.assertEqual(len(wlist), 1)
#       self.assertEqual(expected_warning, str(wlist[0].message))

class TestCollectIntegrations(unittest.TestCase):
   def test_high_level(self):
      if os.path.exists('testcase_insertions.txt'):
         os.unlink('testcase_insertions.txt')
      tripeline.collect_integrations('testcase_starcode.txt', 'testcase.map')
                                     
      self.assertTrue(os.path.exists('testcase_insertions.txt'))
      with open('testcase_insertions.txt') as f:
         lines = [line for line in f if line[0] != '#']
      with open('testcase_expected.txt') as f:
         expected = f.readlines()
      self.assertEqual(len(lines), len(expected))
      for (line1, line2) in zip(lines, expected):
         self.assertEqual(line1, line2)

      os.unlink('testcase_insertions.txt')
         


      
if __name__ == '__main__':
    unittest.main()


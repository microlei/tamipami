import math
import subprocess
import argparse
import logging
import tempfile
import os
import re
import itertools
from collections import Counter
import statistics
from typing import Dict, List, Tuple, Literal

import matplotlib.pyplot as plt
import pandas
from Bio import SeqIO
import scipy.stats
from skbio.stats import composition
from statsmodels.stats.multitest import fdrcorrection

import pytest
from TAMIPAMI.main import *

def test_iterate_kmer():
    expected_2mer = {'AA': 0, 'AC': 0, 'AT': 0, 'AG': 0, 'CA': 0, 'CC': 0, 'CT': 0, 'CG': 0, 'TA': 0, 'TC': 0, 'TT': 0, 'TG': 0, 'GA': 0, 'GC': 0, 'GT': 0, 'GG': 0}
    actual_2mer = iterate_kmer(2)
    assert actual_2mer == expected_2mer, "kmers not generated properly!"

test_dir = tempfile.mkdtemp()
def test_merge():
    merge_reads('tests/test_data/NT-10seqs_R2.fastq', 'test/test_data/NT-10seqs_R1.fastq','out.fastq')
    #count_pam('GGAATCCCTTCTGCAGCACCTGG', 'out.fastq', 4, '3prime')

def test_count_pam():
    result = count_pam('GGAATCCCTTCTGCAGCACCTGG', 'tests/test_data/NT-10seqs_merged_correct.fastq', 4, '3prime')
    assert result['GCGC'] == 9
#!/usr/bin/env python3

import matplotlib.pyplot as plt
plt.style.use('seaborn-poster')
from matplotlib import colors

import numpy as np
import pandas as pd
import cooler
import cooltools
import cooltools.expected
import bioframe

import os
import pickle
import sys

merged = None
bIsEmpty = True
strDir = sys.argv[1]
for strFile in os.listdir(strDir):
   if strFile.endswith(".dat"):
      print(strFile)
      cvd = pd.read_pickle(strDir + '/' + strFile)
      if bIsEmpty:
         merged = cvd
         bIsEmpty = False
      else:
         merged = pd.concat([merged, cvd], axis=0)

print(merged.head(5))
print(merged.tail(5))

merged.to_pickle(sys.argv[2])


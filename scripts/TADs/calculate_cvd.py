#!/usr/bin/env python3

import matplotlib.pyplot as plt
plt.style.use('seaborn-poster')
from matplotlib import colors

import numpy as np
import pandas as pd
import cooler
import cooltools
import cooltools.expected
import sys
import bioframe
#import pickle

strCoolerFile = sys.argv[1]
print(f"Loading the cooler file: {strCoolerFile}")
clr = cooler.Cooler(strCoolerFile)

strRegionsFile = sys.argv[2]
print(f"Loading the regions file: {strRegionsFile}")
regions = pd.read_csv(strRegionsFile, sep='\t', header=None, names=['chrom', 'start', 'end'])

print("Calculating the probabilities")
cvd = cooltools.expected.diagsum(clr=clr, regions=regions, transforms={'balanced': lambda p: p['count']*p['weight1']*p['weight2']})

pd.to_pickle(cvd, sys.argv[3])

#with open(sys.argv[3], 'wb') as hFile:
#	pickle.dump({'cooler': strCoolerFile, 'regions': strRegionsFile, 'cvd': cvd}, hFile)

cvd.to_csv(sys.argv[4], index=False)

print(cvd.head(4))
print(cvd.tail(4))

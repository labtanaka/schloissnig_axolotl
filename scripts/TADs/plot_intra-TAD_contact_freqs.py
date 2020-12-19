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
import re
import pdb

print(f'Reading the pickle file : {sys.argv[1]}')
cvd = pd.read_pickle(sys.argv[1])

if not 'tad_len' in cvd.columns:
   print('Assigning the TAD groups')
   tadgroup = []
   tadlen = []
   for i, row in cvd.iterrows():
      name = row['region']
      m = re.match('.+:([0-9]+)-([0-9]+)', name)
      if m:
         l = int(m.group(2)) - int(m.group(1))
         tadlen.append(l)
      else:
         print(f'Unmatched region: {name}')
         sys.exit(1)
      # Regions
      # <100k
      # 100-250k
      # 250-500k
      # 500-750k
      # 750-1000k
      # 1000-5000k
      # 5000-10000k
      # 10000-15000k
      # 15000-20000k
      # 20000-30000k
      # 30000-50000k
      # >50M
      if l < 100_000:
         tadgroup.append('<100k')
      elif l < 250_000:
         tadgroup.append('100k-250k')
      elif l < 500_000:
         tadgroup.append('250k-500k')
      elif l < 750_000:
         tadgroup.append('500k-750k')
      elif l < 1_000_000:
         tadgroup.append('750k-1M')
      elif l < 5_000_000:
         tadgroup.append('1M-5M')
      elif l < 10_000_000:
         tadgroup.append('5M-10M')
      elif l < 15_000_000:
         tadgroup.append('10M-15M')
      elif l < 20_000_000:
         tadgroup.append('15M-20M')
      elif l < 30_000_000:
         tadgroup.append('20M-30M')
      elif l < 50_000_000:
         tadgroup.append('30M-50M')
      else:
         tadgroup.append('>50M')

   print(f'Saving the data to {sys.argv[1]}')
   cvd['tad_len'] = tadlen
   cvd['tad_group'] = tadgroup
   cvd.to_pickle(sys.argv[1])


print('Preparing the data')
cvd_agg = (cvd.groupby(['diag','tad_group']).agg({'n_valid':'sum','count.sum':'sum','balanced.sum':'sum'}).reset_index())
cvd_agg['s_bp'] = (cvd_agg['diag'] * int(sys.argv[2]))
cvd_agg['count.avg'] = (cvd_agg['count.sum'] / cvd_agg['n_valid'])
cvd_agg['balanced.avg'] = (cvd_agg['balanced.sum'] / cvd_agg['n_valid'])

print('Plotting the data')
f, (ax1, ax2) = plt.subplots(1,2, figsize=(30, 10))
for grp in ['<100k', '100k-250k', '250k-500k', '500k-750k', '750k-1M', '1M-5M', '5M-10M', '10M-15M', '15M-20M', '20M-30M', '30M-50M', '>50M']:
   print(f"  Processing the group '{grp}'")
   cvd_grp = cvd_agg[cvd_agg['tad_group'] == grp]
   tmp = cvd[cvd['tad_group'] == grp]
   nElem = len(tmp['region'].unique())
   ax1.loglog(cvd_grp['s_bp'], cvd_grp['balanced.avg'], label=f'{grp}: N={nElem}')

   cvd_grp = cvd[cvd['tad_group'] == grp]
   if cvd_grp.empty:
      continue
   
   try:
      lb_cvd, lb_slopes, lb_distbins = cooltools.expected.logbin_expected(cvd_grp)
      lb_cvd_agg, lb_slopes_agg = cooltools.expected.combine_binned_expected(lb_cvd, binned_exp_slope=lb_slopes)
      lb_cvd_agg['s_bp'] = lb_cvd_agg['diag.avg'] * int(sys.argv[2])
      lb_slopes_agg['s_bp'] = lb_slopes_agg['diag.avg'] * int(sys.argv[2])
      ax2.semilogx(lb_slopes_agg['s_bp'], lb_slopes_agg['slope'], alpha=0.5, label=f'{grp}')
   except:
      print(f"ERROR: An error occured. The slope of '{grp}' was not calculated")

ax1.set(xlabel='separation, bp', ylabel='Intra TAD contact frequency')
ax1.set_aspect(1.0)
ax1.grid(lw=0.5)
ax1.legend()

ax2.set(xlabel='separation, bp', ylabel='slope')
ax2.set_aspect(1.0)
ax2.grid(lw=0.5)
ax2.legend()

f.savefig(sys.argv[3])

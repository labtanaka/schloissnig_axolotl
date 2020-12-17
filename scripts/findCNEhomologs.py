#!/usr/bin/env python3

import os
import numpy as np
import sys
import math
import pandas as pd
import re
import subprocess


def usage():
    print('USAGE:')
    print(f'   python3 {sys.argv[0]} <hg19_vs_ambMex.matches> <human_enhancers.hg19.fa> <AmexG_v6.DD.corrected.round2.chr.fa>')



if len(sys.argv) < 4:
    usage()
    sys.exit(1)


taskid = 1
with(open(sys.argv[1])) as hFile:
    for line in hFile.readlines():
        if line.startswith('Human locus'):
            continue
        
        locusHsap, cnesHsap, genesHsap, locusAmex, genesAmex = line.strip().split('\t')
        if cnesHsap == 'N/A':
            continue

        with(open('_CNE_hs.list', 'w')) as hQueries:
            for m in re.findall('chr[^:]+:[0-9]+-[0-9]+', cnesHsap):
                print(m.replace(':', '_').replace('-', '_'), file=hQueries)

        print(f'Processing {locusHsap}', file=sys.stderr)
        print(f'   Extracting the human sequences', file=sys.stderr)
        os.system(f'samtools faidx -r _CNE_hs.list {sys.argv[2]} > _tmp_/{taskid}_CNE_hs.fa')
        
        print(f'   Extracting the axolotl sequences', file=sys.stderr)
        os.system(f'samtools faidx {sys.argv[3]} {locusAmex} > _tmp_/{taskid}_amex.fa')

        print(f"lastz_32 _tmp_/{taskid}_amex.fa _tmp_/{taskid}_CNE_hs.fa --format=sam W=5 K=2200 L=3000 Q=HoxD55.q > _tmp_/{taskid}_out.sam")

        taskid += 1

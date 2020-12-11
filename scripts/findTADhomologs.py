#!/usr/bin/env python3

import os
import numpy as np
import sys
import math
import pandas as pd
import re


# Load the tads+genes files for human and axolotl separately 
# Format:
#   chr     start    end      gene1;gene2;gene3
def loadTads(strFilename):
    tads = dict()
    with open(strFilename, 'r') as hFile:
        for strLine in hFile.readlines():
            chrID, start, end, genes = strLine.strip().split('\t')
            tads[f'{chrID}:{start}-{end}'] = set()
            for geneID in genes.split(';'):
                for geneID2 in geneID.split('|'):
                    if geneID2 and not re.match('^LOC[0-9]+$', geneID2):
                        tads[f'{chrID}:{start}-{end}'].add(geneID2)
    return (tads)


def loadCNEs(strFilename):
    pddata = []
    with open(strFilename, 'r') as hFile:
        for strLine in hFile.readlines():
            if strLine.startswith('>'):
                chrID, start, end = (strLine.strip()[1:]).split('_')
                pddata.append([chrID, int(start), int(end)])
    return pd.DataFrame(data=pddata, columns=['Chromosome', 'Start', 'End'])


def getTADsforGene(tads):
    genes = dict()
    for coord in tads:
        for geneID in tads[coord]:
            if genes.get(geneID):
                genes[geneID].append(coord)
            else:
                genes[geneID] = [coord]
    a = list(genes.keys())
    return (genes)


def findMatchingTAD(queryGenes, targetGenes):
    tadlist = []
    for geneID in queryGenes:
        tads = targetGenes.get(geneID)
        if tads:
            for tad in tads:
                tadlist.append(tad)
    # Now find the most common TAD
    nMax = 0
    best = None
    for tad in tadlist:
        n = tadlist.count(tad)
        if n > nMax:
            nMax = n
            best = tad
    return(best) 
    

####### Main #######
if len(sys.argv) < 5:
    print("USAGE", file=sys.stderr)
    print(f"python3 {sys.argv[0]} /path/to/axolotl/tads /path/to/human/tads /path/to/human/CNEs /path/to/the/output", file=sys.stderr)
    sys.exit(1)


print("ARGUMENTS:", file=sys.stderr)
print(f"  Axolotl TADS: {sys.argv[1]}", file=sys.stderr) # '/groups/tanaka/Projects/axolotl-genome/AmexG_v6.0/AmexG_v6.0_DD/work/manuscript/compare_TADs/ambMex60DD.putative.tads+genes.bed'
print(f"  Human TADS:   {sys.argv[2]}", file=sys.stderr) # '/groups/tanaka/Projects/axolotl-genome/AmexG_v6.0/AmexG_v6.0_DD/work/manuscript/compare_TADs/hg19.tads+genes'
print(f"  CNEs:         {sys.argv[3]}", file=sys.stderr) # '/groups/tanaka/Projects/axolotl-genome/AmexG_v6.0/AmexG_v6.0_DD/work/manuscript/rebuttal/CNEs/human_enhancers.hg19.fa'
print(f"  Output:       {sys.argv[4]}", file=sys.stderr) # '/groups/tanaka/Projects/axolotl-genome/AmexG_v6.0/AmexG_v6.0_DD/work/manuscript/compare_TADs/hg19_vs_ambMex.matches'


print(f'Loading axolotl TADs', file=sys.stderr)
ambMex_tads = loadTads(sys.argv[1])
print(f'   Loaded {len(ambMex_tads)} TADs', file=sys.stderr)

print(f'Loading human TADs', file=sys.stderr)
hg19_tads = loadTads(sys.argv[2])
print(f'   Loaded {len(hg19_tads)} TADs', file=sys.stderr)

print(f'Loading human CNEs', file=sys.stderr)
hg19_cnes = loadCNEs(sys.argv[3])
print(f'   Loaded {len(hg19_cnes)} CNEs', file=sys.stderr)

# Iterate over the human TADs and find the axolotl TAD with the most genes present
with open(sys.argv[4], 'w') as hFile:
    print("Human locus\tHuman genes\tHuman CNEs\tAxolotl locus\tAxolotl genes", file=hFile)
    ambMex_genes = getTADsforGene(ambMex_tads)
    for tadCoord in hg19_tads:
        print(f"Processing the locus '{tadCoord}'...", file=sys.stderr)
        bestMatch = findMatchingTAD(hg19_tads[tadCoord], ambMex_genes)
        if bestMatch:
            m = re.search('^([^:]+):([0-9]+)-([0-9]+)$', tadCoord)
            if m:
                cnes = hg19_cnes[hg19_cnes['Start'].between(int(m.group(2)), int(m.group(3))) | hg19_cnes['End'].between(int(m.group(2)), int(m.group(3)))].values.tolist()
                cnes = [f'{x[0]}:{x[1]}-{x[2]}' for x in cnes]
                if len(cnes) == 0:
                    cnes = 'N/A'
                print(f'{tadCoord}\t{hg19_tads[tadCoord]}\t{cnes}\t{bestMatch}\t{ambMex_tads[bestMatch]}', file=hFile)
            else:
                print(f"ERROR: malformatted coordinate '{tadCoord}' ", file=sys.stderr)
                sys.exit(2)

print("Done", file=sys.stderr)

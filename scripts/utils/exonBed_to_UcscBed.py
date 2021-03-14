#!/usr/env/bin python3

import sys

genes = dict()

for line in sys.stdin:
    strChr, strStart, strEnd, strName, strDummy, strStrand = line.strip().split('\t')

    exons = genes.get(strName)
    if not exons:
        exons = []
        genes[strName] = exons
    exons.append([strChr, int(strStart), int(strEnd), strName, strDummy, strStrand])

for strName in genes:
    exons = sorted(genes[strName], key = lambda l: l[1])

    iStart = exons[0][1]
    iEnd = exons[-1][2]

    offsets = []
    lengths = []
    for exon in exons:
        offsets.append(str(exon[1] - iStart))
        lengths.append(str(exon[2] - exon[1] + 1))

    print(f'{strChr}\t{iStart - 1}\t{iEnd}\t{strName}\t{strDummy}\t{strStrand}\t{iStart - 1}\t{iStart - 1}\t0,0,0\t{len(exons)}\t{",".join(lengths)}\t{",".join(offsets)}')
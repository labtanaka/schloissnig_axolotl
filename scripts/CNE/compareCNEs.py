#!/usr/bin/env python3

import os
import sys
import re
import argparse
import pandas as pd


def loadGenes(bedFile):
    data = dict()
    with(open(bedFile, 'r')) as hFile:
        for line in hFile.readlines():
            chrName, start, end, symbol = line.strip().split('\t')
            start = int(start)
            end = int(end)

            # Keep the longest isoform
            if ( not data.get(symbol) or (end - start) > (data[symbol]['end'] - data[symbol]['start']) ):
                for s in symbol.split('|'):
                    s = s.replace(' [nr]', '').replace(' [hs]', '')
                    data[s] = {'chr': chrName, 'start': start, 'end': end}
    return pd.DataFrame(data=[[data[x]['chr'], data[x]['start'], data[x]['end'], x] for x in data], columns=['chr', 'start', 'end', 'symbol'])


def loadTADs(bedFile):
    data = []
    with(open(bedFile, 'r')) as hFile:
        for line in hFile.readlines():
            if line.startswith('#'):
                continue
            chrName, start, end, _ = line.split('\t', 3)
            data.append([chrName, int(start), int(end)])
    return pd.DataFrame(data=data, columns=['chr', 'start', 'end'])


def loadCNEs(samFile):
    data = []
    with(open(samFile, 'r')) as hFile:
        for line in hFile.readlines():
            cne_hg, _, chr_amex, start_amex, _ = line.split('\t', 4)
            chr_hg, start_hg, end_hg = cne_hg.split('_')
            data.append({'chr_hg': chr_hg, 'start_hg': int(start_hg), 'end_hg': int(end_hg), 'chr_amex': chr_amex, 'start_amex': int(start_amex)})
    return data


def main():
    ap = argparse.ArgumentParser(description='This script compares the CNEs and their closest genes within the same TAD between human and axolotl')

    ap.add_argument('--hg19', type=str, required=True, nargs=2, help='paths to the human genes and TADs BED files')
    ap.add_argument('--amex', type=str, required=True, nargs=2, help='path to the axolotl genes and TADs in BED file')
    ap.add_argument('--sam', type=str, required=True, nargs=1, help='lastz SAM file')
    
    try:
        opts = ap.parse_args()
    except:
        ap.print_help()
        sys.exit(0)


    print(f"Reading the human TADs from '{opts.hg19[1]}'", file=sys.stderr)
    tads_hg = loadTADs(opts.hg19[1])
    print(f"Reading the human genes from '{opts.hg19[0]}'", file=sys.stderr)
    genes_hg = loadGenes(opts.hg19[0])
    print(f'  Loaded {genes_hg["chr"].count()} genes in {tads_hg["chr"].count()} TADs', file=sys.stderr)

    print(f"Reading the axolotl TADs from '{opts.amex[1]}'", file=sys.stderr)
    tads_amex = loadTADs(opts.amex[1])
    print(f"Reading the axolotl genes from '{opts.amex[0]}'", file=sys.stderr)
    genes_amex = loadGenes(opts.amex[0])
    print(f'  Loaded {genes_amex["chr"].count()} genes in {tads_amex["chr"].count()} TADs', file=sys.stderr)

    print(f"Reading the CNEs from '{opts.sam[0]}'", file=sys.stderr)
    cnes = loadCNEs(opts.sam[0])
    print(f'  Loaded {len(cnes)} CNEs', file=sys.stderr)

    # Find the closest gene for each CNE in both human and axolotl
    nSkipped = 0
    nValid = 0
    for cne in cnes:
        # First, find the human TAD the current CNE is contained in
        tads = tads_hg[(tads_hg['chr'] == cne['chr_hg']) & 
                       (tads_hg['start'] <= cne['start_hg']) & 
                       (tads_hg['end'] >= cne['end_hg'])]
        if tads['chr'].count() == 0:
            print(f"WARNING: there is no TAD for {cne}. Skipped", file=sys.stderr)
            nSkipped += 1
        else:
            nValid += 1

            # Use the longest TAD
            maxlen = pd.DataFrame.max(tads['end'] - tads['start'])
            _hg_tad = tads[tads['end'] - tads['start'] == maxlen]
            [[chrID, start, end]] = _hg_tad.values

            # Find the human genes within the tad and their axolotl counterparts
            _genes_hg = genes_hg[ (genes_hg['chr'] == chrID) & (genes_hg['start'].between(start, end) | genes_hg['end'].between(start, end)) ]
            _genes_amex = genes_amex[genes_amex['symbol'].isin(_genes_hg['symbol'].values)]

            # Find the TADs in axolotl that contain the genes identified above
            _tads = pd.DataFrame()
            for idx in _genes_amex.index:
                _tads_amex = tads_amex[ (tads_amex['chr'] == _genes_amex['chr'][idx]) & 
                                        (tads_amex['start'] <= _genes_amex['start'][idx]) & 
                                        (tads_amex['end'] >= _genes_amex['end'][idx])]
                # Find the longest TAD for each gene
                maxlen = pd.DataFrame.max(_tads_amex['end'] - _tads_amex['start'])
                _tads_amex = _tads_amex[_tads_amex['end'] - _tads_amex['start'] == maxlen].head(1)
                _tads_amex['symbol'] = _genes_amex['symbol'][idx]
                _tads = pd.concat([_tads, _tads_amex])
            # Take the TAD with the highest number of genes in it.
            gene_counts = _tads.groupby(['chr', 'start', 'end']).count()
            mxN = pd.DataFrame.max(gene_counts['symbol'])
            _amex_tad = gene_counts[gene_counts['symbol'] == mxN].head(1)
            
            # Find the gene closest to the CNE
            mindist = pd.DataFrame.min( pd.DataFrame.abs(genes['start'] - cne['start_hg']) )
            closest = genes[ abs(genes['start'] - cne['start_hg']) == mindist ]

            # Find the homologous axolotl TAD

            
            sys.exit(0)


    print(f"Analyzed {nValid} valid CNEs. Skipped {nSkipped} CNEs", file=sys.stderr)

if __name__ == "__main__":
    main()
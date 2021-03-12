#!/usr/bin/env python3

"""

This script converts the coordinates within the v2 contigs into the coordinates within the
chromosome-scale scaffolds.

"""

import sys
import os
import argparse
import re
import gzip


SPACER='n' * 20


def loadContigSizes(strFilename):
    contigSizes = dict()
    with open(strFilename, "r") as hFile:
        for strLine in hFile:
            arrFields = strLine.split("\t")
            contigSizes[arrFields[0]] = int(arrFields[1])
    return (contigSizes)


def loadScaffoldStructure(strFilename, contigSizes):
    contigs = dict()
    scaffolds = dict()
    offset = 0
    strChr = None
    hFile = None
    if strFilename.endswith('.gz'):
        hFile = gzip.open(strFilename, 'rt')
    else:
        hFile = open(strFilename, "r")
    for strLine in hFile:
        strLine = strLine.strip()
        if strLine.startswith('>'):
            strChr = strLine[1:]
            offset = 0
        else:
            strand, contigID = strLine.split(" ")
            contigs[contigID] = {'start': offset,
                                 'end': offset + contigSizes[contigID],
                                 'strand': strand,
                                 'chr': strChr}
            offset = contigs[contigID]['end'] + len(SPACER)
            scaffolds[strChr] = offset
    hFile.close()
    return (contigs, scaffolds)


def convertCoordinates(strFilename, contigs, scaffolds):
    with open(strFilename, "r") as hFile:
        for strLine in hFile:
            # Skip any commented lines
            if strLine.startswith('#'):
                print(strLine)
                continue

            arrFields = re.split("\t", strLine.strip())

            # If the contig is not a part of any scaffold, output the data unchanged
            props = contigs.get(arrFields[0])
            if not props:
                print(strLine)
                continue

            # There are four possible combinations of the orientations of the contig within
            # the scaffold and the region within the contig:
            # 
            # Contig within scaffold    |   Region within contig    |   Output strand   |   Output start coordinate
            # --------------------------+---------------------------+-------------------+---------------------------
            #           +               |           +               |         +         |       start + offset
            #           +               |           -               |         -         |       start + offset
            #           -               |           +               |         -         |       contig_end - start
            #           -               |           -               |         +         |       contig_end - start
            
            # Thus, if the contig itself is in + orientation within the scaffold, then simply add 
            # the coordinates to the start position of the contig.
            iFrom = None
            iTo = None
            if props['strand'] == '+':
                iFrom = props['start'] + int(arrFields[1])
                iTo = props['start'] + int(arrFields[2])
            else:
                # If the contig if on the opposite strand, the coordinates need to
                # be re-calculated relatively to the end of the contig
                iFrom = props['end'] + 1 - int(arrFields[1])
                iTo  = props['end'] + 1 - int(arrFields[2])

            if iTo < iFrom:
                tmp = iTo
                iTo = iFrom
                iFrom = tmp

            arrFields[0] = props['chr']
            arrFields[1] = str(iFrom)
            arrFields[2] = str(iTo)

            if len(arrFields) >= 6:
                if props['strand'] == arrFields[5]:
                    arrFields[5] = '+'
                else:
                    arrFields[5] = '-'

            print("\t".join(arrFields))


def main():
    ap = argparse.ArgumentParser(description='This script merges RepeatMasker output files for individual contigs into annotation for chromosomes')

    ap.add_argument('--chr', type=str, required=True, nargs=1, help='path to the chromosome structure file. Each line is a contig name, + or - at '+
                                                                    'the beginning of the line indicate the orientation, e.g. + C0001')
    ap.add_argument('--sizes', type=str, required=True, nargs=1, help='path to the chromosome sizes file, e.g. a fai file')
    ap.add_argument('--bed', type=str, required=True, nargs=1, help='path to the input BED file')
    
    try:
        opts = ap.parse_args()
    except:
        ap.print_help()
        sys.exit(0)

    # Load contig sizes
    contigSizes = loadContigSizes(opts.sizes[0])

    # Create scaffold structure
    contigs, scaffolds = loadScaffoldStructure(opts.chr[0], contigSizes)

    # Parse the input *.bed file and convert the coordinates
    convertCoordinates(opts.bed[0], contigs, scaffolds)

if __name__ == "__main__":
    main()
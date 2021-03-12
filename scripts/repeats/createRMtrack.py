#!/usr/bin/env python3

import sys
import os
import argparse
import re


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
    with open(strFilename, "r") as hFile:
        for strLine in hFile:
            strLine = strLine.strip()
            if strLine.startswith('>'):
                strChr = strLine[1:]
                offset = 0
            else:
                strand, contigID = strLine.split(" ")
                contigID = "%s_pilon_pilon" % (contigID)
                contigs[contigID] = {'start': offset,
                                     'end': offset + contigSizes[contigID],
                                     'strand': strand,
                                     'chr': strChr}
                offset = contigs[contigID]['end'] + len(SPACER)
                scaffolds[strChr] = offset
    return (contigs, scaffolds)


def parseOutFile(strFilename, contigs, scaffolds):
    bPrintHeader = True
    with open(strFilename, "r") as hFile:
        for strLine in hFile:
            arrFields = re.split("\s+", strLine.strip())
            if len(arrFields) < 15:
                if bPrintHeader:
                    print(strLine, end="")
                continue
            props = contigs.get(arrFields[4])
            if props:
                bPrintHeader = False
                # If the contig itself is in + orientation within the scaffold, then
                # simply add the coordinates to the start position of the contig.
                iFrom = None
                iTo = None
                if props['strand'] == '+':
                    iFrom = props['start'] + int(arrFields[5])
                    iTo = props['start'] + int(arrFields[6])
                else:
                    # If the contig if on the opposite strand, the coordinates need to
                    # be re-calculated relatively to the end of the contig
                    iFrom = props['end'] + 1 - int(arrFields[6])
                    iTo  = props['end'] + 1 - int(arrFields[5])

                print("%s\t%s\t%s\t%s\t%s\t%d\t%d\t(%d)\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % 
                        (arrFields[0], arrFields[1], arrFields[2], arrFields[3],
                        props['chr'], iFrom, iTo, scaffolds[props['chr']] - iTo, 
                        arrFields[8], arrFields[9], arrFields[10], arrFields[11], arrFields[12], arrFields[13], arrFields[14]))


def main():
    ap = argparse.ArgumentParser(description='This script merges RepeatMasker output files for individual contigs into annotation for chromosomes')

    ap.add_argument('--chr', type=str, required=True, nargs=1, help='path to the chromosome structure file. Each line is a contig name, + or - at '+
                                                                    'the beginning of the line indicate the orientation, e.g. + C0001')
    ap.add_argument('--sizes', type=str, required=True, nargs=1, help='path to the chromosome sizes file, e.g. a fai file')
    ap.add_argument('--rm_out', type=str, required=True, nargs=1, help='path to the RepeatMasker output file')
    
    try:
        opts = ap.parse_args()
    except:
        ap.print_help()
        sys.exit(0)

    # Load contig sizes
    contigSizes = loadContigSizes(opts.sizes[0])

    # Create scaffold structure
    contigs, scaffolds = loadScaffoldStructure(opts.chr[0], contigSizes)

    # Parse the *.out file and write the BED9 file
    parseOutFile(opts.rm_out[0], contigs, scaffolds)

if __name__ == "__main__":
    main()
#!/usr/bin/env python3

import sys
import os
import argparse
import gzip
import re
from time import gmtime, strftime
import pdb

"""
Parses the files downloaded from, for example, from 
  ftp://ftp.ensembl.org/pub/release-102/gff3/lepisosteus_oculatus/Lepisosteus_oculatus.LepOcu1.102.abinitio.gff3.gz
and
  ftp://ftp.ensembl.org/pub/release-102/fasta/lepisosteus_oculatus/pep/Lepisosteus_oculatus.LepOcu1.pep.all.fa.gz
and merges them into a single FASTA file, which contains the protein ID and the position in the genome
"""


def parseGFF3(gff3):
    genes = dict()
    with(gzip.open(gff3, 'rt')) as hFile:
        for line in hFile.readlines():
            if line.startswith('#'):
                continue

            strChrName, _, strType, iStart, iEnd, _, strStrand, _, attr = line.strip().split('\t', 9)
            if strType != 'gene':
                continue

            m = re.search('ID=gene:([^;]+);', attr)
            if m:
                if genes.get(m.group(1)):
                    print(f"   WARN: the gene ID '{m.group(1)}' is duplicated and will be overwritten", file=sys.stderr)
                genes[m.group(1)] = {'chr': strChrName, 'start': int(iStart), 'end': int(iEnd), 'strand': strStrand}
                m2 = re.search(';Name=([^;]+);', attr)
                if m2:
                    genes[m.group(1)]['symbol'] = m2.group(1)

    return genes


def parsePeptides(fasta, genes):
    peptides = dict()
    currentEntry = None
    with(gzip.open(fasta, 'rt')) as hFile:
        for line in hFile.readlines():
            if line.startswith('>'):
                m = re.search('^>([^ ]+).+gene:([^ ]+)\.[0-9]+ transcript:([^ ]+) gene_biotype:protein_coding', line)
                if m:
                    if peptides.get(m.group(1)):
                        print(f"   WARN: the peptide entry '{m.group(1)}' is already present and will be overwritten", file=sys.stderr)
                    peptides[m.group(1)] = {'sequence': '', 'transcript': m.group(3)}
                    ge = genes.get(m.group(2))
                    if not ge:
                        print(f"   WARN: the gene name '{m.group(2)}' was not found in the GFF3 file", file=sys.stderr)
                    else:
                        peptides[m.group(1)]['gene'] = ge
                        peptides[m.group(1)]['gene']['id'] = m.group(2)
                    currentEntry = peptides[m.group(1)]
                else:
                    print(f"   WARN: The entry '{line.strip()}' does not define a valid protein-coding gene. Skipped", file=sys.stderr)
                    currentEntry = None
            elif currentEntry:
                currentEntry['sequence'] += line.strip()
    return peptides


def outputFasta(strFileOut, peptides, strPrefix):
    hFile = None
    if strFileOut:
        hFile = open(strFileOut, 'w')
    for id in peptides:
        ge = peptides[id].get('gene')
        if ge:
            print(f">{strPrefix}_{id} gene_id:{ge['id']};transcript_id:{peptides[id]['transcript']};" + 
                      f"locus:{ge['strand']}{ge['chr']}:{ge['start']}-{ge['end']}\n" + 
                      peptides[id]['sequence'], file=(hFile if hFile else sys.stdout))
    if hFile:
        hFile.close()


def main():
    ap = argparse.ArgumentParser(description='Merges the Ensembl peptide and gff3 files into a single file')

    ap.add_argument('--pep', type=str, required=True, nargs=1, help='path to the FASTA file contaning the peptide sequences')
    ap.add_argument('--gff3', type=str, required=True, nargs=1, help='path to the gff3 file')
    ap.add_argument('--prefix', type=str, required=True, nargs=1, help='prefix to append to the protein IDs')
    ap.add_argument('--out', type=str, required=False, nargs=1, help='path to the output file (default STDOUT')

    try:
        opts = ap.parse_args()
    except:
        ap.print_help()
        sys.exit(0)

    # Check if the files exist
    if not os.path.isfile(opts.pep[0]):
        printLog("The file %s does not exist" % (opts.dna[0]), 0, "error")
        sys.exit(1)
    if not os.path.isfile(opts.gff3[0]):
        printLog("The file %s does not exist" % (opts.blastout[0]), 0, "error")
        sys.exit(1)

    print(f"Reading the GFF3 file '{opts.gff3[0]}'", file=sys.stderr)
    genes = parseGFF3(opts.gff3[0])
    print(f"   Loaded {len(genes)} entries", file=sys.stderr)

    print(f"Reading the peptides file '{opts.pep[0]}'", file=sys.stderr)
    peptides = parsePeptides(opts.pep[0], genes)
    print(f"   Loaded {len(peptides)} entries", file=sys.stderr)

    print(f"Saving the data", file=sys.stderr)
    outputFasta(opts.out[0] if opts.out else None, peptides, opts.prefix[0])


if __name__ == "__main__":
    main()
#!/usr/bin/env python3

import os
import sys
import re


samBody = []

for f in os.listdir(sys.argv[1]):
    if f.endswith('_amex.fa'):
        print(f'Processing the file {f}', file=sys.stderr)

        chrName = None
        iOffset = -1

        with(open(f'{sys.argv[1]}/{f}', 'r')) as hFile:
            line = hFile.readline()
            if not line:
                continue

            m = re.search('^>([^:]+):([0-9]+)', line)
            if not m:
                print(f'   ERROR: the header line is invalid: {line}', file=sys.stderr)
                sys.exit(1)

            chrName = m.group(1)
            iOffset = int(m.group(2))

        with(open(f'{sys.argv[1]}/{f.replace("_amex.fa", "_out.sam")}', 'r')) as hFile:
            for line in hFile.readlines():
                if line.startswith('@'):
                    continue
                fields = line.strip().split('\t')
                if fields[2] != chrName:
                    print(f'   ERROR: the chromosome ID in the SAM file {fields[2]} does not match the chromosome name in the fasta file {chrName}', file=sys.stderr)
                    sys.exit(2)

                fields[3] = str(int(fields[3]))
                print("\t".join(fields))


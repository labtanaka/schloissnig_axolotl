# Convert individual *.out files for contigs into a single *.out file for chromosomes
python3 createRMtrack.py --chr chr.txt --sizes AmexG_v6.DD.corrected.round2.fa.fai --rm_out masked.out > masked.chr.out

# Convert the *.out file to a tab-separated file suitable for the Genome Browser
hgLoadOut -tabFile=repeats.tab -nosplit test masked.chr.out

# Split the *.tab file into separate files: one per repeat class
mkdir rmcls; sort -k12,12 repeats.tab | splitFileByColumn -ending=tab -col=12 -tab stdin rmcls

# For some weird reason, DNA and LINE classes also have a class with the question mark. Paste those together
cd rmcls
cat LINE\?.tab >> LINE.tab
rm LINE\?.tab
cat DNA\?.tab >> DNA.tab
rm DNA\?.tab

# Convert *.tab to *.bed
for FILE in $(find . -name "*.tab"); do echo $FILE; ./toBed6+10.pl $FILE | sort -k1,1 -k2,2n > $(basename $FILE .tab).bed; done

# Convert *.bed to *.bigBed
for FILE in $(find . -name "*.bed"); do echo $FILE; bedToBigBed -tab -type=bed6+10 -as=rmsk16.as.txt $FILE ambMex60DD.chrom.sizes $(basename $FILE .bed).bb; done


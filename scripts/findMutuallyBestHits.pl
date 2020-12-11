#!/usr/bin/env perl

use strict;
use Storable;
use Getopt::Long;

my $strOrg1 = undef;
my $strOrg2 = undef;
my $strOutput = undef;

sub usage {
    print STDERR "USAGE:\nperl $0 --org1 /path/to/blast/org1-org2 --org2 /path/to/blast/org2-org1 --out /path/to/output/file\n";
}

GetOptions('org1=s' => \$strOrg1, 'org2=s' => \$strOrg2, 'out:s' => \$strOutput);

if(!$strOrg1 || !$strOrg2) {
    usage();
    exit(0);
}

print STDERR "Loading organism1 vs organism2: '$strOrg1'\n";
my $refData1 = loadData($strOrg1);
print STDERR "  Loaded " . (scalar keys %{$refData1}) . " blast hits\n";

print STDERR "Loading organism2 vs organism1: '$strOrg2'\n";
my $refData2 = loadData($strOrg2);
print STDERR "  Loaded " . (scalar keys %{$refData2}) . " blast hits\n";

print STDERR "Identifying mutually best hits\n";
open(OUT, ">" . (($strOutput) ? $strOutput : '-'));
foreach my $strID (keys %{$refData1}) {
    my $strHitID = $refData1->{$strID}->{hit};
    if ( ($refData2->{$strHitID}) && ($refData2->{$strHitID}->{hit} eq $strID) ) {
        print OUT sprintf("%s\t%s:%s-%s\t%s\t%s:%s-%s\n", 
                             $strID, $refData2->{$strHitID}->{chr}, $refData2->{$strHitID}->{start}, $refData2->{$strHitID}->{end},
                             $strHitID, $refData1->{$strID}->{chr}, $refData1->{$strID}->{start}, $refData1->{$strID}->{end});
    }
}
close(OUT);


sub loadData {
    my ($strPath) = @_;
    my $refBlastHits = {};
    my $refData = retrieve($strPath);
    foreach my $key (sort keys %{$refData}) {
        # Normally, the query ID is <prefix>_<queryID>. Thus, trim the organism specifix prefix. 
        # Only the axolotl queries are not built using this pattern and have to be handled slightly 
        # differently.
        my $strQuery = undef;
        if (($key =~m/^[a-zA-Z]{6}_(.+)$/) || ($key =~ m/^(AMEX[^;]+)/)) {
            $strQuery = $1;
        } 
        else {
            print STDERR "  WARN: Malformed query: '$key'. Skipped\n";
            next;
        }

        # Only take the top hit. This may not be very clean, but good enough for most
        # cases. I have rearly seen that the best Blast hit is not the first to be output.
        my @arrHits = @{$refData->{$key}};
        my $strHit = $arrHits[0]->{hit};
        my $strHitDef = $arrHits[0]->{hDef};

        if ($strHitDef eq 'No definition line') {  # Axolotl blast database
            if ($strHit =~ m/^(AMEX[^;]+);([^:]+):([0-9]+)-([0-9]+)$/) {
                $refBlastHits->{$strQuery} = {'hit' => $1, 'chr' => $2, 'start' => $3, 'end' => $4};
            }
        } 
        elsif ($strHit =~ m/^[a-zA-Z]{6}_(.+)$/) {
            my $strHitID = $1;
            if ($strHitDef =~ m/^gene_id:([^;]+);transcript_id:([^;]+);locus:([-+])([^:]+):([0-9]+)-([0-9]+)$/) {
                $refBlastHits->{$strQuery} = {'hit' => $strHitID, 'gene' => $1, 'transcript'=> $2, 'strand' => $3, 'chr' => $4, 'start' => $5, 'end' => $6};
            }
        } 
        else {
            print STDERR "  WARN: Malformed hit line: '$strHit'. Skipped\n";
            next;
        }
    }
    return $refBlastHits;
} 
#! usr/bin/perl -w

use strict;

#A program to parse genbank files and extract the sense and antisense CDS names and chromosome locations in order to determine TPMs of the various genes.

#Declaring variables:
my @everything;
my $sense_start = 0;
my $sense_end = 0;
my $sense_product;
my $sense_locus;
my $cds_length = 0;
my $cds_length_total = 0;
my $trnarrna_length = 0;
my $trna_length_total = 0;
my $ncrna_length = 0;
my $ncrna_length_total = 0;
my $rrna_length_total = 0;
my $old_locus = 0;

open INFILE, "/media/shaun/My Passport/Coxiella_RNA_SEQ/sequence.gb" or die "Can't do this stupid thing!\n";
open OUTFILE6, ">/media/shaun/My Passport/Coxiella_RNA_SEQ/CoxiellaNM2_genbank_locus_text" or die "Can't make genbank file\n";

while (<INFILE>) {
   push @everything, $_;
}

close INFILE;

I: foreach my $i (@everything) {
   if ($i =~ /\s+CDS\s+(\d+)..(\d+)/) { #if sense CDS
	$sense_start = $1;
	$sense_end = $2;
	$cds_length = $sense_end - $sense_start + 1;
	$cds_length_total += $cds_length;
   }
   if ($sense_start == 0) {
   if ($i =~ /\s+CDS\s+complement\((\d+)..(\d+)/) { #if antisense CDS
	$sense_start = $1;
	$sense_end = $2;
        print "$1\t$2\n";
	$cds_length = $sense_end - $sense_start + 1;
	$cds_length_total += $cds_length;
   }
   }
   if ($sense_start == 0) {
   if ($i =~ /\s+tRNA\s+(\d+)..(\d+)/) { #if sense tRNA
	$sense_start = $1;
	$sense_end = $2;
	$cds_length = $sense_end - $sense_start + 1;
	$cds_length_total += $cds_length;
   }
   }
   if ($sense_start == 0) {
   if ($i =~ /\s+tRNA\s+complement\((\d+)..(\d+)/) { #if antisense tRNA
	$sense_start = $1;
	$sense_end = $2;
	$cds_length = $sense_end - $sense_start + 1;
	$cds_length_total += $cds_length;
   }
   }
   if ($sense_start == 0) {
   if ($i =~ /\s+rRNA\s+(\d+)..(\d+)/) { #if sense rRNA
	$sense_start = $1;
	$sense_end = $2;
	$cds_length = $sense_end - $sense_start + 1;
	$cds_length_total += $cds_length;
   }
   }
   if ($sense_start == 0) {
   if ($i =~ /\s+rRNA\s+complement\((\d+)..(\d+)/) { #if antisense rRNA
	$sense_start = $1;
	$sense_end = $2;
	$cds_length = $sense_end - $sense_start + 1;
	$cds_length_total += $cds_length;
   }
   }
   if ($sense_start == 0) {
   if ($i =~ /\s+ncRNA\s+(\d+)..(\d+)/) { #if sense ncRNA
	$sense_start = $1;
	$sense_end = $2;
	$cds_length = $sense_end - $sense_start + 1;
	$cds_length_total += $cds_length;
   }
   }
   if ($sense_start == 0) {
   if ($i =~ /\s+ncRNA\s+complement\((\d+)..(\d+)/) { #if antisense ncRNA
	$sense_start = $1;
	$sense_end = $2;
	$cds_length = $sense_end - $sense_start + 1;
	$cds_length_total += $cds_length;
   }
   }
   if ($sense_end > 0) {
	$sense_locus = locus($i);
	if ($sense_locus) {	
#         if ($old_locus gt 0) {
#          if ($old_locus eq $sense_locus) {
#            next I;
#          }
#         }	
#FOR HIGHLIGHTING:
   if ($sense_start < 10) {
	print OUTFILE6 "000000$sense_start\t";
	if ($sense_end < 10) {
	print OUTFILE6 "000000$sense_end\t$sense_locus\n";
        $old_locus = $sense_locus;
        $sense_start = 0;
        $sense_end = 0;
        $sense_locus = 0;
        next I;
	}
	if (($sense_end >= 10) && ($sense_end < 100)) {
	print OUTFILE6 "00000$sense_end\t$sense_locus\n";
        $old_locus = $sense_locus;
        $sense_start = 0;
        $sense_end = 0;
        $sense_locus = 0;
        next I;
	}
	if (($sense_end >= 100) && ($sense_end < 1000)) {
	print OUTFILE6 "0000$sense_end\t$sense_locus\n";
        $old_locus = $sense_locus;
        $sense_start = 0;
        $sense_end = 0;
        $sense_locus = 0;
        next I;
	}
	if (($sense_end >= 1000) && ($sense_end < 10000)) {
	print OUTFILE6 "000$sense_end\t$sense_locus\n";
        $old_locus = $sense_locus;
        $sense_start = 0;
        $sense_end = 0;
        $sense_locus = 0;
        next I;
	}
	if (($sense_end >= 10000) && ($sense_end < 100000)) {
	print OUTFILE6 "00$sense_end\t$sense_locus\n";
        $old_locus = $sense_locus;
        $sense_start = 0;
        $sense_end = 0;
        $sense_locus = 0;
        next I;
	}
	if (($sense_end >= 100000) && ($sense_end < 1000000)) {
	print OUTFILE6 "0$sense_end\t$sense_locus\n";
        $old_locus = $sense_locus;
        $sense_start = 0;
        $sense_end = 0;
        $sense_locus = 0;
        next I;
	}
	if ($sense_end >= 1000000) {
	print OUTFILE6 "$sense_end\t$sense_locus\n";
        $old_locus = $sense_locus;
        $sense_start = 0;
        $sense_end = 0;
        $sense_locus = 0;
        next I;
	}
        $old_locus = $sense_locus;
        $sense_start = 0;
        $sense_end = 0;
        $sense_locus = 0;
        next I;
   }	
   if (($sense_start >= 10) && ($sense_start < 100)) {
	print OUTFILE6 "00000$sense_start\t";
	if (($sense_end >= 10) && ($sense_end < 100)) {
	print OUTFILE6 "00000$sense_end\t$sense_locus\n";
        $old_locus = $sense_locus;
        $sense_start = 0;
        $sense_end = 0;
        $sense_locus = 0;
        next I;
	}
	if (($sense_end >= 100) && ($sense_end < 1000)) {
	print OUTFILE6 "0000$sense_end\t$sense_locus\n";
        $old_locus = $sense_locus;
        $sense_start = 0;
        $sense_end = 0;
        $sense_locus = 0;
        next I;
	}
	if (($sense_end >= 1000) && ($sense_end < 10000)) {
	print OUTFILE6 "000$sense_end\t$sense_locus\n";
        $old_locus = $sense_locus;
        $sense_start = 0;
        $sense_end = 0;
        $sense_locus = 0;
        next I;
	}
	if (($sense_end >= 10000) && ($sense_end < 100000)) {
	print OUTFILE6 "00$sense_end\t$sense_locus\n";
        $old_locus = $sense_locus;
        $sense_start = 0;
        $sense_end = 0;
        $sense_locus = 0;
        next I;
	}
	if (($sense_end >= 100000) && ($sense_end < 1000000)) {
	print OUTFILE6 "0$sense_end\t$sense_locus\n";
        $old_locus = $sense_locus;
        $sense_start = 0;
        $sense_end = 0;
        $sense_locus = 0;
        next I;
	}
	if ($sense_end >= 1000000) {
	print OUTFILE6 "$sense_end\t$sense_locus\n";
        $old_locus = $sense_locus;
        $sense_start = 0;
        $sense_end = 0;
        $sense_locus = 0;
        next I;
	}
        $old_locus = $sense_locus;
        $sense_start = 0;
        $sense_end = 0;
        $sense_locus = 0;
        next I;
   }
   if (($sense_start >= 100) && ($sense_start < 1000)) {
	print OUTFILE6 "0000$sense_start\t";
	if (($sense_end >= 100) && ($sense_end < 1000)) {
	print OUTFILE6 "0000$sense_end\t$sense_locus\n";
        $old_locus = $sense_locus;
        $sense_start = 0;
        $sense_end = 0;
        $sense_locus = 0;
        next I;
	}
	if (($sense_end >= 1000) && ($sense_end < 10000)) {
	print OUTFILE6 "000$sense_end\t$sense_locus\n";
        $old_locus = $sense_locus;
        $sense_start = 0;
        $sense_end = 0;
        $sense_locus = 0;
        next I;
	}
	if (($sense_end >= 10000) && ($sense_end < 100000)) {
	print OUTFILE6 "00$sense_end\t$sense_locus\n";
        $old_locus = $sense_locus;
        $sense_start = 0;
        $sense_end = 0;
        $sense_locus = 0;
        next I;
	}
	if (($sense_end >= 100000) && ($sense_end < 1000000)) {
	print OUTFILE6 "0$sense_end\t$sense_locus\n";
        $old_locus = $sense_locus;
        $sense_start = 0;
        $sense_end = 0;
        $sense_locus = 0;
        next I;
	}
	if ($sense_end >= 1000000) {
	print OUTFILE6 "$sense_end\t$sense_locus\n";
        $old_locus = $sense_locus;
        $sense_start = 0;
        $sense_end = 0;
        $sense_locus = 0;
        next I;
	}
        $old_locus = $sense_locus;
        $sense_start = 0;
        $sense_end = 0;
        $sense_locus = 0;
        next I;
   }
   if (($sense_start >= 1000) && ($sense_start < 10000)) {
	print OUTFILE6 "000$sense_start\t";
	if (($sense_end >= 1000) && ($sense_end < 10000)) {
	print OUTFILE6 "000$sense_end\t$sense_locus\n";
        $old_locus = $sense_locus;
        $sense_start = 0;
        $sense_end = 0;
        $sense_locus = 0;
        next I;
	}
	if (($sense_end >= 10000) && ($sense_end < 100000)) {
	print OUTFILE6 "00$sense_end\t$sense_locus\n";
        $old_locus = $sense_locus;
        $sense_start = 0;
        $sense_end = 0;
        $sense_locus = 0;
        next I;
	}
	if (($sense_end >= 100000) && ($sense_end < 1000000)) {
	print OUTFILE6 "0$sense_end\t$sense_locus\n";
        $old_locus = $sense_locus;
        $sense_start = 0;
        $sense_end = 0;
        $sense_locus = 0;
        next I;
	}
	if ($sense_end >= 1000000) {
	print OUTFILE6 "$sense_end\t$sense_locus\n";
        $old_locus = $sense_locus;
        $sense_start = 0;
        $sense_end = 0;
        $sense_locus = 0;
        next I;
	}
        $old_locus = $sense_locus;
        $sense_start = 0;
        $sense_end = 0;
        $sense_locus = 0;
        next I;
   }
   if (($sense_start >= 10000) && ($sense_start < 100000)) {
	print OUTFILE6 "00$sense_start\t";
	if (($sense_end >= 10000) && ($sense_end < 100000)) {
	print OUTFILE6 "00$sense_end\t$sense_locus\n";
        $old_locus = $sense_locus;
        $sense_start = 0;
        $sense_end = 0;
        $sense_locus = 0;
        next I;
	}
	if (($sense_end >= 100000) && ($sense_end < 1000000)) {
	print OUTFILE6 "0$sense_end\t$sense_locus\n";
        $old_locus = $sense_locus;
        $sense_start = 0;
        $sense_end = 0;
        $sense_locus = 0;
        next I;
	}
	if ($sense_end >= 1000000) {
	print OUTFILE6 "$sense_end\t$sense_locus\n";
        $old_locus = $sense_locus;
        $sense_start = 0;
        $sense_end = 0;
        $sense_locus = 0;
        next I;
	}
        $old_locus = $sense_locus;
        $sense_start = 0;
        $sense_end = 0;
        $sense_locus = 0;
        next I;
   }
   if (($sense_start >= 100000) && ($sense_start < 1000000)) {
	print OUTFILE6 "0$sense_start\t";
	if (($sense_end >= 100000) && ($sense_end < 1000000)) {
	print OUTFILE6 "0$sense_end\t$sense_locus\n";
        $old_locus = $sense_locus;
        $sense_start = 0;
        $sense_end = 0;
        $sense_locus = 0;
        next I;
	}
	if ($sense_end >= 1000000) {
	print OUTFILE6 "$sense_end\t$sense_locus\n";
        $old_locus = $sense_locus;
        $sense_start = 0;
        $sense_end = 0;
        $sense_locus = 0;
        next I;
	}
        $old_locus = $sense_locus;
        $sense_start = 0;
        $sense_end = 0;
        $sense_locus = 0;
        next I;
   }
   else {
	print OUTFILE6 "$sense_start\t$sense_end\t$sense_locus\n";
        $old_locus = $sense_locus;
        $sense_start = 0;
        $sense_end = 0;
        $sense_locus = 0;
        next I;
   }
   $old_locus = $sense_locus;
   $sense_start = 0;
   $sense_end = 0;
   $sense_locus = 0;
   next I;
 }
 }
}

$trna_length_total = 0;
$rrna_length_total = 0;
$ncrna_length_total = 0;
$cds_length_total = 0;


close OUTFILE6;

##########################
######SUBROUTINE##########
##########################

sub product {
foreach my $i (@_) {
my $product;
   if ($i =~ /\s+\/product="(.*?)"/) { #Grab the product for each CDS by searching for the first product in the lines following a CDS
	$product = $1;
   }
   if ($product =~ /(.*?)\(/) {
	$product = $1;
   }
}
}

sub locus {
foreach my $i (@_) {
my $locus;
   if ($i =~ /\s+\/old_locus_tag="(.*?)"/) { #This will grab every CDS for the first product in the lines following a locus tag
	$locus = $1;
	return $locus;
   }
   elsif ($i =~ /\s+\/note="(sRNA.*?)"/) {
	$locus = $1;
	return $locus;
   }
   else {
    next;
   }
}
}



#! /usr/bin/perl -w

#A program to calculate the TPM's and the p value between two sets of independent (LCV and SCV) transcriptome data using Student's independent T test

use strict;

#Declare variables:
my @col;
my $base;
my $read;
my %scv_read_hash;
my @scv_read_arr;
my @bases;
my $old_base = 0;
my $igr_start = 0;
my $igr_end = 0;
my %igr;
my %feature;
my %start_cds;
my %cds;
my $scv_start = 0;
my $scv_end = 0;
my %scv_base;
my %scv_cds;
my $length;
my %read_scv_length;
my %read_scv_sum;
my $scv_total = 0;
my $sum = 0;
my @nucleotide;
my @read_value;
my %scv_igr;
my $igr_sum = 0;
my $length_scv_sum = 0;
my $total_scv_read = 0;
my $avg_scv_read_length = 0;
my $scv_rg_rl = 0;
my $each_scv_t = 0;
my $scv_t = 0;
my $feature_length = 0;
my $scv_sum_reads = 0;
my @old_reads;
my $average_read = 0;
my $read_avg = 0;
my $numerator = 0;
my $denominator = 0;
my $TPM = 0;
my %lcv_read_hash;
my @lcv_read_arr;
my $lcv_start = 0;
my $lcv_end = 0;
my %lcv_base;
my %lcv_cds;
my %read_lcv_length;
my %read_lcv_sum;
my $lcv_total = 0;
my %lcv_igr;
my $length_lcv_sum = 0;
my $total_lcv_read = 0;
my $avg_lcv_read_length = 0;
my $lcv_rg_rl = 0;
my $each_lcv_t = 0;
my $lcv_t = 0;
my $lcv_sum_reads = 0;
my @cds_reads;
my $cds;
my $lcv_squared_dev = 0;
my $scv_squared_dev = 0;
my $added_dev = 0;
my $halved_dev = 0;
my $pooled_dev = 0;
my $sub_means = 0;
my $t_value = 0;
my $two_n = 0;
my $root_n = 0;
my $den = 0;
my %t_values;
my $statistic;
my %scv_read_length;
my %scv_read_deviation;
my %lcv_read_deviation;
my %scv_read_mean;
my %lcv_read_mean;
my $scv_weird_avg;
my $scv_standard_dev;
my @scv_square_read_avg;
my @scv_read_value;
my @scv_nucleotide;
my $scv_length = 0;
my $scv_mean = 0;
my %scv_cds_overlap;
my $scv_sum = 0;
my $scv_square_total = 0;
my $scv_square_sum = 0;
my $scv_minus_avg = 0;
my $scv_squared_avg = 0;
my %scv_base_overlap;

#Open transcript file:
open INFILE, "/media/shaun/My Passport/Coxiella_RNA_SEQ/Converted_Vero_SCV_2_NM2_Ambiguous.txt" or die "Can't open numerical file\n";

while (<INFILE>) {
   chomp; #Remove newline characters
   @col = split (/,/, $_); #Split everything separated by a comma space and save in the @col array
   $base = $col[0]; #Place the position on the chromosome into the variable $base
   $read = $col[1]; #Place the number of transcripts for that position on the chromosome into the $read variable
   $scv_read_hash{$base} = $read; #This will create a temporary hash of all scv transcript reads associated with their position on the chromosome
}

close INFILE;

push @scv_read_arr, {%scv_read_hash}; #This will create an array of hashes

open INFILE2, "/media/shaun/My Passport/Coxiella_RNA_SEQ/CoxiellaNM2_genbank_locus_text" or die "Can't open genebank file\n";

while (<INFILE2>) {
   chomp;
   @bases = split (/\t/, $_);
   if ($old_base ne 0) { #if the variable $old_base has a value
     if ($bases[0] gt $old_base) { #if the $old_base variable is not the same as the $base variable we are currently looking at
	   $igr_start = $old_base + 1; #start our igr right after our previous cds
	   $igr_end = $bases[0] - 1; #end our igr right before our next cds
	   $igr{$igr_start} = $igr_end; #take the region between the end and start of CDS' into our IGR hash
	   $feature{$igr_start} = $igr_end; #this will allow us to calculate the T for each feature!
	   $old_base = 0; #set our variables back to zero
	   $igr_start = 0;
	   $igr_end = 0;
	 }
   }
   $start_cds{$bases[0]} = $bases[1]; #This will create a hash with the starting and ending nt as the CDS name and value
   $feature{$bases[0]} = $bases[1]; #this will allow us to calculate the T for each feature!
   $old_base = $bases[1];
   if ($bases[0] =~ /^0+([1-9]\d*)/) {
    $bases[0] = ($1);
   }
   $cds{$bases[0]} = $bases[2]; #This will create a hash with the starting nt as the key and the CDS name as the value
}

close INFILE2;

#For transcripts, we want to analyze each read value and average the values for each transcript and give the base range for that transcript
#We will define a transcript as a range of values all greater than 0

open OUTFILE3, ">/media/shaun/My Passport/Coxiella_RNA_SEQ/Vero_SCV2_Ambiguous_NM2_TPM" or die "Can't make a sense outfile\n";

#For transcripts:
#To calculate the TPM = (sum of the read x avg length of all reads x 10^6)/(length of the read x (sum of all reads in a feature x average read length) / lenght of the feature)
#To calculate the standard deviation = (x - mean)squared/thatmean and all square rooted

#To grab transcripts that are within CDS:
foreach my $i (0..$#scv_read_arr) { #Sorting the array everything is ultimately stored in
   foreach my $key (sort {$a <=> $b} keys %{$scv_read_arr[$i]} ) { #sorting the hash the nucleotide position and reads are stored in
     print "$key\n";
     foreach my $number (sort {$a <=> $b} keys %start_cds) { #This will look at the starting bases for each CDS in the Coxiella chromosome
       if (($key >= $number) && ($key <= $start_cds{$number})) { #This will only select reads that are between the starting and ending bases for each CDS
	if ($scv_read_arr[$i]{$key} == 0.0) { #Select every transcript read value that equals zero
           $scv_start = $nucleotide[0]; #Give the starting base of the transcript
	   $scv_end = $nucleotide[-1]; #The ending base of the transcript
	   if ($scv_total > 0) {
	    $scv_base{$scv_start} = $scv_end; #Creating a hash with the start as the key and the end as the value
	    $scv_cds{$scv_start} = $cds; #Creating a hash with the start as the key and the cds name as the value
	    $length = $scv_end - $scv_start +1; #This will calculate the length of the read
	    $read_scv_length{$scv_start} = $length; #This will create a hash with the start as the key and the length of the read as the value
	    $read_scv_sum{$scv_start} = $sum; #This will create a hash with the start as the key and the sum of the reads as the value
		push @cds_reads, $key;
	    } #Close the if ($w_total > 0) loop
	    @nucleotide = ();
	    $scv_total = 0;	#Empty all variables to start over fresh for the next transcript
	    $sum = 0;
	} #Close the if ($read_arr[$i]{$key} == 0) loop
	if ($key == $number) { #Select every transcript read value that equals the start of the cds
           $scv_start = $nucleotide[0]; #Give the starting base of the transcript
	   $scv_end = $nucleotide[-1]; #The ending base of the transcript
	   if ($scv_total > 0) {
	    $scv_base{$scv_start} = $scv_end; #Creating a hash with the start as the key and the end as the value
	    $scv_cds{$scv_start} = $cds; #Creating a hash with the start as the key and the standard deviation as the value
	    $length = $scv_end - $scv_start +1; #This will calculate the length of the read
	    $read_scv_length{$scv_start} = $length; #This will create a hash with the start as the key and the length of the read as the value
	    $read_scv_sum{$scv_start} = $sum; #This will create a hash with the start as the key and the sum of the reads as the value
            push @cds_reads, $key;
	   } #Close the if ($w_total > 0) loop
	    @nucleotide = ();
	    $scv_total = 0;	#Empty all variables to start over fresh for the next transcript
		$sum = 0;
	} #Close the if ($key == $number) loop
	if ($scv_read_arr[$i]{$key} > 0.0) {
		push (@nucleotide, $key); #This will add every nucleotide in the read as long as it's value greater than zero and within the CDS and has an R transcripts
		push (@read_value, $scv_read_arr[$i]{$key}); #This will add every read value in this range to the array
		$scv_total++; #This will total the reads within this range
		$cds = $cds{$number};
		$sum += $scv_read_arr[$i]{$key}; #This will add the values of each read
	    }#Close the if ($read_arr[$i]{$key} > 0.0)
          }#Close the foreach my $number (sort keys %start_cds)
	else {
	  next;
	}#Close else loop
	}#Close the if ($antisense loop
   }#Close the foreach my $key loop
}#Close the foreach my $i loop

#To grab transcripts that are within IGR:
foreach my $i (0..$#scv_read_arr) { #Sorting the array everything is ultimately stored in
   BASE: foreach my $key (sort {$a <=> $b} keys %{$scv_read_arr[$i]} ) { #sorting the hash the nucleotide position and reads are stored in
     foreach my $number (sort {$a <=> $b} keys %igr) { #This will look at the starting bases for each CDS in the Coxiella chromosome
       if (($key >= $number) && ($key <= $igr{$number})) { #This will only select reads that are between the starting and ending bases for each CDS
#	 foreach my $i (@cds_reads) { #for every read that's already been classified as an intragenic read:
#	      if ($nucleotide[0] eq $i) { #if it is a read that's also in an IGR
#	       next BASE; #skip to the next read
#	      } #close if loop
#	   }#close foreach loop
	if ($scv_read_arr[$i]{$key} == 0.0) { #Select every transcript read value that equals zero
           $scv_start = $nucleotide[0]; #Give the starting base of the transcript
	   $scv_end = $nucleotide[-1]; #The ending base of the transcript
	   if ($scv_total > 0) {
	    $scv_base{$scv_start} = $scv_end; #Creating a hash with the start as the key and the end as the value
	    $scv_igr{$scv_start} = "IGR".$igr_sum; #Creating a hash with the start as the key and the cds name as the value
	    $length = $scv_end - $scv_start +1; #This will calculate the length of the read
	    $read_scv_length{$scv_start} = $length; #This will create a hash with the start as the key and the length of the read as the value
	    $read_scv_sum{$scv_start} = $sum; #This will create a hash with the start as the key and the sum of the reads as the value
	    } #Close the if ($w_total > 0) loop
	    @nucleotide = ();
	    $scv_total = 0;	#Empty all variables to start over fresh for the next transcript
	    $sum = 0;
	} #Close the if ($read_arr[$i]{$key} == 0) loop
	if ($key == $number) { #Select every transcript read value that equals the start of the cds
           $scv_start = $nucleotide[0]; #Give the starting base of the transcript
	   $scv_end = $nucleotide[-1]; #The ending base of the transcript
	   if ($scv_total > 0) {
	    $scv_base{$scv_start} = $scv_end; #Creating a hash with the start as the key and the end as the value
	    $scv_igr{$scv_start} = "IGR".$igr_sum; #Creating a hash with the start as the key and the standard deviation as the value
	    $length = $scv_end - $scv_start +1; #This will calculate the length of the read
	    $read_scv_length{$scv_start} = $length; #This will create a hash with the start as the key and the length of the read as the value
	    $read_scv_sum{$scv_start} = $sum; #This will create a hash with the start as the key and the sum of the reads as the value
	   } #Close the if ($w_total > 0) loop
	    @nucleotide = ();
	    $scv_total = 0;	#Empty all variables to start over fresh for the next transcript
	    $sum = 0;
	} #Close the if ($key == $number) loop
	if ($scv_read_arr[$i]{$key} > 0.0) {
		push (@nucleotide, $key); #This will add every nucleotide in the read as long as it's value greater than zero and within the CDS and has an SCV transcript
		push (@read_value, $scv_read_arr[$i]{$key}); #This will add every read value in this range to the array
		$scv_total++; #This will total the reads within this range
		$igr_sum += 1;
		$sum += $scv_read_arr[$i]{$key}; #This will add the values of each read
	    }#Close the if ($read_arr[$i]{$key} > 0.0)
          }#Close the foreach my $number (sort keys %start_cds)
	else {
	  next;
	}#Close else loop
	}#Close the if ($antisense loop
   }#Close the foreach my $key loop
}#Close the foreach my $i loop

$scv_total = 0;

foreach my $bases (sort keys %scv_base) { #for each scv read:
   $length_scv_sum += $read_scv_length{$bases}; #add the total length of each read
   $total_scv_read++; #calculate the total number of reads mapped from scv's
}

$avg_scv_read_length = $length_scv_sum / $total_scv_read; #find the average read length of scv transcripts

FEATURE: foreach my $bases (sort {$a <=> $b} keys %feature) { #for each feature (both cds and igrs):
  print "$bases\t$feature{$bases}\n";
  READ: foreach my $reads (sort {$a <=> $b} keys %scv_base) { #for each scv read:
       if (($reads >= $bases) && ($reads <= $feature{$bases})) { #This will only select reads that are between the starting and ending bases for each CDS
	 foreach my $i (@old_reads) { #for previously used reads:
	   if ($i == $reads) { #if the read mapping to the feature have already been used:
	      next READ; #go back up to the next read
	   }#close if loop
	 }#close foreach loop
	 $scv_sum_reads += $read_scv_sum{$reads}; #sum all of the reads mapping to that feature (CDS and IGR's)
	 $feature_length = $feature{$bases} - $bases + 1; #this is the lenght of the feature mapping these reads
         $scv_total += 1;
         push @old_reads, $reads;
     } #close the if (($bases <.......
     else{
      if ($scv_total gt 0) {
	  $scv_rg_rl = $scv_sum_reads * $avg_scv_read_length; #the numerator of T, the sum of all reads mapping to a feature times the average read length
          $each_scv_t = $scv_rg_rl/$feature_length; #the T for each of the features (cds and igr's)
          $scv_t += $each_scv_t; #the T for the SCV transcriptome
	  $feature_length = 0; #empty variables
	  $scv_sum_reads = 0;
          $scv_total = 0;
	}#close if ($scv_total gt 0)
     } #close else{
   } #close the READ loop
} #close the foreach my $bases (sort keys %feature)

      if ($scv_total gt 0) {
	  $scv_rg_rl = $scv_sum_reads * $avg_scv_read_length; #the numerator of T, the sum of all reads mapping to a feature times the average read length
          $each_scv_t = $scv_rg_rl/$feature_length; #the T for each of the features (cds and igr's)
          $scv_t += $each_scv_t; #the T for the SCV transcriptome
	  $feature_length = 0; #empty variables
	  $scv_sum_reads = 0;
          $scv_total = 0;
	}#close if ($scv_total gt 0)
print "The T for Q is $scv_t\n";

print OUTFILE3 "#OUTPUT FROM TPM.pl The T for Q is $scv_t\n#Start\tEnd\tLength\tSum of read\tAverage of read\tTPM#\n";

my $old_nt = 0;

NTS: foreach my $nts (sort {$a <=> $b} keys %scv_base) {
   if ($old_nt ne 0) {
     print OUTFILE3 "$old_nt\t$scv_base{$old_nt}\t$read_scv_length{$old_nt}\t$read_scv_sum{$old_nt}\t$average_read\t$TPM\tIntergenic\n";
     $old_nt = 0;
   }
   $average_read = $read_scv_sum{$nts} / $read_scv_length{$nts}; #the average of the read, just to add to the output
   $read_avg = $read_scv_sum{$nts} * $avg_scv_read_length; #the sum of the read times the average of all reads, part of the numeratore for the TPM
   $numerator = $read_avg * 1000000; #finishing up the numerator for the TPM by mutliplying it by 10^6
   $denominator = $read_scv_length{$nts} * $scv_t; #the denominator for the TPM, the length of the read times the T of the SCV transcriptome
   $TPM = $numerator/$denominator; #determining the TPM for this SCV read
   GENES: foreach my $genes (sort {$a <=> $b} keys %cds) {
     if ($nts eq $genes) {
       print OUTFILE3 "$nts\t$scv_base{$nts}\t$read_scv_length{$nts}\t$read_scv_sum{$nts}\t$average_read\t$TPM\t$cds{$nts}\n";
       $old_nt = 0;
       next NTS; 
     }
     if ($nts ne $genes) {
       $old_nt = $nts;
       next GENES;
     }
   }
} #close the foreach loop

if ($old_nt ne 0) {
 print OUTFILE3 "$old_nt\t$scv_base{$old_nt}\t$read_scv_length{$old_nt}\t$read_scv_sum{$old_nt}\t$average_read\t$TPM\tIntergenic\n";
 $old_nt = 0;
}
@old_reads = ();

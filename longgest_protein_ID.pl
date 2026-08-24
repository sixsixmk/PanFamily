#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

my $input_file = "";
my $output_file = "longest_proteins.fa";

GetOptions(
    "input|i=s"  => \$input_file,
    "output|o=s" => \$output_file,
) or die "Usage: perl $0 -i <input_fasta_file> [-o <output_fasta_file>]\n";

unless ( -f $input_file ) {
    die "ERROR: Input file '$input_file' does not exist or was not provided. Please specify a valid FASTA file with -i parameter.\n";
}

open my $in_fh, "<", $input_file or die "Unable to open input file '$input_file': $!\n";
$/=">";
<$in_fh>;
my %hash;

while(<$in_fh>){
    chomp;
    my($id,$seq)=(split /\n/,$_,2)[0,1];
    $seq=~s/\n//g;
    my($ID1,$ID2)=split /\./,$id;   
    my $length=length($seq);
   
    if(exists  $hash{$ID1}){
        if($length > length($hash{$ID1})) {
            $hash{$ID1}=$seq;           
        } else {
            next;
        }
    } else {
        $hash{$ID1}=$seq;
    }
}

close $in_fh;
open my $out_fh, ">", $output_file or die "Unable to open output file '$output_file': $!\n";

foreach my $_ (keys %hash){
    print $out_fh ">$_\n$hash{$_}\n";
}

close $out_fh;

print "Processing completed. The longest sequences have been saved to '$output_file'.\n";

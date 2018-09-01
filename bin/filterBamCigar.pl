#!/usr/bin/perl
use strict;
use warnings;

my $sam = shift @ARGV;
my %reads;

open(SAM,$sam);
while(<SAM>){
	my $currLine = $_;
	if($currLine =~ m/^@/){
		next;
	}
	my @columns = split("\t",$currLine);
	if($columns[5] =~ m/N/){
		$reads{$columns[0]} = 1;
	}
}
close(SAM);

open(SAM,$sam);
while(<SAM>){
	my $currLine = $_;
        if($currLine =~ m/^@/){
		print $currLine;
                next;
        }
        my @columns = split("\t",$currLine);
	if(exists $reads{$columns[0]}){
		print $currLine;
	}
}
close(SAM);

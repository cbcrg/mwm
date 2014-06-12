#!/usr/bin/env perl

##############################################################################################
### fixCoord.pl                                                                            ###
##############################################################################################
### Script to process text files generated by SMART. Sometimes coordinates were not        ###
### correctly introduced when doing the recording. This script corrects this type of files ###
### example fixCoord.pl                                                                    ###
##############################################################################################
###OPTIONS                                                                                 ###
###                                                                                        ###
##############################################################################################

use strict;
use FileHandle;
use Data::Dumper;

my $file = shift (@ARGV);

my $F=new FileHandle;
open ($F, $file);        

my @fileLines; 

my $coordReached = 0;
my $session = "";

while (<$F>)
	{   
		#Skip lines above coordinates
		if ($_ =~ /^( )+\d/)
			{
				chomp;
#				print $_;
				my $aryLine = [ split('\t', $_) ];
#				print Dumper ($aryLine);
#				last;
				print "$aryLine->[3]\t";
				$aryLine->[3] = $aryLine->[3] * 175 / 1475.25; 
				print "$aryLine->[3]\n";		
			}
		else 
			{
#				print $_;
			}
									
	}
#!/usr/bin/env perl

##############################################################################################
###  all2singleSessionsSmart.pl                                                            ###
##############################################################################################
### Script to process a file that have several sessions of a MWM test into single files    ###
### example all2singleSessionsSmart.pl youngG2AllSessions.txt                              ###
##############################################################################################
###OPTIONS                                                                                 ###
###                                                                                        ###
##############################################################################################

use strict;
use FileHandle;
use Data::Dumper;

my $file = shift (@ARGV);

my $F=new FileHandle;
open($F, $file);        

my @fileLines; 

my $firstFile = 1;
my $session = "";

while (<$F>)
	{
		if ($firstFile == 1 && $_ =~ /TRACKING/)
			{ 				
				$firstFile = 0;
				push (@fileLines, $_);
			}
		elsif ($firstFile == 0 && $_ =~ /TRACKING/)
			{			
				my $file = $session.".txt";
    			my $F= new FileHandle;
    				
				vfopen ($F, ">$file");
						
				foreach my $element (@fileLines)
					{							
    					print $F "$element";
					}
				close ($F);	
					
				@fileLines = ();
				
			}
		else
			{
				if ($_ =~ /File\sName(.*)\\(.*)(\.trs)/)
					{
						$session = $2;
						print STDERR "-------------$session\n";
					}
				push (@fileLines, $_);
			}
		
	}

# Last File
my $LastFile = $session.".txt";
    			my $F= new FileHandle;
    				
				vfopen ($F, ">$LastFile");
						
				foreach my $element (@fileLines)
					{							
    					print $F "$element";
					}
				close ($F);	
	
sub vfopen 
  {
    my $f=shift;
    my $file=shift;

    if (($file =~/^\>/) && !($file =~/^\>\>/ )){open ($f, $file); return $f;}
    elsif (($file =~/^\>\>(.*)/))
      {
	if (!-e $1){	print STDERR "\nERROR: $file does not exist [FATAL]\n";exit(1);}
      }
    elsif (!-e $file){	print STDERR "\nERROR: $file does not exist [FATAL]\n";exit(1);}
   
    open ($f,$file);
    return $f;
  }	
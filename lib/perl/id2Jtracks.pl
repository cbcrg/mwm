#!/usr/bin/env perl

use strict;
use FileHandle;
use Data::Dumper;

my $file = shift (@ARGV);
#print "$file\n";

my $fileIdGen = shift (@ARGV);
my $FIdGen=new FileHandle;
open($FIdGen, $fileIdGen);

my $hGenID = {};
my $firstL = 0;

while (<$FIdGen>)
	{	
		chomp;	
		my $line = $_;
		
#		print "---$line";
		if ($firstL == 0) 
			{
				$firstL = 1; 
				next;
				
			}
			
		my @ary = split (',', $line);
#		print "$ary[0]    $ary[1]\n";
		$hGenID->{$ary[1].$ary[2]}->{$ary[0]} = 1;		
	}

#print Dumper ($hGenID);
#die;
my $F=new FileHandle;
open($F, $file);        

my $first = 0;

my (@aryWtSaline, @aryTSSaline, @aryWtEGCG, @aryTSEGCG);
my $mouseFind = 0;

while (<$F>)
  {
    my $l=$_;
           
    if ($l=~/^Subject Identification:\s(\w+)/)
      { 
      	my $id = $1;
#      	print "$id-----\n";
      	foreach my $mouse (sort {$a<=>$b} keys (%{$hGenID->{"WT1"}}))
      		{
      				if ($id =~ /$mouse/)
      					{
#      						print "$id--- WT";
      						push (@aryWtSaline, $id);
      						$mouseFind = 1;
      						next;
      					}      				
      		}
      		
      	foreach my $mouse (sort {$a<=>$b} keys (%{$hGenID->{"TS1"}}))
      		{
#      				print "mouse->$mouse\n";
      				if ($id =~ /$mouse/)
      					{
#     						print "$id--- TS";
      						push (@aryTSSaline, $id);
      						$mouseFind = 1;
      						next;
      					}
      		}
      		
      	foreach my $mouse (sort {$a<=>$b} keys (%{$hGenID->{"WT2"}}))
      		{
#      				print "mouse->$mouse\n";
      				if ($id =~ /$mouse/)
      					{
#     						print "$id--- WT2\n";
      						push (@aryWtEGCG, $id);
      						$mouseFind = 1;
      						next;
      					}
      		}
      		
      	foreach my $mouse (sort {$a<=>$b} keys (%{$hGenID->{"TS2"}}))
      		{
#      				print "mouse->$mouse\n";
      				if ($id =~ /$mouse/)
      					{
#     						print "$id--- TS";
      						push (@aryTSEGCG, $id);
      						$mouseFind = 1;
      						next;
      					}
      		}		     
        if ($mouseFind == 0) {print STDERR "Mouse not found $id\n";}      

      }
  }

#Printing results WT and TS mice
$first = 1;
foreach my $wtAnimal (@aryWtSaline)
	{
		if ($first == 1) 
			{
				print "WT_Saline\n";
				print "$wtAnimal"; 
				$first=0;
			}
        else 
        	{
        		print ",$wtAnimal"
        	};
	}
	     
print "\n";
$first = 1;
foreach my $TSAnimal (@aryTSSaline)
	{
		if ($first == 1) 
			{
				print "TS_Saline\n";
				print "$TSAnimal"; 
				$first=0;
			}
        else 
        	{
        		print ",$TSAnimal"
        	};
	}   

print "\n";
$first = 1;

foreach my $wtAnimalEGCG (@aryWtEGCG)
	{
		if ($first == 1) 
			{
				print "WT_EGCG\n";
				print "$wtAnimalEGCG"; 
				$first=0;
			}
        else 
        	{
        		print ",$wtAnimalEGCG"
        	};
	}
	     
print "\n";
$first = 1;
foreach my $TSAnimalEGCG (@aryTSEGCG)
	{
		if ($first == 1) 
			{
				print "TSE_EGCG\n";
				print "$TSAnimalEGCG"; 
				$first=0;
			}
        else 
        	{
        		print ",$TSAnimalEGCG";
        	}
	}   
print "\n";


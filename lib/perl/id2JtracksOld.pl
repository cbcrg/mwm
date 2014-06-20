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
				
		if ($firstL == 0) 
			{
				$firstL = 1; 
				next;				
			}
			
		my @ary = split (';', $line);
		$hGenID->{$ary[1].$ary[2]}->{$ary[0]} = 1;		
	}

#print Dumper ($hGenID);

#die;
my $F=new FileHandle;
open($F, $file);        

my $first = 0;

#    WT1                WT2            WT3            WT4             TS1             TS2               TS3             TS4
my (@aryWtSaline_NE, @aryWtSaline_EE, @aryWtEGCG_NE, @aryWtEGCG_EE, @aryTSSaline_NE, @aryTSSaline_EE, @aryTSEGCG_NE, @aryTSEGCG_EE);
my $mouseFind = 0;

while (<$F>)
  {
    my $l=$_;
           
    if ($l=~/^Subject Identification:\s(\w+)/)
#	if ($l=~/^Subject Identification\:/)
      { 
#		print $l;
      	my $id = $1;
      	foreach my $mouse (sort {$a<=>$b} keys (%{$hGenID->{"WT1"}}))
      		{
      				if ($id =~ /$mouse/)
      					{
#      						print "$id--- WT";
      						push (@aryWtSaline_NE, $id);
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
      						push (@aryWtSaline_NE, $id);
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
      						push (@aryWtSaline_EE, $id);
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
      						push (@aryTSSaline_EE, $id);
      						$mouseFind = 1;
      						next;
      					}
      		}
      	
      	foreach my $mouse (sort {$a<=>$b} keys (%{$hGenID->{"WT3"}}))
      		{
      				if ($id =~ /$mouse/)
      					{
#      						print "$id--- WT";
      						push (@aryWtEGCG_NE, $id);
      						$mouseFind = 1;
      						next;
      					}      				
      		}
      		
      	foreach my $mouse (sort {$a<=>$b} keys (%{$hGenID->{"TS3"}}))
      		{
#      				print "mouse->$mouse\n";
      				if ($id =~ /$mouse/)
      					{
#     						print "$id--- TS";
      						push (@aryTSEGCG_NE, $id);
      						$mouseFind = 1;
      						next;
      					}
      		}
      		
      	foreach my $mouse (sort {$a<=>$b} keys (%{$hGenID->{"WT4"}}))
      		{
#      				print "mouse->$mouse\n";
      				if ($id =~ /$mouse/)
      					{
#     						print "$id--- WT\n";
      						push (@aryWtEGCG_EE, $id);
      						$mouseFind = 1;
      						next;
      					}
      		}
      		
      	foreach my $mouse (sort {$a<=>$b} keys (%{$hGenID->{"TS4"}}))
      		{
#      				print "mouse->$mouse\n";
      				if ($id =~ /$mouse/)
      					{
#     						print "$id--- TS";
      						push (@aryTSEGCG_EE, $id);
      						$mouseFind = 1;
      						next;
      					}
      		}		     
      					     
        if ($mouseFind == 0) {print STDERR "Mouse not found $id\n";}      

      }
  }

#Printing results WT and TS mice
$first = 1;

# WT1
foreach my $wtAnimal (@aryWtSaline_NE)
	{
		if ($first == 1) 
			{
				print "WT_Water_NE\n";
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

# TS1
foreach my $TSAnimal (@aryTSSaline_NE)
	{
		if ($first == 1) 
			{
				print "TS_Water_NE\n";
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

# WT2
foreach my $wtAnimalEGCG (@aryWtSaline_EE)
	{
		if ($first == 1) 
			{
				print "WT_Water_EE\n";
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

# TS2
foreach my $TSAnimalEGCG (@aryTSSaline_EE)
	{
		if ($first == 1) 
			{
				print "TS_Water_EE\n";
				print "$TSAnimalEGCG"; 
				$first=0;
			}
        else 
        	{
        		print ",$TSAnimalEGCG";
        	}
	}   

print "\n";
$first = 1;

# WT3
foreach my $wtAnimalEGCG (@aryWtEGCG_NE)
	{
		if ($first == 1) 
			{
				print "WT_EGCG_NE\n";
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

# TS3
foreach my $TSAnimalEGCG (@aryTSEGCG_NE)
	{
		if ($first == 1) 
			{
				print "TSE_EGCG_NE\n";
				print "$TSAnimalEGCG"; 
				$first=0;
			}
        else 
        	{
        		print ",$TSAnimalEGCG";
        	}
	}   
print "\n";
$first = 1;

# WT4
foreach my $wtAnimalEGCG (@aryWtEGCG_EE)
	{
		if ($first == 1) 
			{
				print "WT_EGCG_EE\n";
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

# TS4
foreach my $TSAnimalEGCG (@aryTSEGCG_EE)
	{
		if ($first == 1) 
			{
				print "TSE_EGCG_EE\n";
				print "$TSAnimalEGCG"; 
				$first=0;
			}
        else 
        	{
        		print ",$TSAnimalEGCG";
        	}
	}   
print "\n";
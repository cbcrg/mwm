#!/usr/bin/env perl

use strict;
use FileHandle;
use XML::Simple;
use Data::Dumper;

if ($#ARGV ==-1)
  {    
    print "****************** Description **************************\n";
    print "Generates a csv file from a xml comming from MWM java application\n";
    print "****************** Command Line  ***********************\n";
    print "rhmm.pl -data <data_file>\n";
    print "****************** Flags      **************************\n";
    print "  -data  <file1 file2.. > ........File: input data from file(s).\n";
    print "  -action      <mode> ............Mode: 'convert' Converts xml files into csv.\n";
    print "  -output      <mode>.............Mode: csv or xml.\n";    
    print "****************** END **************************\n\n\n";
    die;	
  }
				
#my $cl=join(" ", @ARGV);
#my @commands=split (/\-+/,$cl);
#my @files = split (" ", shift @commands);

my $param;

$param = &process_param (@ARGV);

my ($data);

#Reads the data
$data = &readData ($data, $param);
print Dumper($data);

#print Dumper($param); #del
if ($data && $param->{action} eq "convert")
  {
		if (!$param->{output} || $param->{output} eq "csv")
			{
				&xml2csv ($data);
			}
		
		elsif ($param->{output} eq "xml")
			{
				;#do a function able to transform csv to xml
			}
		
		else
			{
				print "\n****ERROR: $param->{output} is an unknown pararmeter[FATAL]***\n";
	    	die;
			}
		
  }

#############################################
#                                           #
# FUNCTIONS                                 #
#                                           #
#############################################

#####################
# Parameters
#####################
	
sub process_param
  
  {
    my @arg=@_;
    my $cl=join(" ", @arg);
    
    my @commands=split (/\s\-+/,$cl);
    my $param={};
    
    
    foreach my $c (@commands)
      {
				if (!($c =~ /\S/)){next;}
				$c =~ /(\w+)\s*(.*)\s*/;
				my $k = $1;
				if (!$2) {$param->{$k} = 1;}
				else {$param->{$k} = $2;}
				$param->{$k} =~ s/\s*$//;
      }
    
    return check_parameters ($param);
  }

sub check_parameters 
  
  {
    my $p = shift;
    my $rp = {};
    
    $rp->{data} = 1;
    $rp->{action} = 1;
    $rp->{output} = 1;
        
    foreach my $k (keys (%$p))
      {
				if (!$rp->{$k})
	  			{
	    			print "\n****ERROR: $k is an unknown pararmeter[FATAL]***\n";
	    			die;
	  			}
	  			
				else
	  			{
	    			print "PARAM: -$k ---> [$p->{$k}]\n";
	  			}
      }
    return $p;
	}
	
sub readData 
	
	{
		my $d = shift;
		my $p = shift;
		
		if ($p->{data}) 
      {
				my @fl=split (/\s+/, $p->{data});
				
				return (\@fl);
#				foreach my $ff (@fl)
#	  			{
#	    			if ( -e $ff) 
#	    				{
#	    					$d=generic_undump_data($ff,$d, $p);
#	    				}
#	  			}
      	}	
	}
	
   
sub xml2csv
	{	
		my $ary_files = shift;
		my ($H, $header, $data, $file);
		
		foreach $file (@$ary_files)
			{
				my $xml = new XML::Simple;
				$H = $xml->XMLin($file);				
				$header = $H->{'header'};
				
	      &printHeader ($header, $file);							
			}
		
		foreach $file (@$ary_files)
			{
				my $xml = new XML::Simple;
				$H = $xml->XMLin($file);								
				$data = $H->{'data'}{'record'};
	            
	      &printData ($data, $file);				
			}
	}
	
sub printHeader
	{
		my $H = shift;
		my $f = shift;
		
		my ($k_1, $k_2, $k_3, $v);
				
		foreach $k_1 (sort ({$a cmp $b} keys(%$H)))
			{		
				foreach $k_2 (sort {$a cmp $b} keys (%{$H->{$k_1}}))
					{
												
						if ($k_2 eq "goalPosition")
							{
								foreach $k_3 (keys (%{$H->{$k_1}{$k_2}}))
									{
										$v = $H->{$k_1}{$k_2}{$k_3};
										print "#h;file;$f;$k_1;$k_2;$k_3;$v;\n";
									}
							}
						
						else 
							{
								$v = $H->{$k_1}{$k_2};
								print "#h;file;$f;$k_1;$k_2;$v;\n";		
							}																		 
					}
			}
	}
		
sub printData
	{
		my $H = shift;
		my $f = shift;
		
		my ($k_1, $k_2, $k_3, $v);
		
		foreach $k_1 (sort {$a <=> $b} keys(%$H))
			{								
				print "#d;index;$k_1;";
															
				foreach $k_2 (sort {$b cmp $a} keys(%{$H->{$k_1}}))
					{						
						if ($k_2 eq "time")
							{
								foreach $k_3 (sort {$a cmp $b} keys(%{$H->{$k_1}{$k_2}}))
									{									
										$v = $H->{$k_1}{$k_2}{$k_3};
										$k_3 .="T";
										print "$k_3;$v;"; 
									}
							}
						else
							{							
								foreach $k_3 (sort {$a cmp $b} keys(%{$H->{$k_1}{$k_2}}))
									{									
										$v = $H->{$k_1}{$k_2}{$k_3};
																				
										###Harcode for fast binning!!!
										if ($k_3 eq "orientation")
											{
												$k_3 = "bin";	
											}
										
										if ($v eq "-1")
											{
												$v = 2;	
											}
										###END -- Harcode for fast binning!!!
										
										print "$k_3;$v;"; 
									}												
							}						
					}
					print "file;$f;\n";
				
			}
	}
		
		
		

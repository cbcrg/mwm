#!/usr/bin/env perl

use strict;
use FileHandle;
use XML::Simple;
use Data::Dumper;

my $cl=join(" ", @ARGV);
my @commands=split (/\-+/,$cl);
my @files = split (" ", shift @commands);

&xml2csv (\@files);

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
		
		
		

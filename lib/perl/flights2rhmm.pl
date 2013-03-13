#!/usr/bin/env perl

use strict;
use FileHandle;
use XML::Simple;
use Data::Dumper;

if ($#ARGV ==-1)
  {    
    print "\n****************** Description **************************\n";
    print "Read flights and generates a fiel readable by rhmm.pl\n";
#    print "Bins the data and show statistics\n";
    print "****************** Command Line  ***********************\n";
    print "flights2rhmm.pl -data <data_file>\n";
    print "****************** Flags      **************************\n";
    print "  -data  <file1 file2.. > ........File: input data from file(s).\n"; 
#    print "  -infile_format <mode>.............Mode: csv or xml.\n";   
#    print "  -output        <mode>.............Mode: csv or xml.\n";
#    print "  -add           <mode>.............Mode: direction.\n";
#    print "  -bin           <mode>.............Mode: angle.\n";
#    print "  -action        <mode> ............Mode: 'convert' Converts xml files into csv.\n";
#    print "  ........................................'countBins' Perform bins counts.\n";
#    print "  ........................................'logodd' Calculates bins logodd ratio.\n";
#    print "  -outBins       <mode> ............Mode:  name of the output files containning the bin and the dibin counts.\n";
#    print "  .........................................'no' bins counts are not shown.\n";
#    print "  -outLogodd     <mode> ............Mode:  'screen' logodd ratio results on screen\n";
#    print "  ..........................................name of the output files containning the logodd ratio results.\n";
#    print "  .........................................'no' logodd ratios not shown.\n";
#    print "  -outdata       <mode> ............Mode:  'no' the data is not shown\n";
    #print "  -stats         <mode> ............Mode: 'logodd' Calculates logodd of bins.\n"; #REVIEW IF IS WORTH TO IMPLEMENT A DIFFERENT OPTION (STATS)
    #AND NO INSIDE -action    
    print "****************** END **************************\n\n\n";
    die;  
  }

my $param;
my $BIN=0;

$param = &process_param (@ARGV);

my ($data, $bin, $logodd);

#Reads the data
$data = &readData ($data, $param);

if ($data && $param->{bin})
  {     
    &data2bin ($data, $param); 
  }

if ($data && $param->{outdata} ne "no")
  {  
    &genericWriteData ($data, $param); 
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
    $rp->{outdata} = 1;
    $rp->{readMode} = 1;
    $rp->{bin} = 1;
    $rp->{binField} = 1;
    $rp->{binN} = 1;
    $rp->{binDelta} = 1;
    $rp->{binName} = 1;
    $rp->{binBoundaries} = 1;
    $rp->{format} = 1;
    $rp->{out} = 1;
#    $rp->{output} = 1;

#    $rp->{add} = 1;

#    $rp->{countBin} = 1;
#    $rp->{infile_format} = 1;
#    $rp->{outBins} = 1;
#    $rp->{outLogodd} = 1;
    
    
    foreach my $k (keys (%$p))
      {
        if (!$rp->{$k})
          {
            print STDERR "\n****ERROR: $k is an unknown pararmeter[FATAL]***\n";
            die;
          }
          
        else
          {
            print STDERR "PARAM: -$k ---> [$p->{$k}]\n";
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
                
        foreach my $ff (@fl)
          {                       
            if ( -e $ff) 
              {
                my $F= new FileHandle;
                vfopen ($F,$ff);
                my $l=<$F>;                                
                close ($F);
                                    
                #$d = &generic_undump_data ($ff, $d, $p);
                #Format: flights2rhmm.01
                if ($l =~ /Format: flights2rhmm.01/ || $l=~/^#d/ || $l =~ /Format: rhmm.data.01/)
                  {
                    $d = &csvFile2hash ($ff, $d, $p);
                  }
                elsif ($l!~/^#/)
                  {
                    $d = &vectorFlightsFile2hash ($ff, $d, $p);                    
                  }                               
              }
             
            else 
            	{
            		print STDERR "\nERROR: $ff does not exist [FATAL]\n";
            		exit(1);
            	}
          }
        }
         
      return ($d);
  }
  
sub vectorFlightsFile2hash
  {
    my $file = shift;
    my $d = shift;    
    my $param = shift;
    
    my $F= new FileHandle;
    my $pFlight;
    my $ctr = 1;
    
    vfopen ($F, $file);
    
    while (<$F>) 
      {
        my $l=$_;
        chomp($l);
        
        my @aryL = split (" ",$l);
        
        my ($index, $flight) = 1;
        
        if (!$param->{readMode} || $param->{readMode} eq "join")
          {
            $flight = $aryL[1];
            
            if ($flight == $pFlight)
              {               
                $pFlight = $flight;
                next;
              }
            else
              {                
                $index = $ctr;                
                $pFlight = $flight;
                $ctr++;  
              }  
            
          }
        elsif ($param->{readMode} eq "raw")
          {            
            $index = $aryL[0];
            $flight = $aryL[1];  
          }
        else
          {
            &errorMng ("\nERROR: readMode: $param->{readMode} unknown [FATAL]\n");  
          }            
        
        $d->{$file}{$index}{'1'} = $index;
        $d->{$file}{$index}{'flight'} = $flight;
        $d->{$file}{$index}{'file'} = $file;
           
      }
      
    return ($d);  
  }

sub csvFile2hash
  {
    my $file = shift;
    my $d = shift;    
    my $param = shift;
    
    my $F= new FileHandle; 
    vfopen ($F, $file);   
    
    while (<$F>)
	 {
      my $l=$_;
      chomp($l);
      my $L={};
      #print "$l\n";
      if ( $l=~/#d/)
        {
  	     my @v = ($l=~/([^;]+)/g);
  	     shift @v;	  	     

		 while (@v)
	      {
		    my $key=shift @v;
		    my $value= shift @v;
		    $L->{$key}=$value;		
	      }
        }
	  
	  my $f = $L->{file};
	  my $index = $L->{1};
		 
	  foreach my $k (keys(%$L))
	   {
	     $d->{$f}{$index}{$k} = $L->{$k};
	   }  		 
    }
    close ($F);
	#print Dumper ($d);
	return $d;
  }
   
sub vfopen 
  {
    my $f=shift;
    my $file=shift;

    if (($file =~/^\>/) && !($file =~/^\>\>/ )){open ($f, $file); return $f;}
    elsif (($file =~/^\>\>(.*)/))
      {
  			if (!-e $1){  print STDERR "\nERROR: $file does not exist [FATAL]\n";exit(1);}
      }
    elsif (!-e $file){  print STDERR "\nERROR: $file does not exist [FATAL]\n";exit(1);}
   
    open ($f,$file);
    return $f;
  }

sub data2bin
  {
    my $d = shift;
    my $param = shift;    
    
    my $field = $param->{binField};
    my $nbin = $param->{binN};
    my $delta = $param->{binDelta};
    my $name = $param->{binName};
        
    if (!$field) {$field="flight";}
    if (!$nbin) {$nbin=1;}
    if (!$delta) {$delta=0.02;}
    if (!$name) 
      {
    	$BIN++;
    	$name = "BIN$BIN"."_";
      }
    elsif ($name eq "numeric")
      {
        $name = "";
      }
    else
      {
        $name .= "_";
      }  
      
    foreach my $f (keys(%$d))
      {
    	foreach my $i (keys (%{$d->{$f}}))
    	  {    	            
    	     #if (($S->{$c}{$i}{$field}/$delta < 0) && ($S->{$c}{$i}{$field}/$delta))
    		 if ($d->{$f}{$i}{$field}/$delta < 0)
              {
    		    my $bin = 0;
    		        
    		    if ($d->{$f}{$i}{bin}) 
    		      {
    		        $d->{$f}{$i}{bin}="$d->{$f}{$i}{bin}"."::"."$name"."$bin";
    		      }
    		    else
    		      {
    		        $d->{$f}{$i}{bin}="$name"."$bin";
    		      }
        	  }
        		    
             else
        	  {
        	    my $bin = int ($d->{$f}{$i}{$field}/$delta);
        		        
        		if ( $bin<0) {$bin=0;}
        		else
  		          {
  			       $bin = ($bin >= $nbin)? $nbin:$bin+1;
  		          }
  		        
  		        if ($d->{$f}{$i}{bin}) 
  		          {
  		            $d->{$f}{$i}{bin}="$d->{$f}{$i}{bin}"."::"."$name"."$bin";
  		          }
  		        else
  		          {
  		            $d->{$f}{$i}{bin}="$name"."$bin";
  		          }
        		}         
              }
    	   }
            
    delete ($param->{binField});
    delete ($param->{binN});
    delete ($param->{binDelta});
    delete ($param->{binName});

    
    return ($d);
  }


sub genericWriteData
  {
    my $d = shift;
    my $p = shift;
    
    my $file = (exists ($p->{out}))? $p->{out} : "";
    my $csvExt = ".csv";
    my $RExt = ".R";
    
    if (exists ($p->{format}))
      {       
        if ($p->{format} eq "csv")
          {            
            $file .= $csvExt;
            &writeData2csv ($d, $p, $file);
          }
        elsif ($p->{format} eq "R")
          {
            if ($file) {$file .= $RExt;}
            &writeData2R ($d, $p, $file);
          }
        else
          {
            &errorMng ("\nERROR: Output format: $param->{format} unknown [FATAL]\n");     
          }
      }
    else
      {
         if ($file) {$file .= $csvExt;}
         &writeData2csv ($d, $p, $file); 
      }      
  }  
  
sub writeData2csv
  {
    my $d = shift;
    my $p=shift;
    my $file = shift;
    
    #print Dumper ($H);
    my ($i, $f, $ext, $last) = "";                 
    my ($k_1, $k_2, $v);
    
    
    my $F= new FileHandle;
        
    if (!$file){open ($F, ">-");}
    else {open ($F, ">$file");}
    
    print $F "#comment;Format: flights2rhmm.01\n";    
    
    foreach $f (sort {$a cmp $b} keys(%$d))
      { 
#        $file= $f;
#               
#        if ($file =~ /(\.[^.]+)$/)
#          {
#            $ext = $1;
#            $file =~ s/(\.[^.]+)$//;
#            $file .= ".csv";  
#          }
#        else
#          {
#            $file .= ".csv";
#          }  
                    
        if ($p->{binBoundaries}) {print $F "#d;1;0;bin;BEGIN;file;$f;flight;0;\n";}                      
        foreach $i (sort {$a <=> $b} keys(%{$d->{$f}}))
          {            
            
            print $F "#d;";
            foreach my $k (sort {$a cmp $b} keys(%{$d->{$f}{$i}}))
              {                                                          
                print $F "$k;$d->{$f}{$i}{$k};";
                $last = $i+1;                               
              }
            print $F "\n";  
              
          } 
        if ($p->{binBoundaries}) {print $F "#d;1;", "$last", ";bin;END;file;$f;flight;0;\n";}                                        
        
      }
    close ($F);
  } 
   
 sub writeData2R
  {
    my $d = shift;
    my $p=shift;
    my $file = shift;
    
    my $F= new FileHandle;
        
    if (!$file){open ($F, ">-");}
    else {open ($F, ">$file");}
    
    &data2tableNames ($d, $F);
    &data2tableRec ($d, $F);
    
    close ($F);    
    
  } 

sub data2tableNames
  {
    my $d = shift;
    my $F = shift;
    
	foreach my $f (sort ({$a<=>$b}keys(%$d)))
	  { 
	    foreach my $i (sort {$a<=>$b}keys (%{$d->{$f}}))
	      {
		    my $first=0;
		
		    foreach my $k (sort (keys (%{$d->{$f}{$i}})))
		      {		        
		        if ($first == 0) 
		          {
		            $k = ($k == 1)? "seqOrder" : $k; 
		            print $F "$k"; $first=1;
		          }
		        else {print $F "\t$k";}		     
		      }		
		
	       print $F "\n";
	       last;
	  }	  
	  last;	    
    }	
  }

sub data2tableRec
  {
    my $d = shift;
	my $F = shift;
	  
    foreach my $f (sort ({$a<=>$b}keys(%$d)))
      { 
        foreach my $i (sort {$a<=>$b}keys (%{$d->{$f}}))
          {
	       my $first=0;
	       foreach my $k (sort (keys (%{$d->{$f}{$i}})))
	         {
	           if ($first == 0) 
	             {
		          print $F "$d->{$f}{$i}{$k}"; 
		          $first=1;
	             }
	           else 
	             {
		          print $F "\t$d->{$f}{$i}{$k}";
	             }		     
	         }		
	      
	    print $F "\n";
      }	    
    }
  }
	  
sub errorMng
  {
    my $msg = shift;
    
    $msg = "\nERROR: ".$msg." [FATAL]\n";
    print "$msg";
  }          
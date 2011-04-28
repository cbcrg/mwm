#!/usr/bin/env perl

use strict;
use FileHandle;
use XML::Simple;
use Data::Dumper;

my $HEADER;

if ($#ARGV ==-1)
  {    
    print "\n****************** Description **************************\n";
    print "Generates a csv file from an xml comming from MWM java application\n";
    print "Bins the data and show statistics\n";
    print "****************** Command Line  ***********************\n";
    print "xml2csv.pl -data <data_file>\n";
    print "****************** Flags      **************************\n";
    print "  -data  <file1 file2.. > ........File: input data from file(s).\n"; 
    print "  -infile_format <mode>.............Mode: csv or xml.\n";   
    print "  -output        <mode>.............Mode: csv or xml.\n";
    print "  -bin           <mode>.............Mode: angle.\n";
    print "  -action        <mode> ............Mode: 'convert' Converts xml files into csv.\n";
    print "  ........................................'countBins' Show bin statistics.\n";
    print "  -outbins       <mode> ............Mode:  name of the output files containning the bin and the dibin counts.\n";    
    print "****************** END **************************\n\n\n";
    die;  
  }
        
#my $cl=join(" ", @ARGV);
#my @commands=split (/\-+/,$cl);
#my @files = split (" ", shift @commands);

my $param;

$param = &process_param (@ARGV);

my ($data, $bin);

#Reads the data
$data = &readData ($data, $param);

#Bins data
if ($param->{bin})
  {        
    $data = &data2bin ($data, $param->{bin});        
  }

#if ($data && $param->{action} eq "convert")
#if (!$param->{output} || $param->{output} eq "data")

if ($param->{action} eq "convert")
  {      
    if ($param->{infile_format} eq "xml")
      {            
        &data2printCsv ($data);
      }
    
    elsif ($param->{infile_format} eq "csv")
      {
        ;#do a function able to transform csv to xml
      }
    
    else
      {
        print STDERR "\n****ERROR: file format is unknown [FATAL]***\n";
        die;
      }
    die;
  }

elsif ($param->{action} eq "countBins")
  
  {
    $bin = &checkDataBin ($data, $bin);            
  }
  
#if ($Data && $param->{outdata} ne "no")
##outdata option in rhmm is used to provided a name for the output files

#Dumping results
$param = set_output_name ($param);

if ($data && $param->{outdata} ne "no")
  {  
    if (!$param->{output} || $param->{output} eq "csv")
      {            
        &data2printCsv ($data);
      }
    
    elsif ($param->{output} eq "xml")
      {
        ;#do a function able to transform csv to xml #REVIEW
      }
    
    else
      {
        print STDERR "\n****ERROR: $param->{output} is an unknown pararmeter[FATAL]***\n";
        die;
      }
  }

if ($bin && $param->{outBins} ne "no")
  {  
  	printBins ($bin->{'bin'}, $bin->{'bin_T'}, "$param->{outBins}.bin");
		printBins (&hash2hashCsvLike ($bin->{'dibin'},"1"), &hash2hashCsvLike ($bin->{'dibin_T'}, "0"), "$param->{outBins}.dibin");
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
    $rp->{outdata} = 1;
    $rp->{bin} = 1;
    $rp->{countBin} = 1;
    $rp->{infile_format} = 1;
    $rp->{outBins} = 1;
    $rp->{out} = 1;
    
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

sub set_output_name
	{
		my $param=shift;
		
		if (!$param->{out})
	  	{
	  		$param->{out} = "out";	    
	  	} 
	  
	  if (!$param->{outBins})
	  	{
	  		$param->{outBins} = "$param->{out}";
	  	}	
	  
	  return ($param); 	
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
                $d = &generic_undump_data ($ff, $d, $p);                
              }
          }
        }  
      return ($d);
  }
  
#sub xml2csv
#  {  
#    my $ary_files = shift;
#    my ($H, $header, $data, $file);
#    
#    foreach $file (@$ary_files)
#      {
#        my $xml = new XML::Simple;
#        $H = $xml->XMLin($file);        
#        $header = $H->{'header'};
#        
#        #&printHeader ($header, $file);              
#      }
#    
#    foreach $file (@$ary_files)
#      {
#        my $xml = new XML::Simple;
#        $H = $xml->XMLin($file);                
#        $data = $H->{'data'}{'record'};
#        print Dumper ($data);      
#        #&printData ($data, $file);        
#      }
#  }

sub generic_undump_data
  {
    my $file=shift;
    my $d=shift;
    my $p=shift;
    my $field=shift;
    my $F= new FileHandle;

    vfopen ($F,$file);
    my $l=<$F>;
    close ($F);
    
    if ($l=~/\?xml/)
      { 
        if (exists ($p->{infile_format}) && ($p->{infile_format} ne "xml"))
          {
            print STDERR "***** ERROR: Format provided $p->{infile_format} is not the actual format of file $file [FATAL]\n";
            die;
          }
          
        $p->{infile_format} = "xml";       
        return &xml2hash ($file,$d,$p,$field);
      }
    
    elsif ($l=~/Format: csv/ || $l=~/Format: rhmm.data.01/)
      {
        if (exists ($p->{infile_format}) && ($p->{infile_format} ne "csv"))
          {
            print STDERR "***** ERROR: Format provided $p->{infile_format} is not the actual format of file $file [FATAL]\n";
            die;
          }
          
        $p->{infile_format} = "csv";
        
        return &csv2hash ($file,$d,$p,$field);
      }
      
    else 
      {        
        print STDERR "***** ERROR: Format of $file is unknown [FATAL]\n";
        die;
      }
  }
  
sub xml2hash
  {  
    my $file = shift;
    my $data = shift;    
    my $param = shift;
    
    my ($H, $header, $xml);
    
    $xml = new XML::Simple;
    $H = $xml->XMLin($file);            
    $data->{$file}{'header'} = &hash2hashCsvLike ($H->{'header'}, 0);
    $data->{$file}{'data'} = &hash2hashCsvLike ($H->{'data'}{'record'}, 1);    
    #$data->{$file} = &hash2hashCsvLike ($H);   
    return ($data)
  }

sub csv2hash
  {
    my $file = shift;
    my $data = shift;    
    my $param = shift;
    my $F= new FileHandle;
    my ($d, $file2save);
    
    vfopen ($F, $file);
    
    while (<$F>)
      {
        my $l=$_;
                
        if ($l=~/#h/ || /#comment/)
          {            
            $HEADER.=$l;
          }
          
        elsif ($l=~/#d/)
          {
            chomp($l);
                        
            if ($l =~ m/\;file\;.([A-Za-z0-9_\/\.]+)/) 
              {
                $file2save = $1;                 
              }
            
            my @ary_l = ($l=~/([^;]+)/g);
            shift @ary_l;  #get rid of #d
            my $exp = shift (@ary_l);
            my $index = shift(@ary_l);
                                        
            for (my $a=0; $a<=$#ary_l; $a+=2)
              {
                #$d->{$file}{'data'}{'record'}{$index}{$ary_l[$a]}=$ary_l[$a+1];
                $d->{$file2save}{'data'}{$index}{$ary_l[$a]}=$ary_l[$a+1];
              }            
          }
      } 
    return ($d)
  }

#hasIndex: If the hash has a record number we do not use is in the value we are going to show
#          eg  {1}{animalPosition}{angle}=90, key =animalPosition#angle and not 1#animalPosition#angle
 
sub hash2hashCsvLike
  {
    my $H = shift;
    my $hasIndex = shift;
    
    my ($h, $d, $r, $main_k, $k_1, $k_2, $k_3, $d, $key, $v);
        
    if ($hasIndex)
      {         
        foreach $k_1 (sort {$a <=> $b} keys(%$H))
          {
            foreach $k_2 (sort {$b cmp $a} keys(%{$H->{$k_1}}))
              {  
                foreach $k_3 (sort {$a cmp $b} keys(%{$H->{$k_1}{$k_2}}))
                  {
                    $key= $k_2."#".$k_3;
                    $v = $H->{$k_1}{$k_2}{$k_3};                    
                    $r->{$k_1}{$key} = $v; 
                  }                                  
              }                                                
          }  
        }
      
      else
        {
          foreach $k_1 (sort {$a <=> $b} keys(%$H))
            {
              foreach $k_2 (sort {$b cmp $a} keys(%{$H->{$k_1}}))
                {  
                  if (ref ($H->{$k_1}{$k_2}) eq "HASH")
                    {                                                    
                      foreach $k_3 (sort {$a cmp $b} keys(%{$H->{$k_1}{$k_2}}))
                        {
                          $key= $k_1."#".$k_2."#".$k_3;
                          $v = $H->{$k_1}{$k_2}{$k_3};                          
                          $r->{$key} = $v; 
                        }
                                    
                    }
                  else
                    {                      
                      $key= $k_1."#".$k_2;
                      $v = $H->{$k_1}{$k_2};
                      $r->{$key} = $v; 
                                          
                                    
                    }
                }
            }
        }            
      return ($r);
  }

sub data2printCsv
  {
    my $d = shift;
    my ($f);
    
    if (defined ($HEADER))
      {
        print "$HEADER";
      }
      
    else
      {
        print "#comment;Format: csv\n";
        
        foreach $f (sort ({$a cmp $b} keys(%$d)))
          {        
            printHeader ($data-> {$f}{'header'}, $f);            
          }
      }    
            
    foreach $f (sort ({$a cmp $b} keys(%$d)))
      {
        printData ($data-> {$f}{'data'}, $f);
      } 
  }  

#sub printHeader
#  {
#    my $H = shift;
#    my $f = shift;
#    
#    my ($k_1, $k_2, $k_3, $v);
#    
#    foreach $k_1 (sort ({$a cmp $b} keys(%$H)))
#      {    
#        foreach $k_2 (sort {$a cmp $b} keys (%{$H->{$k_1}}))
#          {
#                        
#            if ($k_2 eq "goalPosition")
#              {
#                foreach $k_3 (keys (%{$H->{$k_1}{$k_2}}))
#                  {
#                    $v = $H->{$k_1}{$k_2}{$k_3};
#                    print "#h;file;$f;$k_1;$k_2;$k_3;$v;\n";
#                  }
#              }
#            
#            else 
#              {
#                $v = $H->{$k_1}{$k_2};
#                print "#h;file;$f;$k_1;$k_2;$v;\n";    
#              }                                     
#          }
#      }
#  }
  
sub printHeader
  {
    my $H = shift;
    my $f = shift;
    
    my ($k_1, $v);
        
    foreach $k_1 (sort {$a <=> $b} keys(%$H))
      {
        $v = $H->{$k_1};        
        print "#h;file;$f;$k_1;$v\n";
      }
  }  

sub printData
  {
    my $H = shift;
    my $f = shift;
    
    my ($k_1, $k_2, $v);
        
    foreach $k_1 (sort {$a <=> $b} keys(%$H))
      {
        print "#d;index;$k_1;";
        foreach $k_2 (sort {$a <=> $b} keys(%{$H->{$k_1}}))
          {
            if ($k_2 eq "file")
              {
                $f = $H->{$k_1}{$k_2};
                next;
              } 
            
            $v = $H->{$k_1}{$k_2};
        
            print "$k_2;$v;";
          }
        
        $f = path2fileName ($f);
        print "file;$f\n";
      }
  }    
  
  
#sub printData
#  {
#    my $H = shift;
#    my $f = shift;
#    
#    my ($k_1, $k_2, $k_3, $v);
#    
#    foreach $k_1 (sort {$a <=> $b} keys(%$H))
#      {                
#        print "#d;index;$k_1;";
#                              
#        foreach $k_2 (sort {$b cmp $a} keys(%{$H->{$k_1}}))
#          {            
#            if ($k_2 eq "time")
#              {
#                foreach $k_3 (sort {$a cmp $b} keys(%{$H->{$k_1}{$k_2}}))
#                  {                  
#                    $v = $H->{$k_1}{$k_2}{$k_3};
#                    $k_3 .="T";
#                    print "$k_3;$v;"; 
#                  }
#              }
#            else
#              {              
#                foreach $k_3 (sort {$a cmp $b} keys(%{$H->{$k_1}{$k_2}}))
#                  {                  
#                    $v = $H->{$k_1}{$k_2}{$k_3};
#                                                        
#                    ###Harcode for fast binning!!!
##                    if ($k_3 eq "orientation")
##                      {
##                        $k_3 = "bin";  
##                      }
##                    
##                    if ($v eq "-1")
##                      {
##                        $v = 2;  
##                      }
#                    ###END -- Harcode for fast binning!!!
#                    
#                    print "$k_3;$v;"; 
#                  }                        
#              }            
#          }
#          print "file;$f;\n";
#        
#      }
#  }

sub data2bin
  {
    my $data = shift;
    my $field = shift; 
    my ($f, $H, $k_1, $k_2, $k_3, $v, $bin, $flag, $pbin); 
    
    foreach $f (sort ({$a cmp $b} keys(%$data)))
      {
        
        $H = $data-> {$f}{'data'};
        
        foreach $k_1 (sort ({$a cmp $b} keys(%$H)))
          {    
            foreach $k_2 (sort {$a cmp $b} keys (%{$H->{$k_1}}))
              {                
                next if ($k_2 !~ /($field)/);                
                
                $flag = 1;
                                
                $v = $H->{$k_1}{$k_2};
                
                if ($field =~ "angle")
                  {
                    SWITCH: 
                      {
                        (($v <= 45 && $v >0) || $v > 320 && $v <=360 ) && do 
                          { 
                            $bin = 1;                            
                            last SWITCH;
                          };
                                                    
                        ($v <= 320 && $v > 225 ) && do 
                          { 
                            $bin = 2;                             
                            last SWITCH;
                          };
                            
                        ($v <= 225 && $v > 135 ) && do 
                          { 
                            $bin = 3;                            
                            last SWITCH;
                          };
                            
                        ($v <= 135 && $v > 45) && do 
                          { 
                            $bin = 4;                            
                            last SWITCH;
                          };
                        
                        ($v == 999.000) && do #Mouse in same position within two consecutives moments
                          {                             
                            $bin = $pbin;                            
                            last SWITCH;
                          };                         
                      }
                    
                    $H->{$k_1}{'postProcessed#bin'} = $bin;
                    $pbin = $bin;
                  }
                
                else
                  {
                    print STDERR "\n****ERROR: Is not possible to bin by $field [FATAL]***\n";
                    die;
                  }
              }
          }
      }
        
      if (!defined ($flag))
        {          
          print STDERR "\n****ERROR: Is not possible to bin by $field, field does not exist [FATAL]***\n"; 
          die; 
        }
      
      return ($data);
  }

##Old version of data2bin 
#sub data2bin
#  {
#    my $data = shift;
#    my $field = shift; 
#    my ($f, $H, $k_1, $k_2, $k_3, $v, $bin); 
#    
#    foreach $f (sort ({$a cmp $b} keys(%$data)))
#      {
#        
#        $H = $data-> {$f}{'data'};
#        print Dumper ($H);#del 
#        
#        foreach $k_1 (sort ({$a cmp $b} keys(%$H)))
#          {    
#            foreach $k_2 (sort {$a cmp $b} keys (%{$H->{$k_1}}))
#              {
#                next if ($k_2 eq "time");
#                
#                foreach $k_3 (keys (%{$H->{$k_1}{$k_2}}))
#                  {
#                    next if ($k_3 ne $field);
#                    #print  "$k_3\n";#del
#                    $v = $H->{$k_1}{$k_2}{$k_3};
#                    
#                    SWITCH: 
#                      {
#                        (($v <= 45 && $v >0) || $v > 320 && $v <=360 ) && do 
#                          { 
#                            $bin = 1; 
#                            last SWITCH;
#                          };
#                                                  
#                        ($v <= 320 && $v > 225 ) && do 
#                          { 
#                            $bin = 2; 
#                            last SWITCH;
#                          };
#                          
#                        ($v <= 225 && $v > 135 ) && do 
#                          { 
#                            $bin = 4; 
#                            last SWITCH;
#                          };
#                          
#                        ($v <= 135 && $v > 45) && do 
#                          { 
#                            $bin = 3; 
#                            last SWITCH;
#                          };
#                                                              
#                      }
#                      #print "$v ---> $bin\n"; 
#                      $H->{$k_1}{$k_2}{'bin'} = $bin;
#                  }                
#              }
#          
#          }
#          $data-> {$f}{'data'} = $H;
#      } 
#      #print Dumper ($data);#del
#      return ($data);
#      
#  }
  

    
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
    
sub checkDataBin
  {
    
    my $d = shift;
    my $bn = shift;
    
    my ($k_1, $k_2, $k_3);
    
    foreach $k_1 (keys(%$d))
      {
        foreach $k_2 (keys(%{$d->{$k_1}}))        
          {  
            next if ($k_2 eq "header");
            
            foreach $k_3 (sort {$a cmp $b} keys(%{$d->{$k_1}{$k_2}}))
              {
                if (exists ($d->{$k_1}{$k_2}{$k_3}{"postProcessed#bin"}))
                  {
                    $bn =  (&countBins ($d));
                    return ($bn);
                  }
                  
                else
                  {
                    print STDERR "\n****ERROR: Data is not binned [FATAL]***\n";
                    die;
                  }
              }
          }
      }
  }
  
sub countBins
  {
    my $d = shift;    
    
    my ($k_1, $k_2, $k_3, $v, $bin, $f, $bin_T, $dibin, $dibin_T, $pv, $h);
    
    foreach $k_1 (keys(%$d))
      {
        foreach $k_2 (keys(%{$d->{$k_1}}))        
          {  
            next if ($k_2 eq "header");
            
            foreach $k_3 (sort {$a cmp $b} keys(%{$d->{$k_1}{$k_2}}))
              {                                  
                  $v = $d->{$k_1}{$k_2}{$k_3}{'postProcessed#bin'};
                  #$v .= "_bin";###REVIEW if it should be key value or only value and add to all of them the same option               
                  
                  $f = $k_1;  
                  
                  $f = path2fileName ($f);

                  $bin->{$f}{$v}++;#del ###REVIEW if it should be key value or only value and add to all of them the same option
                  #$bin->{$f}{'bin'}{$v}++;
                  
                  $bin_T->{$v}++;                  
                  
                  if (!defined ($pv))
                    {
                      $pv = $v;
                      next;
                    } 
                    
                  $dibin->{$f}{$pv}{$v}++;
                  $dibin_T->{$pv}{$v}++;                  
                  $pv = $v;                            
              }
          }        
        }
      
#      print Dumper ($bin);#del
#      print Dumper ($bin_T);#del
#      print Dumper ($dibin);#del
#      print Dumper ($dibin_T);#del
      
      $h->{'bin'} = $bin;
      $h->{'bin_T'} = $bin_T;
      $h->{'dibin'} = $dibin;    
      $h->{'dibin_T'} = $dibin_T;          
      
      $h->{'dibin'} = $dibin;
      
      return ($h);
  }

sub printBins
  {
    my $h = shift;
    my $binList = shift;    
    my $outf = shift;
    
    my ($f, $bin, $v); #, @bin_angle);
    
    my $F= new FileHandle;
	
		if ($outf)
			{
				&vfopen ($F, ">$outf");
			}
		
		else 
			{
				$F=*STDOUT;
			}
	
		#@bin_angle = ("1", "2", "3", "4");
        
    #print "file\tbin1\tbin2\tbin3\tbin4\n";  
    print $F "file\t";
    
    foreach $bin ((sort {$a cmp $b} keys(%$binList)))
    	{
      	print $F "\t$bin";
      } 
    
    print $F "\n";
    
    foreach $f ((sort {$a cmp $b} keys(%$h)))
      {          
        print $F "$f";
                        
        #foreach $bin (@bin_angle)
        foreach $bin ((sort {$a cmp $b} keys(%$binList)))
          {
            print $F "\t";
            $v = $h->{$f}{$bin};
            if (!defined ($v)) 
              {
                $v = 0;
              }
              
            print $F "$v";    
          }
        print $F "\n";  
      }
    
  }

sub path2fileName 
  {
    my $f = shift;
                  
    if ($f =~ /^*.\//) #avoiding names as us/cn/file.act -> file.act  #REVIEW!!!!!
      {      
        my @a = split ("/",$f);        
        $f = pop (@a);      
      } 
    
    return ($f);
    
  }
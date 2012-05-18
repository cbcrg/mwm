#!/usr/bin/env perl

use strict;
use FileHandle;
use XML::Simple;
use Data::Dumper;

my $pi = 3.14159265358979;;

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
    print "  -add           <mode>.............Mode: direction.\n";
    print "  -bin           <mode>.............Mode: angle.\n";
    print "  -action        <mode> ............Mode: 'convert' Converts xml files into csv.\n";
    print "  ........................................'countBins' Perform bins counts.\n";
    print "  ........................................'logodd' Calculates bins logodd ratio.\n";
    print "  -outBins       <mode> ............Mode:  name of the output files containning the bin and the dibin counts.\n";
    print "  .........................................'no' bins counts are not shown.\n";
    print "  -outLogodd     <mode> ............Mode:  'screen' logodd ratio results on screen\n";
    print "  ..........................................name of the output files containning the logodd ratio results.\n";
    print "  .........................................'no' logodd ratios not shown.\n";
    print "  -outdata       <mode> ............Mode:  'no' the data is not shown\n";
    #print "  -stats         <mode> ............Mode: 'logodd' Calculates logodd of bins.\n"; #REVIEW IF IS WORTH TO IMPLEMENT A DIFFERENT OPTION (STATS)
    #AND NO INSIDE -action    
    print "****************** END **************************\n\n\n";
    die;  
  }
        
#my $cl=join(" ", @ARGV);
#my @commands=split (/\-+/,$cl);
#my @files = split (" ", shift @commands);

my $param;

$param = &process_param (@ARGV);

my ($data, $bin, $logodd);

#Reads the data
$data = &readData ($data, $param);

#Adds new fields (fields that are not in xml and should be calculated)
if ($param->{add})
  {        
    $data = &addField2data ($data, $param->{add});      
  }


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

elsif ($param->{action} eq "logodd")
  
  {
  	#no bins no logodds
  	$bin = &checkDataBin ($data, $bin);
  	#&bins2logodd (&hash2hashCsvLike ($bin->{'dibin'},"1"), &hash2hashCsvLike ($bin->{'dibin_T'}, "0"));	  	
  	$logodd = &bins2logodd ($bin, $param->{'outLogodd'});
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
        print STDERR "\n****ERROR: $param->{output} is an unknown parameter[FATAL]***\n";
        die;
      }
  }

if ($bin && $param->{outBins} ne "no")
  {  
  	printBins ($bin->{'bin'}, $bin->{'bin_T'}, "$param->{outBins}.bin");
		printBins (&hash2hashCsvLike ($bin->{'dibin'}, "1"), &hash2hashCsvLike ($bin->{'dibin_T'}, "0"), "$param->{outBins}.dibin");
  }

if ($logodd && ($param->{outLogodd} ne "no" && $param->{outLogodd} ne "screen"))
  {
  	&dumpLogodd ($logodd, $param->{outLogodd});
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
    $rp->{add} = 1;
    $rp->{bin} = 1;
    $rp->{countBin} = 1;
    $rp->{infile_format} = 1;
    $rp->{outBins} = 1;
    $rp->{outLogodd} = 1;
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
	  
#	  if (!$param->{outBins})
#	  	{
#	  		$param->{outBins} = "$param->{out}";
#	  	}
	  		
	  if (!$param->{outBins})
	  	{
	  		$param->{outBins} = "$param->{out}";
	  	}	
	  
	  if (!$param->{outLogodd})
	  	{
	  		$param->{outLogodd} = "$param->{out}";
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
             
            else 
            	{
            		print STDERR "\nERROR: $ff does not exist [FATAL]\n";
            		exit(1);
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

    vfopen ($F, $file);
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
    $file = &path2fileName ($file);            
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
                    $v =~ s/,/\./g; #decimals always point separated                    
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
                          $v =~ s/,/\./g; #decimals always point separated                          
                          $r->{$key} = $v; 
                        }
                                    
                    }
                  else
                    {                      
                      $key= $k_1."#".$k_2;
                      $v = $H->{$k_1}{$k_2};
                      $v =~ s/,/\./g; #decimals always point separated 
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
        foreach $k_2 (sort {$a cmp $b} keys(%{$H->{$k_1}}))
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
    
    foreach $f (sort ({$a cmp $b} keys(%$data)))#Review there is a new function which perform this appart I could use it here checkField 
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
                        (($v <= 45 && $v >0) || $v > 315 && $v <=360 ) && do 
                          { 
                            $bin = 1;                            
                            last SWITCH;
                          };
                         
                        ($v <= 135 && $v > 45) && do 
                          { 
                            $bin = 2;                            
                            last SWITCH;
                          };                                                                           
                            
                        ($v <= 225 && $v > 135 ) && do 
                          { 
                            $bin = 3;                            
                            last SWITCH;
                          };
                            
                        ($v <= 315 && $v > 225 ) && do 
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
                
                elsif ($field =~ "direction")
                	{
                		#&direction2bin ($d);#binDirection
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
    
    my ($k_1, $k_2, $k_3, $v, $bin, $f, $bin_T, $dibin, $dibin_T, $h);
    
    foreach $k_1 (keys(%$d))
      {
      	my ($pv);
      	
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
                  $bin->{$f}{'total'}++;
                  
                  #$bin->{$f}{'bin'}{$v}++;
                  
                  $bin_T->{$v}++;                  
                  
                  if (!defined ($pv))
                    {                      
                      $pv = $v;                      
                      next;
                    } 
                    
                  $dibin->{$f}{$pv}{$v}++;
                  $dibin->{$f}{'total'}{'transitions'}++;
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
    print $F "file";
    
    foreach $bin ((sort {$a cmp $b} keys (%$binList)))
    	{
      	print $F "\t$bin";
      } 
    
    foreach $f ((sort {$a cmp $b} keys (%$h)))
      {          
        if (exists ($h->{$f}{'total'}))
    			{
    				print STDERR "IWH";
    				print $F "\ttotal\n";
    				last;
    			}
    
		    elsif (exists ($h->{$f}{'total#transitions'}))
		    	{
		    		print STDERR "IWH";
		    		print $F "\ttotal\n";
		    		last;
		    	}
		    else
		    	{
		    		print $F "\n";
		    		last;
		    	}
      }
      
    foreach $f ((sort {$a cmp $b} keys (%$h)))
      {                  		    
        print $F "$f";
                        
        #foreach $bin (@bin_angle)
        foreach $bin ((sort {$a cmp $b} keys (%$binList)))
          {
            print $F "\t";
            $v = $h->{$f}{$bin};
            
            if (!defined ($v)) 
              {
                $v = 0;
              }
            
            #printing frecuencies
            elsif (exists ($h->{$f}{'total'}))
            	{
            		#print "$h->{$f}{'total'} ----------------\n";#del
            		$v /= $h->{$f}{'total'};
            	}
            
            elsif (exists ($h->{$f}{'total#transitions'}))
            	{
            		#print "$h->{$f}{'total#transitions'} ----------------\n";#del
            		$v /= $h->{$f}{'total#transitions'};
            	}
              
            printf $F  "%6.3f", $v	;    
          }
        
        if (exists ($h->{$f}{'total'}))
    			{    				
    				print $F "\t$h->{$f}{'total'}";
    			}
    
		    elsif (exists ($h->{$f}{'total#transitions'}))
		    	{
		    		print $F "\t$h->{$f}{'total#transitions'}";		    		
		    	}  
		    	
        print $F "\n";  
      }
    
  }

#			 $h->{'bin'} = $bin;
#      $h->{'bin_T'} = $bin_T;
#      $h->{'dibin'} = $dibin;    
#      $h->{'dibin_T'} = $dibin_T; 

sub bins2logodd
	{
		my $h = shift; 
    my $p = shift;
    
    my ($bin, $binList, $dibins, $id, $b1, $b2, $N_tr, $N_1, $N_2, $fr_tr, $fr_1, $fr_2, $logodd, $H_logodd);
    
    #Using the whole hash which summarizes the bins counts of all files, some files might not have all bins
    $bin = $h->{'bin'};
    $binList = $h->{'bin_T'};  
    $dibins = $h->{'dibin'};
    
    foreach $id ((sort {$a cmp $b} keys (%$dibins)))
    	{
      	if ($p eq "screen")
      		{
      			print "\n\n$id=================\n\n";
      		}
      		
      	foreach my $b1 ((sort {$a <=> $b} keys (%{$binList})))
	  			{	  			
	    			foreach my $b2 ((sort {$a <=> $b} keys (%{$binList})))
	      			{	      				
	      				$N_tr = $dibins->{$id}{$b1}{$b2};	      					      			
	      				$N_1 = $bin->{$id}{$b1};
	      				$N_2 = $bin->{$id}{$b2};	      				
	      				
	      				$fr_tr = $N_tr / $dibins->{$id}{'total'}{'transitions'};
	      				$fr_1 = $N_1 / $bin->{$id}{'total'};
	      				$fr_2 = $N_2 / $bin->{$id}{'total'};
	      				
	      				$logodd = (($fr_1*$fr_2)==0 || $fr_tr==0)? 0 : log (($fr_tr) / ($fr_1*$fr_2));	      				
	      				
	      				$H_logodd->{$id}{$b1}{$b2} = $logodd;   
	      				
	      				if ($p eq "screen")
	      					{
	      						printf "\tID: %10s \t %1s -- %1s : %6.3f (N Trans: %5d)(N Tot Trans: %5d)(FC: %5d)($b1: $N_1, $b2: $N_2)\n", $id, $b1, $b2, $logodd, $dibins->{$id}{$b1}{$b2}, $dibins->{$id}{'total'}{'transitions'}, $bin->{$id}{'total'}; 	      						
	      					}
	      			}
	  			}
      }
      
      return ($H_logodd);
    	
	}

sub dumpLogodd
	{
		my $h = shift;
		my $fileName = shift;
		
		my ($ls);			
		
			
		$h = &hash2hashCsvLike($h, "1");
		
		foreach my $f (keys (%$h))
			{
				$ls = $h-> {$f};
				last; #all of them have all the possible pairs of dibins
			}
		
		printBins ($h, $ls, "$fileName.logodd"); 		
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

sub addField2data
	{
		my $data = shift;
		my $p = shift;
		
		my ($f, $h, $k_1, $window, $angle);
		
		$window = 1; #that means we use the previous and following points
		
		if ($p eq "direction")
			{
				&checkField ($data, "Xabs", "Yabs");
				
				foreach $f (sort ({$a cmp $b} keys(%$data)))
					{						
						$h = $data-> {$f}{'data'};
						
						my (@lastDiffPoint, $flag);
						
						foreach $k_1 (sort ({$a <=> $b} keys(%$h)))
		          { 
	
		          	if ($k_1 < $window)
		          		{		          			      			
		          			$h->{$k_1}{'animalPosition#direction'} = 0;		          			
		          			next;	
		          		}
		          			          			          	
		          	elsif (exists ($h->{$k_1+1}) && exists ($h->{$k_1+1})) #elsif (exists ($h->{$k_1+$window}) && exists ($h->{$k_1+$window})) No tendria que ser asi
		          		{
										
										my @pv = ($h->{$k_1-$window}{'animalPosition#Xabs'}, $h->{$k_1-$window}{'animalPosition#Yabs'});
										
										my @v = ($h->{$k_1}{'animalPosition#Xabs'}, $h->{$k_1}{'animalPosition#Yabs'});
										
										my @fv =  ($h->{$k_1+$window}{'animalPosition#Xabs'}, $h->{$k_1+$window}{'animalPosition#Yabs'});
																				
										if (&aryCompare (\@v, \@pv))#if it remains in the same previous position
											{
												 $h->{$k_1}{'animalPosition#direction'} = 0;												 											
												 next;
											}
                                         	
                                        elsif (&aryCompare (\@v, \@fv))#if it remains in the same following position
										{
											 $h->{$k_1}{'animalPosition#direction'} = 0;
											 
											 #I have to keep the last position that was different to make the calculation when it changes
											 next if ($flag);
											 $flag = 1;												 
											 @lastDiffPoint = @fv;												 
											 next;
										}
																																
										else
											{
																								
												if ($flag) 
													{
														@pv = @lastDiffPoint;
														$flag = 0; 													
													}
																				
												my $a = &substrVec (\@v, \@pv);
												my $b = &substrVec (\@fv, \@v);
												
												my $c = (&scalarDot ($a, $b)) / (&modulus ($a) * &modulus ($b)); 

												$angle = &acos ($c);	
												$angle = &rad2deg ($angle);		
												
												#As the angle to respect the target it angle 0 corresponds to y=0 and the angle grows in anticlockwise sense
												#$angle = ($b->[0] > $a->[0])? 360-$angle : $angle;
												
												#movements to the left negative angles until -180 to the right positive until 180												
												$angle = ($b->[0] > $a->[0])? $angle : -$angle;
																					
												$h->{$k_1}{'animalPosition#direction'} = $angle;												
											}
										
		          		}
		          		
		          	else 
		          		{
										$angle = 0;
                    $h->{$k_1}{'animalPosition#direction'} = $angle;
									}								
		          }	
		          $data-> {$f}{'data'} = $h;	          
					}										
			}	
		return ($data);
	}
	
sub checkField
	{
		
		my $h = shift;
		my @list = @_;
		
		my ($f, $k_1, $k_2, $field);
		
		foreach $field (@list)
			{
				my $flag;
				
				foreach $f (sort ({$a cmp $b} keys(%$data)))
		      {
		        
		        $h = $data-> {$f}{'data'};
		        
		        foreach $k_1 (sort ({$a cmp $b} keys(%$h)))
		          {    
		            foreach $k_2 (sort {$a cmp $b} keys (%{$h->{$k_1}}))
		              {                
		                next if ($k_2 !~ /($field)/);
		                $flag = 1;
		                last;
									}
								
								if ($flag==1) 
									{
										last;
									}
																	
								else 
									{
										print STDERR "\n****ERROR: $field is not present at data [FATAL]***\n";
                    die;
									}													
		          }
		      }
			}
	}
	
sub substrVec
	{
		my $v1 = shift;
		my $v2 = shift;
		my @v3;
		
		if (scalar(@$v1) == scalar (@$v2))
			{
				for (my $i = 0; $i < scalar(@$v1); $i++)
					{
						$v3[$i] = $v1->[$i] - $v2->[$i];
					}
				
				return (\@v3);
			}
		else
			{				
				print STDERR "\n****ERROR: Vectors $v1 - $v2 have different size substraction fail [FATAL]***\n";
				die;
			}
	}

sub scalarDot
	{
		my $v1 = shift;
		my $v2 = shift;
		
		my $sc = 0;
		
			if (scalar(@$v1) == scalar (@$v2))
				{
					for (my $i=0; $i < scalar(@$v1); $i++ )
						{
							$sc += $v1->[$i] * $v2->[$i];		
						}
					
					return ($sc);
				}
			
			else
				{
					print STDERR "\n****ERROR: Vectors have different size scalar product fail [FATAL]***\n";
					die;
				}
	}
	
sub modulus
	{
		my $v = shift;
		
		my $mod = 0;
		
		for (my $i=0; $i < scalar(@$v); $i++ )
			{
				$mod += $v->[$i] ** 2;
			}
		
		return (sqrt ($mod));	
							
		
	}

sub acos 
	{ 
		my $n = shift;
		my $res;
		
		$n = ($n > 1)? 1 : $n; 
		
		$res = atan2 (sqrt (1 - $n * $n), $n);
		
		return ($res);  
	}
	
sub aryCompare
	{
		my $ary1 = shift;
		my $ary2 = shift;
		
		if (scalar (@$ary1) == scalar (@$ary1))
			{
				for (my $i=0; $i < scalar(@$ary1); $i++ )
					{
						if ($ary1->[$i] != $ary2->[$i])
							{
								return (0);
							}
					}
				
				return (1);
			}
		
		else
			{
				return (0);
			}	
	}

	sub rad2deg 
		{ 
			my $alpha = shift;
			my $alpha_deg = shift;
			 
			$alpha_deg = ($alpha/$pi) * 180; 
		}
	
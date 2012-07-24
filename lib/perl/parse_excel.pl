#!/usr/bin/perl -w
use strict;
use Data::Dumper;
use FileHandle;

our ($path,$outpath);
our %mice;
our @types;

my @measures;
my (%xt,%x);
my ($P_zones,$P_quants,$P_measures,$P_xt,$P_x);#pointers


$path="/users/cn/ierb/work/MaraDierssen/data/Ts1Ts2Ts1Cje/MWM/";
#$outpath="/users/cn/ierb/work/MaraDierssen/results/matrices_from_excel/contingency_matrices/equal_binning/";

@types=("Wt","Ts1CjeTs1","Ts1CjeTs2");

@{$mice{"Wt"}}=(173,178,179,181,182,183,184,904,909,914,924,932,937);
@{$mice{"Ts1CjeTs1"}}=(161,171,185,187,193,194,197);
@{$mice{"Ts1CjeTs2"}}=(905,908,910,911,913,916,918,922,923,934,939);

($P_xt,$P_zones,$P_quants)=parse_excel();

#print Dumper ($P_xt);die;#del

%xt=%$P_xt;#$xt{$id}{$quant.$zone}[$day] ~ $xt{$mouse}{$measure}[$day]

($P_x,$P_measures)=trial_averages($P_xt,$P_zones,$P_quants);

%x=%$P_x;#$x{$mouse}{$measure}[$day] #average over trial per mouse, i.e. the information of one day containing 4 trials
@measures=@$P_measures;

#Printing the table with the mean measure of a given variable over trials for each mouse (rows)
#&printTableVar_vs_Mice (\%x);
&printTableByVar_vs_Mice (\%x);

#Printing the table with the mean value of a given measure over trials for each mouse (rows)
#                           V1  V2  V3 ....  
#  mice n
#  mice n+1
#  ...
#The input of function is a %hash with this structure $x{$mouse}{$measure}[$day]
sub printTableVar_vs_Mice
  {
    my ($h, $k_1, $k_2, $genotype, $day);
    $h = shift;    
    #print STDERR "Are you talking to me!\n";#del
    
    #print table header
    for ($day=1;$day<=5;$day++)
      {
      
        foreach $k_1 (sort ({$a <=> $b} keys(%$h)))
          { 
            if ($day == 1)
              {
                print "animal\tgenotype";
              }
                         
            foreach $k_2 (sort {$a cmp $b} keys (%{$h->{$k_1}}))
              {
                #We don't like "(%)" symbol in headers
                $k_2 =~ s/\(%\)/Perc/g;
                print "\t$k_2"."Day"."$day";
              }
              
             #print "\n";
             last;
          }
        } 
    print "\n";
    
    #print values 
    
    #animal
    foreach $k_1 (sort ({$a <=> $b} keys(%$h)))
      {
        
        print "$k_1";
    
        #which genotype has the animal
        #due to the structure of the hash I have to search for each of the genotype whether the animal is present or not
        foreach $genotype (keys (%mice)) 
          {             
            (searchValInAry ($k_1, $mice{$genotype}) == 1)? print "\t$genotype" : next;            
          }         
                          
        for ($day=1;$day<=5;$day++)
          {
            foreach $k_2 (sort {$a cmp $b} keys (%{$h->{$k_1}}))
              {                       
                my $ary = $h->{$k_1}{$k_2};
                print "\t@$ary[$day]"; #I was only printing day four
              }
          }
          
        print "\n"; 
          
      }
    
      
      
      #print "\n";
      
    #die;
  }
  
#Printing the table with the mean value of a given measure over trials for each mouse (rows)
#                           V1  V2  V3 ....  
#  mice n
#  mice n+1
#  ...
#The input of function is a %hash with this structure $x{$mouse}{$measure}[$day]
#Modify to make the plots that they use to see difference across days by groups
#
#   - -
#   * *     
#        - - -   
#        * * *
#                *  
#                - 
sub printTableByVar_vs_Mice
  {
    my ($h, $m, $v, $genotype, $day);
    $h = shift;    
          
    foreach $v (sort {$a cmp $b} keys (%{$h->{'910'}})) #variables
      {                                               
        my $file = $v."_linePlot.tbl";
        my $F= new FileHandle;
        open ($F, ">$file");
        
        print STDERR "$v\n";
        
        #Header
        print $F "animal\tgenotype\tday1\tday2\tday3\tday4\tday5\n";
        
        #animal, i.e.mice
        foreach $m (sort ({$a <=> $b} keys(%$h)))
          { 
            my $ary = $h->{$m}{$v};
            
            #print "$m";
            print $F "$m";
            
            #Printing genotypes
            foreach $genotype (keys (%mice)) 
              {             
                (searchValInAry ($m, $mice{$genotype}) == 1)? print $F "\t$genotype" : next;  #guardarlo para printar          
              } 
               
            for ($day=1;$day<=5;$day++)
              {                               
                print $F "\t@$ary[$day]"; 
              }
              
            print $F "\n";
          }              
        close ($F);
        
      }

  }
  
sub trial_averages{
    
    my ($P_xt,$P_zones,$P_quants)=@_;
    my %xt=%$P_xt;
    my @zones=@$P_zones;
    my @quants=@$P_quants;

    my ($quant,$zone,$day,$type,$mouse,$id,$measure,$trial);
    my @tmp;
    my ($P_x,$P_measures);

    foreach $quant (@quants){
	foreach $zone (@zones){
	    
	    if ($quant eq "Lat.T." && $zone eq "TOTAL"){
		$measure="whishaw";
	    }
	    else{
		$measure=$quant.$zone;
	    }
	    push @measures,$measure;
	    
	    for ($day=1;$day<=5;$day++){
		foreach $type (@types){
		    foreach $mouse (@{$mice{$type}}){ 
		    
			$x{$mouse}{$measure}[$day]=0;
			@tmp=();
			foreach $trial ("N","S","W","E"){
			    $id=$mouse." ".$trial;
			    
			    if (defined $xt{$id}{$measure}[$day]){
				
				push @tmp,$xt{$id}{$measure}[$day];
				
			    }
			    else{
				$xt{$id}{$measure}[$day]="NA";
				#print "undefined $id day$day\n";
			    }
			}
			if (@tmp==0){
			    die "don't have $mouse $day\n";
			}
			else{
			    $x{$mouse}{$measure}[$day]=average(@tmp);
			}
			
		    }
		}
	    }
	}
    }
    return (\%x,\@measures);
}



sub parse_excel{
    
    my ($day,$file,$id,$zone,$i);
    my (@tmp,@variable);
    my (@zones,@quants);
    my %xt;

    for ($day=1;$day<=5;$day++){
	
	$file="A$day";
	
	

	open (IN,$path.$file."/$file".".txt") or die $path.$file."/$file".".txt\n";
	while (<IN>){
	    chomp;
	    @tmp=split(/\t/,$_);
	    if ($_=~/Subject Identification/){
	
		$tmp[1]=~s/\s$//;#delete possible ending in space
		$id=$tmp[1];
		
	    }
	    elsif($_=~/^Sbj.Code/){
		@variable=@tmp;
		
	    }
	    elsif(defined $variable[0] && $_=~/^$id/ && $tmp[1] ne "Zone 31"){
		$zone=$tmp[1];
		$zone=~s/\s//g;
		unless (@zones==15){
		    push @zones,$zone;
		}
		
		for ($i=2;$i<12;$i++){
		    $variable[$i]=~s/\s//g;
		    unless (@quants==10){
			push @quants,$variable[$i];
			
		    }
		    
		    $xt{$id}{$variable[$i].$zone}[$day]=$tmp[$i]*1;
		    
		}
		
	    }	
	    elsif($_=~/Error/){
		
		$_=readline (IN);
		@tmp=split(/\t/,$_);
		$xt{$id}{"whishaw"}[$day]=$tmp[0]*1;
		
	
	    }
		
	}
	close IN;
    }

    return (\%xt,\@zones,\@quants);
    

 
}

		
sub average{
    my @data=@_;
    my $i;
    my $x=0;

    
    for ($i=0;$i<@data;$i++){
	$x+=$data[$i];
    }
    
    $x/=scalar(@data);
    return $x;

}

#This function performes a search of a value in an array
#returns 1 or 0
sub searchValInAry
  {
    my $element = shift;
    my $ary = shift;

    if (grep {$_ eq $element} @$ary) 
      {      
        return (1);
      }
    else
      {
        return (0)
      }  
  }
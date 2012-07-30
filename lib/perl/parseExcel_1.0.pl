#!/usr/bin/perl -w
use strict;
use Data::Dumper;
use FileHandle;
use Cwd; #PWD function

our ($path,$outpath);
our %mice;
our @genotypes;

my @measures;
my (%xt,%x);
my ($P_zones,$P_quants,$P_measures,$P_xt,$P_x);#pointers


#$path="/users/cn/ierb/work/MaraDierssen/data/Ts1Ts2Ts1Cje/MWM/";
#$outpath="/users/cn/ierb/work/MaraDierssen/results/matrices_from_excel/contingency_matrices/equal_binning/";

@genotypes=("Wt","Ts1CjeTs1","Ts1CjeTs2");

@{$mice{"Wt"}}=(173,178,179,181,182,183,184,904,909,914,924,932,937);
@{$mice{"Ts1CjeTs1"}}=(161,171,185,187,193,194,197);
@{$mice{"Ts1CjeTs2"}}=(905,908,910,911,913,916,918,922,923,934,939);

#my %sessionOrder = (
#                    '1' => 'PT',
#                    '2' => 'A1',
#                    '3' => 'A2',
#                    '4' => 'A3',
#                    '5' => 'A4',
#                    '6' => 'A5', 
#                   );

my %sessionOrder = (                    
                    '1' => 'A3', 
                   );                   
                   
($P_xt, $P_zones, $P_quants)= &parseFiles ();

%xt = %$P_xt;

#Check whether empty fields are present
($P_xt) = &checkEmptyFields ($P_xt, $P_zones, $P_quants);

#Average over trials (cardinals point entry into the pool)
($P_x) = &trialAverages ($P_xt);


##%x=%$P_x;#$x{$mouse}{$measure}[$day] #average over trial per mouse, i.e. the information of one day containing 4 trials
##@measures=@$P_measures;

#Printing the table with the mean measure of a given variable over trials for each mouse (rows)

&printAllVar ($P_x);

#Printing a table for each single variabla containing all the values for all sessions
#Classical plot
#&printTableByVar_vs_Mice ($P_x);

###################
#FUNCTIONS
##########

sub parseFiles 
  {
    my $ext = "txt";
    my ($file, $id, $zone, $i, $day, $mouse, $point);
    my (@mousePoint, @tmp,@variable);
    my (@zones,@quants);
    my (%xt, %ids);
    
    my $pwd = getcwd;
    
    opendir DIR, $pwd or die "cannot open dir $pwd: $!";
    
    for (readdir DIR)
      {
        if (-d $_) {next;} #skip directories
        elsif (/^[.]/) {next;} #skip dot-files
        
        #elsif (/(.+)[.]\Q$ext\E$/)
        elsif (/(.+)[.]$ext$/)
          { 
            #$day = $1;
            #print "text file: ", $1;
            $file =  ($_);
                        
            if ($file =~ /(.+)[_](.+)$ext$/)
              {
                $day = $1;
              }
              
            else
              {
                $day = $1;
              }
            
            $day = uc ($day); #forcing upper case
            
            print STDERR "Processing file: ".$file."\n";
                        
            my $F=new FileHandle;
            
            %ids = {};
            
            open ($F, "$file")  or die "FATAL ERROR: Can't open file: $file";
                        
            while (<$F>)
              {
                                
          	    chomp;
          	    @tmp=split(/\t/,$_);
          	    
          	    if ($_=~/Subject Identification/)
          	     {          	
            		$tmp[1]=~s/\s$//;#delete possible ending in space
            		$id=$tmp[1];
            		
            		#Checking whether same id has been already used
            		if (exists ($ids {$id}))
            		  {
            		    print STDERR "FATAL ERROR: Animal Id $id has already been used in file --> $file\n";
            		    die;
            		  }
            		else
            		  {
            		    $ids {$id} = 1;
            		  }  
            		
            		#Checking animal id (format as example: 173 E)
            		if ($id !~ /^\d{3}\s+\w/) 
            		  {
            		    print STDERR "WARNING: Animal Id $id is not standardize file --> $file\n";
            		  }
            		   
            		@mousePoint = split ( " ", $id);          		
            		$mouse = $mousePoint[0];
            		$point = $mousePoint[1];          		          	
          	     }
          	     
          	    elsif($_=~/^Sbj.Code/)
          	     {
          		  @variable=@tmp;          		
          	     }
          	    
          	    #zone 31 is the goal, depends on how much time spent the experimenter before retrieving the animal
          	    elsif(defined $variable[0] && $_=~/^$id/ && $tmp[1] ne "Zone 31")
          	     {
          		  $zone=$tmp[1];
          		  $zone =~ s/\s//g;
          		  $zone = uc ($zone); #forcing upper case
          		  
          		  unless (@zones==15)
          		    {
          		      push @zones,$zone;
          		    }
          		
          		  for ($i=2;$i<12;$i++)
          		    {
            		  $variable[$i]=~s/\s//g;
            		  unless (@quants==10)
            		    {
            			 push @quants,$variable[$i];          			
                        }
            		    
            	      #$xt{$id}{$variable[$i].$zone}{$day} = $tmp[$i] * 1;
            	      $xt{$mouse}{$variable[$i].$zone}{$day}{$point} = $tmp[$i] * 1;          		   
          		    }          		
          	       }	
          	    
          	    elsif($_=~/Error/)
          	     {          	
          		    my $dumpLine =  <$F>;
          		
          		    @tmp=split(/\t/,$_);
          		    $xt{$mouse}{"whishaw"}{$day}{$point} = $tmp[0] * 1;
			
	             }
			                                                            
              }
           
            close ($F);
            	        	                                  
          }
          
          else 
            {
              print STDERR "WARNING: File: $_ has an unsupported format\n"
            }
                
      }
       
    closedir DIR;
    
    return (\%xt,\@zones,\@quants);
  }

sub trialAverages
  {
    #my ($P_xt, $P_zones, $P_quants) = @_;
    my ($P_xt) = shift;
    
    my ($mouse, $var, $session, $point, $trial);
    
    #my ($quant, $zone, $measure, $type, $trial, $id );
    my (@tmp);
    

      
    foreach $mouse (sort ({$a cmp $b} keys (%$P_xt)))
      {        
        foreach $var (keys (%{$P_xt->{$mouse}}))
        
          {            
            foreach $session (keys (%{$P_xt->{$mouse}{$var}}))
			 {
			   foreach $point (keys (%{$P_xt->{$mouse}{$var}{$session}}))
			     {
			       if ($P_xt->{$mouse}{$var}{$session}{$point} ne "NA") 
			         {
			           push @tmp, $P_xt->{$mouse}{$var}{$session}{$point};
			         }
			       else
			         {
			           next;
			         }
			     }
			   
			   if (scalar(@tmp) == 0)
			     {
			       print "culo";
			     }
			   $x {$mouse}{$var}{$session} = average(@tmp); 
			         
               @tmp = ();      
			 }
            
          }
      }
    
    return (\%x);  
    
  }

sub checkEmptyFields
  {
    my ($P_xt, $P_zones, $P_quants) = @_;
    
    my ($quant, $zone, $measure, $type, $mouse,  $var, $k_1, $session, $trial);
    
    #This array will be used to check whether all animals have all the variables defined or not
    foreach $quant (@$P_quants)
      {
        foreach $zone (@$P_zones)
          {
            if ($quant eq "Lat.T." && $zone eq "TOTAL")#This field is always zero, that is why we use it us a dump variable
              {
                $measure = "whishaw"; #wishaw has not a zone
              }
            else
              {
    	        $measure = $quant.$zone;
              }
            
            push @measures,$measure;
          }
         
      }
      
    foreach $type (@genotypes)
      {
        foreach $mouse (@{$mice{$type}})
          {
            foreach $var (@measures)
			  {
                foreach $k_1 (sort ({$a <=> $b} keys(%sessionOrder)))
                  {
                    $session = $sessionOrder {$k_1};
                    
                    foreach $trial ("N","S","W","E")
    		          {
    		            if (!defined $P_xt->{$mouse}{$var}{$session}{$trial}) #session is always defined, each individual file is a session
    			             {
    			               $P_xt->{$mouse}{$var}{$session}{$trial}="NA"; #if it does not exist we define it as empty
    				           print STDERR "WARNING: undefined value for $mouse, $var, $session, $trial\n";
    			             }    			      
    		          }
                  }
			  }
          }
      }
      
    return ($P_xt);  
  }
      
#Printing the table with the mean value of a given measure over trials for each mouse (rows)
#                           V1  V2  V3 ....  
#  mice n
#  mice n+1
#  ...
#The input of function is a %hash with this structure $x{$mouse}{$measure}[$day]
sub printAllVar
  {
    my $h = shift;
    
    my ($mouse, $var, $k_1, $session, $genotype);

    print "animal\tgenotype";
    
    foreach $mouse (sort ({$a <=> $b} keys(%$h)))        
     {
         foreach $var (sort {$a cmp $b} keys (%{$h->{$mouse}}))
          {
            #we keep the order of sessions by using the hash sessionOrder
            foreach $k_1 (sort ({$a <=> $b} keys(%sessionOrder)))
              {
                $session = $sessionOrder{$k_1};
 
                #We don't like "(%)" symbol in headers
                $var =~ s/\(%\)/Perc/g;
                print "\t$var"."Day"."$session";
              }     
          }
       last;   
     }
    print "\n";
    
    #print table values
    foreach $mouse (sort ({$a <=> $b} keys(%$h)))        
     { 
       #which genotype has the animal
       #due to the structure of the hash I have to search for each of the genotype whether the animal is present or not
       foreach $genotype (keys (%mice)) 
        {             
          (searchValInAry ($mouse, $mice{$genotype}) == 1)? print "\t$genotype" : next;            
        }
           
       foreach $var (sort {$a cmp $b} keys (%{$h->{$mouse}}))
        {
          #we keep the order of sessions by using the hash sessionOrder
          foreach $k_1 (sort ({$a <=> $b} keys(%sessionOrder)))
            {
              $session = $sessionOrder{$k_1};
              
              print "\t$h->{$mouse}{$var}{$session}";                  
            }     
        }
        
      print "\n";              
     }
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
    my ($h, $pwd, $newDir, $mouse, $var, $v, $m, $genotype, $k_1, $session);
    $h = shift;    
    
    $newDir = "singleVarTbl";    
    &createGoDir ($newDir);
        
    my @variables;
            
    foreach $mouse (sort ({$a <=> $b} keys(%$h)))        
     {
      foreach $var (sort {$a cmp $b} keys (%{$h->{$mouse}}))
        {
          push @variables, $var;
        }
       last;
     }      
    
    #print Dumper (@variables);#del
    
    foreach $v (@variables) 
      {
        my $file = $v."_linePlot.tbl";
        my $F= new FileHandle;
        open ($F, ">$file");
        
        print STDERR "INFO: Variable being processed: $v\n";
        
        #printing table header
        print $F "mice\tgenotype";
        
        foreach $k_1 (sort ({$a <=> $b} keys(%sessionOrder)))
          {
            $session = $sessionOrder {$k_1};
            print $F "\t$session";
          }
        
        print $F "\n";
        
        #printing table content
        foreach $m (sort ({$a <=> $b} keys(%$h)))
          {
            print $F "$m";
            
            #Printing genotypes
            foreach $genotype (keys (%mice)) 
              {             
                (searchValInAry ($m, $mice{$genotype}) == 1)? print $F "\t$genotype" : next;           
              }
             
            foreach $k_1 (sort ({$a <=> $b} keys(%sessionOrder)))
              {
                $session = $sessionOrder {$k_1};
                
                print $F "\t$h->{$m}{$v}{$session}";
              }
            print $F "\n";   
          }
        close ($F);
      }
                 
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

#This function creates a directory in the current path and enters it
sub createGoDir
  {
    my ($newDir, $pwd, $fullPath);
    
    $newDir = shift;
    $pwd = getcwd;
    $fullPath = $pwd. '/' . $newDir;
    
    #Checking whether dir exists
    if (-d $fullPath) 
      {
        print STDERR "WARNING: Directory $newDir already exists!\n";
      }
    else
      {
        mkdir ($fullPath) or die "FATAL ERROR: Directory $newDir can not be done!";
      }
    
    chdir ($fullPath) or die "FATAL ERROR: Directory $newDir can not be reach!";
    
    return (1);          
  }
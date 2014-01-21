#!/usr/bin/env perl

use strict;
use FileHandle;
use Data::Dumper;

my $file = shift (@ARGV);
#print "$file\n";

my $F=new FileHandle;
    open($F, $file);        

my $first = 0;

while (<$F>)
  {
    my $l=$_;
            
    if ($l=~/^Subject Identification:\s(\w+)/)
      {      
        if ($first == 0) {print "$1"; $first=1;}      
        else {print ",$1"};
      }
  }
     
print "\n";
   




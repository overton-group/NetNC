#!/usr/bin/perl -w

##
## Calculates the Functional Coherence Score for NetNC
##

# Copyright 2016 University of Edinburgh

# @author Ian Overton (i.overton@qub.ac.uk)
#
#    FCS.pl is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#        
#    FCS.pl is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#     
#    You should have received a copy of the GNU General Public License
#    along with FCS.pl in the file LICENSE.txt
#    If not, see <http://www.gnu.org/licenses/>.
#

use strict;


my ($genelist, $netnc, $out) = @ARGV;
die "usage [genelist] [NetNC results] [output]\n" unless $out;

open (I, $netnc) or die "cant open $netnc $!\n"; ;

my (%p, %d); # p-values and degree

while (<I>) {
	chomp();
	my @l = split (/\t/, $_);
	$p{$l[0]} += $l[2];
	$p{$l[1]} += $l[2]; 
	++$d{$l[0]}; 
	++$d{$l[1]}; 
}

close I;


open (A, "$genelist") or die "cant open $genelist $!\n"; 
open (OUT, ">$out") or die "cant open $out $!\n";

while (<A>) { 
	chomp(); 
	my $gene = $1 if (/(\S+)/);  
	if ($p{$gene}) { 
		my $np = $p{$gene} / $d{$gene};
		print OUT "$gene\t$p{$gene}\t$d{$gene}\t$np\n";
	 } else { 
		print OUT "$gene\t0\t0\t0\n";
	}
}

close A;
close OUT;

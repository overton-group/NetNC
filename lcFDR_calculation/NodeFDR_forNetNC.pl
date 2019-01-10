#!/usr/bin/perl -w

##
## This script identifies the minimum pFDR value assigned to a node (gene) by analysing NetNC output files
##

# Copyright 2018 Queen's University Belfast
#
# @author Ian Overton (i.overton@qub.ac.uk)
#
#    NodeFDR_forNetNC.pl is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    
#    NodeFDR_forNetNC.pl is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#   
#   You should have received a copy of the GNU General Public License
#   along with NodeFDR_forNetNC.pl in the file LICENSE.txt
#   If not, see <http://www.gnu.org/licenses/>.
#
# v1.1

use strict;

my ($netnc, $fdr, $fcs) = @ARGV;

die "\nusage: NodeFDR_forNetNC.pl [NNC_Edge_Scores] [NetNC_FDR_file] [NetNC_FCS_file] > NodeFDR_output_file.txt\n\n## This script analyses files produced by the NetNC software:\n## NetNC_EdgeScores has suffix .FG.id1_id2_nlp.txt\n## NetNC_FDR_file has suffix .FDR.txt\n## NetNC_FCS_file has suffix .FCS.txt\n" unless ($fcs);

# Read in the NetNC edge scores, assigning highest NetNC score to each node

open (IN, $netnc) or die "cant open $netnc $!\n";

my (%scor, %fdr, %deg, %fcs);

while (<IN>) { 
        chomp();
        my @line = split (/\s+/, $_);

        my $n1 = $line[0];
        my $n2 = $line[1];
        my $sc = $line[2];

        if ($scor{$n1}) {
                $scor{$n1} = $sc if ($sc > $scor{$n1});

        } else {
                $scor{$n1} = $sc;
        }

        if ($scor{$n2}) {
                $scor{$n2} = $sc if ($sc > $scor{$n2});

        } else {
                $scor{$n2} = $sc;
        }
}

close IN;

# Read in the FDR file
open (FIN, $fdr) or die "cant open $fdr $!\n";

while (<FIN>) {
        chomp();
        my @line = split(/\s+/, $_);
        my $scor = $line[0];
        $fdr{$scor} = $line[1];

}

close FIN;

# Read the FCS data
open (FC, $fcs)  or die "cant open $fcs $!\n";

while (<FC>) { 
        chomp();
        my @line = split (/\s+/, $_);

        my $nod = $line[0];
        $deg{$nod} = $line[2];
        $fcs{$nod} = $line[3];

}

# Assign FDR to each node based on the score and print the results

foreach my $node (keys %fcs) {

        my $a = 1; 
        my $sco = 0;

        if ($scor{$node}) { 
                $a = $fdr{$scor{$node}};
                $sco = $scor{$node};
        }

        print "$node\t$deg{$node}\t$fcs{$node}\t$sco\t$a\n";

}


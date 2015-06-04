# Copyright (C) 2014 Alejandro Ochoa, Singh Research Group, Princeton University
# This file is part of dPUC.
# dPUC is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
# dPUC is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with dPUC.  If not, see <http://www.gnu.org/licenses/>.

package DpucOvsCompact;
our $VERSION = 1.00;
use lib '.';
use NetCC;
use strict;

# this package implements a rather messy bit of code that tries to compact overlap definitions so lp_solve is more efficient
# the idea is to find overlap cliques, which require fewer statements to declare (so they will be much more restrictive and otherwise efficient inside lp_solve)

# OLD description:
# try to compact overlap definitions
# idea is overlaps often appear in cliques, so we can make them a single statement
# this isn't only faster for feasibility verification, but might actually converge faster at the LP relaxation level, at least we hope...
# we won't worry about getting cliques exactly right, let's just get them greedily and hope for the best

sub compactOvs {
    # input is array ref of overlap pairs, each pair itself an array ref that contains two integer indexes (which uniquely identify the domains involved in the overlap)
    # returned is an array ref of compacted overlaps, each element an array ref that contains two or more integer indexes that are cliques, that is, every pair within that set is a disallowed overlap!
    my ($ovs) = @_;
    
    my @ovsCompact; # the data we want
    
    # this might seem backward, but doing it simplifies analysis...
    # load all edges into a double hash net (both ways)
    my %ovsNetBothWays;
    foreach my $ov (@$ovs) {
	$ovsNetBothWays{$ov->[0]}{$ov->[1]} = 1;
	$ovsNetBothWays{$ov->[1]}{$ov->[0]} = 1;
    }
    # now use this code to separate into connected components
    my ($CCs) = NetCC::getCCs(\%ovsNetBothWays);
    # process each connected component separately
    # we use a while loop because we may throw things to the back of the @$CCs list (a regular for/foreach won't work, I think)
    while (@$CCs) { # while there are things to process
	# get the next connected component
	my $CC = shift @$CCs;
	# need to make it both ways for degrees analysis...
	my $CCbothWays = makeAllEdgesBothWays($CC);
	# get number of nodes and edges in CC
	my $numNodes = scalar keys %$CCbothWays;
	my $numEdges = countEdgesInNet($CC); # has to be one way for accurate counts
	# take care of easy cases
	# the number of edges tells us if it's a full clique (covers the most common case, numNodes == 2)
	if ($numEdges*2 == $numNodes*($numNodes-1)) {
	    push @ovsCompact, [keys %$CCbothWays]; # all nodes are the solution clique we wanted.
	}
	else {
	    # else continue, finding greedy cliques will take a bit more work...
	    # find the node with the highest degree...
	    my $nodeBest = getNodeHighestDegree($CCbothWays);
	    # clique is a list, initialize with this node
	    my @clique = ($nodeBest);
	    # greedily grow clique by visiting $nodeBest's neighbors
	    # this is a random process which may have multiple correct answers
	    foreach my $neigh (keys %{$CCbothWays->{$nodeBest}}) {
		my $fullyConnected = 1; # remember if we found disconnected pairs
		my $neighsOfNeigh = $CCbothWays->{$neigh}; # copy down neighbors of this neighbor (net around "clique" candidate)
		foreach my $nodeInClique (@clique) {
		    unless (exists $neighsOfNeigh->{$nodeInClique}) { # this means $neigh and $nodeInClique aren't connected!
			$fullyConnected = 0; last; # stop looking if we know this isn't a good candidate
		    }
		}
		push @clique, $neigh if $fullyConnected; # add $neigh to clique if it was fully-connected to previous clique members
	    }
	    # sort cause we use it sorted later
	    @clique = sort @clique;
	    # clique is now maximal, add to list we want
	    push @ovsCompact, \@clique;
	    # since we already covered the CC==clique case, in this case there are edges leftover that are not in clique
	    my $sizeClique = scalar @clique;
	    # navigate all pairs in clique to delete, in canonical order (nodei lt nodej since @clique is sorted)
	    for (my $i = 0; $i < $sizeClique; $i++) {
		my $nodei = $clique[$i];
		for (my $j = $i+1; $j < $sizeClique; $j++) {
		    my $nodej = $clique[$j];
		    # delete from "one way" net
		    # however, note only keys in the second level are being deleted currently
		    delete $CC->{$nodei}{$nodej};
		}
		# this makes sure we delete keys in the first level IF their second level hash is empty
		# this works because in this order we should have deleted all second-level entries that were to be deleted
		delete $CC->{$nodei} unless scalar keys %{$CC->{$nodei}};
	    }
	    # now we should have a clean $CC (one way net)
	    # we might have more cliques (or trivial pairs) left to be processed, so try to separate into CCs once more and append to queue
	    $CCbothWays = makeAllEdgesBothWays($CC); # need to repeat this on pruned net
	    my ($CCs2) = NetCC::getCCs($CCbothWays); # partition again
	    push @$CCs, @$CCs2; # add new partitioned things to queue (could be just one, we'll worry about that in next iteration)
	    # done with this original $CC!
	}
    }
    # very last thing is to reorder the overlaps so the largest (most restrictive) sets come first
    @ovsCompact = sort { @$b <=> @$a } @ovsCompact;
    
    return \@ovsCompact;
}


### INTERNAL FUNCTIONS

sub makeAllEdgesBothWays {
    # normally we keep the undirected network so that edges satisfy $acc1 lt $acc2, but for searching purposes we want both directions
    # this makes a new graph, so the old graph is not modified
    # this works on multinets as well as regular nets ($value is simply a reference instead of a scalar)
    # assumes input net is undirected.  Inputting a directed graph will produce a graph that is still symmetric in the end, but exactly which of the two values is kept for the undirected edge is decided randomly.
    my $net = shift;
    my %net; # net net with edges both ways
    # navigate input network
    while (my ($acc1, $acc2s) = each %$net) {
	while (my ($acc2, $value) = each %$acc2s) {
	    # transfer values to new network
	    $net{$acc1}{$acc2} = $value;
	    $net{$acc2}{$acc1} = $value;
	}
    }
    return \%net; # return new net
}

sub getNodeHighestDegree {
    # very fast if we only want the top node (random choice for ties)
    my ($netBothWays) = @_;
    # initialize the data we want
    my $nodeBest = '';
    my $degBest = 0;
    # navigate nodes in a random order
    while (my ($node, $netHere) = each %$netBothWays) {
	# degree is simply number of neighbors!
	my $deg = scalar keys %$netHere;
	# overwrite "best" if necessary
	if ($degBest < $deg) {
	    $degBest = $deg;
	    $nodeBest = $node;
	}
    }
    # done, return node we wanted
    return $nodeBest;
}

sub countEdgesInNet {
    # should work just fine with directional nets, or symmetric nets that are one-way
    # net has to be one-way only!
    # this breaks with nets that are both ways, watch out!!!
    # actually divide by two if it's both ways, but IF there are no self edges.  Self edges really mess this up.
    my ($net) = @_;
    my $c = 0;
    while (my ($acc1, $acc2s) = each %$net) {
	$c += scalar keys %$acc2s;
    }
    return $c;
}

1;

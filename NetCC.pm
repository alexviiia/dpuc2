# Copyright (C) 2014 Alejandro Ochoa, Singh Research Group, Princeton University
# This file is part of dPUC.
# dPUC is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
# dPUC is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with dPUC.  If not, see <http://www.gnu.org/licenses/>.

package NetCC;
our $VERSION = 1.00;
use strict;

# package for implementing the search of connected components on a network (in my hash format)
# the main function has a very limited use, it makes sense to keep it separate from other Net code

sub getCCs {
    # identifies all connected components of the given network.  They are found randomly, and sorted in the end by their sizes (largest first)
    # each returned subnet in list of CCs satisfies $acc1 lt $acc2 (i.e. they're one way only)
    my ($netBothWays) = @_;
    # to know when we've retrieved all the connected components
    my @nodes = sort keys %$netBothWays;
    my $numNodes = scalar @nodes;
    # this array keeps the connected components, list of references to nets
    my @CCs;
    my @CCsNodes; # parallel array with set of nodes (sometimes we'd rather have that than the net)
    # this hash maps a reference for the connected component to its node size, so we don't have to recalculate when we're sorting
    my %ccRef2size;
    
    # try to exhaust nodes
    while (@nodes) { # there are still nodes we haven't visited
	# get a random node from the network
	my $query = shift @nodes;
	# get the connected component of this node
	my ($netCC, $nodesCC) = findCCOneNode($netBothWays, $query);
	# store subnet in list of connected components
	push @CCs, $netCC;
	push @CCsNodes, $nodesCC; # nodes-only version (a set, not a list)
	# map reference to size in nodes
	$ccRef2size{$netCC} = keys %$nodesCC;
	# subtract the CC's nodes from the total nodes...
	@nodes = setMinusHash(\@nodes, \%$nodesCC);
    }
    
    # now we've split the entire network into connected components, sort list descending by node size and return
    @CCs = sort { $ccRef2size{$b} <=> $ccRef2size{$a} } @CCs;
    return (\@CCs, \%ccRef2size, \@CCsNodes);
}

### INTERNAL FUNCTIONS

sub findCCOneNode {
    # this does a breath-first search and returns all the subnetwork that is connected to the node
    # assumes the $netBothWays has been processed with "sub makeAllEdgesBothWays", so we can find all of a node's neighbors directly
    # returned subnet satisfies $acc1 lt $acc2 (i.e. it's one way only)
    # we also return a sorted list of all the nodes that are part of the subnet
    my ($netBothWays, $node) = @_;
    my %netCC; # connected component network
    my %visited; # remembers the nodes we've analized
    my %queue = ($node => 1); # nodes to examine, start with our query
    while (keys %queue) {
	# get current node to analyse
	my $query = (keys %queue)[0]; # get any key, we forgo order here
	delete $queue{$query}; # remove from queue
	# mark it as analyzed
	$visited{$query} = 1;
	# get all the neighbors of the query and their values
	while (my ($neighbor, $value) = each %{$netBothWays->{$query}}) {
	    # add to queue if we haven't analyzed it yet
	    $queue{$neighbor} = 1 if !$queue{$neighbor} && !$visited{$neighbor};
	    # add edge to new net (this may overwrite the end coming from the other direction, but it doesn't matter)
	    my ($acc1, $acc2) = sort ($query, $neighbor); # first sort so the network is good
	    $netCC{$acc1}{$acc2} = $value; # now add edge
	}
    }
    # return the subnetwork of interest, and also a list of nodes in this subnetwork (generated incidentally and may be useful for some applications)
    return (\%netCC, \%visited);
}

sub setMinusHash {
    # return an array that contains all the elements in $a that are not in hash $b
    my ($a, $b) = @_;
    my @c;
    foreach my $a0 (@$a) {
	push @c, $a0 unless exists $b->{$a0};
    }
    return @c;
}

1;


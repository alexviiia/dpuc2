# Copyright (C) 2014 Alejandro Ochoa, Singh Research Group, Princeton University
# This file is part of dPUC.
# dPUC is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
# dPUC is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with dPUC.  If not, see <http://www.gnu.org/licenses/>.

package DpucNetScores;
our $VERSION = 1.03;

# 2015-06-05 10:52:14 EDT
# v1.01 - script now parses files that start with list of nodes
# will incrementaly adjust scripts to make proper use of this new info

# 2015-08-05 13:37:44 EDT
# v1.02 - added direct alpha setting for prior (rather than through scale)
# - also obsoleted "net blending" option through dPUC (many simplifications!)
# - minor loop change (avoids using "each") stops newer Perls from complaining about it
# 2015-08-09 10:39:06 EDT
# - same version (haven't posted): added "old" panning option, that behaves like the old dPUC used to (didn't match any panning weight construct)
# 2015-08-11 20:18:41 EDT
# - same version (haven't posted): added "old" scale option, that behaves like the old dPUC used to (setting scale as 1/n^2, regardless of actual scale achieved).
# 2015-08-14 11:24:16 EDT
# - same version (haven't posted): removed all self/trans distinction and parameter that controlled it. the code that handled negative self-context scores also went away.
# - also removed panning option, added "oldMode" boolean that sets old panning and old scale at the same time.
# - also removed gamma option, which probably doesn't make much sense, and would be a pain to benchmark. This drops an internal function completely!
# - overall, cleaned internal counts2scores considerably!

# 2015-10-30 22:04:07 EDT
# v1.03 - removed "scale" prior parametrization in favor of directly setting alpha (which is also assumed to be linear, not log like before)

use lib '.';
use FileGz;
use EncodeIntPair;
use strict;

# for Dpuc only
sub loadNet {
    # this is a wrapper around the file that normally loads the net, includes the typical pre-processing that is necessary for posElim or nonOvNeg
    my ($scoreScale, $accs, $fi, $alpha, $scaleContext, $shiftContext, $cCut, $oldMode) = @_;
    
    # always load raw counts and compute scores on the fly!
    $cCut = undef if $oldMode; # old mode didn't have these things!!! Ensure cuts aren't used there!!!
    my ($net, $nodes) = parseNet($fi, $cCut);
    my $n = @$nodes; # counts2scores needs to know how many families there are!

    # as soon as we can, let's compare the outside @$accs with @nodes hash/set, they should be same or die!
    # this way, we catch a version mismatch as early as possible, though this is technically only between Pfam-A.hmm.dat and the dpucNet (input domain preds could be yet of a different version, so we'll be checking them again later as well.
    die "Fatal: domain family sets disagree between input Pfam-A.hmm.dat and this dpucNet: $fi\nTheir Pfam versions probably disagree!\nNo output was generated\n" unless setIdentity($accs, $nodes);
    
    # construct scores given counts and other parameters
    $net = counts2scores($net, $n, $alpha, $scaleContext, $shiftContext, $oldMode);
    
    # need to make sure this exists before processing nets (necessary for C code and to define $n for on-the-fly scores)
    my $acc2i = set2hashmap($nodes); # now make hash that maps each accession to its index, store as global!
    
    # when doing posElim, scale scores to ints, turn into parallel arrays, store in C space (so fast code access it fast)
    
    # scale to have integers, trying to speed up calculations
    my $netScaleInt = scaleEdgesAndInt($net, $scoreScale);
    # this turns double hash net into a parallel array that is a lot more C-friendly, but not perlish at all
    # additionally, this function removes non-positive scores
    my ($netTs, $netTscores) = mapPfamNetAccs2ints($netScaleInt, $acc2i);
    
    # posElim needs a list of domains that carry context information, also gets stored in C space
    # acc2cc is an hash ref that maps each acc to 1 if it is context carrying, undef otherwise
    # note acc2cc contains NEG too, doesn't affect out analysis at all!
    my $acc2cc = getNodes($net);
    
    # only nonOvNeg uses $net, return it (safe to ignore for posElim), also return the stuff needed for C structures (for posElim)
    return ($net, $acc2i, $netTs, $netTscores, $acc2cc); # stuff that Dpuc2 needs, many for passing to C
}

sub parseNet {
    my ($fi, $cut, $noNodeList) = @_;
    # loads a net from file, with a threshold on edge values
    # default is to assume file has a node list, but it's possible to say there isn't (for backward compatibility?)
    my $isCut = defined $cut; # default is no cut
    my %net;
    my @nodes; # list of nodes which might be present in the file
    my $fhi = FileGz::getInFh($fi);
    # first line is node list... (unless explicitly told otherwise)
    unless ($noNodeList) {
	# read one line, chomp
	$_ = <$fhi>;
	chomp;
	# get list of nodes, split by tabs
	@nodes = split /\t/;
    }
    while (<$fhi>) {
	chomp;
	my ($acc1, $acc2, $value) = split /\t/;
	$net{$acc1}{$acc2} = $value if !$isCut || $value >= $cut; # add only if cutoff is passed
    }
    close $fhi;
    return (\%net, \@nodes); # return these two things
}


### functions only used internally

sub setIdentity {
    # returns true if sets are identical (unordered array values, or if hash: same keys; values are ignored), false otherwise
    # this version assumes keys have no newlines, then comparison can be faster!
    my ($s1, $s2) = @_;
    $s1 = [keys %$s1] if 'HASH' eq ref $s1;
    $s2 = [keys %$s2] if 'HASH' eq ref $s2;
    # turn arrays into strings, such that strings are equal (assuming no newlines in keys) iff sets are equal.
    $s1 = join "\n", sort @$s1;
    $s2 = join "\n", sort @$s2;
    return $s1 eq $s2; # this comparison gives the answer
}

sub counts2scores {
    # copied and simplified from dpucNet04.counts2scores.pl, in order to generate one single net on the fly from a single set of counts
    # no combinatorics here like original script (which was left unchanged), and we dropped some useless features (notably clans map)
    # dropped intermediate panned counts net creation, made sense for combinatorics but not here!
    # on the other hand, made algorithms faster, which matters a lot more when doing things on the fly!
    my ($acc2acc2c, $n, $alpha, $scaleBits, $shiftBits, $oldMode) = @_; # input count structure and params
    # get sum of counts from network, parameter independent!
    my $cTot = sumEdges($acc2acc2c);
    # add counts in both directions (except for diagonal) if we want "old" mode.
    # Note: in old mode $cTot was as it is above (and it would be wrong to compute it after the transformation below)
    $acc2acc2c = makeOldDpucCounts($acc2acc2c) if $oldMode;
    
    # part 2: score params
    my $n2 = $n*$n; # the number of directed pairs appears often!
    # oldMode overrides alpha here!
    $alpha = 1/$n2 if $oldMode; # emulates dPUC 1.0 behavior
    # score super denominator includes denominator of pij probability (cTot + alphaSum), and pi and pj backgrounds, both 1/n
    my $superDenom = log2( $cTot/$n2 + $alpha );
    
    my %acc2acc2s; # scores we want!
    # navigate net
    foreach my $acc1 (keys %$acc2acc2c) {
	my $acc2c = $acc2acc2c->{$acc1};
	while (my ($acc2, $c) = each %$acc2c) {
	    # make score, store in net!
	    $acc2acc2s{$acc1}{$acc2} = log2($c + $alpha) - $superDenom;
	}
    }
    
    # add entries for the uniform negative scores
    $acc2acc2s{NEG}{NEG} = log2($alpha) - $superDenom;
    my $net = \%acc2acc2s; # save as reference, more convenient for following operations...
    
    # if we want to scale or shift context (relative to HMM scores), overwrite net with scaled net, this will affect all downstream scores appropriately
    # default is no scaling or shifting, don't do anything in that case (triggered if scale is false [including undefined and zero], and if it's 1; and if shift is false [undefined or zero])
    # note shifting is done first, so combinations of shift and weights act like "weight*(score-shift)"
    if ($shiftBits) {
	# unlike scaling, shifting can be either positive or negative
	$net = shiftEdges($net, $shiftBits);
    }
    if ($scaleBits && $scaleBits != 1) { # if a true value, and not 1
	# if user provides a string, that should be an error too, won't bother testing it but perl will complain above if warnings are on
	die "Error: context scale ($scaleBits) must be positive" if $scaleBits < 0; # this should be a fatal error
	$net = scaleEdges($net, $scaleBits);
    }
    
    # all done, return the scores net we want!!!
    return $net;
}

sub makeOldDpucCounts {
    # want to make a symmetric count network exactly like the old dPUC used to
    my ($net) = @_; # input is directed count network
    my %netOld; # here's the output we'll construct
    foreach my $acc1 (keys %$net) {
	my $acc2c = $net->{$acc1};
	while (my ($acc2, $c) = each %$acc2c) {
	    $netOld{$acc1}{$acc2} += $c; # add count in the direction it was observed
	    $netOld{$acc2}{$acc1} += $c if $acc1 ne $acc2; # if we're not in the diagonal, add it in the other direction as well!!!
	}
    }
    return \%netOld; # done, return this!!!
}

sub mapPfamNetAccs2ints {
    # this sub takes a standard net with two pfam accessions mapped to a score, and turns it into a parallel array, one with the accessions turned into integers and combined into a matrix flat number ij2m(i,j), and another array with the corresponding score
    # this format is optimal for fast queries in C code, but it's not very perlish at all
    # this code also filters to only keep positive scores
    my ($net, $acc2i) = @_;
    # these are the parallel arrays we want
    my @ts;
    my @tscores;
    while (my ($acc1, $acc2s) = each %$net) {
	next if $acc1 eq 'NEG'; # this special entry won't work here because "NEG" isn't an integer, things would die
	my $accInt1 = $acc2i->{$acc1}; # explicit map
	while (my ($acc2, $s) = each %$acc2s) {
	    next unless $s > 0; # here we implement the positive score filter for posElim downstream
	    my $accInt2 = $acc2i->{$acc2}; # explicit map
	    my $t = EncodeIntPair::ij2m($accInt1, $accInt2); # flatten with ij2m(), with i and j in the same order as in input
	    push @ts, $t; # add flattened index of Pfam accessions pair
	    push @tscores, $s; # add score in parallel
	}
    }
    # unfortunately, the above procedure doesn't produce a list of @ts that is sorted ascending (important for binary searches!), so we have to sort now, then return
    return paraSortNum(\@ts, \@tscores);
}

sub scaleEdgesAndInt {
    # this takes a net with floating point values and scales them by the provided $scoreScale, then turns into integers via rounding
    # net is modified by copy, original remains the same
    # works equally well with directional nets, or one or both ways
    my ($net, $scoreScale) = @_;
    my %net; # new network with scaled ints
    while (my ($acc1, $acc2s) = each %$net) {
	while (my ($acc2, $s) = each %$acc2s) {
	    $net{$acc1}{$acc2} = int $s*$scoreScale+1/2;
	}
    }
    return \%net; # done, return modified network
}

sub getNodes {
    # returns a hash ref that maps every node key to 1 if it was present in this network
    
    # this is a naive method that assumes the net is not both ways (otherwise the keys of the first level are all keys we want).
    my ($net) = @_;
    my %nodes; # this will be the (unique) set of acci's (acc indexes) with context
    while (my ($acc1, $acc2s) = each %$net) {
	$nodes{$acc1} = 1;
	foreach my $acc2 (keys %$acc2s) {
	    $nodes{$acc2} = 1;
	}
    }
    return \%nodes; # return the map we want!
}

sub paraSortNum {
    # this sorts multiple arrays in parallel by sorting the first one numerically ascending
    # adapted from http://www.perlmonks.org/?node_id=350180
    # NOTE: I could probably also try a schwartzian version, which should be easy to write for two arrays (not for n) and it might be faster than this, who knows
    # http://en.wikipedia.org/wiki/Schwartzian_transform
    my ($arrayMaster, @arrayRefs) = @_;
    # make the list of indeces
    my @indexes = (0 .. $#$arrayMaster);
    # this sorts the indexes so that the master array is sorted numerically ascending
    @indexes = sort { $arrayMaster->[$a] <=> $arrayMaster->[$b] } @indexes;
    # finally, this takes care of sorting all the arrays appropriately before returning
    map { [ @$_[ @indexes ] ] } ( $arrayMaster, @arrayRefs );
}

sub shiftEdges {
    # this takes a net with real values and shifts them by the provided $scoreShift
    # net is modified by copy, original remains the same
    # works equally well with directional nets, or one or both ways
    my ($net, $scoreShift) = @_;
    my %net; # new network with scaled ints
    while (my ($acc1, $acc2s) = each %$net) {
	while (my ($acc2, $s) = each %$acc2s) {
	    $net{$acc1}{$acc2} = $s + $scoreShift;
	}
    }
    return \%net; # done, return modified network
}

sub scaleEdges {
    # this takes a net with real values and scales them by the provided $scoreScale
    # net is modified by copy, original remains the same
    # works equally well with directional nets, or one or both ways
    my ($net, $scoreScale) = @_;
    my %net; # new network with scaled ints
    while (my ($acc1, $acc2s) = each %$net) {
	while (my ($acc2, $s) = each %$acc2s) {
	    $net{$acc1}{$acc2} = $s*$scoreScale;
	}
    }
    return \%net; # done, return modified network
}

sub sumEdges {
    # add the values of every edge in net
    # edges will be double counted if net is symmetric and edges are listed both ways, but should be ok one way or directional (will add both directions)
    my $net = shift;
    my $cTot = 0;
    # navigate net
    while (my ($acc1, $acc2c) = each %$net) {
	while (my ($acc2, $c) = each %$acc2c) {
	    $cTot += $c;
	}
    }
    return $cTot;
}

sub set2hashmap {
    # given a set, make a hash that maps each key to an index
    # input can also be an array ref (order is preserved in that case)
    # order is arbitrary here
    my ($set) = @_;
    my %hash;
    my @keys = 'HASH' eq ref $set ? keys %$set : @$set; # assume array ref if not hash ref
    my $l = @keys;
    for (my $i = 0; $i < $l; $i++) {
	$hash{$keys[$i]} = $i;
    }
    return \%hash;
}

sub log2 {
    return log($_[0])/log(2);
}

1;

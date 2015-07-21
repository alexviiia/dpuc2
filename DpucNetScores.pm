# Copyright (C) 2014 Alejandro Ochoa, Singh Research Group, Princeton University
# This file is part of dPUC.
# dPUC is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
# dPUC is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with dPUC.  If not, see <http://www.gnu.org/licenses/>.

package DpucNetScores;
our $VERSION = 1.01;

# 2015-06-05 10:52:14 EDT
# v1.01 - script now parses files that start with list of nodes
# will incrementaly adjust scripts to make proper use of this new info

use lib '.';
use FileGz;
use EncodeIntPair;
use strict;

# for Dpuc only
sub loadNet {
    # this is a wrapper around the file that normally loads the net, includes the typical pre-processing that is necessary for posElim or nonOvNeg
    my ($scoreScale, $accs, $fi, $panning, $alphaSelfToTrans, $scaleExp, $scaleContext, $shiftContext, $cCut, $w, $fi2, $cCut2, $gamma) = @_;
    
    # always load raw counts and compute scores on the fly!
    my ($net, $nodes) = loadContextCountsNet($fi, $cCut, $w, $fi2, $cCut2);
    my $n = @$nodes; # counts2scores needs to know how many families there are!

    # as soon as we can, let's compare the outside @$accs with @nodes hash/set, they should be same or die!
    # this way, we catch a version mismatch as early as possible, though this is technically only between Pfam-A.hmm.dat and the dpucNet (input domain preds could be yet of a different version, so we'll be checking them again later as well.
    die "Fatal: domain family sets disagree between input Pfam-A.hmm.dat and this dpucNet: $fi\nTheir Pfam versions probably disagree!\nNo output was generated\n" unless setIdentity($accs, $nodes);
    
    # construct scores given counts and other parameters
    $net = counts2scores($net, $n, $panning, $scaleExp, $alphaSelfToTrans, $gamma, $scaleContext, $shiftContext);
    
    # need to make sure this exists before processing nets (necessary for C code and to define $n for on-the-fly scores)
    my $acc2i = set2hashmap($nodes); # now make hash that maps each accession to its index, store as global!
    
    # when doing posElim, scale scores to ints, turn into parallel arrays, store in C space (so fast code access it fast)
    
    # scale to have integers, trying to speed up calculations
    my $netScaleInt = scaleEdgesAndInt($net, $scoreScale);
    # this turns double hash net into a parallel array that is a lot more C-friendly, but not perlish at all
    # additionally, this function removes non-positive scores
    my ($netTs, $netTscores) = mapPfamNetAccs2ints($netScaleInt, $acc2i);
    
    # posElim needs a list of domains that carry context information, also gets stored in C space
    # however, all domains have context (self context in particular) if the "countless self" score is positive
    my $negSelfInt = $netScaleInt->{PFNEG}{SELF}; # get scaled version!
    # acc2cc is an hash ref that maps each acc to 1 if it is context carrying, undef otherwise
    # note acc2cc contains PFNEG, SELF, and TRANS too, doesn't affect out analysis at all!
    # if negative self contex is actually positive, make a map that makes every acc context-carrying
    my $acc2cc = $negSelfInt > 0 ? {map { $_ => 1 } @$accs} : getNodes($net);
    
    # only nonOvNeg uses $net, return it (safe to ignore for posElim), also return the stuff needed for C structures (for posElim)
    return ($net, $acc2i, $netTs, $netTscores, $negSelfInt, $acc2cc); # stuff that Dpuc2 needs, many for passing to C
}

# Dpuc and CODD
sub loadContextCountsNet {
    # loads and process counts network (does not apply scores).  This applies to CODD as well (which only uses the topology)
    # assumes fi is defined, cCut doesn't need to be (means no filter)
    my ($fi, $cCut, $w, $fi2, $cCut2) = @_;
    
    my ($net1, $nodes1) = parseNet($fi, $cCut); # loads net and applies counts threshold
    
    if ($w) { 
	# blend with second net
	# assume parameters are the same as for first if undefined
	$fi2 = $fi unless $fi2;
	$cCut2 = $cCut unless defined $cCut2; # but zero shouldn't be overwritten!!!  only undefined!
	my ($net2, $nodes2) = parseNet($fi2, $cCut2);
	# to blend, node lists should match (otherwise there's a version mismatch)
	die "Error: network files have different node lists: $fi, $fi2\n" unless setIdentity($nodes1, $nodes2);
	# blend nets
	$net1 = blendNets([$net1, $net2], [$w, 1-$w]); # apply weight w to net1, (1-w) to net2
    }
    
    return ($net1, $nodes1);
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

sub counts2scores {
    # copied and simplified from dpucNet04.counts2scores.pl, in order to generate one single net on the fly from a single set of counts
    # no combinatorics here like original script (which was left unchanged), and we dropped some useless features (notably clans map)
    # dropped intermediate panned counts net creation, made sense for combinatorics but not here!
    # on the other hand, made algorithms faster, which matters a lot more when doing things on the fly!
    my ($acc2acc2c, $n, $panning, $scaleExp, $alphaSelfToTrans, $gamma, $scaleBits, $shiftBits) = @_; # input count structure and params
    $acc2acc2c = expEdges($acc2acc2c, $gamma) if $gamma && $gamma != 1; # exponentiate edges if necessary (a gamma of zero/false wouldn't be allowed, but it can be used conveniently to mean "don't exponentiate".
    # get sum of counts from network, parameter independent!
    my $cTot = sumEdges($acc2acc2c);
    
    # part 1: panning counts
    my $pan = $panning eq 'inf' ? 0 : (2**-$panning)/2; # get the regular parameter, without the log2 and also divided by two (so max panning is a factor of 1/2, NOT 1)
    my $pan2 = 1-$pan; # the complement of the panning, used very often too
    
    # part 2: score params
    # deduce priors from their ratio and total sum... (math leading to these solutions not shown, but it's really easy!)
    $alphaSelfToTrans = 2**$alphaSelfToTrans; # undo log2 scale!
    my $scale = 2**-$scaleExp; # get actual scale from specified exponent
    my $alphaTrans = $scale*$cTot/(($alphaSelfToTrans + $n-1)*$n);
    my $alphaSelf = $alphaTrans*$alphaSelfToTrans;
    # super denominator (pre-calculated for efficiency) includes
    # - denominator that completes the pij probability (cTot + [forced alphaSum = scale*cTot])
    # - pi and pj backgrounds, which are just 1/n each
    my $superDenom = $cTot*(1+$scale)/($n*$n);
    
    my %acc2acc2s; # scores we want!
    # navigate net
    # pairs observed both ways will be processed twice, but meh, it'll just be overwritten... and it doesn't happen that often to try to prevent it (we'll be checking more often than it's worth).
    while (my ($acc1, $acc2c) = each %$acc2acc2c) {
	while (my ($acc2, $c) = each %$acc2c) {
	    my $isDiag = $acc1 eq $acc2; # remember if we're in diagonal or not
	    
	    # part 1: panning
	    my $cp = $c; # panned count, default unpanned
	    my $cr = 0; # reverse panned count, default zero
	    # in two trivial cases counts don't change (no global panning, or we're in diagonal, which is fixed to panning).  This shouldn't affect reverse count at all, don't change anything there either.  Don't run through code below because I'm afraid math will introduce precision errors.  And this is faster!
	    if ($pan && !$isDiag) { # do a bit more work
		# get reverse count too, or undefined if not available
		my $c2 = $acc2acc2c->{$acc2}{$acc1};
		$c2 = 0 unless defined $c2; # to simplify code, set to zero when unavailable
		# perform panning!  Overwrite previous defaults
		$cp = $pan2*$c + $pan*$c2; # original direction
		$cr = $pan2*$c2 + $pan*$c; # other direction.  Note since $pan>0 here, this is always non-trivial! (because $c is never zero in this case, and neither is $pan2 for allowed inputs, it's always in [1/2,1])
	    }
	    
	    # part 2: scores, store in net
	    $acc2acc2s{$acc1}{$acc2} = log2(($cp + ($isDiag ? $alphaSelf : $alphaTrans))/$superDenom); # regular score always gets saved
	    $acc2acc2s{$acc2}{$acc1} = log2(($cr + ($isDiag ? $alphaSelf : $alphaTrans))/$superDenom) if $cr; # reverse score only if non-trivial (which correlates with $cr != 0 because of "if" above)
	}
    }
    
    # add entries for the uniform negative scores (which are now different between self and trans! so two scores!)
    $acc2acc2s{PFNEG} = {
	'SELF'  => log2($alphaSelf /$superDenom),
	'TRANS' => log2($alphaTrans/$superDenom),
    };
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

sub mapPfamNetAccs2ints {
    # this sub takes a standard net with two pfam accessions mapped to a score, and turns it into a parallel array, one with the accessions turned into integers and combined into a matrix flat number ij2m(i,j), and another array with the corresponding score
    # this format is optimal for fast queries in C code, but it's not very perlish at all
    # this code also filters to only keep positive scores
    my ($net, $acc2i) = @_;
    # these are the parallel arrays we want
    my @ts;
    my @tscores;
    while (my ($acc1, $acc2s) = each %$net) {
	next if $acc1 eq 'PFNEG'; # new entry that only occurs in Dpuc 2.0, obviously wouldn't work because "NEG" isn't an integer, things would act funny if not die
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

sub expEdges {
    # given a network with raw counts, exponentiates each value by the given exponent $gamma
    # has the effect of increasing the frequency of the heaviest edges relative to the weakest ones, which leads to smaller scores for the weakest edges, but notably without setting any of them to the worst negative score reserved for unobserved edges.  Here weakest edges (after log-odds, of course) will be smoothly downweighed to have near zero scores, or slightly negative scores, as desired by adjusting gamma
    my ($net, $gamma) = @_;
    my %net; # new net, so don't edit input structure
    # navigate net
    while (my ($acc1, $acc2c) = each %$net) {
	while (my ($acc2, $c) = each %$acc2c) {
	    $net{$acc1}{$acc2} = $c**$gamma; # here we exponentiate the counts, then place in new network structure
	}
    }
    return \%net; # return new net
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

sub blendNets {
    my ($nets, $ws) = @_;
    # given a list of nets with parallel blending weight vector, will produce a normalized network where each edge is the weighted average of the input edges, after each being normalized
    # really only want pair case, but it's harder to not generalize!
    # weights must be between 0 and 1, won't verify that this is the case, but you might get unexpected results otherwise.
    # preserves directionality of input nets!
    # turns out this stupid loop is best way I could think of producing new net, since edges might be missing in any of the inputs
    my %net; # this is the new network that is the blend of the inputs
    my $n = @$nets;
    for (my $i = 0; $i < $n; $i++) {
	my $net = $nets->[$i];
	# first get normalizing constant, incorporate into weight
	my $w = $ws->[$i]/sumEdges($net);
	# navigate net
	while (my ($acc1, $acc2c) = each %$net) {
	    while (my ($acc2, $c) = each %$acc2c) {
		$net{$acc1}{$acc2} += $c*$w; # add normalized and weighted version of this edge to final blended net
	    }
	}
    }
    return \%net; # done, return desired net
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

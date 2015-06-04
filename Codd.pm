# Copyright (C) 2014 Alejandro Ochoa, Singh Research Group, Princeton University
# This file is part of dPUC.
# dPUC is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
# dPUC is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with dPUC.  If not, see <http://www.gnu.org/licenses/>.

package Codd;
our $VERSION = 1.00;
use lib '.';
use Domains;
use strict;

# implements the original CODD (for Pfam GA only), plus some extensions inspired by dPUC (directionality and negative context, in particular)

# this is our implementation of the CODD algorithm (Terrapon et al 2009), with an added optional negative filter to try to match our dPUC algorithms!  Also with an additional context directionality constraint that the original CODD lacked.
# importantly, this function should produce the same output as the old codd() when a symmetrized $net is input
# key is for sorting/ranking at overlap removal steps, all parts default to sorting by E if key is undefined.  But it's useful to use q-values, for example

# we have a few differences of necessity, including treating every domain separately instead of clustering by families, since our benchmarks can only handle the former.  Their lack of "repeat" recovery might show using this procedure.
# they also have versions where "potential" domains can rescue each other, but using different thresholds than using GA only as validating domains.  Lastly, they can use Interpro to rescule potential domains, but that method wouldn't compete fairly with ours and I don't think it would be that much better.

sub codd {
    # the hits of a single protein as input, along with two nets needed for processing
    my ($hits, $net, $nestingNet, $doNeg, $fCut, $lCut, $key, $doGaFree) = @_;
    $key = 'E' unless defined $key;
    
    if ($doGaFree) {
	# this second version is used by DPUC2 as a prefilter, to have a smaller set of domains without overlaps as input
	# since algorithm is much simplified, let's skip some of the unnecessary steps (compared to the full CODD that uses full GA)
	
	# coddSeparateGa not needed because it's trivial (it's the top ranking domain in input, which is selected inside coddAddNonGa instead
	# coddNegCleanGa not needed because there's only one "GA" domain (the top-ranking one)
	# coddAddNonGa (below) is the workhorse for this version
	# note that no GA hits are passed, all hits are passed as $hitsNonGa (with doGaFree=1, it'll be ok)
	return coddAddNonGa(undef, $hits, $net, $nestingNet, $doNeg, $fCut, $lCut, $key, $doGaFree);
    }
    else {
	# the standard GA-based CODD follows... with our new bells and whistles
	
	# separate GA from non GA, remove overlaps in GA too
	my ($hitsGa, $hitsNonGa) = coddSeparateGa($hits, $nestingNet, $fCut, $lCut, $key);
	return [] unless scalar @$hitsGa; # if we had no GA domains, this method doesn't return anything... bail out early!
	
	# now that we have the GA domains, we need to make sure there are no negative interactions within this set.  This is the first (and maybe the largest) step towards improving our signal-to-noise ratio
	# not done if we aren't using negative context
	$hitsGa = coddNegCleanGa($hitsGa, $net, $key) if $doNeg;
	
	# a similar sub will add the non-GA hits now, and sort by start!  This is the final output
	return coddAddNonGa($hitsGa, $hitsNonGa, $net, $nestingNet, $doNeg, $fCut, $lCut, $key);
    }
}

### INTERNAL FUNCTIONS that implement modular steps in CODD extensions

sub coddSeparateGa {
    # all the flavors of CODD use this function
    # besides separating GA from non-GA, we also remove disallowed overlaps in the GA list (should typically remove overlaps within the same clan, but it's safer to remove all overlaps that Pfam hasn't observed before)
    my ($hits, $nestingNet, $fCut, $lCut, $key) = @_;
    # first thing we do is separate GA from non GA, using the flag we always save per hit in the initial file
    my @hitsGa = ();
    my @hitsNonGa = ();
    foreach my $h (@$hits) {
	if ($h->{GA}) { push @hitsGa, $h; }
	else { push @hitsNonGa, $h; }
    }
    # The GA domains we use are supposed to be the final GA, excluding disallowed overlaps, so let's process like that...
    my $hitsGa = Domains::removeOverlapsByScore(\@hitsGa, $nestingNet, $fCut, $lCut, $key); # remove overlaps as desired, also sorts hits by start
    return ($hitsGa, \@hitsNonGa);
}

sub coddNegCleanGa {
    # this sub takes in a list of GA domains, sorts them by E-value, goes through each one, and removes each if it has any negative context with the previous domains that passed.  This way the output list has no negative context domains.
    # We don't have to worry about overlaps here, that has to be taken care of before, outside this sub.  All remaining overlaps here are allowed.
    my ($hitsGa, $net, $key) = @_;
    my @hitsGa; # this is the new list of GA domains with positive context only.
    # sort input hits by score
    $hitsGa = Domains::sortByCol($hitsGa, $key); # sort by "normalized score"
    # browse input domains
    foreach my $h (@$hitsGa) {
	my $acc = $h->{acc}; # copy variable we use often
	my $haveNegContext = 0; # this boolean reminds us us if we had negative context or not
	# browse the domains that have already passed
	foreach my $h2 (@hitsGa) {
	    # decide if we have positive context (absence of positive context is negative context!), which requires the two domains to be properly ordered
	    # the order is determined by the start sites, which has to be different since we had no overlaps getting here!
	    unless ($h->{start} < $h2->{start} ? $net->{$acc}{$h2->{acc}} : $net->{$h2->{acc}}{$acc}) {
		$haveNegContext = 1; last; # this pair had negative context, mark and stop looking
	    }
	}
	push @hitsGa, $h unless $haveNegContext; # add to list unless it had negative context with anything
    }
    # done! return the GA hits that passed!
    return \@hitsGa;
}

sub coddAddNonGa {
    # this sub takes in a list of GA domains, and non-GA domains, sorts the last by E-value, goes through each one, passes first filter if it doesn't overlap with domains that have already passed, second filter if there is positive context with anything, and third filter if it doesn't have any negative context with the previous domains that passed (if requested with $doNeg).  This way the output list has all positive and no negative context domains.
    # We do have to worry about overlaps here, for the non-GA domains
    my ($hitsGa, $hitsNonGa, $net, $nestingNet, $doNeg, $fCut, $lCut, $key, $doGaFree) = @_;
    my $haveNestingNet = defined $nestingNet && scalar keys %$nestingNet; # a boolean we use...
    # sort input non-GA hits by score.  This is usually a very long list
    $hitsNonGa = Domains::sortByCol($hitsNonGa, $key); # sort by "normalized score"
    # this is the new list of domains that pass, initialize it by copying GA domains
    # but note that the "doGaFree" mode has no $hitsGa (since they're only used in the next line, it's ok, they can even be undefined!), instead the first domain that apsses is the top-ranking domain in $hitsNonGa (which in this mode is all domains)!  Also note, this domain is removed from $hitsNonGa so we don't go through it again!
    my @hitsPass = $doGaFree ? (shift @$hitsNonGa) : @$hitsGa;
    # browse input non-GA domains
    foreach my $h (@$hitsNonGa) {
	# copy values we call repeatedly
	my ($s, $e, $acc) = @{$h}{qw(start end acc)};
	# there are really three filters: no overlaps, at least one positive context pair, and no negative pairs.  We can do all tests in the same loop!  (it might look uglier on the code, but I think it makes more sense, and is it faster?)
	# booleans that keep track of things
	my $haveOverlaps = 0; # Filter 1: keep track of overlaps (all overlaps are bad, either with GA or non-GA that have passed)
	my $havePosScore = 0; # Filter 2: have at least one positive score (with GA domains only!  This is how CODD does it, and might be important to a lower FP rate.  Unless we've turned on the flag that allows all passing domains to be validating domains)
	my $haveNegScore = 0; # Filter 3: have at least one negative score (all negative scores are bad, either with GA or non-GA that have passed)
	# navigate domains that have passed to find conflicts (filters 1,3) or support (filter 2)
	foreach my $h2 (@hitsPass) {
	    my ($s2, $e2, $acc2) = @{$h2}{qw(start end acc)};
	    # Filter 1: check for overlaps first
	    # $h overlaps with a region that has already passed
	    if ( Domains::overlap($s, $e, $s2, $e2, $fCut, $lCut) ) {
		# overlap isn't allowed
		# if we don't have a nesting net, we won't bother looking for accessions
		unless ( $haveNestingNet && $nestingNet->{$acc}{$acc2} ) {
		    $haveOverlaps = 1; last;
		}
	    }
	    # if we continue, then there was no disallowed overlap
	    # decide if we have positive context, which requires the two domains to be properly ordered
	    # the order is determined by the start sites, which has to be different since we had no overlaps getting here!
	    # (with permissive overlaps it's a bit more subtle, but two regions with the same start are necessarily nested, which is always a disallowed overlap, so even in that case the below code/reasoning works)
	    if ($s < $s2 ? $net->{$acc}{$acc2} : $net->{$acc2}{$acc}) {
		# this pair has positive context (will only count if the "validating" domain was indeed GA, unless we've turned on the flag that allows all passing domains to be validating domains), mark and move on.  We need to keep checking because overlaps might cancel this positive contribution.
		$havePosScore = 1 if $h2->{GA};
	    }
	    elsif ($doNeg) { $haveNegScore = 1; last; } # this pair had negative context, mark and stop looking (but only if $doNeg, we definitely want to continue otherwise, no "last")
	}
	# we have seen all relevant pairs, decide whether this domain passes or not
	push @hitsPass, $h if $havePosScore && !$haveOverlaps && !($doNeg && $haveNegScore); # add to list if we had at least one positive context edge, no overlaps, and no negative context with anything (if requested via $doNeg)
    }
    # done! return our new list of GA domains plus the non-GA that passed the filters!
    # also sort by start because this is the final output
    return Domains::sortByCol(\@hitsPass, 'start');
}



1;

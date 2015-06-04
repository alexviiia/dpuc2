# Copyright (C) 2014 Alejandro Ochoa, Singh Research Group, Princeton University
# This file is part of DomStratStats and dPUC.
# DomStratStats and dPUC are free software: you can redistribute them and/or modify them under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
# DomStratStats and dPUC are distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with DomStratStats and dPUC.  If not, see <http://www.gnu.org/licenses/>.

package Domains;
our $VERSION = 1.00;
use strict;

# general functions for my domain structure
# I frequently use "Hits" instead of "Domains", just some old nomenclature from who-knows-where (Pfam or BLAST?)

sub sortByCol {
    # a general method for sorting by the value of a column $col
    # always assumes we're scoring numbers and that we want ascending (i.e. best E-values, p-values, and q-values are smaller)
    my ($hits, $col) = @_;
    return [] unless $hits && @$hits; # bail out early if hits are undefined (return defined empty array ref), so programs don't complain with missing entries
    
    # let's check for undefined values, these are always a problem...
    # finding them should aways be fatal, we're doing something wrong if we're trying to sort undefined values
    foreach my $h (@$hits) {
	die "Fatal: in Domains::sortByCol, wanted to sort undefined value in column $col\n" unless defined $h->{$col};
    }
    
    # simply sort and return a new array reference (original array reference is not modified)
    return [sort { $a->{$col} <=> $b->{$col} } @$hits];
}

sub removeOverlapsByListOrder {
    # this sub takes a list of hits and removes all overlaps, keeping whatever came first in the list
    # the idea is that the hits have been sorted outside this sub with a particular priority in mind
    # we allow some overlaps via a nesting net, which maps accessions to the boolean of interest
    my ($hits, $nestingNet, $fCut, $lCut) = @_;
    return [] unless defined $hits && @$hits; # return an empty list if $hits is undefined or an empty list
    my $haveNestingNet = defined $nestingNet && keys %$nestingNet; # a boolean we use...
    my @hitsNew; # hits that pass
    foreach my $h (@$hits) {
	my ($s, $e, $acc) = @{$h}{qw(start end acc)}; # get coordinates and acc
	my $isOverlapping = 0;
	foreach my $h2 (@hitsNew) {
	    # $h overlaps with a region that has already passed
	    if ( overlap($s, $e, $h2->{start}, $h2->{end}, $fCut, $lCut) ) {
		# overlap isn't allowed
		# if we don't have a nesting net, we won't bother looking for accessions
		unless ( $haveNestingNet && $nestingNet->{$acc}{$h2->{acc}} ) {
		    $isOverlapping = 1; last;
		}
	    }
	}
	# will only pass if it didn't overlap with another domain that ranked higher
	push @hitsNew, $h unless $isOverlapping;
    }
    return \@hitsNew;
}

sub removeOverlapsByScore {
    # I tried to be more general than Pfam.  This removes all overlaps, not only the ones of the same accession, unless they're allowed as defined by a hash.  The priority is set by the E-value (essentially a normalized score).  This should be useful for Supfam and Smart processing as well!!!  Again, this is more general than that processing, since it assumes no GA filtering (can do with predictions in general), and the nestingNet allows special overlaps (a feature that only Pfam has, in a way, but may be defined for Supfam as well)
    my ($hits, $nestingNet, $fCut, $lCut, $key) = @_;
    $key = 'E' unless defined $key; # E-values are maintained as default for backward-compatibility, but can set thresholds on q instead, or other options.
    $hits = sortByCol($hits, $key); # sort by "normalized score"
    $hits = removeOverlapsByListOrder($hits, $nestingNet, $fCut, $lCut); # remove overlaps by score
    return sortByCol($hits, 'start'); # re-sort by start, return
}

sub overlap {
    # determines whether domains overlap using thresholds on the allowed overlap size in terms of both amino acids and fraction relative to the smallest domain
    # if f cut is false the function behaves like old overlap, and is nearly as fast too!  So this can behave as a drop-in substitute for the less general function (when called with old params and no additional thresholds)!
    # tried to optimize relative to the version that always returns l,f and doesn't set thresholds.  Figured can take shortcuts if we don't care about the exact values.
    # comparison is l<=lCut and f<=fCut are good, unlike Orengo's (l off by 1 cause it's an int, f is practically the same cause it's continuous)
    # for efficiency, we assume nesting is never acceptable.  (except in known exceptions, i.e. nesting net, but that's handled outside this function)  This is really ok since goal is to emulate Orengo's method, and it doesn't make sense biologically otherwise anyway.
    # for convenience, setting lCut to false (zero or undef) should mean we don't want to use it (behaves like lCut=inf).  Idea is, if we don't want to use fractional threshold, we can set that to 1 (the logical max, although we always mark that case as overlap for code simplicity, read "nesting" above), but there isn't a logical max for length in amino acids.  On the other hand, if we don't want to tolerate any amount of overlap (like the old, strict overlap tester), we can set fCut=0 and lCut to anything.  This way we can be absolutely strick and also absolutely permissive for both kinds of thresholds
    my ($s1, $e1, $s2, $e2, $fCut, $lCut) = @_; # start1, end1, start2, end2
    ($s1, $e1, $s2, $e2) = ($s2, $e2, $s1, $e1) if $s1 > $s2; # force regions to be ordered so that $s1 <= $s2
    if ($s2 < $e1) { # not <= cause ends are exclusive (a la perl, but not like hmmer).
	return 1 unless $fCut; # fCut==0 (or just false in general) is strict case, return quickly if that's the case.
	# there are really only two cases for overlap... e1<e2 means hits aren't nested, otherwise second one (in new order) is nested in first (and quantities to threshold over are more trivial)
	if ($e1 < $e2) {
	    my $l = $e1-$s2; # in this case (not nesting) this is the size of the overlap in amino acids
	    # get lengths, unfortunately we need to bother with this here, but maybe we can define vars with global scope to save a bit of time there
	    my $l1 = $e1-$s1;
	    my $l2 = $e2-$s2;
	    my $lRef = $l1 < $l2 ? $l1 : $l2; # this saves the min
	    # set threshold and return
	    # in first case, have an overlap regardless of the $lCut
	    # in second case set threshold only if lCut is true ( first special case lCut==0 (or false in general) behaves like lCut==inf... decide only using fCut!)
	    return ($l/$lRef > $fCut || ($lCut && $l > $lCut) );
	}
	else { return 1; } # if e1 >= e2, f==1 guaranteed, and that's never good biologically, so let's stop it there.  If we wanted to actually apply lCut, we'd have to figure out the length of the smallest region, but again let's not waste time doing that...
    }
    else { return 0; } # no overlap
}

sub filterPfamClans {
    # for removing overlapping hits of the same clan, keeps only the one with the best E-value
    # can also be used for Superfamily to remove overlapping hits of the same fold, provided the right mapping is passed in $acc2clan
    
    # remember not to use this on Superfamily hits, only PFAM
    # this script will not check that all entries are valid PFAM, that would make it slower (not sure if significantly so though)
    my ($hits, $acc2clan, $fCut, $lCut) = @_;
    # sort hits by E-value
    my $hitsSorted = sortByCol($hits, 'E');
    my @hitsNew;
    foreach my $h (@$hitsSorted) {
	my $acc = $h->{acc};
	my $clan = $acc2clan->{$acc};
	unless (defined $clan) { push @hitsNew, $h; next; } # hit automatically passes if it doesn not belong to a clan
	# if passes, $h belongs to a clan
	my ($s1, $e1) = @{$h}{qw(start end)}; # get hit1's range with hash ref slice
	my $isOvC = 0;
	foreach my $h2 (@hitsNew) {
	    my $acc2 = $h2->{acc};
	    my $clan2 = $acc2clan->{$acc2};
	    next unless defined $clan2; # proceed to comparison if this also belongs to a clan
	    if ( overlap($s1, $e1, $h2->{start}, $h2->{end}, $fCut, $lCut) ) { # $h overlaps with a region that has already passed
		if ($clan eq $clan2) { $isOvC = 1; last; } # if clans are also the same, toss this $h
	    }
	}
	push @hitsNew, $h unless $isOvC; # will only pass if it didn't overlap with another domain of the same clan that ranked higher
    }
    return sortByCol(\@hitsNew, 'start'); # sort for print now, in case this is the end of the line
}

sub removeOverlapsByAccPopularity {
    # NOTE: maybe a better way of doing this is using sequence scores!
    # idea is to first count number of domains per family in GA, to get a rough idea of what families have the most "trustworthy" repeats
    # in the second round we look at the raw data (not just GA) and sort by these counts (breaking ties by E-value as is traditional), then we remove overlaps by list order
    # the output solutions should be more homogeneous than regular GA, but we could be introducing noise too
    # it's like assigning highest context to self, but we're not performing dPUC directly, just kinda guessing what should have the best context
    # this guess is then used by dPUC when we don't reach an optimal solution using lp_solve (so it's just a backup plan, not the main plan)
    # it should be fed back into lp_solve in case negative context remains
    # this solution is also compared to the dPUC solution using regular removeOverlapsByScore(), which is another sensible simplified solution
    my ($hits, $nestingNet, $fCut, $lCut, $key) = @_;
    $key = 'E' unless defined $key; # E-values are maintained as default for backward-compatibility, but can set thresholds on q instead, or other options.
    # get GA approximation from raw input hits
    # I say approximation because we don't use "clans" as we should, we remove all overlaps in GA instead except when explicitly allowed
    # it should be the same in theory, but the structure is much simpler than having to keep track of clans
    my $ga = filterBooleanCol($hits, 'GA'); # this applies the GA thresholds we've pre-calculated
    $ga = removeOverlapsByScore($ga, $nestingNet, $fCut, $lCut, $key);
    # now count accs in GA
    my %acc2c;
    foreach my $h (@$ga) { $acc2c{$h->{acc}}++; }
    # make sure non-GA accs get zeros as counts (to keep perl from complaining later at "sort")
    foreach my $h (@$hits) { $acc2c{$h->{acc}} = 0 unless $acc2c{$h->{acc}}; }
    # now that we have these counts, use them to rank input raw domains
    # acc counts are descending, but $key is ascending
    @$hits = sort { $acc2c{$b->{acc}} <=> $acc2c{$a->{acc}} || $a->{$key} <=> $b->{$key} } @$hits;
    # the rest is easy given these functions...
    $hits = removeOverlapsByListOrder($hits, $nestingNet, $fCut, $lCut); # remove overlaps by score
    return sortByCol($hits, 'start'); # re-sort by start, return
}

sub overlapNetBinarySearch {
    # this calculates an overlap networt (both ways) from scratch, doesn't assume anything about the starting data, and shouldn't alter it either
    # if available, uses an "allowed" net (like the Pfam nesting net) to omit overlaps from net.  Right now the net has to be exactly the Pfam nesting net because it uses accessions rather than hit references to do the mapping!!!  It's not worth doing the mapping outside in a different structure, the map is a one-time thing.  Keep in mind for future applications though!
    # uses the binary search subs to make this very fast
    my ($hits, $nestingNet, $fCut, $lCut) = @_;
    
    # at least once a protein sent an empty hits array.  Let's just return an empty net in that case (it should not be used anyway)
    my $numHits = scalar @$hits;
    return {} unless $numHits;
    
    my %overlapNet; # the structure we want, for searching in dPUC (but hit refs are stringified in keys)
    my @pairs; # a different structure, for looping rather than searching (and hit refs are preserved, not stringified)
    my $haveNestingNet = defined $nestingNet; # a boolean that tells us quickly if we have a net or not
    
    # initialize the data (sorting and lowStart mapping)
    # we overwrite the local $hits reference with an array ref that is sorted differently, but the original input array (which you see outside this function) is unchanged!
    ($hits, my $hit2lowStart) = overlapBinarySearchInit($hits);
    
    # now start searching, ascending in array, and avoid rechecking old overlaps by adjusting the 'lower' range
    my $upper = $numHits-1; # precalculate the upper range, should speed things up a bit
    for (my $i = 0; $i < $numHits-1; $i++) { # note the last hit is never the query, since it's always in other query's results
	my $hi = $hits->[$i]; # the query
	my $lower = $i+1; # this skips the query and the things before it
	# these are the things that overlap our query
	my $hitsOv = overlapBinarySearch($hi, $hits, $hit2lowStart, $lower, $upper, $fCut, $lCut);
	# now we need to add the results to the net, both ways
	my $nestingNetHere = $nestingNet->{$hi->{acc}} if $haveNestingNet; # copy down structure to reduce hash calls
	foreach my $hj (@$hitsOv) {
	    # skip if overlap is allowed
	    next if $haveNestingNet && $nestingNetHere->{$hj->{acc}};
	    $overlapNet{$hi}{$hj} = 1;
	    $overlapNet{$hj}{$hi} = 1;
	    push @pairs, [$hi, $hj]; # add pair that overlaps to list.  Note this is actual reference, not string version of reference
	}
    }
    
    # done, return the network!
    return (\%overlapNet, \@pairs);
}

sub overlapBinarySearchInit {
    # sorts the list the way we need it for queries, and also gets the hit2lowStart auxiliary array.  Doesn't alter the original array
    my ($hits, $s) = @_;
    $s = '' unless defined $s; # suffix, to use startAln and endAln instead ($s='Aln' in that case)
    
    # at least once, it happened that a protein sent an empty hits array, but the initial $lowStart chokes on that!  Let's just return empty arrays instead of attempting to make things work
    my $numHits = scalar @$hits;
    return ([], []) unless $numHits;
    
    # the very first thing to do is to sort the hits by end first, then by start to break ties
    my @hits = sort { $a->{'end'.$s} <=> $b->{'end'.$s} || $a->{'start'.$s} <=> $b->{'start'.$s} } @$hits;
    
    # now we create the array that tells us what is the lowest start position remaining relative to the current position
    my @hit2lowStart;
    my $lowStart = $hits[$numHits-1]{'start'.$s}; # the last entry initializes this variable
    # navigate it backwards
    for (my $i = $numHits-1; $i >= 0; $i--) {
	my $starti = $hits[$i]{'start'.$s}; # get this start
	$lowStart = $starti if $lowStart > $starti; # update the current lowest start
	$hit2lowStart[$i] = $lowStart;
    }
    
    return (\@hits, \@hit2lowStart);
}

sub overlapBinarySearch {
    # this is a very efficient search for overlaps, using a binary-search strategy
    # it assumes the query list is already sorted by end, then breaking ties by start (both ascending)
    # we also definitely need an array that for each position marks the lowest remaining start position
    # lower and upper are not required, but can be set to reduce searching ranges (if we've traversed the array up to down, for example).  They are in terms of the array indeces, not in terms of start or end coordinates!
    
    # binary search code based on General::binarySearch, itself based on/stolen from code: http://staff.washington.edu/jon/dsa-perl/bsearch
    # the rest of the algorithm was informed by http://stackoverflow.com/questions/303591/a-range-intersection-algorithm-better-than-on
    # particularly, this entry:
    #   1. Use binary search to find the first range with an End value of >= TestRange.Start
    #   2. Iterator starting at the binary search until you find an Start > TestRange.End:
    #      2a. If the range if the current range is within the TestRange, add it to your result.
    # I realized I needed to edit the structure to add one more array, one that maps each position to the lowest start position remaining (can be built in O(n), not bad compared to sorting)
    
    my ($hitQuery, $hits, $hit2lowStart, $lower, $upperIni, $fCut, $lCut, $s) = @_;
    $s = '' unless defined $s; # suffix, to use startAln and endAln instead ($s='Aln' in that case)
    # the query hit itself matters very little compared to just knowing the start and end, so we extract them now
    my $startQuery = $hitQuery->{'start'.$s};
    my $endQuery = $hitQuery->{'end'.$s};
    my $l1 = $endQuery-$startQuery; # might need it for permissive overlap thresholds.  Doesn't hurt to have it here even if we don't use it, it's a lot cheaper than the binary search
    
    # set these two default values
    $lower = 0 unless defined $lower;
    $upperIni = scalar(@$hits)-1 unless defined $upperIni;
    
    # @$hits are sorted by end first, so it's easy to check if any domains satisfy end > query.start by looking at the last one
    # if this isn't the case, then there's no matches!
    # use upperIni for this, most useful if that's defined
    return [] if $hits->[$upperIni]{'end'.$s} <= $startQuery;
    
    # from now on there's at least one domain with end > query.start, let's find the least such domain!
    # the first part is using a binary search to find the first range with end > query.start (this accounts for end-exclusive, the way we have it)
    my $upper = $upperIni; # we'll remember what the initial upper value was, since $upper gets updated as the search progresses
    my $i; # index of probe
    while ($lower < $upper) {
	$i = int(($lower + $upper)/2); # choose the index right at the middle of our range
	my $end = $hits->[$i]{'end'.$s}; # get the current end (remember this is end-exclusive)
	if ($end > $startQuery) { $upper = $i; } # this means we might have it, but we might need to look a bit lower
	else { $lower = $i+1; } # we definitely need to look higher.  The +1 guarantees that rounding down won't get us stuck in infinite iterations, right?
    }
    # since we're not looking for an exact number (many ranges can have the same end point, and we could have "hits" with identical ranges), we're only done when the range has been reduced to a point
    if ($lower == $upper) { $i = $lower; } # this is our starting point
    elsif ($upper == $upperIni && $lower == $upper+1) { return []; } # this is an anomalous case, happens when our query's start was higher than any of the ends present, so no hits were found!
    else { die "Unexpected error: lower $lower != upper $upper (upperIni $upperIni)!"; } # I don't think this ever happens
    
    # the second part is to iterate up until we find a hit with lowStart >= query.end (accounts for end-exclusive setup again)
    # we need to *check* and collect hits until this happens
    my @hitsOv;
    # we get the first start to start the loop (this also verifies that we got anything at all, failure to satisfy the condition means nothing overlapped our hit)
    my $lowStart = $hit2lowStart->[$i];
    while ($lowStart < $endQuery) {
	# this hit might be promising, but we need to check it too
	my $hi = $hits->[$i];
	my ($s2,$e2) = @{$hi}{'start'.$s, 'end'.$s};
	if ($s2 < $endQuery) { # this hit was good, in strict overlap sense
	    # but the story isn't over if we're using permissive overlap thresholds...
	    # we reimplement overlap here with specific shortcuts so we save a bit more time that we might need
	    if ($fCut) { # overlaps are strict if fCut == 0 (or any other false value)
		# will need to compute smallest of both ranges
		my $l2 = $e2 - $s2;
		my $lRef = $l1 < $l2 ? $l1 : $l2; # this saves the min
		# after this point, the identities of h1,h2 shouldn't matter, we can overwrite s2,e2 cause we don't use them anymore, but we shouldn't overwrite the query's range!
		my ($s1,$e1) = ($startQuery, $endQuery); # copy into variables that we might reorder
		($s1, $e1, $s2, $e2) = ($s2, $e2, $s1, $e1) if $s1 > $s2; # force regions to be ordered so that $s1 <= $s2
		# there are really only two cases for overlap... e1<e2 means hits aren't nested, otherwise second one (in new order) is nested in first (and quantities to threshold over are more trivial)
		if ($e1 < $e2) {
		    my $l = $e1-$s2; # in this case (not nesting) this is the size of the overlap in amino acids
		    # set thresholds and decide
		    # in first case, have an overlap regardless of the $lCut
		    # in second case set threshold only if lCut is true ( first special case lCut==0 (or false in general) behaves like lCut==inf... decide only using fCut!)
		    if ($l/$lRef > $fCut || ($lCut && $l > $lCut) ) { push @hitsOv, $hi; }
		    # else not an overlap
		}
		else { push @hitsOv, $hi; } # if e1 >= e2, f==1 guaranteed, and that's never good biologically, so let's stop it there and call it a disallowed overlap.  If we wanted to actually apply lCut, we'd have to figure out the length of the smallest region, but again let's not waste time doing that...
	    }
	    else { push @hitsOv, $hi; } # passes if strict overlaps
	}
	# look at next hit
	$i++;
	last if $i > $upperIni; # check that we don't go over bounds (happens if the very last range overlapped with our query)
	$lowStart = $hit2lowStart->[$i];
    }
    # done! return 
    return \@hitsOv;
}

sub ga {
    # applies gathering thresolds filter (from scratch, does not rely on pre-compted labels)
    my ($hits, $acc2ds2ga) = @_;
    my @hitsPass; # new hits we want
    foreach my $h (@$hits) {
	# calculate whether this hit passes GA or not (ignores overlaps from clans and modes, doesn't matter to us at this point)
	my $ds2ga = $acc2ds2ga->{$h->{acc}}; # copy down structure
	push @hitsPass, $h if
	    $h->{scoreSeq} >= $ds2ga->{s} && $h->{score} >= $ds2ga->{d};
    }
    return \@hitsPass;
}

sub labelGa {
    # creates the GA column to pre-compute if things pass the thresholds or not
    # ignores overlaps from clans and modes
    # order unchanged, hits modified by reference
    my ($hits, $acc2ds2ga) = @_;
    foreach my $h (@$hits) {
	my $ds2ga = $acc2ds2ga->{$h->{acc}}; # copy down structure
	$h->{GA} = ($h->{scoreSeq} >= $ds2ga->{s} && $h->{score} >= $ds2ga->{d}) ? 1 : 0; # passed gathering cutoff or not (boolean)
    }
}

sub filterBooleanCol {
    # is a generalization of filtering using any of the Pfam cutoffs we have precalculated into the columns (GA,TC,NC).
    # if the file is missing the required column, nothing passes
    my ($hits, $col) = @_; # col is GA,TC,NC, or something else that can be treated as a boolean
    my @hitsNew;
    foreach my $h (@$hits) {
	push @hitsNew, $h if $h->{$col};
    }
    return \@hitsNew;
}

sub filterEOrGa {
    # a filter for both CODD and Dpuc, so that GA always stays, and everything else has to satisfy an E-value cutoff
    my ($hits, $Ecut, $gaCol, $key) = @_;
    $key = 'E' unless defined $key;
    $gaCol = 'GA' unless defined $gaCol; # default is to actually use GA column, but may use a different one (i.e. 'Dpuc') to define trusted domains differently
    my @hitsNew;
    foreach my $h (@$hits) {
	push @hitsNew, $h if $h->{$gaCol} || $h->{$key} <= $Ecut;
    }
    return \@hitsNew;
}

1;

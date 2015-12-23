# Copyright (C) 2014 Alejandro Ochoa, Singh Research Group, Princeton University
# This file is part of dPUC.
# dPUC is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
# dPUC is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with dPUC.  If not, see <http://www.gnu.org/licenses/>.

package Dpuc;
our $VERSION = 2.06;
use lib '.';
use Domains;
use DpucPosElim;
use DpucLpSolve;
#use DpucLpSolveOld; # needed for -old! Not included in public release.
use DpucNetScores;
use DpucOvsCompact;
use EncodeIntPair;
use strict;

# 2015-08-05 13:52:55 EDT
# v2.03 - besides major version bump, removed explicit inputs to loadNet, which simply passes all to DpucNetScores::loadNet
# 2015-08-12 21:15:10 EDT
# - same version (haven't posted yet): removed removeOvsCodd and removeOvsPRank option!
# 2015-08-14 11:24:16 EDT
# - same version (haven't posted): removed all self/trans distinction and parameter that controlled it. the code that handled negative self-context scores also went away.
# 2015-10-30 22:04:07 EDT
# v2.04 - major version bump (reflects change in DpucNetScores prior parametrization)
# 2015-11-13 13:01:09 EST - v2.05
# - a few days ago added score halving option (kept old versions so we don't break -old)
# - today added additional options so sequence threshold gets moved to objective (and no actual separate thresholds are set!). Domain thresholds are set exactly, but sequence threshold is approximated.  The ratio of F=Td/Ts is used to consider when sequence thresholds are actually necessary (see global var below for threshold, which can be changed).
# 2015-12-23 16:10:57 EST - v2.06
# - major version bump (reflects change in DpucNet parser)

# package constants
my $scoreScale = 1000; # in posElim, we multiply context scores internally by this much, then turn into ints (so additions are faster), but we also have to keep this scale in mind when we combine context and HMM scores!  This loses some precision, but I think it's acceptable for the sake of speed, and also given the uncertainty involved in the actual calculation of context scores.
# useful for debugging the lp_solve end, not otherwise...
my $lpSolveVerbose = 0; # print extra messages
my $lpSolveDoSolve = 1; # skip solving if zero (can still get a lot of model stats this way, saving time)
my $timeoutDef = 1; # default is to give up with lp_solve after 10 seconds
my $EcutDef = 1e-2; # in terms of p-values (new setup), the least stringent one we can use (most stringent HMMER3 filter) is 1e-2 with new filters
our $doNewHitRefs = 1; # may edit this outside to change behavior as desired
our $newHalving = 1; # new change so scores are internally counted differently (same sometimes, halved some other times), and fixes sequence thresholds
our $cutTsd = 0.9; # decides if we don't use sequence threshold by thresholding F=Td/Ts

sub loadNet {
    # this is a wrapper around the file that normally loads the net, includes the typical pre-processing that is necessary for posElim or nonOvNeg

    my ($net, $acc2i, $netTs, $netTscores, $acc2cc) = DpucNetScores::loadNet($scoreScale, @_); # pass the same input here, plus this list of accs
    
    # copy perl structures into global and native C structures, should speed a few things up
    DpucPosElim::copyNetTscoresToC($netTs, $netTscores);
    
    # return this $netData array, contains an array ref of things we want (no point in presenting them separately on the outside, only dPUC needs them internally!
    return [$net, $acc2i, $acc2cc];
}

sub getInfoKeys {
    # this sub gets the list of keys that posElim and nonOvNeg %info hashes have (when combined)
    # they're returned sorted
    # getting them from these dummy function calls (both with trivial hits cases) ensures key lists are always up to date
    # - no need to keep a separate list that has to be manually ensured to match what these functions return
    # - disadvantage is we can't sort keys manually in a way that make sense
    my (undef, $info) = posElim([]);
    my (undef, $infoLp) = lpNonOvNeg([]);
    mergeHashes($info, $infoLp); # consolidate info hashes!
    return sort keys %$info; # return list directly (instead of as ref, as I usually do), it's small so no big deal
}

sub dpuc {
    my ($hits, $netData, $nestingNet, $acc2ds2t, $acc2ts, $timeout, $Ecuts, $foLps, $fCut, $lCut) = @_;
    # this sub efficiently computes dpuc predictions over a range of E-value cutoffs
    # the idea is to factor and skip work in these cases
    # - posElim preds are always nested by E-value, so we start from the previously-solved problem instead of starting from scratch
    # - nonOvNeg preds are a bit more complicated (not just nested), but solutions are guaranteed to be the same (when optimal) when no solution domains were eliminated by E-value cutoffs
    # - other data can be shared between the two and across E-values
    # ASSUMPTIONS: @$Ecuts has to be sorted descending (1000,900,...1,0.1,...), this program doesn't check for that but depends on that to be true
    # returns array of $hits (nonOvNeg solutions) that are parallel to the $Ecuts array
    # $foLps is optional, but if is specified, it should be an array of LP output files parallel to @$Ecuts, it will be ignored otherwise
    
    # OR, input $Ecuts can be scalar, for simplicity, but in that case the return is different...
    # output will be single data ($hitsDpuc, $hitsPE, $info) instead of three arrays each with a single member (again, only when input $Ecuts was scalar)
    # similarly, if defined, $foLps will be expected to be a scalar (a single file) instead of an array (or be ignored)
    # in the case of a single E-cut, this is as performant as it would be without this extra structure, so this is the main api for that case too
    
    # only confusing bit is then: if we pass an @$Ecuts array with a single ecut, we return the parallel arrays anyway (since it's more natural to think user expects this array setup in this case, to not complicate their side of things)
    
    my ($net, $acc2i, $acc2cc) = @$netData; # expand $netData
    
    # defaults
    $timeout = $timeoutDef unless defined $timeout;
    $Ecuts = $EcutDef unless defined $Ecuts; # set as scalar, not list, which as said above changes the return data structure, watch out!
    
    my $oneEcut = 0; # boolean, false by default
    if (!ref $Ecuts) { # if Ecuts is a scalar
	$oneEcut = 1; # make this boolean true
	$Ecuts = [$Ecuts]; # turn input into an array ref with the single original Ecut scalar in array
	$foLps = [$foLps] if defined $foLps && !ref $foLps; # similar processing if defined and scalar
    }
    # so after this step, $Ecuts is always an array, and so is $foLps if it was defined (loops don't have to think about this too much anymore)
    
    my @cut2hitsDpuc;
    my @cut2hitsPosElim; # for debugging/random interest
    my @cut2info; # for debugging/random interest
    my $maxE; # used to determine when rerunning nonOvNeg (actually dPUC as whole) is necessary...
    my $hitsRaw = $hits; # remember all input hits, for "dpucPlusRepeats" extension (not needed otherwise)
    my $hitsDpuc; # copy reference to previous solution when nonOvNeg is not redone
    
    # do pre-elimination step here, only once for all Ecuts, necessary for both posElim and nonOvNeg
    ($hits, my $hit2passSeq, my $errStr) = preElim($hits, $acc2ds2t, $acc2ts, $nestingNet, $fCut, $lCut, $acc2i, $acc2cc, $net);

    # if there was an error, return with it immediately
    return (undef, undef, undef, $errStr) if $errStr;
    
    my ($overlapNet, $scoresNet); # compute once, remember later, undef at first is fine
    
    # navigate cuts.  Start with the least stringent cutoffs, so we can analyze progressively
    foreach my $Ecut (@$Ecuts) {
	my $foLp = shift @$foLps if defined $foLps; # get next file in list, if list is defined (otherwise $foLp should be undefined)
	my $info = {}; # set info as empty hash, to distinguish when dPUC wasn't rerun
	# analyze old solution (with more permissive E-value cutoff), if available, to decide whether rerunning dPUC is necessary or not for more stringent E-value cutoff
	# only process if domains in solutions will be eliminated (otherwise solution is guaranteed to be the same as before)
	if (!defined $maxE || $maxE > $Ecut) {
	    # run the dPUC positive elimination, replace old hits because these are subsets as we advance @cuts
	    # it will compute the $overlapNet the first time, and remember it subsequent times
	    ($hits, $info, $overlapNet) = posElim($hits, $nestingNet, $acc2ts, $hit2passSeq, $Ecut, $overlapNet, $fCut, $lCut);
	    # remember what overlaps are forbidden, and store the final scores to make a single hash call per pair.
	    # only do first time
	    $scoresNet = getScoresNet_nonOvNeg($hits, $net, $overlapNet) unless defined $scoresNet;
	    # run the dPUC lp_solve general optimization, these are kept separate because they don't have an obvious relation as we advance @cuts
	    ($hitsDpuc, my $infoLp) = lpNonOvNeg($hits, $scoresNet, $nestingNet, $acc2ts, $hit2passSeq, $timeout, $foLp, $fCut, $lCut);
	    mergeHashes($info, $infoLp); # consolidate info hashes!
	    $maxE = 0; # identify worst E-value in solution to avoid rerunning dPUC next time if possible, initialize with zero (the best possible E-value), which covers the trivial no-hits case!
#	    if ($info->{lpStatus} eq 'OPTIMAL' || $info->{lpStatus} eq 'TRIVIAL') {
	    foreach my $h (@$hitsDpuc) {
		# go through data to find worst E-value, only consider non-GA domains, the only ones posElim can eliminate (otherwise lpNonOvNeg will be fed the same problem again)!
		$maxE = $h->{E} if !$h->{GA} && $maxE < $h->{E};
	    }
#	    }
#	    else { $maxE = $Ecut; } # if not optimal, setting maxE to current cut guarantees that dPUC will be run next time, as desired
	}
	push @cut2hitsDpuc, $hitsDpuc; # store in larger structure, regardless of lpStatus
	push @cut2hitsPosElim, $hits;
	push @cut2info, $info;
    }
    
    # done, return desired list of predictions
    if ($oneEcut) { return ($cut2hitsDpuc[0], $cut2hitsPosElim[0], $cut2info[0]); } # but if we passed a scalar E-value, return single data instead of parallel arrays
    else { return (\@cut2hitsDpuc, \@cut2hitsPosElim, \@cut2info); } # return parallel arrays
}

sub preElim {
    # this sub applies a few one-time-only processing steps to data
    my ($hits, $acc2ds2t, $acc2ts, $nestingNet, $fCut, $lCut, $acc2i, $acc2cc, $net) = @_;

    # first thing to do is map accs to integers, because this helps catch version mismatch errors (and we want to stop as soon as possible when that happens)
    # we'd like to check before filters have been applied! (filters reduce the change that we'll catch the error sooner)
    foreach my $h (@$hits) {
	$h->{acci} = $acc2i->{$h->{acc}}; # create acci column, things that interact with C need it
	# NOTE: now the %$acc2i map is only for valid nodes in network, so undefined values must be fatal!!!
	unless (defined $h->{acci}) {
	    my $errStr = "FATAL: domain family $h->{acc} was not part of context network!  The HMM DB version of the domain predictions and the dpucNet files probably differs.\n";
	    return (undef, undef, $errStr); # return right away
	}
    }

    # label domains that pass thresholds in the absence of context (this allows for shortcuts)
    $newHalving ? labelGaDpuc($hits, $acc2ds2t, $acc2ts) : Domains::labelGa($hits, $acc2ds2t); # calculate whether this hit passes GA or not, needed for filterCC
    
    # remove non-GA domains that are not "context-carrying" (domains that don't have positive context with anything, including themselves)
    # this filter is super quick and independent of E-value thresholds, and doesn't do anything when applied multiple times, so might as well put it out here
    $hits = filterCC($hits, $acc2cc);

    # this is done earlier in labelGaDpuc, so only have to do if it wasn't used (if $newHalving was false)
    unless ($newHalving) {
	# now permanently associate domain thresholds to domains, and store normalized score and acci too
	# this modifies the hits outside dpuc() but by adding columns that no other script uses (doesn't modify standard columns!), so it should have no side-effects
	foreach my $h (@$hits) {
	    my $td = $acc2ds2t->{$h->{acc}}{d}; # get domain threshold
	    $h->{td} = $td; # store threshold in hit as another column/key (needed by lp_solve unfortunately, when computing sequence threshold to un-normalize domain scores)
	    $h->{scoreNorm} = $h->{score} - $td; # compute and score normalized score too!
	}
    }
    
    # determine which sequence thresholds need to be evaluated, to try to save time (primarily a posElim shortcut; shouldn't be used in nonOvNeg)
    # this boolean depends only on which hits are GA (which are never eliminated by posElim) and which pfam thresholds are non-trivial.  Therefore, its results are valid for all Ecuts for a given protein
    my $hit2passSeq = getPassSeq($hits, $acc2ts);
    
    # yey, done!
    return ($hits, $hit2passSeq);
}

sub lpNonOvNeg {
    # this is the generalized sub that runs lpSolve with nonOvNeg on a set of hits
    # only funny param is $foLp, if provided, we only print the original formulation (not the timeout simplifications), just so you know!
    my ($hits, $scoresNet, $nestingNet, $acc2ts, $hit2passSeq, $timeout, $foLp, $fCut, $lCut) = @_;
    
    my $numHitsIn = scalar @$hits;
    
    # this handy sub transforms a single @$hits problem, given the general scoring structures, into the lp_solve problem, which gets solved
    my ($hitsPass, $info) = strucLp($hits, $scoresNet, $acc2ts, $hit2passSeq, $timeout, $foLp);
    
    # if we had an optimal solution, return!
    if ($info->{lpStatus} eq 'OPTIMAL' || $info->{lpStatus} eq 'TRIVIAL') {
	return ($hitsPass, $info);
    }
    # else continue...
    # this should only be the suboptimal or timeout case, but don't assume other errors can't occur
    
    # let's try to reasonably simplify the problem without just returning GA...
    # there are two ways of doing it
    # let's remove all overlaps while keeping the top scoring hits, then run lp_solve again
    my $hitsNonOv = Domains::removeOverlapsByScore($hits, $nestingNet, $fCut, $lCut); # remove overlaps as desired, also sorts hits by start
    # another way of removing overlaps, this time favoring homogeneity (which should be a better bet when we have tons of repeats with overlapping clan members)
    my $hitsNonOvPop = Domains::removeOverlapsByAccPopularity($hits, $nestingNet, $fCut, $lCut);
    my ($hitsPassNonOv, $infoNonOv); # we might not run these things, keep undefined in that case...
    my ($hitsPassNonOvPop, $infoNonOvPop); # ditto
    
    # try if we actually eliminated stuff
    # in one important case this is not true, so we don't keep wasting our time
    if (@$hitsNonOv < $numHitsIn){
	($hitsPassNonOv, $infoNonOv) = strucLp($hitsNonOv, $scoresNet, $acc2ts, $hit2passSeq, $timeout);
	$info->{lpTime} += $infoNonOv->{lpTime}; # always add new runtime to old runtime
    }
    if (@$hitsNonOvPop < $numHitsIn){
	($hitsPassNonOvPop, $infoNonOvPop) = strucLp($hitsNonOvPop, $scoresNet, $acc2ts, $hit2passSeq, $timeout);
	$info->{lpTime} += $infoNonOvPop->{lpTime}; # always add new runtime to old runtime
    }
    # decide if the second simplified version is better than the first...
    # first: only consider if the second one is defined
    # second: either first one is undefined, or both are defined but second score is better...
    if (defined $infoNonOvPop && (!defined $infoNonOv || $infoNonOvPop->{score} > $infoNonOv->{score}) ) {
	# replace first with second in these cases
	($hitsPassNonOv, $infoNonOv) = ($hitsPassNonOvPop, $infoNonOvPop);
    }
    # if this problem was solved at any capacity (optimally or suboptimally), return that
    # we can safely assume that if we couldn't solve any of the simplified problems, that there's no way we could have solved the original problem either
    if (defined $infoNonOv && ($infoNonOv->{lpStatus} eq 'OPTIMAL' || $infoNonOv->{lpStatus} eq 'TRIVIAL' || $infoNonOv->{lpStatus} eq 'SUB-OPTIMAL')) {
	# if original lpStatus was suboptimal, then we have two plausible solutions, the original and the simplified one...
	# decide which one to return based on total score
	if ($info->{lpStatus} eq 'SUB-OPTIMAL' && $info->{score} >= $infoNonOv->{score}) {
	    return ($hitsPass, $info); # return original one, pretend second thing wasn't run except for time (already updated)
	}
	# else continue...
	# if original timed out, or if it was suboptimal but the new solution was better, return new sol
	$infoNonOv->{lpTime} = $info->{lpTime}; # put cumulative time back on this thing (only original hash had such cumulative time)
	$infoNonOv->{lpStatus} = $info->{lpStatus}.','.$infoNonOv->{lpStatus}; # lpStatus is now the combination of the original and the second lpStatus, just to be informative
	$infoNonOv->{lpNonZeroes} = $info->{lpNonZeroes}; # however, we won't increase the size of the problem, leave original size
	# should probably add lpNodesProc and lpIter, but meh, not too useful if not debugging/optimizing, and in almost all cases this simplified problem won't happen anyway, so meh again
	return ($hitsPassNonOv, $infoNonOv);
    }
    # else continue...
    $hits = $hitsNonOv; # use these for further simplification below...
    my $infoNonOvLpStatus = $infoNonOv->{lpStatus}; # might use this below...

    # I think this doesn't happen anymore, will have to test in hard problems (newFdr? newFdrOrthoAlns? change in dpucNet params?)
    # remember to get back to this... will keep old version for now...
    ################
    # if we still couldn't solve this problem (double timeout?), return top non-overlapping GA, even if it's empty
    # but we still want to process cause negative context will bump down some GA...
    my $hitsNonOvGa = Domains::filterBooleanCol($hits, 'GA'); # hits were already sorted by start, and this operation does not change sorting, so this is it
    my ($hitsPassNonOvGa, $infoNonOvGa); # we might not run this thing, keep undefined in that case...
    # try if we actually eliminated stuff
    # in one important case this is not true, so we don't keep wasting our time
    if (@$hitsNonOvGa < $numHitsIn){
	($hitsPassNonOvGa, $infoNonOvGa) = strucLp($hitsNonOvGa, $scoresNet, $acc2ts, $hit2passSeq, $timeout);
	$infoNonOvGa->{lpTime} += $info->{lpTime}; # always add new runtime to old runtime, but assing to new info hash this time (old hash definitely won't be returned here)
	$infoNonOvGa->{lpStatus} = join ',', grep { defined } $info->{lpStatus}, $infoNonOvLpStatus, $infoNonOvGa->{lpStatus}; # lpStatus is now the combination of the original and the subsequent lpStatuss, just to be informative
	$infoNonOvGa->{lpNonZeroes} = $info->{lpNonZeroes}; # however, we won't increase the size of the problem, leave original size
	# return whatever we got if optimal or suboptimal
	if ($infoNonOvGa->{lpStatus} =~ /(OPTIMAL|TRIVIAL)$/) { # comparison is more awkward because statuses are mixed by this point, but meh, ending should be this to pass
	    return ($hitsPassNonOvGa, $infoNonOvGa);
	}
    }
    # I'll be damned, but if this didn't get solved either, just return GA without negative context processing (should be exactly GA, old regular non-context output, ignoring clan filtering [so disallowing overlaps except those explicitly allowed, instead of the oposite, which is to disallow overlaps explicitly by clans])
    # Need to wrapUp here, in all other cases strucLp() does that but not here!
    $hitsNonOvGa = wrapUp_nonOvNeg($hitsNonOvGa, $scoresNet); # write final context scores into domain objects (hashes)
    $info->{score} = getScoreManually($hitsNonOvGa); # have to set total score manually to the right value
    return ($hitsNonOvGa, $info); # only info that made it this far was original info, return that
}

sub strucLp {
    # this sub takes a particular hits structure, and a few general data structures, and returns the text corresponding to the lp_solve problem to solve, computing all the necessary intermediate steps
    # also detects and returns trivial problems to spare us calling lp_solve
    my ($hits, $scoresNet, $acc2ts, $hit2passSeq, $timeout, $foLp) = @_;
    
    # somewhat often there are no input domains, bail out early in that case!
    # equally trivial is the 1-domain case, there's no context to figure things out for... (but score is the domain's score)
    my $numHitsIn = scalar @$hits;
    # construct the info hash to return in that case
    my %infoTrivial = (
	'lpStatus' => 'TRIVIAL',
	'lpTime' => 0, 
	'lpNonZeroes' => 0, 
	'score' => 0, # will be updated below if we have 1 domain rather than zero...
	'lpNodesProc' => 0, 
	'lpIter' => 0,
	'lpIn' => $numHitsIn,
	'lpOut' => $numHitsIn, # in these trivial case number of domains in and out are the same (again, 0 or 1, or all)
	'lpCols' => 0,
	'lpRows' => 0,
	);
    unless ($numHitsIn > 1) {
	if ($numHitsIn) { # need to do a bit more
	    my $h = $hits->[0]; # get the one hit
	    $infoTrivial{score} = $h->{scoreNorm}; # total score is simply normalized score of the one domain, no context)
	    # also fill in columns that context data should have, in this case it's also trivial and doesn't differ from non-context in the regular columns of qw(score scoreSeq), so it's ok to modify directly instead of creating a copy
	    # not having this properly defined can mess up downstream analyses, and at best have perl complain a lot about undefined values
	    # (scoreSeq is never used directly, but before eliminations could have had a larger value than just the single surviving domain's score, so here we make sure they agree, as wrapUp() would have done)
	    $h->{scoreHmm} = $h->{score};
	    $h->{scoreSeq} = $h->{score};
	    $h->{scoreHmmSeq}  = $h->{score};
	    $h->{scoreContext} = 0;
	    $h->{scoreContextSeq} = 0;
	}
	return ($hits, \%infoTrivial);
    }
    
    # next "trivial" case is no overlaps and all positive context
    # contingent on posElim having been run already (to ensure thresholds were satisfied).  In that case the positive score is actually the real score, and it's obviously maximized by including everything.
    # one huge Pf protein is so large lpSolve doesn't like it, but satisfies this condition and is therefore actually easy
    my $hasOvs = 0; # booleans that keep track of what we want to know...
    my $hasNegs = 0;
    # navigate all domain pairs
    for (my $i = 0; $i < $numHitsIn; $i++) {
	my $hi = $hits->[$i];
	my $scoresNetHere = $scoresNet->{$hi};
	for (my $j = $i+1; $j < $numHitsIn; $j++) {
	    my $hj = $hits->[$j];
	    # if we found a score...
	    my $Cij = $scoresNetHere->{$hj};
	    if (defined $Cij) {
		$hasNegs = 1 if $Cij < 0; # mark if we found a negative score
	    }
	    # a missing score is necessarily a disallowed overlap, mark overlaps boolean
	    else { $hasOvs = 1; }
	    last if $hasOvs && $hasNegs; # stop looking if both booleans have been set.  this only breaks inner loop
	}
	last if $hasOvs && $hasNegs; # stop looking if both booleans have been set.  this breaks outer loop
    }
    if (!$hasOvs && !$hasNegs) {
	# this is another trivial case, wrap up and return what we have
	$hits = wrapUp_nonOvNeg($hits, $scoresNet); # write final context scores into domain objects (hashes)
	$infoTrivial{score} = getScoreManually($hits); # have to set total score manually to the right value
	return ($hits, \%infoTrivial);
    }
    
    # not trivial, have to use lp_solve...
    # restructure and send to lp_solve via C binding, get answers!
#    my ($hitsPass, $info) = DpucLpSolveOld::strucLpForC($hits, $scoresNet, $acc2ts, $hit2passSeq, $timeout, $foLp); # old version can be accessed this way! unless $newHalving
    my ($hitsPass, $info) = strucLpForC($hits, $scoresNet, $acc2ts, $timeout, $foLp);
    
    # write final context scores into domain objects (hashes), return!
    $hitsPass = wrapUp_nonOvNeg($hitsPass, $scoresNet);
    return ($hitsPass, $info);
}

sub getScoreManually {
    # normally lp_solve sets the total score of the solution, but it has to be done manually in trivial cases
    # NOTE: do this after wrapUp, so the domain "scores" already include context
    my ($hits) = @_;
    my $score = 0; # the sum of scores we want
    foreach my $h (@$hits) {
	$score += $h->{score} - $h->{td}; # score (including context) minus domain threshold
    }
    return $score; # done, return the score
}

sub strucLpForC {
    # we want to write the objective function, but we get the edge-defining constraints and declarations trivially, so might as well get them
    # the scores net already excludes scores between overlapping domains, so we won't waste text writing that, and we won't define those edges that cannot be used!
    # 2015-11-21: NOTE this is new version, which assumes newHalving=1 throughout... old version was preserved in DpucLpSolveOld.pm
    my ($hits, $scoresNet, $acc2ts, $timeout, $foLp) = @_;
    
    # not sure what Inline::C does if $foLp is undefined, better set it to an empty string
    $foLp = '' unless defined $foLp;
    
    # we use this all the time
    my $n = scalar @$hits;
    
    # first we need to map domains to variable indexes (which otherwise have no distinctive names)
    # since we also declare edges and families as variables, we need an easy way of mapping to final indexes without a particularly clean structure in mind.  Let's textify variables in the old way and map to contiguous indexes via this hash
    my @varNames = (); # list of names!
    my @obj = (0); # old setup had this defined automatically by Si, but new setup aims to phase that out and set objective more directly!!! To avoid casting problems, let's define the first element as zero (lp_solve ignores it)
    # domain states Di are 1..$n, the only answers we're interested in
    my %acc2is; # lastly, get all indexes per family (for sequence threshold)
    for (my $i = 0; $i < $n; $i++) {
	$varNames[$i+1] = 'x'.$i; # add to map, 1-based
	my $hi = $hits->[$i]; # get hit #i
	# in new setup, let's populate the objective function
	# using 1-based notation here, and storing objective non-sparsely!!!
	$obj[$i+1] = $hi->{scoreNorm}; # normalized score (Hi-Ti)
	# try to organize "sequence" threshold stuff
	my $acc = $hi->{acc};
	# add if $ts is there!
	push @{$acc2is{$acc}}, $i+1 if $acc2ts->{$acc}; # add to list of domains in family by index, 1-based
    }
    
    # start indexes after Di's and Si's (for edges first down here, then families)
    my $iNext = $n+1; # remember last index not used yet, for growing list
    
    # now let's add edges that have scores (other edges aren't allowed by overlap, they won't be defined)
    my @edgeDefs = (); # sparse edge definitions
    my @ovs = (); # sparse overlap list
    for (my $i = 0; $i < $n; $i++) {
	my $hi = $hits->[$i];
	my $scoresNetHere = $scoresNet->{$hi};
	for (my $j = $i+1; $j < $n; $j++) {
	    my $hj = $hits->[$j];
	    my $Cij = $scoresNetHere->{$hj};
	    if (defined $Cij) {
		# there's a score, so we need an edge variable
		$varNames[$iNext] = 'x'.$i.'_'.$j;
		
		# lastly, add an edge definition constraint, which is extremely sparse...
		push @edgeDefs, [$i+1, $j+1, $iNext]; # indexes have to be 1-based, not 0-based as in map, but $iNext is already good!
		#		my @sparseRow = (1,1,-2); # always the same for all sparse rows, will set this in C directly
		# constraint type and bounds are also set directly in C
		# in new setup, we put scores in objective and nowhere else
		# note no halving is needed here, since now we don't double-count things!!!
		$obj[$iNext] = $Cij;
		$iNext++; # increment for next round
	    }
	    # a missing score in the new setup is forcibly a disallowed overlap, write constraint
	    else { push @ovs, [$i+1, $j+1]; } # indexes also 1-based
	}
    }
    
    # use bulky outside function to try to compact overlap definitions, which will make the LP faster!
    # nothing to consider if overlaps have already been removed! (pass an empty array ref instead)
    my $ovsCompact = DpucOvsCompact::compactOvs(\@ovs);
    
    # lastly, define the "sequence" thresholds, each of which requires three statements
    my @seqDefs = (); # this is an unusual structure but it minimizes what we send to C while being complete
    while (my ($acc, $is) = each %acc2is) {
	$varNames[$iNext] = $acc; # assign latest index to this family
	push @seqDefs, [$iNext, @$is]; # need family index and list of family members, we sort this out in C
	$obj[$iNext] = -$acc2ts->{$acc}; # add threshold to objective, we don't need it anywhere else!!!
	$iNext++; # increase for next round
    }

    # done! time to send to C code! get answers!
    no warnings 'uninitialized'; # fixes a complaint that the below function call has when @seqDefs are empty arrays.  The exact warning is "Use of uninitialized value in subroutine entry at Dpuc.pm line xxx."
    my ($sol, $lpStatus, $lpTime, $numRows, $lpNonZeroes, $score, $lpNodesProc, $lpIter) = DpucLpSolve::sendToLpSolve($n, $iNext-1, \@obj, \@edgeDefs, $ovsCompact, \@seqDefs, $timeout, \@varNames, $lpSolveVerbose, $lpSolveDoSolve, $foLp);
    # we got things that make sense in C, but we can clean it up a bit for Perl...
    $lpStatus =~ s/ solution//; # this is a stupid string, remove this useless tail
    $lpStatus = 'TIMEOUT' if $lpStatus =~ /timeout/;
    # now let's return the subset of hits we want
    my @hitsPass;
    if ($lpStatus eq 'OPTIMAL' || $lpStatus eq 'SUB-OPTIMAL') { # solution stays empty otherwise
	for (my $i = 0; $i < $n; $i++) {
	    push @hitsPass, $hits->[$i] if int $sol->[$i]+1/2; # sol_i should be either 1 or zero, acts as boolean!  Round to closest int just in case though, sometimes it's not exactly zero and passes as "true" but shouldn't
	}
    }
    # now let's make a hash that passes all the secondary info (the non-hits essentially)
    my %info = (
	'lpStatus' => $lpStatus, 
	'lpTime' => $lpTime, 
	'lpNonZeroes' => $lpNonZeroes, 
	'score' => $score, 
	'lpNodesProc' => $lpNodesProc, 
	'lpIter' => $lpIter,
	'lpIn' => $n, # number of domains that went into the problem
	'lpOut' => scalar @hitsPass, # and number of domains out
	'lpCols' => $iNext-1, # number of variables in problem (matrix columns)
	'lpRows' => $numRows, # number of constraints in problem (matrix rows)
	);
    # awesome, return hits sorted naturally!
    return (Domains::sortByCol(\@hitsPass, 'start'), \%info);
}

sub posElim {
    # this sub does everything that must be done for posElim, while receiving as inputs all that can be factored out, and the problem itself ($hits)
    my ($hits, $nestingNet, $acc2ts, $hit2passSeq, $Ecut, $overlapNet, $fCut, $lCut) = @_;
    
    # remove hits that are neither GA nor satisfy our E-value cutoff (same as CODD)
    # this overwrites the internal $hits array ref, but the external $hits is unmodified!
    $hits = Domains::filterEOrGa($hits, $Ecut);
    
    # make info hash we might return right away...
    my %info = (
	'peIter' => 0,
	'peIn' => scalar @$hits, # total number of starting hits, E-value eliminations shouldn't count!
	'peGa' => 0,
	'pePassNoGa' => 0,
	'peNoPass' => 0,
	);
    
    # take care of trivial cases first, instead of initializing C structures and all the fuss...
    # no domains is super easy...
    if ($info{peIn} == 0) { return ($hits, \%info); }
    # one domain has no context, so it only passes if it's GA
    elsif ($info{peIn} == 1) {
	if ($hits->[0]{GA}) {
	    $info{peGa} = 1; # update info hash to say one GA domain was there
	    return ($hits, \%info);
	}
	else {
	    $info{peNoPass} = 1; # instead, update info hash to say one domain was eliminated
	    return ([], \%info); # make sure empty array is returned, and not original hits array
	}
    }
    # else keep going!
    
    my $numPass = $info{peIn}; # main number we need to keep for loops
    $info{peGa} = countGa($hits); # just for stats, use separate function
    
    # an expensive part of DpucPosElim::getScoresNet_posElim is mapping hits to their accessions, so we can get their scores
    # let's do that mapping in O(n) instead of in the O(n^2) loop
    my @hitis2acc;
    for (my $i = 0; $i < $numPass; $i++) {
	$hitis2acc[$i] = $hits->[$i]{acci}; # store in final structure, hopefully as a number internally
    }
    # now actually get scores, store them in C structure that Perl can't see directly
    DpucPosElim::getScoresNet_posElim(\@hitis2acc);
    
    # need to initialize @$hitis with all the indexes in $hits.  Do as array reference for easy overwritting with elimMaxScore.
    my $hitis = [0..$numPass-1];
    
    # initialize loop vars
    my $iter = 0;
    my ($numPassNow); # undefined is ok
    
    # very first time we run posElim we will have too many hits to eliminate, and it's not worth computing the $overlapNet until the set is much smaller
    # however, subsequent runs (for more stringent E-values) can take the previously computed $overlapNet and skip directly to that part of the elimination...
    unless (defined $overlapNet) {
	# do massive eliminations by setting cutoffs at the domain and sequence level
	$iter++;
	$hitis = elimMaxScore($hits, $hitis, $acc2ts, $hit2passSeq);
	$numPassNow = scalar @$hitis;
	while ($numPass != $numPassNow) { # this condition means we haven't converged
	    $numPass = $numPassNow; # update counts
	    # run again
	    $iter++;
	    $hitis = elimMaxScore($hits, $hitis, $acc2ts, $hit2passSeq);
	    $numPassNow = scalar @$hitis;
	}
	
	# now that total hits are largely reduced, include overlap information into eliminations (too expensive to do at first)
	my $hitsLeft = [@{$hits}[@$hitis]]; # this list includes only hits that have passed.  This reduced list requires overlap checking.
	($overlapNet) = Domains::overlapNetBinarySearch($hitsLeft, $nestingNet, $fCut, $lCut); # get which hits actually overlap and aren't allowed to, mapping from reference
    }
    
    # need to do this always, to use $overlapNet info
    zeroScoresOvPairs($hits, $hitis, $overlapNet); # update C's $scoresNet half matrix by zeroing the scores that aren't allowed due to overlaps
    
    # do massive eliminations by setting cutoffs at the domain and sequence level
    $iter++;
    $hitis = elimMaxScore($hits, $hitis, $acc2ts, $hit2passSeq);
    $numPassNow = scalar @$hitis;
    while ($numPass != $numPassNow) { # this condition means we haven't converged
	$numPass = $numPassNow; # update counts
	# run again
	$iter++;
	$hitis = elimMaxScore($hits, $hitis, $acc2ts, $hit2passSeq);
	$numPassNow = scalar @$hitis;
    }
    
    # done, turn passed indeces into a finalized domain list, to return
    $hits = [@{$hits}[@$hitis]]; # this updates hits to include only the hits that have passed.
    
    # done, we got the posElim solution to the input set of hits!
    # update stats in existing %info hash
    # note $iter and $numPass were updated in loops above so couldn't be assigned until now
    $info{peIter} = $iter;
    $info{pePassNoGa} = $numPass - $info{peGa}; # NOTE: Dpuc < 2 used to return $numPass to the exclusion of $peGa, but Dpuc 2.0's $numPass include $peGa !  For compatibility, we subtract $peGa from the new $numPass to produce the same old result
    $info{peNoPass} = $info{peIn} - $numPass; # spell out the meaning of this number
    # return $overlapNet too, since hit refs don't change when needed in lpNonOvNeg!
    return ($hits, \%info, $overlapNet);
}

# actually this had to be split due to $scoresNet differences (posElim uses double array, nonOvNeg uses double hash)
# since, wrapUp_posElim has been deprecated, so this is actually the only version that exists, but we have to note that the inputs are nonOvNeg specific!
sub wrapUp_nonOvNeg {
    my ($hits, $scoresNet) = @_;
    # calculate all context scores
    my $hit2contextScore = getSumScoreEdgesFromNode_nonOvNeg($hits, $scoresNet);
    # now get total domain and sequence scores, as well as HMM seq scores (cause we don't trust the ones in the structure, might not hold after eliminations)
    my %hit2scoreDom;
    my %acc2scoreSeq;
    my %acc2scoreHmmSeq;
    foreach my $h (@$hits) {
	my $contextScore = $hit2contextScore->{$h}; # may be undefined if sum was zero
	my ($scoreDom, $acc) = @{$h}{qw(score acc)};
	$acc2scoreHmmSeq{$acc} += $scoreDom;
	$scoreDom += $contextScore/2 if $contextScore; # note here actual context scores (sum of edges) are shared between domains, so we count half toward each (so sum makes sense)
	$hit2scoreDom{$h} = $scoreDom;
	$acc2scoreSeq{$acc} += $scoreDom;
    }
    # now make new hits with the modified scores (original hit scores will not be touched!)
    my @hitsNew;
    foreach my $h (@$hits) {
	my $acc = $h->{acc};
	# get total score
	my $scoreTotalDom = sprintf '%0.1f', $hit2scoreDom{$h};
	my $scoreTotalSeq = sprintf '%0.1f', $acc2scoreSeq{$acc};
	# get HMM score only (old score)
	my $scoreHmmDom = $h->{score};
	my $scoreHmmSeq = sprintf '%0.1f', $acc2scoreHmmSeq{$acc};
	# infer the context score
	my $scoreContextDom = sprintf '%0.1f', $scoreTotalDom-$scoreHmmDom;
	my $scoreContextSeq = sprintf '%0.1f', $scoreTotalSeq-$scoreHmmSeq;
	# now store these numbers again.  Note score and scoreSeq will now refer to the total scores
	my $h2 = $doNewHitRefs ? {%$h} : $h; # we might want to make a new copy of hits to edit without changing originals (default), or we might want to overwrite original hits (in rare applications it is more useful)
#	my %hNew = %$h; # this should copy the hits, so we can edit them without changing the originals!
	$h2->{score} = $scoreTotalDom;
	$h2->{scoreSeq} = $scoreTotalSeq;
	$h2->{scoreHmm} = $scoreHmmDom;
	$h2->{scoreHmmSeq} = $scoreHmmSeq;
	$h2->{scoreContext} = $scoreContextDom;
	$h2->{scoreContextSeq} = $scoreContextSeq;
	push @hitsNew, $h2;
    }
    # return new hits with modified scores
    return \@hitsNew;
}

# this is for posElim only (eliminations works with the positive scores assumption, and besides nonOvNeg only eliminates through lp_solve, outside of Perl)
sub elimMaxScore {
    my ($hits, $hitis, $acc2ts, $hit2passSeq) = @_;
    # go through hitis (which passed previously) and eliminate those that have a partner from those that don't, and those whose partners may raise their scores enough to pass the threshold
    my $numHitis = scalar @$hitis;
    # calculate all context scores
    my $hitk2contextScore = DpucPosElim::getSumScoreEdgesFromNode_posElim($hitis);
    # now get sequence scores
    my %acc2scoreSeq;
    # calculate the sum of scores for each family that needs the sequence thresholds
    for (my $k = 0; $k < $numHitis; $k++) {
	# for now we need to get the original hit structure, to use the hit2blah mappings, but we should consider permanently shifting into indeces for all code!
	my $h = $hits->[$hitis->[$k]];
	# decide whether a family needs a sequence score by looking at the passSeq boolean (family may pass automatically, so why bother with the sequence score)
	next if $hit2passSeq->{$h};
	# get all these quantities of the hit at the same time (is it faster?)
###	my ($score, $acc) = @{$h}{qw(score acc)}; # weird hash slice ### old version
	my ($score, $acc) = @{$h}{qw(scoreNorm acc)}; # weird hash slice ### $newHalving==1 version
	$acc2scoreSeq{$acc} += $score + $hitk2contextScore->[$k]/$scoreScale; # Also need to reverse context score int scaling!
    }
    # now for each family that didn't pass automatically, remember if the sequence score passed this threshold or not.
    # this passing is true for all members of the family, so it makes sense to calculate it separately from the individual domain thresholds (this way we're factoring that calculation out!)
    my %acc2passSeq; # remember if the sequence score passed the sequence threshold
    while (my ($acc, $scoreSeq) = each %acc2scoreSeq) {
	# get precalculated sequence threshold
	my $ts = $acc2ts->{$acc};
	# now we check if the threshold is passed by the sequence score
	$acc2passSeq{$acc} = 1 if $scoreSeq >= $ts;
    }
    # now browse all the domains, decide whether they still pass or not.
    my @hitisPass; # has partners and max score allows domain to pass threshold
    for (my $k = 0; $k < $numHitis; $k++) {
	my $i = $hitis->[$k]; # get $i (index relative to all hits)
	# for now we need to get the original hit structure, to use the hit2blah mappings, but we should consider permanently shifting into indeces for all code!
	my $h = $hits->[$i];
	# hits with positive scores always pass automatically
	if ($h->{GA}) { push @hitisPass, $i; }
	else {
	    # if this had no partners, no point in looking at thresholds (since this isn't already GA, then it automatically doesn't pass)
	    my $contextScore = $hitk2contextScore->[$k];
	    if ($contextScore) {
		# test if it passed the threshold here or not
		my $scoreDom = $h->{scoreNorm} + $contextScore/$scoreScale; # add self HMM score (normalized by threshold) and max context score, reverse context score int scaling
		if ($scoreDom >= 0) { # has to be non-negative (because it's already normalized)
		    # domain threshold was passed, now check sequence threshold.
		    # the first part is true if the domain automatically passes the sequence threshold (due to many reasons).
		    # the second part is true if the family of this domain passed the sequence threshold (which was only evaluated for domains that didn't pass automatically, so we need to put it second)
		    if ($hit2passSeq->{$h} || $acc2passSeq{$h->{acc}}) {
			push @hitisPass, $i; # this passed
		    }
		}
		# else max score wasn't good enough
	    }
	}
    }
    return \@hitisPass;
}

# this is shared by both posElim and nonOvNeg
sub getPassSeq {
    # this function determines which seqThreshs are necessary, a shortcut that might save us posElim and nonOvNeg time
    my ($hits, $acc2ts) = @_;
    
    # make this structure that keeps track of which families pass GA, and we'll mark them and additional family members as automatically passing all sequence thresholds.
    # this is guaranteed as long as there's a single domain (trivial) or all domains have positive scores
    # (copied from another part of the code...): an assumption breaks down when using negative context, but I doubt that it will matter in practice, let's just trust GAs when it comes to the sequence threshold.
    # anomalous cases might not be eliminated when they should, but this filter is strong enough that we won't care about this problem
    my %acc2pass;
    # record which accessions pass either because one of their domains is GA, or there is no sequence threshold to impose
    # this pass sets a single passing value per acc
    foreach my $h (@$hits) {
	my $acc = $h->{acc};
	$acc2pass{$acc} = 1 if $h->{GA} || !defined $acc2ts->{$acc};
    }
    # lastly, permanently map automatic passing to each hit
    my %hit2passSeq; # they pass sequence threshold as a family (no need to calculate sequence score)
    foreach my $h (@$hits) {
	$hit2passSeq{$h} = 1 if $acc2pass{$h->{acc}};
    }
    return \%hit2passSeq;
}

# the only reason we have a nonOvNeg version is that this takes the scoresNet that is a double hash, too painful to rewrite the rest of nonOvNeg to use the double array system that posElim uses.  This is used specifically for the wrapUp_nonOvNeg, and nothing else.
sub getSumScoreEdgesFromNode_nonOvNeg {
    # returns a hash that maps each hit to its total context score (should be halved in some contexts since edges are shared)
    my ($hits, $scoresNet) = @_;
    my %hit2score; # this is the hash we want
    foreach my $h (@$hits) {
	my $scoresNetThisHit = $scoresNet->{$h}; # copy down structure, should reduce the number of hash calls
	my $sumContextScore = 0; # this is the context score it gets from all possible interactions
	foreach my $h2 (@$hits) {
	    my $edgeScore = $scoresNetThisHit->{$h2}; # score of this pair
	    $sumContextScore += $edgeScore if defined $edgeScore; # will be undefined if this wasn't a context pair
	}
	$hit2score{$h} = $sumContextScore if $sumContextScore; # store the total score if non-zero (saves a bit of memory and hash lookup time?)
    }
    return \%hit2score;
}

# this is the nonOvNeg version (uses negative scores, special rules apply to them!).
# A lot of the code is shared with *_posElim, but the tighter integration of each for their particular needs is worth for the speed differences I anticipate.
sub getScoresNet_nonOvNeg {
    # this sub combines the overlap disallowed rule, with the context net that uses accessions, to make a net that maps each hit pair (as refs, not Pfam accessions) to the scores they should map to.  This should significantly speed up iterations (without affecting single-iteration runs a whole lot)
    # providing an empty overlap hash should have no effect in processing, I only hope perl doesn't complain
    my ($hits, $net, $overlapNet) = @_;
    my %scoresNet; # the structure we want, both ways
    # since most edges are probably negative, save this score from the network to save on hash calls
    my $neg = $net->{NEG}{NEG};
    # we can halve the number of operations by not repeating pairs in loop!
    my $numHits = scalar @$hits;
    for (my $i = 0; $i < $numHits; $i++) {
	my $hi = $hits->[$i];
	my $acci = $hi->{acc};
	my $netThisAcc = $net->{$acci}; # copy down structure, shoud reduce the number of hash calls
	for (my $j = $i+1; $j < $numHits; $j++) {
	    my $hj = $hits->[$j];
	    # check for overlap and ask the nesting net if it should be allowed (precalculated)
	    if (!$overlapNet->{$hi}{$hj}) {
		### This is where the posElim-nonOvNeg difference lies, changed for Dpuc2
		# here for Dpuc11 we used to force counts with negative scores to zero scores, we don't do that anymore (but all counts get positive scores now anyway, we're simply not checking anymore)
		# this part also changed from Dpuc11, negative scores were allowed to take arbitrary p_i,p_j background values, but now they're explicitly uniform.
		# lastly, Dpuc11 allowed for some domains to not participate in context at all (positive or negative), in the case that they were only present in eliminated architectures.  Now we make exceptions for those architectures, so all domains have context.  Missing positive context implies negative context!
		my $accj = $hj->{acc};
		my $edgeScore = $netThisAcc->{$accj}; # score of this pair
		$edgeScore = $neg unless defined $edgeScore;
		### end posElim-nonOvNeg difference
		# save both ways!
		$scoresNet{$hi}{$hj} = $edgeScore;
		$scoresNet{$hj}{$hi} = $edgeScore;
	    }
	}
    }
    return \%scoresNet;
}

# for posElim only, since it uses the half matrix $scoresNet (not the double hash net that nonOvNeg uses!)
sub zeroScoresOvPairs {
    # this sub updates the $scoresNet half matrix by zeroing the scores that aren't allowed due to overlaps
    my ($hits, $hitis, $overlapNet) = @_;
    
    # C-like initialize all vars outside of potentially tight loop
    my ($i,$j,$k,$l,$ti) = (0,0,0,0,0);
    my ($hi, $hj) = ({},{});
    
    # browse the @$hitis pairs with triangle indeces of indeces
    my $numHitis = scalar @$hitis;
    for ($k = 1; $k < $numHitis; $k++) { # starts at k=1 cause l<k, so k=0 doesn't have an l!
	$i = $hitis->[$k];
	$hi = $hits->[$i];
	for ($l = 0; $l < $k; $l++) {
	    $j = $hitis->[$l];
	    $hj = $hits->[$j];
	    # zero score if we have a disallowed overlap
	    if ($overlapNet->{$hi}{$hj}) { # works cause net is both ways
		# zero score by calling the C function, which interacts with the C data directly
		DpucPosElim::zeroScore(EncodeIntPair::t($i, $j));
	    }
	}
    }
    # done, nothing to return since $scoresNet was altered in C space
}

# quick filter independent of E-value cutoffs, used once in dpuc before posElim
sub filterCC {
    # ask that non-positive domains are "context-carrying", otherwise they don't pass
    my ($hits, $acc2cc) = @_;
    my @hitsPass;
    foreach my $h (@$hits) {
	# GA pass automatically
	# else query map $acc2cc whether an accession is context-carrying
	if ($h->{GA} || $acc2cc->{$h->{acc}}) {
	    push @hitsPass, $h;
	}
    }
    return \@hitsPass;
}

sub countGa {
    # as simple as it sounds, just return number of positive normalized score (old:GA) domains
    # not important for Dpuc, just stats collection
    my $hits = shift;
    my $c = 0;
    foreach my $h (@$hits) { $c++ if $h->{GA}; }
    return $c;
}

sub getPfamThresholdsSeq {
    # returns %acc2ts hash excluding the unnecessary thresholds!
    # NOTE: new behavior ($newHalving==1) sets a filter on ratio of thresholds (to prevent using sequence thresholds that are actually very close to the domain thresholds, so they're practically redundant; omitting them is much more efficient), and stores the delta of thresholds (Ts-Td, which is what we actually have in the objective function; Ts by itself is never used anymore!)
    # input are the pfam GA thresholds
    my ($acc2ds2t) = @_;
    my %acc2ts;
    while (my ($acc, $ds2t) = each %$acc2ds2t) {
	# get thresholds of interest (used more under new rules)
	my $Ts = $ds2t->{s};
	my $Td = $ds2t->{d};
	# this sequence threshold has to be applied, since it's higher than the domain threshold (so it's not trivially satisfied)
	if ($Ts > $Td) {
	    if ($newHalving) {
		# add *delta* of thresholds, and only if sequence threshold is sufficiently large relative to domain threshold
		$acc2ts{$acc} = $Ts-$Td if $Td/$Ts < $cutTsd;
	    } else {
		$acc2ts{$acc} = $Ts; # always add these cases (old behavior)
	    }
	}
    }
    return \%acc2ts;
}

sub labelGaDpuc {
    # 2015-11-13 14:09:58 EST
    # like Domains::labelGa(), but here we use a modified sequence threshold that is more compatible with our new objective function... always sets domain thresholds (as regular GA), but sequence thresholds are different (in non-context experiments, this worked well)
    # GA precomputes whether domains pass thresholds without context, so here we make that consistent with our change!
    # 
    # OLD description that still applies:
    # creates the GA column to pre-compute if things pass the thresholds or not
    # ignores overlaps from clans and modes
    # order unchanged, hits modified by reference
    my ($hits, $acc2ds2t, $acc2ts) = @_;
    die "Fatal: labelGaDpuc() should only be used if newHalving is true!!!" unless $newHalving; # a sanity check
    my %acc2scoreSeq; # increment these for non-trivial cases!
    
    # first look looks at domain thresholds and collects data for sequence threshold, set later
    # integral to the new system are the normalized scores, so let's compute them here too!
    foreach my $h (@$hits) {
	my $acc = $h->{acc};
	my $td = $acc2ds2t->{$acc}{d}; # get domain threshold
	$h->{td} = $td; # store threshold in hit as another column/key (needed by lp_solve unfortunately, when computing sequence threshold to un-normalize domain scores, and for summaries of trivial cases)
	$h->{scoreNorm} = $h->{score} - $td; # compute and score normalized score too!
	# mark if domain threshold passed (always needed in the new setup)
	$h->{GA} = $h->{scoreNorm} >= 0 ? 1 : 0; # passed threshold or not (boolean)
	if ($acc2ts->{$acc}) { # only collect things when this threshold is defined
	    # if there is a sequence threshold to look at, we need to collect the sum of normalized scores...
	    $acc2scoreSeq{$acc} += $h->{scoreNorm};
	}
    }

    # let's look at sequence threshold per family,
    my %acc2passSeq; # booleans that say if family passed or not (really a set, since only values of 1 need to be stored)
    while (my ($acc, $scoreSeq) = each %acc2scoreSeq) {
	# the sum of scores needs to be compared to the difference Ts-Td (that's the approximation we derived)
	$acc2passSeq{$acc} = 1 if $scoreSeq >= $acc2ts->{$acc};
    }
    
    # navigate domains again, and relabel if it had passed domain threshold but not sequence...
    foreach my $h (@$hits) {
	next unless $h->{GA}; # only need to re-evaluate things that passed domain threshold
	my $acc = $h->{acc};
	next unless $acc2ts->{$acc}; # also only need to re-evaluate things that needed sequence thresholds
	$h->{GA} = 0 unless $acc2passSeq{$acc}; # mark if sequence threshold wasn't passed! (stays 1 otherwise)
    }

    # this complicated labeling is finally complete!!!
    # nothing to return, since hits were edited by reference
}

sub mergeHashes {
    # this particular sub takes the second hash and transfers all the entries into the first hash, overwritting any repeated keys.  It does so in a memory efficient way (one key at the time), and doesn't need to return anything since it only takes in hash references
    # for the moment I don't need a fancier method, but many come to mind.  This is not fancy at all.
    my ($hash1, $hash2) = @_;
    while (my ($key, $val) = each %$hash2) {
	$hash1->{$key} = $val;
    }
}

1;

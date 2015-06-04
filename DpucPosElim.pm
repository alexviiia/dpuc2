# Copyright (C) 2014 Alejandro Ochoa, Singh Research Group, Princeton University
# This file is part of dPUC.
# dPUC is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
# dPUC is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with dPUC.  If not, see <http://www.gnu.org/licenses/>.

package DpucPosElim;
our $VERSION = 2.00;
use strict;

use Inline C => 'DATA';

1;

__DATA__
__C__

// global variables, which live while dPUC is needed (across proteins).  These two arrays map pfam accession pairs to the positive scores
// so they shouldn't be freed until dPUC is no longer needed (which is usually until the end, so meh).  These aren't considered memory leaks as far as I'm concerned.
int numNetTs = 0;
int* netTs = NULL;
int* netTscores = NULL;

// this holds the temporary Ts and Tscores arrays used per protein
// it's good practice to free "ts" and "tscores" between proteins, for cleanliness (I end up doing it before a new protein is started).  Since the same structures are reused (i.e. pointers not reset unless freed), this also doesn't constitute but a minor memory leak (gets worse if very large problems are considered towards the middle).
int numTsMax = 0; // this should signal us when it's time to grow.  Set to 0 so nothing is actually allocated if posElim isn't used (i.e. nonOvNeg), realloc should behave like malloc the first time around anyway (since pointers below are set to NULL initially too)
int stepGrowTs = 1024*1024; // this may be tuned somehow, not sure how much it will matter.  An array of 1024*1024 ints will take up 4 MB (cause each int takes 4 bytes), this isn't much at all!!!
int numTs = 0; // this tells us how big the data in the top two arrays is, for the purposes of loops, but isn't the total memory allocated, that's the number above
int* ts = NULL; // it's important to initialize these as NULL cause they are called with realloc only.  Nothing else uses realloc, so meh
int* tscores = NULL;
int negSelf = 0; // DPUC 2.0's directed network has a new series of priors that in some cases allow self interactions to be positive even in the absense of counts.  We need this value for posElim!

// copy this value from Perl space to C space
void setNegSelf (int negSelfIn) {
  negSelf = negSelfIn;
}

// these are necessary to avoid memory leaks (if process does Context in the beginning and decides to stop and do other things later)!!!  Need to call from Perl when done to free that memory!!!
void freeNetTscoresInC () {
  // don't do anything unless they had already been defined.  Better for nonOvNeg
  if (numNetTs != 0) {
    numNetTs = 0;
    free(netTs);
    free(netTscores);
    netTs = NULL;
    netTscores = NULL;
  }
}
// this one should be run after each protein is done!
void freeTscoresInC () {
  if (numTsMax != 0) {
    numTs = 0;
    numTsMax = 0;
    free(ts);
    free(tscores);
    ts = NULL;
    tscores = NULL;
  }
}
// shortcut to free all the structures, to finish program or bail out for errors
// NOTE: currently unused!
void freeAllInC () {
  freeNetTscoresInC();
  freeTscoresInC();
}

void growTscores () {
  // update size
  numTsMax += stepGrowTs;
  // for ts, attempt to grow the memory into a temporary pointer
  void *tmp;
  tmp = realloc(ts, numTsMax * sizeof(int));
  if (tmp != NULL) { ts = tmp; } // success, can overwrite old pointer
  else {
    // die with a message
    croak("Error: tried to grow ts array but realloc failed!\n");
  }
  // now try again for tscores
  tmp = realloc(tscores, numTsMax * sizeof(int));
  if (tmp != NULL) { tscores = tmp; } // success, can overwrite old pointer
  else {
    // die with a message
    croak("Error: tried to grow tscores array but realloc failed!\n");
  }
  // success!  nothing to return
}

// this sub, when called from perl, transfers the parallel arrays into native C arrays, and other methods access these directly so no further perl2c conversions are needed
void copyNetTscoresToC (AV* netTs_av, AV* netTscores_av) {
  // clean up previous scores if they had been defined!
  freeNetTscoresInC();
  // get length of arrays, will be the same for netTs and netTscores, save as global
  numNetTs = av_len(netTs_av) + 1;
  // initialize native versions, now that we know their size
  netTs = malloc(numNetTs * sizeof(int));
  netTscores = malloc(numNetTs * sizeof(int));
  // copy old values over
  int ti; // iterator
  for (ti = 0; ti < numNetTs; ti++) {
    // do not check for null values, transfer right away
    netTs[ti] = (int)SvIV(* av_fetch(netTs_av, ti, 0));
    netTscores[ti] = (int)SvIV(* av_fetch(netTscores_av, ti, 0));
  }
}

// this function zeroes the score with index t (preflatened, t(i,j)), which can be easily called from perl!!!
void zeroScore (int t) {
  // continue if there are ts/tscores to search through.  Searching is otherwise a waste of time and we'll have to error check the output to make sure we don't cause a segfault!!!
  if (numTs != 0) {
    // first find the index in the ts array
    int ti = binarySearchNumMin(t, ts, 0, numTs-1);
    // zero the entry if it was found
    // have to check that we actually found the entry, otherwise my func returns the least upper bound.
    if (ts[ti] == t) { 
      tscores[ti] = 0;
    }
  }
}

// this is the posElim version (uses positive scores only)
void getScoresNet_posElim (AV* hitis2acc_av) {
  // this sub turns the scores net that maps from Pfam accessions, into a half matrix that references each pair score by the domain pair index (from the input @$hits array).
  // I decided to separate the overlap net into a far more efficient process, once the @hitsPass list (in the outside of this sub) is short enough that it makes sense.
  // this used to store a full matrix, but now we only store half, with the largest index going first!  Diagonal is always zero so it isn't stored either!
  // factoring out variable allocation from the loops makes a difference in the edge of the high-performance spectrum...
  int i, j, acci, accj, netT, netTi;
    
  // it's always faster to translate perl arrays into native C arrays before using them in n^2 loops like the one below
  int numHitis = av_len(hitis2acc_av)+1;
  int hitis2acc[numHitis];
  // copy old values over to the new native array
  for (i = 0; i < numHitis; i++) { // I decided to reuse the iterator i instead of creating another variable!
    hitis2acc[i] = (int)SvIV(* av_fetch(hitis2acc_av, i, 0));
  }
  
  // we want to precalculate which hits are overlapping and shouldn't (doing it in every iteration takes a long time)
  // the structure we want is a parallel array, with the matrix indeces (i,j) i>j (flattened into a single number t) on the first array, and the corresponding scores on the other array.  The @ts are sorted naturally, since navigating the loop as we do, $t is always increasing so this is automatically satisfied!!!
  // we can halve the number of operations by not repeating pairs in loop!
  int t = 1; // starts at 1 cause we skipped i=j=0 (since k counts the diagonal too but here we don't)
  freeTscoresInC(); // to prevent side-effects or bugs, erase previous structure (we're supposed to overwrite the parts we want, but who knows if there are bugs, this is safer (is it slower?).  It will also release more memory after a particularly big problem has been found (it'd be kept till the end otherwise).  Note that this also blanks numTs if needed.
  //  numTs = 0; // blank this count, since we're starting all over
  // sanity check
  if (numTs != 0) {
    // die with a message
    croak("Error: presumably numTsMax is zero but numTs isn't!\n");
  }
  for (i = 1; i < numHitis; i++) { // start at i=1 cause j<i, so i=0 doesn't have a j!
    acci = hitis2acc[i];
    // I decided not to store some dumb intermediate variables (map $i to the hit structure, then to it's accession, then copy down structure)
    for (j = 0; j < i; j++) {
      accj = hitis2acc[j];
      // flatten pair indeces into a single number
      // DPUC < 2 used to encode netT pairs as ij2t, which ignores directionality, but DPUC 2.0 uses the ij2m encoding instead
      // ij2m function reimplemented from General.pm
      // note acci and accj are reversed compared to General.pm because they're navigated backwards (implicit ij2t map needs j<i but ij2m in General.pm needs acci,accj to reflect ordering in prot [i<j])
      if (accj >= acci) { netT = accj*accj + acci; }
      else { netT = (acci+1)*(acci+1) - accj - 1; }
      // find the score through binary search
      netTi = binarySearchNumMin(netT, netTs, 0, numNetTs-1); // get the index of netT in netTs
      if (netTs[netTi] == netT) { // have to do this check cause when the number wasn't found, the least upper bound is returned instead!
	// this means the index was found!!  Get score.
	// save only one way (i>j) via implicit ij2t map
	// note that t goes in the same order as we traverse i,j in this loop, so we don't have to calculate it, but we need to be careful to add 1 when we should have traversed the diagonal (cause we don't use the diagonal)
	// before adding, check that the current size can be saved, otherwise a segfault will occurr...
	if (numTs >= numTsMax) { growTscores(); } // we have to grow the array more...
	// now proceed, if we still have an error then something terrible happened!  But growTscores() should exit the program and cleanup if we had a problem anyway.
	ts[numTs] = t;
	tscores[numTs] = netTscores[netTi]; // actually get the score from the parallel array
	numTs++; // now we officially have a new element!
      }
      else if (negSelf > 0 && acci == accj) {
	// score wasn't found, but we still have a positive score if negSelf > 0 and acci == accj
	// go through the same drill as above...
	if (numTs >= numTsMax) { growTscores(); }
	ts[numTs] = t;
	tscores[numTs] = negSelf; // main difference here is using negSelf instead of the explicit score from the network
	numTs++;
      }
      // increment t for next j
      t++;
    }
    // increment t for the diagonal that we skip otherwise
    t++;
  }
}

// this new version uses ints instead of floats!!!  Proper scaling is done outside of this function
// for posElim only, nonOvNeg never uses a total sum like this (since we write each edge as a different variable)
SV* getSumScoreEdgesFromNode_posElim (AV* hitis_av) { // yey, the local typemap handles AV* !!!
  // returns an array that maps each hit index to its total context score (should be halved in some contexts since edges are shared)
  // note that, when accessing pair scores, we take advantage that, by construction, and maintained by the elimination, the hitis array is always sorted ascending!
    
  // factoring out variable allocation from the tight loops makes a difference in the edge of the high-performance spectrum...
  int k, l, i, ti, tri, t_want, s;
    
  // get length of incoming data
  int numHitis = av_len(hitis_av) + 1;
    
  // this is the mapping array we want, in terms of the indeces of the hitis indeces (so it is not a sparse array)
  // make native array first, and copy the values to an AV* when we are done
  int hitk2score[numHitis];
  // make sure array is initialized with zeroes, so the += code works
  for (k = 0; k < numHitis; k++) {
    hitk2score[k] = 0;
  }
    
  // it's odd, but sometimes the case that none of the hits have context.  I had a segfault problem originally, but regardless, it makes sense to halt extra processing and skip to returning the hitk2score with zeroes (perlified into an AV*) if we know there's no context in the problem!
  if (numTs != 0) {
    // create new array that has native ints
    int hitis[numHitis];
    // copy old values over to the new native array
    for (k = 0; k < numHitis; k++) {
      // do not check for null values, transfer right away
      hitis[k] = (int)SvIV(* av_fetch(hitis_av, k, 0));
    }
        
    // traverse pairs as in triangle matrix
    // note that in this loop ti=j+tri is going to increase monotonically, but not as ti++.  This should make it really easy to navigate ks[] in parallel, without binary search or anything similar
    ti = 0; // this stores the smallest index that is larger than the index of the last used score
    int t = ts[0]; // this remembers the t(i,j) that has not been used (so we save some array calls in comparisons below)
    for (k = 1; k < numHitis; k++) { // again, starts from k=1 cause l<k, so k=0 has no partner
      i = hitis[k]; // get i, the actual domain index, not the same as k (the index of the index, I know, confusing)
      tri = i*(i+1)/2; // the triangle version of i, to be used to generate t(i,j) efficiently
      for (l = 0; l < k; l++) {
	// get score of pair (i,j are already sorted correctly, j<i), then add it to both i's and j's total count
	// this is the index we want, note that j=hitis[l] but we do not store that variable
	t_want = hitis[l] + tri;
	// if this is the same as the next index in the scores array, we use that score and then increase the index for next time
	if (t == t_want) {
	  s = tscores[ti]; // this is the score we want to add
	  hitk2score[k] += s; // add score to native array to both entries
	  hitk2score[l] += s;
	  ti++; // now we have to compare to next entry, so update t index (important cause we use it to retrieve scores!)
	  t = ts[ti]; // also update the t for next round of comparisons
	}
	// if this is larger than the next index in the scores array (ti), we need to find a more appropriate ki
	// this may happen anywhere in the problem, although we expect it more for the beginning.  Not sure how often we expect it compared to equality
	// to be fast, we should use a binary search strategy, but we want it to return the t that either matches t_want, or the one next up in the list if nothing matched
	else if (t < t_want) {
	  // this returns the index ti of the least upper bound of t_want in ts, so that t >= t_want, and then we're in the same case as above
	  // note that we start looking with the old ti as the lower bound, not zero, hopefully it will make things a bit faster
	  // the upper bound is the last index of the matrix, can't be fancier there
	  ti = binarySearchNumMin(t_want, ts, ti, numTs-1);
	  t = ts[ti]; // update t with the right ti
	  // now we ask again if we got exactly the index we wanted or not, and repeat those steps accordingly
	  if (t == t_want) {
	    s = tscores[ti]; // this is the score we want to add
	    hitk2score[k] += s; // add score to native array to both entries
	    hitk2score[l] += s;
	    ti++; // now we have to compare to next entry, so update t index (important cause we use it to retrieve scores!)
	    t = ts[ti]; // also update the t for next round of comparisons
	  }
	  // else, the answer has to be forcibly t > t_want due to binary search step, meaning this pair did not have a score, and we have to wait for next round
	}
	// if this is smaller than the next index in the scores array, then this does not have a score, and we simply need to wait for next round without doing anything
      }
    }
  }
    
  // done with the math, now have to turn hitk2score into a perl AV*
  AV* hitk2score_av;
  hitk2score_av = newAV();
  for (k = 0; k < numHitis; k++) {
    av_store(hitk2score_av, k, newSViv(hitk2score[k]) );
  }
    
  return newRV_noinc((SV*)hitk2score_av); // typemap turns AV* into SV* automatically.  However, decided to call the reference creator myself so there aren't any reference count problems (no memory leaks)
}


int binarySearchNumMin (int x, int *a, int l, int u) {
  // search array of integers @$a for $x.
  // return index of match, or of the least upper bound otherwise
  // IMPORTANT: assumes list @$a is already sorted ascending, it doesn't work otherwise!!!
  // we don't sort here cause it would run the sort at every single call, but we can probably do it much less often than that...
  // based on code stolen from
  // http://staff.washington.edu/jon/dsa-perl/bsearch
  // actually this c-code was altered enough that now it looks more like the "Single comparison per iteration" in this article
  // http://en.wikipedia.org/wiki/Binary_search
  // typically l = 0, u = length -1, but they are really lower, upper end of search interval, so they can be used more wisely sometimes
    
  // bail out in trivial cases...
  // if u==-1, return -1.  This is common if u=length-1 and length==0, as computed outside...  No sound solution is possible...
  if (u == -1) { return -1; }
    
  // general code
  int i; // index of probe
  while (l < u) {
    i = l+(u-l)/2; // choose the index right at the middle of our range, c rounds automatically.  This version (instead of [u+l]/2) is more numerically stable.
    // note that since we want the least upper bound, upper has to include the probe i always when it needs to be updated!
    if (a[i] < x) { l = i+1; } // update bounds accordingly (next line too), get number in array only once!
    else { u = i; }
  }
  return u; // either returns the match index, or the index of the least upper bound!
}

# Copyright (C) 2014 Alejandro Ochoa, Singh Research Group, Princeton University
# This file is part of dPUC.
# dPUC is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
# dPUC is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with dPUC.  If not, see <http://www.gnu.org/licenses/>.

package DpucLpSolve;
our $VERSION = 1.00;
use strict;

use Inline C => 'Config',
#    'ENABLE' => 'AUTOWRAP',
    'INC' => '-I/usr/include/lpsolve', # needed to find lp_lib.h for our C code
    'LIBS' => "-llpsolve55 ",
    # paths that work on cetus
#    'INC' => '-I/Genomics/grid/users/ochoa/bin/lp_solve', # needed to find lp_lib.h for our C code
#    'LD' => 'gcc -Wl,-rpath=/Genomics/grid/users/ochoa/bin/lp_solve/lpsolve55/bin/ux64',
#    'LIBS' => "-L/Genomics/grid/users/ochoa/bin/lp_solve/lpsolve55/bin/ux64 -llpsolve55",
    ;
use Inline C => 'DATA';

1;

__DATA__
__C__

// conflicting definitions between lp_solve and perl :(

// /usr/lib64/perl5/CORE/handy.h
//#define TRUE (1)
//#define FALSE (0)
// /usr/lib64/perl5/CORE/perl.h
//#define STATIC static
// /usr/lib64/perl5/CORE/pp.h
//#define NORMAL PL_op->op_next

// /usr/include/lpsolve/lp_types.h
//#define FALSE                    0
//#define TRUE                     1
//#if 0
//  #define STATIC static
//#else
//  #define STATIC
//#endif
// /usr/include/lpsolve/lp_lib.h
//#define NORMAL                   4

// strategy is to undefine everything conflicting before entering the lp_solve code
#undef TRUE
#undef FALSE
#undef STATIC
#undef NORMAL
// lp_solve then re-defines these four things...
#include "lp_lib.h"
// we probably don't have to reset true/false, they're the same anyway, but we'll do it anyway
// we definitely have to reset static and normal for perl xs code to be happy
#undef TRUE
#undef FALSE
#undef STATIC
#undef NORMAL
#define TRUE (1)
#define FALSE (0)
#define STATIC static
#define NORMAL PL_op->op_next

void sendToLpSolve (int n, int numCols, AV* edgeDefs, AV* ovs, AV* SiDefs, AV* seqThreshs, AV* seqDefs, long timeout, AV* varNames, int verbose, int doSolve, char* foLp) {
  lprec *lp = NULL;
  int* colno = NULL; // define sparse rows this way
  REAL* row = NULL;
  char* lpStatus = NULL; // return value, should be set by lp_solve if successfull, or by us if there's a malloc error
  int i,j; // iterators
  REAL lpTime = 0.0; // time lp_solve ran for, init to zero when errors occur
  int lpNonZeroes = 0; // another potentially useful stat
  REAL score = 0.0; // value of objective function (i.e. total score of solution)
  int lpNodesProc = 0; // B&B nodes processed
  long long iter = 0; // B&B iterations + LP relaxations
  AV* sol = newAV(); // initialize solution as empty perl array
  
  // we need the number of constraints (ignoring objective function)
  // SiDefs = n
  // seq stuff = 3*@seqDefs, actually it's more complicated now...
  // edgeDefs are also 3 statements per edge
  int numEdgeDefs = av_len(edgeDefs)+1;
  int numOvs = av_len(ovs)+1;
  int numSeqDefs = av_len(seqDefs)+1;

  // counting the number of seqDefs is now a bit trickier, but a quick loop determines that size...
  // this only counts the strong seqDefs, there's one more weak statement that has its own loop (and there's exactly numSeqDefs of those)
  int numSeqDefsStrong = 0;
  for (i = 0; i < numSeqDefs; i++) {
    // this weird thing gets a single definition, a perl array with a bunch of integers
    AV* seqDef = (AV*)SvRV(* av_fetch(seqDefs, i, 0));
    // this is how many statements this particular seqDef has (total-1, remember av_len already includes this -1)
    numSeqDefsStrong += av_len(seqDef);
  }
  
  // need length of input data to compute numRows...
  int numRows = n + 3*numEdgeDefs + numOvs + 2*numSeqDefs + numSeqDefsStrong;
  
  // make colno/row large enough to fit an entire row
  colno = malloc(numCols * sizeof(int));
  row = malloc(numCols * sizeof(REAL));
  // also create a new LP model
  lp = make_lp(numRows, numCols);
  // make sure we don't continue (except to cleanup) unless these structures were allocated correctly
  if ((colno == NULL) || (row == NULL) || (lp == NULL)) {
    fprintf(stderr, "Unable to create new LP model: Malloc error\n"); // warn user!
    lpStatus = 'MALLOCFAIL'; // return value if we had problems
  }
  else {
    if (verbose || (foLp != NULL && *foLp != '\0')) {
      // set names for each variable, useful for debugging
      for (i = 1; i <= numCols; i++) {
	set_col_name(lp, i, (char*)SvPV_nolen(* av_fetch(varNames, i, 0)));
      }
    }
    
    // set binary types for certain columns
    // 1..n are Di's, all binary
    for (i = 1; i <= n; i++) {
      set_binary(lp, i, TRUE);
    }
    // n+1..2n are Si's, NOT binary
    // 2n+1 till end are Eij's and acc's, again all binary
    for (i = 2*n+1; i <= numCols; i++) {
      set_binary(lp, i, TRUE);
    }
    // note Si's default type and range are perfect
    // in particular, Si >= 0 by default, which is the domain threshold we wish to enforce!  No need to write it out separately
  
    /* Model created */
    // this makes model building slightly faster
    set_add_rowmode(lp, TRUE);
  
    /* first! set the objective function */
    for (i = 0; i < n; i++) { // there are n Si's
      colno[i] = i+1+n; // the 1-based index of Si
      row[i] = 1.0; // the coefficient is always 1
    }
    set_obj_fnex(lp, n, row, colno); // actually set objective function in *lp, only look at the first n entries of row/colno!
  
    /* now add the constraints */
    int offset = 1; // initial index offset is just +1 to be 1-based
  
    // define Si's in terms of Di's and Eij's
    for (i = 0; i < n; i++) {
      // this weird thing gets a single definition, a perl array with two other arrays that define the sparse row
      AV* SiDef = (AV*)SvRV(* av_fetch(SiDefs, i, 0));
      AV* SiDefColno = (AV*)SvRV(* av_fetch(SiDef, 0, 0));
      AV* SiDefRow = (AV*)SvRV(* av_fetch(SiDef, 1, 0));
      // get length of sparse row to iterate over it
      int numSiDef = av_len(SiDefColno)+1;
      // transfer numbers from perl to C
      for (j = 0; j < numSiDef; j++) {
	colno[j] = (int)SvIV(* av_fetch(SiDefColno, j, 0));
	row[j] = (REAL)SvNV(* av_fetch(SiDefRow, j, 0));
      }
      set_rowex(lp, offset+i, numSiDef, row, colno); // add constraint!
      set_constr_type(lp, offset+i, EQ); // this is an "equal" constraint
      set_rh(lp, offset+i, 0.0); // and it's equal to zero
    }
    offset += n; // now offset also includes the previously defined Si definitions
    
    // do overlap constraints
    for (i = 0; i < numOvs; i++) {
      // this weird thing gets an "overlap", a perl array with two or more ints
      AV* ov = (AV*)SvRV(* av_fetch(ovs, i, 0));
      int numOv = av_len(ov)+1;
      for (j = 0; j < numOv; j++) {
	colno[j] = (int)SvIV(* av_fetch(ov, j, 0));
	row[j] = 1.0; // since we can have more than two OVs, we need to make sure the 1.0 coefficients extend as far as we need them
      }
      set_rowex(lp, offset+i, numOv, row, colno); // add constraint!
      // default constraint type is LE, won't change that
      set_rh(lp, offset+i, 1.0); // upper bound is 1
    }
    offset += numOvs; // now offset also includes the previously defined overlap constraints
    
    // edge definitions, part 1 (most restrictive)
    // the two coefficients are the same for these two statements, only the indexes change
    row[0] = -1.0; row[1] = 1.0;
    for (i = 0; i < numEdgeDefs; i++) {
      // this weird thing gets an "edge", a perl array with three ints
      AV* edge = (AV*)SvRV(* av_fetch(edgeDefs, i, 0));
      // this section is easier if we give the entries names
      colno[1] = (int)SvIV(* av_fetch(edge, 2, 0)); // set the index of the thing AND is assigned to
      // now create two statements...
      // since two constraints arise from this data, the output indexes take this alternating form "offset+2*i" and "offset+2*i+1"
      for (j = 0; j < 2; j++) {
	colno[0] = (int)SvIV(* av_fetch(edge, j, 0));
	set_rowex(lp, offset+2*i+j, 2, row, colno); // add constraint!
	// default constraint type is LE, won't change that
	set_rh(lp, offset+2*i+j, 0.0); // and it's <= zero
      }
    }
    offset += 2*numEdgeDefs; // now offset also includes the previously defined edge constraints
    
    // edge definitions, part 2 (least restrictive)
    // these coefficients are all the same for these statements, only the indexes change
    row[0] = 1.0; row[1] = 1.0; row[2] = -1.0;
    for (i = 0; i < numEdgeDefs; i++) {
      // this weird thing gets an "edge", a perl array with three ints
      AV* edge = (AV*)SvRV(* av_fetch(edgeDefs, i, 0));
      for (j = 0; j < 3; j++) {
	colno[j] = (int)SvIV(* av_fetch(edge, j, 0));
      }
      set_rowex(lp, offset+i, 3, row, colno); // add constraint!
      // default constraint type is LE, won't change that
      set_rh(lp, offset+i, 1.0); // upper bound is 1
    }
    offset += numEdgeDefs; // now offset also includes the previously defined edge constraints
    
    // now enforce sequence thresholds, code is very similar to above because colno/row are both fully precalculated
    for (i = 0; i < numSeqDefs; i++) {
      // this weird thing gets a single definition, a perl array with two other arrays that define the sparse row
      AV* seqThresh = (AV*)SvRV(* av_fetch(seqThreshs, i, 0));
      AV* seqThreshColno = (AV*)SvRV(* av_fetch(seqThresh, 0, 0));
      AV* seqThreshRow = (AV*)SvRV(* av_fetch(seqThresh, 1, 0));
      // get length of sparse row to iterate over it
      int numSeqThresh = av_len(seqThreshColno)+1;
      // transfer numbers from perl to C
      for (j = 0; j < numSeqThresh; j++) {
	colno[j] = (int)SvIV(* av_fetch(seqThreshColno, j, 0));
	row[j] = (REAL)SvNV(* av_fetch(seqThreshRow, j, 0));
      }
      set_rowex(lp, offset+i, numSeqThresh, row, colno); // add constraint!
      set_constr_type(lp, offset+i, GE); // this is an "greater or equal" constraint
      set_rh(lp, offset+i, 0.0); // and it's >= zero
    }
    offset += numSeqDefs; // now offset also includes the previously defined sequence thresholds
  
    // last things left are the sequence definitions (in terms of the Di's)
    // this portion sets all the strong constraints (numSeqDef-1 per seqDef)
    // the two coefficients are the same for these statements, only the indexes change
    row[0] = 1.0; row[1] = -1.0;
    for (i = 0; i < numSeqDefs; i++) {
      // this weird thing gets a single definition, a perl array with a bunch of integers
      AV* seqDef = (AV*)SvRV(* av_fetch(seqDefs, i, 0));
      // get length of data to iterate over it
      int numSeqDef = av_len(seqDef)+1;
      // the first integer is special (sequence index), same for all constraints (put in second position)
      colno[1] = (int)SvIV(* av_fetch(seqDef, 0, 0)); // index of family/seq variable
      // navigate all domains now, put their indexes in the first position
      for (j = 1; j < numSeqDef; j++) { // start from non-special entries
	colno[0] = (int)SvIV(* av_fetch(seqDef, j, 0));
	// can't use offset+i*k+l because numSeqDef is different for each seqDef (not a rectangular matrix), but the formula below should be accurate enough
	set_rowex(lp, offset+j-1, 2, row, colno); // add constraint!
	// default constraint type is LE, won't change that
	set_rh(lp, offset+j-1, 0.0); // and it's <= zero
      }
      offset += numSeqDef-1; // now offset also includes the previously defined strong sequence definition for a single seqDef
    }
    
    // last things left are the sequence definitions (in terms of the Di's)
    // this portion sets the weak constraints (one per seqDef)
    for (i = 0; i < numSeqDefs; i++) {
      // this weird thing gets a single definition, a perl array with a bunch of integers
      AV* seqDef = (AV*)SvRV(* av_fetch(seqDefs, i, 0));
      // get length of data to iterate over it
      int numSeqDef = av_len(seqDef)+1;
      // transfer numbers from perl to C
      for (j = 0; j < numSeqDef; j++) {
	colno[j] = (int)SvIV(* av_fetch(seqDef, j, 0));
	row[j] = -1.0; // all coefficients are -1 (except for first one, which we set below)
      }
      row[0] = 1.0; // set first element for first constraint to what it should be
      set_rowex(lp, offset+i, numSeqDef, row, colno); // add constraint!
      // default constraint type is LE, won't change that
      set_rh(lp, offset+i, 0.0); // and it's <= zero
    }
    offset += numSeqDefs; // now offset also includes the previously defined weak sequence definitions
  
    // finally done constructing model
    // need to reset this before solving
    set_add_rowmode(lp, FALSE);
  
    // minimizing is the default optimization, change that
    set_maxim(lp);
    
    // since we have a sum of Si's and each Si > 0, we can safely say solution is > 0, not sure if this will make things easier for B&B algorithm
    //set_obj_bound(lp, 0.0);
    
    // messing with B&B settings
    set_bb_floorfirst(lp, BRANCH_FLOOR);
    //set_bb_floorfirst(lp, BRANCH_CEILING);
    set_bb_rule(lp, NODE_AUTOORDER);
    
    //set_scaling(lp, SCALE_GEOMETRIC|SCALE_EQUILIBRATE); // S,S , yey, maybe all the "unbounded" problems have to do with int scaling??????? nodes/iter = (2.74, 2.56)
    set_scaling(lp, SCALE_RANGE|SCALE_EQUILIBRATE);      // S,S, same score/time but better/same nodes/iter (2.78, 2.55)
    
    // set timeout!
    set_timeout(lp, timeout);
    
    // write LP model to output file if provided
    if (foLp != NULL && *foLp != '\0')
      write_lp(lp, foLp);
    
    // use this line if you want to see the model created!
    if (verbose)
      write_LP(lp, stdout);
    
    // set verbosity, let's almost mute unless we're debugging
    if (!verbose)
      set_verbose(lp, CRITICAL); // SEVERE
    //else
    //set_verbose(lp, FULL); // extreme debugging
    
    
    // get size of problem (non-zero elements of matrix)
    lpNonZeroes = get_nonzeros(lp);
    
    // SOLVE!
    if (doSolve) {
      int ret = solve(lp);
      // see how long lp_solve thinks this solution took to solve, to compare to wall time
      lpTime = time_elapsed(lp);
      lpStatus = get_statustext(lp, ret); // get text version of "ret", which is a lot more useful to us humans
      score = get_objective(lp); // get value of objective function
      lpNodesProc = get_total_nodes(lp); // number of int nodes processed for B&B
      iter = get_total_iter(lp); // get number of iterations for LP relaxations + B&B
      // we have solutions to look at if ret is 0 or 1
      if ((ret == 0) || (ret == 1)) {
	// read solution
	get_variables(lp, row);
	// add to perl array to return
	for (i = 0; i < numCols; i++) {
	  av_push(sol, newSVnv(row[i]));
	}
	// print solution if we want that much detail...
	if (verbose) {
	  for (i = 0; i < numCols; i++) {
	    printf("%s: %f\n", get_col_name(lp, i+1), row[i]);
	  }
	}
      }
    }
    else {
      lpStatus = "OPTIMAL"; // fool loops into thinking this worked! (would be MALLOCFAIL otherwise)
    }
  }
  
  // since we return two things into perl, we need to deal with this stack...
  Inline_Stack_Vars;
  Inline_Stack_Reset;
  Inline_Stack_Push(sv_2mortal(newRV_noinc((SV*) sol))); // add the solution array
  Inline_Stack_Push(sv_2mortal(newSVpv(lpStatus, 0))); // add string return value
  Inline_Stack_Push(sv_2mortal(newSVnv(lpTime))); // add runtime
  Inline_Stack_Push(sv_2mortal(newSViv(numRows))); // add number of rows of problem
  Inline_Stack_Push(sv_2mortal(newSViv(lpNonZeroes))); // add size of problem
  Inline_Stack_Push(sv_2mortal(newSVnv(score))); // add value of objective function
  Inline_Stack_Push(sv_2mortal(newSViv(lpNodesProc))); // add number of int nodes processed by B&B
  Inline_Stack_Push(sv_2mortal(newSViv(iter))); // add number B&B + LP relaxation iterations
  Inline_Stack_Done; // say we're done

  // clean up regardless of lpStatus
  if (row != NULL)
    free(row);
  if (colno != NULL)
    free(colno);
  if (lp != NULL)
    delete_lp(lp);
}


# Copyright (C) 2014 Alejandro Ochoa, Singh Research Group, Princeton University
# This file is part of DomStratStats and dPUC.
# DomStratStats and dPUC are free software: you can redistribute them and/or modify them under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
# DomStratStats and dPUC are distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with DomStratStats and dPUC.  If not, see <http://www.gnu.org/licenses/>.

package Hmmer3ScanTab;
our $VERSION = 1.03;
use lib '.';
use FileGz;
use strict;

# version 1.03 2014-11-06 delta
# - runHmmscan has new optional parameter pCut for final p-value threshold.  Default stays at 1e-4 for public code, but for benchmarks need it to be 1e-2 (for DSS paper), so now that option is available.  Old code should run just as before, although I haven't made any explicit tests!
# - not published yet, just internal for now

# Note: column where insertion happens is defined here!  A global constant!
my $numColsDef = 23;

# functions for parsing, editing, and printing HMMER3 domain tabular outputs from hmmscan

sub runHmmscan { 
    my ($hmmscan, $pfamA, $fiSeq, $fo, $pCut) = @_;
    
    # most interesting processing is treatment of heuristic filters (dPUC doesn't get enough with the defaults).  Below the defaults are shown for reference.  We should change them if we aim to be more permissive!
    #  --F1 <x> : MSV threshold: promote hits w/ P <= F1  [0.02]
    #  --F2 <x> : Vit threshold: promote hits w/ P <= F2  [1e-3]
    #  --F3 <x> : Fwd threshold: promote hits w/ P <= F3  [1e-5]
    
    # parameters
    $pCut = 1e-4 unless defined $pCut; # (F3) a separate analysis showed that this covers reasonable q-values and local FDRs with no problem.  Let's not compute and/or store more preds than needed!
    my $pCut2 = 1e-1; # F1 and F2 more permissive thresholds for HMMER3 heuristic filters
    
    # I won't validate input, the plan is just to let hmmscan complain
    
    $fiSeq = FileGz::realFile($fiSeq); # input may be compressed but specified without .gz extension, this will make it work transparently!
    
    # actually run hmmscan
    my @system = ($hmmscan, '--F1', $pCut2, '--F2', $pCut2, '--F3', $pCut, '-Z', 1, '--domZ', 1, '-E', $pCut, '--domE', $pCut, '-o', '/dev/null', '--domtblout', $fo, $pfamA, $fiSeq);
    print "Command running:\n> @system\n";
    my $ex = system @system;
    die "Error: the previous command returned $ex!\n" if $ex;
}

sub parseInit {
    # this is a version specifically optimized for the hmmscan version, which processes one protein at the time.  In that case, there's no memory issues, so there's need for protein filtering or compression or any column loss.  Assumes an input filehandle is the only case accepted too.
    # USAGE: my ($fhi, $header, $protNext, $hitNext) = Hmmer3ScanTab::parseInit($fi);
    my ($fi, $extra) = @_;
    my $header = ''; # stores lines of header
    my $nextLine; # input to next protein's processing
    # start reading file
    my $fhi = FileGz::getInFh($fi);
    while (<$fhi>){
	if (/^\#/) { $header .= $_; } # store header, will regurgitate onto dpuc2's output
	else { $nextLine = $_; last; }
    }
    my ($protNext, $hitNext) = _parseOneLine($nextLine, $extra); # just parse this thing for next round
    return ($fhi, $header, $protNext, $hitNext);
}

sub parseProt {
    # this is a version specifically optimized for the hmmscan version, which processes one protein at the time.  In that case, there's no memory issues, so there's need for protein filtering or compression or any column loss.  Assumes an input filehandle is the only case accepted too.
    ### NOTE: hits might be unsorted, dpuc2 doesn't care
    # USAGE: (my $prot, my $hits, $protNext, $hitNext) = Hmmer3ScanTab::parseProt($fhi, $protNext, $hitNext);
    my ($fhi, $prot, $hitPrev, $extra) = @_;
    my @hits = ($hitPrev); # this is the main structure we want to return, add last previously-parse line
    # start reading file
    my ($protNext, $hitNext); # to return or next round
    while (<$fhi>) {
	next if /^\#/; # skip commented lines, more proteins might follow (from concatenated HMMER3 outputs) or it might just be the end of file, but let's not assume either way.
	($protNext, $hitNext) = _parseOneLine($_, $extra);
	last if $protNext ne $prot; # stop processing when protein changes
	# if we're still in the same protein, add this hit to the list of hits
	push @hits, $hitNext;
	$hitNext = ''; # also clear these two things so that next time we can distinguish a new protein from the end of the file
	$protNext = '';
    }
    return ($prot, \@hits, $protNext, $hitNext); # done, return hits and next line's processed data (or empty otherwise)
}

# internal function, not needed by end users
sub _parseOneLine {
    # this keeps the entire original line from the HMMER3 input.  It is memory-inefficient but we assume this is used for hmmscan so only one protein at the time is analyzed, so we never load a whole proteome (or worse, uniref) using this!!!
    my ($line, $extra) = @_; # only input is line to be parsed, and a string that specifies if we should have additional columns or not (for our newest DomStratStats outputs)
    my $origLine = $line; # copy a version that is unedited, for final output (includes newline)
    chomp $line;
    my @fields = split / +/, $line; # split with spaces, will work for everything except the last description, which does have spaces, but we don't actually care about it
    my $prot = $fields[3]; # get protein's name
    # where we store the parsed data, let's go minimal since everything is preserved in "origLine"
    # start coordinates are always reduced by one so we have a 0-based, end exclusive range system
    my %hit = (
	'origLine' => $origLine,
	'acc' => $fields[1],
	'start' => $fields[17] - 1,
	'end' => $fields[18],
	'name' => $fields[0],
	'start2' => $fields[15] - 1,
	'end2' => $fields[16],
	'score' => $fields[13],
	'E' => $fields[12],
	'scoreSeq' => $fields[7],
	'ESeq' => $fields[6],
	);
    $hit{acc} =~ s/\.\d+$//; # remove version number for accession (shouldn't affect non-Pfams)
    # add extras as needed
    if ($extra) { # in rawest data this is undefined, so get that out of the way
	if ($extra eq 'domStratStats') {
	    $hit{q} = $fields[22];
	    $hit{lFDR} = $fields[23];
	    $hit{FDRl} = $fields[24];
	}
	elsif ($extra eq 'tieredStratStats') {
	    $hit{qSeq} = $fields[22];
	    $hit{qDomCond} = $fields[23];
	}
	else { die "Error: Hmmer3ScanTab::_parseOneLine had extra='$extra' unmatched"; } # if this happens, it is probably a bug!
    }
    return ($prot, \%hit); # this is all we need
}

sub parseAllCut {
    # unlike general functions that load one prot at the time to be memory efficient, this loads it all into memory (it is more practical sometimes), and sets arbitrary significance thresholds if needed
    # undefined threshold $cut will set no filter
    # so for now assumes that $cut is maximum value to keep (but in the future could provide a variable that indicates if we want inequality reversed)
    # note that undefined $key behaves like $key='E' because of eCuts()
    my ($fi, $cut, $key, $extra) = @_;
    
    my %prot2hits; # the structure we want
    
    # now read the full files entry-by-entry
    my ($fhi, $header, $protNext, $hitNext) = parseInit($fi, $extra);
    while ($protNext) {
	# read next protein
	(my $prot, my $hits, $protNext, $hitNext) = parseProt($fhi, $protNext, $hitNext, $extra);
	# set threshold and overwrite
	$hits = eCuts($hits, $cut, $key) if defined $cut; # if there was no threshold provided, do not try to filter!
	# store if there were hits (always an array ref, so this won't complain)
	$prot2hits{$prot} = $hits if @$hits;
    }
    
    # done with input file
    close $fhi;
    # return structure we wanted
    return \%prot2hits;
}

sub eCuts {
    # simple implementation of E-value thresholds.
    # This doesn't remove overlaps.
    my ($hits, $Ecut, $key) = @_;
    $key = 'E' unless defined $key; # E-values are maintained as default for backward-compatibility, but can set thresholds on q instead, or other options.
    my @hitsPass; # new hits we want
    foreach my $h (@$hits) {
	# apply threshold to the "domain" E-value
	push @hitsPass, $h if $h->{$key} <= $Ecut; 
    }
    return \@hitsPass;
}

sub parsePValues {
    # this function focuses in parsing the p-value distributions, so most things aren't parsed, to be more efficient
    # have to keep track of p-value counts per HMM (by name), of additional counts for m adjustment, and of sequence p-values too (for tiered analysis only)
    # so there are two forms
    # 1) the simplest is the domain-only analysis which only needs domain p-values and the m adjustment
    # 2) the tiered analysis doesn't need the domain p-values directly or their m adjustment, but it does need sequence p-values and domain p-values tiered by sequence
    # I tried to optimize the code so we don't gather data we don't need, which varies depending on the mode.  This might be faster, but I'm particularly worried about memory usage in large files.  We'd never use both modes simultaneously, they're similar but incompatible.
    # note a list of files (an array ref) is also accepted, and will return p-values for all inputs combined!
    # NOTE: used to use accessions, but names are equally unique in Pfam and they are better for Supfam (Supfam accs are not unique)
    # USAGE DOMAIN: my ($name2pDom2c, $name2mDom) = ParseHmmer3ScanTab::parsePValues($fis);
    # USAGE TIERED: my ($name2pSeq2c, $name2pSeq2pDom2c) = ParseHmmer3ScanTab::parsePValues($fis, $doTiered);
    my ($fis, $doTiered) = @_;
    # the data we want (combined across files!), but depending on mode some of it will not be gathered!
    my %name2mDom; # dom only
    my %name2pDom2c; # dom only
    my %name2pSeq2c; # tiered only
    my %name2pSeq2pDom2c; # tiered only
    # generalize to file lists
    my @fis = 'ARRAY' eq ref $fis ? @$fis : ($fis); # get file list, or make a list with the one file otherwise
    foreach my $fi (@fis) {
	# some temp vars, these store per-protein data
	my $protLast = ''; # input to next protein's processing
	my %name2c; # dom only, keep counts of names per prot, to adjust m per name
	my %name2pSeq; # tiered only, get unique pSeq maps per prot, each family should be counted (per prot) once at this level!
	# start reading file
	my $fhi = FileGz::getInFh($fi);
	while (<$fhi>){
	    next if /^\#/;
	    chomp;
	    my @fields = split / +/; # split with spaces, will work for everything except the last description, which does have spaces, but we don't actually care about it
	    my $prot = $fields[3]; # get protein's name
	    if ($prot ne $protLast) {
		# changed protein, have to do some things
		if ($doTiered) {
		    while (my ($name, $pSeq) = each %name2pSeq) {	
			$name2pSeq2c{$name}{$pSeq}++; # add a count for this sequence p-value and this name to master hash
		    }
		    %name2pSeq = (); # clear this hash for next round
		}
		else {
		    # store adjustments to m if necessary
		    while (my ($name, $c) = each %name2c) {
			$name2mDom{$name} += $c-1 if $c > 1; # these are the necessary adjustments (when $c > 1 only)
		    }
		    %name2c = (); # clear this hash for next round
		}
		$protLast = $prot; # update for next round
	    }
	    # after transfering things from last protein, now we process for the new current prot!
	    my $name = $fields[0]; # use name (for Supfam in particular!)
	    my $pDom = $fields[12]; # E
	    if ($doTiered) {
		my $pSeq = $fields[6]; # ESeq
		$name2pSeq2pDom2c{$name}{$pSeq}{$pDom}++; # domain p-value histogram stratified by sequence p-values!
		$name2pSeq{$name} = $pSeq; # store sequence p-value (overwrite ok cause it should be the same value everytime we see it)
	    }
	    else {
		$name2pDom2c{$name}{$pDom}++; # add count for this name-pDom pair
	    }
	}
	close $fhi; # close input!
	# when we are done, we haven't stored the last protein's data, let's do that now!
	if ($doTiered) {
	    while (my ($name, $pSeq) = each %name2pSeq) {	
		$name2pSeq2c{$name}{$pSeq}++; # add a count for this sequence p-value and this name to master hash
	    }
	}
	else {
	    # store adjustments to m if necessary
	    while (my ($name, $c) = each %name2c) {
		$name2mDom{$name} += $c-1 if $c > 1; # these are the necessary adjustments (when $c > 1 only)
	    }
	}
	# no need to clear temp vars anymore!
    }
    ### return data!!!  Behaves differently depending on mode!
    if ($doTiered) { return (\%name2pSeq2c, \%name2pSeq2pDom2c); }
    else { return (\%name2pDom2c, \%name2mDom); }
}

#my $offset = 180; # hardcoded, since HMMER3 format is fixed space.
#sub addCols {
#    # takes a domain's string and adds columns (as a single substring) to the original string by adding them just before the "description of target"
#    # have to handle some anomalous cases unfortunately!
#    my ($origLine, $newStr) = @_;
#    # let's handle the special case in which a line is not long enough (I think only actual case is first line of header, meh)
#    my $lengthOrigLine = length $origLine;
#    if ($lengthOrigLine < $offset) {
#	$origLine .= ' ' x (1+$offset - $lengthOrigLine); # this adds just enough whitespace, plus one more to make next offset2 procedure work for this case too
#    }
#    # now, we might have to shift offset if it isn't on a space character (this means something in the HMMER3 printing was larger than it should have been, I personally encountered it on protein PF3D7_1038400 versus domain NPR, it has a huge sequence bias value of 11330.7 when normal bias values look like 32.3 at the most
#    # in general, increase offsets until a space character is found (which is usually right away, at $offset2==0)
#    my $offset2 = 0; # do not edit global offset, this is a local offset that varies per line
#    while (' ' ne substr $origLine, $offset+$offset2, 1) { $offset2++; }
#    # we're finally ready to insert the string we want!
#    # note input gets padded with spaces, so input doesn't need to have that.  This edits input by reference, but it was copied on the way in so the outside value should be unedited, I think
#    substr $origLine, $offset+$offset2, 0, ' '.$newStr;
#    return $origLine; # this is the edited version now, return
#}

sub posOfCol {
    my ($str, $nCol) = @_;
    $nCol = $numColsDef unless defined $nCol; # use default column number to insert in!
    # finds the first character position of a given column in a string assumed to be separated by variable amounts of whitespace.
    # since this is only for HMMER3, whitespace is always single spaces, and not other whitespaces in general
    # this is necessary since long IDs, or scores that are too long, can shift columns, so we can't rely on fixed spacing.  However, we can rely on space-separated columns, so we just have to count those.  In my tests I verified that if input ID has spaces then HMMER3 only uses the first word as the ID, so there are no spaces where they shouldn't be!
    my $pos = -1; # tells me where the last whitespace position was...
    for (my $i=0; $i < $nCol-1; $i++) {
	# main procedure is iterating the search for spaces, this way we count the columns
	if ($str =~ / +/g) {
	    # moved ahead, keep going... but because of block-local annoyances, must record end of range now
	    $pos = $+[0];
	}
	else { # warn if we didn't find enough columns before the end!
	    die "Error in Hmmer3ScanTab::posOfCol(): only ".($i+1)." columns were found in sequence below but $nCol were desired!\n$str\n";
	}
    }
    return $pos;
}

sub addCols {
    my ($origLine, $newStr, $pos) = @_;
    # like addCols but here the position of insertion has been previously determined, so only the splicing happens
    # separated from posOfCol because for headers we only determine this position for the last of the three header lines, it's impossible to determine it from the first one in isolation
    
    # in this case we determine insertion position from input string, most common case
    # but we can't always do that so sometimes we provide a defined $pos that allows us to skip this case
    # use default column number (it has been ommited)!
    $pos = posOfCol($origLine) unless defined $pos;
    
    # let's handle the special case in which a line is not long enough (I think only actual case is first line of header)
    my $lengthOrigLine = length $origLine;
    if ($lengthOrigLine < $pos) {
	$origLine .= ' ' x ($pos - $lengthOrigLine); # this adds just enough whitespace
    }
    
    # note input gets padded with spaces, so input doesn't need to have that.  This edits input by reference, but it was copied on the way in so the outside value should be unedited, I think
    substr $origLine, $pos, 0, $newStr.' ';
    return $origLine; # this is the edited version now, return
}

sub addColsHeader {
    # wrapper for additional header processing
    my ($header, $newStr1, $newStr2, $newStr3) = @_;
    my $length = length $newStr1; # get length of first string
    die "Error: input lines have different lengths!" unless $length == length $newStr2 && $length == length $newStr3; # just make sure these are the same, or output will appear misaligned
    my @header = split /\n/, $header; # split to get separate lines
    # unfortunately only the third header line can be used to count columns, so let's do that first
    my $pos = posOfCol($header[2]); # use default column number, same as usual
    # generously use the standard column insertion function to add strings to each, with the exception that here the insertion position is provided (necessary for first two header lines)
    $header[0] = addCols($header[0], $newStr1, $pos);
    $header[1] = addCols($header[1], $newStr2, $pos);
    $header[2] = addCols($header[2], $newStr3, $pos);
    # done, just put lines back together to return a string
    return join "\n", @header, '';
}

sub printInit {
    # USAGE: my $fho = Hmmer3ScanTab::printInit($fo, $header, $comp);
    my ($fo, $header, $comp) = @_; # format is predetermined here, can't specify it because you change it
    my $fho = FileGz::getOutFh($fo, $comp);
    print $fho $header; # print the header that the HMMER3 parser provided (has newlines already)
    return $fho; # only the file handle is needed for subsequent processing (header isn't needed)
}

sub printProt {
    # super easy hit printer that assumes original HMMER3 line was kept (like my _parseOneLine does)
    # USAGE: Hmmer3ScanTab::printProt($fho, $hits);
    my ($fho, $hits) = @_;
    foreach my $h (@$hits) {
	print $fho $h->{origLine}; # this repeats original input in same format (nothing changed at all, includes newline)
    }
}

1;

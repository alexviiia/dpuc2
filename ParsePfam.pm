# Copyright (C) 2014 Alejandro Ochoa, Singh Research Group, Princeton University
# This file is part of DomStratStats and dPUC.
# DomStratStats and dPUC are free software: you can redistribute them and/or modify them under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
# DomStratStats and dPUC are distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with DomStratStats and dPUC.  If not, see <http://www.gnu.org/licenses/>.

package ParsePfam;
our $VERSION = 1.00;

use lib '.';
use FileGz;
use strict;

# define hash structures, let's make these globals cause they're used so often.  However, they're undefined until parsePfamDat is called!  Will be overwritten if we call parsePfamDat with different versions, so be careful how this is used!!!
our ($acc2name, $acc2de, $acc2ds2ga, $acc2type, $acc2length, $acc2clan, $nestingNet);
# same but for pfamC file
our ($clan2name, $clan2de, $clan2cc, $clan2scop, $clan2cath, $clan2merops);

sub dat {
    my ($fi) = @_; # input file
    # a parser of the useful data contained in a small DAT file.  Doesn't exist for Pfam 23 and older.
    # usage
    # my ($acc2name, $acc2de, $acc2ds2ga, $acc2length, $acc2clan, $acc2type, $nestingNet) = ParsePfam::dat($fiDat);
    
    # the actual new structures we want
    my (%acc2name, %acc2de, %acc2ds2ga, %acc2length, %acc2clan, %acc2type);
    my %acc2nesting; # maps accession to a list of accessions that are allowed to overlap (nesting).  This one is temporary for parsing, but we always want it in the form of a network instead, so that's returned instead.
    # a temporary variable that holds the current accession, so everything else can be mapped
    my $acc = '';
    my $name = ''; # name comes before accession, stupid, but we have to remember it and defer mapping of this one by one line...
    
    # open file, start reading
    my $fhi = FileGz::getInFh($fi);
    while (<$fhi>) {
	if ( my ($type, $s) = /^\#=GF (\w+)\s+(.+)$/ ) { # things we're parsing for... many lines are not interesting
	    if ($type eq 'ID') { $name = $s; }
	    elsif ($type eq 'AC') {
		$acc = $s; # read new accession
		$acc =~ s/\.\d+$//; # remove version digits
		$acc2name{$acc} = $name; # name gets mapped here
	    }
	    elsif ($type eq 'DE') { $acc2de{$acc} = $s; }
	    elsif ($type eq 'TP') { $acc2type{$acc} = $s; }
	    elsif ($type eq 'ML') { $acc2length{$acc} = $s; }
	    elsif ($type eq 'CL') { $acc2clan{$acc} = $s; }
	    elsif ($type eq 'NE') { push @{$acc2nesting{$acc}}, $s; }
	    elsif ($type eq 'GA') {
		# split string into the two values it contains, the sequence and the domain thresholds
		# store them directly in the hash we want via this complicated hash ref slice
		@{$acc2ds2ga{$acc}}{qw(s d)} = split /; */, $s;
	    }
	}
    }
    close $fhi;
    
    # annoying, nesting is in terms of domain names, but we want them in accessions
    # first step is to get name2acc mapping
    my %name2acc = reverse %acc2name; # this works as long as the names are unique!
    # now go through each value and map
    foreach my $acc (keys %acc2nesting) {
	@{$acc2nesting{$acc}} = map { $name2acc{$_} } @{$acc2nesting{$acc}};
    }
    # turn into network, the only useful form of this thing
    $nestingNet = _getNestingNet(\%acc2nesting, \%acc2clan); # this will overwrite the global, but not returned references we may have loaded earlier
    # copy new data to global references
    ($acc2name, $acc2de, $acc2ds2ga, $acc2length, $acc2clan, $acc2type) = (\%acc2name, \%acc2de, \%acc2ds2ga, \%acc2length, \%acc2clan, \%acc2type);
    
    # done, return structures!
    return (\%acc2name, \%acc2de, \%acc2ds2ga, \%acc2length, \%acc2clan, \%acc2type, $nestingNet); # , \%acc2nesting
}

sub clans {
    # a parser of the useful data contained in a small Pfam-C file.  It should be about as fast to load the data directly from here as it would be to parse it into hashes and load those instead as needed, but I think this saves me scripting time by not having to define so many random files.  Of course, this doesn't work for Pfam 23 and older (the file format changed).
    # usage
    # my ($clan2name, $clan2de, $clan2cc, $clan2scop, $clan2cath, $clan2merops) = ParsePfam::clans($fiCla);
    my ($fi) = @_; # input file
    
    # define hash structures
    my %clan2name;
    my %clan2de; # short description
    my %clan2cc; # long comments
    # DR can point to several databases, including (examples):
    my %clan2scop;   # DR   SCOP; 49944;       # 263 entries
    my %clan2cath;   # DR   CATH; 3.30.60.30;  # 139 entries
    my %clan2merops; # DR   MEROPS; I1;        # 6 entries, a peptidase database
    
    # a temporary variable that holds the current accession, so everything else can be mapped
    my $acc = ''; # current clan accession
    my $name = ''; # current name, comes before accession sadly
    
    # open file, start reading
    my $fhi = FileGz::getInFh($fi);
    while (<$fhi>) {
	if ( my ($type, $s) = /^\#=GF (\w+)\s+(.+)$/ ) { # things we're parsing for... many lines are not interesting
	    if ($type eq 'ID') { $name = $s; }
	    elsif ($type eq 'AC') {
		$acc = $s; # read new accession
		$acc =~ s/\.\d+$//; # remove version digits
		$clan2name{$acc} = $name; # name gets mapped here
	    }
	    elsif ($type eq 'DE') { $clan2de{$acc} = $s; }
	    elsif ($type eq 'CC') { $clan2cc{$acc} .= $s.' '; } # comments may occupy multiple lines, so concatenate (do we need to add spaces?)
	    # MB gets ignored, redundat with PfamADat data
	    elsif ($type eq 'DR') {
		my ($db, $dbAcc) = split /;\s*/, $s; # parse the two entries (database name and accession).  Removes all semicolons
		if ($db eq 'SCOP') { push @{$clan2scop{$acc}}, $dbAcc; }
		if ($db eq 'CATH') { push @{$clan2cath{$acc}}, $dbAcc; }
		if ($db eq 'MEROPS') { push @{$clan2merops{$acc}}, $dbAcc; }
	    }
	}
    }
    close $fhi;
    
    # copy new data to global references
    ($clan2name, $clan2de, $clan2cc, $clan2scop, $clan2cath, $clan2merops) = (\%clan2name, \%clan2de, \%clan2cc, \%clan2scop, \%clan2cath, \%clan2merops);
    
    # done, return structures!
    return (\%clan2name, \%clan2de, \%clan2cc, \%clan2scop, \%clan2cath, \%clan2merops);
}


### experimental function (for some benchmarks in which newer clans seemed to make an improvement, but not for general use)

sub datExt {
    # loads pfamDat for old data too (last one loaded, so it sticks around)!  Extends old data in place! (so rest of code doesn't have to change!)
    # an analysis suggests keeping old clans, except gains (old undef, new def) and differences (both def) will improve old clans.  Ignore losses and dead families (keep them the same in old).
    my ($fi1, $fi2) = @_;
    # get old and new data
    # get new first because we want to use the old data exclusively in the end, let that one hang around!
    dat($fi2);
    my $acc2clan2 = $acc2clan; # copy reference, will be overwritten next
    dat($fi1); # these references stick around
    
    # navigate only old accs (the ones we have data for), but note it's full list, the didn't have to have had clans before!
    foreach my $acc (keys %$acc2name) {
	# get newest clan annotations, could be undef
	my $clan2 = $acc2clan2->{$acc};
	if ($clan2) { # this means we ignore losses as well as things that didn't have clans in both versions
	    $acc2clan->{$acc} = $clan2; # overwrite old clan with new... this means we either keep things the same, or we made additions, or we made a change, all cases are acceptably dealt this way.  Let's not waste time comparing values then!
	}
    }
    # done, new extended clans are ready to be used out of the box!
}


### internal function

sub _getNestingNet {
    # simply turns a map between nesting families into a bidirectional network.  Excludes nestings between families of the same clan (read below).
    # only used internally by dat()
    # 2012-09-13 14:16:38 EDT
    # update: discovered some nesting pairs are of the same clan.  While in a broader setting this might be acceptable, in out setup (which only uses ranges and not full HMM alignments) a true nesting is indistinguishable from a bad/redundant overlap between members of the same clan.  To be consistent with other code, to keep reasonable assumptions, let's remove these same-clan pairs from here.
    my ($acc2nesting, $acc2clan) = @_;
    my %net;
    # navigate keys
    while (my ($acc1, $partners) = each %$acc2nesting) {
	my $clan1 = $acc2clan->{$acc1};
	foreach my $acc2 (@$partners) { # @acc2
	    my $clan2 = $acc2clan->{$acc2};
	    next if $clan1 && $clan2 && $clan1 eq $clan2; # ignore nesting between two accs of the same clan!
	    # add to net both ways
	    $net{$acc1}{$acc2} = 1;
	    $net{$acc2}{$acc1} = 1;
	}
    }
    return \%net;
}


1;

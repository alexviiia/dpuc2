# Copyright (C) 2014 Alejandro Ochoa, Singh Research Group, Princeton University
# This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.

package FileGz;
our $VERSION = 1.01;
#use PerlIO::gzip;
use File::Copy (); # core module, only used for File::Copy::move
use strict;

# 1.01, unpublished so far
# 2014-11-24
# - myRename added option $ifExists, so it returns quietly instead of dying if file to rename does not exist (mostly used in my personal code only)
# - mvBak added to move fi to fi~ handling compressed extension adequately (so fi.gz becomes fi.gz~ instead of fi~.gz).  
# - both myRename and mvBak return destination file (real one) if there was a move.
# - both myRename and mvBak die if base rename() failed for some other reason than the ones we had already covered.
# - for maximum functionality, both now use File::Copy::move instead of rename, the latter of which doesn't work across system boundaries.
# - myRenameReal is shared part used by myRename and mvBak, it complains but doesn't perform checks directly.  Unlike the other two, this one also works on directories!
# - myUnlink doesn't complain anymore about removing files that evaluate to false ('', 0, or undef)

# Main goal: handle files whether they are compressed or not, there's no need to specify ahead of time if some assumptions hold.  Mostly built for gzip, but in some cases, other compression formats (bzip2) and (uncompressed) stdin/stdout are also handled!  Lastly, catching and dying due to errors is important to me.
# Main assumption: only one version of a file (plain or compressed) exists, so it should be obvious which one to open.  If this isn't so, then your code might misbehave, particularly if the versions differ.  This assumption is always checked in the background for you, and code will warn if it is not met (but it will always be non-fatal).
# I included other related file-open operations even when they do not or cannot use compression, just for consistency
# History: based/copied on some code from the Keating Lab at MIT, circa 2006.  I think it was Gevorg Grigoryan's code.

### TO DO
# flag to return undef instead of dying? (for sitivs)

sub getInFh {
    my ($fi, $layer) = @_;
    $layer = '' unless defined $layer; # default that won't make things die
    ($fi, my $format) = solveInput($fi); # check that there's only one version (fatal otherwise), also preserve format
    my $fhi;
    if ($format eq '-') {
	$fhi = \*STDIN; # useful hack
    }
    elsif ($format eq '') {
	open($fhi, '<'.$layer, $fi) || die "Could not open for reading $fi: $!";
    }
    elsif ($format eq 'gzip') {
	open($fhi, '-|'.$layer, "gunzip -c $fi") || die "Could not open for reading $fi: $!";
#	open($fhi, '<:gzip', $fi) || die "Could not open for reading $fi: $!";
    }
    elsif ($format eq 'bzip2') {
	open($fhi, '-|'.$layer, "bzip2 -dc $fi") || die "Could not open for reading $fi: $!";
    }
    else { die "Compression format not recognized: $format\n"; }
    return $fhi;
}

sub getOutFh {
    my ($fo, $format, $layer) = @_;
    $layer = '' unless defined $layer; # default that won't make things die
    ($fo, $format) = solveOutput($fo, $format); # do some checks on this input, clean it up too
    return $fo if $format eq 'GLOB'; # if input looks like a filehandle, return it without doing anything else (this is very convenient!)
    my $fho;
    if ($format eq '-') {
	$fho = \*STDOUT; # useful hack.  Note STDOUT is never gzipped, even if that option is passed!
    }
    elsif ($format eq 'gzip') {
	open($fho, '|-'.$layer, "gzip -c > $fo") || die "Could not open for writing $fo: $!";
#	open($fho, '>:gzip', $fo) || die "Could not open for writing $fo: $!";
    }
    elsif ($format eq 'bzip2') {
	open($fho, '|-'.$layer, "bzip2 -zcf > $fo") || die "Could not open for writing $fo: $!";
    }
    elsif ($format eq '') {
	open($fho, '>'.$layer, $fo) || die "Could not open for writing $fo: $!";
    }
    else { die "Compression format not recognized: $format\n"; }
    return $fho;
}

sub getAppFh {
    # this only works with uncompressed files...
    my ($fo, $layer) = @_;
    $layer = '' unless defined $layer; # default that won't make things die
    die "Undefined or empty file name passed!" if !defined $fo || $fo eq '';
    die "Could not open $fo: append only works on plain files" if $fo =~ /\.gz$/i;
    open(my $fho, '>>'.$layer, $fo) || die "Could not open $fo: $!";
    return $fho;
}

sub getPipeFh {
    my ($fo, $layer) = @_;
    $layer = '' unless defined $layer; # default that won't make things die
    die "Undefined or empty program name passed!" if !defined $fo || $fo eq '';
    open(my $fho, '|-'.$layer, $fo) || die "Could not open $fo: $!";
    return $fho;
}

sub myUnlink {
    # handles compression extensions transparently (e.g. will delete readme.txt.gz if readme.txt was an input and readme.txt does not exist, but readme.txt.gz does exist).  Non fatal.
    # also provides informative warning messages when files are not found or could not be deleted
    # @_ is the list of files to delete
    foreach my $fi (@_) {
	if ( $fi = realFile($fi) ) { # if defined, then it was found
	    unlink($fi) || warn "Note: file $fi was found, but could not be deleted: $!";
	}
	else { warn "Note: file $fi was not found, although it was scheduled for deletion!" if $fi; } # if it was defined, sometimes it isn't!
    }
}

sub myRename {
    # move fi to fo, or their compressed versions
    my ($fi, $fo, $ifExists) = @_;
    # first let's figure out what the real input file is (its real path) and its format
    ($fi, my $format) = solveInput($fi, 1); # non-fatal version, keep format!
    unless ($fi) {
	# if file to move does not exists, there are two possible outcomes
	if ($ifExists) { return; } # if this option is set, return quietly without doing anything else
	else {
	    die "Fatal: file to move $fi does not exist!\n"; # should stop if this is missing!  Seems like a big deal to me!
	}
    }
    # now get output in the same format, but depending on what $fo was, there could be funny behaviors (because $fo's extension takes precedence over $format).
    ($fo, my $format2) = solveOutput($fo, $format);
    die "Fatal: input and output files $fi, $fo, in myRename are not of the same format!\n" unless $format eq $format2;
    # if everything was good, let's just rename.  Errors might happen here too, catch them!
    myRenameReal($fi, $fo);
    return $fo; # return this file, it may be of interest!
}

sub mvBak {
    # move fi to a backup file fi~, handling the compressed form adequately!
    # (so fi.gz becomes fi.gz~ instead of fi~.gz)
    # if output exists, it always gets overwritten (this is expected for backup files), so we don't even check
    my ($fi, $ifExists) = @_;
    # get real file, but it may not exist
    my $fiReal = realFile($fi);
    unless ($fiReal) {
	# if file to move does not exists, there are two possible outcomes
	if ($ifExists) { return; } # if this option is set, return quietly without doing anything else
	else {
	    die "Fatal: file to move $fi does not exist!\n"; # should stop if this is missing!  Seems like a big deal to me!
	}
    }
    # simple move here... real file gets tilde added
    my $fo = $fiReal.'~';
    # if everything was good, let's just rename.  Errors might happen here too, catch them!
    myRenameReal($fiReal, $fo);
    return $fo; # return this file, it may be of interest!
}

sub myRenameReal {
    # this version moves only real files (so if they're compressed, they must include the right extensions).  It complains if things don't work.  It also works on directories!
    # this is separated from myRename so it can be shared with mvBak, and also be used on its own on directories
    # note here files or dirs aren't directly tested for existence, although File::Copy will probably complain if there were problems.
    my ($fi, $fo) = @_;
    File::Copy::move($fi, $fo) || die "Fatal: could not rename file $fi to $fo: $!";
}

sub realFile {
    # for a given file, returns itself if it exists, or the compressed version if that exists, or undefined if neither exist
    # never fatal!
    my ($fi) = @_;
    ($fi) = solveInput($fi, 1); # now a wrapper for solveInput, but non-fatal.  Throw away format.
    return $fi;
}

### INTERNAL FUNCTIONS, really the workhorses behind this code

sub solveInput {
    # looks for all versions of a given file, supposedly with the intent of reading it
    # warns when we find multiple (possibly conflicting) versions
    # dies, unless asked not to, if inputs are missing or no versions to read were found
    # NOTE: only returns the one pair of file-format, or the first one if multiple versions were found
    my ($fi, $nonFatal) = @_;
    # First special case, if there wasn't a valid input...
    if (!defined $fi || $fi eq '') {
	if ($nonFatal) { return (undef, undef); }
	else { die "Fatal: undefined or empty file name passed!"; }
    }
    # Second special case, stdin gets preserved (so this works with getInFh())
    if ($fi eq 'STDIN' || $fi eq '-') {
	return ('-','-'); # in this case there's only this one "file", there are no alternative version to look for, and note format is also '-'
    }
    # now for regular files...
    # first try to remove compressed extensions, if they're present.  It'll be easier to look for all versions later
    $fi =~ s/\.gz$//i;
    $fi =~ s/\.bz2$//i;
    # start scanning for some possibilities
    my @files; # array of arrays of file-format pairs
    if (-f $fi) { push @files, [$fi, '']; }
    if (-f $fi.'.gz') { push @files, [$fi.'.gz', 'gzip']; }
    if (-f $fi.'.GZ') { push @files, [$fi.'.GZ', 'gzip']; }
    if (-f $fi.'.bz2') { push @files, [$fi.'.bz2', 'bzip2']; }
    if (-f $fi.'.BZ2') { push @files, [$fi.'.BZ2', 'bzip2']; }
    # time to warn or die if needed
    if (@files>1) {
	# for multiple files (a serious conflict), let's warn in non-fatal case
	my @fis = map { $_->[0] } @files; # flatten double array to single array of files only
	warn "Note: multiple versions of file $fi were found: @fis\n";
    }
    if (@files == 0) {
	if ($nonFatal) { return (undef, undef); } # returned undefined quietly
	else { die "Fatal: no versions of file $fi were found!\n"; }
    }
    else {
	return @{$files[0]}; # return "the one" file-format pair.
    }
}

sub solveOutput {
    # from this given command, figure out what the actual output will be, and the final format (which might have been implicit in input's extension)
    # also, make sure output doesn't have multiple existing versions and that, if we're overwritting, that we're overwritting the only copy
    # as solveInput for the input version, this one will do plenty of complaining if needed
    my ($fo, $format) = @_;
    
    # these special cases need no further processing
    return ($fo, 'GLOB') if ref($fo) eq 'GLOB'; # if input looks like a filehandle, return it without doing anything else (this is very convenient!)
    die "Fatal: Passed a non-GLOB reference instead of a scalar or GLOB!" if ref $fo; # complain if it's any other kind of reference
    die "Fatal: Undefined or empty file name passed!" if !defined $fo || $fo eq '';
    if ($fo eq 'STDOUT' || $fo eq '-') {
	return ('-','-');
    }
    
    # now for regular files, which need some checking
    # infer from extension, ignore and correct? input $format value.  In these cases $fo gets preserved
    if ($fo =~ /\.gz$/i) { $format = 'gzip'; }
    elsif ($fo =~ /\.bz2$/i) { $format = 'bzip2'; }
    # use $format if it appeared to be a regular file input.  In these cases $format gets preserved.
    elsif (!defined $format) { $format = ''; } # will have a plain file
    elsif ($format eq 'gzip') { $fo .= '.gz'; }
    elsif ($format eq 'bzip2') { $fo .= '.bz2'; }
    
    # now that we've figured out what the real output file will be, let's make sure there aren't multiple versions lying around
    if (my $foReal = realFile($fo)) { # non fatal, because it's ok if nothing exists.  But if multiple versions already exist, they'll only appear as a warning.
	warn "Note: File to write to ($fo) will not overwrite existing version ($foReal)!\n" unless $fo eq $foReal;
    }
    return ($fo, $format); # done, return this for actually opening file handles, or something
}

1;

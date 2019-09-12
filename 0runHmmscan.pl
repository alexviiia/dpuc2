# Copyright (C) 2014-2019 Alejandro Ochoa, Singh Research Group, Princeton University
# This file is part of DomStratStats and dPUC.
# DomStratStats and dPUC are free software: you can redistribute them and/or modify them under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
# DomStratStats and dPUC are distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with DomStratStats and dPUC.  If not, see <http://www.gnu.org/licenses/>.

# core Perl modules
use FindBin ();
# local modules
use lib '.';
use Hmmer3ScanTab;
use strict;

# this hack is to announce the right package, without actually having to change the code...
# if everything is missing, nothing is printed
# if both are present, two lines are printed
# so this handles all cases
my $version_strings = '';
# in each case add initial comment and final newline
if (eval "use DomStratStats; 1") {
    $version_strings .= '# ' . DomStratStats::version_string() . "\n";
}
if (eval "use Dpuc; 1") {
    $version_strings .= '# ' . Dpuc::version_string() . "\n";
    # this gets rid of an annoying warning that only happens in this setting
    Inline->init();
}

# get input, ask politely for it otherwise
my ($hmmscan, $pfamA, $fiSeq, $fo, $file_stdout) = @ARGV;

unless ($fo) {
    print "# $FindBin::Script: Get domain predictions from your protein sequences
$version_strings# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

Usage: perl -w $FindBin::Script <hmmscan> <Pfam-A.hmm> <FASTA input> <output table> \\
         [<file stdout>]

The required inputs are
    <hmmscan>      the path to the HMMER3 hmmscan executable.
    <Pfam-A.hmm>   the path to your HMM library of choice (in HMMER3 format).
    <FASTA input>  the FASTA sequence file, may be compressed with gzip.
    <output table> the hmmscan output is plain text table delimited by whitespace (always
                   uncompressed).

The optional input is
    <file stdout>  the file to which hmmscan's standard output goes, including alignments
                   (default /dev/null)

This script changes hmmscan parameters that are important for dPUC and DomStratStats to
work.  In particular, it forces outputs to report p-values in place of E-values, and it 
relaxes the heuristic p-value filters to allow more predictions through.  This script sets 
the most stringent p-value threshold at 1e-4.
";
    exit 0;
}

# these stay undefined (sets default behaviors)
my $pCut;
my $singleThread;

# this sub does all the magic
Hmmer3ScanTab::runHmmscan($hmmscan, $pfamA, $fiSeq, $fo, $pCut, $singleThread, $file_stdout);

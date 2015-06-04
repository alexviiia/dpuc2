# Copyright (C) 2014 Alejandro Ochoa, Singh Research Group, Princeton University
# This file is part of DomStratStats and dPUC.
# DomStratStats and dPUC are free software: you can redistribute them and/or modify them under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
# DomStratStats and dPUC are distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with DomStratStats and dPUC.  If not, see <http://www.gnu.org/licenses/>.

my $VERSION = '1.01';
use lib '.';
use Hmmer3ScanTab;
use strict;

# get input, ask politely for it otherwise
my ($hmmscan, $pfamA, $fiSeq, $fo) = @ARGV;

die "# $0 $VERSION - Get domain predictions from your protein sequences
# DomStratStats 1.xx, viiia.org/domStratStats
# dPUC 2.xx, viiia.org/dpuc2
# Alejandro Ochoa, John Storey, Manuel Llinás, and Mona Singh.
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

Usage: perl -w $0 <hmmscan> <Pfam-A.hmm> <FASTA input> <output table>

The required inputs are
    <hmmscan>      the path to the HMMER3 hmmscan executable.
    <Pfam-A.hmm>   the path to your HMM library of choice (in HMMER3 format).
    <FASTA input>  the FASTA sequence file, may be compressed with gzip.
    <output table> the hmmscan output is plain text table delimited by whitespace (always
                   uncompressed).

See the online manual for more info.
" unless $fo;

# this sub does all the magic
Hmmer3ScanTab::runHmmscan($hmmscan, $pfamA, $fiSeq, $fo);

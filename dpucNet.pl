# Copyright (C) 2014-2019 Alejandro Ochoa, Singh Research Group, Princeton University
# This file is part of dPUC.
# dPUC is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
# dPUC is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with dPUC.  If not, see <http://www.gnu.org/licenses/>.

# core Perl modules
use FindBin ();
# local modules
use lib '.';
use Dpuc; # just for version_string :(
use DpucNet;
use strict;

# constants
my $comp = 'gzip'; # compress output
my $verbose = 1;

my ($fi, $fo) = @ARGV;

unless ($fo) {
    print "# $FindBin::Script: Extract the context count network from Pfam predictions
# " . Dpuc::version_string() . "
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

Usage: perl -w $FindBin::Script <Pfam-A.full.uniprot> <dpucNet output>

The required inputs are
    <Pfam-A.full.uniprot>  Input \"full\" domain prediction file from Pfam (29-latest).
                           Use Pfam-A.full for Pfam 23-28.
    <dpucNet output>       Directed context network of domain family pair counts.

Input file may be compressed with gzip, and may be specified with or without the .gz 
extension.  Output file will be automatically compressed with gzip.

Run once per Pfam release.  This Pfam-specific script counts every ordered domain family 
pair observed in the full predictions of standard Pfam releases.  The code is optimized to
parse such large files quickly and store domain predictions with minimal memory usage.
However, for the latest Pfam releases it remains a considerable amount of time and memory.
Most users should skip this step, downloading instead the desired precomputed dPUC context 
networks file from:

    https://github.com/alexviiia/dpuc2-data

This script is provided for completeness and transparency.
";
    exit 0;
}

# this does all the magic!
DpucNet::makeNet($fi, $fo, $comp, $verbose);

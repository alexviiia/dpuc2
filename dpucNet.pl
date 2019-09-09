# Copyright (C) 2014-2019 Alejandro Ochoa, Singh Research Group, Princeton University
# This file is part of dPUC.
# dPUC is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
# dPUC is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with dPUC.  If not, see <http://www.gnu.org/licenses/>.

my $VERSION = '1.03';
use lib '.';
use DpucNet;
use strict;

# constants
my $comp = 'gzip'; # compress output
my $verbose = 1;

my ($fi, $fo) = @ARGV;

unless ($fo) {
    print "# $0 $VERSION - Extract the context count network from Pfam predictions
# dPUC       ".(sprintf '%0.2f', $DpucNet::VERSION)." - https://github.com/alexviiia/dpuc2
# Alejandro Ochoa, John Storey, Manuel Llinás, and Mona Singh.
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

Usage: perl -w $0 <Pfam-A.full.uniprot> <dpucNet output>

The required inputs are
    <Pfam-A.full.uniprot>  Input \"full\" domain prediction file from Pfam.
    <dpucNet output>       Directed context network of domain family pair counts.

Input file may be compressed with gzip, and may be specified with or without the .gz 
extension.  Output file will be automatically compressed with gzip.

See the online manual for more info.
";
    exit 0;
}

# this does all the magic!
DpucNet::makeNet($fi, $fo, $comp, $verbose);

# Copyright (C) 2014 Alejandro Ochoa, Singh Research Group, Princeton University
# This file is part of dPUC.
# dPUC is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
# dPUC is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with dPUC.  If not, see <http://www.gnu.org/licenses/>.

my $VERSION = '1.00';
use lib '.';
use Hmmer3ScanTab;
use Dpuc;
use ParsePfam;
use strict;

# this script accepts as input a hmmer3 hmmscan output file, and outputs a table file with the predictions that dpuc2 keeps with these parameters
my ($fiPfamADat, $fiCounts, $fi, $fo) = @ARGV;

die "# $0 $VERSION - Produce dPUC domain predictions from raw hmmscan data
# dPUC ".(sprintf '%0.2f', $Dpuc::VERSION).", viiia.org/dpuc2
# Alejandro Ochoa, John Storey, Manuel Llinás, and Mona Singh.
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

Usage: perl -w $0 <Pfam-A.hmm.dat> <dpucNet> <input table> <output table>

The required inputs are
    <Pfam-A.hmm.dat>  Pfam annotation file, contains GA thresholds and nesting network.
    <dpucNet>         Directed context network of domain family pair counts.
    <input table>     The output from hmmscan (previous script).
    <output table>    Input with most domains removed by dPUC. Format is identical to input.

It is your responsability to ensure that the same Pfam version is used for Pfam-A.hmm.dat, 
the dPUC context network, and the input table.

Input file may be compressed with gzip, and may be specified with or without the .gz 
extension.  Output file will be automatically compressed with gzip.

See the online manual for more info.
" unless $fo;

# dpuc2 parameter defaults (all hardcoded for now, for simplicity)
my $cut = '1e-4';
my $panning = 'inf'; # inf=fully directional scoring, did best in my dissertation's benchmarks
my $alphaSelfToTrans = 0;
my $scaleExp = 23;
my $timeout = 1; # default lp_solve timeout in seconds
my $scaleContext; # default undefined (translates to no scaling, or 1)
my $shiftContext; # default undefined (translates to no shifting, or 0)
my ($fCut, $lCut) = qw(.5 40); # new permissive overlap params
my $cCutNet = 2; # minimum count to filter in input dpucNet
my $removeOvsPRank; # boolean, default false
my $removeOvsCodd; # boolean, default false

# constants
my $comp = 'gzip';

# get Pfam GA thresholds and nesting net
ParsePfam::dat($fiPfamADat); # populates $ParsePfam::acc2ds2ga and $ParsePfam::nestingNet
# this wrapper choses the right file, given the parameters, and does all necessary processing
my $net = Dpuc::loadNet($fiCounts, $panning, $alphaSelfToTrans, $scaleExp, $scaleContext, $shiftContext, $cCutNet);
# get necessary sequence thresholds
my $acc2ts = Dpuc::getPfamThresholdsSeq($ParsePfam::acc2ds2ga);

# now read the full files entry-by-entry
my ($fhi, $header, $protNext, $hitNext) = Hmmer3ScanTab::parseInit($fi);
# prepare output
my $fho = Hmmer3ScanTab::printInit($fo, $header, $comp);
while ($protNext) {
    # read next protein
    (my $prot, my $hits, $protNext, $hitNext) = Hmmer3ScanTab::parseProt($fhi, $protNext, $hitNext);
    # apply Dpuc, overwriting hits
    ($hits) = Dpuc::dpuc($hits, $net, $ParsePfam::nestingNet, $ParsePfam::acc2ds2ga, $acc2ts, $timeout, $cut, undef, $fCut, $lCut, $removeOvsPRank, $removeOvsCodd);
    # print to output
    Hmmer3ScanTab::printProt($fho, $hits);
}

# done with Hits input and output files
close $fhi;
close $fho;

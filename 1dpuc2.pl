# Copyright (C) 2014 Alejandro Ochoa, Singh Research Group, Princeton University
# This file is part of dPUC.
# dPUC is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
# dPUC is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with dPUC.  If not, see <http://www.gnu.org/licenses/>.

my $VERSION = '1.01';
use lib '.';
use Hmmer3ScanTab;
use Dpuc;
use ParsePfam;
use strict;

# this script accepts as input a hmmer3 hmmscan output file, and outputs a table file with the predictions that dpuc2 keeps with these parameters
my ($fiPfamADat, $fiCounts, $fi, $fo0, @cuts) = @ARGV;

# 2015-07-27 14:31:03 EDT - v1.01
# - Now multiple p-value thresholds for candidate domains can be set as an option (on the command line).  The script used to have a single threshold hardcoded, and did not handle outputs for multiple thresholds.
# - Removed comment about Pfam versions having to match (since now they are explicitly checked, user isn't trusted!)

die "# $0 $VERSION - Produce dPUC domain predictions from raw hmmscan data
# dPUC ".(sprintf '%0.2f', $Dpuc::VERSION).", viiia.org/dpuc2
# Alejandro Ochoa, John Storey, Manuel Llinás, and Mona Singh.
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

Usage: perl -w $0 <Pfam-A.hmm.dat> <dpucNet> <input table> <output table> [<p-value cuts...>]

The required inputs are
    <Pfam-A.hmm.dat>   Pfam annotation file, contains GA thresholds and nesting network.
    <dpucNet>          Directed context network of domain family pair counts.
    <input table>      The output from hmmscan (previous script).
    <output table>     Input with most domains removed by dPUC. Format is identical to input.
                       File path will be used as the base to generate multiple files, one for 
                       each p-value thresholds, if more than one threshold is requested.

The following parameters are optional.
    <p-value cuts...>  The p-value thresholds for candidate domains [default 1e-4]. If it is 
                       a single threshold, the output table is output to the path provided.
                       If multiple thresholds are provided, the output table has these 
                       thresholds inserted as text just before the file extension, creating 
                       separate outputs for each threshold.

If more than one p-value threshold is provided on candidate domains:
Note that dPUC predictions are not necessarily nested as candidate domain thresholds are varied, 
since dPUC solves a complicated pairwise optimization problem, so these files are not redundant 
and are not obtained by filtering the largest file. Each output is as if dPUC had been run anew
on the given p-value threshold, although the code does this efficiently by factoring out shared 
steps between thresholds.

Input file may be compressed with gzip, and may be specified with or without the .gz 
extension.  Output file will be automatically compressed with gzip.

See the online manual for more info.
" unless $fo0;

# set default threshold if none were provided
@cuts = ('1e-4') unless @cuts;

# dpuc2 parameter defaults (all hardcoded for now, for simplicity)
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

# pre-processing: make sure p-values are listed descending!!!
@cuts = sort { $b <=> $a } @cuts;

# pre-processing: create output file paths for the case of multiple thresholds
my @fos = ($fo0); # singleton case
if (@cuts > 1) {
    if ($fo0 =~ /^(.*)(\.\w+(?:\.gz)?)$/) {
	# if here, we sucessfully extracted base and extension (with compression added further at the end)
	@fos = map { $1.'.p'.$_.$2 } @cuts;
    } else {
	# if there's no extension or something else weird is going on, just append cuts to the end
	@fos = map { $fo0.'.p'.$_ } @cuts;
    }
#    print "Outputs are: \n".join("\n", @fos)."\n";
}

# get Pfam GA thresholds and nesting net
ParsePfam::dat($fiPfamADat); # populates $ParsePfam::acc2ds2ga and $ParsePfam::nestingNet
my @accs = keys %$ParsePfam::acc2name; # this helps check that dpucNet and pfamDat agree (they won't if versions are different)
# this wrapper choses the right file, given the parameters, and does all necessary processing
my $net = Dpuc::loadNet(\@accs, $fiCounts, $panning, $alphaSelfToTrans, $scaleExp, $scaleContext, $shiftContext, $cCutNet);
# get necessary sequence thresholds
my $acc2ts = Dpuc::getPfamThresholdsSeq($ParsePfam::acc2ds2ga);

# now read the full files entry-by-entry
my ($fhi, $header, $protNext, $hitNext) = Hmmer3ScanTab::parseInit($fi);

# prepare outputs
my @fhos;
foreach my $fo (@fos) {
    my $fho = Hmmer3ScanTab::printInit($fo, $header, $comp);
    push @fhos, $fho;
}

# main dPUC processing, efficiently creating curves if necessary (for more than one p-value threshold)
while ($protNext) {
    # read next protein
    (my $prot, my $hits, $protNext, $hitNext) = Hmmer3ScanTab::parseProt($fhi, $protNext, $hitNext);
    # apply Dpuc, obtaining an array of predictions corresponding to each p-value threshold used
    my ($cut2hits, undef, undef, $errStr) = Dpuc::dpuc($hits, $net, $ParsePfam::nestingNet, $ParsePfam::acc2ds2ga, $acc2ts, $timeout, \@cuts, undef, $fCut, $lCut, $removeOvsPRank, $removeOvsCodd);
    # errors need to be included in every output!
    if ($errStr) {
	foreach my $fho (@fhos) {
	    print $fho $errStr; # include error message in output!
	    close $fho; # make sure file is closed correctly
	}
	die $errStr; # die with the same message for STDERR
    }
    # print predictions to outputs
    for (my $i = 0; $i < @cuts; $i++) {
	Hmmer3ScanTab::printProt($fhos[$i], $cut2hits->[$i]); # this sends the correct prediction to the correct file
    }
}

# done with Hits input and output files
close $fhi;
foreach my $fho (@fhos) {
    close $fho;
}

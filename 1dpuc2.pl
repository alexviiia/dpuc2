# Copyright (C) 2014 Alejandro Ochoa, Singh Research Group, Princeton University
# This file is part of dPUC.
# dPUC is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
# dPUC is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with dPUC.  If not, see <http://www.gnu.org/licenses/>.

my $VERSION = '1.04';
use Getopt::Long (); # Core perl package, should always be available!
use lib '.';
use Hmmer3ScanTab;
use Dpuc;
use ParsePfam;
use strict;

# 2015-07-27 14:31:03 EDT - v1.01
# - Now multiple p-value thresholds for candidate domains can be set as an option (on the command line).  The script used to have a single threshold hardcoded, and did not handle outputs for multiple thresholds.
# - Removed comment about Pfam versions having to match (since now they are explicitly checked, user isn't trusted!)

# 2015-08-12 21:18:19 EDT - v1.02
# - removed removeOvsCodd and removeOvsPRank option
# 2015-09-04 14:20:56 EDT - v1.02 still (not published yet)
# - changed default scale parameter from 23 to 3, which is more sensible.
# - added Getopt::Long, to set a whole bunch of params on command line

# 2015-10-30 22:31:01 EDT - v1.03
# - internally, alpha (pseudocount) is set directly rather than through scale or as logarithm. Here we change the parameter passing to abide to this change.

# 2016-06-28 19:28:39 EDT - v1.04
# - fixed a "--pvalues" bug, where if we added "--pvalues 1e-2", the output was supposed to be (and now is) only for that threshold, but the bug added it to the default of 1e-4, resulting in two outputs (for 1e-2 and 1e-4).


# clean script name
my ($scriptName) = $0 =~ /(\w+\.pl)$/;

# parameters accessible as options
my $cutDef = '1e-4'; # default p-value threshold for candidate domains
my @cuts = (); # list of p-value thresholds for candidate domains
my $alpha = 1e-2; # new scale default, changed for version 2.04!
my ($fCut, $lCut) = qw(.50 40); # new permissive overlap params
my $cCutNet = 2; # minimum count to filter in input dpucNet
my $timeout = 1; # default lp_solve timeout in seconds
my $comp = 'gzip';

my $usage = "# $scriptName $VERSION - Produce dPUC domain predictions from raw hmmscan data
# dPUC ".(sprintf '%0.2f', $Dpuc::VERSION).", viiia.org/dpuc2
# Alejandro Ochoa, John Storey, Manuel Llinás, and Mona Singh.
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

Usage: perl -w $scriptName [-options] <Pfam-A.hmm.dat> <dpucNet> <input table> <output table>

The required inputs are
    <Pfam-A.hmm.dat>   Pfam annotation file, contains GA thresholds and nesting network.
    <dpucNet>          Directed context network of domain family pair counts.
    <input table>      The output from hmmscan (previous script).
    <output table>     Input with most domains removed by dPUC. Format is identical to input.
                       File path will be used as the base to generate multiple files, one for 
                       each p-value thresholds, if more than one threshold is requested.

Options:
    --pvalues <x...>   The p-value thresholds for candidate domains  [$cutDef]
                       If multiple thresholds are provided, the output table has these 
                       thresholds inserted as text just before the file extension, creating 
                       separate outputs for each threshold.
    --alpha <x>        Pseudocount (sym-Dirichlet param) for context scores  [$alpha]
    --cutNet <x>       Keep context network counts >= cutNet  [$cCutNet]
    --fCut <x>         Permissive overlap max proportion  [$fCut]
    --lCut <x>         Permissive overlap max length in amino acids  [$lCut]
    --timeout <x>      lp_solve timeout in seconds  [$timeout]
    --noGzip           Outputs will not be compressed (default is gzip compression).

If more than one p-value threshold is provided on candidate domains:
Note that dPUC predictions are not necessarily nested as candidate domain thresholds are varied, 
since dPUC solves a complicated pairwise optimization problem, so these files are not redundant 
and are not obtained by filtering the largest file. Each output is as if dPUC had been run anew
on the given p-value threshold, although the code does this efficiently by factoring out shared 
steps between thresholds.

Input file may be compressed with gzip, and may be specified with or without the .gz 
extension.  Output file will be automatically compressed with gzip unless --noGzip is used.

See the online manual for more info.
";

# try to get parameters from command line
Getopt::Long::GetOptions(
    'pvalues=f{1,}' => \@cuts,
    'alpha=f' => \$alpha,
    'cutNet=f' => \$cCutNet,
    'fCut=f' => \$fCut,
    'lCut=f' => \$lCut,
    'timeout=f' => \$timeout,
    'noGzip' => sub { $comp = '' },
    ) or die $usage;

# this script accepts as input a hmmer3 hmmscan output file, and outputs a table file with the predictions that dpuc2 keeps with these parameters
my ($fiPfamADat, $fiCounts, $fi, $fo0) = @ARGV;
die $usage unless $fo0;

# dpuc2 hardcoded parameters (untested, users shouldn't use them!)
my $scaleContext; # default undefined (translates to no scaling, or 1)
my $shiftContext; # default undefined (translates to no shifting, or 0)

# pre-processing: make sure p-values are listed descending!!!
if (@cuts) {
    @cuts = sort { $b <=> $a } @cuts;
} else {
    @cuts = ($cutDef); # use default threshold as sole threshold if none were specified!
}

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
my $net = Dpuc::loadNet(\@accs, $fiCounts, $alpha, $scaleContext, $shiftContext, $cCutNet);
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
    my ($cut2hits, undef, undef, $errStr) = Dpuc::dpuc($hits, $net, $ParsePfam::nestingNet, $ParsePfam::acc2ds2ga, $acc2ts, $timeout, \@cuts, undef, $fCut, $lCut);
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

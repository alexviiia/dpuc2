# Copyright (C) 2014-2019 Alejandro Ochoa, Singh Research Group, Princeton University
# This file is part of dPUC.
# dPUC is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
# dPUC is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with dPUC.  If not, see <http://www.gnu.org/licenses/>.

package EncodeIntPair;
our $VERSION = 1.01;
use strict;

### map symmetric pair w/ diagonal

sub t {
    # maps $i, $j into $t in a triangle matrix of any size!
    # maps contiguously if we do: i in 0..total, j in 0..i (allowing diagonal), good for loops so we don't have to calculate this, but we do need to calculate it to map back
    # in current form, assumes i>=j
    my ($i, $j) = @_;
    return $j + $i*($i+1)/2;
}

# for completeness, should add reverse map of t()
# don't have it in perl anymore, but do I have it in C?
# it's possible I never had it

### map symmetric pair w/o diagonal

sub ij2sn { # ij to Symmetric Non-diagonal
    # maps $i, $j into $t in a triangle matrix of any size!
    # maps contiguously if we do: i in 0..total, j in 0..i-1 (excludes diagonal), good for loops so we don't have to calculate this, but we do need to calculate it to map back
    # in current form, assumes i>=j
    my ($i, $j) = @_;
    return $j + $i*($i-1)/2;
}

sub sn2ij {
    # reverse map of ij2sn
    my ($k) = @_;
    ## I had to resort to the quadratic equation to get $i from $k:
    my $i = int( (1 + sqrt(1+8*$k))/2 ); # correct $i is rounded-down solution!
    my $j = $k - $i*($i-1)/2; # j is just the remainder!
    return ($i, $j); # this is it!
}


### may asymmetric pair

sub ij2m {
    # like &t, but maps i,j into a full matrix of any size!
    # see &ij2m_test for a graphical representation of this mapping
    my ($i, $j) = @_;
    return $i > $j ? $i**2+$j : ($j+1)**2-$i-1;
}

sub m2ij {
    # reverse map of ij2m, needed to figure out what a pair was from the encoded version
    my ($m) = @_;
    # size of matrix is bound by the square root of this number (this is exactly the value of the larger of i and j)
    my $d = int sqrt $m; # round sqrt down
    my $r = $m - $d**2; # find the difference between the input and the square matrix fit (which is smaller because it was rounded *down* to the nearest square), so this $r is positive or zero
    if ($r <= $d) { return ($d, $r); } # if this difference is smaller than the dimension, these two values correspond exactly to i and j
    else { return ($d*2-$r, $d); } # in the other case, j is the largest of the two, and i is the difference along the length of both sides (so d*2 instead of just d; no -1,etc because this way boundaries are met correctly)
}

#############
### TESTS ###
#############

### symmetric w/ diagonal 

sub t_test {
    # code that makes a nice STDOUT-printed matrix that illustrates the ij2m mapping
    # input a matrix size to know how far to take the example
    my $d = shift;
    print join("\t", '', 0..$d-1)."\n";
    for (my $i = 0; $i < $d; $i++) { # y-axis
	my @row = ($i);
	for (my $j = 0; $j <= $i; $j++) { # x-axis
	    push @row, t($i, $j);
	}
	print join("\t", @row)."\n";
    }
}

### symmetric w/o diagonal

sub ij2sn_test {
    # code that makes a nice STDOUT-printed matrix that illustrates the ij2m mapping
    # input a matrix size to know how far to take the example
    my $d = shift;
    print join("\t", '', 0..$d-1)."\n";
    for (my $i = 0; $i < $d; $i++) { # y-axis
	my @row = ($i);
	for (my $j = 0; $j < $i; $j++) { # x-axis
	    push @row, ij2sn($i, $j);
	}
	print join("\t", @row)."\n";
    }
}

sub ij2sn_testRev {
    # code that makes a nice STDOUT-printed matrix that illustrates the ij2m mapping
    # input a matrix size to know how far to take the example
    my $d = shift;
    print join("\t", '', 0..$d-1)."\n";
    for (my $i = 0; $i < $d; $i++) { # y-axis
	my @row = ($i);
	for (my $j = 0; $j < $i; $j++) { # x-axis
	    push @row, join ',', sn2ij(ij2sn($i, $j));
	}
	print join("\t", @row)."\n";
    }
}

### asymmetric w/ diagonal

sub ij2m_test {
    # code that makes a nice STDOUT-printed matrix that illustrates the ij2m mapping
    # input a matrix size to know how far to take the example
    my $d = shift;
    print join("\t", '', 0..$d-1)."\n";
    for (my $j = 0; $j < $d; $j++) { # y-axis
	my @row = ($j);
	for (my $i = 0; $i < $d; $i++) { # x-axis
	    push @row, ij2m($i, $j);
	}
	print join("\t", @row)."\n";
    }
}

sub ij2m_testRev {
    # code that makes a nice STDOUT-printed matrix that illustrates the ij2m mapping
    # input a matrix size to know how far to take the example
    my $d = shift;
    print join("\t", '', 0..$d-1)."\n";
    for (my $j = 0; $j < $d; $j++) { # y-axis
	my @row = ($j);
	for (my $i = 0; $i < $d; $i++) { # x-axis
	    push @row, join ',', m2ij(ij2m($i, $j));
	}
	print join("\t", @row)."\n";
    }
}

1;

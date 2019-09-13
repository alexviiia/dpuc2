<img src="logo.png" alt="dPUC2" align="right" />

# dPUC2: Domain Prediction Using Context

The dPUC (Domain Prediction Using Context) software improves domain prediction by considering every domain in the context of other domains, rather than independently as standard approaches do. 
Our framework maximized the probability of the data under an approximation, which reduces to a pairwise context problem, and we have shown that our probabilistic method is indeed more powerful than competing methods.

Here you can download our code, and learn how to use it.
This page is the software's manual.

## Version 2 brief history

The dPUC 2.xx series is a major update from dPUC 1.0. 
dPUC2 works with HMMER 3 and Pfam 24 onward, and is much faster than before.
If you were a dPUC 1.0 user, note that inputs and outputs are completely different; the new setup is simpler due to the changes that HMMER3 and the new Pfams have brought. 
Lastly, there are improvements in the context network parametrization, most notably the addition of directed context. 
See `Changes` file for more information.


# Installation

dPUC2 installation is a bit complex, broken down into the following parts below:

- Install software dependencies (easy on Fedora and Ubuntu Linux)
- Clone this repository
- Download Pfam HMM database files
- Download precomputed context network file

## Dependencies

dPUC2 contains code written in C that requires the `lp_solve` library headers to compile.
dPUC2 requires:

- Perl 5
- `gzip`
- The [Inline::C](http://search.cpan.org/~etj/Inline-C-0.62/lib/Inline/C.pod) Perl package
- The [HMMER3](http://hmmer.janelia.org/) protein hidden Markov model software
- The [lp_solve 5.5](http://lpsolve.sourceforge.net/5.5/) integer linear programming library
- Git versioning control software (optional, to clone repository)

Perl and gzip are part of practically all Linux distributions.
Below are commands for installing the rest of these dependencies on two common Linux distributions: Fedora and Ubuntu.

### Fedora Linux

On a terminal, type:
```bash
sudo dnf install perl-Inline-C lpsolve-devel hmmer git
```
This command should also work on Red Hat Linux too (untested).

### Ubuntu Linux

On a terminal, type:
```bash
sudo apt install libinline-c-perl liblpsolve55-dev hmmer git
sudo ln -s /usr/lib/lp_solve/liblpsolve55.so /usr/lib/liblpsolve55.so
```
This worked on Ubuntu 19.04.
This command may also work on other Debian Linux systems (untested).

If the above doesn't prevent an error about linking to the lp_solve library, other users reported these additional steps helping:

- Adding `/usr/lib/lp_solve/` to the LD library path in `~/.bashrc`
- Creating a file in `/etc/ld.so.conf.d` called `lpsolve.conf` with the content: `/usr/lib/lp_solve`

Thank you to Philipp Meister for helping me figure out these Ubuntu instructions.

### Other systems

Besides the above dependencies, dPUC2 has to be able to find `lp_lib.h`.
`DpucLpSolve.pm` has a hardcoded path for `/usr/include/lpsolve`, which works for Fedora and Ubuntu.
On other systems, you can run this command to make sure `lp_lib.h` is in the right place:
```bash
find / -name lp_lib.h
```
If the result includes `/usr/include/lpsolve/lp_lib.h`, then no changes to `DpucLpSolve.pm` are needed.
Otherwise, open `DpucLpSolve.pm` with a text editor and replace `/usr/include/lpsolve` with the directory that contains `lp_lib.h`.

Email me if you enconter any issues installing dPUC2 and its dependencies, especially if you were able to solve them, so I can update these instructions.


## Cloning repository

You can download a ZIP copy of this repository on GitHub.
However, cloning the repository has many advantages.
You can clone it with this command:
```bash
git clone https://github.com/alexviiia/dpuc2.git
```
If there are any updates in the future, you can simply update your copy of the repository by typing, inside the `dpuc2` directory:
```bash
git pull
```

## Running scripts from arbitrary directories

Each script can be run easily from the directory that contains it (and its `*.pm` Perl module dependencies).
For example, if you cloned this repository onto the local directory `dpuc2/`, on a terminal you can run one of the help messages by typing:
```bash
cd dpuc2/ 
perl -w 1dpuc2.pl # <ARGS>... 
```

To run the scripts from other directories, you have to specify the directory containing the `*.pm` Perl modules using the `-I` option for `perl`, like so:
```bash
cd ..
perl -Idpuc2/ -w dpuc2/1dpuc2.pl # <ARGS>...
```
Note that the home directory shortcut `~` doesn't work with `-I`, but you can use `$HOME` instead.
So if the code is in `~/dpuc2/`, then this will work:
```bash
perl -I$HOME/dpuc2/ -w ~/dpuc2/1dpuc2.pl # <ARGS>...
```

## Pfam database

For the HMM database, you will need to download `Pfam-A.hmm.gz` and `Pfam-A.hmm.dat.gz` from [Pfam](ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/). 
Use HMMER3's hmmpress to prepare database for searching.
A minimal set of commands is:
```bash
# download files
wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz
wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.dat.gz
# uncompress Pfam-A.hmm, while keeping compressed version around
gunzip -k Pfam-A.hmm.gz
# press the database
hmmpress Pfam-A.hmm
# remove temporary uncompressed copy (no longer needed)
rm Pfam-A.hmm
```

## Precomputed dPUC2 context network file

Lastly, download a precomputed context network file that corresponds with the Pfam version you want to use, such as `dpucNet.pfam32.txt.gz`, from the [dpuc2-data](https://github.com/alexviiia/dpuc2-data) repository. 
A quick download command is (change the pfam version accordingly):
```bash
wget https://github.com/alexviiia/dpuc2-data/raw/master/dpucNet.pfam32.txt.gz
```
Code to recompute this file from Pfam is also provided (slow and uses lots of memory, not recommended).

# Running the scripts

## Sample files

This repository contains the following sample files:

Sample input:

- [`sample/sample.fa`](sample/sample.fa): 10 random protein sequences from the *Plasmodium falciparum* proteome.

Sample outputs:

- [`sample/sample.txt`](sample/sample.txt): HMMER3/Pfam32 domain predictions using weak filters (produced by `0runHmmscan.pl` as below).
- [`sample/sample.dpuc.txt`](sample/sample.dpuc.txt): Dpuc2 domain predictions using the previous sample as input (produced by `1dpuc2.pl` as in the first command below).

## Synopsis of scripts

All scripts give detailed usage instructions when executed without arguments.
The following commands can be run on the sample file placed in the same directory as the code and called from that location.
The input file may be compressed with `gzip` and may be specified with or without the GZ extension.
The output file is compressed with `gzip`, whether the outputs indicate a GZ extension or not.
So the commands as they are below produce and will work entirely with compressed files, without a single GZ extension indicated.

Do this the first time only, to compile the C portion of the code or to troubleshoot dependency problems. 
If successful, it will display the "usage" message.
```bash
perl -w 1dpuc2.pl
```

Create the dPUC context count network for your Pfam version (skip this step by downloading above the precomputed dPUC network you need).
```bash
perl -w dpucNet.pl Pfam-A.full.uniprot dpucNet.txt 
```

Produce domain predictions with HMMER3 (provides hmmscan) with weak filters. 
This is the slowest step.
```bash
perl -w 0runHmmscan.pl hmmscan Pfam-A.hmm sample.fa sample.txt 
```

(Optional) compress output, our scripts will read it either way 
```bash
gzip sample.txt 
```

Produce dPUC predictions with default parameters 
```bash
perl -w 1dpuc2.pl Pfam-A.hmm.dat dpucNet.txt sample.txt sample.dpuc.txt 
```

Produce dPUC predictions for more than one p-value threshold for candiate domains (produces one output for each threshold)
```bash
perl -w 1dpuc2.pl Pfam-A.hmm.dat dpucNet.txt sample.txt sample.dpuc.txt --pvalues 1e-1 1e-4 1e-9
```


## Using `dpucNet.pl`: Extract the context count network from Pfam predictions

This is the help message you get by running the script without arguments:
```bash
perl -w dpucNet.pl
```
```
# dpucNet.pl: Extract the context count network from Pfam predictions
# dPUC 2.08 - https://github.com/alexviiia/dpuc2
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

Usage: perl -w dpucNet.pl <Pfam-A.full.uniprot> <dpucNet output>

The required inputs are
    <Pfam-A.full.uniprot>  Input "full" domain prediction file from Pfam (29-latest).
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
```

Let me further insist that, despite optimizations, this script consumes prohibitive amounts of memory (and time, to a lesser extent), so most regular users won't even be able to complete a run.
First, downloading the input Pfam file `Pfam-A.full.uniprot.gz` takes a very long time (16 GB for Pfam 32).
Then, for Pfam 32, 15G of memory are needed for the script to run.

The memory problems are due to Pfam-A.full.uniprot listing domains grouped by family, whereas these family pair counts can only be calculated after predictions are grouped by protein, so all the data must be in memory before it can be rearranged.

The output file starts with one line listing all Pfam accession codes observed (sorted), followed by one line for every pair of Pfam families (also sorted), with the count of observed pairs.
For Pfam 32, this file looks a bit like this (abbreviated):
```
PF00001 PF00002 PF00003 PF00004 PF00005 ... PF18867 PF18868 PF18869 PF18870 PF18871
PF00001 PF00001 1992
PF00001 PF00002 1
PF00001 PF00003 2
PF00001 PF00004 2
PF00001 PF00005 1
...
PF18870 PF18870 1
PF18871 PF04326 3
PF18871 PF04471 3
PF18871 PF18735 10
PF18871 PF18871 5
```


## Using `0runHmmscan.pl`: Get domain predictions from your protein sequences

This is the help message you get by running the script without arguments:
```bash
perl -w 0runHmmscan.pl
```
```
# 0runHmmscan.pl: Get domain predictions from your protein sequences
# dPUC 2.08 - https://github.com/alexviiia/dpuc2
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

Usage: perl -w 0runHmmscan.pl <hmmscan> <Pfam-A.hmm> <FASTA input> <output table> \
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
```

Note that the p-value threshold of 1e-4 here matches the default domain candidate p-value threshold that dPUC uses, and which worked best in our benchmarks.

Run on the sample file provided here, which is a standard FASTA file that looks like this (abbreviated):
```
>MAL8P1.134
MSSIAKKTQYNIKVDIHEVK...
>PFA0015c
MALKKGVINESKLSARNVLE...
...
```
The main output is a large table that looks like this (abbreviated and split into 2 columns):
```
#                                                                            --- full sequence --- 
# target name        accession   tlen query name           accession   qlen   E-value  score  bias 
#------------------- ---------- ----- -------------------- ---------- ----- --------- ------ ----- 
C2                   PF00168.30   103 MAL8P1.134           -           1704   2.3e-39  120.2   0.8 
C2                   PF00168.30   103 MAL8P1.134           -           1704   2.3e-39  120.2   0.8 
C2                   PF00168.30   103 MAL8P1.134           -           1704   2.3e-39  120.2   0.8 
C2                   PF00168.30   103 MAL8P1.134           -           1704   2.3e-39  120.2   0.8 
C2                   PF00168.30   103 MAL8P1.134           -           1704   2.3e-39  120.2   0.8 
Ferlin_C             PF16165.5    154 MAL8P1.134           -           1704   1.9e-08   20.6   0.0 
ATS                  PF15445.6    446 PFA0015c             -           1327  3.1e-194  632.5  39.3 
Duffy_binding        PF05424.11   182 PFA0015c             -           1327   7.5e-67  211.0  19.8 
Duffy_binding        PF05424.11   182 PFA0015c             -           1327   7.5e-67  211.0  19.8 
NTS                  PF15447.6     37 PFA0015c             -           1327   5.5e-14   38.8   0.2 
ArsR                 PF09824.9    159 PFA0015c             -           1327   3.5e-05    9.6  15.0 
...
```
```
-------------- this domain -------------   hmm coord   ali coord   env coord
  #  of  c-Evalue  i-Evalue  score  bias  from    to  from    to  from    to  acc description of target
--- --- --------- --------- ------ ----- ----- ----- ----- ----- ----- ----- ---- ---------------------
  1   5   3.6e-10   3.6e-10   26.4   0.7     2    90    11   105    10   118 0.80 C2 domain
  2   5     7e-08     7e-08   19.0   0.0     3    57   177   234   175   279 0.72 C2 domain
  3   5   3.5e-10   3.5e-10   26.4   0.0     7    91   670   761   666   773 0.77 C2 domain
  4   5   5.6e-08   5.6e-08   19.3   0.0     4    95  1333  1434  1330  1441 0.78 C2 domain
  5   5   1.4e-07   1.4e-07   18.0   0.0     3   102  1494  1599  1492  1600 0.82 C2 domain
  1   1   1.9e-08   1.9e-08   20.6   0.0    63   143  1616  1702  1587  1704 0.82 Ferlin C-terminus
  1   1  3.1e-194  3.1e-194  632.5  39.3     2   446   895  1327   894  1327 0.95 acidic terminal segme...
  1   2     2e-45     2e-45  141.2   0.7     1   156   120   267   120   289 0.91 Duffy binding domain
  2   2   4.8e-25   4.8e-25   74.8   5.2     1   162   528   657   528   684 0.84 Duffy binding domain
  1   1   5.5e-14   5.5e-14   38.8   0.2     1    37    14    48    14    48 0.98 N-terminal segments o...
  1   1   2.8e-05   2.8e-05    9.9   1.0    68   143   756   832   743   842 0.86 ArsR transcriptional ...
...
```
Although the file still calls the various columns by "E-value", they are all are p-values when output by `0runHmmscan.pl`.

## Using `1dpuc2.pl`: Produce dPUC domain predictions from raw hmmscan data

This is the help message you get by running the script without arguments:
```bash
perl -w 1dpuc2.pl
```
```
# 1dpuc2.pl: Produce dPUC domain predictions from raw hmmscan data
# dPUC 2.08 - https://github.com/alexviiia/dpuc2
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

Usage: perl -w 1dpuc2.pl [-options] <Pfam-A.hmm.dat> <dpucNet> <input table> <output table>

The required inputs are
    <Pfam-A.hmm.dat>   Pfam annotation file, contains GA thresholds and nesting network.
    <dpucNet>          Directed context network of domain family pair counts.
    <input table>      The output from hmmscan (previous script).
    <output table>     Input with most domains removed by dPUC. Format is identical to input.
                       File path will be used as the base to generate multiple files, one for 
                       each p-value thresholds, if more than one threshold is requested.

Options:
    --pvalues <x...>   The p-value thresholds for candidate domains  [1e-4]
                       If multiple thresholds are provided, the output table has these 
                       thresholds inserted as text just before the file extension, creating 
                       separate outputs for each threshold.
    --alpha <x>        Pseudocount (sym-Dirichlet param) for context scores  [0.01]
    --cutNet <x>       Keep context network counts >= cutNet  [2]
    --fCut <x>         Permissive overlap max proportion  [.50]
    --lCut <x>         Permissive overlap max length in amino acids  [40]
    --timeout <x>      lp_solve timeout in seconds  [1]
    --noGzip           Outputs will not be compressed (default is gzip compression).

If more than one p-value threshold is provided on candidate domains:
Note that dPUC predictions are not necessarily nested as candidate domain thresholds are varied, 
since dPUC solves a complicated pairwise optimization problem, so these files are not redundant 
and are not obtained by filtering the largest file. Each output is as if dPUC had been run anew
on the given p-value threshold, although the code does this efficiently by factoring out shared 
steps between thresholds.

Input file may be compressed with gzip, and may be specified with or without the .gz 
extension.  Output file will be automatically compressed with gzip unless --noGzip is used.
```

# Original Perl module descriptions

Since my source code is not self-documented, here's a brief description of what each module does.

| File              | Description                                                              |
|-------------------|--------------------------------------------------------------------------|
| FileGz.pm         | Handles normal and compressed files transparently.                       |
| ParsePfam.pm      | Parses Pfam-A.hmm.dat and other Pfam-specific files.                     |
| Domains.pm        | Processes domains, particularly overlaps.                                |
| Hmmer3ScanTab.pm  | Runs and parses hmmscan outputs.                                         |
| Dpuc.pm           | Main dPUC package connects various context domain prediction strategies. |
| DpucPosElim.pm    | C code that solves the dPUC "positive elimination".                      |
| DpucLpSolve.pm    | Perl-C glue for `lp_solve` library.                                      |
| DpucNet.pm        | Counts domain family pairs in Pfam-A.full.                               |
| DpucNetScores.pm  | Turns context pair counts into bitscores for domain prediction.          |
| DpucOvsCompact.pm | Compacts domain overlap LP definitions by finding cliques.               |
| NetCC.pm          | Finds connected components in a network.                                 |
| EncodeIntPair.pm  | Maps non-negative integer pairs into single integers.                    |


# Software license

This code is released under the GNU GPLv3 (GNU General Public License version 3).
See [LICENSE](LICENSE).

# Compatibility

This code has been tested on

- Perl 5.18, 5.20, 5.22, 5.28
- HMMER 3.0, 3.1b1, 3.1b2, 3.2.1
- Pfam 25, 27-32
- Inline::C 0.53, 0.62, 0.78
- lp_solve 5.5.2.0

# Citations

2017-04-12.
Alejandro Ochoa, Mona Singh.
Domain prediction with probabilistic directional context.
Bioinf 33(16) 2471-2478. 
[Article](http://dx.doi.org/10.1093/bioinformatics/btx221),
[bioRxiv](http://biorxiv.org/content/early/2016/12/14/094284) 2016-12-14.

2015-11-17.
Alejandro Ochoa, John D Storey, Manuel Llinás, and Mona Singh.
Beyond the E-value: stratified statistics for protein domain prediction.
PLoS Comput Biol. 11 e1004509.
[Article](http://dx.doi.org/10.1371/journal.pcbi.1004509),
[arXiv](http://arxiv.org/abs/1409.6384) 2014-09-23.

2011-03-31.
Alejandro Ochoa, Manuel Llinás, and Mona Singh.
Using context to improve protein domain identification.
BMC Bioinformatics, 12:90.
[Article](http://dx.doi.org/10.1186/1471-2105-12-90).


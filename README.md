dPUC 2
===

Domain Prediction Using Context.
Extend your Pfam predictions without loss of precision using domain context!

About dPUC 2
===

The dPUC (Domain Prediction Using Context) software improves domain prediction by considering every domain in the context of other domains, rather than independently as standard approaches do. Our specific framework maximized the probability of the data assuming a naive Bayes approximation, which reduces to a pairwise context problem, and we have shown that our probabilistic method is indeed more powerful than competing methods.

The dPUC 2.xx series is a major update from dPUC 1.0. For the average user, what matters is that now dPUC works with HMMER 3 and Pfam 24 onward, and dPUC is much faster than before. If you were a dPUC 1.0 user, you should know that inputs and outputs are completely different; the new setup is simpler due to the changes that HMMER3 and the new Pfams have brought. Lastly, there are improvements in the context network parametrization, most notably the addition of directed context. See the release notes for more information.

Installation
===

You'll need Perl 5, the [Inline::C](http://search.cpan.org/~etj/Inline-C-0.62/lib/Inline/C.pod) Perl package, [HMMER3](http://hmmer.janelia.org/), the [lp_solve 5.5 library](http://lpsolve.sourceforge.net/5.5/) and gzip installed.

For the HMM database, you will need to download Pfam-A.hmm.gz and Pfam-A.hmm.dat.gz from [Pfam](ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/). Use HMMER3's hmmpress to prepare database for searching.

Lastly, you may want to download a precomputed context network file, such as dpucNet.pfam27.txt.gz, from [viiia.org](http://viiia.org/dpuc2/?l=en-us). But code to recompute this file from Pfam is also provided (slow, not recommended).

Synopsis of scripts
===

Do this the first time only, to compile the C portion of the code or to troubleshoot dependency problems. If successful, it will display the "usage" message.
```
perl -w 1dpuc2.pl
```

Create the dPUC context count network for your Pfam version (skip this step by downloading above the precomputed dPUC network you need).
```
perl -w dpucNet.pl Pfam-A.full dpucNet.txt 
```

Produce domain predictions with HMMER3 (provides hmmscan) with weak filters. This is the slowest step.
```
perl -w 0runHmmscan.pl hmmscan Pfam-A.hmm sample.fa sample.txt 
```

(Optional) compress output, our scripts will read it either way 
```
gzip sample.txt 
```

Produce dPUC predictions with default parameters 
```
perl -w 1dpuc2.pl Pfam-A.hmm.dat dpucNet.txt sample.txt sample.dpuc.txt
```

More details
===

All scripts give detailed usage instructions when executed without arguments.  I have more detailed recommendations, usage examples (code snippets and test sample inputs and outputs) at my personal website, viiia.org, in [English](http://viiia.org/dpuc2/?l=en-us) and [Spanish](http://viiia.org/dpuc2/).


Citation
===

Alejandro Ochoa, Manuel Llinás, and Mona Singh. Using context to improve protein domain identification. BMC Bioinformatics 2011, 12:90. 2011-03-31. [Pubmed](http://www.ncbi.nlm.nih.gov/pubmed/21453511). [Article](http://dx.doi.org/10.1186/1471-2105-12-90).

Alejandro Ochoa, John Storey, Manuel Llinás, and Mona Singh.  "Beyond the E-value: stratified statistics for protein domain prediction." Submitted. [arXiv preprint](http://arxiv.org/abs/1409.6384).


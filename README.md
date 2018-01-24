# hmmfreq
## HMM approach to reconstruct pooled haplotypes in known concentration from short sequencing reads

---

Prequisite:

MGLM-master; package from https://github.com/Hua-Zhou/MGLM

SamToMsa; package from https://github.com/thomaskf/SamToMsa

---

##  How to cite


---

##  Documentation

Below are the recommended values for the parameters of the main function, hmmfreq_dirmn2D(). 

* rdpair, first output file from converting a SAM file to an alignment file with SamToMsa
* rdsingle, second output file from converting a SAM file to an alignment file with SamToMsa
* siz, analysis window size in nucleotides (50)
* ovlap, overlap between consecutive windows as a percentage (0.99)
* prewinsiz, analysis window size for the initial estimate of coverage (400)
* outfile, path to write a fasta output file, set to 0 to not write the output (0)
* smoothing, a Fourier transform local smoothing of the coverage can be used by setting this option to 1 at the cost of slowing down the reconstruction (0)
* smoothing_factor, threshold for local smoothing (2)
* prewinshape, shape of the analysis window for local coverage estimation ('triangle_around')
* verbose, print information messages and some plots of the results (0)
* error_count, percentage of error that is allowed in the coverage estimates (0.05)
* substate, optionnally limit the number of states to observations (0)
* cnt_per_state, record the coverage for each state (1)
* transit_error, number of mutations threshold to forbid transitions between windows, must be defined in accordance with the overlap length ovlap (1)
* freq_limit, lower limit for coverage of the local windows to keep a subsequence (0.01)
* bothend, if 1 run the HMM from both ends and generate two sequences for each haplotype (0)
* true_freq, value for the expected frequencies in the pool ([.625 .25 .125])

##  Example

Assuem the short reads have been mapped to a areference sequence and the resulting SAM file has been converted using SamToMsa tool (outputs of SamToMsa being file.paired.txt and file.single.txt).

```matlab
addpath('../hmmfreq');
addpath(genpath('../MGLM-master'));
options = struct('rdpair', 'file.paired.txt',...
     'rdsingle', 'file.single.txt',...
     'siz', 50, 'ovlap', 0.99, 'prewinsiz', 400, 'outfile', 0, 'smoothing', 0, 'smoothing_factor', 2,...
     'prewinshape', 'triangle_around', 'verbose', 0, 'error_count', 0.05, 'substate', 0, 'cnt_per_state', 1,...
     'transit_error', 1, 'freq_limit', 0.01, 'bothend', 0, 'true_freq', target_freq) ;
haplo = hmmfreq_dirmn2D(options) ;
fastawrite('reconstructed_haplotypes.fasta',haplo) ; % using the Matlab Bioinformatics Toolbox
```

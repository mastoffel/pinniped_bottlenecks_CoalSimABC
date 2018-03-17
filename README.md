# Analysis workflow part 1 for "Recent demographic histories and genetic diversity across pinnipeds are shaped by anthropogenic interactions and mediated by ecology and life-history"

## Coalescent-simulations and ABC analysis.

The scripts contained in this repository were used to simulate genetic data under 
a bottleneck and a neutral demographic scenario using `strataG` as an interface
to [fastsimcoal2](http://cmpg.unibe.ch/software/fastsimcoal2/). These data
were then used for ABC analyses across all 29 pinniped species. 


## Prerequisites

(1) These script should be run on a multi-core machine, optimally using around 20
or more cores. However, everything can run quickly for testing purposes based
on a small number of simulations (say 1000 instead of 20000000)

(2) Install fastsimcoal2 and check using the [strataG](https://cran.r-project.org/web/packages/strataG/index.html) vignette
that the package can access fastsimcoal.

(3) Several other packages need to be installed, which are mentioned in the scripts.
Among them is a small package specifically written for the analysis of this paper,
the `sealABC` package, which can be install from GitHub with:

```
devtools::install_github("mastoffel/sealABC")
```

## Sequence of scripts

1) 1_sim_msats_simcoal.R    # simulates genetic data under the two scenarios
2) 2_comp_abc_1_simcoal.R   # First part of the ABC analysis
3) 3_comp_abc_2_simcoal.R   # Second part of the ABC analysis
4) 4_comp_abc_3_save.R      # Saving ABC results

## Analysis workflow part 2

The rest of the analyses can be found here:
https://github.com/mastoffel/pinniped_bottlenecks



# scTreeSim

Simulator for developmental processes on lineage trees.

``` text
         ┌── o
     ┌───┤
     │   └── o
 ───┤
     │   ┌── o
     └───┤
         └── o
```

> [!WARNING]
> **Work in progress.** scTreeSim is under active development and has not been released.
> Currently, **only the tree simulation is validated** (see [Tree simulation](#tree-simulation)).
> All other modules (see [Experimental modules](#experimental-modules-not-yet-validated)) are not yet fully tested or documented. 
> Their interfaces may change, and their output should not be relied on.
> A fully documented and validated release, accompanied by a dedicated publication, is planned.

## Installation

The development version (`main`) can be installed from GitHub:

``` r
# install.packages("remotes")
remotes::install_github("sccevo/scTreeSim")
```

To reproduce the tree simulations of the publication, install the version tagged `v0.0.1`:

``` r
remotes::install_github("sccevo/scTreeSim@v0.0.1")
```

## Tree simulation

In **scTreeSim**, phylogenetic trees are simulated under an age-dependent branching process with incomplete sampling (for details, see [Mulberry & Pilarski et al., 2026](https://doi.org/10.1093/molbev/msag080)). Note: The simulations in the publication were run with version [`v0.0.1`](https://github.com/sccevo/scTreeSim/releases/tag/v0.0.1) of this package, which was functional but not yet CRAN-conformant.

The user can specify

1. number $n$ of extant tips, or
2. time of origin $t_{or}$.

The simulation starts with a single cell. Cell lifetimes are sampled from the Gamma distribution with shape $k > 0$ and scale $\theta > 0$. 
At the end of its lifetime, each cell divides with probability $(1-d)$ or dies with probability $d$. 
The process continues until

1.  $\frac{n}{\rho}$ cells are alive at a time, where $\rho$ is the sampling probability of extant cells, or
2.  the time $t_{or}$ has passed.

Next, the lifetimes of surviving cells are censored to the present. 
Optionally, dead cells are pruned and extant cells are sampled uniformly with probability $\rho$.

``` r
library(scTreeSim)

# tree with 50 sampled extant tips, 20% sampling, 10% death probability
# (returns NULL with a message if the process goes extinct; simply re-run or change the seed)
tree <- sim_adb_ntaxa_samp(ntaxa = 50, shape = 3, mean_lifetime = 1,
                           death_prob = 0.1, sampling_prob = 0.2)
```

See `vignette("simulate_tree")` for more examples, including time-of-origin simulations and multi-type trees.

## Experimental modules (not yet validated)

The following features build on a simulated tree. They are included for transparency and early testing, but are **not yet validated**:

| Module | Main functions | Vignette |
|---|---|---|
| Gene expression | `sim_expression_ou()`, `sim_expression_matrix()`, `sim_expression_counts()` | `simulate_tree` |
| Dynamic lineage barcodes | `sim_barcode_seq()`, `sim_barcode_nonseq()`, `sim_barcode_generic()` | `simulate_barcodes` |
| Binary traits | `sim_binary()` | `simulate_binary` |
| Single-nucleotide variants | `sim_gt16_nt_diploid()` | `simulate_SNV` |

## Citation

If you use the tree simulation, please cite: 

Nicola Mulberry, Julia Pilarski, Jana Dinger, Tanja Stadler, Age-dependent Phylodynamics with Application to Single-cell Lineage Trees, Molecular Biology and Evolution, Volume 43, Issue 4, April 2026, msag080, https://doi.org/10.1093/molbev/msag080

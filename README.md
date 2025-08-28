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

## Tree simulation

In **scTreeSim**, phylogenetic trees are simulated under an age-dependent branching process with incomplete sampling. The user can specify

(1) number $n$ of extant tips, or
(2) time of origin $t_{or}$.

The simulation starts with a single cell. Cell lifetimes are sampled from the Gamma distribution with shape $k > 0$ and scale $\theta > 0$. At the end of its lifetime, each cell divides with probability $(1-d)$ or dies with probability $d$. The process continues until

1.  $\frac{n}{\rho}$ cells are alive at a time, where $\rho$ is the sampling probability of extant cells, or

2.  the time $t_{or}$ has passed.

Next, the lifetimes of surviving cells are censored to the present. Optionally, dead cells are pruned and extant cells are sampled uniformly with probability $\rho$.

# Pneumococcus NFDS modelling

This is a private repo for MAGPIE's pneumococcus modelling.


## Reproducing paper results

Each paper lives in `paper/<paper>/` and has its own pipeline, Docker image and README.
To build all results for one paper, run from the repo root:

```bash
docker compose build <paper>
docker compose run --rm <paper>
```

Results are written to `projects/<paper>/results/`.

| Paper | Folder |
|---|---|
| South Africa NFDS vaccine paper | [`projects/sa_nfds_paper`](projects/sa_nfds_paper) |


## Data
The data/ folder has USA-Mass, SouthAfrica and Invasiveness folders.
The Invasiveness folder contains the serotype-specific invasiveness values from Nick Croucher's
2022 paper (Disease progression rates... ).
The other two contain the data for Massachusetts and South Africa.

## Code
The code/ folder contains scripts that group and reweight the original data to create
"the NFDS model" (for South Africa).

## An NFDS Model

### Model specification

We refer to the information needed to define the NFDS dynamics (for a location) as
"the NFDS model" or "the model" (for that location).

For now, we refer to the gpsc-st groups as "strains".

The information needed to simulate an NFDS model forward in time is the model, and some parameters. The model consists of:

1. The matrix G. $G_{jl}$ is the fraction of isolates in strain $j$ that have the gene at locus $l$.
+ The columns of G are in order of most stable loci first.
+ The rows of G are "strains": for example, a Serotype (ST) and GPSC combination.
2. The serotype for each strain (this is also in the strain name, given in the rownames of G).
3. el (the equilibrium frequencies).
4. x0 (the starting values of the populations for each strain).

For the South Africa data (published in Bellman et al.),
this is in the form of a list in a file called samodel.rds, in the data/ folder. This object also includes `delta`, the "Delta statistic" (from Corander et al 2017). Low values of delta correspond to small changes in frequency between pre- and post-vaccination populations.

### Model parameters
The parameters consist of:

1. The vaccine composition (which strains are in the vaccine).
2. The strength of vaccination ($v \in [0,1]$, a constant; if a strain is in the vaccine, its growth is inhibited, by a amount governed by this parameter).
3. The strength of the NFDS effect, $\sigma$ (positive, usually in $[0,1]$ but not necessarily).
4. The portion $p$ of the $L$ loci that are under strong NFDS. We order the loci by increasing delta statistic; the first $pL$ of them are considered to be under strong NFDS.
5. The NFDS weights, $w$. We'll use $w_{high}$ and $w_{low}$. If a gene is under "strong NFDS", its locus has the higher weight, and otherwise, it has the lower weight. The lower weight may be 0.
6. A migration parameter $m$, allowing strain loss and gain from the rest of the world.

## Model simulation

There is a file called simulate_ode.R which does what it says: simulates the ode version of the NFDS dynamics.

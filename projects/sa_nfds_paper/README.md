# South Africa NFDS vaccine paper

The below workflow description is WIP of the current state of this...

To start an interactive R session for the project run :

```bash
docker compose run --rm sa_nfds_paper R
```


## Create the results

With Docker (recommended):

```bash
docker compose build sa_nfds_paper
docker compose run --rm sa_nfds_paper
```

Without Docker, from the repo root in R:

```r
Sys.setenv(TAR_PROJECT = "sa_nfds_paper")
targets::tar_make()
```

## Outputs

- `results/plots/`: figures
- `results/computed_data/`: tables and data behind the figures

## Workflow

See the [pipeline diagram](pipeline.md) for how the targets fit together.

After adding or removing targets, regenerate it with
`Rscript projects/sa_nfds_paper/update_pipeline_diagram.R`.

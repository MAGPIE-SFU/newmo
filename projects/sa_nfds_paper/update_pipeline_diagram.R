# Regenerates pipeline.md. Run from the repo root after adding or removing targets:
#   Rscript projects/sa_nfds_paper/update_pipeline_diagram.R
Sys.setenv(TAR_PROJECT = "sa_nfds_paper")
writeLines(
  c("```mermaid", targets::tar_mermaid(targets_only = TRUE), "```"),
  "projects/sa_nfds_paper/pipeline.md"
)

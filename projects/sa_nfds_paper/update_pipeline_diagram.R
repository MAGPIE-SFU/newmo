# Regenerates pipeline.md. Run from the repo root after adding or removing targets:
#   Rscript projects/sa_nfds_paper/update_pipeline_diagram.R
Sys.setenv(TAR_PROJECT = "sa_nfds_paper")
writeLines(
  c("```mermaid", targets::tar_mermaid(targets_only = TRUE), "```"),
  "projects/sa_nfds_paper/pipeline.md"
)

# The below generates the target worflow network as a pdf
# net <- targets::tar_visnetwork(
#   targets_only     = TRUE,   # leave out the functions; much more readable
#   level_separation = 250     # horizontal space between the columns of targets
# )
# net$width  <- "100%"
# net$height <- "95vh"         # 95% of the browser window's height

# htmlwidgets::saveWidget(net, "projects/sa_nfds_paper/pipeline.html",
#                         selfcontained = TRUE, title = "sa_nfds_paper pipeline")

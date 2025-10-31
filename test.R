##############################
## KSM worker access diag.  ##
##############################

require("doFuture")           # foreach via futures
require("future.batchtools")  # Slurm via batchtools
require("foreach")            # foreach frontend
require("future")             # future core (for safety)
require("fs")                 # filesystem utilities
require("parallel")           # detectCores, etc.

# where the *driver* sees ksm (may be "")
ksm_lib_hint <- dirname(system.file(package = "ksm"))

# make sure our user lib is in path on the driver
Sys.setenv(R_LIBS_USER = "~/R/library")
.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))

cat("Driver .libPaths():\n"); print(.libPaths())
cat("Driver sees ksm_lib_hint: '", ksm_lib_hint, "'\n", sep = "")

##############################
## Packages to load on worker
##############################
libraries_to_load <- c(
  "doFuture",        # registers foreach adaptor
  "future",          # future core
  "future.batchtools",
  "foreach",
  "parallel",
  "fs"
  # NOTE: intentionally NOT loading "ksm" here
)

##############################
## Exports to workers
##############################
vars_to_export <- c(
  "ksm_lib_hint"
)

##################
## Hyper-params ##
##################
RR <- 1:4  # a few test tasks

##################
## Slurm resources
##################
cores_per_node <- 1
resources_list <- list(
  cpus_per_task = cores_per_node,
  mem = "2G",
  walltime = "00:05:00",
  nodes = 1
)

###############################
## Plan: Slurm + multisession
###############################
# IMPORTANT: Your Slurm template should export R_LIBS_USER before Rscript:
#   export R_LIBS_USER="$HOME/R/library"
#   export R_LIBS="$R_LIBS_USER"

myslurm <- tweak(
  batchtools_slurm,
  template  = "batchtools.slurm.iid2d.tmpl",  # use your existing template
  resources = resources_list
)

registerDoFuture()
plan(list(myslurm, multisession))

#####################
## Worker diagnostic
#####################
path <- getwd()
start_time <- Sys.time()

res <- foreach(
  r = RR,
  .combine  = "rbind",
  .export   = vars_to_export,
  .packages = libraries_to_load
) %dopar% {
  # --- worker context ---
  host <- tryCatch(system("hostname", intern = TRUE), error = function(e) NA_character_)
  r_ver <- paste(R.version$major, R.version$minor)
  r_libs_user <- Sys.getenv("R_LIBS_USER", unset = NA_character_)
  lp_before <- .libPaths()

  # Try load ksm the "normal" way
  have_ns1 <- requireNamespace("ksm", quietly = TRUE)
  ver1 <- if (have_ns1) as.character(utils::packageVersion("ksm")) else NA_character_
  path1 <- if (have_ns1) system.file(package = "ksm") else NA_character_

  # If not found, try with driver's hint (prepend its dir to .libPaths)
  if (!have_ns1 && nzchar(ksm_lib_hint)) {
    .libPaths(c(ksm_lib_hint, .libPaths()))
  end_if <- NULL
  }
  have_ns2 <- requireNamespace("ksm", quietly = TRUE)
  ver2 <- if (have_ns2) as.character(utils::packageVersion("ksm")) else NA_character_
  path2 <- if (have_ns2) system.file(package = "ksm") else NA_character_

  lp_after <- .libPaths()
  search_path <- paste(search(), collapse = " | ")

  data.frame(
    r = r,
    host = host,
    r_version = r_ver,
    R_LIBS_USER = r_libs_user,
    libpaths_before = paste(lp_before, collapse = ";"),
    libpaths_after  = paste(lp_after,  collapse = ";"),
    ksm_found_initial = have_ns1,
    ksm_version_initial = ver1,
    ksm_path_initial = path1,
    ksm_found_after_hint = have_ns2,
    ksm_version_after_hint = ver2,
    ksm_path_after_hint = path2,
    search_path = search_path,
    stringsAsFactors = FALSE
  )
}

# Save & print
outfile <- file.path(path, "ksm_worker_diag.csv")
write.csv(res, outfile, row.names = FALSE)
print(res)

elapsed_time_minutes <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
cat("Elapsed time (min):", round(elapsed_time_minutes, 2), "\n")
cat("Wrote: ", outfile, "\n", sep = "")

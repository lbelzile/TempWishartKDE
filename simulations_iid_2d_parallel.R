# Simulation study
#
# 1. generate from six distributions @f a total of @n samples of positive definite symmetric matrices of size @d
# 2. evaluate the optimal bandwidth using 'ksm:::bandwidth_optim' using one of the five available methods @m
# 3. Calculate the integrated square error

require("cubature")           # adaptive multidimensional integration (e.g., Cuhre)
require("doFuture")           # registers future-based backends for foreach (%dopar%)
require("fs")                 # cross-platform filesystem ops (paths, files, dirs)
require("future.batchtools")  # runs futures via batchtools (Slurm templates/resources)
require("ksm")                # SPD kernel density tools (kdens_symmat, integrate_spd, bandwidth_optim)
require("parallel")           # base R parallelism (PSOCK clusters, detectCores)

# Define the list of libraries to load on each cluster node

libraries_to_load <- c(
  "cubature",
  "doFuture",
  "fs",
  "future.batchtools",
  "parallel"
)

# Define the list of variables/functions to export to the worker nodes

vars_to_export <- c(
  "combo",
  "cores_per_node",
  "criteria",
  "ISE",
  "IAE",
  "kernels",
  "libraries_to_load",
  "models",
  "nobs",
  "path",
  "replication",
  "resources_list",
  "RR",
  "setup_parallel_cluster",
  "time_min"
)

# Sets up a parallel cluster, loads necessary libraries, and exports required variables globally

setup_parallel_cluster <- function() {
  num_cores <<- detectCores()
  cl <<- makeCluster(num_cores)
  
  # Export the list of libraries to the worker nodes
  clusterExport(cl, varlist = "libraries_to_load")
  
  # Load necessary libraries on each cluster node
  invisible(clusterEvalQ(cl, {
    lapply(libraries_to_load, library, character.only = TRUE)
  }))
  
  # Export all necessary objects, functions, and parameters to the worker nodes
  clusterExport(cl, varlist = vars_to_export)
  
  return(cl) # Return the cluster object
}

# Initialize all variables in the list as NULL except vars_to_export and setup_parallel_cluster

invisible(
  lapply(
    vars_to_export[
      !(vars_to_export %in% c("vars_to_export", "setup_parallel_cluster"))
    ],
    function(x) assign(x, NULL, envir = .GlobalEnv)
  )
)

# Hyper-parameters

nobs <- c(125L)
models <- 1:6
combo <- 1:5
kernels <- c("smnorm", "smlnorm", "Wishart")
criteria <- c("lscv", "lcv")
RR <- 1:128 # replications

#' Integrated squared error of kernel density estimator for symmetric matrices
#'
#' Given a sample \code{x} from a target density, and the optimal bandwidth, compute the integrated squared error via numerical integration.
#' @param S matrix at which to evaluate
#' @param x sample of symmetric matrix observations from which to build the kernel density kernel
#' @param bandwidth double for the bandwidth of the kernel
#' @param model number of the target function from which points are drawn.
#' @param kernel string, one of \code{Wishart}, \code{smlnorm} (log-Gaussian) or \code{smnorm} (Gaussian).
#' @param return integrated squared error
#' @param ... additional parameters, currently ignored

ISE <- function(
    S,
    x,
    bandwidth,
    model = 1:6,
    kernel = c("Wishart", "smlnorm", "smnorm"),
    ...
) {
  #model <- match.arg(model, choices = 1:6)
  model <- as.integer(match.arg(as.character(model), choices = as.character(1:6)))
  kernel <- match.arg(kernel)
  dim <- ncol(x)
  # Compute the squared difference between the kernel and the target density
  (ksm::kdens_symmat(
    x = S,
    xs = x,
    b = bandwidth,
    kernel = kernel,
    log = FALSE
  ) -
      ksm::simu_fdens(S, model = model, d = dim))^2
}

#' Integrated absolute error of kernel density estimator for symmetric matrices
#'
#' Given a sample \code{x} from a target density, and the optimal bandwidth, compute the integrated squared error via numerical integration.
#' @param S matrix at which to evaluate
#' @param x sample of symmetric matrix observations from which to build the kernel density kernel
#' @param bandwidth double for the bandwidth of the kernel
#' @param model number of the target function from which points are drawn.
#' @param kernel string, one of \code{Wishart}, \code{smlnorm} (log-Gaussian) or \code{smnorm} (Gaussian).
#' @param return integrated squared error
#' @param ... additional parameters, currently ignored

IAE <- function(
    S,
    x,
    bandwidth,
    model = 1:6,
    kernel = c("Wishart", "smlnorm", "smnorm"),
    ...
) {
  #model <- match.arg(model, choices = 1:6)
  model <- as.integer(match.arg(as.character(model), choices = as.character(1:6)))
  kernel <- match.arg(kernel)
  dim <- ncol(x)
  # Compute the absolute difference between the kernel and the target density
  abs(
    ksm::kdens_symmat(
      x = S,
      xs = x,
      b = bandwidth,
      kernel = kernel,
      log = FALSE
    ) -
      ksm::simu_fdens(S, model = model, d = dim)
  )
}

##################
## Set the path ##
##################

path <- getwd()
# setwd(path)

##############################
## Parallelization on nodes ##
##############################

cores_per_node <- 64 # number of cores for each node in the super-computer

resources_list <- list(
  cpus_per_task = cores_per_node,
  mem = "240G",
  walltime = "24:00:00",
  nodes = 1
  # Omit 'partition' to let SLURM choose
)

###############################
## Main code (exact version) ##
###############################

.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))

# Disable the check for random number generation misuse in doFuture
options(doFuture.rng.onMisuse = "ignore")

# Register the doFuture parallel backend
registerDoFuture()

# Tweak the batchtools_slurm with the custom template and resources
myslurm <- tweak(
  batchtools_slurm,
  template = "batchtools.slurm.iid2d.tmpl",
  resources = resources_list
)

# Set the plan for future
plan(list(myslurm, multisession))

# Create empty data frames to store the results
raw_results <- data.frame(
  replication = integer(),
  nobs = integer(),
  model = integer(),
  combo = integer(),
  kernel = integer(),
  criterion = character(),
  logISE = numeric(),
  # logIAE = numeric(),
  bandwidth = numeric(),
  time_min = numeric(),
  stringsAsFactors = FALSE
)

# Capture the start time
start_time <- Sys.time()

# Split RR into node-sized blocks
RR_blocks <- split(RR, ceiling(seq_along(RR) / cores_per_node))

# For each over blocks (one Slurm job per block)
res <- foreach(
  block = RR_blocks,
  .combine = "rbind",
  .export = vars_to_export,
  .packages = libraries_to_load
) %dopar% {
  
  # spin up a per-node PSOCK cluster using your existing helper
  cl <- setup_parallel_cluster()
  
  # prepare (i,j) grid once
  combinations <- expand.grid(
    i = seq_along(nobs),
    j = seq_along(models),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  # make the grid visible on workers (minimal)
  parallel::clusterExport(cl, varlist = "combinations", envir = environment())
  
  block_results <- data.frame(
    replication = integer(),
    nobs = integer(),
    model = integer(),
    combo = integer(),
    kernel = character(),
    criterion = character(),
    logISE = numeric(),
    # logIAE = numeric(),
    bandwidth = numeric(),
    time_min = numeric(),
    stringsAsFactors = FALSE
  )
  
  # each r in this block runs in parallel across the node's cores; (i,j) loop runs inside each worker
  r_list <- parallel::parLapply(cl, X = block, fun = function(r_local) {
    set.seed(r_local)
    
    out_rows_all <- vector("list", nrow(combinations))
    for (idx in seq_len(nrow(combinations))) {
      i <- combinations$i[[idx]]
      j <- combinations$j[[idx]]
      
      xs <- ksm::simu_rdens(n = nobs[i], model = models[j], d = 2L)
      
      local_rows <- vector("list", length(combo))
      for (k in combo) {
	t0 <- Sys.time() # <-- start timing

        band <- ksm::bandwidth_optim(
          x = xs,
          criterion = criteria[k %% 2L + 1L],
          kernel = kernels[k %% 3L + 1L]
        )
        ISE_int <- try(
          ksm::integrate_spd(
            f = function(S) {
              ISE(
                S,
                x = xs,
                kernel = kernels[k %% 3L + 1L],
                bandwidth = band,
                model = j
              )
            },
            dim = 2L,
            method = "cuhre",
            lb = 1e-8,
            ub = Inf,
            neval = 1e7L
          ),
          silent = TRUE
        )
        if (!inherits(ISE_int, "try-error")) {
          log_ISE <- ISE_int$integral
        } else {
          log_ISE <- NA
        }
        # IAE_int <- try(
        #   ksm::integrate_spd(
        #     f = function(S) {
        #       IAE(
        #         S,
        #         x = xs,
        #         kernel = kernels[k %% 3L + 1L],
        #         bandwidth = band,
        #         model = j
        #       )
        #     },
        #     dim = 2L,
        #     method = "cuhre",
        #     lb = 1e-8,
        #     ub = Inf,
        #     neval = 1e7L
        #   ),
        #   silent = TRUE
        # )
        # if (!inherits(IAE_int, "try-error")) {
        #   log_IAE <- IAE_int$integral
        # } else {
        #   log_IAE <- NA
        # }

	elapsed_min <- as.numeric(difftime(Sys.time(), t0, units = "mins")) # <-- stop timing

        local_rows[[k]] <- data.frame(
	  replication = r_local,
          nobs = nobs[i],
          model = models[j],
	  combo = combo[k],
          kernel = kernels[k %% 3L + 1L],
          criterion = criteria[k %% 2L + 1L],
          logISE = log_ISE,
          # logIAE = log_IAE,
          bandwidth = band,
	  time_min = elapsed_min,
          stringsAsFactors = FALSE
        )
      }
      out_rows_all[[idx]] <- do.call(rbind, local_rows)
    }
    
    do.call(rbind, out_rows_all)
  })
  
  block_results <- rbind(block_results, do.call(rbind, r_list))
  
  # tear down the per-node cluster
  try(parallel::stopCluster(cl), silent = TRUE)
  
  # return results for this block
  block_results
}

# Combine results from all nodes
raw_results <- res

# Stop parallel execution
plan(sequential)

# Calculate the duration in minutes
elapsed_time_minutes <- as.numeric(difftime(
  Sys.time(),
  start_time,
  units = "mins"
))
print(paste("Elapsed time:", round(elapsed_time_minutes, 2), "minutes"))

# Save the raw results to a CSV file
raw_output_file <- file.path(path, "raw_ISE_iid_2d.csv")
write.csv(raw_results, raw_output_file, row.names = FALSE)

print("Raw results saved to raw_ISE_iid_2d.csv")

#########################
## Process the results ##
#########################

# Build a summary table by (nobs, model, kernel, criterion)

# Guard against accidental factor coercion
raw_results$kernel    <- as.character(raw_results$kernel)
raw_results$criterion <- as.character(raw_results$criterion)

# Split by grouping keys
.grp <- interaction(raw_results$nobs, raw_results$model,
                    raw_results$kernel, raw_results$criterion,
                    drop = TRUE)

summ_list <- lapply(split(raw_results, .grp), function(df) {
  data.frame(
    nobs      = df$nobs[1],
    model     = df$model[1],
    kernel    = df$kernel[1],
    criterion = df$criterion[1],
    # ISE stats
    mean_ISE   = mean(df$logISE, na.rm = TRUE),
    sd_ISE     = sd(df$logISE, na.rm = TRUE),
    median_ISE = median(df$logISE, na.rm = TRUE),
    IQR_ISE    = IQR(df$logISE, na.rm = TRUE),
    # # IAE stats
    # mean_IAE   = mean(df$logIAE, na.rm = TRUE),
    # sd_IAE     = sd(df$logIAE, na.rm = TRUE),
    # median_IAE = median(df$logIAE, na.rm = TRUE),
    # IQR_IAE    = IQR(df$logIAE, na.rm = TRUE),
    mean_time_min = mean(df$time_min, na.rm = TRUE),
    # # Bandwidth stats (useful to see selection variability)
    # mean_bandwidth   = mean(df$bandwidth, na.rm = TRUE),
    # sd_bandwidth     = sd(df$bandwidth, na.rm = TRUE),
    # median_bandwidth = median(df$bandwidth, na.rm = TRUE),
    # IQR_bandwidth    = IQR(df$bandwidth, na.rm = TRUE),
    # Count of successful runs (non-NA ISE)
    # n_eff_ise = sum(!is.na(df$logISE)),
    # n_eff_iae = sum(!is.na(df$logIAE)),
    stringsAsFactors = FALSE
  )
})

summary_results <- do.call(rbind, summ_list)
row.names(summary_results) <- NULL

# Save the summary to CSV next to the raw file
summary_output_file <- file.path(path, "summary_iid_2d.csv")
write.csv(summary_results, summary_output_file, row.names = FALSE)

print("Summary results saved to summary_iid_2d.csv")


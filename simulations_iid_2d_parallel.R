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
  #"ksm",
  "parallel"
)

# Define the list of variables/functions to export to the worker nodes

vars_to_export <- c(
  "ISE",
  "IAE",
  "cores_per_node",
  "libraries_to_load",
  "path",
  "resources_list",
  "RR",
  "setup_parallel_cluster",
  "vars_to_export",
  "nobs",
  "models",
  "combo",
  "kernels",
  "criteria",
  "ksm_lib"
)

# Sets up a parallel cluster, loads necessary libraries, and exports required variables globally

setup_parallel_cluster <- function() {
  num_cores <<- detectCores() - 1
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
models <- 1:1
combo <- 1:5
kernels <- c("smnorm", "smlnorm", "Wishart")
criteria <- c("lscv", "lcv")
RR <- 1:2

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
  walltime = "14:00:00",
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
  nobs = integer(),
  model = integer(),
  kernel = integer(),
  criterion = character(),
  logISE = numeric(),
  bandwidth = numeric(),
  stringsAsFactors = FALSE
)

# Capture the start time
start_time <- Sys.time()

# Parallel loop over the replications (RR), each node processes one set of RR values
res <- foreach(
  r = RR,
  .combine = "rbind",
  .export = vars_to_export,
  .packages = libraries_to_load
) %dopar%
  {
    # Set a unique seed for each node (replication)
    set.seed(r)

    # Set library paths within each worker node
    .libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))

    local_raw_results <- data.frame(
      nobs = integer(),
      model = integer(),
      kernel = integer(),
      criterion = character(),
      logISE = numeric(),
      bandwidth = numeric(),
      stringsAsFactors = FALSE
    )

    for (i in seq_along(nobs)) {
      for (j in seq_along(models)) {
        # set.seed(b + 2025 * i)
        xs <- ksm::simu_rdens(n = nobs[i], model = models[j], d = 2L)
        for (k in seq_along(combo)) {
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
          IAE_int <- try(
            ksm::integrate_spd(
              f = function(S) {
                IAE(
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
          if (!inherits(IAE_int, "try-error")) {
            log_IAE <- IAE_int$integral
          } else {
            log_IAE <- NA
          }
          local_raw_results <- rbind(
            local_raw_results,
            data.frame(
              nobs = nobs[i],
              model = models[j],
              kernel = kernels[k %% 3L + 1L],
              criterion = criteria[k %% 2L + 1L],
              logISE = log_ISE,
              logIAE = log_IAE,
              bandwidth = band
            )
          )
        }
      }
    }

    # Return the raw results for this replication
    return(local_raw_results)
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


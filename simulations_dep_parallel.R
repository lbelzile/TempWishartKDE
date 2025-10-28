require(ksm)

libraries_to_load <- c(
  "cubature",
  "doFuture",
  "future",
  "future.batchtools",
  "ksm",
  "parallel"
)


#' Simulation study serial dependence simulation from WAR(1)
#' @param n sample size
#' @param Mmod model for autoregressive coefficients matrix
#' @param Smod covariance matrix of VAR components
#' @param K degrees of freedom
#' @return a cube of dimension 2 by 2 by \code{n}.
simu_rWAR <- function(
  n,
  Mmod = c("M1", "M2", "M3"),
  Smod = c("S1", "S2", "S3"),
  K = 1L
) {
  Mmod <- match.arg(Mmod)
  Smod <- match.arg(Smod)
  n <- as.integer(n)
  stopifnot(n >= 1)
  if (Mmod == "M1") {
    M <- matrix(
      c(0.9, 0, 1, 0),
      byrow = TRUE,
      nrow = 2,
      ncol = 2
    )
  } else if (Mmod == "M2") {
    M <- matrix(
      c(0.3, -0.3, -0.3, 0.3),
      byrow = TRUE,
      nrow = 2,
      ncol = 2
    )
  } else {
    M <- diag(rep(0.5, 2))
  }
  if (Smod == "S1") {
    S <- diag(2)
  } else if (Smod == "S2") {
    S <- matrix(0.5, 2, 2) +
      diag(rep(0.5, 2))
  } else {
    S <- matrix(
      c(1, 0.9, 0.9, 1),
      byrow = TRUE,
      nrow = 2,
      ncol = 2
    )
  }
  ksm::rWAR(
    n = n,
    M = M,
    S = S,
    K = K
  )
}
#' Integrated squared error of kernel density estimator for symmetric matrices
#'
#' Given a sample \code{x} from a target density, and the optimal bandwidth, compute the integrated squared error via numerical integration.
#' @param x sample of symmetric matrix observations from which to build the kernel density kernel
#' @param bandwidth double for the bandwidth of the kernel
#' @param fdens density of the target function from which points are drawn.
#' @param kernel string, one of \code{Wishart}, \code{smlnorm} (log-Gaussian) or \code{smnorm} (Gaussian).
#' @param tol double for tolerance of numerical integral
#' @param lb lower bound for integration range of eigenvalues
#' @param ub upper bound for integration range of eigenvalues
#' @param return double with value of the integrated squared error
#' @param ... additional parameters, currently ignored
logISE <- function(
  x,
  bandwidth,
  fdens,
  kernel = c("Wishart", "smlnorm", "smnorm"),
  tol = 1e-3,
  lb = 0,
  ub = Inf,
  neval = 1e6L,
  method = c("suave", "hcubature"),
  ...
) {
  method <- match.arg(
    method[1],
    choices = c("hcubature", "pcubature", "cuhre", "divonne", "suave", "vegas")
  )
  kernel <- match.arg(kernel)
  dim <- ncol(x)
  stopifnot(
    dim %in% c(2L, 3L),
    length(lb) == 1L,
    length(ub) == 1L,
    lb >= 0,
    ub > lb,
    tol > 0,
    is.numeric(bandwidth) & length(bandwidth) == 1L
  )
  integrand <- function(vars, kernel) {
    # Construct SPD matrix S
    if (length(vars) == 3L) {
      S <- ksm::rotation_scaling(vars[1], vars[2:3])
      jacobian_value <- abs(vars[2] - vars[3]) / 4
    } else {
      S <- ksm::rotation_scaling(vars[1:3], vars[4:6])
      jacobian_value <- abs(vars[4] - vars[5]) *
        abs(vars[4] - vars[6]) *
        abs(vars[5] - vars[6]) *
        sin(vars[2]) /
        24
    }
    S <- array(S, dim = c(nrow(S), ncol(S), 1))

    # Compute the squared difference between the kernel and the target density
    diff_squared <- (ksm::kdens_symmat(
      x = S,
      xs = x,
      b = bandwidth,
      kernel = kernel,
      log = FALSE
    ) -
      fdens(S))^2

    # Return the product of diff_squared and the Jacobian factor
    return(diff_squared * jacobian_value)
  }
  if (dim == 2L) {
    # Define the limits of integration
    lower_limit <- c(0, rep(lb, length.out = dim))
    upper_limit <- c(2 * pi, rep(ub, length.out = dim))
  } else if (dim == 3L) {
    lower_limit <- c(rep(0, 3), rep(lb, length.out = dim))
    upper_limit <- c(2 * pi, pi, 2 * pi, rep(ub, length.out = dim))
  }
  result <- cubature::cubintegrate(
    f = integrand,
    lower = lower_limit,
    upper = upper_limit,
    method = "suave",
    relTol = 1e-3,
    maxEval = neval,
    kernel = kernel
  )

  # print(result)
  # Return the result of the integral
  return(log(result$integral))
}


vars_to_export <- c(
  "logISE",
  "simu_rWAR"
  #TODO figure out what needs to be added here, beyond locally defined functions
)

# Sets up a parallel cluster, loads necessary libraries, and exports required variables globally

setup_parallel_cluster <- function() {
  num_cores <<- parallel::detectCores() - 1
  cl <<- parallel::makeCluster(num_cores)

  # Export the list of libraries to the worker nodes
  parallel::clusterExport(cl, varlist = "libraries_to_load")

  # Load necessary libraries on each cluster node
  invisible(parallel::clusterEvalQ(cl, {
    lapply(libraries_to_load, library, character.only = TRUE)
  }))

  # Export all necessary objects, functions, and parameters to the worker nodes
  list_vars <- ls()
  parallel::clusterExport(cl, varlist = ls())

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

##################
## Set the path ##
##################

path <- getwd()
# setwd(path)

################
## Parameters ##
################

# Variables for the simulation
nrep <- 1000L


##############################
## Parallelization on nodes ##
##############################

resources_list <- list(
  cpus_per_task = cores_per_node,
  mem = "240G",
  walltime = "60:00:00",
  nodes = 1
  # Omit 'partition' to let SLURM choose
)


# Variables for the simulation
nrep <- 1000L
cores_per_node <- 63 # number of cores for each node in the super-computer
nobs <- c(125L, 250L)
Mmod <- c("M1", "M2", "M3")
Smod <- c("S1", "S2", "S3")
kernels <- c("Wishart", "smlnorm")


###############################
## Main code (exact version) ##
###############################

.libPaths("~/R/library")

# Disable the check for random number generation misuse in doFuture
options(doFuture.rng.onMisuse = "ignore")

# Register the doFuture parallel backend
doFuture::registerDoFuture()

# Tweak the batchtools_slurm with the custom template and resources
myslurm <- future::tweak(
  batchtools_slurm,
  template = "batchtools.slurm.dependent.tmpl",
  resources = resources_list
)

# Set the plan for future
future::plan(list(myslurm, multisession))

# Create empty data frames to store the results
raw_results <- data.frame(
  nobs = integer(),
  Mmod = integer(),
  Smod = integer(),
  kernel = character(),
  logISE = numeric(),
  bandwidth = numeric(),
  timing = numeric()
)

# Capture the start time
start_time <- Sys.time()

# Parallel loop over the replications (RR), each node processes one set of RR values
res <- foreach::foreach(
  r = RR,
  .combine = "rbind",
  .export = vars_to_export,
  .packages = libraries_to_load
) %dopar%
  {
    # Set a unique seed for each node (replication)
    set.seed(r)
    nobs <- c(125L, 250L)
    Mmod <- c("M1", "M2", "M3")
    Smod <- c("S1", "S2", "S3")
    kernels <- c("Wishart", "smlnorm")

    # Set library paths within each worker node
    .libPaths("~/R/library")

    local_raw_results <- data.frame(
      nobs = integer(),
      Mmod = integer(),
      Smod = integer(),
      kernel = character(),
      logISE = numeric(),
      bandwidth = numeric(),
      timing = numeric()
    )
    for (i in seq_along(nobs)) {
      for (j in seq_along(Mmod)) {
        for (k in seq_along(Smod)) {
          xs <- simu_rWAR(
            n = nobs[i],
            Mmod = Mmod[j],
            Smod = Smod[k]
          )
          for (l in seq_along(kernels)) {
            start_time <- Sys.time()
            band <- ksm:::bandwidth_optim(
              data = xs,
              criterion = "lscv",
              h = ceiling(nobs[i]),
              kernel = kernels[l]
            )

            log_ISE <- try(
              logISE(
                x = xs,
                bandwidth = band,
                fdens = function(x) {
                  simu_fdens(x, model = model, 2)
                },
                kernel = kernels[k],
                lb = 0,
                ub = Inf,
                tol = 1e-3
              ),
              silent = TRUE
            )
            if (inherits(log_ISE, "try-error")) {
              log_ISE <- NA
            }
            timing <- as.numeric(difftime(
              Sys.time(),
              start_time,
              units = "seconds"
            ))
            raw_results <- rbind(
              raw_results,
              data.frame(
                nobs = nobs[i],
                Mmod = j,
                Smod = k,
                kernel = kernels[l],
                logISE = log_ISE,
                bandwidth = band,
                timing = timing
              )
            )
          }
        }
      }
    }

    # Return the raw results for this replication
    return(local_raw_results)
  }

# Combine results from all nodes
raw_results <- res

# Stop parallel execution
future::plan(sequential)

# Calculate the duration in minutes
elapsed_time_minutes <- as.numeric(difftime(
  Sys.time(),
  start_time,
  units = "mins"
))
print(paste("Elapsed time:", round(elapsed_time_minutes, 2), "minutes"))

# Save the raw results to a CSV file
raw_output_file <- file.path(path, "raw_ISE_dep.csv")
write.csv(raw_results, raw_output_file, row.names = FALSE)

print("Raw results saved to raw_ISE_dep.csv")

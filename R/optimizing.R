

# ========================================================================*
# Helper functions ----
# ========================================================================*

# Get the value of the objective at the best set of parameters found.
get_val <- function(.op, .pkg) {
    if (.pkg %in% c("stats", "nloptr")) {
        val <- .op$value
    } else {
        val <- .op$fval
    }
    return(val)
}
# Get convergence code
get_converge_code <- function(.op, .pkg) {
    if (.pkg == "minqa") {
        conv_code <- .op$ierr
    } else if (.pkg == "nloptr") {
        conv_code <- .op$convergence
    } else {
        conv_code <- .op$convergence
    }
    return(conv_code)
}
# Get whether an optimization has converged.
get_converged <- function(.op, .pkg) {
    if (.pkg == "minqa") {
        conv <- .op$ierr == 0L
    } else if (.pkg == "nloptr") {
        conv <- .op$convergence > 0L
    } else {
        conv <- .op$convergence == 0L
    }
    return(conv)
}
# Combine arguments and run optimizer
run_optim <- function(.optim, .pars, .fn, .control, .fn_args) {
    .args <- c(list(par = .pars, fn = .fn, control = .control), .fn_args)
    # Adjust arguments across packages:
    if (packageName(environment(.optim)) == "nloptr") {
        if (!is.null(.args[["control"]][["maxit"]])) {
            .args[["control"]][["maxeval"]] <- .args[["control"]][["maxit"]]
            .args[["control"]][["maxit"]] <- NULL
        }
        if (!is.null(.args[["control"]][["reltol"]])) {
            .args[["control"]][["xtol_rel"]] <- .args[["control"]][["reltol"]]
            .args[["control"]][["reltol"]] <- NULL
        }
        .args[["x0"]] <- .pars
        .args[["par"]] <- NULL
    }
    if (packageName(environment(.optim)) == "minqa") {
        if (!is.null(.args[["control"]][["maxit"]])) {
            .args[["control"]][["maxfun"]] <- .args[["control"]][["maxit"]]
            .args[["control"]][["maxit"]] <- NULL
        }
        # No obvious equivalent in `minqa`, so just remove this:
        .args[["control"]][["reltol"]] <- NULL
    }
    op <- do.call(.optim, .args)
    return(op)
}



# ========================================================================*
# Main function ----
# ========================================================================*





#' Winnowing optimization
#'
#' @details
#' All output files are tables.
#' The output for the box evaluations has the box number (column `box`),
#' parameter values (named based on bounds arguments or `par1`, `par2`, ...),
#' and output from objective function at those values (`val`).
#' The output for optimizations include the
#' starting parameter values (with `_start` suffix),
#' ending parameter values (with `_end` suffix),
#' convergence code for optimization (`conv_code`), and
#' output from objective function ending parameter values (`val`).
#'
#'
#'
#' @param fn Objective function to minimize.
#' @param lower_bounds Lower bounds of boxes for each parameter.
#'     If this vector has names, these names are used for column names in
#'     intermediate output tables (argument `file_names` below).
#'     If names are present, they must match those in `upper_bounds`,
#'     but don't need to be in the same order.
#' @param upper_bounds Upper bounds of boxes for each parameter.
#'     See `lower_bounds` for providing names for this vector.
#' @param fn_args List containing other arguments to use for `fn`.
#'     Defaults to `list()`.
#' @param n_bevals Number of evaluations of `fn` per box. Defaults to `100L`.
#' @param n_boxes Number of boxes. Defaults to `1000L`.
#' @param n_outputs A numeric vector indicating the number of best-fitting
#'     optimization outputs from each optimization step to pass to the next one.
#'     Because the first step is the optimizations for the best evaluation per
#'     box, the first item in this vector must be `<= n_boxes`.
#'     Because each optimization step is taking a number of fits from the
#'     previous step, each successive item in this vector must be less than
#'     or equal to the previous.
#'     The length of this vector must be `>= 1` and
#'     must match the lengths of `controls`,
#'     `optimizers`, and, if requested, `file_names`.
#'     Defaults `c(100L, 20L, 3L)`.
#' @param controls A list where each item is a named list.
#'     Each list in the top-most list contains arguments to use for
#'     the `control` argument for the optimization that occurs for that step.
#'     The first step is the optimization for each box.
#'     The length of this vector must be `>= 1` and  must match the lengths
#'     of `n_outputs`, `optimizers`, and, if requested, `file_names`.
#'     The default argument for this is
#'     `list(list(maxit = 100, reltol = 1e-4), list(maxit = 500, reltol = 1e-6), list(maxit = 1000, reltol = 1e-8))`.
#'      This results in three optimization steps that get increasingly
#'      polished.
#' @param optimizers A list of optimizer functions to use for each optimization step.
#'     The length of this vector must be `>= 1` and  must match the lengths
#'     of `n_outputs`, `controls`, and, if requested, `file_names`.
#'     The default argument for this is `c(optim, optim, optim)`.
#'     This results in three optimization steps using `stats::optim`.
#' @param file_names A character vector specifying the file name(s)
#'     where to save intermediate output.
#'     The first output is from box evaluations, and the subsequent ones
#'     are from all the optimizations except for the last one.
#'     To output all the results from the last optimizer step,
#'     change the last item of the `n_outputs` argument.
#'     If provided, this vector must be the same length as
#'     `n_outputs`, `controls`, and `optimizers`.
#'     If you want some output to be written and others not to be, just
#'     set the file names for the step(s) you don't want written to `NA`.
#'     For example, for the default 3 optimizations, if
#'     `file_names = c(NA, "file1.txt", "file2.txt")`, the first
#'     two optimization steps will be written to files, but the box
#'     evaluations will not be written.
#'     File names that end with `.txt` will be tab-delimited,
#'     and those that end with `.csv` will be comma-delimited.
#'     See `Details` above for info on the output from these files.
#'     Note that an error will trigger if attempting to overwrite an existing
#'     file unless the `overwrite` argument is set to `TRUE`.
#'     If `NULL` (the default), no output is written.
#' @param overwrite A single logical for whether to allow overwriting files
#'     for intermediate output (the `file_names` argument to this function).
#'     Defaults to `FALSE`.
#' @param na_stop Single logical for whether to return the matrix of initial
#'     evaluations in each box if it contains `NA`s.
#'     If `FALSE`, the optimizer ignores these values and continues on
#'     (unless they're all `NA`).
#'     Defaults to `FALSE`.
#'
#' @importFrom stats optim
#' @importFrom tools file_ext
#'
#' @returns A list containing `tail(n_outputs, 1)` object(s) of the class
#'     returned by the last optimization step.
#'
#'
#' @examples
#' # "Continuous location planning problem with Manhattan metric"
#' # From https://ds-pl-r-book.netlify.app/optimization-in-r.html
#' fn <- function(loc, a, x, y) sum(a * (abs(x - loc[1]) + abs(y - loc[2]) ) )
#' n <- 100
#' a.vec <- sample(1:100, size = n)    # sample weights for each point/customer
#' x.vec <- rnorm(n)                   # sample x coordinates
#' y.vec <- rnorm(n)                   # sample y coordinates
#'
#' res <- winnowing_optim(fn, lower_bounds = rep(-1, 2), upper_bounds = rep(1, 2),
#'                        fn_args = list(a = a.vec, x = x.vec , y = y.vec))
#'
#' @export
#'
winnowing_optim <- function(fn,
                            lower_bounds,
                            upper_bounds,
                            fn_args = list(),
                            n_bevals = 100L,
                            n_boxes = 1000L,
                            n_outputs = c(100L, 20L, 3L),
                            controls = list(list(maxit = 100, reltol = 1e-4),
                                            list(maxit = 500, reltol = 1e-6),
                                            list(maxit = 1000, reltol = 1e-8)),
                            optimizers = c(optim, optim, optim),
                            file_names = NULL,
                            overwrite = FALSE,
                            na_stop = FALSE) {

    # Type and length checks:
    stopifnot(is.function(fn))
    stopifnot(is.numeric(lower_bounds) && is.numeric(upper_bounds))
    stopifnot(length(lower_bounds) == length(upper_bounds))
    stopifnot(is.list(fn_args))
    stopifnot(length(n_bevals) == 1L && as.integer(n_bevals) == n_bevals)
    stopifnot(length(n_boxes) == 1L && as.integer(n_boxes) == n_boxes)
    stopifnot(is.list(controls))
    stopifnot(length(controls) >= 1)
    stopifnot(all(sapply(controls, is.list)))
    stopifnot(is.list(optimizers))
    stopifnot(length(optimizers) == length(controls))
    stopifnot(all(sapply(optimizers, is.function)))
    stopifnot(is.numeric(n_outputs))
    stopifnot(length(n_outputs) == length(controls))
    stopifnot(min(n_outputs) < n_boxes)
    stopifnot(all(diff(n_outputs) < 0))
    stopifnot(length(overwrite) == 1L && inherits(overwrite, "logical"))
    stopifnot(length(na_stop) == 1L && inherits(na_stop, "logical"))

    # Value checks:
    stopifnot(all(lower_bounds < upper_bounds))
    stopifnot(n_bevals > 0L)
    stopifnot(n_boxes > 0L)
    stopifnot(all(n_outputs > 0L))

    # Number of optimization steps (NOT including box evaluations):
    n_optims <- length(controls)

    # Checking for names in one or more `_bounds` objects:
    if ((!is.null(names(lower_bounds)) && is.null(names(upper_bounds))) ||
        (is.null(names(lower_bounds)) && !is.null(names(upper_bounds)))) {
        stop("If `lower_bounds` is named, `upper_bounds` should be too, ",
             "and vice versa")
    }
    # If they are named, make sure they're the same and in the same order,
    # and use these for `par_names` (used for output) if present
    if (!is.null(names(lower_bounds))) {
        par_names <- names(lower_bounds)
        if (!identical(sort(par_names), sort(names(upper_bounds)))) {
            stop("Names for `lower_bounds` and `upper_bounds` must be the same ",
                 "(order doesn't matter)")
        }
        upper_bounds <- upper_bounds[par_names]
    } else par_names <- paste0("par", 1:length(lower_bounds))

    # Check validity of file names:
    if (is.null(file_names)) {
        file_names <- rep(NA_character_, n_optims)
    } else stopifnot(is.character(file_names) && length(file_names) == length(controls))
    for (i in 1:length(file_names)) {
        file_i <- file_names[[i]]
        if (!is.na(file_i)) {
            if (length(file_i) != 1 || !is.character(file_i)) {
                stop(sprintf("The item `file_names[[%i]]` ('%s') is not NA or a single string",
                             i, file_i))
            }
            if (! file_ext(file_i) %in% c("txt", "csv")) {
                stop(sprintf("Extension for `file_names[[%i]]` ('%s') is not 'txt' or 'csv'",
                             i, file_i))
            }
            if (!dir.exists(dirname(file_i))) {
                stop(sprintf("Directory for `file_names[[%i]]` ('%s') does not exist",
                             i, file_i))
            }
            if (!overwrite && file.exists(file_i)) {
                stop(paste0("File for `file_names[[", i, "]]` ('", file_i,
                            "') already exists and `overwrite` is `FALSE`"))
            }
        }
    }

    pkgs <- lapply(optimizers, \(o) packageName(environment(o)))
    if (!all(pkgs %in% c("stats", "minqa", "nloptr"))) {
        stop(paste("this function only programmed for stats::optim and",
                   "optimizers from minqa and nloptr packages"))
    }

    mids <- (upper_bounds + lower_bounds) / 2
    steps <- (upper_bounds - mids) / n_boxes
    n_pars <- length(mids)

    # ---------------------------------------------------------*
    # Box evaluations (potentially writing output):
    # ---------------------------------------------------------*
    # Objective function evaluations for box `i`
    bevals_i <- matrix(0.0, n_bevals, n_pars+2)
    colnames(bevals_i) <- c("box", par_names, "val")
    # Best set of parameters for each box:
    best_box_bevals <- matrix(0.0, n_boxes, n_pars+1)
    colnames(best_box_bevals) <- c(par_names, "val")
    for (i in 1:n_boxes) {
        for (j in 1:n_bevals) {
            .pars <- runif(n_pars, mids - steps * i, mids + steps * i)
            .val <- do.call(fn, c(list(.pars), fn_args))
            bevals_i[j,] <- c(i, .pars, .val)
        }
        if (any(is.na(bevals_i[,"val"])) && na_stop) {
            warning(sprintf(paste("\nThere were NAs in the %ith box.",
                                  "The matrix of parameter values and",
                                  "evaluations is being returned."), i))
            return(bevals_i)
        }
        if (all(is.na(bevals_i[,"val"]))) {
            warning(sprintf(paste("\nThe %ith box was all NAs!",
                                  "The matrix of parameter values and",
                                  "evaluations is being returned."), i))
            return(bevals_i)
        }
        if (!is.na(file_names[[1]])) {
            write.table(bevals_i, file_names[[1]], append = i > 1, quote = FALSE,
                        sep = ifelse(file_ext(file_names[[1]]) == "txt", "\t", ","),
                        row.names = FALSE, col.names = i == 1)
        }
        best_idx <- which(bevals_i[,"val"] == min(bevals_i[,"val"], na.rm = TRUE))[[1]]
        best_box_bevals[i,] <- bevals_i[best_idx, c(par_names, "val")]
    }
    # sort with lowest (best fit) at the top:
    best_box_bevals <- best_box_bevals[order(best_box_bevals[,"val"]),]

    # ---------------------------------------------------------*
    # Optimizations and optionally write files
    # ---------------------------------------------------------*
    # Number of optimizations for each step:
    n_optims_vec <- c(n_boxes, head(n_outputs, -1))
    # List to store outputs from optimizations
    optim_outs <- rep(list(NA), n_optims)
    # Now do optimizations themselves:
    for (k in 1:length(optim_outs)) {
        n <- n_optims_vec[[k]]
        # Initial box optimizations take a different input:
        if (k == 1) {
            op_inputs <- best_box_bevals[,par_names]
        } else {
            op_inputs <- optim_outs[[k-1]][,paste0(par_names, "_end")]
        }
        # Last optimizer outputs different object type:
        if (k == n_optims) {
            op_out <- \(op, pars) return(op)
            op_bind_sort <- function(x) {
                vals <- sapply(x, \(y) get_val(y, pkgs[[k]]))
                sort_idx <- order(vals)
                x <- x[sort_idx]
                return(x)
            }
        } else {
            op_out <- function(op, pars) {
                row <- rbind(c(pars, op$par, get_converge_code(op, pkgs[[k]]),
                               get_val(op, pkgs[[k]])))
                colnames(row) <- c(paste0(par_names, "_start"),
                                   paste0(par_names, "_end"),
                                   "conv_code", "val")
                return(row)

            }
            op_bind_sort <- function(x) {
                x <- do.call(rbind, x)
                x <- x[order(x[,"val"]),]
                return(x)
            }
        }
        optim_outs[[k]] <- lapply(1:n, \(i) {
            pars <- unname(op_inputs[i,])
            op <- run_optim(optimizers[[k]], pars, fn, controls[[k]], fn_args)
            return(op_out(op, pars))
        }) |>
            op_bind_sort()
        # Optionally write file if not on the final optimization:
        if (k < n_optims && !is.na(file_names[[k+1]])) {# `+1` bc first file is for box evaluations
            write.table(optim_outs[[k]], file_names[[k+1]],
                        quote = FALSE, row.names = FALSE,
                        sep = ifelse(file_ext(file_names[[k+1]]) == "txt", "\t", ","))
        }
    }

    best_ops <- optim_outs[[n_optims]][1:tail(n_outputs, 1)]
    not_conv <- sapply(best_ops, \(o) !get_converged(o, tail(pkgs, 1)))

    for (i in which(not_conv)) {
        warning("Final optimization ", i, " did not converged.")
    }

    return(best_ops)

}




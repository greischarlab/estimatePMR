

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
#' @param fn Objective function to minimize.
#' @param lower_bounds Lower bounds of boxes for each parameter.
#'     If this vector has names, these names are used for column names in
#'     intermediate output tables (arguments starting with `file_` below).
#'     If names are present, they must match those in `upper_bounds`,
#'     but don't need to be in the same order.
#' @param upper_bounds Upper bounds of boxes for each parameter.
#'     See `lower_bounds` for providing names for this vector.
#' @param fn_args List containing other arguments to use for `fn`.
#'     Defaults to `list()`.
#' @param box_control List containing arguments to use for the `control`
#'     argument for the optimization that occurs for each box.
#'     Defaults to `list(maxit = 100, reltol = 1e-4)`
#'     (these names are adjusted if using an optimizer other than `optim`).
#' @param fine_control List containing arguments to use for the `control`
#'     argument for the fine optimizations.
#'     Defaults to `list(maxit = 500, reltol = 1e-6)`
#'     (these names are adjusted if using an optimizer other than `optim`).
#' @param polished_control List containing arguments to use for the `control`
#'     argument for the polished optimizations.
#'     Defaults to `list(maxit = 1000, reltol = 1e-8)`
#'     (these names are adjusted if using an optimizer other than `optim`).
#' @param box_optim Function to use for the optimization for each box.
#'     Must be `stats::optim` or a function from packages `minqa` or `nloptr`.
#'     Defaults to `stats::optim`.
#' @param fine_optim Function to use for the fine optimizations.
#'     Must be `stats::optim` or a function from packages `minqa` or `nloptr`.
#'     Defaults to `stats::optim`.
#' @param polished_optim Function to use for the polished optimizations.
#'     Must be `stats::optim` or a function from packages `minqa` or `nloptr`.
#'     Defaults to `stats::optim`.
#' @param n_bevals Number of evaluations of `fn` per box. Defaults to `100L`.
#' @param n_boxes Number of boxes. Defaults to `1000L`.
#' @param n_fine Number of fine optimizations.
#'     Must be `< n_boxes` and `> n_polished`.
#'     Defaults to `100L`.
#' @param n_polished Number of polished optimizations.
#'     Must be `< n_fine`.
#'     Defaults to `20L`.
#' @param n_outputs Number of top output object(s) to return.
#'     Must be `< n_polished`. Defaults to `3L`.
#' @param file_bevals Single string specifying the filename where to save the
#'     output from the box evaluations step.
#'     The output is a table with the box number (column `box`),
#'     parameter values (named based on bounds arguments or `par1`, `par2`, ...),
#'     and output from objective function at those values (`val`).
#'     If the file ends with `.txt`, the output will be tab-delimited,
#'     and if the file ends with `.csv`, the output will be comma-delimited.
#'     Other extensions are not allowed.
#'     An error will trigger if attempting to overwrite an existing file
#'     unless the `overwrite` argument is set to `TRUE`.
#'     If `NULL` (the default), no output is written.
#' @param file_boxes Single string specifying the filename where to save the
#'     output from the box evaluations step.
#'     The output is a table with the starting parameter values (with `_start` suffix),
#'     ending parameter values (with `_end` suffix),
#'     convergence code for optimization (`conv_code`), and
#'     output from objective function ending parameter values (`val`).
#'     See description of argument `file_bevals` above for file extensions and
#'     overwriting.
#'     If `NULL` (the default), no output is written.
#' @param file_fine Single string specifying the filename where to save the
#'     output from the box evaluations step.
#'     The output is a table with the same columns as for `file_boxes`
#'     (see above).
#'     See description of argument `file_bevals` above for file extensions and
#'     overwriting.
#'     If `NULL` (the default), no output is written.
#' @param overwrite A single logical for whether to allow overwriting files
#'     for intermediate output (arguments to this function starting
#'     with `file_`). Defaults to `FALSE`.
#' @param na_stop Single logical for whether to return the matrix of initial
#'     evaluations in each box if it contains `NA`s.
#'     If `FALSE`, the optimizer ignores these values and continues on
#'     (unless they're all `NA`).
#'     Defaults to `FALSE`.
#'
#' @importFrom stats optim
#' @importFrom tools file_ext
#'
#' @returns A list containing `n_outputs` object(s) of the class returned
#'     by `polished_optim`.
#'
#'
#' @example
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
                            box_control = list(maxit = 100, reltol = 1e-4),
                            fine_control = list(maxit = 500, reltol = 1e-6),
                            polished_control = list(maxit = 1000, reltol = 1e-8),
                            box_optim = optim,
                            fine_optim = optim,
                            polished_optim = optim,
                            n_bevals = 100L,
                            n_boxes = 1000L,
                            n_fine = 100L,
                            n_polished = 20L,
                            n_outputs = 3L,
                            file_bevals = NULL,
                            file_boxes = NULL,
                            file_fine = NULL,
                            overwrite = FALSE,
                            na_stop = FALSE) {

    # Type and length checks:
    stopifnot(is.function(fn))
    stopifnot(is.numeric(lower_bounds) && is.numeric(upper_bounds))
    stopifnot(length(lower_bounds) == length(upper_bounds))
    stopifnot(is.list(fn_args))
    stopifnot(is.list(box_control))
    stopifnot(is.list(fine_control))
    stopifnot(is.list(polished_control))
    stopifnot(is.function(box_optim))
    stopifnot(is.function(fine_optim))
    stopifnot(is.function(polished_optim))
    stopifnot(length(n_bevals) == 1L && as.integer(n_bevals) == n_bevals)
    stopifnot(length(n_boxes) == 1L && as.integer(n_boxes) == n_boxes)
    stopifnot(length(n_fine) == 1L && as.integer(n_fine) == n_fine)
    stopifnot(length(n_polished) == 1L && as.integer(n_polished) == n_polished)
    stopifnot(length(n_outputs) == 1L && as.integer(n_outputs) == n_outputs)
    stopifnot(length(na_stop) == 1L && inherits(na_stop, "logical"))
    stopifnot(length(overwrite) == 1L && inherits(overwrite, "logical"))

    # Value checks:
    stopifnot(all(lower_bounds < upper_bounds))
    stopifnot(n_bevals > 0L)
    stopifnot(n_boxes > 0L)
    stopifnot(n_fine > 0L && n_fine < n_boxes)
    stopifnot(n_polished > 0L && n_polished < n_fine)
    stopifnot(n_outputs > 0L && n_outputs < n_polished)

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

    # Check validity of intermediate file names:
    interm_files <- list("file_bevals" = file_bevals, "file_boxes" = file_boxes,
                         "file_fine" = file_fine)
    for (n in names(interm_files)) {
        if (!is.null(interm_files[[n]])) {
            if (length(interm_files[[n]]) != 1 || !is.character(interm_files[[n]])) {
                stop(sprintf("Argument `%s` should be NULL or a single string", n))
            }
            if (! file_ext(interm_files[[n]]) %in% c("txt", "csv")) {
                stop(sprintf("Extension for argument `%s` should be 'txt' or 'csv'", n))
            }
            if (!dir.exists(dirname(interm_files[[n]]))) {
                stop(sprintf("Directory '%s' does not exist for argument `%s`",
                             dirname(interm_files[[n]]), n))
            }
            if (!overwrite && file.exists(interm_files[[n]])) {
                stop(paste0("File for argument `", n, "` ('", interm_files[[n]],
                            "') already exists and `overwrite` is `FALSE`"))
            }
        }
    }


    pkgs <- list(box = packageName(environment(box_optim)),
                       fine = packageName(environment(fine_optim)),
                       polished = packageName(environment(polished_optim)))

    if (!all(pkgs %in% c("stats", "minqa", "nloptr"))) {
        stop(paste("this function only programmed for stats::optim and",
                   "optimizers from minqa and nloptr packages"))
    }

    mids <- (upper_bounds + lower_bounds) / 2
    steps <- (upper_bounds - mids) / n_boxes
    n_pars <- length(mids)

    box_optims <- matrix(0.0, n_boxes, 2*n_pars+2)
    colnames(box_optims) <- c(paste0(par_names, "_start"),
                             paste0(par_names, "_end"), "conv_code", "val")
    evals_i <- matrix(0.0, n_bevals, n_pars+2)
    colnames(evals_i) <- c("box", par_names, "val")

    for (i in 1:n_boxes) {
        for (j in 1:n_bevals) {
            .pars <- runif(n_pars, mids - steps * i, mids + steps * i)
            .val <- do.call(fn, c(list(.pars), fn_args))
            evals_i[j,] <- c(i, .pars, .val)
        }
        if (any(is.na(evals_i[,"val"])) && na_stop) {
            warning(sprintf(paste("\nThere were NAs in the %ith box.",
                                  "The matrix of parameter values and",
                                  "evaluations is being returned."), i))
            return(evals_i)
        }
        if (all(is.na(evals_i[,"val"]))) {
            warning(sprintf(paste("\nThe %ith box was all NAs!",
                                  "The matrix of parameter values and",
                                  "evaluations is being returned."), i))
            return(evals_i)
        }
        if (!is.null(file_bevals)) {
            write.table(evals_i, file_bevals, append = i > 1, quote = FALSE,
                        sep = ifelse(file_ext(file_bevals) == "txt", "\t", ","),
                        row.names = FALSE, col.names = i == 1)
        }

        best_idx <- which(evals_i[,"val"] == min(evals_i[,"val"], na.rm = TRUE))[[1]]
        best_pars <- unname(evals_i[best_idx, par_names])
        op <- run_optim(box_optim, best_pars, fn, box_control, fn_args)

        box_optims[i,] <- c(best_pars, op$par, get_converge_code(op, pkgs$box),
                            get_val(op, pkgs$box))
    }

    # sort with lowest at the top:
    box_optims <- box_optims[order(box_optims[,"val"]),]
    # Output if requested:
    if (!is.null(file_boxes)) {
        write.table(box_optims, file_boxes, quote = FALSE, row.names = FALSE,
                    sep = ifelse(file_ext(file_boxes) == "txt", "\t", ","))
    }

    fines <- matrix(0.0, n_fine, ncol(box_optims))
    colnames(fines) <- colnames(box_optims)
    for (i in 1:n_fine) {
        pars <- unname(box_optims[i,paste0(par_names, "_end")])
        op <- run_optim(fine_optim, pars, fn, fine_control, fn_args)
        fines[i,] <- c(pars, op$par, get_converge_code(op, pkgs$fine),
                       get_val(op, pkgs$fine))
    }
    # sort with lowest at the top:
    fines <- fines[order(fines[,"val"]),]
    # Output if requested:
    if (!is.null(file_fine)) {
        write.table(fines, file_fine, quote = FALSE, row.names = FALSE,
                    sep = ifelse(file_ext(file_fine) == "txt", "\t", ","))
    }

    # keep all optimizer output in this case:
    polished_op <- lapply(1:n_polished, \(i) {
        pars <- unname(fines[i,paste0(par_names, "_end")])
        op <- run_optim(polished_optim, pars, fn, polished_control, fn_args)
        return(op)
    })
    polished_vals <- sapply(polished_op, \(x) get_val(x, pkgs$polished))
    sort_idx <- order(polished_vals)
    polished_op <- polished_op[sort_idx]
    polished_vals <- polished_vals[sort_idx]

    best_ops <- polished_op[1:n_outputs]
    not_conv <- sapply(best_ops, \(o) !get_converged(o, pkgs$polished))

    for (i in which(not_conv)) {
        warning("Final optimization ", i, " did not converged.")
    }

    return(best_ops)

}




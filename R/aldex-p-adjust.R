##' Adjust p-values jointly across fitted coefficients and contrasts
##'
##' Recomputes multiple-testing adjustments from an existing ALDEx3 fit,
##' without refitting or drawing new samples. Each group pools all features
##' of its named coefficients or stored contrasts into one adjustment family.
##'
##' @param object An object returned by \code{aldex}, retaining the per-draw
##'   \code{estimate}, \code{std.error}, and \code{df} arrays for coefficients
##'   and any requested contrasts. These arrays are returned by default,
##'   including during streaming. Composition and scale draws are not needed.
##' @param groups A character vector of coefficient or contrast names for one
##'   family, or a list of character vectors for several disjoint families.
##'   Unlisted rows each retain a separate family across all features.
##'   \code{NULL} restores separate adjustment for every row. Use exact row
##'   names. \strong{Prefixes are not normally needed}: a name that occurs
##'   only once among all coefficient and contrast rows can be used directly,
##'   for example \code{groups = c("age_effect", "group_effect")}. If both
##'   components have a row named \code{"age"}, use \code{"coefficient:age"}
##'   or \code{"contrast:age"} to select the intended row. Either prefix may
##'   also be used when the name is unique. Qualified names take precedence
##'   over unqualified names. A row cannot occur twice.
##'   Intercepts are included only when explicitly named in a group.
##' @param method A method accepted by \code{\link[stats]{p.adjust}}. By
##'   default, uses \code{object$p.adjust.method}, or \code{"BH"} if absent.
##' @param onesided \code{FALSE} for two-sided inference, or \code{"lower"}
##'   or \code{"upper"} for a directional test. By default, uses
##'   \code{object$onesided}, or \code{FALSE} if absent. Fits made before
##'   version 1.6.1 did not record either setting: explicitly supply
##'   \code{method} and \code{onesided} if those fits used nondefault settings.
##' @details Directional probabilities are reconstructed as
##'   \code{pt(estimate / std.error, df)} using both tails directly. For
##'   two-sided inference, each tail is doubled and capped at one. Adjustment
##'   pools the feature-by-row hypotheses in a group separately within each
##'   draw and direction, before averaging over draws and selecting the
##'   smaller directional aggregate. This preserves ALDEx3's aggregation
##'   procedure; applying \code{p.adjust} to the final \code{p.val} matrix
##'   would generally give different results.
##'
##'   Groups define the multiplicity families; separate groups do not provide
##'   a joint adjustment for their combined discovery list. Choose families
##'   from the scientific questions rather than from observed significance.
##'   The assumptions of the selected adjustment method still apply, and
##'   grouping does not establish a new FDR guarantee for Monte Carlo
##'   aggregated probabilities. Only contrasts already present in the fit
##'   can be included; this function does not construct new contrasts.
##' @section Choosing groups:
##'   Find the available coefficient names with \code{rownames(fit$estimate)}
##'   and stored contrast names with \code{rownames(fit$contrasts$matrix)}.
##'   Use these exact names, including factor level suffixes where present.
##'
##'   \code{groups = c("age", "group")} puts all features for both named
##'   rows in one family. With D features, this is \code{2 * D} hypotheses per draw
##'   and direction. A list creates separate families: for example,
##'   \code{groups = list(c("age", "group"), c("early", "late"))}
##'   creates two families of \code{2 * D} hypotheses each. Rows not named in either
##'   family are still adjusted individually across their D features.
##'   These examples assume the names are unique across both components.
##'   A family may contain coefficients, contrasts, or both; prefixes only
##'   identify which row is meant and do not determine family membership.
##' @return A modified \code{aldex} object, with \code{p.val} and
##'   \code{p.val.adj} matrices for coefficients and stored contrasts.
##'   Raw probabilities change only if \code{onesided} changes. Estimates,
##'   standard errors, degrees of freedom, and draws are preserved.
##'   \code{p.adjust.method}, \code{onesided}, and \code{p.adjust.groups}
##'   record the settings; the latter is a list of fully qualified row names,
##'   including unlisted singleton families. Assign the result to keep it.
##' @seealso \code{\link{aldex}}, \code{\link{summary.aldex}}
##' @examples
##' set.seed(16)
##' Y <- matrix(rpois(4 * 12, 50), 4, 12)
##' metadata <- data.frame(age = seq_len(12), group = rep(c(0, 1), each = 6))
##' C <- rbind(age_effect = c(0, 1, 0), group_effect = c(0, 0, 1))
##' colnames(C) <- colnames(model.matrix(~ age + group, metadata))
##' fit <- aldex(Y, ~ age + group, metadata, nsample = 8,
##'              scale = clr.sm, gamma = 0.2, contrasts = C)
##' rownames(fit$estimate)         # "(Intercept)", "age", "group"
##' rownames(fit$contrasts$matrix) # "age_effect", "group_effect"
##' # Unique contrast names need no prefix. Pool all four features from
##' # each of these two contrasts into one family of eight hypotheses.
##' joint <- aldex.p.adjust(fit, groups = c("age_effect", "group_effect"))
##' summary(joint, component = "contrasts")
##' # A family can also mix coefficients and stored contrasts.
##' mixed <- aldex.p.adjust(fit, groups = c("age", "group_effect"))
##' # Separate families for coefficients and contrasts.
##' separate <- aldex.p.adjust(fit, groups = list(
##'   c("age", "group"), c("age_effect", "group_effect")))
##' # If a contrast shares a coefficient name, identify which one you mean.
##' C.shared <- C
##' rownames(C.shared)[1] <- "age"
##' fit.shared <- aldex(Y, ~ age + group, metadata, nsample = 8,
##'                     scale = clr.sm, gamma = 0.2, contrasts = C.shared)
##' # "age" alone is now ambiguous. Select the contrast named "age":
##' joint.shared <- aldex.p.adjust(fit.shared,
##'   groups = c("contrast:age", "group_effect"))
##' # Select the coefficient instead, grouped with the other contrast:
##' mixed.shared <- aldex.p.adjust(fit.shared,
##'   groups = c("coefficient:age", "group_effect"))
##' @details Not yet available for test="ssrv.bootstrap".
##' @export
aldex.p.adjust <- function(object, groups=NULL, method=NULL, onesided=NULL) {
  .aldex_guard_ssrv(object, "aldex.p.adjust")
  if (!inherits(object, "aldex") || !is.list(object)) {
    stop("object must be an aldex fit")
  }
  if (is.null(method)) method <- object$p.adjust.method
  if (is.null(method)) method <- "BH"
  if (!is.character(method) || length(method) != 1L || is.na(method)) {
    stop("method must be a single p.adjust method")
  }
  method <- match.arg(method, stats::p.adjust.methods)
  if (is.null(onesided)) onesided <- object$onesided
  if (is.null(onesided)) onesided <- FALSE
  if (!identical(onesided, FALSE) && !identical(onesided, "lower") &&
      !identical(onesided, "upper")) {
    stop("onesided must be FALSE, 'lower', or 'upper'")
  }

  validate <- function(values, component) {
    fields <- c("estimate", "std.error", "df")
    for (field in fields) {
      x <- values[[field]]
      if (is.null(x)) {
        stop(component, " must retain per-draw ", field,
             "; include 'estimate', 'std.error', and 'df' in return.pars when fitting")
      }
      if (!is.numeric(x) || length(dim(x)) != 3L || any(dim(x) < 1L)) {
        stop(component, "$", field, " must be a numeric row x feature x draw array")
      }
      if (!identical(dim(x), dim(values$estimate)) ||
          !identical(dimnames(x), dimnames(values$estimate))) {
        stop(component, " estimate, std.error, and df arrays must have matching dimensions and names")
      }
    }
    rn <- rownames(values$estimate)
    if (is.null(rn) || anyNA(rn) || any(!nzchar(rn)) || anyDuplicated(rn)) {
      stop(component, " must have unique nonempty row names")
    }
    if (any(values$std.error < 0, na.rm=TRUE) ||
        any(values$df <= 0, na.rm=TRUE)) {
      stop(component, " must have nonnegative std.error and positive df")
    }
    rn
  }
  coefficient.names <- validate(object, "coefficients")
  contrast.names <- character()
  if (!is.null(object$contrasts)) {
    contrast.names <- validate(object$contrasts, "contrasts")
    if (!identical(dim(object$estimate)[2:3], dim(object$contrasts$estimate)[2:3]) ||
        !identical(dimnames(object$estimate)[2:3],
                   dimnames(object$contrasts$estimate)[2:3])) {
      stop("coefficients and contrasts must have matching feature and draw dimensions and names")
    }
  }
  resolved <- .aldex_resolve_p_adjust_groups(groups, coefficient.names, contrast.names)
  labels <- c(paste0("coefficient:", coefficient.names),
              if (length(contrast.names)) paste0("contrast:", contrast.names))
  if (is.null(resolved)) {
    resolved <- list(indices=as.list(seq_along(labels)), labels=as.list(labels))
  }

  tails <- function(values) {
    statistic <- values$estimate / values$std.error
    list(lower=array(stats::pt(statistic, values$df, lower.tail=TRUE),
                     dim(statistic), dimnames(statistic)),
         upper=array(stats::pt(statistic, values$df, lower.tail=FALSE),
                     dim(statistic), dimnames(statistic)))
  }
  probabilities <- tails(object)
  if (length(contrast.names)) {
    contrast.probabilities <- tails(object$contrasts)
    probabilities <- list(
      lower=abind::abind(probabilities$lower, contrast.probabilities$lower, along=1),
      upper=abind::abind(probabilities$upper, contrast.probabilities$upper, along=1))
  }
  adjusted <- aldex.pvals(probabilities$lower, probabilities$upper, method,
                          onesided, groups=resolved$indices)
  dimnames(adjusted$p.res) <- dimnames(adjusted$p.adj.res) <-
    list(c(coefficient.names, contrast.names), colnames(object$estimate))
  coefficient.rows <- seq_along(coefficient.names)
  object$p.val <- adjusted$p.res[coefficient.rows, , drop=FALSE]
  object$p.val.adj <- adjusted$p.adj.res[coefficient.rows, , drop=FALSE]
  if (length(contrast.names)) {
    contrast.rows <- length(coefficient.names) + seq_along(contrast.names)
    object$contrasts$p.val <- adjusted$p.res[contrast.rows, , drop=FALSE]
    object$contrasts$p.val.adj <- adjusted$p.adj.res[contrast.rows, , drop=FALSE]
  }
  object$p.adjust.method <- method
  object$onesided <- onesided
  object$p.adjust.groups <- resolved$labels
  object
}

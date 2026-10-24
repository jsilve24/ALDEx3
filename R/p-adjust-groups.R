.aldex_resolve_p_adjust_groups <- function(groups, coefficient.names,
                                           contrast.names=character()) {
  if (is.null(groups)) return(NULL)

  coefficient.names <- .aldex_validate_p_adjust_group_names(
    coefficient.names, "coefficient.names"
  )
  contrast.names <- .aldex_validate_p_adjust_group_names(
    contrast.names, "contrast.names"
  )
  canonical <- c(
    if (length(coefficient.names)) paste0("coefficient:", coefficient.names),
    if (length(contrast.names)) paste0("contrast:", contrast.names)
  )
  unqualified <- c(coefficient.names, contrast.names)

  if (is.character(groups) && is.null(dim(groups))) {
    groups <- list(groups)
  } else if (!is.list(groups) || !length(groups)) {
    stop("`groups` must be a character vector or a nonempty list of character vectors")
  }

  groups <- lapply(groups, .aldex_validate_p_adjust_group_tokens)
  used <- integer()
  indices <- lapply(groups, function(group) {
    resolved <- unname(vapply(group, function(token) {
      canonical.match <- which(canonical == token)
      if (length(canonical.match)) return(canonical.match)

      if (identical(token, "coefficient:") || identical(token, "contrast:")) {
        stop("malformed p-adjustment group token: `", token, "`")
      }
      unqualified.match <- which(unqualified == token)
      if (!length(unqualified.match)) {
        stop("unknown p-adjustment group token: `", token, "`")
      }
      if (length(unqualified.match) > 1L) {
        stop("ambiguous p-adjustment group token: `", token,
             "`; use `coefficient:", token, "` or `contrast:", token, "`")
      }
      unqualified.match
    }, integer(1)))
    if (anyDuplicated(resolved) || any(resolved %in% used)) {
      stop("p-adjustment groups cannot repeat a coefficient or contrast")
    }
    used <<- c(used, resolved)
    resolved
  })

  remaining <- setdiff(seq_along(canonical), used)
  indices <- c(indices, lapply(remaining, function(index) as.integer(index)))
  list(indices=indices,
       labels=lapply(indices, function(index) canonical[index]))
}

.aldex_validate_p_adjust_group_names <- function(names, argument) {
  if (!is.character(names) || !is.null(dim(names)) || anyNA(names) ||
      any(!nzchar(names)) || anyDuplicated(names)) {
    stop("`", argument,
         "` must be a character vector of unique, nonempty names without missing values")
  }
  names
}

.aldex_validate_p_adjust_group_tokens <- function(group) {
  if (!is.character(group) || !is.null(dim(group)) || !length(group) ||
      anyNA(group) || any(!nzchar(group))) {
    stop("each p-adjustment group must be a nonempty character vector without missing or empty tokens")
  }
  group
}

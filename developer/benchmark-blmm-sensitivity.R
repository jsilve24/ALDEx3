#!/usr/bin/env Rscript
# Compare two installed ALDEx3 trees without loading their same-named DLLs together.
# Rscript developer/benchmark-blmm-sensitivity.R OLD_LIBRARY NEW_LIBRARY INPUTS_RDS NEW_OUTPUT [REPS=3]
# Build BOTH trees with R CMD INSTALL --preclean and matching release flags first.
# Reusing pkgload/debug objects can silently compare -O0 against optimized code.
args <- commandArgs(TRUE)
self <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1L])
fields <- c("estimate", "std.error", "df", "p.lower", "p.upper", "random.eff")
clock <- function() proc.time()[[3L]]
write_table <- function(rows, file) {
  write.csv(if (length(rows)) do.call(rbind, rows) else data.frame(), file, row.names = FALSE)
}
compiler_provenance <- function(lib) {
  files <- list.files(file.path(lib, "ALDEx3", "libs"), recursive = TRUE, full.names = TRUE)
  dlls <- files[basename(files) == paste0("ALDEx3", .Platform$dynlib.ext)]
  producers <- character()
  readelf <- Sys.which("readelf")
  if (nzchar(readelf)) for (dll in dlls) {
    info <- suppressWarnings(system2(readelf, c("--debug-dump=info", "--dwarf-depth=1", shQuote(dll)),
                                     stdout = TRUE, stderr = FALSE))
    if (is.null(attr(info, "status"))) producers <- c(producers,
      sub("^.*DW_AT_producer\\s*:\\s*(\\([^)]*\\):\\s*)?", "",
          grep("DW_AT_producer", info, value = TRUE), perl = TRUE))
  }
  producers <- unique(producers)
  flags <- regmatches(producers, gregexpr("(?<!\\S)-O(?:[0-9]+|g|s|z|fast)?(?!\\S)",
                                         producers, perl = TRUE))
  effective <- vapply(flags, function(x) if (length(x)) tail(x, 1L) else NA_character_, character(1))
  effective[!is.na(effective) & effective == "-O"] <- "-O1"
  verified <- length(effective) > 0L && !anyNA(effective)
  list(dlls = dlls, md5 = tools::md5sum(dlls), producers = producers,
       effective_optimization = sort(unique(effective)), verified = verified,
       status = if (verified) "DLL DW_AT_producer optimization flags available" else
         "DLL compiler flags unavailable: requires manual compiler verification")
}
compiler_text <- function(info) c(info$status,
  paste("Effective optimization:", paste(info$effective_optimization, collapse = ", ")),
  paste(names(info$md5), unname(info$md5)), info$producers)
validate_inputs <- function(fixtures) {
  if (!is.list(fixtures) || !length(fixtures)) stop("Inputs must be a nonempty fixture list")
  if (is.null(names(fixtures))) names(fixtures) <- paste0("fixture_", seq_along(fixtures))
  if (anyNA(names(fixtures)) || any(!nzchar(names(fixtures))) || anyDuplicated(names(fixtures)))
    stop("Fixture names must be nonempty and unique")
  for (z in fixtures) {
    if (!is.list(z) || !is.array(z$logW) || length(dim(z$logW)) != 3L ||
        any(dim(z$logW) < 1L) || !all(is.finite(z$logW)) ||
        !is.data.frame(z$data) || nrow(z$data) != dim(z$logW)[1L] ||
        !inherits(z$formula, "formula")) stop("Invalid fixture")
  }
  fixtures
}
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1",
           VECLIB_MAXIMUM_THREADS = "1", BLIS_NUM_THREADS = "1")

if (length(args) && identical(args[1L], "--worker")) {
  if (length(args) != 6L) stop("Invalid worker arguments")
  lib <- args[2L]; input <- args[3L]; output <- args[4L]
  engine <- args[5L]; repetition <- as.integer(args[6L])
  compiler <- readRDS(file.path(dirname(input), "compiler-provenance.rds"))[[engine]]
  if (!identical(compiler$md5, tools::md5sum(compiler$dlls))) stop("DLL changed after compiler preflight")
  writeLines(compiler_text(compiler), sub("\\.rds$", "-compiler.txt", output))
  t0 <- clock()
  library("ALDEx3", lib.loc = lib, character.only = TRUE)
  ns <- asNamespace("ALDEx3")
  options(ALDEx3.blmm.adaptive.diagnostics = TRUE)
  loading_seconds <- clock() - t0
  package_path <- normalizePath(find.package("ALDEx3"))
  package_files <- list.files(package_path, recursive = TRUE, full.names = TRUE,
                               all.files = TRUE, no.. = TRUE)
  package_files <- package_files[!dir.exists(package_files)]
  checksums <- tools::md5sum(package_files)
  t0 <- clock(); fixtures <- validate_inputs(readRDS(input)); input_seconds <- clock() - t0
  warnings <- errors <- timings <- phases <- answers <- snapshots <- list()
  if (repetition == 1L) for (name in names(fixtures)) {
    z <- fixtures[[name]]
    snapshots[[name]] <- tryCatch(withCallingHandlers({
      parsed <- ns$blmm_parse_formula(z$formula, z$data, backend = "auto")
      theta <- ifelse(is.finite(parsed$lower), .7, -.2)
      Y <- matrix(z$logW[, 1L, ], nrow = dim(z$logW)[1L])
      objective <- ns$blmm_make_adfun(parsed$X, Y, parsed$basis,
        integer(length(theta)), theta)
      lapply(list(interior = theta, zero = numeric(length(theta))), function(location)
        list(location = location, objective = objective$fn(location),
             gradient = objective$gr(location), Hessian = objective$he(location)))
    }, warning = function(w) {
      warnings[[length(warnings) + 1L]] <<- data.frame(fixture = name,
        rep = repetition, engine, phase = "snapshot", message = conditionMessage(w))
      invokeRestart("muffleWarning")
    }), error = function(e) {
      errors[[length(errors) + 1L]] <<- data.frame(fixture = name,
        rep = repetition, engine, phase = "snapshot", message = conditionMessage(e))
      e
    })
  }
  engine_run <- function(z) {
    fit <- ns$blmm(z$logW, z$formula, z$data, n.cores = 1L,
                   inference = "adaptive", contrasts = z$contrasts)
    # Explicit extraction is part of the timed operation; preserve the raw fit too.
    arrays <- fit[intersect(c(fields, paste0("contrast.", fields)), names(fit))]
    tfit <- clock()
    prefixes <- c(coefficients = "", if (!is.null(z$contrasts)) c(contrasts = "contrast."))
    final <- lapply(prefixes, function(pre)
      ns$aldex.pvals(arrays[[paste0(pre, "p.lower")]],
                     arrays[[paste0(pre, "p.upper")]], "BH", FALSE))
    list(raw = fit, fit = arrays, final = final,
         diagnostics = attr(fit, "blmm.adaptive.diagnostics"), aggregation_started = tfit)
  }
  run <- function(z, fixture, phase) {
    gc(); started <- clock()
    value <- tryCatch(withCallingHandlers(engine_run(z), warning = function(w) {
      warnings[[length(warnings) + 1L]] <<- data.frame(fixture, rep = repetition,
        engine, phase, message = conditionMessage(w))
      invokeRestart("muffleWarning")
    }), error = function(e) e)
    ended <- clock(); failed <- inherits(value, "error")
    if (failed) errors[[length(errors) + 1L]] <<- data.frame(fixture, rep = repetition,
      engine, phase, message = conditionMessage(value))
    row <- data.frame(fixture, rep = repetition, engine, phase,
                      seconds = ended - started, failed)
    if (phase == "fit") timings[[length(timings) + 1L]] <<- row else
      phases[[length(phases) + 1L]] <<- row
    if (!failed) {
      row$engine_seconds <- value$aggregation_started - started
      row$aggregation_seconds <- ended - value$aggregation_started
      for (part in c("engine", "aggregation")) {
        detail <- row[, c("fixture", "rep", "engine", "phase", "seconds", "failed")]
        detail$phase <- paste(phase, part, sep = "_")
        detail$seconds <- row[[paste0(part, "_seconds")]]
        phases[[length(phases) + 1L]] <<- detail
      }
      value$timing <- row; value$aggregation_started <- NULL
    }
    value
  }
  # Warm every design before timing any fixture, excluding all warmup work.
  for (name in names(fixtures)) {
    z <- fixtures[[name]]
    z$logW <- z$logW[, 1L, seq_len(min(2L, dim(z$logW)[3L])), drop = FALSE]
    invisible(run(z, name, "warmup"))
  }
  for (name in names(fixtures)) {
    message(engine, " repetition ", repetition, ": ", name)
    answers[[name]] <- run(fixtures[[name]], name, "fit")
  }
  if (!identical(checksums, tools::md5sum(package_files)))
    stop("Installed package changed during worker execution")
  saveRDS(list(answers = answers, snapshots = snapshots, timings = timings, phases = phases,
    warnings = warnings, errors = errors, provenance = list(engine = engine,
      library = normalizePath(lib), package_path = package_path,
      description = packageDescription("ALDEx3"), files = data.frame(
        file = substring(names(checksums), nchar(package_path) + 2L), md5 = unname(checksums)),
      session = capture.output(sessionInfo()), input_md5 = unname(tools::md5sum(input)),
      loading_seconds = loading_seconds, input_seconds = input_seconds,
      compiler = compiler,
      thread_environment = Sys.getenv(c("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS",
                                        "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "BLIS_NUM_THREADS")))), output)
  quit(status = 0L)
}

if (length(args) < 4L || length(args) > 5L)
  stop("Usage: benchmark-blmm-sensitivity.R OLD_LIBRARY NEW_LIBRARY INPUTS_RDS NEW_OUTPUT [REPS=3]")
reps <- if (length(args) == 5L) suppressWarnings(as.integer(args[5L])) else 3L
if (length(reps) != 1L || is.na(reps) || reps < 1L) stop("REPS must be a positive integer")
libraries <- setNames(vapply(args[1:2], normalizePath, character(1), mustWork = TRUE), c("old", "new"))
for (lib in libraries) if (!file.exists(file.path(lib, "ALDEx3", "DESCRIPTION")))
  stop("Each library must contain an installed ALDEx3 package")
input <- normalizePath(args[3L], mustWork = TRUE)
fixtures <- validate_inputs(readRDS(input))
out <- args[4L]
if (file.exists(out)) stop("Output already exists; use a new output directory")
dir.create(out, recursive = TRUE); out <- normalizePath(out)
dir.create(file.path(out, "workers"))
invisible(file.copy(input, file.path(out, "inputs.rds")))
input <- file.path(out, "inputs.rds")
invisible(file.copy(self, file.path(out, "benchmark-blmm-sensitivity.R")))
writeLines(unname(tools::md5sum(input)), file.path(out, "inputs-md5.txt"))
write.csv(do.call(rbind, lapply(names(fixtures), function(name) {
  z <- fixtures[[name]]
  data.frame(fixture = name, N = dim(z$logW)[1L], D = dim(z$logW)[2L],
    S = dim(z$logW)[3L], formula = paste(deparse(z$formula), collapse = " "))
})), file.path(out, "inputs.csv"), row.names = FALSE)
if (file.exists("/proc/cpuinfo")) invisible(file.copy("/proc/cpuinfo", file.path(out, "cpuinfo.txt")))
compiler <- lapply(libraries, compiler_provenance)
saveRDS(compiler, file.path(out, "compiler-provenance.rds"))
for (engine in names(compiler))
  writeLines(compiler_text(compiler[[engine]]), file.path(out, paste0(engine, "-compiler.txt")))
if (all(vapply(compiler, function(x) x$verified, logical(1)))) {
  if (!identical(compiler$old$effective_optimization, compiler$new$effective_optimization))
    stop("DLL optimization flags differ; rebuild both packages with matching release flags and --preclean")
} else message("DLL compiler flags unavailable: requires manual compiler verification; see *-compiler.txt")
if (identical(Sys.getenv("ALDEX3_BENCH_INPUTS_ONLY"), "1")) {
  message("Inputs validated and copied; no fits run: ", out)
  quit(status = 0L)
}
results <- list(); timings <- phases <- warnings <- errors <- list()
for (rep in seq_len(reps)) {
  order <- if (rep %% 2L) c("old", "new") else c("new", "old")
  for (position in seq_along(order)) {
    engine <- order[position]
    stem <- file.path(out, "workers", paste0(engine, "-", rep))
    command <- c("--vanilla", normalizePath(self), "--worker", libraries[[engine]],
                  input, paste0(stem, ".rds"), engine, rep)
    status <- system2(file.path(R.home("bin"), "Rscript"), vapply(command, shQuote, character(1)),
                       stdout = paste0(stem, ".log"), stderr = paste0(stem, ".log"))
    if (status != 0L || !file.exists(paste0(stem, ".rds")))
      stop("Worker failed; inspect ", stem, ".log")
    result <- readRDS(paste0(stem, ".rds"))
    results[[paste(engine, rep, sep = "-")]] <- result
    for (row in result$timings) { row$position <- position; timings[[length(timings) + 1L]] <- row }
    phases <- c(phases, result$phases); warnings <- c(warnings, result$warnings); errors <- c(errors, result$errors)
    write.csv(result$provenance$files, paste0(stem, "-package-md5.csv"), row.names = FALSE)
    writeLines(result$provenance$session, paste0(stem, "-session.txt"))
  }
}
comparisons <- rescue <- list()
compare_arrays <- function(actual, reference, fixture, rep, component, quantity, scale = NULL, decisions = FALSE) {
  if (!identical(dim(actual), dim(reference)) || !identical(dimnames(actual), dimnames(reference)))
    stop("Mismatched result dimensions or labels: ", fixture, " ", quantity)
  ok <- is.finite(actual) & is.finite(reference); differences <- abs(actual - reference)
  statistic <- function(values, fun) if (length(values)) fun(values) else NA_real_
  comparisons[[length(comparisons) + 1L]] <<- data.frame(fixture, rep, component, quantity,
    pairs = length(ok), finite_pairs = sum(ok), old_nonfinite = sum(!is.finite(reference)),
    new_nonfinite = sum(!is.finite(actual)), median_abs = statistic(differences[ok], median),
    max_abs = statistic(differences[ok], max), max_scaled = if (is.null(scale)) NA_real_ else
      statistic(differences[ok] / pmax(abs(scale[ok]), .Machine$double.eps), max),
    decisions_disagree = if (decisions) sum((actual[ok] <= .05) != (reference[ok] <= .05)) else NA_integer_)
}
for (fixture in names(fixtures)) {
  before <- results[["old-1"]]$snapshots[[fixture]]
  after <- results[["new-1"]]$snapshots[[fixture]]
  if (inherits(before, "error") || inherits(after, "error")) next
  for (state in names(before)) for (quantity in c("objective", "gradient", "Hessian"))
    compare_arrays(after[[state]][[quantity]], before[[state]][[quantity]], fixture,
      1L, paste0("observed_REML_", state), quantity)
}
for (rep in seq_len(reps)) for (fixture in names(fixtures)) {
  before <- results[[paste0("old-", rep)]]$answers[[fixture]]
  after <- results[[paste0("new-", rep)]]$answers[[fixture]]
  if (inherits(before, "error") || inherits(after, "error")) next
  for (component in names(before$final)) {
    pre <- if (component == "contrasts") "contrast." else ""
    for (field in setdiff(fields, "random.eff")) {
      key <- paste0(pre, field)
      scale <- if (field == "estimate") before$fit[[paste0(pre, "std.error")]] else
        if (field %in% c("std.error", "df")) before$fit[[key]] else NULL
      compare_arrays(after$fit[[key]], before$fit[[key]], fixture, rep, component, field, scale)
    }
    for (field in c("p.res", "p.adj.res"))
      compare_arrays(after$final[[component]][[field]], before$final[[component]][[field]],
        fixture, rep, component, field, decisions = TRUE)
  }
  compare_arrays(after$fit$random.eff, before$fit$random.eff, fixture, rep,
                   "random_covariance", "random.eff")
  if (length(before$diagnostics) != length(after$diagnostics)) stop("Mismatched diagnostics")
  for (feature in seq_along(before$diagnostics)) {
    a <- before$diagnostics[[feature]]$diagnostics
    b <- after$diagnostics[[feature]]$diagnostics
    if (!identical(a$draw, b$draw)) stop("Mismatched diagnostic draw indices")
    rescue[[length(rescue) + 1L]] <- data.frame(fixture, rep, feature, draw = a$draw,
      old_rescued = a$rescued, new_rescued = b$rescued, old_reason = a$reason,
      new_reason = b$reason, status_changed = if (is.null(a$status) || is.null(b$status))
        rep(NA, nrow(a)) else a$status != b$status,
      rescue_changed = a$rescued != b$rescued, reason_changed = a$reason != b$reason)
  }
}
write_table(timings, file.path(out, "timing.csv")); write_table(phases, file.path(out, "phases.csv"))
write_table(warnings, file.path(out, "warnings.csv")); write_table(errors, file.path(out, "errors.csv"))
write_table(comparisons, file.path(out, "comparison.csv")); write_table(rescue, file.path(out, "rescue-comparison.csv"))
good <- do.call(rbind, timings); good <- good[!good$failed, , drop = FALSE]
if (nrow(good)) {
  summary <- aggregate(seconds ~ fixture + engine, good,
    function(x) c(median = median(x), min = min(x), max = max(x)))
  write.csv(summary, file.path(out, "timing-summary.csv"), row.names = FALSE)
  paired <- merge(good[good$engine == "old", c("fixture", "rep", "seconds")],
                  good[good$engine == "new", c("fixture", "rep", "seconds")],
                  by = c("fixture", "rep"), suffixes = c("_old", "_new"))
  paired$speedup <- paired$seconds_old / paired$seconds_new
  write.csv(paired, file.path(out, "paired-timing.csv"), row.names = FALSE)
}
writeLines(c("Each old/new repetition runs in a separate R subprocess, sequentially in alternating order.",
  "Times include default adaptive BLMM fitting, inference, extraction and directional BH aggregation.",
  "Loading, input reads, per-design warmup, gc, comparisons and serialization are excluded.",
  "Untimed repetition-one snapshots compare observed REML objective/gradient/Hessian at interior and zero covariance.",
  "Threads are limited through environment variables before worker startup; n.cores=1.",
  "Workers retain raw fits, full adaptive diagnostics, warnings/errors, package hashes and session information.",
  "Compiler provenance comes from DLL DW_AT_producer strings when available; optimization mismatches stop the run.",
  "Missing compiler debug metadata requires manual verification of matching release/preclean builds.",
  "Inputs are public or synthetic only if the supplied fixtures establish that provenance.",
  "Other workloads on the shared machine must be paused externally."), file.path(out, "README.txt"))
message("Finished: ", out)

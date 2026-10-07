# MANUSFIL: RESULTAT FÖR H1–H4 I BERLIN-PRESENTATIONEN
# ------------------------------------------------------------
# 1. Kör den uppdaterade 01.clean_all.R när clean_all.rds behöver skapas.
# 2. Kör HELA denna fil med Source. Du behöver inte köra 04. analyses.R först.
# Samma datamapp som i 01 och 04 används. Ändra vid behov före körningen:
# options(buratti.data_dir = "C:/sökväg/till/dina/data")
#
# Output sparas automatiskt i undermappen manus_output i datamappen:
#   output_04_manus.txt: samtliga resultat, varningar och sessionInfo()
#   *.csv: numeriska tabeller med n, skattningar, SE, CI och p-värden
#   manus_results.rds: tabeller, modeller och bootstrapresultat för återbruk
# Samma outputfiler ersätts vid omkörning. clean_all.rds ändras inte.
#
# Avgränsning: H1 korrelationer och faktormodeller; H2 bootstrapjämförelse;
# H3 domänskillnader och accuracy-analyser; H4 kohortskillnader inom domän.
# Deskriptiv statistik, confidence-bortfall och reliabilitet hör till dessa.
# Separata könsanalyser, gamma och analyser för Paper 2 ingår inte.
# Kohortanalysen kallades H5 i den tidigare 04-filen, men är H4 i presentationen.
#
# Måtten läses färdigberäknade från 01. Bias använder som tidigare separata
# tillgängliga medelvärden för confidence och accuracy; matchad bias införs inte.
# Domänspecifik tiouppgifts-accuracy används i H3:s sensitivitet. Total accuracy
# används i H2 och inomdomänregressionerna. Befintliga mått skrivs inte över.
# Sex01 och Cohort01 är numeriska 0/1-kovariater, som i 04. Kohort 0=1948, 1=1953.
#
# Paket som redan används i 04. Installera vid behov i konsolen (en gång):
# install.packages(c("lmerTest", "emmeans", "lavaan"))
#
# Tolkningsstöd: vanliga CI är 95 %, RMSEA-intervall är 90 %.
# Korrelationernas CI är ojusterade; p justeras med Holm inom varje mått.
# H3:s domänkontraster använder Tukey. H4 behåller en kohortkontrast per domän:
# Holm inom en enda kontrast ändrar inte p och justerar inte över tre domäner.
# Regressionskoefficienternas p-värden och CI är ojusterade. Linjära regressioner
# använder t-intervall; logistiska regressioner använder Wald-intervall för b/OR.
# Bootstrappercentiler för H2 är ett CI. Split-half-percentiler visar däremot
# variationen mellan uppdelningar och är INTE ett gemensamt 95-procentigt CI.
# ICC(C,1) mellan femuppgiftshalvor är inte fulltestreliabilitet eller retest.
#
# API-stöd: https://lavaan.ugent.be/tutorial/inspect.html
# https://rvlenth.github.io/emmeans/articles/confidence-intervals.html

.manus_domains <- function() {
  list(
    Antonyms = list(prefix = "atonyms", n_answered = "atonym_n_answered",
      n_conf = "atonym_n_conf", confidence = paste0("m_con_", 1:10),
      correct = paste0("mCat_", c(24,31,35,37,40,41,42,43,44,45), "_correct"),
      all_correct = paste0("mCat_", 1:45, "_correct")),
    `Metal Folding` = list(prefix = "metal", n_answered = "metal_n_answered",
      n_conf = "metal_n_conf", confidence = paste0("p_con_", 1:10),
      correct = paste0("p_Q", c(10,12,16,19,20,21,24,26,27,28), "_correct"),
      all_correct = paste0("p_Q", 1:40, "_correct")),
    `Number Series` = list(prefix = "number", n_answered = "number_n_answered",
      n_conf = "number_n_conf", confidence = paste0("tcat_t", c(2,10,11,13,14,15,17,18,19,22)),
      correct = sprintf("Tal%02d_correct", c(2,10,11,13,14,15,17,18,19,22)),
      all_correct = sprintf("Tal%02d_correct", 1:40))
  )
}

.manus_describe <- function(data, domains) {
  measures <- c("acc", "conf", "bias", "abs_bias", "slope", "total_acc")
  groups <- list(All = seq_len(nrow(data)), `Cohort 1948` = which(data$Cohort01 == 0),
                 `Cohort 1953` = which(data$Cohort01 == 1))
  rows <- list()
  for (group in names(groups)) for (domain in names(domains)) for (measure in measures) {
    x <- data[[paste0(domains[[domain]]$prefix, "_", measure)]][groups[[group]]]
    x <- x[is.finite(x)]; n <- length(x)
    m <- if (n) mean(x) else NA_real_
    s <- if (n > 1L) sd(x) else NA_real_
    se <- s / sqrt(n)
    half <- if (n > 1L) qt(.975, n - 1L) * se else NA_real_
    rows[[length(rows) + 1L]] <- data.frame(Group = group, Domain = domain,
      Measure = measure, n_group = length(groups[[group]]), n = n,
      n_missing = length(groups[[group]]) - n, M = m, SD = s, SE = se,
      lower95 = m - half, upper95 = m + half,
      min = if (n) min(x) else NA_real_, max = if (n) max(x) else NA_real_)
  }
  do.call(rbind, rows)
}

.manus_correlations <- function(data, domains) {
  measures <- c("total_acc", "conf", "bias", "abs_bias", "slope")
  pairs <- combn(seq_along(domains), 2L)
  rows <- list()
  for (measure in measures) {
    current <- list()
    for (i in seq_len(ncol(pairs))) {
      a <- pairs[1L, i]; b <- pairs[2L, i]
      x <- data[[paste0(domains[[a]]$prefix, "_", measure)]]
      y <- data[[paste0(domains[[b]]$prefix, "_", measure)]]
      ok <- is.finite(x) & is.finite(y); n <- sum(ok)
      r <- p <- lower <- upper <- NA_real_; status <- "OK"
      if (n < 3L || sd(x[ok]) == 0 || sd(y[ok]) == 0) {
        status <- "Not estimable: insufficient pairs or no variance"
      } else {
        fit <- cor.test(x[ok], y[ok], method = "pearson", conf.level = .95)
        r <- unname(fit$estimate); p <- fit$p.value
        if (length(fit$conf.int) == 2L) { lower <- fit$conf.int[1L]; upper <- fit$conf.int[2L] }
      }
      current[[i]] <- data.frame(Measure = measure, Domain_A = names(domains)[a],
        Domain_B = names(domains)[b], n = n, r = r, lower95 = lower,
        upper95 = upper, p_raw = p, Status = status)
    }
    current <- do.call(rbind, current)
    current$p_Holm <- p.adjust(current$p_raw, method = "holm")
    rows[[measure]] <- current
  }
  do.call(rbind, rows)
}

.manus_bootstrap <- function(data, domains, B = 5000L, seed = 12345L) {
  conf_vars <- vapply(domains, function(x) paste0(x$prefix, "_conf"), character(1))
  acc_vars <- vapply(domains, function(x) paste0(x$prefix, "_total_acc"), character(1))
  mean_z_diff <- function(x) {
    conf <- cor(x[conf_vars], use = "pairwise.complete.obs")
    acc <- cor(x[acc_vars], use = "pairwise.complete.obs")
    mean(atanh(conf[lower.tri(conf)])) - mean(atanh(acc[lower.tri(acc)]))
  }
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
  on.exit(if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv) else
    if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
      rm(".Random.seed", envir = .GlobalEnv), add = TRUE)
  set.seed(seed)
  observed <- mean_z_diff(data)
  draws <- replicate(B, mean_z_diff(data[sample(seq_len(nrow(data)), replace = TRUE), , drop = FALSE]))
  ok <- is.finite(draws)
  ci <- if (any(ok)) unname(quantile(draws[ok], c(.025, .975))) else c(NA_real_, NA_real_)
  if (!is.finite(observed) || !all(ok)) warning("H2: some bootstrap estimates are not finite; inspect B_estimable and the data.")
  list(summary = data.frame(n_rows_resampled = nrow(data), B_requested = B,
    B_estimable = sum(ok), seed = seed, mean_Fisher_z_conf_minus_acc = observed,
    lower95 = ci[1L], upper95 = ci[2L]), draws = draws, RNG_kind = RNGkind())
}

.manus_long <- function(data, domains) {
  result <- do.call(rbind, lapply(names(domains), function(domain) {
    prefix <- domains[[domain]]$prefix
    data.frame(id = data$id, Sex01 = data$Sex01, Cohort01 = data$Cohort01,
      Domain = domain, Accuracy = data[[paste0(prefix, "_acc")]],
      Bias = data[[paste0(prefix, "_bias")]],
      Absolute_bias = data[[paste0(prefix, "_abs_bias")]],
      Slope = data[[paste0(prefix, "_slope")]])
  }))
  result$Domain <- factor(result$Domain, levels = names(domains))
  rownames(result) <- NULL
  result
}

# Dessa counts hämtas från exakt de modellrader som användes, inte från df totalt.
.manus_model_n <- function(model, label, input_n) {
  mf <- model.frame(model)
  mixed <- inherits(model, "merMod")
  used_ids <- if (mixed) mf$id else attr(model, "manus_ids")
  data.frame(Model = label, n_input_observations = input_n,
    n_used_observations = nobs(model), n_excluded_observations = input_n - nobs(model),
    n_unique_participants = length(unique(used_ids)))
}

.manus_regression <- function(data, formula, label, logistic = FALSE) {
  vars <- all.vars(formula)
  used <- data[complete.cases(data[vars]), , drop = FALSE]
  fit <- if (logistic) glm(formula, data = used, family = binomial(), na.action = na.omit) else
    lm(formula, data = used, na.action = na.omit)
  attr(fit, "manus_ids") <- used$id
  s <- summary(fit); cf <- s$coefficients
  critical <- if (logistic) qnorm(.975) else qt(.975, df.residual(fit))
  result <- data.frame(Model = label, Term = rownames(cf), b = cf[, 1L],
    SE = cf[, 2L], df = if (logistic) Inf else df.residual(fit),
    statistic = cf[, 3L], test = if (logistic) "z" else "t", p = cf[, 4L],
    lower95 = cf[, 1L] - critical * cf[, 2L],
    upper95 = cf[, 1L] + critical * cf[, 2L], row.names = NULL)
  if (logistic) {
    result$OR <- exp(result$b); result$OR_lower95 <- exp(result$lower95)
    result$OR_upper95 <- exp(result$upper95)
  }
  list(model = fit, coefficients = result, n = .manus_model_n(fit, label, nrow(data)),
    fit = data.frame(Model = label, n = nobs(fit), df_residual = df.residual(fit),
      R2 = if (logistic) NA_real_ else s$r.squared,
      adjusted_R2 = if (logistic) NA_real_ else s$adj.r.squared))
}

.manus_reliability <- function(data, domains) {
  # Samma råa alfa och ICC(C,1) som i 04, enbart hela urvalet och måtten i manuset.
  choices <- combn(2:10, 4L)
  partitions <- lapply(seq_len(ncol(choices)), function(i) {
    a <- c(1L, choices[, i]); list(A = a, B = setdiff(1:10, a))
  })
  row_mean <- function(x) { v <- rowMeans(x, na.rm = TRUE); v[is.nan(v)] <- NA_real_; v }
  half_scores <- function(conf, correct, indices) {
    conf <- conf[, indices, drop = FALSE]; correct <- correct[, indices, drop = FALSE]
    bias <- row_mean(conf) - row_mean(correct)
    right <- wrong <- conf
    right[is.na(correct) | correct != 1] <- NA_real_
    wrong[is.na(correct) | correct != 0] <- NA_real_
    cbind(bias = bias, abs_bias = abs(bias), slope = row_mean(right) - row_mean(wrong))
  }
  alpha_rows <- split_rows <- list()
  for (domain in names(domains)) {
    cat("Reliability: ", domain, " — 126 balanced splits\n", sep = "")
    spec <- domains[[domain]]
    conf <- as.matrix(data[spec$confidence]) / 100
    correct <- as.matrix(data[spec$correct])
    total <- as.matrix(data[spec$all_correct]); total[is.na(total)] <- 0
    inputs <- list(total_acc = total, acc = correct, conf = conf)
    for (measure in names(inputs)) {
      x <- inputs[[measure]]
      eligible <- if (measure == "total_acc") is.finite(data[[paste0(spec$prefix, "_total_acc")]]) else rep(TRUE, nrow(x))
      complete <- eligible & rowSums(is.finite(x)) == ncol(x)
      used <- x[complete, , drop = FALSE]
      alpha <- NA_real_; n_constant <- NA_integer_; status <- "OK"
      if (nrow(used) < 3L) status <- "Fewer than three complete participants" else {
        item_var <- apply(used, 2L, var); n_constant <- sum(item_var == 0)
        sum_var <- var(rowSums(used))
        if (!is.finite(sum_var) || sum_var == 0) status <- "No total-score variance" else {
          k <- ncol(used); alpha <- k / (k - 1) * (1 - sum(item_var) / sum_var)
        }
      }
      alpha_rows[[length(alpha_rows) + 1L]] <- data.frame(Domain = domain, Measure = measure,
        n_items = ncol(x), n_input = nrow(x), n_eligible = sum(eligible),
        n_used = sum(complete), n_constant_items = n_constant, Alpha_raw = alpha, Status = status)
    }
    for (split in seq_along(partitions)) {
      a <- half_scores(conf, correct, partitions[[split]]$A)
      b <- half_scores(conf, correct, partitions[[split]]$B)
      for (measure in colnames(a)) {
        ok <- is.finite(a[, measure]) & is.finite(b[, measure]); n <- sum(ok)
        icc <- lower <- upper <- NA_real_; status <- "OK"
        if (n < 3L) status <- "Fewer than three participants with both halves estimable" else {
          x <- cbind(a[ok, measure], b[ok, measure])
          ms_person <- 2 * var(rowMeans(x))
          residual <- sweep(sweep(x, 1L, rowMeans(x)), 2L, colMeans(x)) + mean(x)
          ms_error <- sum(residual^2) / (n - 1)
          if (ms_person + ms_error == 0) status <- "No participant or residual variance" else {
            icc <- (ms_person - ms_error) / (ms_person + ms_error)
            if (ms_error == 0) lower <- upper <- 1 else {
              f <- ms_person / ms_error
              f_lower <- f / qf(.975, n - 1, n - 1); f_upper <- f * qf(.975, n - 1, n - 1)
              lower <- (f_lower - 1) / (f_lower + 1); upper <- (f_upper - 1) / (f_upper + 1)
            }
          }
        }
        split_rows[[length(split_rows) + 1L]] <- data.frame(Domain = domain,
          Measure = measure, Split = split, items_A = paste(partitions[[split]]$A, collapse = ","),
          items_B = paste(partitions[[split]]$B, collapse = ","), n_input = nrow(data),
          n_A = sum(is.finite(a[, measure])), n_B = sum(is.finite(b[, measure])),
          n_complete = n, ICC_C_half = icc, lower95 = lower, upper95 = upper, Status = status)
      }
    }
  }
  details <- do.call(rbind, split_rows)
  keys <- unique(details[c("Domain", "Measure")])
  summaries <- lapply(seq_len(nrow(keys)), function(i) {
    x <- details[details$Domain == keys$Domain[i] & details$Measure == keys$Measure[i], ]
    v <- x$ICC_C_half[is.finite(x$ICC_C_half)]
    data.frame(keys[i, ], n_splits = nrow(x), n_splits_estimable = length(v),
      n_complete_min = min(x$n_complete), n_complete_median = median(x$n_complete),
      n_complete_max = max(x$n_complete), ICC_half_median = if (length(v)) median(v) else NA_real_,
      split_q025 = if (length(v)) unname(quantile(v, .025)) else NA_real_,
      split_q975 = if (length(v)) unname(quantile(v, .975)) else NA_real_)
  })
  list(alpha = do.call(rbind, alpha_rows), split_summary = do.call(rbind, summaries),
       split_details = details, partitions = partitions)
}

.run_manuscript <- function(data_dir = getOption("buratti.data_dir", "/safe/data/Buratti"),
                            output_dir = file.path(data_dir, "manus_output")) {
  packages <- c("lmerTest", "emmeans", "lavaan")
  missing_packages <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_packages)) stop("Installera först följande paket: ", paste(missing_packages, collapse = ", "), call. = FALSE)
  data_path <- file.path(data_dir, "clean_all.rds")
  if (!file.exists(data_path)) stop("clean_all.rds saknas i ", data_dir,
    ". Kör 01.clean_all.R eller ange options(buratti.data_dir = ...).", call. = FALSE)
  data <- readRDS(data_path)
  domains <- .manus_domains()
  measures <- c("acc", "conf", "bias", "abs_bias", "slope", "total_acc")
  required <- c("id", "Cohort", "RSSEX", "Cohort01", "Sex01",
    unlist(lapply(domains, function(x) c(x$n_answered, x$n_conf, x$confidence, x$correct,
      x$all_correct, paste0(x$prefix, "_", measures))), use.names = FALSE))
  missing <- setdiff(required, names(data))
  if (length(missing)) stop("Kolumner saknas i clean_all.rds: ", paste(missing, collapse = ", "),
    ". Kör den uppdaterade 01.clean_all.R först.", call. = FALSE)
  if (anyNA(data$id) || anyDuplicated(data$id)) stop("id måste vara unikt och finnas för varje deltagarrad.", call. = FALSE)
  if (any(!is.na(data$Cohort) & !data$Cohort %in% c(1948, 1953)) ||
      any(!is.na(data$RSSEX) & !data$RSSEX %in% c(1, 2))) stop("Oväntade kohort- eller könskoder; kontrollera data.", call. = FALSE)
  expected_cohort <- ifelse(data$Cohort == 1948, 0, 1)
  expected_sex <- ifelse(data$RSSEX == 1, 0, 1)
  if (!isTRUE(all.equal(as.numeric(data$Cohort01), expected_cohort, check.attributes = FALSE)) ||
      !isTRUE(all.equal(as.numeric(data$Sex01), expected_sex, check.attributes = FALSE)))
    stop("Cohort01 eller Sex01 stämmer inte med kodningen i 01. Kör den uppdaterade 01-filen.", call. = FALSE)
  numeric_cols <- setdiff(required, c("id", "Cohort", "RSSEX"))
  if (!all(vapply(data[numeric_cols], is.numeric, logical(1)))) stop("Analys- och uppgiftskolumnerna måste vara numeriska.", call. = FALSE)
  if (any(vapply(data[numeric_cols], function(x) any(is.infinite(x)), logical(1)))) stop("Oändliga värden finns i analysdata; kontrollera dem först.", call. = FALSE)
  # En lokal kopia används. Inget sparas tillbaka till clean_all.rds.
  data$Cohort01 <- as.numeric(data$Cohort01); data$Sex01 <- as.numeric(data$Sex01)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(output_dir)) stop("Outputmappen kunde inte skapas: ", output_dir, call. = FALSE)
  log_file <- file.path(output_dir, "output_04_manus.txt")
  connection <- file(log_file, open = "wt", encoding = "UTF-8")
  previous_depth <- sink.number()
  on.exit({ while (sink.number() > previous_depth) sink(); close(connection) }, add = TRUE)
  sink(connection, split = TRUE)
  tables <- models <- model_counts <- list()
  report <- function(name, value, print_result = TRUE) {
    value <- as.data.frame(value)
    tables[[name]] <<- value
    write.csv(value, file.path(output_dir, paste0(name, ".csv")), row.names = FALSE, fileEncoding = "UTF-8", na = "NA")
    if (print_result) { cat("\n", name, "\n", sep = ""); print(value, row.names = FALSE, digits = 6) }
    flush(connection)
    invisible(value)
  }
  section <- function(title) { cat("\n\n", strrep("=", 72), "\n", title, "\n", sep = ""); flush(connection) }
  register_regression <- function(result, label) {
    models[[label]] <<- result$model; model_counts[[label]] <<- result$n
    report(paste0(label, "_coefficients"), result$coefficients)
    report(paste0(label, "_fit"), result$fit)
  }
  mixed_model <- function(formula, label, cohort_comparison = FALSE) {
    fit <- lmerTest::lmer(formula, data = long, REML = TRUE, na.action = na.omit)
    models[[label]] <<- fit
    model_counts[[label]] <<- .manus_model_n(fit, label, nrow(long))
    report(paste0(label, "_n"), model_counts[[label]])
    mf <- model.frame(fit)
    groups <- unique(mf[c("Domain", "Cohort01")])
    counts <- lapply(seq_len(nrow(groups)), function(i) {
      x <- mf[mf$Domain == groups$Domain[i] & mf$Cohort01 == groups$Cohort01[i], , drop = FALSE]
      data.frame(Domain = groups$Domain[i], Cohort = if (groups$Cohort01[i] == 0) 1948 else 1953,
        n_observations = nrow(x), n_unique_participants = length(unique(x$id)))
    })
    report(paste0(label, "_n_by_domain_cohort"), do.call(rbind, counts))
    omnibus <- anova(fit, type = 3, ddf = "Satterthwaite")
    report(paste0(label, "_anova"), cbind(Term = rownames(omnibus), as.data.frame(omnibus)))
    # Asymptotiska EMM/kontraster motsvarar outputen från den tidigare 04-körningen.
    # N för varje domän/kohort redovisas separat ovan; EMM är modellskattade medel.
    if (cohort_comparison) {
      emm <- emmeans::emmeans(fit, ~ Cohort01 | Domain, at = list(Cohort01 = c(0, 1)), lmer.df = "asymptotic")
      means <- as.data.frame(confint(emm, level = .95, adjust = "none"))
      means$Cohort <- ifelse(means$Cohort01 == 0, 1948, 1953)
      contrast <- emmeans::contrast(emm, "pairwise", adjust = "holm")
      differences <- as.data.frame(summary(contrast, infer = c(TRUE, TRUE), level = .95, adjust = "none"))
      # En kontrast i varje by-grupp: Holm-p är identiskt med ojusterat p här.
      differences$p_adjustment <- "Holm within domain (one contrast; no adjustment across domains)"
      differences$CI_adjustment <- "none"
    } else {
      emm <- emmeans::emmeans(fit, ~ Domain, lmer.df = "asymptotic")
      means <- as.data.frame(confint(emm, level = .95, adjust = "none"))
      contrast <- emmeans::contrast(emm, "pairwise", adjust = "tukey")
      differences <- as.data.frame(summary(contrast, infer = c(TRUE, TRUE), level = .95, adjust = "tukey"))
      differences$p_adjustment <- differences$CI_adjustment <- "Tukey over three domain contrasts"
    }
    report(paste0(label, "_EMM"), means)
    report(paste0(label, "_contrasts"), differences)
    report(paste0(label, "_diagnostics"), data.frame(Model = label,
      singular_fit = lme4::isSingular(fit), convergence_messages = paste(fit@optinfo$conv$lme4$messages, collapse = "; ")))
    invisible(fit)
  }
  withCallingHandlers({
    cat("Manuscript analyses H1–H4 started: ", format(Sys.time()), "\n", sep = "")
    cat("Data: ", normalizePath(data_path, winslash = "/"), "\n", sep = "")
    cat("Output: ", normalizePath(output_dir, winslash = "/"), "\n", sep = "")
    cat("N = unique participants; n_observations can include repeated domains.\n")
    cat("95% CI unless labelled RMSEA (90%). No derived measures are replaced.\n")
    section("PARTICIPANTS AND DESCRIPTIVE STATISTICS — TABLES 1 AND 7")
    report("participants_total", data.frame(n_rows = nrow(data), n_unique_participants = length(unique(data$id))))
    report("participants_cohort_sex", as.data.frame(table(Cohort = data$Cohort, RSSEX = data$RSSEX, useNA = "ifany")))
    report("descriptives", .manus_describe(data, domains))
    section("CONFIDENCE AVAILABILITY — SUPPLEMENTARY TABLE S1")
    cat("Logistic regression: unadjusted coefficient p-values and Wald 95% CIs for b and odds ratios.\n")
    availability <- list()
    for (domain in names(domains)) {
      spec <- domains[[domain]]
      eligible <- !is.na(data[[spec$n_answered]]) & data[[spec$n_answered]] >= 1
      indicator <- rep(NA_real_, nrow(data))
      indicator[eligible] <- as.numeric(data[[spec$n_conf]][eligible] >= 1)
      # Varje domän använder sina egna svar/confidence. Ingen kohort filtreras bort.
      data[[paste0(spec$prefix, "_has_conf")]] <- indicator
      for (group in c("All", "1948", "1953")) {
        selected <- if (group == "All") rep(TRUE, nrow(data)) else !is.na(data$Cohort) & data$Cohort == as.numeric(group)
        ne <- sum(eligible & selected); nc <- sum(indicator[selected] == 1, na.rm = TRUE)
        availability[[length(availability) + 1L]] <- data.frame(Group = group, Domain = domain,
          n_total = sum(selected), n_eligible = ne, n_no_eligible_accuracy = sum(selected) - ne,
          n_confidence_available = nc, n_without_confidence = sum(indicator[selected] == 0, na.rm = TRUE),
          n_eligible_availability_unknown = sum(eligible & selected & is.na(indicator)),
          pct_available = if (ne) 100 * nc / ne else NA_real_)
      }
    }
    report("S1_confidence_availability", do.call(rbind, availability))
    indicators <- paste0(vapply(domains, function(x) x$prefix, character(1)), "_has_conf")
    report("S1_confidence_overlap", data.frame(n_total = nrow(data),
      n_with_confidence_all_three = sum(rowSums(data[indicators] == 1, na.rm = TRUE) == 3)))
    for (domain in names(domains)) {
      prefix <- domains[[domain]]$prefix; label <- paste0("S1_availability_", prefix)
      formula <- as.formula(paste0(prefix, "_has_conf ~ ", prefix, "_total_acc + Sex01 + Cohort01"))
      register_regression(.manus_regression(data, formula, label, logistic = TRUE), label)
    }
    section("H1/H2 CROSS-DOMAIN CORRELATIONS — TABLE 2")
    correlations <- .manus_correlations(data, domains)
    report("H1_H2_correlations", correlations)
    report("H1_H2_mean_r", do.call(rbind, lapply(unique(correlations$Measure), function(measure) {
      v <- correlations$r[correlations$Measure == measure]
      data.frame(Measure = measure, n_pairs_estimable = sum(is.finite(v)),
        arithmetic_mean_r = if (all(is.finite(v))) mean(v) else NA_real_)
    })))
    cat("Pair-specific n above. Holm p over three pairs per measure; unadjusted Fisher 95% CIs.\n")
    section("H2 CONFIDENCE VERSUS TOTAL ACCURACY — PARTICIPANT BOOTSTRAP")
    cat("5,000 resamples; same row draw for confidence and accuracy; pairwise missingness.\n")
    bootstrap <- .manus_bootstrap(data, domains)
    report("H2_bootstrap", bootstrap$summary)
    cat("Difference = mean Fisher z(confidence) minus mean Fisher z(total accuracy).\n")
    cat("The three correlations have different pair-specific n; they are not treated as independent in resampling.\n")
    long <- .manus_long(data, domains)
    section("H3 DOMAIN DIFFERENCES — TABLE 3")
    cat("Primary models: sex/cohort adjustment, random participant intercept; REML, Type III Satterthwaite F.\n")
    for (outcome in c("Bias", "Absolute_bias", "Slope")) {
      mixed_model(as.formula(paste0(outcome, " ~ Domain + Sex01 + Cohort01 + (1 | id)")), paste0("H3_", outcome))
    }
    section("H3 SENSITIVITY: DOMAIN-SPECIFIC TEN-ITEM ACCURACY")
    cat("Accuracy comes from each row's domain. Accuracy is part of bias; interpret these as sensitivity models.\n")
    for (outcome in c("Bias", "Absolute_bias", "Slope")) {
      mixed_model(as.formula(paste0(outcome, " ~ Domain + Accuracy + Sex01 + Cohort01 + (1 | id)")), paste0("H3_accuracy_adjusted_", outcome))
    }
    section("H3 TOTAL ACCURACY AND MONITORING WITHIN DOMAINS — TABLE 4")
    cat("Linear regression: unadjusted coefficient p-values and t-based 95% CIs, with residual df and R-squared.\n")
    for (measure in c("abs_bias", "slope")) for (domain in names(domains)) {
      prefix <- domains[[domain]]$prefix; label <- paste0("H3_performance_", prefix, "_", measure)
      formula <- as.formula(paste0(prefix, "_", measure, " ~ ", prefix, "_total_acc + Sex01 + Cohort01"))
      register_regression(.manus_regression(data, formula, label), label)
    }
    section("H4 COHORT CONTRASTS WITHIN EACH DOMAIN — TABLES 5 AND 6")
    cat("1948 minus 1953; sex adjusted; accuracy is not included in these models.\n")
    cat("One contrast per domain: no multiplicity adjustment across the three domains, as in 04.\n")
    for (outcome in c("Absolute_bias", "Slope")) {
      mixed_model(as.formula(paste0(outcome, " ~ Domain * Cohort01 + Sex01 + (1 | id)")), paste0("H4_", outcome), cohort_comparison = TRUE)
    }
    section("H4 PERFORMANCE/CONFIDENCE CONTEXT — TABLE 7")
    for (measure in c("total_acc", "conf")) for (domain in names(domains)) {
      prefix <- domains[[domain]]$prefix; label <- paste0("H4_context_", prefix, "_", measure)
      formula <- as.formula(paste0(prefix, "_", measure, " ~ Cohort01 + Sex01"))
      register_regression(.manus_regression(data, formula, label), label)
    }
    section("H1 FACTOR MODELS — SUPPLEMENTARY TABLE S2")
    cat("Three existing specifications, ML/FIML. A three-indicator saturated model is not used as evidence of fit.\n")
    specifications <- list(
      S2_abs_bias_two_factors = paste(
        "metacog_abs =~ atonyms_abs_bias + metal_abs_bias + number_abs_bias",
        "cog_g =~ atonyms_total_acc + metal_total_acc + number_total_acc", "metacog_abs ~~ cog_g", sep = "\n"),
      S2_slope_two_factors = paste(
        "metacog_slope =~ atonyms_slope + metal_slope + number_slope",
        "cog_g =~ atonyms_total_acc + metal_total_acc + number_total_acc", "metacog_slope ~~ cog_g", sep = "\n"),
      S2_slope_one_factor = paste("general_factor =~ atonyms_slope + metal_slope + number_slope +",
        "atonyms_total_acc + metal_total_acc + number_total_acc"))
    for (label in names(specifications)) {
      cat("\n", label, "\n", sep = "")
      fit <- lavaan::cfa(specifications[[label]], data = data, estimator = "ML", missing = "fiml")
      models[[label]] <- fit
      converged <- lavaan::lavInspect(fit, "converged")
      report(paste0(label, "_status"), data.frame(Model = label,
        n_used_participants = sum(lavaan::lavInspect(fit, "nobs")), converged = converged))
      vars <- lavaan::lavNames(fit, "ov")
      report(paste0(label, "_indicator_n"), data.frame(Variable = vars,
        n_observed = vapply(data[vars], function(x) sum(is.finite(x)), integer(1))))
      if (!converged) { warning(label, ": did not converge; no fit indices interpreted."); next }
      print(summary(fit, fit.measures = TRUE, standardized = TRUE))
      indices <- lavaan::fitMeasures(fit)
      requested <- c("chisq", "df", "pvalue", "cfi", "tli", "rmsea", "rmsea.ci.lower", "rmsea.ci.upper", "srmr", "aic", "bic",
        "cfi.robust", "tli.robust", "rmsea.robust", "rmsea.ci.lower.robust", "rmsea.ci.upper.robust")
      report(paste0(label, "_fit"), data.frame(Index = requested,
        Value = as.numeric(indices[requested])))
      report(paste0(label, "_parameters"), lavaan::parameterEstimates(fit, standardized = TRUE, ci = TRUE, level = .95))
      report(paste0(label, "_standardized"), lavaan::standardizedSolution(fit, ci = TRUE, level = .95))
    }
    section("INTERNAL CONSISTENCY AND SPLIT-HALF — SUPPLEMENTARY TABLE S3")
    cat("Alpha: full-test missing items are zero as in total accuracy; ten-item accuracy/confidence use complete cases.\n")
    cat("Derived scores: all 126 five-plus-five splits, consistency ICC(C,1), no Spearman–Brown conversion.\n")
    reliability <- .manus_reliability(data, domains)
    report("S3_alpha", reliability$alpha)
    report("S3_split_half", reliability$split_summary)
    report("S3_split_half_details", reliability$split_details, print_result = FALSE)
    cat("Split percentiles are variation over splits, not a 95% CI. Each split's own CI is in S3_split_half_details.csv.\n")
    section("MODEL SAMPLE SIZES")
    report("all_model_n", do.call(rbind, model_counts))
    cat("Mixed models: counts include repeated domain observations; participants are counted separately.\n")
    cat("CFA participant counts are in the corresponding S2_*_status tables (FIML, not complete-case n).\n")
    section("REPRODUCIBILITY")
    print(sessionInfo())
    cat("Global contrasts option: "); print(getOption("contrasts"))
    results <- list(tables = tables, models = models, bootstrap = bootstrap,
      reliability_partitions = reliability$partitions, data_path = normalizePath(data_path, winslash = "/"),
      input_md5 = unname(tools::md5sum(data_path)), session_info = sessionInfo(),
      contrasts = getOption("contrasts"), created = Sys.time())
    saveRDS(results, file.path(output_dir, "manus_results.rds"))
    cat("\nManuscript analyses completed: ", format(Sys.time()), "\n", sep = "")
    cat("Output saved to: ", normalizePath(log_file, winslash = "/"), "\n", sep = "")
    flush(connection)
    invisible(results)
  }, warning = function(w) { writeLines(paste0("WARNING: ", conditionMessage(w)), connection); flush(connection) },
  message = function(m) { writeLines(paste0("MESSAGE: ", conditionMessage(m)), connection); flush(connection) },
  error = function(e) { writeLines(paste0("ERROR: ", conditionMessage(e)), connection); flush(connection) })
}

# Source kör alla avsnitt en gång och skapar ett eget objekt för manusresultaten.
manus_results <- .run_manuscript()

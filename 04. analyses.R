# 04. analyses.R — försiktigt städad version
# Beräkningar, modellinställningar och körordning är bevarade.
# Kända problem är markerade med GRANSKA och beskrivna i README.txt.
# Samma datamapp som i 01.clean_all.R. Anpassa via options() vid behov:
# options(buratti.data_dir = "C:/sökväg/till/dina/data")
# Automatisk paketinstallation och rm(list = ls()) har tagits bort.
# Installera vid behov i konsolen (en gång):
# install.packages(c("dplyr", "psych", "effectsize", "lavaan",
#                    "lmerTest", "emmeans"))

required_packages <- c("dplyr", "psych", "effectsize", "lavaan", "lmerTest", "emmeans")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0L) {
  stop("Installera först följande paket: ", paste(missing_packages, collapse = ", "), call. = FALSE)
}

data_dir <- getOption("buratti.data_dir", "/safe/data/Buratti")
df <- readRDS(file.path(data_dir, "clean_all.rds"))

# Kontrollera kopplingen innan de första analyserna körs.
required_columns <- c(
  "Sex01", "Cohort01", "atonym_n_conf", "metal_n_conf", "number_n_conf",
  "atonym_n_answered", "metal_n_answered", "number_n_answered"
)
missing_columns <- setdiff(required_columns, names(df))
if (length(missing_columns) > 0L) {
  stop(
    "Variabler saknas i clean_all.rds: ", paste(missing_columns, collapse = ", "),
    ". Kör först den uppdaterade 01.clean_all.R i samma datamapp.",
    call. = FALSE
  )
}

c("atonyms_total_acc",
  "metal_total_acc",
  "number_total_acc"
) %in% names(df)

# ------------------------------------------------------------
# Kontroll av bias när confidence och accuracy har olika uppgiftsunderlag
# ------------------------------------------------------------
# Nuvarande bias = tillgänglig medelconfidence/100 minus tillgänglig accuracy.
# Matchad bias använder ENBART uppgifter där både confidence och rätt/fel finns.
# När samma uppgifter ingår blir måtten identiska, bortsett från avrundningsfel.
# Jämförelsen visar känsligheten för bortfall; matchad bias är inte ett facit.
# Befintliga mått, totalpoäng, deltagarfilter, slope och gamma ändras inte.
# Resultaten sparas separat i calibration_check. Inga kolumner i df skrivs över.

check_calibration_overlap <- function(data) {
  domains <- list(
    Antonyms = list(
      prefix = "atonyms",
      confidence = paste0("m_con_", 1:10),
      correct = paste0("mCat_", c(24, 31, 35, 37, 40, 41, 42, 43, 44, 45), "_correct")
    ),
    `Metal Folding` = list(
      prefix = "metal",
      confidence = paste0("p_con_", 1:10),
      correct = paste0("p_Q", c(10, 12, 16, 19, 20, 21, 24, 26, 27, 28), "_correct")
    ),
    `Number Series` = list(
      prefix = "number",
      confidence = paste0("tcat_t", c(2, 10, 11, 13, 14, 15, 17, 18, 19, 22)),
      correct = sprintf("Tal%02d_correct", c(2, 10, 11, 13, 14, 15, 17, 18, 19, 22))
    )
  )
  required <- c("Cohort", unlist(lapply(domains, function(spec) {
    c(spec$confidence, spec$correct, paste0(spec$prefix, c("_bias", "_abs_bias")))
  }), use.names = FALSE))
  missing <- setdiff(required, names(data))
  if (length(missing)) {
    stop("Kalibreringskontrollen saknar kolumner: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  # En person utan matchade uppgifter får NA, inte ett värde beräknat på tomma data.
  safe_row_mean <- function(x) {
    result <- rowMeans(x, na.rm = TRUE)
    result[is.nan(result)] <- NA_real_
    result
  }
  people <- list()
  items <- list()
  for (domain in names(domains)) {
    spec <- domains[[domain]]
    confidence <- as.matrix(data[spec$confidence]) / 100
    correct <- as.matrix(data[spec$correct])
    paired <- !is.na(confidence) & !is.na(correct)
    answered_without_conf <- !is.na(correct) & is.na(confidence)
    conf_without_answer <- !is.na(confidence) & is.na(correct)
    matched_conf <- confidence
    matched_acc <- correct
    matched_conf[!paired] <- NA_real_
    matched_acc[!paired] <- NA_real_
    matched_bias <- safe_row_mean(matched_conf) - safe_row_mean(matched_acc)
    current_bias <- data[[paste0(spec$prefix, "_bias")]]
    current_abs_bias <- data[[paste0(spec$prefix, "_abs_bias")]]
    people[[domain]] <- data.frame(
      Row = seq_len(nrow(data)),
      id = if ("id" %in% names(data)) data$id else seq_len(nrow(data)),
      Domain = domain, Cohort = data$Cohort,
      n_answered = rowSums(!is.na(correct)),
      n_conf = rowSums(!is.na(confidence)),
      n_pairs = rowSums(paired),
      n_answered_without_conf = rowSums(answered_without_conf),
      n_conf_without_answer = rowSums(conf_without_answer),
      bias_current = current_bias,
      bias_matched = matched_bias,
      abs_bias_current = current_abs_bias,
      abs_bias_matched = abs(matched_bias),
      # Delta i procentenheter: positivt = högre värde med matchade uppgifter.
      delta_bias_pp = 100 * (matched_bias - current_bias),
      delta_abs_bias_pp = 100 * (abs(matched_bias) - current_abs_bias)
    )
    # Uppgiftsvis bortfall per kohort kan visa t.ex. om mCat_43_correct/m_con_8
    # har systematiskt olika underlag mellan kohorter. Saknad kohort visas som NA.
    cohort_items <- lapply(unique(data$Cohort), function(cohort) {
      rows <- if (is.na(cohort)) is.na(data$Cohort) else !is.na(data$Cohort) & data$Cohort == cohort
      data.frame(
        Domain = domain, Cohort = cohort,
        correct_item = spec$correct, confidence_item = spec$confidence,
        n_participants = sum(rows),
        n_answered = colSums(!is.na(correct[rows, , drop = FALSE])),
        n_pairs = colSums(paired[rows, , drop = FALSE]),
        n_answered_without_conf = colSums(answered_without_conf[rows, , drop = FALSE]),
        n_conf_without_answer = colSums(conf_without_answer[rows, , drop = FALSE]),
        row.names = NULL
      )
    })
    items[[domain]] <- do.call(rbind, cohort_items)
  }
  by_person <- do.call(rbind, people)
  rownames(by_person) <- NULL
  summarize <- function(x) {
    # Alla förändringsmått använder samma personer med båda definitioner tillgängliga.
    # n_comparable kan vara mindre än n_current om inga matchade uppgifter finns.
    ok <- is.finite(x$bias_current) & is.finite(x$bias_matched) &
      is.finite(x$abs_bias_current) & is.finite(x$abs_bias_matched)
    delta <- x$delta_bias_pp[ok]
    abs_delta <- x$delta_abs_bias_pp[ok]
    mismatch <- (x$n_answered_without_conf + x$n_conf_without_answer)[ok] > 0
    data.frame(
      n_participants = nrow(x),
      n_current = sum(is.finite(x$bias_current)),
      n_matched = sum(is.finite(x$bias_matched)),
      n_comparable = sum(ok),
      n_different_item_sets = sum(mismatch),
      pct_different_item_sets = if (length(delta)) 100 * mean(mismatch) else NA_real_,
      # Positiva och negativa förändringar kan ta ut varandra i mean_delta_bias_pp.
      mean_delta_bias_pp = if (length(delta)) mean(delta) else NA_real_,
      # Medel av förändringens storlek; detta är INTE förändringen i absolut bias.
      mean_size_delta_bias_pp = if (length(delta)) mean(abs(delta)) else NA_real_,
      p95_size_delta_bias_pp = if (length(delta)) unname(quantile(abs(delta), 0.95)) else NA_real_,
      max_size_delta_bias_pp = if (length(delta)) max(abs(delta)) else NA_real_,
      # Absolut bias måste jämföras separat eftersom absolutbeloppet är icke-linjärt.
      mean_delta_abs_bias_pp = if (length(abs_delta)) mean(abs_delta) else NA_real_,
      mean_size_delta_abs_bias_pp = if (length(abs_delta)) mean(abs(abs_delta)) else NA_real_
    )
  }
  by_domain <- do.call(rbind, lapply(names(domains), function(domain) {
    cbind(Domain = domain, summarize(by_person[by_person$Domain == domain, ]))
  }))
  by_cohort <- do.call(rbind, lapply(names(domains), function(domain) {
    do.call(rbind, lapply(unique(data$Cohort), function(cohort) {
      rows <- by_person$Domain == domain &
        if (is.na(cohort)) is.na(by_person$Cohort) else !is.na(by_person$Cohort) & by_person$Cohort == cohort
      cbind(Domain = domain, Cohort = cohort, summarize(by_person[rows, ]))
    }))
  }))
  list(by_domain = by_domain, by_cohort = by_cohort,
       by_person = by_person, by_item_cohort = do.call(rbind, items))
}

calibration_check <- check_calibration_overlap(df)
# Börja med dessa tabeller. Förändringar anges i PROCENTENHETER, inte procent.
# Noll/liten skillnad visar att bortfallshanteringen har liten effekt på själva måttet.
# Om skillnaderna är större eller kohortspecifika bör även modeller jämföras separat.
# Dessa tabeller visar inte i sig om regressionsresultat eller slutsatser ändras.
print(calibration_check$by_domain, row.names = FALSE)
print(calibration_check$by_cohort, row.names = FALSE)
# Mer detaljer vid behov (kör valfri rad interaktivt):
# View(calibration_check$by_person)
# View(calibration_check$by_item_cohort)

# ------------------------------------------------------------
# Intern konsistens och split-half-konsistens INOM varje domän
# ------------------------------------------------------------
# Intern reliabilitet undersöks med alfa för accuracy/confidence och
# upprepad split-half för bias, absolut bias, slope och gamma.
# Inga befintliga mått skrivs över. Alfa beräknas inte för kalibreringsmåtten.
# Alfa skattar intern konsistens för total accuracy, tiouppgifts-accuracy och
# medelconfidence. För total accuracy räknas saknade uppgifter som 0, precis
# som i din totalpoäng med fast nämnare; bara personer med tillgänglig
# totalpoäng ingår. För tiouppgifts-accuracy och confidence används kompletta
# fall på respektive tio kolumner. Därför kan urvalen skilja sig åt.
# Alfa är inte ett test av endimensionalitet eller stabilitet över tid.
#
# För bias, absolut bias, slope och gamma delas samma tio uppgifter i 5+5.
# Bias i varje halva använder separata tillgängliga medelvärden för confidence
# och accuracy, precis som dina nuvarande mått. Matchad bias införs inte här.
# Måttet räknas om i vardera halvan och halvorna jämförs med ICC(C,1):
# tvåvägs konsistens för en halvtestskattning. Nivåskillnader mellan halvorna
# tillåts. Detta görs för ALLA 126 unika balanserade uppdelningar av tio items.
# Samma uppdelning används för alla personer; uppgift 1 ligger alltid i A,
# vilket undviker att samma uppdelning räknas två gånger med omvända halvor.
#
# ICC avser HALVTESTEN med fem uppgifter, inte reliabiliteten för hela testet.
# Ingen automatisk Spearman–Brown-korrigering görs för dessa härledda mått:
# gamma, slope och absolut bias är inte enkla summor av likvärdiga items.
# Median och 2,5/97,5-percentiler sammanfattar variationen mellan uppdelningar;
# percentilerna är INTE ett 95 % konfidensintervall. split_details innehåller
# ett separat 95 % ICC-intervall per beräkningsbar uppdelning.
#
# Fem uppgifter per halva kan ge odefinierad gamma/slope, särskilt om alla
# tillgängliga svar är rätt/fel eller confidence saknas. Inga värden imputeras.
# n_complete visar personer med ett ändligt mått i BÅDA halvor; samma personer
# ingår i båda kolumnerna för varje ICC. Urvalet kan variera mellan uppdelningar.
# Resultaten visas totalt och per kohort; låg ICC kan avspegla få uppgifter,
# selektivt bortfall och verklig variation mellan uppgifter. Detta mäter
# split-half-konsistens för de härledda måtten. Den ska rapporteras tillsammans
# med antal beräkningsbara splits/personer och variationen över uppdelningar,
# inte som ett direkt reliabilitetsmått för hela tiouppgiftstestet.
# Upprepade testtillfällen eller samband mellan domäner analyseras inte här.
# Läsning: https://doi.org/10.3758/s13423-021-01948-3

check_within_test_reliability <- function(data) {
  domains <- list(
    Antonyms = list(prefix = "atonyms", confidence = paste0("m_con_", 1:10),
      correct = paste0("mCat_", c(24,31,35,37,40,41,42,43,44,45), "_correct"),
      all_correct = paste0("mCat_", 1:45, "_correct")),
    `Metal Folding` = list(prefix = "metal", confidence = paste0("p_con_", 1:10),
      correct = paste0("p_Q", c(10,12,16,19,20,21,24,26,27,28), "_correct"),
      all_correct = paste0("p_Q", 1:40, "_correct")),
    `Number Series` = list(prefix = "number", confidence = paste0("tcat_t", c(2,10,11,13,14,15,17,18,19,22)),
      correct = sprintf("Tal%02d_correct", c(2,10,11,13,14,15,17,18,19,22)),
      all_correct = sprintf("Tal%02d_correct", 1:40))
  )
  required <- c("Cohort", unlist(lapply(domains, function(spec) {
    c(spec$confidence, spec$correct, spec$all_correct, paste0(spec$prefix, "_total_acc"))
  }), use.names = FALSE))
  missing <- setdiff(required, names(data))
  if (length(missing)) stop("Reliabilitetskontrollen saknar kolumner: ", paste(missing, collapse = ", "), call. = FALSE)
  groups <- list(`Alla kohorter` = seq_len(nrow(data)))
  for (cohort in unique(data$Cohort)) {
    label <- if (is.na(cohort)) "Kohort saknas" else paste("Kohort", cohort)
    groups[[label]] <- which(if (is.na(cohort)) is.na(data$Cohort) else data$Cohort == cohort)
  }
  partitions <- lapply(seq_len(ncol(combn(2:10, 4))), function(i) {
    a <- c(1L, combn(2:10, 4)[, i])
    list(A = a, B = setdiff(1:10, a))
  })
  safe_mean <- function(x) {
    result <- mean(x, na.rm = TRUE)
    if (is.nan(result)) NA_real_ else result
  }
  score <- function(conf, correct) {
    bias <- safe_mean(conf) - safe_mean(correct)
    slope <- safe_mean(conf[correct == 1]) - safe_mean(conf[correct == 0])
    ok <- is.finite(conf) & is.finite(correct)
    x <- conf[ok]; y <- correct[ok]
    gamma <- NA_real_
    if (length(x) >= 2L) {
      pairs <- combn(seq_along(x), 2L)
      products <- (x[pairs[1, ]] - x[pairs[2, ]]) * (y[pairs[1, ]] - y[pairs[2, ]])
      C <- sum(products > 0); D <- sum(products < 0)
      if (C + D > 0L) gamma <- (C - D) / (C + D)
    }
    c(bias = bias, abs_bias = abs(bias), slope = slope, gamma = gamma)
  }
  half_scores <- function(conf, correct, indices) {
    t(vapply(seq_len(nrow(conf)), function(i) score(conf[i, indices], correct[i, indices]), numeric(4)))
  }
  alpha_rows <- list()
  split_rows <- list()
  for (group in names(groups)) {
    group_data <- data[groups[[group]], , drop = FALSE]
    for (domain in names(domains)) {
      spec <- domains[[domain]]
      conf <- as.matrix(group_data[spec$confidence]) / 100
      correct <- as.matrix(group_data[spec$correct])
      total <- as.matrix(group_data[spec$all_correct])
      total[is.na(total)] <- 0
      alpha_inputs <- list(total_acc = total, acc = correct, conf = conf)
      for (measure in names(alpha_inputs)) {
        x <- alpha_inputs[[measure]]
        eligible <- if (measure == "total_acc") is.finite(group_data[[paste0(spec$prefix, "_total_acc")]]) else rep(TRUE, nrow(x))
        complete <- eligible & rowSums(is.finite(x)) == ncol(x)
        used <- x[complete, , drop = FALSE]
        alpha <- NA_real_
        status <- "OK"
        n_constant <- NA_integer_
        if (nrow(used) < 3L) {
          status <- "För få kompletta personer (minst 3 krävs här)"
        } else {
          item_var <- apply(used, 2, var)
          n_constant <- sum(item_var == 0)
          sum_var <- var(rowSums(used))
          if (!is.finite(sum_var) || sum_var == 0) {
            status <- "Ingen variation i summapoängen"
          } else {
            # Rått Cronbach-alfa; binära items ger samma koefficient som KR-20.
            # Items utan varians behålls för att motsvara den angivna poängen.
            k <- ncol(used)
            alpha <- k / (k - 1) * (1 - sum(item_var) / sum_var)
          }
        }
        alpha_rows[[length(alpha_rows) + 1L]] <- data.frame(
          Group = group, Domain = domain, Measure = measure, n_items = ncol(x),
          n_group = nrow(x), n_eligible = sum(eligible), n_used = sum(complete),
          n_constant_items = n_constant, Alpha_raw = alpha, Status = status)
      }
      for (split in seq_along(partitions)) {
        a <- half_scores(conf, correct, partitions[[split]]$A)
        b <- half_scores(conf, correct, partitions[[split]]$B)
        for (measure in colnames(a)) {
          ok <- is.finite(a[, measure]) & is.finite(b[, measure])
          n <- sum(ok)
          icc <- lower <- upper <- NA_real_
          status <- "OK"
          if (n < 3L) {
            status <- "För få personer med beräkningsbara mått i båda halvor"
          } else {
            x <- cbind(a[ok, measure], b[ok, measure])
            # ICC(C,1) = (MS_person - MS_error)/(MS_person + MS_error), k = 2.
            ms_person <- 2 * var(rowMeans(x))
            residual <- sweep(sweep(x, 1, rowMeans(x)), 2, colMeans(x)) + mean(x)
            ms_error <- sum(residual^2) / (n - 1)
            if (ms_person + ms_error == 0) {
              status <- "ICC odefinierad: ingen person- eller residualvariation"
            } else {
              icc <- (ms_person - ms_error) / (ms_person + ms_error)
              if (ms_error == 0) {
                lower <- upper <- 1
              } else {
                # ANOVA/F-baserat 95 % intervall för just denna uppdelning.
                f <- ms_person / ms_error
                f_lower <- f / qf(.975, n - 1, n - 1)
                f_upper <- f * qf(.975, n - 1, n - 1)
                lower <- (f_lower - 1) / (f_lower + 1)
                upper <- (f_upper - 1) / (f_upper + 1)
              }
            }
          }
          split_rows[[length(split_rows) + 1L]] <- data.frame(
            Group = group, Domain = domain, Measure = measure, Split = split,
            n_group = nrow(group_data), n_A = sum(is.finite(a[, measure])),
            n_B = sum(is.finite(b[, measure])), n_complete = n,
            ICC_C_half = icc, lower95 = lower, upper95 = upper, Status = status)
        }
      }
    }
  }
  details <- do.call(rbind, split_rows)
  keys <- unique(details[c("Group", "Domain", "Measure")])
  summaries <- lapply(seq_len(nrow(keys)), function(i) {
    key <- keys[i, ]
    x <- details[details$Group == key$Group & details$Domain == key$Domain & details$Measure == key$Measure, ]
    values <- x$ICC_C_half[is.finite(x$ICC_C_half)]
    data.frame(key, n_splits = nrow(x), n_splits_estimable = length(values),
      n_complete_min = min(x$n_complete), n_complete_median = median(x$n_complete),
      n_complete_max = max(x$n_complete),
      ICC_half_median = if (length(values)) median(values) else NA_real_,
      split_q025 = if (length(values)) unname(quantile(values, .025)) else NA_real_,
      split_q975 = if (length(values)) unname(quantile(values, .975)) else NA_real_)
  })
  list(alpha = do.call(rbind, alpha_rows), split_summary = do.call(rbind, summaries),
       split_details = details, partitions = partitions)
}

within_reliability <- check_within_test_reliability(df)
cat("\nIntern konsistens: alfa för accuracy och confidence inom varje domän\n")
print(within_reliability$alpha, row.names = FALSE)
cat("\nSplit-half-konsistens: ICC mellan femuppgiftshalvor inom varje domän\n")
cat("Median/percentiler avser uppdelningarna, inte fulltestreliabilitet eller ett gemensamt konfidensintervall.\n")
print(within_reliability$split_summary, row.names = FALSE)
# Detaljer vid behov: View(within_reliability$split_details)
# En uppdelnings uppgiftspositioner: within_reliability$partitions[[1]]
# Uppgiftspositionerna avser de tio domänspecifika kolumnnamnen ovan.


####describe
desc <- psych::describe(df[,c(
      "atonyms_acc",
      "atonyms_conf",
      "atonyms_bias",
      "atonyms_abs_bias",
      "atonyms_gamma",
      "atonyms_slope",
      "atonyms_total_acc",
      "metal_acc",
      "metal_conf",
      "metal_bias",
      "metal_abs_bias",
      "metal_gamma",
      "metal_slope",
      "metal_total_acc",
      "number_acc",
      "number_conf",
      "number_bias",
      "number_abs_bias",
      "number_gamma",
      "number_slope",
      "number_total_acc"
)])

round(desc[, c("n", "mean", "sd", "min", "max")], 2)

library(dplyr)

desc <- df %>%
group_by(Sex01, Cohort01) %>%
summarise(
  across(
    c(
      atonyms_acc,
      atonyms_conf,
      atonyms_bias,
      atonyms_abs_bias,
      atonyms_slope,
      atonyms_total_acc,
      metal_acc,
      metal_conf,
      metal_bias,
      metal_abs_bias,
      metal_slope,
      metal_total_acc,
      number_acc,
      number_conf,
      number_bias,
      number_abs_bias,
      number_slope,
      number_total_acc
    ), list(
      n = ~sum(!is.na(.)),
      M = ~mean(., na.rm = TRUE),
      SD = ~sd(., na.rm = TRUE)
  )),.groups = "drop")

desc

desc <- df %>%
group_by(Cohort01) %>%
summarise(
  across(
    c(
      atonyms_acc,
      atonyms_conf,
      atonyms_bias,
      atonyms_abs_bias,
      atonyms_slope,
      atonyms_total_acc,
      metal_acc,
      metal_conf,
      metal_bias,
      metal_abs_bias,
      metal_slope,
      metal_total_acc,
      number_acc,
      number_conf,
      number_bias,
      number_abs_bias,
      number_slope,
      number_total_acc
    ), list(
      n = ~sum(!is.na(.)),
      M = ~mean(., na.rm = TRUE),
      SD = ~sd(., na.rm = TRUE)
  )),.groups = "drop")

print(desc, width = Inf)

#######################
sum(!is.na(df$atonyms_conf))

sum(!is.na(df$atonyms_total_acc) & is.na(df$atonyms_conf))

sum(!is.na(df$metal_conf))

sum(!is.na(df$metal_total_acc) & is.na(df$metal_conf))

sum(!is.na(df$number_conf))

sum(!is.na(df$number_total_acc) & is.na(df$number_conf))

#########################################################
table(df$atonym_n_answered >=1, useNA = "ifany")

table(df$atonym_n_conf >=1, useNA = "ifany")

table(
  answered = df$atonym_n_answered >=1,
  confidence = df$atonym_n_conf >=1,
  useNA = "ifany"
)

library(effectsize)

df$atonyms_has_conf <- NA

df$atonyms_has_conf[df$atonym_n_answered >= 1] <-
as.numeric(df$atonym_n_conf[df$atonym_n_answered >=1] >=1)

table(df$atonyms_has_conf, useNA = "ifany")

######################################
df$metal_has_conf <- NA

df$metal_has_conf[df$metal_n_answered >= 1] <-
as.numeric(df$metal_n_conf[df$metal_n_answered >=1] >=1)

table(df$metal_has_conf, useNA = "ifany")

#################################################
df$number_has_conf <- NA

df$number_has_conf[df$number_n_answered >= 1] <-
as.numeric(df$number_n_conf[df$number_n_answered >=1] >=1)

table(df$number_has_conf, useNA = "ifany")

####################################################
t.test(atonyms_total_acc ~ atonyms_has_conf, data = df)

cohens_d(atonyms_total_acc ~ atonyms_has_conf, data = df)

df %>%
group_by(atonyms_has_conf) %>%
summarise(
  n = n(),
  M = mean(atonyms_total_acc, na.rm = TRUE),
  SD = sd(atonyms_total_acc, na.rm = TRUE)
)

# Behåll indikatorn ovan: NA för personer utan besvarade Metal-uppgifter.
# Samma urvalsprincip gäller därmed för alla tre testen.

t.test(metal_total_acc ~ metal_has_conf, data = df)

cohens_d(metal_total_acc ~ metal_has_conf, data = df)

df %>%
group_by(metal_has_conf) %>%
summarise(
  n = n(),
  M = mean(metal_total_acc, na.rm = TRUE),
  SD = sd(metal_total_acc, na.rm = TRUE)
)

########################## Vad ger bortfallet? 
df$Cohort01 <- ifelse(df$Cohort == 1948, 0,1)

table(df$Cohort, df$Cohort01)

df$Sex01 <- ifelse(df$RSSEX == 1, 0,1)

table(df$RSSEX, df$Sex01)

summary(glm(atonyms_has_conf ~ atonyms_total_acc + Sex01 + Cohort01, data = df, family = binomial))

summary(glm(metal_has_conf ~ metal_total_acc + Sex01 + Cohort01, data = df, family = binomial))

########## Hur många har/har inte confidence per kohort?
table(df$Cohort, df$atonyms_has_conf)

prop.table(table(df$Cohort, df$atonyms_has_conf),1)

table(df$Cohort, df$metal_has_conf)

prop.table(table(df$Cohort, df$metal_has_conf),1)

########## Medelprestationen bland dem som har respektive saknar confidence
aggregate(atonyms_total_acc ~ atonyms_has_conf, data=df, mean, na.rm=TRUE )

aggregate(metal_total_acc ~ metal_has_conf, data=df, mean, na.rm=TRUE )

########## Är de samma personer som saknar confidence data
table(df$atonyms_has_conf, df$metal_has_conf)

prop.table(table(df$atonyms_has_conf, df$metal_has_conf))

####### kolla överlapp mellan metal items och konfidensbedömningar
metal_items <- paste0("p_Q", c(10,12,16,19,20,21,24,26,27,28), "_correct")

df$metal_n_answered <- rowSums(!is.na(df[metal_items]))

metal_conf_items <- paste0("p_con_", 1:10)

df$metal_n_conf <- rowSums(!is.na(df[metal_conf_items]))

table(df$metal_n_answered, df$metal_n_conf)

###### Missing confidence among participants ansvwering all 10 metal items ######
metal_complete_items <- df [df$metal_n_answered == 10, ]

table(metal_complete_items$metal_n_conf)

prop.table(table(metal_complete_items$metal_n_conf))

####### kolla överlapp mellan metal items och konfidensbedömningar
metal_items <- paste0("p_Q", c(10,12,16,19,20,21,24,26,27,28), "_correct")

df$metal_n_answered <- rowSums(!is.na(df[metal_items]))

metal_conf_items <- paste0("p_con_", 1:10)

df$metal_n_conf <- rowSums(!is.na(df[metal_conf_items]))

table(df$metal_n_answered, df$metal_n_conf)

####### kolla överlapp mellan atonym items och konfidensbedömningar
atonym_items <- paste0("mCat_", c(24,31,35,37,40,41,42,43,44,45), "_correct")

df$atonym_n_answered <- rowSums(!is.na(df[atonym_items]))

atonym_conf_items <- paste0("m_con_", 1:10)

df$atonym_n_conf <- rowSums(!is.na(df[atonym_conf_items]))

table(df$atonym_n_answered, df$atonym_n_conf)

##################################################
m_conf <- df[, paste0("m_con_", 1:10), drop = FALSE]

colSums(is.na(m_conf[df$atonym_n_answered == 10, , drop = FALSE]))

##################################Missing confidence for atonyms item 43
atonym_complete <- df[df$atonym_n_answered == 10, ]

table(atonym_complete$Cohort, is.na(atonym_complete$m_con_8))

prop.table(table(atonym_complete$Cohort, is.na(atonym_complete$m_con_8)), 1)

######ny bortfallsanalys där vi jämför de som har skattat med de som kunde skatta
df$atonym_has_conf_corrected <- NA

df$atonym_has_conf_corrected[df$atonym_n_answered >=1] <-
ifelse(
  df$atonym_n_conf[df$atonym_n_answered >=1] >=1,
  1, #har minste en confidence rating
  0 #kunde skatta men har ingen confidence rating
)

table(df$atonym_has_conf_corrected, useNA ="ifany")

prop.table(table(df$atonym_has_conf_corrected))

df$metal_has_conf_corrected <- NA

df$metal_has_conf_corrected[df$metal_n_answered >=1] <-
ifelse(
  df$metal_n_conf[df$metal_n_answered >=1] >=1,
  1, #har minste en confidence rating
  0 #kunde skatta men har ingen confidence rating
)

table(df$metal_has_conf_corrected, useNA ="ifany")

prop.table(table(df$metal_has_conf_corrected))

###############################Kolla deltta 
number_items <- paste0("Tal",
  sprintf("%02d", c(2,10,11,13,14,15,17,18,19,22)),
  "_correct")

conf_vars <- c("tcat_t2", "tcat_t10", "tcat_t11",
  "tcat_t13", "tcat_t14", "tcat_t15",
  "tcat_t17", "tcat_t18", "tcat_t19", "tcat_t22" )

df$number_n_answered <- rowSums(!is.na(df[number_items]))

df$number_n_conf <- rowSums(!is.na(df[conf_vars]))

table(df$number_n_answered, df$number_n_conf)

df$number_has_conf_corrected <- NA

df$number_has_conf_corrected[df$number_n_answered >=1] <-
ifelse(
  df$number_n_conf[df$number_n_answered >=1] >=1,
  1, #har minste en confidence rating
  0 #kunde skatta men har ingen confidence rating
)

#######bortfallsanalyser mellan de som kunde skatta konfidens men inte gjorde de och de som skattade
summary(glm(atonym_has_conf_corrected ~ atonyms_total_acc +Sex01 + Cohort01,
    data= df,
    family = binomial
))

summary(glm(metal_has_conf_corrected ~ metal_total_acc +Sex01 + Cohort01,
    data= df,
    family = binomial
))

summary(glm(number_has_conf_corrected ~ number_total_acc +Sex01 + Cohort01,
    data= df,
    family = binomial
))

table(df$number_has_conf_corrected, useNA = "ifany")

########## huvudanalyser
summary(lm(atonyms_abs_bias ~ atonyms_total_acc + Sex01 + Cohort01, data = df))

summary(lm(metal_abs_bias ~ metal_total_acc + Sex01 + Cohort01, data = df))

summary(lm(number_abs_bias ~ number_total_acc + Sex01 + Cohort01, data = df))

summary(lm(atonyms_slope ~ atonyms_total_acc + Sex01 + Cohort01, data = df))

summary(lm(metal_slope ~ metal_total_acc + Sex01 + Cohort01, data = df))

summary(lm(number_slope ~ number_total_acc + Sex01 + Cohort01, data = df))

#könseffekter
summary(lm(atonyms_total_acc ~ Sex01, data = df))

summary(lm(metal_total_acc ~ Sex01, data = df))

summary(lm(number_total_acc ~ Sex01, data = df))

summary(lm(atonyms_acc ~ Sex01, data = df))

summary(lm(metal_acc ~ Sex01, data = df))

summary(lm(number_acc ~ Sex01, data = df))

summary(lm(atonyms_conf ~ Sex01, data = df))

summary(lm(metal_conf ~ Sex01, data = df))

summary(lm(number_conf ~ Sex01, data = df))

summary(lm(atonyms_abs_bias ~ Sex01, data = df))

summary(lm(metal_abs_bias ~ Sex01, data = df))

summary(lm(number_abs_bias ~ Sex01, data = df))

summary(lm(atonyms_slope ~ Sex01, data = df))

summary(lm(metal_slope ~ Sex01, data = df))

summary(lm(number_slope ~ Sex01, data = df))

#Cohort
summary(lm(atonyms_total_acc ~ Cohort01, data = df))

summary(lm(metal_total_acc ~ Cohort01, data = df))

summary(lm(number_total_acc ~ Cohort01, data = df))

summary(lm(atonyms_acc ~ Cohort01, data = df))

summary(lm(metal_acc ~ Cohort01, data = df))

summary(lm(number_acc ~ Cohort01, data = df))

summary(lm(atonyms_conf ~ Cohort01, data = df))

summary(lm(metal_conf ~ Cohort01, data = df))

summary(lm(number_conf ~ Cohort01, data = df))

summary(lm(atonyms_abs_bias ~ Cohort01, data = df))

summary(lm(metal_abs_bias ~ Cohort01, data = df))

summary(lm(number_abs_bias ~ Cohort01, data = df))

summary(lm(atonyms_slope ~ Cohort01, data = df))

summary(lm(metal_slope ~ Cohort01, data = df))

summary(lm(number_slope ~ Cohort01, data = df))

summary(lm(atonyms_total_acc ~ Cohort01 + Sex01, data = df))

summary(lm(metal_total_acc ~ Cohort01 + Sex01, data = df))

summary(lm(number_total_acc ~ Cohort01 + Sex01, data = df))

summary(lm(atonyms_conf ~ Cohort01 + Sex01, data = df))

summary(lm(metal_conf ~ Cohort01 + Sex01, data = df))

summary(lm(number_conf ~ Cohort01 + Sex01, data = df))

library(lavaan)

#step 1 en metakognitiv faktor för absolute bias
model_abs <- 'metacog_abs =~ atonyms_abs_bias + metal_abs_bias + number_abs_bias'

fit_abs <- cfa(model_abs, data = df, missing = "fiml")

summary(fit_abs, fit.measures = TRUE, standardized = TRUE)

####Step 2
model_twofactor <- 'metacog_abs =~ atonyms_abs_bias + metal_abs_bias + number_abs_bias

cog_g =~ atonyms_total_acc + metal_total_acc+ number_total_acc

metacog_abs ~~ cog_g'

fit_twofactor <- cfa(model_twofactor, data = df, missing = "fiml")

summary(fit_twofactor, fit.measures = TRUE, standardized = TRUE)

############# Confidence
model_conf <- 'conf_f =~ atonyms_conf + metal_conf + number_conf

cog_g =~ atonyms_total_acc + metal_total_acc+ number_total_acc

conf_f ~~ cog_g'

fit_conf <- cfa(model_conf, data = df, missing = "fiml")

summary(fit_conf, fit.measures = TRUE, standardized = TRUE)

############# Slope
model_slope <- 'metacog_slope =~ atonyms_slope + metal_slope + number_slope

cog_g =~ atonyms_total_acc + metal_total_acc+ number_total_acc

metacog_slope ~~ cog_g'

fit_slope <- cfa(model_slope, data = df, missing = "fiml")

summary(fit_slope, fit.measures = TRUE, standardized = TRUE)

#Slope 1 factor with cognitive ability
model_slope_1f <- 'general_factor =~ atonyms_slope + metal_slope + number_slope + atonyms_total_acc + metal_total_acc+ number_total_acc
'

fit_slope_1f <- cfa(model_slope_1f, data = df, missing = "fiml")

summary(fit_slope_1f, fit.measures = TRUE, standardized = TRUE)

#korrelationer
library(psych)

slope_corr<- corr.test(df[,c("atonyms_slope", "metal_slope", "number_slope")])

slope_corr$r

slope_corr$p

conf_corr<- corr.test(df[,c("atonyms_conf", "metal_conf", "number_conf")])

conf_corr$r

conf_corr$p

abs_bias_corr<- corr.test(df[,c("atonyms_abs_bias", "metal_abs_bias", "number_abs_bias")])

abs_bias_corr$r

abs_bias_corr$p

total_acc_corr<- corr.test(df[,c("atonyms_total_acc", "metal_total_acc", "number_total_acc")])

total_acc_corr$r

total_acc_corr$p

###########################################
set.seed(12345)

#Variabler
conf_var <- c("atonyms_conf", "metal_conf", "number_conf")

acc_var <- c("atonyms_total_acc", "metal_total_acc", "number_total_acc")

mean_z_diff <- function(data){
  r_conf <- cor(
    data[, conf_var],
    use = "pairwise.complete.obs"
  )
  r_acc <- cor(
    data[,acc_var],
    use = "pairwise.complete.obs"
  )
  conf_r <- r_conf[lower.tri(r_conf)]
  acc_r <- r_acc[lower.tri(r_acc)]
  conf_z <- atanh(conf_r)
  acc_z <- atanh(acc_r)
  mean(conf_z) - mean(acc_z)
}

observed_diff <- mean_z_diff(df)

B <- 5000

boot_diff <- replicate(B, {
    rows <- sample(seq_len(nrow(df)), replace = TRUE)
    mean_z_diff(df[rows, ])
})

observed_diff

quantile(
  boot_diff,
  probs = c(.025, .975),
  na.rm = TRUE
)

mean(boot_diff > 0, na.rm = TRUE)

######################## Kolla huvudanalyser
plot(df$number_total_acc, df$number_gamma)

##################################Hur många har confidenskattningar på alla tre testen
Confidence_combinations <- as.data.frame(
  table(
    Atonym = df$atonym_has_conf_corrected,
    Metal = df$metal_has_conf_corrected,
    Numbers = df$number_has_conf_corrected,
    useNA = "ifany"
))

Confidence_combinations

has_number <- df$number_has_conf_corrected ==1

has_metal <- df$metal_has_conf_corrected ==1

has_atonym <- df$atonym_has_conf_corrected ==1

data.frame(
  Group = c(
    "Metal but not Numbers",
    "Antonym but not Numbers",
    "Metal AND Antonyms but not Numbers"
  ),
  N = c(
    sum(has_metal & !has_number, na.rm = TRUE),
    sum(has_atonym & !has_number, na.rm = TRUE),
    sum(has_metal & has_atonym & !has_number, na.rm = TRUE)
))

df$Cohort01 <- ifelse(df$Cohort == 1948, 0,1)

table(df$Cohort, df$Cohort01)

df$Sex01 <- ifelse(df$RSSEX == 1, 0,1)

table(df$RSSEX, df$Sex01)

library(lmerTest)

library(emmeans)

# ------------------------------------------------------------
# H3: Domain sensitivity
# Long-format data
# ------------------------------------------------------------
h3_long <- bind_rows(
  df %>% transmute(
    id,
    Sex01,
    Cohort01,
    Domain = "Antonyms",
    Bias = atonyms_bias,
    Absolute_bias = atonyms_abs_bias,
    Slope = atonyms_slope
  ),
  df %>% transmute(
    id,
    Sex01,
    Cohort01,
    Domain = "Metal Folding",
    Bias = metal_bias,
    Absolute_bias = metal_abs_bias,
    Slope = metal_slope
  ),
  df %>% transmute(
    id,
    Sex01,
    Cohort01,
    Domain = "Number Series",
    Bias = number_bias,
    Absolute_bias = number_abs_bias,
    Slope = number_slope
  )
)

# Ange ordningen på domänerna
h3_long$Domain <- factor(
  h3_long$Domain,
  levels = c("Antonyms", "Metal Folding", "Number Series")
)

# ------------------------------------------------------------
# H3a: Bias
# ------------------------------------------------------------
model_bias <- lmer(
  Bias ~ Domain + Sex01 + Cohort01 + (1 | id),
  data = h3_long,
  na.action = na.omit
)

anova(model_bias)

emmeans(
  model_bias,
  pairwise ~ Domain,
  adjust = "tukey"
)

# ------------------------------------------------------------
# H3b: Absolute bias
# ------------------------------------------------------------
model_absbias <- lmer(
  Absolute_bias ~ Domain + Sex01 + Cohort01 + (1 | id),
  data = h3_long,
  na.action = na.omit
)

anova(model_absbias)

emmeans(
  model_absbias,
  pairwise ~ Domain,
  adjust = "tukey"
)

# ------------------------------------------------------------
# H3c: Discrimination (Slope)
# ------------------------------------------------------------
model_slope <- lmer(
  Slope ~ Domain + Sex01 + Cohort01 + (1 | id),
  data = h3_long,
  na.action = na.omit
)

anova(model_slope)

emmeans(
  model_slope,
  pairwise ~ Domain,
  adjust = "tukey"
)

"Sex01" %in%(df)

df$Cohort01 <- ifelse(df$Cohort == 1948, 0,1)

table(df$Cohort, df$Cohort01)

df$Sex01 <- ifelse(df$RSSEX == 1, 0,1)

table(df$RSSEX, df$Sex01)

# ------------------------------------------------------------
# H3: Domain sensitivity
# Long-format data
# ------------------------------------------------------------
# Accuracy hämtas från respektive domäns tio confidence-bedömda uppgifter.
h3_long <- bind_rows(
  df %>% transmute(
    id,
    Sex01,
    Cohort01,
    Domain = "Antonyms",
    Accuracy = atonyms_acc,
    Bias = atonyms_bias,
    Absolute_bias = atonyms_abs_bias,
    Slope = atonyms_slope
  ),
  df %>% transmute(
    id,
    Sex01,
    Cohort01,
    Domain = "Metal Folding",
    Accuracy = metal_acc,
    Bias = metal_bias,
    Absolute_bias = metal_abs_bias,
    Slope = metal_slope
  ),
  df %>% transmute(
    id,
    Sex01,
    Cohort01,
    Domain = "Number Series",
    Accuracy = number_acc,
    Bias = number_bias,
    Absolute_bias = number_abs_bias,
    Slope = number_slope
  )
)

# Ange ordningen på domänerna
h3_long$Domain <- factor(
  h3_long$Domain,
  levels = c("Antonyms", "Metal Folding", "Number Series")
)

# ------------------------------------------------------------
# H3a: Bias
# ------------------------------------------------------------
model_bias_difficulty <- lmer(
  Bias ~ Domain + Accuracy + Sex01 + Cohort01 + (1 | id),
  data = h3_long,
  na.action = na.omit
)

anova(model_bias_difficulty)

emmeans(
  model_bias_difficulty,
  pairwise ~ Domain,
  adjust = "tukey"
)

# ------------------------------------------------------------
# H3b: Absolute bias
# ------------------------------------------------------------
model_absbias_difficulty <- lmer(
  Absolute_bias ~ Domain + Accuracy + Sex01 + Cohort01 + (1 | id),
  data = h3_long,
  na.action = na.omit
)

anova(model_absbias_difficulty)

emmeans(
  model_absbias_difficulty,
  pairwise ~ Domain,
  adjust = "tukey"
)

# ------------------------------------------------------------
# H3c: Discrimination (Slope)
# ------------------------------------------------------------
model_slope_difficulty <- lmer(
  Slope ~ Domain + Accuracy + Sex01 + Cohort01 + (1 | id),
  data = h3_long,
  na.action = na.omit
)

anova(model_slope_difficulty)

emmeans(
  model_slope_difficulty,
  pairwise ~ Domain,
  adjust = "tukey"
)

# ------------------------------------------------------------
# H5: Create long-format data and define Domain
# ------------------------------------------------------------
h5_long <- bind_rows(
  df %>% transmute(
    id,
    Sex01,
    Cohort01,
    Domain = "Antonyms",
    Absolute_bias = atonyms_abs_bias,
    Slope = atonyms_slope
  ),
  df %>% transmute(
    id,
    Sex01,
    Cohort01,
    Domain = "Metal Folding",
    Absolute_bias = metal_abs_bias,
    Slope = metal_slope
  ),
  df %>% transmute(
    id,
    Sex01,
    Cohort01,
    Domain = "Number Series",
    Absolute_bias = number_abs_bias,
    Slope = number_slope
  )
)

# Define Domain as a factor
h5_long$Domain <- factor(
  h5_long$Domain,
  levels = c("Antonyms", "Metal Folding", "Number Series")
)

# Check
table(h5_long$Domain, useNA = "ifany")

# ------------------------------------------------------------
# H5a: Absolute bias
# ------------------------------------------------------------
model_h5_absbias <- lmer(
  Absolute_bias ~ Domain * Cohort01 + Sex01 + (1 | id),
  data = h5_long,
  na.action = na.omit
)

anova(model_h5_absbias)

# Cohort differences within each domain
emmeans(
  model_h5_absbias,
  pairwise ~ Cohort01 | Domain,
  adjust = "holm"
)

# ------------------------------------------------------------
# H5b: Discrimination (Slope)
# ------------------------------------------------------------
model_h5_slope <- lmer(
  Slope ~ Domain * Cohort01 + Sex01 + (1 | id),
  data = h5_long,
  na.action = na.omit
)

anova(model_h5_slope)

# Cohort differences within each domain
emmeans(
  model_h5_slope,
  pairwise ~ Cohort01 | Domain,
  adjust = "holm"
)


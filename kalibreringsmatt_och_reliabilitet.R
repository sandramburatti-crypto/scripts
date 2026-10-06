# KALIBRERINGSMÅTT, ALFA OCH SPLIT-HALF — FRISTÅENDE R-FIL
# Bygger på de uppdaterade 01.clean_all.R och 04. analyses.R.
# Kan köras direkt på samma Finallyready.rds; inga andra skript eller paket krävs.
# Filen innehåller bara rengöring för måtten, måttberäkning och intern reliabilitet.
# Den skriver inte över rådata eller clean_all.rds.

# KÖRNING
# 1. Lägg filen där du vill och öppna den i R/RStudio.
# 2. Om rådata ligger på annan plats än /safe/data/Buratti, kör i R-konsolen:
#    options(buratti.data_dir = "C:/sökväg/till/dina/data")
# 3. Kör hela filen med Source eller source(..., encoding = "UTF-8").
#    Den läser Finallyready.rds och sparar resultaten i calibration_results/
#    under datamappen. Inställningen behöver anges igen i en ny R-session.
#    Resultatfiler med samma namn i den undermappen ersätts vid omkörning.

# RÅDATA OCH UPPGIFTSORDNING
# Rådata ska ha en rad per person och befintliga *_correct-kolumner med 0/1.
# Koden räknar inte ut facit. Den sätter *_correct till NA när svar saknas.
# Confidence förutsätts vara 0–100. Tomma/blankstegssvar behandlas som saknade.
# Dessa tio uppgifter paras i exakt den ordning som anges:
# Antonyms: mCat_{24,31,35,37,40,41,42,43,44,45}_correct <-> m_con_1:10.
# Metal: p_Q{10,12,16,19,20,21,24,26,27,28}_correct <-> p_con_1:10.
# Numbers: Tal{02,10,11,13,14,15,17,18,19,22}_correct <-> tcat_t{2,10,...,22}.
# Samtliga råsvar och rätt/fel-kolumner från de ursprungliga testen behövs:
# mCat_1:45, p_Q1:40 och Tal01a/Tal01b till Tal40a/Tal40b, med *_correct.
# Cohort används för reliabilitet per kohort; id följer med om kolumnen finns.

# MÅTTENS DEFINITIONER (samma som i de uppdaterade filerna)
# *_conf: medel av tillgänglig confidence på de tio uppgifterna, skala 0–100.
# *_conf_prop: samma medel omräknat till 0–1.
# *_acc: andel rätt bland besvarade av dessa tio uppgifter, skala 0–1.
# *_bias: *_conf_prop minus *_acc. Positivt = överconfidence, negativt = under.
#   Medelvärdena använder var för sig sina tillgängliga uppgifter, som i originalet.
# *_abs_bias: abs(*_bias), alltså absolutbeloppet EFTER medelvärdesbildningen.
#   Detta är inte medelvärdet av uppgiftsvisa absoluta differenser.
# *_slope: medelconfidence/100 vid rätt svar minus samma medel vid fel svar.
#   Kräver confidence vid minst ett rätt och minst ett fel svar; annars NA.
# *_gamma: Goodman–Kruskal gamma = (C-D)/(C+D), där C/D är samstämmiga/
#   motstridiga uppgiftspar med både confidence och rätt/fel. Lika värden
#   räknas inte. Inga jämförbara par ger NA. Intervallet är -1 till 1.
# *_total_acc: antal rätt i hela testet / 45 (Antonyms) eller / 40 (övriga).
#   Saknade uppgifter bidrar med 0. Ingen besvarad av de tio aktuella ger NA.
# *_n_answered och *_n_conf: antal besvarade respektive confidence-bedömda
#   uppgifter bland de tio. För Antonyms heter räknarna atonym_n_* (singular).

# 1. Läs rådata och kontrollera att nödvändiga kolumner finns
data_dir <- getOption("buratti.data_dir", "/safe/data/Buratti")
input_file <- file.path(data_dir, "Finallyready.rds")
output_dir <- file.path(data_dir, "calibration_results")
df <- readRDS(input_file)
if (!is.data.frame(df) || nrow(df) == 0L) stop("Rådata måste vara en data.frame med minst en person.")
required_raw_columns <- c(
  "Cohort", paste0("p_Q", 1:40), paste0("p_Q", 1:40, "_correct"),
  paste0("p_con_", 1:10), paste0("mCat_", 1:45),
  paste0("mCat_", 1:45, "_correct"), paste0("m_con_", 1:10),
  sprintf("Tal%02da", 1:40), sprintf("Tal%02db", 1:40),
  sprintf("Tal%02d_correct", 1:40),
  paste0("tcat_t", c(2,10,11,13,14,15,17,18,19,22))
)
missing_raw_columns <- setdiff(required_raw_columns, names(df))
if (length(missing_raw_columns)) {
  stop("Rådata saknar kolumner: ", paste(missing_raw_columns, collapse = ", "), call. = FALSE)
}
n_raw <- nrow(df)

# 2. Metal Folding: rengör svar, markera saknade rätt/fel och beräkna mått

df[paste0("p_Q", 1:40)] <- lapply(df[paste0("p_Q", 1:40)], function (x) replace(x, trimws(x) == "", NA))

for (i in 1:40) df[[paste0("p_Q", i, "_correct")]][is.na(df[[paste0("p_Q", i )]])] <- NA

df[paste0("p_con_", 1:10)] <- lapply(df[paste0("p_con_", 1:10)], function (x) replace(x, trimws(x) == "", NA))

df[paste0("p_con_", 1:10)] <- lapply(df[paste0("p_con_", 1:10)], as.numeric)

df$metal_n_conf <- rowSums(!is.na(df[paste0("p_con_", 1:10)]))

df$metal_conf <- rowMeans(df[paste0("p_con_", 1:10)], na.rm = TRUE)

items <- c(10, 12, 16, 19, 20, 21, 24, 26, 27, 28)

df$metal_acc <- rowMeans(df[paste0("p_Q", items, "_correct")], na.rm = TRUE)

metal_items <- paste0("p_Q", c(10, 12, 16, 19, 20, 21, 24, 26, 27, 28), "_correct")

df$metal_bias <- (df$metal_conf / 100) - df$metal_acc

df$metal_conf_prop <- df$metal_conf / 100

df$metal_abs_bias <- abs(df$metal_conf_prop - df$metal_acc)

df$metal_abs_bias <- as.numeric(df$metal_abs_bias)

conf <- as.matrix(df[, paste0("p_con_", 1:10)]) / 100

corr <- as.matrix(df[, paste0("p_Q", items, "_correct")])

df$metal_slope <- sapply(1:nrow(df), function(i) {
    mean_correct <- mean(conf[i, corr[i, ] == 1], na.rm = TRUE)
    mean_incorrect <- mean(conf[i, corr[i, ] == 0], na.rm = TRUE)
    mean_correct - mean_incorrect
})

metal_gamma <- function(x, y) {
  stopifnot(length(x) == length(y))
  ok <- !is.na(x) & !is.na(y)
  x <- x[ok]
  y <- y[ok]
  n <- length(x)
  if (n < 2L) return(NA_real_)

  C <- D <- 0L
  for (i in seq_len(n - 1L)) {
    for (j in seq.int(i + 1L, n)) {
      product <- (x[i] - x[j]) * (y[i] - y[j])
      if (product > 0) C <- C + 1L
      if (product < 0) D <- D + 1L
    }
  }

  if (C + D == 0L) return(NA_real_)
  (C - D) / (C + D)
}

df$metal_gamma <- sapply(1:nrow(df), function(i) metal_gamma(conf[i, ], corr [i, ]))

df$metal_n_answered <- rowSums(!is.na(df[metal_items]))

df$metal_total_acc <- rowSums(df [, paste0("p_Q", 1:40, "_correct")], na.rm =TRUE)/40

df$metal_total_acc[df$metal_n_answered ==0]<- NA

items_number <- c(2, 10, 11, 13, 14, 15, 17, 18, 19, 22)

# 3. Number Series: samma princip, båda delsvaren krävs för ett besvarat item

conf_vars <- c("tcat_t2", "tcat_t10", "tcat_t11",
  "tcat_t13", "tcat_t14", "tcat_t15",
  "tcat_t17", "tcat_t18", "tcat_t19", "tcat_t22")

number_items <- paste0("Tal",
  sprintf("%02d", c(2,10,11,13,14,15,17,18,19,22)),
  "_correct")

tal_items <- sprintf("Tal%02d", 1:40)

tal_response_vars <- as.vector(rbind(
    paste0(tal_items, "a"),
    paste0(tal_items, "b")
))

df[tal_response_vars] <- lapply(df[tal_response_vars], function(x) {
    replace(x, trimws(x) == "", NA)
})

for (i in 1:40) {
  item <- sprintf("Tal%02d", i)
  a_var <- paste0(item, "a")
  b_var <- paste0(item, "b")
  correct_var <- paste0(item, "_correct")
  df[[correct_var]][is.na(df[[a_var]])|is.na(df[[b_var]])] <- NA
}

df[conf_vars] <- lapply(df[conf_vars], function(x) {
    x <- trimws(x)
    x[x== ""] <- NA
    as.numeric(x)
})

df$number_n_answered <- rowSums(!is.na(df[number_items]))

df$number_n_conf <- rowSums(!is.na(df[conf_vars]))

df$number_conf_answered <- df$number_n_conf

df$number_conf <- rowMeans(df[conf_vars], na.rm = TRUE)

df$number_conf[is.nan(df$number_conf)] <- NA

df$number_acc <- rowMeans(df[paste0("Tal", sprintf("%02d", items_number), "_correct")], na.rm = TRUE)

df$number_acc[is.nan(df$number_acc)] <- NA

df$number_bias <- (df$number_conf / 100) - df$number_acc

df$number_conf_prop <- df$number_conf / 100

df$number_abs_bias <- abs(df$number_conf_prop - df$number_acc)

df$number_abs_bias <- as.numeric(df$number_abs_bias)

conf <- as.matrix(df[conf_vars]) / 100

corr <- as.matrix(df[paste0("Tal", sprintf("%02d", items_number), "_correct")])

df$number_slope <- sapply (1:nrow(df), function(i)
  mean(conf[i, corr[i, ] == 1], na.rm = TRUE) -
  mean(conf[i, corr[i, ] == 0], na.rm = TRUE))

number_gamma <- function(x, y) {
  stopifnot(length(x) == length(y))
  ok <- !is.na(x) & !is.na(y)
  x <- x[ok]
  y <- y[ok]
  n <- length(x)
  if (n < 2L) return(NA_real_)

  C <- D <- 0L
  for (i in seq_len(n - 1L)) {
    for (j in seq.int(i + 1L, n)) {
      product <- (x[i] - x[j]) * (y[i] - y[j])
      if (product > 0) C <- C + 1L
      if (product < 0) D <- D + 1L
    }
  }

  if (C + D == 0L) return(NA_real_)
  (C - D) / (C + D)
}

df$number_gamma <- sapply(1:nrow(df), function(i) number_gamma(conf[i, ], corr [i, ]))

df$number_n_answered <- rowSums(!is.na(df[number_items]))

df$number_total_acc <- rowSums(df[, paste0("Tal",sprintf("%02d", 1:40), "_correct")], na.rm =TRUE)/40

df$number_total_acc[df$number_n_answered == 0] <- NA

# 4. Antonyms: samma princip; uppgifterna väljs med kolumnnamn

df[paste0("mCat_", 1:45)] <- lapply(df[paste0("mCat_", 1:45)], function (x) replace(x, trimws(x) == "", NA))

for (i in 1:45) df[[paste0("mCat_", i, "_correct")]][is.na(df[[paste0("mCat_", i )]])] <- NA

df[paste0("m_con_", 1:10)] <- lapply(df[paste0("m_con_", 1:10)], function (x) replace(x, trimws(x) == "", NA))

df[paste0("m_con_", 1:10)] <- lapply(df[paste0("m_con_", 1:10)], as.numeric)

df$atonym_n_conf <- rowSums(!is.na(df[paste0("m_con_", 1:10)]))

df$atonyms_conf <- rowMeans(df[paste0("m_con_", 1:10)], na.rm = TRUE)

items <- c(24, 31, 35, 37, 40, 41, 42, 43, 44, 45)

df$atonyms_acc <- rowMeans(df[paste0("mCat_", items, "_correct")], na.rm = TRUE)

atonym_items <- c(24, 31, 35, 37, 40, 41, 42, 43, 44, 45)

df$atonyms_bias <- (df$atonyms_conf / 100) - df$atonyms_acc

df$atonyms_conf_prop <- df$atonyms_conf / 100

df$atonyms_abs_bias <- abs(df$atonyms_conf_prop - df$atonyms_acc)

df$atonyms_abs_bias <- as.numeric(df$atonyms_abs_bias)

conf <- as.matrix(df[, paste0("m_con_", 1:10)]) / 100

corr <- as.matrix(df[, paste0("mCat_", items, "_correct")])

df$atonyms_slope <- sapply(1:nrow(df), function(i) {
    mean_correct <- mean(conf[i, corr[i, ] == 1], na.rm = TRUE)
    mean_incorrect <- mean(conf[i, corr[i, ] == 0], na.rm = TRUE)
    mean_correct - mean_incorrect
})

atonyms_gamma <- function(x, y) {
  stopifnot(length(x) == length(y))
  ok <- !is.na(x) & !is.na(y)
  x <- x[ok]
  y <- y[ok]
  n <- length(x)
  if (n < 2L) return(NA_real_)

  C <- D <- 0L
  for (i in seq_len(n - 1L)) {
    for (j in seq.int(i + 1L, n)) {
      product <- (x[i] - x[j]) * (y[i] - y[j])
      if (product > 0) C <- C + 1L
      if (product < 0) D <- D + 1L
    }
  }

  if (C + D == 0L) return(NA_real_)
  (C - D) / (C + D)
}

df$atonyms_gamma <- sapply(1:nrow(df), function(i) atonyms_gamma(conf[i, ], corr [i, ]))

df$atonym_n_answered <- rowSums(
  !is.na(df[paste0("mCat_", atonym_items, "_correct")])
)

df$atonyms_total_acc <- rowSums(df [,paste0("mCat_", 1:45, "_correct")], na.rm =TRUE)/45

df$atonyms_total_acc[df$atonym_n_answered == 0] <- NA

# 5. Samma deltagarfilter som i 01.clean_all.R:
# behåll personer med minst ett besvarat rätt/fel-item bland de tio utvalda
# uppgifterna i minst en domän. Antal confidence-bedömningar räknas separat.

df <- df[
  df$number_n_answered >= 1 |
  df$metal_n_answered >= 1 |
  df$atonym_n_answered >= 1,
]

# 6. Saknade/odefinierade beräknade mått sparas konsekvent som NA

derived_measure_vars <- unlist(lapply(c("atonyms", "metal", "number"), function(domain) {
  paste0(domain, c("_conf", "_conf_prop", "_acc", "_bias", "_abs_bias",
                  "_slope", "_gamma", "_total_acc"))
}), use.names = FALSE)

df[derived_measure_vars] <- lapply(df[derived_measure_vars], function(x) {
  x[is.nan(x)] <- NA_real_
  x
})

# ------------------------------------------------------------
# 7. Intern konsistens och split-half-konsistens INOM varje domän
# ------------------------------------------------------------
# Intern reliabilitet undersöks med alfa för accuracy/confidence och
# upprepad split-half för bias, absolut bias, slope och gamma.
# Inga befintliga mått skrivs över. Alfa beräknas inte för kalibreringsmåtten.
# Alfa skattar intern konsistens för total accuracy, tiouppgifts-accuracy och
# medelconfidence. För total accuracy räknas saknade uppgifter som 0, precis
# som i totalpoängen med fast nämnare; bara personer med tillgänglig
# totalpoäng ingår. För tiouppgifts-accuracy och confidence används kompletta
# fall på respektive tio kolumner. Därför kan urvalen skilja sig åt.
# Alfa är inte ett test av endimensionalitet eller stabilitet över tid.
#
# För bias, absolut bias, slope och gamma delas samma tio uppgifter i 5+5.
# Bias i varje halva använder separata tillgängliga medelvärden för confidence
# och accuracy, precis som i måtten ovan.
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

# 8. Beräkna reliabilitet och spara mått och tabeller
if (nrow(df) == 0L) stop("Ingen person återstår efter deltagarfiltret; inga resultat sparas.")
within_reliability <- check_within_test_reliability(df)
measure_names <- c(
  derived_measure_vars, "atonym_n_answered", "atonym_n_conf",
  "metal_n_answered", "metal_n_conf", "number_n_answered", "number_n_conf",
  "number_conf_answered"
)
# Utdata innehåller id/Cohort och beräknade mått; råsvar exporteras inte.
calibration_measures <- df[c(intersect(c("id", "Cohort"), names(df)), measure_names)]
calibration_results <- list(
  measures = calibration_measures,
  alpha = within_reliability$alpha,
  split_summary = within_reliability$split_summary,
  split_details = within_reliability$split_details,
  partitions = within_reliability$partitions,
  n_raw = n_raw, n_retained = nrow(df)
)
if (!dir.exists(output_dir) && !dir.create(output_dir, recursive = TRUE)) {
  stop("Kunde inte skapa resultatmappen: ", output_dir)
}
saveRDS(calibration_results, file.path(output_dir, "calibration_results.rds"))
write.csv(calibration_measures, file.path(output_dir, "calibration_measures.csv"),
          row.names = FALSE, na = "", fileEncoding = "UTF-8")
write.csv(within_reliability$alpha, file.path(output_dir, "alpha.csv"),
          row.names = FALSE, na = "", fileEncoding = "UTF-8")
write.csv(within_reliability$split_summary, file.path(output_dir, "split_half_summary.csv"),
          row.names = FALSE, na = "", fileEncoding = "UTF-8")
write.csv(within_reliability$split_details, file.path(output_dir, "split_half_details.csv"),
          row.names = FALSE, na = "", fileEncoding = "UTF-8")
cat("\nPersoner i rådata:", n_raw, "; efter samma deltagarfilter som originalet:", nrow(df), "\n")
cat("\nAlfa: accuracy och confidence inom varje domän\n")
print(within_reliability$alpha, row.names = FALSE)
cat("\nSplit-half-konsistens: ICC mellan femuppgiftshalvor\n")
cat("Median/percentiler avser uppdelningarna, inte fulltestreliabilitet.\n")
print(within_reliability$split_summary, row.names = FALSE)
cat("\nResultat sparade i:", normalizePath(output_dir, winslash = "/"), "\n")
# Efter körningen finns calibration_measures och within_reliability även i R.
# Samlat sparat resultat kan läsas med readRDS(".../calibration_results.rds").

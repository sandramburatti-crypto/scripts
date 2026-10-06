# 01.clean_all.R — försiktigt städad version
# Gamma, antonymers kolumnval och kopplingen till analysfilen är rättade.
# Kända problem är markerade med GRANSKA och beskrivna i README.txt.
# Båda filerna använder samma datamapp. Anpassa via options() vid behov:
# options(buratti.data_dir = "C:/sökväg/till/dina/data")
# Kör först denna fil och därefter 04. analyses.R.

#### do Finallyready into metal file
data_dir <- getOption("buratti.data_dir", "/safe/data/Buratti")
df <- readRDS(file.path(data_dir, "Finallyready.rds"))

##### Fixing number
##### Setting Missingness to NA in p_QX_correct
df[paste0("p_Q", 1:40)] <- lapply(df[paste0("p_Q", 1:40)], function (x) replace(x, trimws(x) == "", NA))

for (i in 1:40) df[[paste0("p_Q", i, "_correct")]][is.na(df[[paste0("p_Q", i )]])] <- NA

###### check if correct
df[1:10, c("p_Q20", "p_Q20_correct")]

###### Setting missingness to NA in P_con 
df[paste0("p_con_", 1:10)] <- lapply(df[paste0("p_con_", 1:10)], function (x) replace(x, trimws(x) == "", NA))

df[paste0("p_con_", 1:10)] <- lapply(df[paste0("p_con_", 1:10)], as.numeric)

# Räkna confidence efter att tomma strängar har satts till NA.
df$metal_n_conf <- rowSums(!is.na(df[paste0("p_con_", 1:10)]))

#########Compute mean metal confidence 
df$metal_conf <- rowMeans(df[paste0("p_con_", 1:10)], na.rm = TRUE)

######Check
df[1:10, c(paste0("p_con_", 1:10), "metal_conf")]

########  calculcating accuracy for confidence judge items. 
items <- c(10, 12, 16, 19, 20, 21, 24, 26, 27, 28)

df$metal_acc <- rowMeans(df[paste0("p_Q", items, "_correct")], na.rm = TRUE)

metal_items <- paste0("p_Q", c(10, 12, 16, 19, 20, 21, 24, 26, 27, 28), "_correct")

######Check
df[1:10, c(paste0("p_Q", items, "_correct"), "metal_acc")]

df_valid <- df[!is.na(df$metal_acc),]

nrow(df_valid)

###### Summary and histogram. 
summary (df$metal_acc)

hist (df_valid$metal_acc)

######### Create number bias, absolut bias and slope
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

mean(df$metal_slope, na.rm =TRUE)

summary(df$metal_slope)

hist(df$metal_slope)

###### Create gamma
# Goodman–Kruskal gamma: jämför alla par med både confidence och rätt/fel.
# Lika värden räknas inte. Utan jämförbara par returneras NA.
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

summary(df$metal_gamma)

hist(df$metal_gamma)

#######Create cognitive ability for atonyms
df$metal_n_answered <- rowSums(!is.na(df[metal_items]))

df$metal_total_acc <- rowSums(df [, paste0("p_Q", 1:40, "_correct")], na.rm =TRUE)/40

df$metal_total_acc[df$metal_n_answered ==0]<- NA

summary(df$metal_total_acc)

par(mar = c(4,4,2,1))

hist(df$metal_total_acc)

summary(df$metal_abs_bias)

par(mar = c(4,4,2,1))

hist(df$metal_abs_bias)

##### do Finallyready into NUMBERS
##### Define items that were confidence judged
{ conf_vars <- c("tcat_t2", "tcat_t10", "tcat_t11",
    "tcat_t13", "tcat_t14", "tcat_t15",
    "tcat_t17", "tcat_t18", "tcat_t19", "tcat_t22" )
}

items_number <- c(2, 10, 11, 13, 14, 15, 17, 18, 19, 22)

conf_vars <- c("tcat_t2", "tcat_t10", "tcat_t11",
  "tcat_t13", "tcat_t14", "tcat_t15",
  "tcat_t17", "tcat_t18", "tcat_t19", "tcat_t22")

number_items <- paste0("Tal",
  sprintf("%02d", c(2,10,11,13,14,15,17,18,19,22)),
  "_correct")

##### Clean response variables Tal01a to Tal40b
tal_items <- sprintf("Tal%02d", 1:40)

tal_response_vars <- as.vector(rbind(
    paste0(tal_items, "a"),
    paste0(tal_items, "b")
))

df[tal_response_vars] <- lapply(df[tal_response_vars], function(x) {
    replace(x, trimws(x) == "", NA)
})

######### Set Talxx_correct to NA when either a or b responses are missing. 
for (i in 1:40) {
  item <- sprintf("Tal%02d", i)
  a_var <- paste0(item, "a")
  b_var <- paste0(item, "b")
  correct_var <- paste0(item, "_correct")
  df[[correct_var]][is.na(df[[a_var]])|is.na(df[[b_var]])] <- NA
}

###### check if correct
df[20:40, c("Tal20a", "Tal20b", "Tal20_correct")]

###### Clean confidence variables for tcat02 tp tcat_22
df[conf_vars] <- lapply(df[conf_vars], function(x) {
    x <- trimws(x)
    x[x== ""] <- NA
    as.numeric(x)
})

# Båda räknarna beräknas efter rengöring av svar och confidence.
df$number_n_answered <- rowSums(!is.na(df[number_items]))
df$number_n_conf <- rowSums(!is.na(df[conf_vars]))
# Behåll det tidigare namnet som synonym för samma rengjorda antal.
df$number_conf_answered <- df$number_n_conf

cat("Efter rengöring av Number Series (före deltagarfiltret)\n")
cat("Rows:", nrow(df), "\n")
print(table(df$Cohort, useNA = "ifany"))
print(table(df$number_conf_answered, useNA = "ifany"))
print(table(answered = df$number_n_answered >= 1,
            confidence = df$number_n_conf >= 1))

######### Create number conficence mean
df$number_conf <- rowMeans(df[conf_vars], na.rm = TRUE)

df$number_conf[is.nan(df$number_conf)] <- NA

######Check
df[1:10, c(conf_vars, "number_conf")]

########  Create number accuracy  proportion correct of answered e.g. 1/5= 0.20
df$number_acc <- rowMeans(df[paste0("Tal", sprintf("%02d", items_number), "_correct")], na.rm = TRUE)

df$number_acc[is.nan(df$number_acc)] <- NA

######Check
df[1:10, c(paste0("Tal", sprintf("%02d", items_number), "_correct"), "number_acc")]

######### Create number bias, absolut bias and slope
df$number_bias <- (df$number_conf / 100) - df$number_acc

df$number_conf_prop <- df$number_conf / 100

df$number_abs_bias <- abs(df$number_conf_prop - df$number_acc)

df$number_abs_bias <- as.numeric(df$number_abs_bias)

conf <- as.matrix(df[conf_vars]) / 100

corr <- as.matrix(df[paste0("Tal", sprintf("%02d", items_number), "_correct")])

df$number_slope <- sapply (1:nrow(df), function(i)
  mean(conf[i, corr[i, ] == 1], na.rm = TRUE) -
  mean(conf[i, corr[i, ] == 0], na.rm = TRUE))

mean(df$number_slope, na.rm =TRUE)

###### Create gamma
# Goodman–Kruskal gamma: jämför alla par med både confidence och rätt/fel.
# Lika värden räknas inte. Utan jämförbara par returneras NA.
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

summary(df$number_gamma)

hist(df$number_gamma)

#######Create cognitive ability for number
df$number_n_answered <- rowSums(!is.na(df[number_items]))

number_mat <- as.matrix(df[paste0("Tal", sprintf("%02d", items_number), "_correct")])

number_mat <- number_mat * 1

df$number_total_acc<- rowSums(number_mat, na.rm = TRUE) / 40

df$number_total_acc <- rowSums(df[, paste0("Tal",sprintf("%02d", 1:40), "_correct")], na.rm =TRUE)/40

df$number_total_acc[df$number_n_answered == 0] <- NA

summary(df$number_total_acc)

par(mar = c(4,4,2,1))

hist(df$number_total_acc)

##### do Finallyready into confidencefile
##### Fixing atonyms
##### Setting Missingness to NA in p_QX_correct
df[paste0("mCat_", 1:45)] <- lapply(df[paste0("mCat_", 1:45)], function (x) replace(x, trimws(x) == "", NA))

for (i in 1:45) df[[paste0("mCat_", i, "_correct")]][is.na(df[[paste0("mCat_", i )]])] <- NA

###### check if correct
df[1:50, c("mCat_43", "mCat_43_correct")]

###### Setting missingness to NA in m_con_1 
df[paste0("m_con_", 1:10)] <- lapply(df[paste0("m_con_", 1:10)], function (x) replace(x, trimws(x) == "", NA))

df[paste0("m_con_", 1:10)] <- lapply(df[paste0("m_con_", 1:10)], as.numeric)

df$atonym_n_conf <- rowSums(!is.na(df[paste0("m_con_", 1:10)]))

#########Compute mean confidence for atonyms
df$atonyms_conf <- rowMeans(df[paste0("m_con_", 1:10)], na.rm = TRUE)

######Check
df[1:10, c(paste0("p_con_", 1:10), "atonyms_conf")]

########  Create number accuracy  proportion correct of answered e.g. 1/5= 0.20
items <- c(24, 31, 35, 37, 40, 41, 42, 43, 44, 45)

df$atonyms_acc <- rowMeans(df[paste0("mCat_", items, "_correct")], na.rm = TRUE)

atonym_items <- c(24, 31, 35, 37, 40, 41, 42, 43, 44, 45)

######Check
df[1:10, c(paste0("mCat_", items, "_correct"), "atonyms_acc")]

df_valid <- df[!is.na(df$atonyms_acc),]

nrow(df_valid)

######### Create number bias, absolut bias and slope
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

##################################
c(
  N_bias = sum(!is.na(df$atonyms_bias)),
  mean_bias = mean(df$atonyms_bias, na.rm = TRUE),
  sd_bias = sd(df$atonyms_bias, na.rm = TRUE),
  N_abs_bias = sum(!is.na(df$atonyms_abs_bias)),
  mean_abs_bias = mean(df$atonyms_abs_bias, na.rm = TRUE),
  sd_abs_bias = sd(df$atonyms_abs_bias, na.rm = TRUE)
)

c(
  N_conf = sum(!is.na(df$metal_conf)),
  mean_conf = mean(df$metal_conf, na.rm = TRUE),
  sd_conf = sd(df$metal_conf, na.rm = TRUE))

mean(df$atonyms_slope, na.rm =TRUE)

summary(df$atonyms_slope)

hist(df$atonyms_slope)

###### Create gamma
# Goodman–Kruskal gamma: jämför alla par med både confidence och rätt/fel.
# Lika värden räknas inte. Utan jämförbara par returneras NA.
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

summary(df$atonyms_gamma)

hist(df$atonyms_gamma)

#######Create cognitive ability for atonyms
# Räkna besvarade uppgifter bland de tio med confidence-bedömningar.
# Välj rätt/fel-kolumner med namn så att kolumnordningen inte påverkar antalet.
df$atonym_n_answered <- rowSums(
  !is.na(df[paste0("mCat_", atonym_items, "_correct")])
)

df$atonyms_total_acc <- rowSums(df [,paste0("mCat_", 1:45, "_correct")], na.rm =TRUE)/45

df$atonyms_total_acc[df$atonym_n_answered == 0] <- NA

summary(df$atonyms_total_acc)

par(mar = c(4,4,2,1))

hist(df$atonyms_total_acc)

summary(df$number_abs_bias)

par(mar = c(4,4,2,1))

hist(df$atonyms_abs_bias)

# Kodningarna är samma som i analysfilens original.
# Cohort: 1948 = 0, andra icke-missing värden = 1.
# RSSEX: 1 = 0, andra icke-missing värden = 1.
df$Cohort01 <- ifelse(df$Cohort == 1948, 0, 1)
df$Sex01 <- ifelse(df$RSSEX == 1, 0, 1)

########lägg på filter
df <- df[
  df$number_n_answered >= 1 |
  df$metal_n_answered >= 1 |
  df$atonym_n_answered >= 1,
]

# check
nrow(df)

table(df$Cohort, useNA = "ifany")

# ------------------------------------------------------------
# Använd NA konsekvent för saknade eller odefinierade härledda mått.
# rowMeans()/mean() kan ge NaN om inga observationer finns kvar efter
# borttagning av saknade värden, t.ex. slope utan korrekta/felaktiga svar.
# Omvandlingen ändrar endast NaN; giltiga värden och befintliga NA bevaras.
# Rådata och deltagarurval ändras inte. Inf/-Inf omvandlas inte här.
# ------------------------------------------------------------
derived_measure_vars <- unlist(lapply(c("atonyms", "metal", "number"), function(domain) {
  paste0(domain, c("_conf", "_conf_prop", "_acc", "_bias", "_abs_bias",
                  "_slope", "_gamma", "_total_acc"))
}), use.names = FALSE)

df[derived_measure_vars] <- lapply(df[derived_measure_vars], function(x) {
  x[is.nan(x)] <- NA_real_
  x
})

######### SAVE cleaned dataset
saveRDS(df, file.path(data_dir, "super_clean_all.rds"))

                           capture.output(
  source(file.choose(), echo = TRUE, encoding = "UTF-8"),
  file = "output_01_cleaned.txt"

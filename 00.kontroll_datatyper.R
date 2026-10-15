# KONTROLL AV DATAKOPPLING OCH DATATYPER, version 2026-10-15
# Kör hela filen med Source. Inga extra paket behövs.
# Läser samma clean_all.rds som 04.manus.R och sparar en kort textrapport.
# Inga analyser görs och clean_all.rds ändras inte.
# Vid behov, ange datamappen i konsolen före Source:
# options(buratti.data_dir = "C:/sökväg/till/dina/data")

message("KONTROLLVERSION 2026-10-15")
kontroll_data_dir <- getOption("buratti.data_dir", "/safe/data/Buratti")
kontroll_data_file <- file.path(kontroll_data_dir, "clean_all.rds")
if (!file.exists(kontroll_data_file)) {
  stop("clean_all.rds saknas i ", kontroll_data_dir, ". Ange samma datamapp som i 01.")
}
kontroll_data <- readRDS(kontroll_data_file)

# Exakt de kolumner som 04.manus använder i sin numeriska kontroll.
kontroll_measures <- c("acc", "conf", "bias", "abs_bias", "slope", "total_acc")
kontroll_correct <- c(
  paste0("mCat_", 1:45, "_correct"),
  paste0("p_Q", 1:40, "_correct"),
  sprintf("Tal%02d_correct", 1:40)
)
kontroll_required <- c(
  "Sex01", "Cohort01", "atonym_n_answered", "atonym_n_conf",
  "metal_n_answered", "metal_n_conf", "number_n_answered", "number_n_conf",
  paste0("atonyms_", kontroll_measures),
  paste0("metal_", kontroll_measures),
  paste0("number_", kontroll_measures),
  paste0("m_con_", 1:10), paste0("p_con_", 1:10),
  paste0("tcat_t", c(2, 10, 11, 13, 14, 15, 17, 18, 19, 22)),
  kontroll_correct
)
kontroll_missing <- setdiff(kontroll_required, names(kontroll_data))
kontroll_present <- intersect(kontroll_required, names(kontroll_data))
kontroll_bad <- kontroll_present[!vapply(kontroll_data[kontroll_present], is.numeric, logical(1))]

# Kolumnnamn och datatyp, utan några enskilda deltagares svar.
kontroll_table <- data.frame(
  Column = kontroll_bad,
  Type = vapply(kontroll_data[kontroll_bad], function(x) paste(class(x), collapse = "/"), character(1)),
  n_nonmissing = vapply(kontroll_data[kontroll_bad], function(x) sum(!is.na(x)), integer(1)),
  logical_correct = kontroll_bad %in% kontroll_correct &
    vapply(kontroll_data[kontroll_bad], is.logical, logical(1)),
  row.names = NULL
)
kontroll_report <- capture.output({
  cat("DATATYPKONTROLL 2026-10-15\n")
  cat("Datafil: ", normalizePath(kontroll_data_file, winslash = "/"), "\n", sep = "")
  cat("Datafil senast sparad: ", format(file.info(kontroll_data_file)$mtime), "\n", sep = "")
  cat("Antal deltagarrader: ", nrow(kontroll_data), "\n", sep = "")
  cat("Antal kolumner som kontrolleras: ", length(kontroll_required), "\n\n", sep = "")

  if (length(kontroll_missing) > 0) {
    cat("KOLUMNER SOM SAKNAS:\n")
    cat(paste(kontroll_missing, collapse = "\n"), "\n\n", sep = "")
  }
  if (nrow(kontroll_table) == 0) {
    cat("Alla befintliga kolumner i denna kontroll är numeriska.\n")
  } else {
    cat("KOLUMNER SOM STOPPAR DEN ÄLDRE MANUSFILENS DATATYPKONTROLL:\n")
    print(kontroll_table, row.names = FALSE, max = 10000)
    cat("\nType anger R:s datatyp. logical betyder TRUE/FALSE/NA.\n")
    cat("logical_correct=TRUE anger rätt/fel som den nya manusfilen läser som 1/0/NA.\n")
  }
  kontroll_remaining <- kontroll_table[!kontroll_table$logical_correct, , drop = FALSE]
  if (length(kontroll_missing) == 0 && nrow(kontroll_remaining) == 0) {
    cat("\nDen uppdaterade manusfilen kan passera just denna datatypkontroll.\n")
  } else if (nrow(kontroll_remaining) > 0) {
    cat("\nDESSA DATATYPER BEHÖVER FORTFARANDE UTREDAS:\n")
    print(kontroll_remaining, row.names = FALSE, max = 10000)
  }
  cat("\nKontrollen läser datafilen på sökvägen ovan.\n")
  cat("Den kör inga statistiska modeller och kontrollerar inte alla möjliga datafel.\n")
})

cat(paste(kontroll_report, collapse = "\n"), "\n", sep = "")
kontroll_report_file <- file.path(kontroll_data_dir, "kontroll_datatyper_2026-10-15.txt")
writeLines(kontroll_report, kontroll_report_file, useBytes = TRUE)
cat("\nRapport sparad: ", normalizePath(kontroll_report_file, winslash = "/"), "\n", sep = "")

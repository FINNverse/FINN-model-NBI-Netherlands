#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## FINN-model-NBI-Netherlands | calibration/02_extract_nbi_tables.R
##
## WHAT   Exports the tables the preparation needs from the NBI-7 public Access
##        database (01_download_nbi_raw.R) to CSV: the tree, plot-definition and
##        plot-measurement tables of the three cycles (MFV, NBI-6, NBI-7), the
##        two lists of permanent plots remeasured between cycles, and the
##        lookup tables for species, tree class, harvest status and plot
##        appearance.
## WHY    The .accdb can only be read through mdbtools; doing it once here keeps
##        every later step in plain data.table on CSV, and the row counts written
##        to extract_log.csv let a rerun be checked against this vintage.
## NEEDS  mdbtools (`mdb-tables`, `mdb-export`) on the PATH (Homebrew / apt: mdbtools).
## INPUT  <DATA_DIR>/raw/probos/unz/NBI7_openbare database_definitief_UpdateJuli2024.accdb
## OUTPUT <DATA_DIR>/raw/csv/<table>.csv (one per table), <DATA_DIR>/raw/csv/extract_log.csv
## RUN    Rscript calibration/02_extract_nbi_tables.R   (from the repository root, ~1 min)
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
suppressPackageStartupMessages(library(data.table))

DATA_DIR <- file.path("..", "..", "FINN-zoo-data", "nbi")

PARAMS <- list(
  database = file.path(DATA_DIR, "raw", "probos", "unz", "NBI7_openbare database_definitief_UpdateJuli2024.accdb"),
  tables = c(
    ## measurements: one row per tree (incl. removed / dead trees reported back on permanent plots)
    "data_MFV_boomgegevens", "data_NBI6_boomgegevens", "data_NBI7_boomgegevens",
    ## plot definition (municipality, province, permanent flag, 1 km INSPIRE cell)
    "data_MFV_plotdefinitie", "data_NBI6_plotdefinitie", "data_NBI7_plotdefinitie",
    ## plot visit (date, radius, forest area of the circle, appearance / stand type, stand variables)
    "data_MFV_plotmetingen", "data_NBI6_plotmetingen", "data_NBI7_plotmetingen",
    ## the plots whose later measurement is a true remeasurement of the earlier one
    "data_permanent_MFV_NBI6", "data_permanent_NBI6_NBI7",
    ## code lists
    "lookup_MFV_Boomsoort", "lookup_NBI6_Boomsoort", "lookup_NBI7_Boomsoort", "lookup_Boomklasse",
    "lookup_Oogststatus", "lookup_NBI7_Verschijningsvorm", "lookup_Positie", "lookup_Deelgebied2Bos"),
  date_format = "%Y-%m-%d"     # mdb-export's -D: ISO dates instead of the US default
)

if (Sys.which("mdb-export") == "") stop("mdbtools not found on the PATH (brew install mdbtools / apt install mdbtools)")
if (!file.exists(PARAMS$database)) stop("database not found: ", PARAMS$database, " (run 01_download_nbi_raw.R first)")
out_dir <- file.path(DATA_DIR, "raw", "csv"); dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

available <- system2("mdb-tables", c("-1", shQuote(PARAMS$database)), stdout = TRUE)
missing <- setdiff(PARAMS$tables, available)
if (length(missing)) stop("tables not in the database: ", paste(missing, collapse = ", "))

log <- rbindlist(lapply(PARAMS$tables, function(tb) {
  out <- file.path(out_dir, paste0(tb, ".csv"))
  ## -b strip drops binary junk; quotes keep text fields with embedded commas intact
  status <- system2("mdb-export", c("-D", shQuote(PARAMS$date_format), "-b", "strip", shQuote(PARAMS$database), shQuote(tb)), stdout = out)
  if (status != 0L) stop("mdb-export failed for ", tb)
  dt <- fread(out, encoding = "Latin-1", showProgress = FALSE)
  message(sprintf("%-32s %7d rows, %3d columns", tb, nrow(dt), ncol(dt)))
  data.table(table = tb, rows = nrow(dt), columns = ncol(dt), file = basename(out))
}))
log[, database := basename(PARAMS$database)]
fwrite(log, file.path(out_dir, "extract_log.csv"))
message("done: ", out_dir, " (", nrow(log), " tables)")

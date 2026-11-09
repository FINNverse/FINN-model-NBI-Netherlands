#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## FINN-model-NBI-Netherlands | calibration/01_download_nbi_raw.R
##
## WHAT   Downloads the public database of the Dutch national forest inventory
##        (Zevende Nederlandse Bosinventarisatie, NBI-7, with the MFV 2001-2005
##        and NBI-6 2012-2013 tables included) and its documentation from
##        Stichting Probos, verifies the SHA-256 recorded on the day the model
##        was built, unpacks the archive and writes raw/DOWNLOAD.md and the
##        raw-data section of calibration/PROVENANCE.md.
## WHY    The zoo model must be reproducible from an official public source.
##        Probos replaces the archive in place when errata are applied (the
##        README inside lists updates up to July 2024; the file on the server
##        was last modified 2026-11-05), so the checksum pins the vintage: a
##        different file stops here instead of silently preparing other data.
## SOURCE https://www.probos.nl/publicaties/overige/1094-bosinventarisaties
##        (the page WUR's NBI project site points to for the "openbare
##        database"). No account, no click-through. Terms (README.txt in the
##        archive): use permitted with correct attribution; a publication based
##        on NBI-7 data must acknowledge the Dutch ministry (LVVN) and cite the
##        main report (Schelhaas et al. 2022, WOt-rapport 142). Plot locations
##        are given as 1 km INSPIRE grid cells only.
## INPUT  none (network)
## OUTPUT <DATA_DIR>/raw/probos/NBI-7_openbare_database.zip + unz/ (the .accdb and README.txt)
##        <DATA_DIR>/raw/probos/*.pdf (database description, metadata NBI-7 / NBI-6)
##        <DATA_DIR>/raw/SHA256SUMS.txt, <DATA_DIR>/raw/DOWNLOAD.md, calibration/PROVENANCE.md
## RUN    Rscript calibration/01_download_nbi_raw.R   (from the repository root)
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=

## Where the data live, outside git: ../../FINN-zoo-data/nbi relative to the
## repository root (the zoo convention). Change this one line for another layout.
DATA_DIR <- file.path("..", "..", "FINN-zoo-data", "nbi")

PARAMS <- list(
  base_url = "https://www.probos.nl",
  files = list(
    database = list(
      path   = "/images/NBI-7_openbare_database.zip",
      dest   = file.path(DATA_DIR, "raw", "probos", "NBI-7_openbare_database.zip"),
      ## SHA-256 of the archive downloaded 2026-11-09 (server file dated 2026-11-05;
      ## inside: NBI7_openbare database_definitief_UpdateJuli2024.accdb, README dated 2024-07-30)
      sha256 = "3d71e28e485651e5a7d8207286d38588dd7b9ed4f22638bff74d1a0adc5f4960",
      unzip  = TRUE),
    description = list(                     # report appendix 1: what the tables are
      path   = "/images/pdf/rapporten/20230417_Bijlage_1_Beschrijving_database.pdf",
      dest   = file.path(DATA_DIR, "raw", "probos", "20230417_Bijlage_1_Beschrijving_database.pdf"),
      sha256 = "348ae6b15d8fa5134cd9bca589278e39110da79d9533a1852d5c8ef4364a24bf", unzip = FALSE),
    metadata7 = list(
      path   = "/images/pdf/rapporten/20230417_Metadata_NBI_7_bestand.pdf",
      dest   = file.path(DATA_DIR, "raw", "probos", "20230417_Metadata_NBI_7_bestand.pdf"),
      sha256 = "8b349f2e0309100ed1ffeaad5cf06dbbfd49973f52075d5018accb8fb6b0c085", unzip = FALSE),
    metadata6 = list(
      path   = "/images/pdf/overig/NBI6_Metadata_bestand.pdf",
      dest   = file.path(DATA_DIR, "raw", "probos", "NBI6_Metadata_bestand.pdf"),
      sha256 = "4f38981a58764a41cd5c615588987d6bcab7d65d72551634844cacf309cc8d8b", unzip = FALSE)
  ),
  licence  = paste("No named licence. Probos README.txt (2024-07-30): use of the data is permitted with correct",
                   "attribution; a publication based on NBI-7 data must include an acknowledgement of the Dutch",
                   "Ministry of Agriculture, Fisheries, Food Security and Nature (LVVN) and cite the main report."),
  citation = paste("Schelhaas, M.J., Teeuwen, S., Oldenburger, J., Beerkens, G., Velema, G., Kremers, J., Lerink, B.,",
                   "Paulo, M.J., Schoonderwoerd, H., Daamen, W., Dolstra, F., Lusink, M., van Tongeren, K., Scholten, T.,",
                   "Pruijsten, I., Voncken, F., Clerkx, A.P.P.M. (2022). Zevende Nederlandse Bosinventarisatie: methoden",
                   "en resultaten. WOt-rapport 142, Wageningen. Data: Stichting Probos, Database Zevende Nederlandse",
                   "Bosinventarisatie (NBI-7), https://www.probos.nl/publicaties/overige/1094-bosinventarisaties,",
                   "accessed 2026-11-09. Acknowledgement: Ministerie van Landbouw, Visserij, Voedselzekerheid en Natuur (LVVN)."),
  timeout_s = 1800,
  ## the Probos server answers a bare R download agent with 403; a browser-like agent is accepted
  user_agent = "Mozilla/5.0 (compatible; FINN-model-NBI-Netherlands calibration script)"
)

options(timeout = PARAMS$timeout_s, HTTPUserAgent = PARAMS$user_agent)
sha256 <- function(path) unname(tools::sha256sum(path))

lines <- c("# NBI raw data: download record", "",
           sprintf("Written %s by calibration/01_download_nbi_raw.R.", format(Sys.time(), "%Y-%m-%d %H:%M %Z")), "",
           paste("Terms:", PARAMS$licence), "",
           paste("Citation:", PARAMS$citation), "",
           "| file | URL | SHA-256 | size (bytes) | downloaded |", "|---|---|---|---|---|")
sums <- character(0)
for (nm in names(PARAMS$files)) {
  f <- PARAMS$files[[nm]]; url <- paste0(PARAMS$base_url, f$path)
  dir.create(dirname(f$dest), showWarnings = FALSE, recursive = TRUE)
  fresh <- FALSE
  if (!file.exists(f$dest)) {
    message("downloading ", url)
    download.file(url, f$dest, mode = "wb", quiet = TRUE, headers = c("User-Agent" = PARAMS$user_agent))
    fresh <- TRUE
  } else message("present: ", f$dest)
  got <- sha256(f$dest)
  if (!identical(got, f$sha256))
    stop(sprintf("checksum mismatch for %s:\n  expected %s\n  got      %s\nProbos has replaced the file; record the new checksum and re-check the preparation.",
                 basename(f$dest), f$sha256, got))
  if (isTRUE(f$unzip)) {
    unz_dir <- file.path(dirname(f$dest), "unz")
    if (!dir.exists(unz_dir) || !length(list.files(unz_dir))) { dir.create(unz_dir, showWarnings = FALSE); unzip(f$dest, exdir = unz_dir) }
  }
  sums <- c(sums, sprintf("%s  probos/%s", got, basename(f$dest)))
  lines <- c(lines, sprintf("| %s | %s | %s | %d | %s |", basename(f$dest), url, got, file.size(f$dest),
                            if (fresh) format(Sys.Date()) else "2026-11-09 (kept)"))
}
writeLines(sums, file.path(DATA_DIR, "raw", "SHA256SUMS.txt"))
writeLines(lines, file.path(DATA_DIR, "raw", "DOWNLOAD.md"))

## the committed provenance record: raw-data section (03_prepare adds the climate section)
prov <- c("# Provenance of every raw layer", "",
          "Pinned vintage: Probos archive `NBI-7_openbare_database.zip` as served on 2026-11-09 (server file date 2026-11-05;",
          "inside: `NBI7_openbare database_definitief_UpdateJuli2024.accdb`, README dated 2024-07-30, errata to July 2024).",
          "Probos replaces the archive in place when errata are applied; a different checksum is another vintage and needs a",
          "new prepared-data version.", "",
          "## NBI public database (01_download_nbi_raw.R)", "",
          paste("Source: Stichting Probos on behalf of Wageningen Environmental Research and the Dutch ministry (LVVN),",
                "https://www.probos.nl/publicaties/overige/1094-bosinventarisaties. The MS Access database holds the plot and tree",
                "measurements of the MFV (2001-2005), NBI-6 (2012-2013, as corrected during NBI-7) and NBI-7 (2017-2021), with",
                "lookup tables, the per-tree estimates (volume, biomass) and the lists of permanent plots remeasured between cycles."), "",
          paste("Terms:", PARAMS$licence), "",
          paste("Citation:", PARAMS$citation), "",
          "| file | URL | SHA-256 | size (bytes) | downloaded |", "|---|---|---|---|---|",
          lines[(which(lines == "|---|---|---|---|---|") + 1):length(lines)], "")
writeLines(prov, file.path("calibration", "PROVENANCE.md"))
message("done: ", file.path(DATA_DIR, "raw"), "; calibration/PROVENANCE.md (raw section)")

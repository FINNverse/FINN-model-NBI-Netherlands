# Provenance of every raw layer

Pinned vintage: Probos archive `NBI-7_openbare_database.zip` as served on 2026-11-09 (server file date 2026-11-05;
inside: `NBI7_openbare database_definitief_UpdateJuli2024.accdb`, README dated 2024-07-30, errata to July 2024).
Probos replaces the archive in place when errata are applied; a different checksum is another vintage and needs a
new prepared-data version.

## NBI public database (01_download_nbi_raw.R)

Source: Stichting Probos on behalf of Wageningen Environmental Research and the Dutch ministry (LVVN), https://www.probos.nl/publicaties/overige/1094-bosinventarisaties. The MS Access database holds the plot and tree measurements of the MFV (2001-2005), NBI-6 (2012-2013, as corrected during NBI-7) and NBI-7 (2017-2021), with lookup tables, the per-tree estimates (volume, biomass) and the lists of permanent plots remeasured between cycles.

Terms: No named licence. Probos README.txt (2024-07-30): use of the data is permitted with correct attribution; a publication based on NBI-7 data must include an acknowledgement of the Dutch Ministry of Agriculture, Fisheries, Food Security and Nature (LVVN) and cite the main report.

Citation: Schelhaas, M.J., Teeuwen, S., Oldenburger, J., Beerkens, G., Velema, G., Kremers, J., Lerink, B., Paulo, M.J., Schoonderwoerd, H., Daamen, W., Dolstra, F., Lusink, M., van Tongeren, K., Scholten, T., Pruijsten, I., Voncken, F., Clerkx, A.P.P.M. (2022). Zevende Nederlandse Bosinventarisatie: methoden en resultaten. WOt-rapport 142, Wageningen. Data: Stichting Probos, Database Zevende Nederlandse Bosinventarisatie (NBI-7), https://www.probos.nl/publicaties/overige/1094-bosinventarisaties, accessed 2026-11-09. Acknowledgement: Ministerie van Landbouw, Visserij, Voedselzekerheid en Natuur (LVVN).

| file | URL | SHA-256 | size (bytes) | downloaded |
|---|---|---|---|---|
| NBI-7_openbare_database.zip | https://www.probos.nl/images/NBI-7_openbare_database.zip | 3d71e28e485651e5a7d8207286d38588dd7b9ed4f22638bff74d1a0adc5f4960 | 16368482 | 2026-11-09 (kept) |
| 20230417_Bijlage_1_Beschrijving_database.pdf | https://www.probos.nl/images/pdf/rapporten/20230417_Bijlage_1_Beschrijving_database.pdf | 348ae6b15d8fa5134cd9bca589278e39110da79d9533a1852d5c8ef4364a24bf | 87736 | 2026-11-09 (kept) |
| 20230417_Metadata_NBI_7_bestand.pdf | https://www.probos.nl/images/pdf/rapporten/20230417_Metadata_NBI_7_bestand.pdf | 8b349f2e0309100ed1ffeaad5cf06dbbfd49973f52075d5018accb8fb6b0c085 | 23263 | 2026-11-09 (kept) |
| NBI6_Metadata_bestand.pdf | https://www.probos.nl/images/pdf/overig/NBI6_Metadata_bestand.pdf | 4f38981a58764a41cd5c615588987d6bcab7d65d72551634844cacf309cc8d8b | 10435 | 2026-11-09 (kept) |

## Climate: WorldClim 2.1 (03_prepare_pairs_v1.R)

Source: https://www.worldclim.org/data/worldclim21.html (Fick & Hijmans 2017, doi:10.1002/joc.5086), 1970-2000 normals,
30 arc-seconds, files wc2.1_30s_bio.zip and wc2.1_30s_elev.zip, licence CC BY-SA 4.0; values extracted at the centres of
the plots' 1 km INSPIRE cells (EPSG:3035) on the date of PREPARATION.md.

| variable | file | SHA-256 | size (bytes) |
|---|---|---|---|
| tavg_c | wc2.1_30s_bio_1.tif | 382c991cea4a20b42e345fcd5a897788feee05bc0d613dec032ff0cad1f4c375 | 753795566 |
| prec_mm | wc2.1_30s_bio_12.tif | 6df0ef2f57f900733a569fcc4c8debd27b490978c5d82c2970485f0e66abbf28 | 487101349 |
| temp_season | wc2.1_30s_bio_4.tif | 722f65e8fef0ebdb7237fa265ad75bc88c6a2867a696506bdc8828bbdb094d37 | 1196313251 |
| prec_season | wc2.1_30s_bio_15.tif | 9bfde1c593b72b707f24b52d32da292f29a77c5cc5d7bbbb81832df7440296c9 | 1006627399 |
| tmax_warm_c | wc2.1_30s_bio_5.tif | dd770c326664dcd16d8799fdfe28adb65f5e054c7c3692452e1382eba56a37c3 | 378786176 |
| tmin_cold_c | wc2.1_30s_bio_6.tif | fb5532167aa42fed22f03e552a65d1055acd35a946b9c114c903fb5db1218c82 | 398541228 |
| elev_m | wc2.1_30s_elev.tif | 1d57ce9fd3584cea6a227396bc75ef88522fd5346e8487c17370885c0499716e | 342195708 |


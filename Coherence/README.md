# Coherence rasters (not stored in Git)

Copy the 17 SNAP coherence GeoTIFFs here (8 PRE, 1 event-straddling, 8 POST),
then list them in `sar_coherence/config.py`, for example:

    Pre1 = r"Coherence\Coherence_Pre-Event-01_YYYYMMDD_YYYYMMDD_split_Orb_Stack_esd_coh_deb_TC.tif"

Relative paths are read from the repository root. `*.tif` files are ignored by Git.

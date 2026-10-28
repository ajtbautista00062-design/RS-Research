# Active research implementation

Use research_17_rasters/sar_coherence as the authoritative research implementation.

Active entry points:
- run_analysis.py: 17-raster ITSA analysis
- run_rri.py: independent RRI post-processing
- explain_csv.py: CSV_Explained.md generation

Read config.py, RRI_METHOD.md and SPATIAL_ENHANCEMENT.md before changing code.

Older version directories, visualization backups, retired files and standalone
redraw scripts are historical references. Do not use them as implementation
sources or restore them without an explicit user request.

Some older helper modules remain required imports of the current entry points.
Do not remove them merely because their names appear historical.

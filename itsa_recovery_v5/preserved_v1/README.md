# Current correction and boundary settings
The input boundary is the read-only Codes/Boundary/INFANTA_BOUNDARY.kml explicitly
designated by the researcher. All outputs and code edits remain in Simulation.
The reference raster grid and extent are preserved; every pixel outside the KML
polygon union (including its interior rings) is NoData. The boundary is an
administrative extent, not a built-up-area mask.

After alignment and masking, before ANY graph, ITSA or recovery calculation:
small enclosed holes are filled from the nearest finite aligned coherence pixel.
Defaults: maximum missing component area 9 pixels, maximum donor distance 2 pixels.
8-connected missing components touching the outer coverage/window edge are preserved.
A sufficient read halo prevents tile seams. Existing valid coherence is unchanged.
Large gaps, edge-connected background and missing coverage remain missing.
Thresholds are operational spatial settings, not statistical corrections.
No temporal interpolation is performed. Original source TIFFs are unchanged.
These replacements are spatial estimates and may affect trends/classes; compare
FILL_SMALL_HOLES=False for sensitivity, especially with untagged zero background.

Infanta_Coherence_Correction_Flags_<dates>.tif records one band per pair:
0 missing/outside, 1 original valid aligned input, 2 nearest-neighbor estimate.
The observed coherence bands in the ITSA TIFF now contain corrected aligned inputs;
their names retain acquisition dates, but filled cells must not be called observations.
Metadata and audit counts record the correction. Provenance is verified on export.

QGIS sidecar QML files automatically style recovery TIFFs:
- Classes: saturated categorical palette, nearest-neighbor display; labels 1..5.
- Fraction: continuous interpolated color ramp and bilinear display.
Keep QML beside the corresponding TIFF; if necessary use Properties > Style >
Load Style. The separate continuous fraction map is visually smooth. Categorical
pixel codes remain discrete; colored patches in that map are pixel/class boundaries.
PNG/PDF/SVG continuous previews use display-only interpolation, not analytical
resampling. Their display range -1..2 does not clip the actual numerical TIFF.

# Real-data ITSA pilot and baseline recovery

Run from Simulation:
    & .\.venv\Scripts\python.exe .\sar_coherence\run_analysis.py

Outputs:
outputs/Time Series Graphs/series_<dates>_<run>/
outputs/Recovery Map Raster File/series_<dates>_<run>/
outputs/series_<dates>_<run>/  (logs, manifest, summary, run_state.json)

Only use runs whose run_state.json says complete. Old output folders are not deleted.
The actual directory contains 3 distinct PRE, 1 event-straddling, 3 distinct POST files.
Repeated paths are never independent observations. Filename pair dates drive ordering.
Original rasters are read-only. Future additional distinct PRE/POST pairs are supported.

PRIMARY ITSA
y = b0 + b1*t + b2*POST + b3*t*POST, t in elapsed days from EVENT_ORIGIN.
The event-straddling pair is excluded from the primary regression and residuals.
It is displayed as an unfilled diamond. Six fitted observations, four coefficients,
and two residual degrees of freedom make the current pilot exploratory.
Pixel models require all PRE/POST observations; missing event data do not prevent fitting.
Observed bands retain their own support. Count bands retain incomplete PRE/POST counts.
R-squared is undefined for constant series and is exported as NoData.
The first seven bands are actual observations; following bands are named model outputs.
The event-pulse coefficient from the older implementation has been removed.

REGIONAL GRAPHS
Primary fit uses seeded spatial-sample medians on fixed seven-date common support,
not exact whole-raster medians. Sample size and valid-pixel counts are exported.
The available-support graph is descriptive only; changing coverage can bias trends.
IQR is spatial dispersion, not a confidence interval. No temporal values are fabricated.
PRE regression: solid blue. POST: dashed blue. Counterfactual: dashed gray.
Black dashed line: configured event time. Pair intervals are shown below the graph.
PNG 600 DPI, PDF and SVG, plus residual plot, coefficient and fit CSVs.
Classical OLS standard errors, 95% intervals and p-values are exported using NumPy/SciPy
OLS algebra. They assume independent errors, are unstable in this short pilot, and
do not establish significance or causal effects. No autocorrelation/seasonal model.
The study includes shared acquisitions and unverified acquisition comparability.

RASTERS
Combined observation/ITSA Float32 GeoTIFF plus standalone named model rasters:
pre_level_at_event, pre_slope_per_day, post_level_change_at_event,
post_slope_change_per_day, post_slope_per_day, final_minus_pretrend_prediction,
valid_model_observation_count, residual_rmse, model_r_squared.
The final pre-trend deviation is a model diagnostic, NOT the recovery reference.
Model rasters are verified against an independent least-squares solution on real
pixel series when the user runs the pipeline.

BASELINE RECOVERY (researcher-selected definition retained)
B = median valid PRE (minimum 3); S = 1.4826*median(abs(PRE-B)); L = B-E.
Detectable event loss: L >= max(0.10, 2*max(S,0.02)).
Endpoint P = median of the final 3 chronological POST observations, all required valid.
R = (P-E)/L, unclipped, exported only for assessed detectable-loss pixels.
Recovered requires both final POST observations >= B-0.05.
Operational defaults are configurable and require sensitivity/field validation.
This classification is independent of the ITSA counterfactual and is not reconstruction.
Class TIFF: UInt8, NoData 0, palette embedded:
1 recovered (green), 2 partial (amber), 3 not recovered/further degraded (red-orange),
4 no detectable event loss (blue), 5 insufficient POST support after detected loss (gray).
Partial: return criterion fails and endpoint > event. Otherwise class 3.
Class 4 needs a valid baseline/event but does not require POST data.
Missing baseline/event or outside mask -> NoData. Class 5 is not NoData.
Float comparison tolerance 1e-7 only protects representation at thresholds.
Fraction TIFF: Float32, NoData -9999 for classes 0/4/5.
Grid, transform, extent and compression preserved. No coherence baseline cutoff.
ZERO_IS_INVALID=True remains an unverified inherited footprint assumption.

QGIS
Load Infanta_Baseline_Recovery_Classes_<dates>.tif, not the ITSA multiband composite.
Properties > Symbology > Paletted/Unique values > Band 1. Assign labels 1..5 above.
0 is NoData/transparent. Use nearest-neighbor display/resampling for class data.
Add polygon KML with transparent fill and black outline.
Zonal histogram on class Band 1 gives counts. Never average class codes.
Recovered share among assessed loss pixels = N1/(N1+N2+N3); missing if denominator 0.
Report class 4, class 5 and NoData separately.
Run Zonal statistics on the fraction TIFF for numerical recovery summaries.
Set Study_Area in config.py to an existing polygon KML under Simulation to export CSVs.
CSV shares are 0..1; multiply by 100 for percentages. Boundary rule: pixel centers
unless ALL_TOUCHED is explicitly changed. Pixels are spatially dependent.

Map PNG/PDF/SVG is a nearest-neighbor downsampled preview of the actual class raster.
Use the full-resolution GeoTIFF for analysis. The preview has a legend and CRS axes.
No simulation was executed during this code update; runtime/render validation awaits
the user's run. Historical test scripts outside tests/ may refer to old schemas.

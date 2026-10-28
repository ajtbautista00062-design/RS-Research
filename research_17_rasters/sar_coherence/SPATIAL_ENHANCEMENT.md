# 17-raster ITSA spatial enhancement and barangay exports

Implementation date: 29 October 2026. No research workflow or input-raster scan was run during this revision. No original raster, KML, previous output, or Simulation-folder code was modified.

## What the source audit establishes

The existing four-parameter, actual-date segmented ITSA fit remains the analytical model. PRE observations estimate the fixed median baseline and PRE trend. The event-straddling observation remains supporting disturbance evidence; it is displayed but not incorrectly assigned to either pure-phase regression segment. All eight PRE and eight POST observations are available to the fit, subject to their existing validity rules.

Class 5 can arise from fewer than six valid PRE/POST observations, unstable predictions, the existing 0.40 baseline screen, missing late observations, contradictory event evidence, model/observation disagreement, or unresolved effect-size margins. Wide standard errors alone are already quality flags rather than automatic rejection gates. Class 4 requires the implemented no-loss evidence conditions; it is not the general fallback.

The former spatial regularizer merged only adjacent objects already assigned the same class. It could simplify object partitions without increasing Classes 1–3. Fixed PRE-feature bins could also split observations close to a bin boundary. This revision adds PRE-compatible merging across those bin boundaries and a separately exported, model-supported neighborhood refinement.

These are source-code findings, not measurements of why the current real raster has its class frequencies. No empirical dominance, vegetation cause, or accuracy improvement is claimed. The existing adaptive loss rule is retained:

`loss threshold = max(EVENT_LOSS_MINIMUM, PRE_RESIDUAL_NOISE_MULTIPLIER * max(PRE residual MAD sigma, ITSA_NOISE_FLOOR))`.

Defaults remain 0.03, 1.0, and 0.02 respectively. Neither the 0.40 baseline screen nor disturbance thresholds were lowered without calibration evidence. Existing gate counts, exclusion reports and threshold sweeps remain available in future runs. Low coherence is never used as a vegetation label.

## Analytical and spatial rules

The model remains `Y(t)=beta0+beta1*t+beta2*P+beta3*t*P`, with days centered on `EVENT_ORIGIN`. The fixed reference is `B=median(valid PRE)`. Modeled POST coherence is `(beta0+beta2)+(beta1+beta3)*t`; deficit is `B-modeled POST`; progress is the change between modeled first and last supported POST observations. The retained recovery fraction is progress divided by the modeled initial baseline deficit, where the existing disturbance eligibility and positive-denominator rules permit it. Values remain unclipped.

1. Raw pixel ITSA outputs are always retained. `Original_ITSA_Recovery_Classification.tif` records the selected existing analysis mode (default object ITSA), while `Raw_Pixel_ITSA_Classes.tif` is the independent pixel reference. “Original” here means before neighbor refinement, not a reproduction of an older run with different region-size settings.
2. Segmentation uses PRE median, PRE variability, PRE trend and independently supplied barriers. Adjacent PRE fragments can merge if their **full** feature ranges stay within `PRE_FEATURE_BIN_WIDTHS`. No event/POST recovery label enters segmentation. Fixed physical analysis domains bound memory; domain boundaries can still split real landscapes.
3. Existing same-class pooled-object refitting remains available. Objects with evidence scores at or above 0.80 are now protected from that merge as well.
4. Neighbor refinement starts from high-confidence Classes 1–3 in the original selected-mode map. Candidates must independently be Class 1, 2 or 3 in the raw pixel ITSA assessment, with score at least 0.80. A raw Class 4 or genuinely unresolved raw Class 5 cannot acquire recovery status by borrowing a neighbor's label. This conservative requirement is deliberate.
5. A candidate and an accepted neighboring seed must match in proposed class, external stratum and modeled temporal signatures. Four-connected paths cannot traverse NoData, unknown/known land-cover transitions or supplied barriers. Synchronous growth is bounded by `NEIGHBOR_RADIUS_M` (45 m), with a sufficient halo for consistent processing at I/O boundaries.
6. Original high-confidence Classes 1–4 are protected. No erosion, dilation, unrestricted closing, class-frequency balancing, threshold reduction or coherence-value interpolation is performed.
7. Changed pixels use their **own raw fitted metrics and quality/reason codes** in enhanced metric outputs. Original models and their metrics remain separate. The event loss, deficit and progress therefore remain demonstrably tied to fitted ITSA trajectories.

### Evidence score and acceptance

For a classified pixel/object with a classified reason code:

```
C = min(valid_PRE_count/8, valid_POST_count/8)
S = max(adaptive_loss_threshold, recovery_tolerance, minimum_improvement)
U = UNCERTAINTY_MULTIPLIER * max(level_SE, progress_SE, final_prediction_SE)
model_evidence = min(C, 1/(1+U/S))
```

Nonfinite or unresolved classifications receive score zero. Scores are clipped to [0,1]. They are conditional diagnostic heuristics, **not calibrated probabilities** and not a correction for temporal/spatial dependence.

For each growth link, signature similarity is `clip(1-max_abs_signature_difference/NEIGHBOR_SIGNATURE_TOLERANCE,0,1)`. Signatures are baseline, fitted level loss, modeled initial/final coherence, modeled progress and final deficit. POST slope is represented through the trajectory endpoints/progress over the same observation schedule; differences in temporal support also reduce completeness. Proximity is `1-0.2*connected_path_distance/radius`. Neighborhood agreement is the proportion of active same-stratum four-neighbors that support the candidate class and similarity criterion. The accepted reassignment score is the minimum of model evidence, similarity, proximity and agreement, and must be at least 0.80. Agreement must also meet `NEIGHBOR_MIN_AGREEMENT=0.80`.

Because this is a minimum rather than a weighted average, strong proximity cannot compensate for weak temporal evidence. The effect of these provisional spatial tolerances must be assessed against independent evidence in a later authorized run. Broader Classes 1–3 are possible only where the independent pixel evidence was obscured by a regional label; this is not a promise of greater area.

## Configuration

Edit `sar_coherence/config.py`. `Pre1`–`Pre8`, `Event_Straddler`, and `Post1`–`Post8` remain intact. The older `During` alias remains supported for the same single event input. Empty raster paths are not replaced with synthetic data.

| Setting | Default | Implemented meaning |
|---|---:|---|
| `MIN_MAPPING_AREA_M2` | 4000 | Flag small analytical objects and report small contiguous regions; never delete a legitimate small recovery region |
| `TARGET_REGION_AREA_M2` | 30000 | Preferred PRE grouping scale; incompatible features can stop growth much earlier |
| `MAX_REGION_AREA_M2` | 90000 | Maximum analysis-domain/object area; not a cap on connected same-class map regions |
| `MIN_PER_DATE_VALID_COVERAGE` | 0.80 | Object-date medians require this valid fraction; barangay dates below it remain missing in the aggregate fit |
| `EVIDENCE_PROTECTION_SCORE` | 0.80 | Minimum reassignment score and protection threshold; accepted configuration range 0.80–1.00 |
| `NEIGHBOR_RADIUS_M` | 45 | Maximum connected path length, not unrestricted nearest-neighbor filling |
| `NEIGHBOR_SIGNATURE_TOLERANCE` | 0.05 | Signature scale in coherence units; 0.80 score requires considerably closer agreement |
| `NEIGHBOR_MIN_AGREEMENT` | 0.80 | Matching active-neighbor fraction |
| `ENABLE_SPATIAL_REGION_GROWTH` | True | PRE-compatible fragment merging and bounded evidence-guided growth |
| `ENABLE_NEIGHBOR_CLASS_ASSIGNMENT` | True | Enables the separate neighborhood refinement |
| `ENABLE_BUILTUP_PRIORITY` | True | Uses external built-up strata to prevent grouping/growth across land-cover transitions; no preferential lowering of scientific thresholds |
| `PRESERVE_ORIGINAL_CLASSIFICATION` | True | Required; original selected-mode and raw pixel maps remain |
| `EXPORT_SPATIALLY_ENHANCED_CLASSIFICATION` | True | If false, the enhanced comparison product is an unchanged copy |

Older internal names `MIN_REGION_AREA_M2` and `OBJECT_MIN_VALID_FRACTION` are synchronized to the public settings at run validation. Area uses the projected CRS linear-unit conversion and actual affine pixel dimensions. The reference must retain the existing north-up projected-grid requirements. Genuine same-class components can exceed 90,000 m² across several analytical objects; the connected-region report counts them across all object and I/O boundaries.

### Built-up mask

`BUILTUP_MASK` accepts an exact-grid binary TIFF, polygon KML, or polygon GeoJSON. `BUILTUP_MASK_NOTES` must identify provenance, date and interpretation. GeoJSON coordinates use the explicitly configured CRS (default EPSG:4326). TIFF valid values must be 0 or 1; masked values remain unknown. For polygon masks, the researcher must justify that outside-polygons means non-built-up within the study coverage. A municipal or barangay boundary is not a built-up mask.

Without a supplied reliable layer, the study-wide analysis remains active. `Builtup_Status.tif` records -1 unknown; no urban/vegetation class is invented. Known built-up (1), known non-built-up (0) and unknown (-1) have separate area statistics and cannot mix in spatial grouping when built-up priority is enabled. Unknown land cover is explicitly recorded in logs and the manifest.

## Boundary audit and barangay processing

Read-only XML inspection found **36 polygon Placemark features with 36 unique `ADM4_PCODE` values**, using **`ADM4_EN`** for names, in `Boundary/INFANTA_BOUNDARY.kml`. Examples include Abiawin, Agos-agos, Dinahican, Poblacion 1 (Barangay 1) and Tudturan. This is an individual-barangay layer, not merely one municipal outline. The embedded metadata refer to a 2023 boundary version; confirm suitability for the 2022 study.

Full geometric validity and projected overlap checks are implemented but could not be executed in the current Python environment. Therefore 36 is the detected polygon-feature count, **not a claim of 36 geometrically validated polygons**. Runtime audit records invalid/missing features and overlap areas. Invalid geometries are never silently repaired.

Barangay names and codes are retained in CSV metadata and sanitized for Windows directory names. Stable code-sorted, pixel-center ownership assigns a shared boundary/overlap pixel once, preventing double counting. Barangay exports require `ALL_TOUCHED=False` for consistent inclusion. Overlaps are reported for researcher review. This is a deterministic allocation convention, not a substitute for fixing an invalid administrative layer.

For each intersecting barangay, the workflow:

- Clips existing analytical raster results to the grid-aligned bounding window and masks all non-owned centers. It does not refit the pixel/object classification models.
- Exports the requested nine raster products plus existing model/quality products. Original grid, CRS, descriptions, tags, units, class codes and NoData values are retained.
- Computes exact per-date count, mean, median and spatial quartiles from the aligned valid coherence values. Temporary disk-backed values permit exact quantiles without loading all full rasters into memory. Temporary scratch files are removed after use.
- Retains all 17 date rows, both acquisition dates, phase, valid fraction and eligibility. Missing or under-80%-coverage dates are not imputed. A separate descriptive ITSA fit to eligible barangay medians is explicitly labeled; it does not replace the classification model. Means and spatial IQR are displayed; spatial IQR is not a model confidence interval.
- Exports class counts/areas/valid percentages, transition matrices, clipped analysis-object counts, exact connected-region areas, indicator summaries, evidence summaries and built-up/non-built-up/unknown statistics.
- Skips polygons without owned raster centers and reports the reason. Partial raster coverage does not fabricate coverage outside the reference extent. The per-date fraction denominator is owned pixel centers **within the reference raster**, not the full off-raster polygon area.

Class 5 retains the scientifically broader label “Indeterminate.” Insufficient POST support is one reason, but missing PRE support and contradictory evidence also occur. Using “insufficient POST support” as the sole label would misrepresent the implemented class semantics; consult the exported reason rasters/CSV.

## Output organization and provenance

Default output root is `Codes`. By default, a common timestamped run folder beneath each of the three roots preserves previous outputs:

```
Codes/Rasters/<run_id>/*.tif
Codes/ITSA Graphs/<run_id>/*.{png,pdf,svg}
Codes/CSV files/<run_id>/*.csv
Codes/Rasters/<run_id>/Barangay/<actual_name_code>/*.tif
Codes/ITSA Graphs/<run_id>/Barangay/<actual_name_code>/*.{png,pdf,svg}
Codes/CSV files/<run_id>/Barangay/<actual_name_code>/*.csv
```

Set `OUTPUT_VERSIONED_RUNS=False` only for empty output folders if the exact unversioned layout is required. The preflight refuses nonempty target folders; no previous outputs are deleted or overwritten. JSON provenance and logs live alongside CSV reports. All original ITSA raster/CSV outputs remain, and the existing study ITSA figure is retained. Barangay `Coherence_Time_Series`, `ITSA_Fitted_Trends`, and `Pre_Post_Recovery_Trends` filenames currently export the same comprehensive time-plot layout, not three different statistical models.

Primary additions:

| Raster | Meaning |
|---|---|
| `Original_ITSA_Recovery_Classification` | Original selected-mode classification before neighbor refinement; Byte, NoData 0 |
| `Enhanced_ITSA_Recovery_Classification` | Evidence-guided result; same codes 0–5 |
| `Spatial_Recovery_Regions` | PRE-based analytical object IDs, UInt32, NoData 0; not class codes |
| `ITSA_Disturbance_Magnitude` | Fitted `-beta2`, coherence units |
| `ITSA_Recovery_Magnitude` | Modeled POST improvement, coherence units |
| `ITSA_Recovery_Fraction` | Unclipped model-derived fraction, dimensionless |
| `Evidence_Confidence` | Original model evidence for unchanged pixels; accepted link score for changes, dimensionless |
| `Classification_Changes` | Three Byte bands: original class, enhanced class, reason (0 unchanged, 1 supported spatial reassignment); NoData 255 |
| `Valid_Data_Coverage` | 17 Byte bands; 1 observed, 0 missing, 255 outside analysis domain; pair metadata per band |

Continuous rasters use Float32 and -9999 NoData. Extra raw/enhanced multiband metrics, quality/reason rasters, ownership, external strata, built-up status and aligned coherence observations support traceable exports without rereading/reprocessing original input rasters. This increases disk use; processing remains windowed. Source coherence is never filled.

`Pixel_Reassignments.csv` records row, column, old/new class, evidence score and reason for every change. `threshold_adjusted` is false: this revision does not covertly lower disturbance thresholds. `Spatial_Region_Statistics.csv` and `Connected_Region_Counts.csv` describe actual four-connected enhanced class regions, distinct from bounded analytical objects. The 4,000 m² small-region flag, within-20%-of-30,000 m² target flag and above-90,000 m² flag are reporting conventions, not proof of recovery. Evidence medians in the summary are rounded to a documented 0.001 score histogram; per-date coherence medians are exact.

The categorical palette is consistently embedded in raster color tables and comparison graphs: 1 `#167C80`, 2 `#F4C95D`, 3 `#C73A49`, 4 `#D8DDE2`, 5 `#8C78A3`; NoData transparent. Use QGIS Paletted/Unique values for the class rasters and continuous rendering for the magnitude/fraction rasters. Never average categorical class codes.

## Modified files and functions

All paths below are relative to `Codes/sar_coherence`.

| File | Changes / principal functions |
|---|---|
| `config.py` | Requested spatial defaults, masks, barangay fields, evidence controls and versioned output routing |
| `sarflow/research_pipeline.py` | Integrates land-cover strata, ownership, preserved raw evidence, enhancement, clipping, coverage summaries and output routing; `validate_settings`, `run` |
| `sarflow/spatial_objects.py` | `grow_pre_regions` adds PRE-compatible merging; `regularize` protects high-evidence objects |
| `sarflow/spatial_enhancement.py` (new) | `evidence_score`, `refine`, `shift`: independent candidate support and bounded connected growth |
| `sarflow/analysis_boundaries.py` (new) | `load_barangays`, `owners`, `safe_name`, `BuiltupLayer`, `combine_strata` |
| `sarflow/output_layout.py` (new) | `create_layout`, `OutputLayout`: three roots and overwrite protection |
| `sarflow/enhancement_reports.py` (new) | `enhance_outputs`, `Comparison`, `connected_region_report`, `render_comparison` |
| `sarflow/barangay_outputs.py` (new) | `export_barangays`, `export_zone`, `DiskValues`, `graph_series` |
| `sarflow/recovery_io.py` | Additional typed outputs, pair metadata, units and requested palette; existing product definitions retained |
| `sarflow/spatial.py` | Uses actual `ADM4_EN` when available for existing KML zonal labels |
| `sarflow/itsa_figures.py` | CSV routing to table folder, retained plot plus fitted-trend alias, aggregate model table |
| `tests/test_spatial_enhancement.py` (new) | Small deterministic numerical, boundary, routing, coverage and component tests |

The core `itsa.py` and `recovery.py` equations and original class thresholds are unchanged. The 17 editable path variables and legacy observed-workflow dispatch remain intact. No code in `github-publish` or `Simulation` was changed.

## Verification and remaining limits

- Python AST syntax parsing passed for all 38 Python files; no analysis module was executed by the syntax check.
- The new small test module: **11 passed, 7 skipped**. Passed checks cover classes 1–5 from synthetic 17-date trajectories, score behavior, model-supported changes, raw unresolved/no-loss protection, high-confidence preservation, NoData and land-cover barriers, signature disagreement, insufficient support, halo consistency, and safe output directory creation.
- Three SciPy-dependent checks (coverage aggregation, physical object sizes, PRE merging) and four GIS-dependent checks (KML geometry, ownership, disk quantiles, global components) could not run in this environment.
- Attempting the existing combined test suite failed at import: the installed SciPy lacks `scipy._external.packaging_version`. Rasterio, Shapely, PyProj, Matplotlib and DefusedXML are also unavailable in this interpreter. No dependencies were installed or altered, and no successful end-to-end GIS verification is claimed.
- No research rasters, figures or class-frequency reports were generated. No full-dataset or seven-observation run occurred. Actual class-area changes, external accuracy, source zero semantics, acquisition comparability, overlap effects and the suitability of the spatial/evidence defaults require a later authorized validation run.
- Grid-aligned fixed segmentation domains remain a scale assumption. Enhancement can operate across I/O/object boundaries, but does not eliminate all landscape discontinuities introduced by analytical object boundaries. Small high-evidence features remain protected.
- The KML is not a built-up mask. Without an independent mask, built-up prioritization cannot distinguish vegetation from infrastructure.

## Run later, after configuration and runtime verification

Use a Python environment containing the project's requirements and a functioning SciPy/GDAL stack. First run the lightweight test suite there. Fill all 17 actual paths and verify the comparability/NoData settings; do not set confirmation flags merely to bypass validation.

```powershell
Set-Location "C:\Users\Arlan\OneDrive\Desktop\Jhonjhon\Remote Sensing Research\Codes"
if ((Get-Location).Path -ne "C:\Users\Arlan\OneDrive\Desktop\Jhonjhon\Remote Sensing Research\Codes") { throw "Wrong project directory" }
python -B .\sar_coherence\run_analysis.py
```

This command was not executed during implementation. It intentionally does not delete old outputs.

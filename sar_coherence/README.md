# Sentinel-1 coherence ITSA / RRI recovery workflow (Infanta, Quezon)

Typhoon Noru (Karding), September 2022. Inputs are terrain-corrected 12-day
Sentinel-1 coherence GeoTIFFs from SNAP. The workflow fits an interrupted time
series (ITSA) to every pixel, PRE-only object and multi-scale zone, classifies
recovery, computes the Recovery Rate Index (RRI), and exports municipal and
per-barangay rasters, CSV tables, graphs and QGIS styles.

## 1. Setup

```powershell
python -m venv .venv
.venv\Scripts\python.exe -m pip install -r sar_coherence\requirements.txt
```

## 2. Configure (edit only `config.py`)

- `Pre1` … `Pre8`: PRE-event coherence rasters (leave unused slots empty; ≥3 required, 8 recommended)
- `Event_Straddler` (or `During`): the pair that spans the typhoon (required)
- `Post1` … `Post8`: POST-event coherence rasters (≥3 required, 8 recommended)
- `PAIR_DATES`: leave blank to read `YYYYMMDD_YYYYMMDD` from each file name
- `ACQUISITION_COMPARABILITY_CONFIRMED = True` and `COMPARABILITY_NOTES` (track,
  orbit, polarization, bursts, processing) — the run refuses to start without them
- `Output_Folder`: runs go to `<Output_Folder>\Results\<run name>`; keep it short
  (Windows 260-character path limit; the run checks this before starting)

All rasters, the boundary KML and the output folder must be inside `PROJECT_ROOT`
(the folder that contains `sar_coherence`).

## 3. Run

```powershell
cd sar_coherence
$env:PYTHONPATH = (Get-Location).Path
..\.venv\Scripts\python.exe -m unittest discover -s tests -p "test_*.py"
..\.venv\Scripts\python.exe run_study.py --run-name run01
```

`run_study.py` performs the full study and then the zonal product.
`run_zonal.py --run run01 --replace` re-creates only the zonal files of an
existing run. `run_analysis.py` is an alias of `run_study.py`.

## 4. Outputs (`Results\<run>\`)

```
A. RRI
  Raster Files\Overall | Per Barangay\<name_code>
  Graphs\Overall       | Per Barangay\<name_code>
  CSV\Overall          | Per Barangay\<name_code>
B. ITSA
  Overall              \ Raster Files | ITSA Graph | CSV
  Per Barangay\<name_code>\ Raster Files | ITSA Graph | CSV
Validation, run_manifest.json, processing.log, run_state.json
```

| Product | Where |
|---|---|
| Pixel RRI value, interpretation classes, validity status | `A. RRI` (`ITSA_RRI`, `RRI_*`) |
| **Zonal RRI** value, interpretation classes, validity status, class statistics, map | `A. RRI` (`Zonal_RRI*`) |
| Object ITSA classes, raw pixel classes, `Original_*` classes (historical 0.40 rule), reasons, quality flags, metrics, aligned stack | `B. ITSA` |
| **Zonal recovery classes** (primary map), disturbance evidence, Indeterminate reasons, quality flags, zone size used, zone ID; class areas (ha), target assessment, transitions, map | `B. ITSA` (`Zonal_*`) |
| Zone polygons with attributes (WGS84) | `B. ITSA\Overall\Raster Files\Zonal_Recovery_Zones.geojson` |
| Zonal sensitivity, zone attributes, all-barangay target table, method notes | `B. ITSA\Overall\CSV\Zonal_*` |
| Fixed-population original→corrected transitions, gate outcomes, RRI transitions | `Validation` |

Class rasters carry an embedded colour table and a same-name `.qml` style, so
QGIS shows the legend automatically.

Recovery classes: 0 NoData, 1 Coherence recovered, 2 Partially recovered,
3 Not recovered / further degraded, 4 No detectable event-associated loss,
5 Indeterminate.

## 5. Method summary

**ITSA model.** Segmented OLS on actual pair dates centred on landfall
(`EVENT_ORIGIN`): `y = β0 + β1 t + β2 POST + β3 t·POST`. The event-straddling
pair is excluded from the fit and used as observed disturbance evidence.
Baseline B is the median of valid PRE pairs; the adaptive loss threshold is
`H = max(0.03, 1.0 × max(σ_res, 0.02))`.

**Observability gate (changed).** The historical rule required B ≥ 0.40, which
left ~97% of Infanta Indeterminate (vegetated land has 12-day coherence ≈ 0.10–0.30).
`BASELINE_ELIGIBILITY_MODE = "signal_evidence"` instead requires
B > tolerance + H + margin and B > 3 × residual noise. The 0.40 rule is still
exported as `Original_*` products and quality bit 64.

**Pixel / object ITSA.** Unchanged decision rules from the previous version,
applied to pixels and to fixed PRE-only homogeneous objects.

**Multi-scale zonal ITSA (new, `sarflow/zonal_recovery.py`).** Pixel coherence is
too noisy for stable decisions, so maps were speckled. Square zones aligned to
the grid (70 m base, doubling to 140, 280, 560 m; split by barangay) are
classified from the median coherence of their member pixels. Each pixel takes
the class of the finest zone that resolves it. Decision tree per zone:

- Indeterminate: insufficient PRE/POST/event support, unstable fit, or failed observability gate
- Disturbed: event-pair loss B − E ≥ H, or modeled level loss ≥ H confirmed by early-POST deficit ≥ H
- Disturbed → late median ≥ B − 0.05: **Recovered**; late − trough ≥ 0.03: **Partial**; else **Not recovered**
- Not disturbed but early and late POST deficits persist: **Not recovered** (post-event decline)
- Otherwise: **No detectable event-associated loss**

Zonal RRI = (late − trough) / (B − trough) for disturbed zones.

**Simulation result (7 rasters).** Zonal Indeterminate fell from 99.0% to 0.5%
of study pixels (99.5% reduction); classes 1–3 rose from 5,228 to 373,533 pixels.
Every barangay exceeded a 91% reduction, and all sensitivity scenarios
(42/98 m base zones, tolerance 0.03/0.07, loss floor 0.05, no barangay split)
stayed above 98.9%.

## 6. Limitations

- Event-pair decorrelation can come from rain, wind or soil moisture rather than
  damage; quality bit 1 marks zones whose only disturbance evidence is the event
  pair and that were already near baseline in the first POST pair.
- Low coherence magnitudes are biased upward by the estimator; this is not corrected.
- Zones mix land covers within their extent; zone size is exported per pixel.
- ITSA standard errors are descriptive; shared acquisitions correlate errors.
- No field validation: outputs are exploratory, not accuracy-assessed maps.
- RRI and ITSA classes measure recovery differently (RRI is a fraction of the
  loss; ITSA "Recovered" allows 0.05 below baseline), so they can disagree.

## 7. Code layout

```
config.py              user settings (17 input slots)
run_study.py           full workflow entry point
run_zonal.py           zonal product only, for an existing run
sarflow/
  study_inputs.py      config validation, raster alignment
  study_pipeline.py    windowed processing, run manifest
  itsa.py              segmented OLS
  recovery.py          pixel/object classifier (signal observability gate)
  baseline_recovery.py historical 0.40 classifier, kept for comparison
  spatial_objects.py   PRE-only object segmentation
  study_objects.py     object aggregation
  rri.py               Recovery Rate Index
  zonal_recovery.py    multi-scale zonal classification and exports
  study_exports.py, study_tables.py, study_figures.py, study_rasters.py,
  study_comparison.py, study_verify.py, study_layout.py   outputs and checks
  analysis_boundaries.py, spatial.py, temporal.py, recovery_io.py,
  spatial_enhancement.py                                  shared helpers
tests/                 unit and end-to-end tests (python -m unittest discover -s tests)
```

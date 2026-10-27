# ITSA baseline recovery V3

Original V1/V2 sources, input rasters and previous results remain preserved. Edit this folder's config.py. Primary analysis uses unfilled observations. Spatial correction is a separate sensitivity. Read METHOD.md and TECHNICAL_REPORT.md for exact rules and measured results.

## Run from PowerShell in Simulation

```powershell
& .\.venv\Scripts\python.exe -B .\itsa_recovery_v3\run_analysis.py --stage pilot
```

To preserve and verify an archive, delete the active previous Corrected_outputs folder, then launch V3:

```powershell
.\Regenerate_Corrected_Outputs_v3.ps1
```

The reset script checks workspace paths, copies previous outputs into Output_Archive, SHA256-verifies every file, then deletes only the active Corrected_outputs copy. It has been parsed, not executed.

New results go into a timestamped itsa_v3_pilot_* folder under Corrected_outputs, containing Time Series Graphs, Recovery Map Raster File, Sensitivity and Reports. Accept only source runs whose run_state.json says complete. No new dependencies are needed.

Tests only:

```powershell
& .\.venv\Scripts\python.exe -B -m unittest discover -s .\itsa_recovery_v3\tests -p "test_*.py" -v
```

Existing-output diagnostics completed in Corrected_outputs/Class4_Diagnostics_v3_20261027/V3_Final_Assessment. The earlier V3_Assessment is an exploratory coupled-threshold experiment. The final assessment keeps event tolerance fixed while sweeping modeled-loss thresholds. No full original-source workflow or 17-observation run was executed for this revision.

Full mode requires a completed validated V3 pilot, documented zero semantics, verified acquisition comparability and ALLOW_FULL_DATASET. A cached diagnostic assessment is not a validated original-source V3 pilot.

## QGIS

Load Recovery Map Raster File/ITSA_Recovery_Classes.tif, Band 1, Paletted/Unique values. Use its matching QML: 1 green, 2 amber, 3 vermilion, 4 blue, 5 gray, 0 transparent. Keep categorical resampling at nearest neighbor; do not smooth class codes.

Continuous indicators such as baseline_deficit_final and modeled_recovery_progress use Singleband pseudocolor. Positive deficit means below baseline; positive progress means movement toward baseline. Values are coherence units, not recovery percentages.

Use Zonal histogram for class shares, never means of class codes. Use Zonal statistics for continuous indicators, reporting valid counts. Inspect quality_flags and exclusion_reason. Threshold_Sensitivity.tif: 1 stable, 2 changes across candidates, 0 outside. Candidate disturbance-evidence maps: 1 modeled loss, 2 bounded absence, 3 ambiguous, 4 rejected/conflicting, 0 outside; these are not recovery classes.

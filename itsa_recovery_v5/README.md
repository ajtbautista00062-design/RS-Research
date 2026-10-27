# V5 operational ITSA baseline recovery

This is the current source revision. V1-V4 and previous outputs remain preserved. No V5 simulation or tests were executed during this coding task. Read OPERATIONAL_METHOD.md.

Edit itsa_recovery_v5/config.py. Defaults: CLASSIFICATION_MODE="object_itsa", CLASSIFICATION_EVIDENCE_MODE="operational". Pixel, object and regularized modes remain available. "guarded" selects historical trajectory gates for comparison on the same fit. Quality and rejection outputs are mandatory audit products, so disabling those two export settings raises an explicit error.

Later, from PowerShell in Simulation:

```powershell
.\Run_Operational_ITSA.ps1
```

Or:

```powershell
& .\.venv\Scripts\python.exe -B .\itsa_recovery_v5\run_analysis.py --stage pilot
```

This will run software tests before the actual pilot and create a NEW timestamped Corrected_outputs/itsa_v5_* folder. Previous outputs are not deleted. The older V3/V4 launchers still select older implementations. Full mode requires a completed validated V5 pilot in the same classification/evidence mode, resolved source semantics and acquisition comparability.

Tests only, for later authorization:

```powershell
& .\.venv\Scripts\python.exe -B -m unittest discover -s .\itsa_recovery_v5\tests -p "test_*.py" -v
```

Primary results: Recovery Map Raster File/ITSA_Recovery_Classes.tif, separate ITSA_Evidence_Status and uncertainty/count bands, region statistics, gate/rejection counts, distributions, pixel/object agreement and threshold tables. Pixel ITSA Reference retains its operational and same-fit guarded maps. ENABLE_THRESHOLD_SENSITIVITY_REPORT controls future sensitivity work.

QGIS: primary classes use Paletted/Unique values (matching QML): 1 green, 2 amber, 3 vermilion, 4 blue, 5 gray, 0 transparent. Evidence status: 0 unresolved, 1 guarded-rule agreement (conditional), 2 provisional. Display uncertainty alongside classification; do not average categorical codes. Crossing status distinguishes observed bracketed timing from censored or unobservable timing; a Class-1 endpoint does not guarantee a duration estimate.

All thresholds, spatial masks, pair comparability and recovery interpretation require independent calibration. No increase in Classes 1-3 has been measured for V5.

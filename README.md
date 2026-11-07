# RS-Research: Sentinel-1 coherence recovery (Infanta, Quezon)

ITSA and RRI recovery assessment after Typhoon Noru (Karding), September 2022,
with the multi-scale zonal recovery classification. Full documentation:
[`sar_coherence/README.md`](sar_coherence/README.md).

## Quick start (Windows, PowerShell)

Clone to a SHORT path (Windows limits paths to 260 characters), e.g. `C:\RS`:

```powershell
git clone https://github.com/ajtbautista00062-design/RS-Research.git C:\RS
cd C:\RS
python -m venv .venv
.venv\Scripts\python.exe -m pip install -r sar_coherence\requirements.txt
```

Edit `sar_coherence\config.py`:

1. Put the full path of each raster in `Pre1`…`Pre8`, `Event_Straddler`, `Post1`…`Post8`
   (rasters may be on any drive; relative paths are read from the repository root).
2. Set `ACQUISITION_COMPARABILITY_CONFIRMED = True` and describe track, orbit,
   polarization, bursts and processing in `COMPARABILITY_NOTES`.

Then run:

```powershell
cd C:\RS\sar_coherence
$env:PYTHONPATH = (Get-Location).Path
..\.venv\Scripts\python.exe -m unittest discover -s tests -p "test_*.py"
..\.venv\Scripts\python.exe run_study.py --run-name run17
```

Results are written to `C:\RS\Results\run17\` (`A. RRI`, `B. ITSA`, `Validation`).
The primary map is `B. ITSA\Overall\Raster Files\Zonal_Recovery_Classes.tif`, and
`B. ITSA\Overall\CSV\Zonal_Target_Assessment.csv` reports the Indeterminate
reduction and the recovery-class increase against the original rules.

Repository layout:

```
Boundary/        Infanta barangay polygons (KML) used for study mask and per-barangay outputs
Coherence/       optional place for the rasters (*.tif is not committed)
sar_coherence/   code, config.py, tests
Results/         run outputs (not committed)
```

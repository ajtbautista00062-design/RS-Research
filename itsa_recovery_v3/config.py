"""V3 settings: all generated results stay under Simulation/Corrected_outputs."""
from pathlib import Path
from types import SimpleNamespace

PROJECT_ROOT=Path(__file__).resolve().parent.parent
INPUT_DIRECTORY=PROJECT_ROOT/"Coherence"
Output_Folder=str(PROJECT_ROOT/"Corrected_outputs")
Study_Area=str(PROJECT_ROOT.parent/"Codes"/"Boundary"/"INFANTA_BOUNDARY.kml")
EVENT_START="2022-09-24T16:00:00"
EVENT_END="2022-09-26T16:00:00"
EVENT_ORIGIN="2022-09-25T09:30:00"
REFERENCE_DATE="end"
REFERENCE_RASTER=""
BAND=1
RESAMPLING="bilinear"
COHERENCE_SCALE=1.
COHERENCE_OFFSET=0.
ALL_TOUCHED=False
WINDOW_SIZE=256
SAMPLE_PIXELS=20000
RANDOM_SEED=20220925
FIGURE_FORMATS=("png","pdf","svg")
FIGURE_DPI=600
NODATA=-9999.
# No claim that untagged source zeros have been resolved.
# retain: use GDAL masks; zero remains valid inside them.
# background: remove zeros ONLY after a documented source/footprint justification.
ZERO_SEMANTICS="unresolved"  # unresolved / valid / background
ZERO_SEMANTICS_JUSTIFICATION=""
ACQUISITION_COMPARABILITY_CONFIRMED=False
COMPARABILITY_NOTES="Orbit, polarization, bursts and processing remain unverified."
# Runtime stage is selected explicitly; never automatically execute the full dataset.
EXPECTED_PILOT_COUNTS=(3,1,3)
EXPECTED_FULL_COUNTS=(8,1,8)
ALLOW_FULL_DATASET=False  # change only after reviewing completed pilot validations
RUN_SENSITIVITY=True
DETECTION_SWEEP=(.01,.02,.03,.05,.10)
RECALIBRATE_FROM_PRE=True
CALIBRATED_LOSS_MIN=.045887350471781156
CALIBRATION_SAMPLE_PIXELS=50000
REAL_PIXEL_CHECK_LIMIT=100
# Discovery compatibility with the preserved helpers (no legacy classifier is enabled).
ENABLE_BASELINE_CLASSIFICATION=False
ZERO_IS_INVALID=False
FILL_SMALL_HOLES=False
FILL_MAX_HOLE_PIXELS=9
FILL_RADIUS_PIXELS=2

def analysis_settings():
    return SimpleNamespace(
        MIN_PRE=3,MIN_POST=3,REQUIRE_COMPLETE=True,CONFIRM_POST=2,
        LOSS_MIN=CALIBRATED_LOSS_MIN,LOSS_SPREAD_MULTIPLIER=2.,PRE_SPREAD_FLOOR=.02,
        DETECTION_ALPHA=.05,EVENT_EVIDENCE_MIN=CALIBRATED_LOSS_MIN,LOSS_CALIBRATION_VALIDATED=False,
        CALIBRATION_METHOD="PRE-only detrended residual median",
        RETURN_TOLERANCE=.05,MIN_PROGRESS=.03,
        MAX_RMSE=.10,MAX_PRE_SPREAD=.15,MAX_PREDICTION_EXCESS=.05,
        ERROR_FLOOR=.005,BASELINE_GUARD_FLOOR=.005,GUARD_MULTIPLIER=2.,
        USE_UNCERTAINTY_GUARD=True,MAX_CONDITION=1e6,
        ZERO_UNRESOLVED=ZERO_SEMANTICS=="unresolved",
        COMPARABILITY_UNVERIFIED=not ACQUISITION_COMPARABILITY_CONFIRMED)

def reader_settings(exclude_zero=False,filled=False):
    values={k:v for k,v in globals().items() if k.isupper() or k in ("Output_Folder","Study_Area")}
    values["ZERO_IS_INVALID"]=exclude_zero
    values["FILL_SMALL_HOLES"]=filled
    return SimpleNamespace(**values)

"""V5: operational ITSA classes with separate statistical evidence products."""
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
RECOVERY_MODEL="itsa"
RECOVERY_REFERENCE="pre_median"
EVENT_LOSS_METHOD="model_and_observation"
CLASSIFICATION_EVIDENCE_MODE="operational"  # operational / guarded (V4 comparison)
EXPORT_QUALITY_FLAGS=True  # required audit product; disabling is rejected explicitly
EXPORT_REJECTION_REASONS=True  # required by spatial QA and rejection diagnostics
ENABLE_THRESHOLD_SENSITIVITY_REPORT=True
EVENT_LOSS_MINIMUM="pre_noise"  # auto PRE calibration, or a documented numeric floor
EARLY_POST_LOSS_MINIMUM=.03
LOSS_NOISE_MULTIPLIER=.50
RECOVERY_MIN_IMPROVEMENT=.03
BASELINE_RETURN_TOLERANCE=.05
EFFECT_DECISION_MARGIN=.005  # near-threshold effects remain indeterminate
OBSERVATION_CONFIRMATION_TOLERANCE=.02
MAX_OBSERVATION_DISAGREEMENT=.10
MODEL_RMSE_FAILURE=.20  # hard failure; RMSE>.10 remains a separate quality warning
MIN_PRE_OBSERVATIONS=3
MIN_POST_OBSERVATIONS=3
CONFIRM_POST_OBSERVATIONS=2
REQUIRE_COMPLETE_SERIES=False  # all valid observations used; no fabricated dates
REQUIRE_EVENT_OBSERVATION=False  # missing event requires early modeled AND observed loss
# Discovery compatibility with the preserved helpers (no legacy classifier is enabled).
ENABLE_BASELINE_CLASSIFICATION=False
ZERO_IS_INVALID=False
FILL_SMALL_HOLES=False
FILL_MAX_HOLE_PIXELS=9
FILL_RADIUS_PIXELS=2

# Preferred analytical map. Pixel evidence remains independently available.
CLASSIFICATION_MODE="object_itsa"  # pixel_itsa / object_itsa / spatially_regularized_itsa
ENABLE_SPATIAL_SEGMENTATION=True
ENABLE_SPATIAL_REGULARIZATION=True  # used only in spatially_regularized_itsa mode
MIN_REGION_AREA_M2=2000.0  # provisional ~10 pixels at 14 m; not a damage-size claim
TARGET_REGION_AREA_M2=10000.0
MAX_REGION_AREA_M2=40000.0
PRE_FEATURE_TOLERANCES=(.12,.06,.06)  # local PRE median, temporal MAD, spatial MAD
PRE_NEIGHBORHOOD_RADIUS=1
OBJECT_MIN_VALID_PIXELS=3
OBJECT_MIN_VALID_FRACTION=.80
OBJECT_MAX_SPATIAL_MAD=.10  # heterogeneity screen, not a recovery threshold
OBJECT_HETEROGENEITY_POLICY="quality_only"  # quality_only / strict
MMU_POLICY="flag_only"  # flag_only / defer_low_confidence (historical comparison)
CLASS_CONFIDENCE_THRESHOLD=.80  # operational evidence score, NOT probability
REGULARIZATION_STRENGTH=.50  # [0,1] scales PRE similarity tolerances for pair merges
PRESERVE_HIGH_CONFIDENCE_CLASSES=True
PRESERVE_NODATA_BOUNDARIES=True
EXPORT_PIXEL_LEVEL_REFERENCE=True
MAX_SPATIAL_OBJECTS=200000
# Optional EXACT-GRID categorical rasters; no interpolation of boundary classes.
# Positive equal IDs permit adjacency; zero/masked cells are barriers.
SPATIAL_BARRIER_RASTER=""
VERIFIED_BUILTUP_MASK=""
BUILTUP_MASK_JUSTIFICATION=""  # mandatory if VERIFIED_BUILTUP_MASK is configured
# Municipality boundary alone does not establish built-up extent or river barriers.

def analysis_settings():
    loss=CALIBRATED_LOSS_MIN if EVENT_LOSS_MINIMUM=="pre_noise" else float(EVENT_LOSS_MINIMUM)
    return SimpleNamespace(
        EVIDENCE_MODE=CLASSIFICATION_EVIDENCE_MODE,
        EARLY_LOSS_MIN=EARLY_POST_LOSS_MINIMUM,LOSS_NOISE_MULTIPLIER=LOSS_NOISE_MULTIPLIER,
        EFFECT_MARGIN=EFFECT_DECISION_MARGIN,OBS_CONFIRM_TOLERANCE=OBSERVATION_CONFIRMATION_TOLERANCE,
        MAX_OBSERVATION_DISAGREEMENT=MAX_OBSERVATION_DISAGREEMENT,MODEL_RMSE_FAILURE=MODEL_RMSE_FAILURE,
        REQUIRE_EVENT=REQUIRE_EVENT_OBSERVATION,
        MIN_PRE=MIN_PRE_OBSERVATIONS,MIN_POST=MIN_POST_OBSERVATIONS,
        REQUIRE_COMPLETE=REQUIRE_COMPLETE_SERIES,CONFIRM_POST=CONFIRM_POST_OBSERVATIONS,
        LOSS_MIN=loss,LOSS_SPREAD_MULTIPLIER=2.,PRE_SPREAD_FLOOR=.02,
        DETECTION_ALPHA=.05,EVENT_EVIDENCE_MIN=CALIBRATED_LOSS_MIN,LOSS_CALIBRATION_VALIDATED=False,
        CALIBRATION_METHOD="PRE-only detrended residual median",
        RETURN_TOLERANCE=BASELINE_RETURN_TOLERANCE,MIN_PROGRESS=RECOVERY_MIN_IMPROVEMENT,
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

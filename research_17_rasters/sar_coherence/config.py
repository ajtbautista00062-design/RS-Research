"""Edit ONLY this file for a research run. Empty dates are intentional."""
from pathlib import Path

# PRE-DISASTER COHERENCE RASTERS
Pre1 = r""
Pre2 = r""
Pre3 = r""
Pre4 = r""
Pre5 = r""
Pre6 = r""
Pre7 = r""
Pre8 = r""
# DURING-DISASTER / EVENT-ASSOCIATED RASTER
During = r""
# Preferred name for the same event-straddling input; During remains compatible.
Event_Straddler = r""
# POST-DISASTER COHERENCE RASTERS
Post1 = r""
Post2 = r""
Post3 = r""
Post4 = r""
Post5 = r""
Post6 = r""
Post7 = r""
Post8 = r""
# OUTPUT DIRECTORY: must be inside PROJECT_ROOT; use a NEW folder for each run.
Output_Folder = r""
# OPTIONAL STUDY AREA / ZONAL POLYGONS: polygon KML, WGS84 lon,lat.
Study_Area = r""

PROJECT_ROOT = Path(__file__).resolve().parent.parent
PATHS = {
    **{f"Pre{i}": globals()[f"Pre{i}"] for i in range(1, 9)},
    "During": Event_Straddler or During,
    **{f"Post{i}": globals()[f"Post{i}"] for i in range(1, 9)},
}
# (master, slave), ISO YYYY-MM-DD or ISO timestamp; all UTC, timezone-naive.
# Master can be later than slave: code preserves roles and uses min/max dates.
PAIR_DATES = {name: ("", "") for name in PATHS}
# Replace entries, for example PAIR_DATES["Pre1"] = ("actual date", "actual date").
# No example acquisition dates are supplied because they must come from your data.

# Event window in UTC: calendar-day envelope for 25-26 September in Philippines.
# Refine with local impact evidence; event origin is configurable, NOT an image date.
EVENT_START = "2022-09-24T16:00:00"
EVENT_END = "2022-09-26T16:00:00"
EVENT_ORIGIN = "2022-09-25T09:30:00"  # PAGASA Burdeos landfall; proxy for local impact.
REFERENCE_DATE = "end"  # end recommended; midpoint is a sensitivity scenario.

# Must be confirmed by inspecting SLC/SNAP processing records, not inferred from TIFF.
ACQUISITION_COMPARABILITY_CONFIRMED = False
COMPARABILITY_NOTES = ""  # track, orbit direction, polarization, bursts, window, processing
MAX_PAIR_BASELINE_RATIO = 1.25
ALLOW_BASELINE_MISMATCH = False  # requires substantive rationale in notes; see guide.

MIN_PRE = 6
BASELINE_MIN = 0.40
BASELINE_MAD_MAX = 0.15
NOISE_FLOOR = 0.02
LOSS_ABSOLUTE = 0.10
LOSS_MAD_MULTIPLIER = 2.0
RECOVERY_TOLERANCE = 0.05  # absolute coherence units, baseline minus tolerance
ENDSTATE_POST_PAIRS = 3  # median of the last three post-event coherence pairs
RECOVERY_CONFIRMATIONS = 2  # last two consecutive post pairs must meet baseline threshold
MIN_POST = 8  # conservative default: all eight post observations required
CONFIRMATIONS = 2
MIN_CONFIRM_SPAN_DAYS = 0.0
MAX_CONFIRM_GAP_DAYS = 60.0  # assess against actual sampling; no silent gap bridging
ALLOW_OVERLAPPING_CONFIRMATION_PAIRS = False  # touching endpoints allowed, not independent

BAND = 1
RESAMPLING = "bilinear"  # use average for deliberate coarsening; nearest sensitivity
COHERENCE_SCALE = 1.0  # AFTER TIFF scale/offset; use 0.0001 only if metadata lacks scaling
COHERENCE_OFFSET = 0.0
ALL_TOUCHED = False  # center rule; true is a boundary sensitivity, not fractional area
WINDOW_SIZE = 256
SAMPLE_PIXELS = 20000  # deterministic uniform priority sample of complete-case pixels
RANDOM_SEED = 20220925
SENSITIVITY_TOLERANCES = (0.03, 0.05, 0.07)
SENSITIVITY_CONFIRMATIONS = (2, 3)
FIGURE_FORMATS = ("png", "svg", "pdf")
FIGURE_DPI = 350
INTERACTIVE_HTML = True  # skips with a logged message if Plotly is unavailable
NODATA = -9999.0

# ============================================================
# ITSA DISTURBANCE AND RECOVERY CLASSIFICATION PARAMETERS
# ============================================================
RECOVERY_MODEL = "itsa"  # "legacy_observed" preserves the earlier workflow
CLASSIFICATION_MODE = "object_itsa"  # pixel_itsa / object_itsa / spatially_regularized_itsa
RECOVERY_REFERENCE = "pre_median"  # fixed observed PRE median, not a moving counterfactual
USE_FITTED_ITSA_FOR_RECOVERY = True
USE_ADAPTIVE_DISTURBANCE_THRESHOLD = True
EVENT_LOSS_MINIMUM = 0.03  # provisional effect-size floor; calibrate independently
PRE_RESIDUAL_NOISE_MULTIPLIER = 1.0
ITSA_NOISE_FLOOR = 0.02
EVENT_SUPPORT_MINIMUM = 0.03
MIN_RECOVERY_IMPROVEMENT = 0.03
EFFECT_DECISION_MARGIN = 0.005
ITSA_BASELINE_MINIMUM = BASELINE_MIN  # retain existing 0.40 observability screen
ITSA_MIN_PRE = 6  # uses ALL valid PRE dates, not just six
ITSA_MIN_POST = 6  # uses ALL valid POST dates, not just the endpoints
LATE_POST_CONFIRMATIONS = 2
MIN_LATE_CONFIRM_SPAN_DAYS = 0.0
OBSERVATION_DISAGREEMENT_TOLERANCE = 0.10
MAX_ITSA_RMSE = 0.20  # hard model inadequacy screen, not a significance test
MAX_DESIGN_CONDITION = 1e6
UNCERTAINTY_MULTIPLIER = 2.0  # SE screen only, NOT a calibrated confidence interval
ZERO_IS_INVALID = False  # metadata masks still apply; zero is a valid magnitude
SOURCE_VALIDITY_CONFIRMED = False
SOURCE_VALIDITY_NOTES = ""  # explain declared NoData, zeros and any prior spatial filling

# PRE-only, four-connected regions in physically sized, fixed analysis cells.
# Cells are independent of I/O tiles. Minimum area flags objects; never dilates them.
TARGET_REGION_AREA_M2 = 10000.0
MIN_REGION_AREA_M2 = 2000.0
MAX_REGION_AREA_M2 = 40000.0
PRE_FEATURE_BIN_WIDTHS = (0.10, 0.05, 0.10)  # PRE median, MAD sigma, PRE fitted total change
OBJECT_MIN_VALID_FRACTION = 0.80
OBJECT_MIN_VALID_PIXELS = 3  # capped at the actual region member count
ANCILLARY_REGION_BARRIERS = r""  # optional exact-grid positive integer labels; 0=excluded
ENABLE_SPATIAL_REGULARIZATION = True
REGULARIZATION_MIN_NEIGHBORS = 2
REGULARIZATION_FEATURE_TOLERANCE = 0.05
EXPORT_RAW_CLASSIFICATION = True
EXPORT_REGULARIZED_CLASSIFICATION = True
ENABLE_THRESHOLD_SENSITIVITY_REPORT = True
ITSA_LOSS_SENSITIVITY = (0.01, 0.02, 0.03, 0.05, 0.10)

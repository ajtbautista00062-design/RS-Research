"""Edit ONLY this file for a research run.

Fill in up to 8 PRE, exactly 1 event-straddling, and up to 8 POST coherence
rasters (the full study uses 8 + 1 + 8 = 17). Leave unused Pre/Post slots empty
(at least 3 PRE and 3 POST are required).

Raster paths may be absolute (any drive, e.g. r"D:\\SAR\\Coherence_Pre-Event-01_20220601_20220613.tif")
or relative to PROJECT_ROOT (the repository folder that contains sar_coherence/).
The boundary KML and the output folder must lie inside PROJECT_ROOT.
"""
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
# DURING-DISASTER / EVENT-STRADDLING RASTER (pair spanning the typhoon)
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

PROJECT_ROOT = Path(__file__).resolve().parent.parent
# OUTPUT FOLDER: runs are written to <Output_Folder>/Results/<run name>.
# Keep this path short: Windows paths are limited to 260 characters.
Output_Folder = str(PROJECT_ROOT)
RESULTS_ROOT = Path(Output_Folder).resolve()
VERSION_ROOT = RESULTS_ROOT / "Results"
# STUDY AREA / BARANGAY POLYGONS: polygon KML, WGS84 lon,lat, with name/code fields.
Study_Area = str(PROJECT_ROOT / "Boundary" / "INFANTA_BOUNDARY.kml")
BARANGAY_BOUNDARY = Study_Area
BARANGAY_NAME_FIELD = "ADM4_EN"
BARANGAY_CODE_FIELD = "ADM4_PCODE"
ENABLE_BARANGAY_OUTPUTS = True

PATHS = {
    **{f"Pre{i}": globals()[f"Pre{i}"] for i in range(1, 9)},
    "During": Event_Straddler or During,
    **{f"Post{i}": globals()[f"Post{i}"] for i in range(1, 9)},
}
# (master, slave), ISO YYYY-MM-DD or ISO timestamp; all UTC, timezone-naive.
# Blank entries are read from one YYYYMMDD_YYYYMMDD pair in the file name.
PAIR_DATES = {name: ("", "") for name in PATHS}
# Example: PAIR_DATES["Pre1"] = ("2022-05-04", "2022-05-16")
MIN_PRE_INPUTS = 3   # fewest PRE rasters accepted (8 recommended)
MIN_POST_INPUTS = 3  # fewest POST rasters accepted (8 recommended)

# Event window in UTC: calendar-day envelope for 25-26 September in Philippines.
EVENT_START = "2022-09-24T16:00:00"
EVENT_END = "2022-09-26T16:00:00"
EVENT_ORIGIN = "2022-09-25T09:30:00"  # PAGASA Burdeos landfall; proxy for local impact.
REFERENCE_DATE = "end"  # end recommended; midpoint is a sensitivity scenario.

# Must be confirmed by inspecting SLC/SNAP processing records, not inferred from TIFF.
ACQUISITION_COMPARABILITY_CONFIRMED = False
COMPARABILITY_NOTES = ""  # track, orbit direction, polarization, bursts, window, processing
REQUIRE_COMPARABILITY_CONFIRMATION = True  # run refuses to start until confirmed above
MAX_PAIR_BASELINE_RATIO = 1.25
ALLOW_BASELINE_MISMATCH = False  # requires substantive rationale in notes

# Raster reading
BAND = 1
RESAMPLING = "bilinear"
COHERENCE_SCALE = 1.0  # AFTER TIFF scale/offset; use 0.0001 only if metadata lacks scaling
COHERENCE_OFFSET = 0.0
WINDOW_SIZE = 256
NODATA = -9999.0
ZERO_IS_INVALID = False  # metadata masks still apply; zero is a valid magnitude
SOURCE_VALIDITY_CONFIRMED = False
SOURCE_VALIDITY_NOTES = ""  # explain declared NoData, zeros and any prior spatial filling
FIGURE_FORMATS = ("png", "pdf", "svg")
FIGURE_DPI = 300

# ============================================================
# ITSA DISTURBANCE AND RECOVERY CLASSIFICATION PARAMETERS
# ============================================================
RECOVERY_MODEL = "itsa"
CLASSIFICATION_MODE = "object_itsa"  # PRE-only objects; raw pixel products exported too
USE_ADAPTIVE_DISTURBANCE_THRESHOLD = True
EVENT_LOSS_MINIMUM = 0.03  # provisional effect-size floor; calibrate independently
PRE_RESIDUAL_NOISE_MULTIPLIER = 1.0
ITSA_NOISE_FLOOR = 0.02
EVENT_SUPPORT_MINIMUM = 0.03
MIN_RECOVERY_IMPROVEMENT = 0.03
EFFECT_DECISION_MARGIN = 0.005
RECOVERY_TOLERANCE = 0.05  # absolute coherence units, baseline minus tolerance
ITSA_MIN_PRE = 6   # valid PRE dates per series (capped at the number of PRE inputs)
ITSA_MIN_POST = 6  # valid POST dates per series (capped at the number of POST inputs)
LATE_POST_CONFIRMATIONS = 2
MIN_LATE_CONFIRM_SPAN_DAYS = 0.0
MAX_CONFIRM_GAP_DAYS = 60.0
ALLOW_OVERLAPPING_CONFIRMATION_PAIRS = False  # touching endpoints allowed, not independent
OBSERVATION_DISAGREEMENT_TOLERANCE = 0.10
MAX_ITSA_RMSE = 0.20  # hard model inadequacy screen, not a significance test
MAX_DESIGN_CONDITION = 1e6
UNCERTAINTY_MULTIPLIER = 2.0  # SE screen only, NOT a calibrated confidence interval

# Baseline observability. The historical 0.40 absolute screen left ~97% of
# Infanta Indeterminate. "signal_evidence" instead requires the PRE baseline to
# exceed tolerance + loss threshold + margin AND OBSERVABILITY_NOISE_MULTIPLIER x
# residual noise. The 0.40 rule is still exported as the "Original_*" products.
BASELINE_ELIGIBILITY_MODE = "signal_evidence"  # or "legacy_absolute"
OBSERVABILITY_NOISE_MULTIPLIER = 3.0
ITSA_BASELINE_MINIMUM = 0.40  # historical absolute observability screen

# Legacy settings still read by shared validation helpers.
BASELINE_MIN = ITSA_BASELINE_MINIMUM
BASELINE_MAD_MAX = 0.15
NOISE_FLOOR = 0.02
LOSS_ABSOLUTE = 0.10
LOSS_MAD_MULTIPLIER = 2.0
ENDSTATE_POST_PAIRS = 3
RECOVERY_CONFIRMATIONS = 2
MIN_PRE = 6
MIN_POST = 6
CONFIRMATIONS = 2
MIN_CONFIRM_SPAN_DAYS = 0.0
SAMPLE_PIXELS = 20000

# PRE-only, four-connected regions in physically sized, fixed analysis cells.
MIN_MAPPING_AREA_M2 = 4000.0
TARGET_REGION_AREA_M2 = 30000.0
MAX_REGION_AREA_M2 = 90000.0
MIN_REGION_AREA_M2 = MIN_MAPPING_AREA_M2
PRE_FEATURE_BIN_WIDTHS = (0.10, 0.05, 0.10)  # PRE median, MAD sigma, PRE fitted total change
MIN_PER_DATE_VALID_COVERAGE = 0.80
OBJECT_MIN_VALID_FRACTION = MIN_PER_DATE_VALID_COVERAGE
OBJECT_MIN_VALID_PIXELS = 3  # capped at the actual region member count

# Spatial propagation is disabled: no neighbor voting, filling or region growth.
ENABLE_SPATIAL_REGULARIZATION = False
ENABLE_SPATIAL_REGION_GROWTH = False
ENABLE_NEIGHBOR_CLASS_ASSIGNMENT = False
EXPORT_REGULARIZED_CLASSIFICATION = False
REGULARIZATION_MIN_NEIGHBORS = 2
REGULARIZATION_FEATURE_TOLERANCE = 0.05
EVIDENCE_PROTECTION_SCORE = 0.80  # normalized heuristic, NOT a probability
NEIGHBOR_RADIUS_M = 45.0
NEIGHBOR_MIN_AGREEMENT = 0.80
NEIGHBOR_SIGNATURE_TOLERANCE = 0.05
BUILTUP_MASK = r""  # optional verified exact-grid binary TIFF, polygon KML or GeoJSON
BUILTUP_MASK_NOTES = ""
BUILTUP_GEOJSON_CRS = "EPSG:4326"

# ============================================================
# MULTI-SCALE ZONAL RECOVERY (sarflow/zonal_recovery.py)
# ============================================================
ENABLE_ZONAL_RECOVERY = True
ZONAL_BASE_CELLS = 5            # 5 x 14 m pixels = ~70 m base zone
ZONAL_LEVELS = 4                # 70, 140, 280, 560 m nested zones
ZONAL_MIN_VALID_PIXELS = 9      # per date; smaller zones defer to the parent zone
ZONAL_MIN_VALID_FRACTION = 0.5  # per date, of the zone's study pixels
ZONAL_RESPECT_BARANGAYS = True  # zones never mix two barangays
ZONAL_SENSITIVITY_BASE_CELLS = (3, 7)

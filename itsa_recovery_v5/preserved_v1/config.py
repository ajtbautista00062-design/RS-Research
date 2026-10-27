"""Configuration for filename-driven, exploratory ITSA coherence processing."""
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
INPUT_DIRECTORY = PROJECT_ROOT / "Coherence"
Output_Folder = str(PROJECT_ROOT / "outputs")  # dated runs within the two output folders; audit folders also retained
Study_Area = r"C:\Users\Arlan\OneDrive\Desktop\Jhonjhon\Remote Sensing Research\Codes\Boundary\INFANTA_BOUNDARY.kml"
# Municipality boundary limits coverage; it does not identify built-up land cover.
FILL_SMALL_HOLES = True
FILL_MAX_HOLE_PIXELS = 9   # area limit; larger missing regions are preserved
FILL_RADIUS_PIXELS = 2     # maximum donor distance on the aligned reference grid


EVENT_START = "2022-09-24T16:00:00"
EVENT_END = "2022-09-26T16:00:00"
EVENT_ORIGIN = "2022-09-25T09:30:00"
REFERENCE_DATE = "end"  # end or midpoint; filename dates are acquisition intervals
REFERENCE_RASTER = r""  # empty: use event-straddling grid, with no canvas expansion
ACQUISITION_COMPARABILITY_CONFIRMED = False
COMPARABILITY_NOTES = "Orbit, polarization, bursts and coherence processing require verification."

BAND = 1
RESAMPLING = "bilinear"
COHERENCE_SCALE = 1.0  # applied after source metadata scale/offset
COHERENCE_OFFSET = 0.0
# Inherited footprint assumption, NOT established by nodata=None:
# True may exclude real zero coherence; False may include untagged background.
# Verify footprint/source processing records before changing it.
ZERO_IS_INVALID = True
ALL_TOUCHED = False
WINDOW_SIZE = 256
SAMPLE_PIXELS = 20000
RANDOM_SEED = 20220925
FIGURE_FORMATS = ("png", "pdf", "svg")
FIGURE_DPI = 600
NODATA = -9999.0

# Baseline recovery is enabled by the researcher's selected recovery definition.
# Numeric thresholds are operational defaults from the earlier discussion:
# validate sensitivity before thesis inference. Comparability is NOT assumed confirmed.
ENABLE_BASELINE_CLASSIFICATION = True
CLASS_MIN_PRE = 3
CLASS_MIN_POST = 3
LOSS_MIN = 0.10
LOSS_SPREAD_MULTIPLIER = 2.0
PRE_SPREAD_FLOOR = 0.02
RETURN_TOLERANCE = 0.05
END_POST_COUNT = 3          # require valid values for the last three chronological pairs
CONFIRM_POST_COUNT = 2      # both final pairs must reach baseline minus tolerance

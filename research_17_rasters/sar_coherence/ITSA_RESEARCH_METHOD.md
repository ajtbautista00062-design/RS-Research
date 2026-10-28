# Research-code ITSA recovery revision

## What inspection established

The research code in `Codes/sar_coherence` was not the simulation V5 code. Its
`metrics.calculate_fraction` had no regression. It assigned class 4 whenever an
eligible PRE median minus the event-pair coherence failed
`max(0.10, 2*max(PRE MAD sigma, 0.02))`. The baseline screen also required
coherence >=0.40 and PRE MAD sigma <=0.15. Classes 1–3 required all final three
POST values; class 4 did not. This asymmetric support and a fixed 0.10 loss
floor can restrict recovery-eligible coverage. Whether they explain any actual
map's dominance is **unmeasured**: configured paths are empty and no rasters were
scanned or processed for this revision. No claims about measured improvements
or the latest simulation map follow from this source inspection.

The original `metrics.py`, journal portfolio and observed workflow remain
available using `RECOVERY_MODEL="legacy_observed"`. The default `itsa` path also
exports that classifier as a clearly labeled comparison, not an ITSA result.
`Simulation` was not modified.

## Temporal model and recovery target

All 17 editable inputs remain: Pre1–Pre8, Event_Straddler (During alias), Post1–Post8.
Dates may be explicitly supplied as master/slave UTC pairs. When both are blank,
exactly one YYYYMMDD_YYYYMMDD filename pair is required. Filename order does not
prove master/slave roles: verify roles in the source products. Dates are sorted,
duplicate paths/pairs rejected, and phase support checked against EVENT_START/END.
Keep orbit, polarization, coherence estimation, temporal baselines and geometry
comparable. A filename alone cannot establish acquisition comparability.

Time `t` is actual elapsed days from EVENT_ORIGIN. Pure PRE and POST pairs fit
`Y(t)=b0+b1*t+b2*P+b3*t*P+error`. Numerical time scaling is undone in the
reported coefficients. All finite observations are used; defaults require at
least six PRE and six POST, not a truncated six-point subset. The event-straddling
pair supports or contradicts disturbance detection and is drawn on the graph;
it spans both regimes and is not assigned to either pure-phase regression.
The standard model has 12 residual degrees of freedom with 8+8 valid dates.
There is no automatic pilot execution or substitution for missing full-data inputs.

`B=median(all valid PRE coherence)` is the fixed recovery target.
`M(t)=(b0+b2)+(b1+b3)*t` is the fitted POST trajectory.
`L=-b2` is the level loss at time zero, **not** the first-POST baseline deficit.
`D(t)=B-M(t)` and `progress=M(last valid POST)-M(first valid POST)`.
`counterfactual=b0+b1*t` is disturbance context, not a replacement target.
`fraction=progress/D(first POST)` is exported, unclipped, only where credible
modeled disturbance exists and the denominator is positive. It differs from the
legacy event-endpoint fraction; do not compare those as interchangeable metrics.

Candidate crossing `t_cross=(B-tolerance-b0-b2)/(b1+b3)` is reported as days since
EVENT_ORIGIN only for class 1, a positive slope, and a crossing bracketed by valid
observed POST dates. An already-returned first POST can support class 1 but has
unresolved recovery timing; crossing stays NoData. This is a modeled threshold
crossing, not a measured building reconstruction date or precise event duration.

## Detection and class precedence

The operational loss cutoff is
`H=max(EVENT_LOSS_MINIMUM, multiplier*max(PRE detrended MAD sigma, noise floor))`.
Noise uses residuals around the fitted PRE trend, separating background trend
from residual variability. The initial floor 0.03, residual multiplier 1.0 and
noise floor 0.02 are **provisional operational settings**, not calibrated SAR
constants or significance thresholds. Existing 0.40 baseline observability is
retained. No floor was inferred from class coverage or independently validated.

1. Class 0: outside mask/barrier exclusion or no finite temporal observations.
2. Class 5: insufficient temporal support, low baseline, nonfinite/rank/condition
   failure, RMSE >0.20, POST predictions beyond [-0.05,1.05], major contradiction,
   near-threshold effect, or unresolved recovery. Reasons are exported separately.
3. Credible loss requires `L>=H+0.005` and either event-pair loss >=0.03 OR both
   modeled and observed early-POST baseline deficits >=0.03. A strong event-pair
   coherence gain contradicting modeled loss stays class 5. Indicators need not
   all pass the same threshold. Missing event data never creates evidence.
4. Class 4 requires `L<=H-0.005`, no supporting event loss, and both modeled and
   observed early deficits <0.025. Unsupported/contradictory level changes and
   margins remain class 5. Nonsignificance alone does not establish class 4.
5. Class 1: credible loss, both latest modeled values and both latest actual
   observations >=B-0.05, valid late temporal support, and no major late residual.
6. Class 2: credible loss, final modeled deficit >0.055, late observed deficit
   >0.045, modeled improvement >0.035, observed late-minus-early improvement
   >=0.025, and late model/observation agreement within 0.10.
7. Class 3: credible loss, the same supported persistent deficit, modeled
   improvement <=0.025, observed improvement <=0.035, and late agreement.
8. Remaining ambiguous recovery trajectories are class 5. A positive slope alone
   cannot create class 1. Late confirmation cannot bridge missing dates, prohibited
   overlapping acquisition intervals, or gaps exceeding MAX_CONFIRM_GAP_DAYS.

Steps describing class 5 are failure conditions; usable series then pass the
mutually exclusive class 4 / credible-loss branches and final recovery branches.
Defaults and margins are all configurable in `config.py`; calibrate before thesis
inference. Sensitivity includes floors 0.01/0.02/0.03/0.05/0.10 and return tolerances
0.03/0.05/0.07 on fixed object geometry. It is not a search for balanced colors.

OLS level, progress and final prediction SEs derive from the same design covariance.
SE screens and valid-date counts are quality indicators, not significance gates.
The multiplier 2 is not claimed as an exact confidence level. Errors may be
autocorrelated through shared SLCs; independence-based SEs may be optimistic.
Spatial aggregation increases spatial support, not temporal sample size. Separate
quality flags disclose uncertainty and unresolved source semantics. A classified
series with a quality warning is provisional; class 5 remains for unresolved evidence.

## Spatial objects and bounded processing

Object mode is the preferred map. Segmentation uses four-connected components
within fixed spatial cells, grouped by PRE median, PRE MAD sigma and PRE fitted
total change bins. **No POST or class values define the objects.** PRE-only trend
features are computed independently of POST missingness. Optional exact-grid
positive-integer ancillary labels prohibit crossing known boundaries; zero excludes.
A municipal KML is not a verified built-up mask or automatic river/coast barrier.

Initial min/target/max areas are 2,000/10,000/40,000 m² (0.2/1/4 ha), physical
scales to be reviewed against settlement/building sizes and native resolution.
At ~14 m pixels, target cells contain about 49 pixels; smaller heterogeneous
components are retained. Cells and maximum domains are rounded down to whole
pixel dimensions; domains are multiples of cell size. Actual area uses projected
CRS unit conversion. Components can be small; minimum area flags them rather than
erasing or expanding them. Maximum-domain edges can split real objects. Test
domain size, anchoring and feature-bin sensitivity before geographic interpretation.
This conservative grid-constrained region method is not SLIC or a validated building
segmentation. Adjacent same-class objects can already form contiguous map zones.

For every object/date, exact member-pixel median, valid count, fraction and MAD
are retained. Default support is >=80% valid members and >=min(3, actual member
count). Missing object dates remain missing. Aggregated observations get their
own ITSA fit with exactly the same class semantics. Insufficient PRE support for
segmentation leaves observable pixels indeterminate, not automatically no-loss.

Optional regularization merges adjacent **same-class**, low-confidence objects
only when baseline, loss, progress and final deficit differ by <=0.05, ancillary
labels agree, and pooled area is <=40,000 m². It reaggregates the actual member
observations and refits ITSA; the pooled result must retain the class. Small
high-confidence objects and classes 0/5 are never overwritten. Merging is bounded
to fixed physical domains, independent of I/O chunk size. It does not relabel
disagreeing objects by majority vote. Consequently categorical regularization
change may be identically zero: partition simplification is not improved accuracy.

Only one I/O window of 17 values is held at a time. Chunk lengths are multiples
of physical domains so objects do not split at arbitrary I/O seams; IDs are
reference-grid minimum member indices+1. There is no full-scene stack in RAM,
spatial gap filling, dilation, morphology-based coverage growth or post hoc class
frequency manipulation. Numerous small objects and sensitivity refits may still
make full runs slow; runtime/performance has not been measured.

## Source validity and output products

Metadata masks, per-band scales/offsets and exact reference-grid VRTs are retained.
Valid zero magnitude is kept unless independently documented as missing. Setting
ZERO_IS_INVALID requires source-validity notes and confirmation. If missing zeros
are undeclared, nearest resampling is required to avoid mixing them through a
bilinear warp. Existing declared NoData is still honored. No values are filled.
Unverified source semantics leave all operational results flagged provisional.

Each future run creates a timestamped subfolder below Output_Folder, with
`Time Series Graphs` and `Recovery Map Raster File`; originals and earlier runs
are not overwritten. Maps use the reference CRS, transform and dimensions,
DEFLATE/tiled compression and explicit NoData:

- ITSA_Recovery_Classes: chosen mode's primary Byte map, codes 1–5, NoData 0.
- Raw_Pixel_ITSA_Classes; Object_ITSA_Classes; optional Spatially_Regularized_ITSA_Classes.
- Legacy_Observed_Classes: original observed-only comparison, codes 1–4.
- Spatial_Object_ID and Regularized_Object_ID: UInt32, NoData 0.
- ITSA_Metrics: Float32 bands for baseline, level loss, event loss, POST slope,
  initial/final model, deficit, progress, fraction, supported crossing, threshold,
  fit errors/SEs, observation counts/df and all four coefficients. NoData -9999.
- ITSA_Quality_Flags: UInt16 bitmask, outside-mask 65535. Bits 1/2 uncertainty,
  4 below MMU, 8 incomplete temporal support, 16 source semantics unresolved.
- ITSA_Exclusion_Reasons: Byte reason codes in metadata, outside code 11;
  reserved NoData 255. Structural missingness and insufficient support are distinct.
- Pixel_Object_Classification_Change and Regularization_Classification_Change:
  0 unchanged / 1 changed / NoData 255 outside study mask.
- Region CSVs preserve IDs, coefficients, area, class/quality, all per-date medians,
  counts/fractions/MAD. Class-area tables include the entire inside-boundary
  denominator, including structural missingness. Object counts have object units.
  Areas represent included grid-cell areas, not exact fractional polygon coverage;
  small/narrow polygons require a boundary-rule sensitivity assessment.
- Gate accepted/rejected counts are independent gate totals, not sequential
  attrition; terminal rejection reasons are disjoint. Quality counts can overlap.
- Exact binned distributions and under/overflow counts; loss/return sensitivity
  counts and changed shares; legacy/pixel and pixel/object transition matrices.
- Zonal class counts/shares and eligible fraction summaries use the configured
  polygon pixel-center rule (or ALL_TOUCHED sensitivity). Class shares exclude
  code 0; class 5 is reported. Means/std/min/max are continuous fraction summaries
  among classes 1–3 only. No class-code mean is calculated. Quantiles require a
  later explicit extraction; pixels/objects are not independent temporal replicates.

The ITSA graph uses fixed common-pixel support, actual pair intervals, monochrome
symbols, separate fitted trends, counterfactual, fixed baseline, return tolerance
and event origin. Spatial IQR bars are **not** confidence intervals. Model error
SEs are available separately; graph fit is descriptive spatial-median analysis,
not the model used for every individual region. PNG 350 DPI and PDF/SVG are default.

## Execute later; validation status

Edit paths, pair dates (or valid filenames), Study_Area, Output_Folder and acquisition
comparability first. Existing empty placeholders intentionally prevent execution.
Keep all input paths inside this designated Codes workspace. Set source-validity
notes only from actual processing evidence. Install existing requirements in your
own environment if needed. From the project root:

```powershell
Set-Location -LiteralPath "C:\Users\Arlan\OneDrive\Desktop\Jhonjhon\Remote Sensing Research\Codes"
if ((Get-Location).Path -ne "C:\Users\Arlan\OneDrive\Desktop\Jhonjhon\Remote Sensing Research\Codes") { throw "Wrong project" }
python -B .\sar_coherence\run_analysis.py
```

29 synthetic tests are written in tests/test_itsa_research.py, covering controlled
recovery outcomes, OLS against independent segment fits, all eight PRE observations,
missing data, zero values, event exclusion, observation confirmation, crossing
support, PRE segmentation, barriers, global IDs, spatial support, MMU flagging and
high-confidence/indeterminate preservation. They were **not executed**. Only AST
syntax checks were performed. No simulation, full raster scan, raster processing,
graphs or analytical results were generated. Real-pixel coefficient checks,
runtime checks, acquisition comparability, coverage, threshold calibration and
independent physical-recovery validation remain outstanding.

After separate authorization, run tests before data processing; then use a new
output directory, independently verify selected real-pixel fits and inside-boundary
coverage, assess scale/bin/threshold sensitivity and validate against held-out
field/optical evidence. A 3+3 pilot remains exploratory and cannot establish full
17-observation accuracy. Code syntax is not numerical or scientific validation.

## QGIS

Open ITSA_Recovery_Classes.tif and use Paletted/Unique values: 1 teal #009E73,
2 amber #E6AB02, 3 vermilion #D55E00, 4 blue #0072B2, 5 gray #999999;
0 transparent. Embedded palettes supply colors. Use nearest-neighbor display
for categorical maps; bilinear interpolation invents class mixtures. For a smooth
continuous display, duplicate ITSA_Metrics and select a named progress/deficit/
fraction band using Singleband pseudocolor. Keep numeric values unclipped.
Use Zonal histogram for class counts and selected continuous bands for Zonal
statistics; do not average class codes. Municipal coverage is not built-up coverage.

## Method references and their scope

[Bernal et al., segmented ITS tutorial](https://pmc.ncbi.nlm.nih.gov/articles/PMC5407170/)
supports explicit impact-model specification and attention to confounding,
autocorrelation and seasonal structure; it does not validate this SAR threshold
scheme. [SciPy connected-component documentation](https://docs.scipy.org/doc/scipy/reference/generated/scipy.ndimage.label.html)
documents the connectivity operation. No cited source is claimed to endorse the
particular provisional recovery cutoffs, mapping areas or observed coherence as
proof of physical reconstruction.

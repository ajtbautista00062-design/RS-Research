# V5: operational ITSA baseline recovery with separate uncertainty

## Status

This is a source-only revision. No tests, simulation, raster scans, regression
experiments, plots or analytical exports were executed. Historical V3/V4 results
are not V5 results. V1-V4 sources, original inputs and previous outputs remain
unchanged. V5 retains pixel, PRE-defined object, and optional regularized modes.

## Source-level audit of excessive indeterminate decisions

The previous recovery.assess required guarded modeled loss, guarded return and
guarded improvement. With six regression observations and four parameters there
are two residual degrees of freedom. Wide t-based intervals could consequently
make an otherwise interpretable fitted trajectory indeterminate. This was a
confidence-bound gateway rather than an explicit p-value test, but it had a
similar operational effect. Insufficient numerical rank/df and genuine missing
support remain distinct from wide statistical uncertainty.

Other identified gates were: complete PRE/POST support even for longer series;
mandatory event pixels; missing late confirmations; RMSE/PRE variability screens;
contradictory observations; and converting a supported endpoint return to Class 5
merely because its crossing time was unobservable. At object level, requiring
three spatial pixels made smaller objects impossible to assess, any single-date
high dispersion forced rejection, and optional MMU handling deferred small
low-confidence classes automatically.

These are mechanisms identified in source, not newly measured causes or counts.
No claim about actual V5 class frequencies can be made before an authorized run.

## Shared fit and fixed baseline

recovery.assess performs one segmented fit and retains the historical guarded
decision. operational_recovery then evaluates the same fitted coefficients using
effect-size and observation rules. CLASSIFICATION_EVIDENCE_MODE selects
operational (default) or guarded. Both use actual pair-reference days from the
same event origin and the same spatial grid.

Y(t)=b0+b1*t+b2*P+b3*t*P+error.
B=median(valid PRE); M(t)=(b0+b2)+(b1+b3)*t.
D(t)=B-M(t); improvement=M(last POST)-M(first valid POST).
Modeled interruption loss L=-b2 differs from observed event loss B-E and the
early fixed-baseline deficit B-M(first valid POST). These differences are
explicitly checked rather than treating all three as the same measurement.

The event-straddling pair spans both phases. It remains visible and supplies
supporting/conflicting evidence, but is not assigned wholly to a regression
segment. Object medians are actual observations on fixed PRE-defined geometry.
Spatial pixel counts never increase the temporal sample size or residual df.

## Effect-size detection and observational support

The loss effect threshold is H=max(configured loss floor,0.5*PRE detrended noise).
The default floor uses PRE-only noise calibration inherited from V3, approximately
0.0459 for the earlier data; future runs recalibrate it. This proxy is provisional,
not a validated detection limit or false-alarm guarantee. It is not reduced to
create recovery pixels. Independent control areas/acquisition comparability must
validate pixel and object thresholds separately.

With margin m=0.005, modeled loss requires L>=H+m. It must be corroborated by
either event loss >= event-evidence tolerance OR both modeled and observed early
baseline deficits >=0.03. All indicators need not exceed the same threshold.
A clearly opposing event increase remains a contradiction. When event data are
missing, both early indicators are required and quality records the missing pair.
Set REQUIRE_EVENT_OBSERVATION=True to require the pair explicitly.

No-loss requires L<=H-m, no positive disturbance evidence, and both early
deficits below their threshold. A level-loss estimate close to H is Class 5,
not a fallback Class 4. A baseline deficit following a PRE decline without a
corresponding interruption step is unresolved, not automatically typhoon damage.

## Mutually exclusive outcome hierarchy

1. Outside mask or all observations structurally missing: Class 0.
2. Insufficient temporal support, missing required confirmation, unavailable fit,
   nonfinite/physically invalid predictions, hard model failure or unresolved
   all-zero semantics: Class 5 with an explicit reason.
3. Disturbance indicators contradict or the effect lies within its operational
   margin: Class 5. Supported absence under the rule above: Class 4.
4. Credible modeled disturbance enters recovery assessment:
   - Class 1: final two fitted values and actual observations reach B-0.05.
   - Class 2: fitted endpoint remains below B-0.05-m, modeled improvement exceeds
     0.03+m, and late observations support both a remaining deficit and improvement.
   - Class 3: same supported remaining deficit, modeled improvement <=0.03-m,
     with observed behavior consistent with negligible improvement or decline.
5. Unresolved improvement near its margin, a model/observation return conflict,
   or large late residuals: Class 5. Positive slope alone never assigns recovery.

Late below-baseline confirmation uses the late median and final observation,
with an explicit 0.02 observation tolerance. A transient earlier high value does
not automatically establish sustained recovery. MAX_OBSERVATION_DISAGREEMENT
defaults to 0.10 coherence units. RMSE>0.10 remains a quality warning; hard failure
uses a separately configurable 0.20. Physical prediction and rank/condition
screens remain. These operational choices require scientific calibration.

Minimum PRE/POST support remains three each and final consecutive confirmations
remain two. REQUIRE_COMPLETE_SERIES=False permits missing earlier observations
in longer series while using every available valid date. The seven-observation
pilot still needs all three PRE and three POST to meet the minimum. Missing late
confirmation is not bypassed by extrapolating a fitted endpoint.

## Recovery state versus recovery duration

A supported late return does not depend on estimating an exact crossing date.
Duration is exported only for an observed below-to-above bracket and a modeled
crossing within the observed POST interval. When baseline is already reached by
the first POST observation, Class 1 can describe that observed state, while
duration remains NoData and crossing_status identifies censored timing. This
does not assert that the crossing occurred inside the observed POST interval.
After-period crossings are candidates, not observed recovery durations.

## Statistical uncertainty is independent

The existing coefficient covariance, PRE noise floor, standard errors and t-based
guards are retained. Wide intervals flag provisional outcomes rather than
automatically rejecting them. evidence_status is 0 unresolved, 1 operational
outcome agrees with the guarded rule (conditional), or 2 provisional. Grade 1 is
not independent scientific confirmation; source, calibration and dependence
warnings remain in quality_flags. No calibrated class probabilities are claimed.

New separate products include valid PRE/POST counts, temporal fit n and df,
RMSE/condition, loss threshold, loss/progress guard widths, early deficit,
event agreement, observed confirmation count, crossing status and evidence status.
ITSA_Evidence_Status is Byte; continuous evidence bands retain explicit units.
Guarded_ITSA_Recovery_Classes preserves the stricter same-fit comparison.
input_provenance still distinguishes observations from spatial estimates. Primary
data are unfilled. Sensitivity exports record actual estimated-observation counts.

## Spatial objects and conservative generalization

PRE-only segmentation geometry is unchanged. Object aggregation retains per-date
median, valid counts/fractions and MAD. The minimum spatial count is capped at
actual membership, so a one-pixel object is not required to have three pixels.
The 80% valid fraction still applies. Missing temporal dates remain NaN.

OBJECT_HETEROGENEITY_POLICY=quality_only records high dispersion as provisional
regional evidence, not per-pixel reconstruction. strict preserves the historical
rejection for comparison. Strong contradictory pixel evidence still forces an
object to Class 5. Region medians cannot prove uniform member-pixel recovery.

MMU_POLICY=flag_only retains a supported small class and marks it provisional;
defer_low_confidence reproduces the optional earlier deferral. Regularization
still requires adjacent PRE-similar same-class objects and a confirming pooled
ITSA fit. It never promotes Class 5 or changes a protected high-confidence class.
Independent masks/barriers remain optional exact-grid inputs; the municipal KML
alone is not verified built-up cover. Unknown boundaries cannot be guaranteed.

## Future diagnostics, not generated during this task

Pixel reports contain class frequencies, guarded comparisons, gate/rejection
counts, confidence categories and sampled loss/slope/deficit/improvement
distributions. Object reports include class counts by objects and member pixels,
spatial rejection reasons, dispersion flags, indicator distributions, conditional
guard support/provisional counts, and fixed-object loss/return threshold sweeps.
The primary map count denominator includes the analysis domain and explicitly
reports Class 0. Object member counts use a separate labeled denominator.
pixel_object_agreement.csv records the full class contingency table.

Object sensitivity holds PRE geometry and the primary pixel conflict reference
fixed while refitting actual regional observations. It is not an independent
accuracy assessment. Pixel threshold sensitivity and object sensitivity are
reported separately. ENABLE_THRESHOLD_SENSITIVITY_REPORT controls future sweeps.
Do not average classification/evidence codes in zonal statistics; numeric
continuous fields are separated from codes by their units.

## Source and validation status

Changed modules: recovery/config, new operational_support/operational_recovery,
object_statistics/spatial_regularization, diagnostics/reporting interfaces,
raster_outputs/object_products, object_pipeline and run_analysis. V1-V4 are
preserved. New synthetic specifications cover nonsignificant yet interpretable
outcomes, effect-margin ambiguity, contradictions, missing event/late data,
state-versus-timing, incomplete longer series, tiny object support, heterogeneity
and MMU policies. Historical guarded tests select guarded mode explicitly.

Tests were written but NOT executed. Static syntax review is not numerical or
scientific validation. No V5 rasters or graphs exist from this task. Resolve
source zero semantics, pair comparability, independent recovery evidence,
spatial mixing, temporal dependence and threshold calibration before authoritative
inference. Full 17-observation execution remains explicitly gated.

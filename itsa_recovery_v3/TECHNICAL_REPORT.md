# ITSA baseline recovery V3: measured investigation and implementation

## Scope and preservation

This revision is confined to Simulation/itsa_recovery_v3 and its new diagnostic
outputs and launch script. Original sar_coherence and itsa_recovery_v2 sources,
input coherence rasters and previous results are preserved. Source hash manifests
provide the preservation check. No dependency changes were made.

Diagnostics use the saved unfilled seven-band stack from the completed V2 run
itsa_v2_pilot_20220721_20221106_20261027T144056520776Z. This is an assessment of
existing real observations, not a new original-source processing run or a
17-observation analysis. No synthetic data were substituted for real inputs.
Synthetic sequences are used only by software tests.

## What actually dominated the previous map

The configured KML boundary contains 876,543 pixel centers on the output grid.
Independent raster reads and reproduction of the old classifier agreed exactly
(zero mismatches). Class 5, not Class 4, dominates the whole boundary. Class 4
dominates the smaller set of assessed Classes 1-4. Display appearance is not a
class-frequency measurement.

| Class | V2 pixels | Revised provisional pixels | Revised boundary share |
|---|---:|---:|---:|
| 0 NoData | 0 | 0 | 0.0000% |
| 1 Recovered | 157 | 1,242 | 0.1417% |
| 2 Partial | 106 | 72 | 0.0082% |
| 3 Not recovered / degraded | 158 | 246 | 0.0281% |
| 4 No detectable loss | 199,348 | 57,955 | 6.6118% |
| 5 Indeterminate | 676,774 | 817,028 | 93.2103% |

V3 changes 153,895 classifications relative to V2. This is a methodological
transition, not evidence that the new labels are independently correct.

## Measured causes and rule audit

Counts below use every boundary pixel as denominator. Conditions overlap; they
are not an additive waterfall. Accepted plus rejected equals 876,543 per row.

| V2 condition | Accepted | Rejected |
|---|---:|---:|
| Finite three PRE / event / three POST, each | 876,543 | 0 |
| Stable eligible fit | 848,586 | 27,957 |
| Raw event loss above old threshold | 96,963 | 779,580 |
| Raw modeled level loss above old threshold | 159,061 | 717,482 |
| Guarded event loss above old threshold | 53,220 | 823,323 |
| Guarded modeled loss above old threshold | 20,719 | 855,824 |
| Both guarded losses | 2,327 | 874,216 |
| Both upper loss bounds below threshold (Class 4) | 199,348 | 677,195 |
| Absolute 0.10 floor controls threshold | 544,931 | 331,612 |
| PRE spread controls threshold | 331,612 | 544,931 |
| All seven observations exactly zero | 10,496 | 866,047 |

The dual guarded-loss gate restricted disturbance eligibility. However, changing
the absolute floor alone cannot resolve weak temporal evidence. In the fixed-seed
spatial sample, median fitted level loss is approximately -0.0129, median observed
event loss 0.0123, and median old level-loss guard 0.2427 coherence units. Some
pixels have little loss and others have very uncertain loss. Neither is proof
that the typhoon caused no physical damage.

All cached primary inputs are finite within the boundary, so this particular
comparison is not explained by missing PRE/event/POST pixels. Source zero
semantics and acquisition comparability remain unresolved. The boundary is not
independently validated as a built-up-only mask.

## Model, detection and recovery rules

Time t is actual days from the configured event origin, with the configured
pair-end reference. Fit Y=b0+b1*t+b2*P+b3*t*P+error on PRE and POST. The
event-straddling pair is displayed and assessed separately because its temporal
support spans both phases; assigning it wholly to either segment is misleading.
Three PRE plus three POST leave only two residual degrees of freedom.

The fixed target is B=median(observed PRE). Fitted POST M(t)=(b0+b2)+(b1+b3)*t.
Deficit D(t)=B-M(t); progress=D(first POST)-D(last POST); level loss L=-b2.
The target is the observed PRE condition, not the extrapolated PRE counterfactual.

The provisional loss floor H=0.04588735047 comes from the median PRE-only
detrended residual noise proxy in a fixed-seed sample of 50,000 pixels from
865,463 PRE controls. It is bounded below by ERROR_FLOOR=0.005. Calibration uses
neither POST outcomes nor class balance. Three PRE observations give one PRE-fit
residual degree of freedom. This noise proxy is NOT a validated detection limit
or false-positive rate; independent controls are still required.

The fitted coefficient variance is at least PRE-only residual variance and the
error floor variance. With a one-sided t multiplier at alpha=0.05 and
gL=tcrit*SE(b2), credible modeled loss requires L-gL>=H; bounded absence requires
L+gL<H. The interval between these decisions is indeterminate. These operational
guards are not calibrated joint confidence intervals and assume simplified OLS
errors despite shared acquisitions and possible temporal dependence.

Event support uses a separate fixed tolerance He=0.04588735047 and
gB=max(0.005,2*S/sqrt(nPRE)), where S=1.4826*MAD(PRE). Event opposition is
(B-E)+gB<-He. Positive support is (B-E)-gB>=He. A neutral pair does not veto
credible modeled loss. Against modeled absence, event loss exceeding max(H,He)
after its guard is conflicting evidence. Event tolerance stays fixed in sweeps.

- Class 0: outside boundary or all observations structurally unavailable.
- Class 5: insufficient PRE/event/POST support, unstable/unavailable model,
  unresolved all-zero sequence, uncertain loss, conflicting indicators or return
  unsupported within the observed POST interval.
- Class 4: sufficiently supported, bounded modeled absence with no contrary
  event evidence. A nonsignificant coefficient is never sufficient alone.
- Class 1: credible modeled loss; final two fitted POST values after guards and
  final two actual observations reach B-0.05; crossing is supported inside POST.
- Class 2: credible loss; guarded final prediction and final two observations
  remain below B-0.05; lower guarded modeled progress exceeds 0.03.
- Class 3: same supported degradation, but upper guarded progress is at most
  0.03. Positive slope alone cannot produce a recovered label.

All remaining cases are Class 5. Primary support is complete PRE/POST, at least
three of each. RMSE>0.10, PRE spread>0.15 or POST predictions outside [-0.05,1.05]
are instability screens. Thresholds are configurable operational assumptions.
Candidate crossing=(B-0.05-b0-b2)/(b1+b3) for positive POST slope. Duration is
exported only for an observed below-to-above bracket within the POST interval,
in days from event origin. No extrapolated duration is reported as observed.

## Final separated-threshold sweep

| Modeled loss floor | Recovered | Partial | Degraded | No loss | Indeterminate |
|---|---:|---:|---:|---:|---:|
| 0.01 | 1,518 | 72 | 410 | 41,868 | 832,675 |
| 0.02 | 1,456 | 72 | 360 | 45,798 | 828,857 |
| 0.03 | 1,373 | 72 | 309 | 50,093 | 824,696 |
| 0.05 | 1,219 | 72 | 235 | 60,934 | 814,083 |
| 0.10 | 762 | 67 | 131 | 100,804 | 774,779 |

59,949 pixels (6.8393% of the boundary) change category somewhere across these
candidates. A higher loss floor can increase bounded no-loss decisions; the
class counts need not vary monotonically as recovery counts. The final sweep is
V3_Final_Assessment; V3_Assessment is an earlier coupled-tolerance experiment.
No ground truth was supplied, so agreement and false-detection rates cannot be
estimated. No threshold was selected to create a balanced map.

## Spatial sensitivity and validation

Primary results are unfilled. Reusing the saved V2 corrected branch and its
provenance isolates filling under the zero-excluded policy. Filling changes
280 categories: 6 NoData to indeterminate, 8 indeterminate to recovered,
3 to degraded and 263 to no loss. This demonstrates that filling affects
analytical classifications; it remains a separately labeled sensitivity.
Zero-excluded unfilled has 12,833 NoData pixels versus zero in the retained-zero
primary. This does not resolve which source convention is correct.

All 34 software tests pass, including complete/partial/no recovery, degradation,
neutral/conflicting event evidence, uncertainty, model dependence, 17-observation
support, zero handling, graph rendering and synthetic GeoTIFF round trip.
100 real-pixel regressions and predictions agree with independent calculations.
Final class counts were independently reread from the exported raster. The
archive-and-run script passed PowerShell parsing but was not executed.

## Products and remaining limitations

Class4_Diagnostics_v3_20261027 contains exact V2 condition counts, sampled
distributions, PNG/PDF/SVG diagnostics and V3_Final_Assessment. The latter holds
the class raster, fixed PRE baseline, fitted POST endpoints, level/event loss,
POST slope, final deficit, progress, supported duration/brackets, detection
boundary, uncertainty guards, quality and exclusion rasters, threshold maps and
transition tables. They retain the reference CRS/grid and explicit NoData.

Quality flags identify unresolved zeros, unverified comparability, sparse support,
unvalidated calibration, conflicts and spatial estimates. Class 5 is valid
uncertainty, not missing-data transparency. Class 0 is transparent NoData.
Do not average class codes. Summarize class counts/proportions and continuous
metrics separately, with spatial dependence and valid counts reported.

Before authoritative 17-observation analysis: verify source masks/zeros, orbit,
polarization, temporal/perpendicular baselines and processing comparability;
validate detection and recovery criteria against independent controls/evidence;
review event origin and pair-reference sensitivity; then complete the V3
original-source pilot and its checks. Coherence recovery alone does not establish
reconstruction. The full-data gate remains explicit and disabled by default.

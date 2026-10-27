# V3 decision rules


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


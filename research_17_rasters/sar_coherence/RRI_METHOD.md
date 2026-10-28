# Independent ITSA final-observation RRI

This optional post-processing command reads a completed 17-observation research run. Existing classifiers, input configuration, region growth, graphs and barangay products are unchanged. It does not refit any pixel or barangay ITSA model. Five Python files were added: `rri_config.py`, `run_rri.py`, `sarflow/rri.py`, `sarflow/rri_pipeline.py`, and `tests/test_rri.py`, plus this Markdown documentation. No existing Python file was modified for RRI.

## Definition and event support

Time is measured in actual days centered on the saved intervention origin. With the existing segmented coefficients:

- Cbaseline = beta0, fitted PRE counterfactual at t=0.
- Cevent = beta0 + beta2, modeled immediate POST coherence at t=0.
- Cfinal = beta0 + beta2 + (beta1+beta3)*t_final, using the final chronological POST reference date in the saved manifest.
- Disturbance = -beta2; recovery magnitude = Cfinal-Cevent.
- RRI = (Cfinal-Cevent)/(-beta2).

This is a **final-observation adaptation**, not conventional maximum-post-recovery RRI. It differs from the existing progress/first-POST-deficit fraction and from the existing fixed observed PRE-median recovery target. Coherence recovery does not establish physical reconstruction.

The event-straddling pair measures change over an interval. It is not an instantaneous measurement at t=0. As in the existing ITSA implementation, it supports/contradicts disturbance detection and appears in the recorded 17 observation dates, while pure-phase OLS uses the eight PRE/eight POST observations. The immediate POST intercept is an extrapolated estimate at the intervention; it may be uncertain if observations are distant. Neither the observed event pair nor enhanced class labels replace this modeled estimate in RRI.

## Eligibility and codebooks

Inputs are the saved **raw pixel** ITSA coefficients, raw class and reason rasters, adaptive threshold, original selected-mode classification and coverage stack. Raw classes 1–3 establish the existing operational disturbance and temporal-confirmation requirements. Class 4 is explicitly undisturbed; unresolved Class 5 is not automatically undisturbed. Requiring a resolved existing raw class is conservative and can exclude disturbed but unresolved trajectories. RRI values are not derived from class codes: class codes only gate evidence; numerical estimates come from coefficients.

The denominator must exceed both MIN_DENOMINATOR (0.01 coherence units) and the saved pixel adaptive disturbance threshold. Modeled baseline, immediate POST and final coherence must be finite and within [0,1]. Missing/structurally unavailable pixels remain NoData. No minimum mapping area or neighbor filling is applied. The score-based spatial protections are preserved by leaving spatial classification untouched and by refusing to generate RRI from enhanced labels.

Continuous RRI, disturbance and recovery rasters: Float32, NoData -9999, exact source grid and CRS, tiled/DEFLATE. Values below zero and above one are preserved. Interpretation/status: Byte, NoData 0.

Interpretation codes: 1 further degradation (<-0.05); 2 no measurable recovery (within 0.05 of zero); 3 partial recovery (above 0.05 and below 0.95); 4 modeled-baseline return (within 0.05 of one); 5 exceeding baseline (>1.05). Floating-point roundoff is handled at tolerance boundaries. These dimensionless tolerances are configurable in `rri_config.py`, provisional and subject to sensitivity analysis. They are not confidence intervals or formal significance tests.

Status codes: 0 outside/structurally unavailable; 1 valid; 2 no detected disturbance; 3 insufficient observations; 4 invalid/near-zero denominator; 5 invalid/nonfinite/out-of-range fitted result; 6 unresolved/contradictory temporal evidence. The extra status 6 avoids misrepresenting ambiguous evidence as absence of disturbance.

## Outputs and statistics

All products go exclusively to `Codes/RII/Raster Files` and `Codes/RII/CSV`, with `Barangay/<actual name and code>` children. The command refuses nonempty output roots to preserve previous results. Archive previous RII results before another run; no automatic deletion is performed.

Five requested study rasters and per-barangay clips are generated. CSV outputs include study statistics, category distribution, statistics by original selected-mode recovery class, processing statuses, barangay summary, and a classification-agreement table. Population standard deviation uses divisor N. Exact medians/quartiles use temporary disk-backed values; no synthetic observations are used. Valid area uses actual projected pixel dimensions converted to square meters/hectares. Percentage RRI >=1 uses the exact continuous value, independently of the tolerance category.

The boundary was previously inspected: 36 polygon features with names `ADM4_EN`, codes `ADM4_PCODE`. Runtime polygon validation and code-sorted pixel-center ownership reuse the existing boundary functions. No arbitrary barangays are created. Each barangay reuses the same RRI raster; its model is not refitted. Every acquisition must have at least 80% valid coverage within owned raster centers for an eligible regional summary. Below this support, mean/median/range/quartiles/std and >=1 percentage are suppressed; counts, areas, status and coverage diagnostics remain visible. The denominator excludes off-raster portions of the polygon; inspect incomplete geographic coverage separately.

Disagreement tables compare RRI interpretation against original Classes 1–3 and explicitly flag differing outcomes. Differences are expected because the existing classifier uses the observed PRE median, sustained late confirmation and trajectory evidence, whereas RRI uses fitted intervention counterfactual and the final modeled observation. Disagreement is not itself a validation error.

## Run later

Use a working project environment containing all listed geospatial dependencies. A completed enhanced research run must already contain `Raw_ITSA_Metrics.tif`, raw classes/reasons, coverage and original-class rasters. Its matching `CSV files/<run_id>/run_manifest.json` supplies dates/settings. Older runs missing these artifacts must first be generated by the existing main workflow.

`python -B .\sar_coherence\run_rri.py --source ".\Rasters\<completed_run_id>"`

Optionally pass `--manifest` with the matching manifest path. Never mix coefficients and a manifest from different runs. No full processing was executed during implementation. Small numerical tests and AST syntax checks were run; GIS raster/CSV integration remains unverified in the current environment because its geospatial dependencies are unavailable and its SciPy installation is broken.

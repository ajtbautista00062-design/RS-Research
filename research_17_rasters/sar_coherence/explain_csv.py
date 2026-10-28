"""Generate CSV_Explained.md from output filenames and actual CSV headers only.

Standard library only. Does not import analysis code, read raster pixels, change
CSV files or execute ITSA/RRI. Rerun after generating new analysis outputs.
"""
import argparse
import csv
from pathlib import Path
import re

FILES = {
 'class_area_statistics.csv':'Class counts, physical areas and study-mask percentages for each analytical product.',
 'recovery_class_area_statistics.csv':'Original/enhanced class areas, separated by external built-up/non-built-up/unknown strata.',
 'raw_vs_enhanced_comparison.csv':'Original versus enhanced class-area comparison; currently repeats the class-area table with the same schema.',
 'class_transition_matrix.csv':'Counts and areas for each original-to-enhanced class combination.',
 'pixel_reassignments.csv':'One record per spatially reassigned pixel; documents its own ITSA-supported change.',
 'spatial_region_statistics.csv':'Four-connected enhanced class components, including components spanning multiple analytical objects.',
 'connected_region_counts.csv':'Component counts by recovery class and mapping-area flags.',
 'analysis_object_statistics.csv':'Areas of the original analysis objects intersecting this barangay.',
 'evidence_score_summary.csv':'Summary of normalized evidence scores; these are not calibrated probabilities.',
 'coherence_time_series.csv':'Per-date available-pixel coherence statistics for the named barangay; spatial support can vary by date.',
 'valid_data_coverage.csv':'Valid pixel counts/fractions for each acquisition pair.',
 'itsa_regression_results.csv':'Descriptive aggregate ITSA model indicators. This aggregate fit is distinct from the pixel/object classification models.',
 'itsa_indicator_statistics.csv':'Descriptive statistics for each named modeled indicator within a barangay.',
 'itsa_plot_observations.csv':'Median coherence and spatial interquartile range from the common-pixel graph sample.',
 'raw_region_statistics.csv':'PRE-segmented analytical objects, their regional observations, model indicators and classes.',
 'regularized_region_statistics.csv':'Pooled/refitted analytical objects after optional same-class spatial regularization.',
 'object_class_counts.csv':'Number of analytical objects in each recovery category and processing mode.',
 'classification_gates.csv':'Accepted/rejected pixel counts for each named classification condition; gates are not mutually exclusive.',
 'operational_evidence_counts.csv':'Counts of provisional decisions, decisions without quality warnings and indeterminate outcomes.',
 'exclusion_reasons.csv':'Pixel counts for each classification exclusion/rejection reason.',
 'quality_flags.csv':'Counts for individual quality bits; a pixel can contribute to several bits.',
 'threshold_sensitivity.csv':'Class counts and changes under alternative loss floors and return tolerances on fixed object geometry.',
 'pixel_to_object.csv':'Transition matrix: row = original pixel class, numeric column = object class; entries are pixel counts.',
 'legacy_to_pixel.csv':'Transition matrix: row = legacy observed class, numeric column = ITSA pixel class; entries are pixel counts.',
 'zonal_recovery_statistics.csv':'Zone-level class shares and continuous recovery-fraction summaries among loss-eligible classes.',
 'raster_inventory.csv':'Source raster/grid metadata used for alignment and input auditing.',
 'barangay_export_status.csv':'Status of each named barangay export, including skipped/off-raster cases.',
 'rri_study_area_statistics.csv':'Statistics of valid continuous final-observation ITSA RRI over the study area.',
 'rri_statistics.csv':'Barangay RRI statistics; unsupported region summaries are suppressed when per-date coverage is below 80%.',
 'rri_class_distribution.csv':'Counts, areas and valid-RRI-area percentages for the separate RRI interpretation categories.',
 'rri_by_recovery_class.csv':'Continuous RRI summaries grouped by original selected-mode ITSA recovery class.',
 'rri_barangay_summary.csv':'Combined RRI summary and coverage eligibility for each exported barangay.',
 'rri_processing_diagnostics.csv':'Pixel counts by RRI validity status, including missing support and invalid model/denominator cases.',
 'rri_classification_agreement.csv':'Cross-tabulation of existing recovery classes and RRI interpretation; differences are not accuracy errors.',
 'rri_date_coverage.csv':'Per-acquisition valid-data coverage used to assess barangay RRI summary support.',
}

COLUMNS = {
 'zone':'Study area or named analytical zone/barangay.', 'barangay':'Original barangay name from the boundary attributes.',
 'barangay_code':'Administrative barangay identifier from ADM4_PCODE.',
 'name':'Observation or zone name; interpret with the filename and other identifying columns.',
 'code':'Recovery class code, except in barangay time-series tables where it is the administrative barangay identifier.',
 'label':'Text description of the associated class/reason code.', 'category':'Separate RRI interpretation code (see codebook below).',
 'interpretation':'RRI category label, or summary-support description in RRI barangay summary tables.',
 'recovery_class':'Existing ITSA recovery category; not an RRI interpretation code.',
 'original_recovery_class':'Existing selected-mode ITSA recovery class used for RRI comparison.',
 'rri_interpretation':'RRI interpretation category compared with the existing class.',
 'original':'Class before spatial refinement.', 'enhanced':'Class after spatial refinement.',
 'original_class':'Class before spatial refinement.', 'enhanced_class':'Class after spatial refinement.',
 'class':'Recovery category of the connected component.',
 'product':'Analytical product (pixel/object/regularized or original/enhanced).',
 'mode':'Analytical processing mode.', 'stage':'Object-processing stage (object or regularized).',
 'stratum':'External land-cover grouping: all, builtup, nonbuiltup or unknown.',
 'pixels':'Number of pixels in the row category or condition.', 'objects':'Number of analytical objects.',
 'count':'Number of valid numerical observations entering the summary.',
 'area_m2':'Pixel/component area in square meters, calculated from the projected grid.',
 'area_ha':'Area in hectares (square meters divided by 10,000).',
 'valid_area_m2':'Area represented by valid indicator pixels, in square meters.',
 'valid_area_ha':'Area represented by valid indicator pixels, in hectares.',
 'percent_valid':'Percentage of valid pixels in the relevant product/stratum (ITSA classes 1–5 or valid RRI pixels, according to filename).',
 'percent_of_study_pixels':'Percentage of all study-mask pixel centers, including structural NoData.',
 'percent_study':'Percentage of nonzero-status RRI study pixels; undefined for outside/structurally unavailable status.',
 'percent_rri_ge_1':'Percentage of valid RRI pixels whose continuous value is at least 1; not the tolerance-class percentage.',
 'denominator':'Explicit population used for the reported proportion.',
 'study_pixels':'Study-mask pixel count; in RRI diagnostics it excludes status 0.',
 'valid_pixels':'Valid pixels for the observation/zone; consult the file-specific definition.',
 'owned_pixels':'Pixel centers allocated to this barangay within the reference raster extent.',
 'total_owned_pixels':'All pixel centers entering this RRI summary, including status 0.',
 'valid_fraction':'Fraction (0–1) of owned/study pixels valid for the acquisition.',
 'minimum_date_valid_fraction':'Lowest valid-pixel fraction across the 17 acquisition dates.',
 'eligible':'Whether the acquisition meets the 80% valid-coverage requirement for the regional aggregate fit.',
 'coverage_eligible':'Whether every required date meets the 80% coverage rule for a supported regional RRI summary.',
 'observation':'Configured acquisition-pair identifier.', 'observation_index':'Zero-based index in the saved chronological 17-pair sequence.',
 'phase':'PRE, event-straddling, or POST analytical phase.',
 'master':'First/master acquisition date as recorded; filename order alone does not independently prove SLC roles.',
 'slave':'Second/slave acquisition date as recorded.', 'reference':'Representative pair timestamp used by ITSA (saved end or midpoint convention).',
 'date':'Acquisition-pair reference date used for the plot.',
 'support':'Description of the spatial sample/cohort contributing observations.',
 'scope':'Description of the aggregate model/sample being summarized.',
 'model_scope':'Distinguishes a separate regional aggregate fit from pixel/object classification models.',
 'metric':'Indicator name whose values are summarized in this row; defines the units of statistical columns.',
 'region_id':'Stable identifier of the PRE-based analytical object.',
 'component_id':'Identifier of a four-connected class component; distinct from an analytical object ID.',
 'pixels_in_barangay':'Number of member pixels of this object allocated to the barangay.',
 'area_in_barangay_m2':'Area of the object intersection with the barangay, in square meters.',
 'connected_regions':'Number of four-connected class components.',
 'below_mmu':'Whether/count of components below the 4,000 m² minimum mapping area; no automatic deletion.',
 'near_target_20pct':'Whether component area is within ±20% of the configured target region area.',
 'near_target':'Count of components within ±20% of target area.',
 'exceeds_object_max':'Whether a connected class component exceeds the analysis-object area limit; such a component may be legitimate.',
 'above_object_max':'Count of connected components above the analysis-object area limit.',
 'row':'Zero-based raster row of the reassigned pixel.', 'column':'Zero-based raster column of the reassigned pixel.',
 'evidence_score':'Accepted spatial-reassignment evidence score, 0–1; not a probability.',
 'median_rounded_0_001':'Median evidence score approximated using a histogram rounded to 0.001.',
 'probability_calibrated':'Whether the evidence score has probability calibration; currently false.',
 'quality_flags':'Bitwise quality-warning mask; several warnings can be present simultaneously.',
 'quality_bit':'Individual warning bit being counted.',
 'reason':'Exclusion code in model/region reports; explanatory text in Pixel_Reassignments.csv.',
 'threshold_adjusted':'Whether a reassignment changed the detection threshold; currently false.',
 'gate':'Named classification condition; accepted and rejected counts are reported separately.',
 'accepted':'True means the named gate passed; false means it failed.',
 'evidence':'Operational evidence grouping: provisional, no quality warning, or indeterminate.',
 'loss_floor':'Candidate minimum meaningful coherence loss used in a sensitivity scenario (coherence units).',
 'return_tolerance':'Candidate absolute baseline-return tolerance (coherence units).',
 'changed_share':'Fraction (0–1) of fixed-cohort object member pixels changing class under the scenario.',
 'lower':'Lower histogram-bin edge in the named indicator units; may be negative infinity.',
 'upper':'Upper histogram-bin edge in the named indicator units; may be positive infinity.',
 'status':'RRI validity code in RRI diagnostics; export status text in Barangay_Export_Status.csv.',
 'meaning':'Description of the RRI validity/status code.', 'folder':'Sanitized output directory name.',
 'disagreement':'Flag for differing ITSA/RRI interpretations under the documented comparison convention; not proof of error.',
 'note':'Context or interpretation limitation associated with the row.',
 'baseline':'Median of valid PRE coherence observations (fixed recovery target).',
 'modeled_level_loss':'Negative fitted ITSA level change, -beta2, in coherence units.',
 'observed_event_loss':'Observed PRE median minus event-straddling coherence (coherence units).',
 'post_slope':'Fitted POST slope beta1+beta3 (coherence/day).',
 'modeled_initial':'Fitted coherence at the first supported POST reference date.',
 'modeled_final':'Fitted coherence at the last supported POST reference date.',
 'baseline_deficit':'Fixed observed PRE median minus modeled final POST coherence.',
 'modeled_progress':'Modeled final minus modeled initial POST coherence.',
 'recovery_fraction':'Existing modeled progress / modeled first-POST baseline deficit; dimensionless and unclipped; different from independent RRI.',
 'crossing_days':'Supported modeled baseline-tolerance crossing time, in days since the intervention origin.',
 'adaptive_loss_threshold':'Pixel/object adaptive minimum disturbance threshold, in coherence units.',
 'rmse':'Regression root mean square residual (coherence units).',
 'level_se':'Conditional OLS standard error of beta2 (coherence units).',
 'progress_se':'Conditional OLS standard error of modeled POST improvement (coherence units).',
 'final_prediction_se':'Conditional OLS standard error of fitted final POST mean coherence.',
 'pre_noise':'Robust PRE regression residual spread (MAD scaled by 1.4826), in coherence units.',
 'pre_count':'Number of valid PRE observations used at this analytical location.',
 'post_count':'Number of valid POST observations used at this analytical location.',
 'residual_df':'Residual regression degrees of freedom: valid pure-phase observations minus four parameters.',
 'design_condition':'Condition number of the time-scaled regression design; larger values indicate poorer numerical stability.',
 'beta0':'Fitted PRE coherence at intervention t=0 (counterfactual intercept).',
 'beta1':'PRE trend (coherence/day).', 'beta2':'POST-minus-PRE fitted level change at intervention (coherence units).',
 'beta3':'POST-minus-PRE slope change (coherence/day).',
 'path':'Input raster filesystem path.', 'crs':'Raster coordinate reference system.',
 'resolution':'Raster pixel spacing in CRS units.', 'width':'Raster width in pixels.', 'height':'Raster height in pixels.',
 'bounds':'Raster extent in CRS coordinates.', 'transform':'Affine pixel-to-map transformation.',
 'nodata':'Declared raster NoData value.', 'dtype':'Stored raster numerical type.',
 'scale':'Source raster scale factor.', 'offset':'Source raster additive offset.',
 'grid_matches_reference':'Whether source CRS, affine transform and dimensions match the reference.',
 'eligible_pixels':'Pixel count in loss-eligible ITSA classes 1–3.',
 'fraction_valid_pixels':'Number of finite fraction values entering zonal fraction statistics.',
 'recovered_among_loss_eligible':'Class 1 count / combined counts of classes 1–3 (fraction 0–1).',
}


def column_description(column, filename):
    name=filename.lower()
    if column.isdigit() and name in ('pixel_to_object.csv','legacy_to_pixel.csv'):
        return f'Pixel count entering destination recovery class {column}; source class is the row index.'
    if not column or column.startswith('Unnamed:'):
        return 'Exported row index; in a transition matrix this is the source recovery class.'
    if column in COLUMNS:
        return COLUMNS[column]
    statistics={'mean':'Arithmetic mean','median':'Median (50th percentile)','q25':'25th percentile',
                'q75':'75th percentile','minimum':'Minimum','maximum':'Maximum',
                'std_population':'Population standard deviation (divisor N)',
                'spatial_q25':'Spatial 25th percentile','spatial_q75':'Spatial 75th percentile'}
    if column in statistics:
        quantity='continuous RRI (dimensionless)' if name.startswith('rri_') else (
            'evidence score (dimensionless)' if name=='evidence_score_summary.csv' else
            'the indicator named in metric (its units apply)' if name=='itsa_indicator_statistics.csv' else
            'valid observed coherence (dimensionless)')
        return statistics[column]+' of '+quantity+'. Plot-sample quantiles use the saved common-pixel sample where applicable.'
    match=re.fullmatch(r'(.+)_(median|valid_count|valid_fraction|mad_sigma)',column)
    if match:
        meanings={'median':'median valid coherence','valid_count':'valid member-pixel count',
                  'valid_fraction':'valid member-pixel fraction (0–1)','mad_sigma':'coherence MAD scaled by 1.4826'}
        return f'For acquisition {match[1]}: {meanings[match[2]]}.'
    match=re.fullmatch(r'class_(\d+)_(count|share_valid)',column)
    if match:
        return f'Class {match[1]} '+('pixel count.' if match[2]=='count' else 'share of valid classes 1–5 (0–1); NoData share undefined.')
    if column in ('mean_fraction','std_fraction','minimum_fraction','maximum_fraction'):
        return column.split('_')[0]+' of finite recovery-fraction values among loss-eligible classes 1–3 (dimensionless).'
    return 'Unmapped column: inspect the producing code before interpretation; no meaning is assumed.'


def file_description(name):
    lower=name.lower()
    if lower.startswith('distribution_'):
        return 'Histogram of '+lower.removeprefix('distribution_').removesuffix('.csv')+'; units follow that indicator.'
    return FILES.get(lower,'Unregistered filename: columns below come from the actual CSV header; verify its producer.')


def escape(value):
    return str(value).replace('|','\\|').replace('\n',' ').replace('\r',' ')


def generate(root, folders=None):
    root=Path(root).resolve()
    folders=folders or [root/'CSV files',root/'RII'/'CSV']
    found=set()
    for folder in folders:
        folder=Path(folder).resolve()
        if not folder.is_relative_to(root):
            raise ValueError('CSV search directories must be inside the Codes project.')
        if folder.is_dir():
            found.update(p.resolve() for p in folder.rglob('*.csv') if p.resolve().is_relative_to(root))
    lines=['# CSV files explained','',
           'Generated from filenames and actual CSV headers. CSV values and all ITSA/RRI calculations are unchanged.',
           'Blank/NaN entries indicate unavailable or suppressed results, not zero. Counts are not independent statistical replicates.',
           'ITSA classes: 0 NoData; 1 recovered; 2 partial; 3 not recovered/degraded; 4 no detectable loss; 5 indeterminate.',
           'RRI interpretation: 0 NoData; 1 further degradation; 2 approximately zero recovery; 3 partial; 4 approximately baseline; 5 above baseline.',
           'RRI status: 0 outside/structurally unavailable; 1 valid; 2 no detectable disturbance; 3 insufficient observations; 4 invalid denominator; 5 invalid model; 6 unresolved evidence.',
           'RRI uses fitted intervention counterfactual and final observation. It differs from the existing fixed-PRE-median recovery fraction.',
           'Interpretation tolerances and processing settings belong to the saved run/configuration; labels do not establish physical reconstruction.','',
           f'CSV files found: {len(found)}. Rerun this generator after new outputs are created.','']
    for path in sorted(found):
        lines.extend(['## '+escape(path.relative_to(root).as_posix()),'',file_description(path.name),''])
        try:
            with path.open(encoding='utf-8-sig',newline='') as stream:
                header=next(csv.reader(stream),[])
            if not header:
                lines.extend(['Empty file: no header is available; column definitions cannot be verified.',''])
                continue
            lines.extend(['| Column | Meaning |','|---|---|'])
            lines.extend(f'| {escape(col) or "(unnamed index)"} | {escape(column_description(col,path.name))} |' for col in header)
            lines.append('')
        except (OSError,UnicodeError,csv.Error) as exc:
            lines.extend(['Could not read header: '+escape(exc),''])
    if not found:
        lines.extend(['## Expected output filenames','',
                      'No CSV outputs were found in the selected folders. The following is a filename guide, not evidence that these files were generated.',
                      'Actual per-file column tables will be added when outputs exist and this command is rerun.','',
                      '| Expected filename | Purpose |','|---|---|'])
        lines.extend(f'| {name} | {escape(description)} |' for name,description in sorted(FILES.items()))
        lines.extend(['','## Column reference for future outputs','',
                      'These definitions come from the current exporters. They are not a claim that every file contains every listed column.',
                      'After CSV generation, the per-file tables above will list only each file’s actual columns.','',
                      '| Column | Meaning |','|---|---|'])
        lines.extend(f'| {escape(name)} | {escape(description)} |' for name,description in sorted(COLUMNS.items()))
        lines.extend(f'| {name} | Arithmetic/distribution statistic of RRI, coherence, evidence score, or the indicator named by metric; units follow that quantity. |'
                     for name in ('mean','median','q25','q75','minimum','maximum','std_population'))
    destination=root/'CSV_Explained.md'
    destination.write_text('\n'.join(lines)+'\n',encoding='utf-8')
    return destination


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--csv-dir',action='append',type=Path,help='Optional project-local CSV directory; repeat for multiple folders')
    args=parser.parse_args()
    print(generate(Path(__file__).resolve().parent.parent,args.csv_dir))

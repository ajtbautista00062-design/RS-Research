"""Bounded-memory ITSA research path; original observed workflow remains selectable."""
from contextlib import ExitStack
from datetime import datetime, timezone
import hashlib
import json
import logging
from pathlib import Path
import numpy as np
import pandas as pd
import rasterio
from rasterio.windows import Window, transform as window_transform
from .temporal import validate, local_path, day_offsets
from .spatial import inventory, aligned, read_stack, read_zones, mask_for, windows
from .metrics import calculate_fraction
from .recovery import assess, Assessment, CLASSES, QUALITY, REASONS
from .recovery_io import create_outputs, write, class_rows, METRICS
from .spatial_objects import (cell_sizes, pre_change_only, segment, region_assessments,
                              regularize, paint)
from .itsa_diagnostics import Diagnostics, write_records
from .sampling import Sample
from .itsa_zonal import ZonalSummary
from .output_layout import create_layout
from .analysis_boundaries import load_barangays, owners, BuiltupLayer, combine_strata
from .spatial_enhancement import evidence_score


def validate_settings(c):
    # Public settings are authoritative; older function interfaces retain aliases.
    c.MIN_REGION_AREA_M2 = c.MIN_MAPPING_AREA_M2
    c.OBJECT_MIN_VALID_FRACTION = c.MIN_PER_DATE_VALID_COVERAGE
    if not c.PRESERVE_ORIGINAL_CLASSIFICATION:
        raise ValueError("Original classification must be retained for comparison.")
    if c.ENABLE_BARANGAY_OUTPUTS and c.ALL_TOUCHED:
        raise ValueError("Barangay exports require ALL_TOUCHED=False for consistent pixel-center ownership.")
    if not .8 <= c.EVIDENCE_PROTECTION_SCORE <= 1:
        raise ValueError("Evidence protection must be between 0.80 and 1.00.")
    if not 0 < c.NEIGHBOR_MIN_AGREEMENT <= 1:
        raise ValueError("Neighbor agreement must be in (0,1].")
    if (not np.isfinite(c.NEIGHBOR_RADIUS_M) or not np.isfinite(c.NEIGHBOR_SIGNATURE_TOLERANCE)
            or c.NEIGHBOR_RADIUS_M <= 0 or c.NEIGHBOR_SIGNATURE_TOLERANCE <= 0):
        raise ValueError("Neighbor distance and signature tolerance must be positive.")
    if c.CLASSIFICATION_MODE not in ("pixel_itsa", "object_itsa", "spatially_regularized_itsa"):
        raise ValueError("Unsupported CLASSIFICATION_MODE")
    if not c.USE_FITTED_ITSA_FOR_RECOVERY or c.RECOVERY_REFERENCE != "pre_median":
        raise ValueError("ITSA mode requires fitted ITSA and fixed PRE median recovery target.")
    positive = ("EVENT_LOSS_MINIMUM", "ITSA_NOISE_FLOOR", "MIN_RECOVERY_IMPROVEMENT",
                "TARGET_REGION_AREA_M2", "MIN_REGION_AREA_M2", "MAX_REGION_AREA_M2",
                "MAX_DESIGN_CONDITION", "MAX_ITSA_RMSE", "EVENT_SUPPORT_MINIMUM",
                "OBSERVATION_DISAGREEMENT_TOLERANCE", "REGULARIZATION_FEATURE_TOLERANCE")
    for key in positive:
        value = getattr(c, key)
        if value is None or not np.isfinite(value) or value <= 0:
            raise ValueError(f"{key} must be a finite positive number.")
    if not 3 <= c.ITSA_MIN_PRE <= 8 or not 3 <= c.ITSA_MIN_POST <= 8:
        raise ValueError("ITSA requires 3..8 actual PRE/POST observations per fit.")
    if not 2 <= c.LATE_POST_CONFIRMATIONS <= c.ITSA_MIN_POST:
        raise ValueError("Require two or more consecutive late POST confirmations.")
    if not 0 < c.OBJECT_MIN_VALID_FRACTION <= 1 or c.OBJECT_MIN_VALID_PIXELS < 1:
        raise ValueError("Invalid regional valid-pixel support settings.")
    if not c.MIN_REGION_AREA_M2 <= c.TARGET_REGION_AREA_M2 <= c.MAX_REGION_AREA_M2:
        raise ValueError("Require min <= target <= maximum region area.")
    if len(c.PRE_FEATURE_BIN_WIDTHS) != 3 or any(v <= 0 for v in c.PRE_FEATURE_BIN_WIDTHS):
        raise ValueError("Three positive PRE feature bin widths are required.")
    if c.EFFECT_DECISION_MARGIN < 0 or c.PRE_RESIDUAL_NOISE_MULTIPLIER < 0:
        raise ValueError("Effect margin/noise multiplier cannot be negative.")
    if c.REGULARIZATION_MIN_NEIGHBORS < 1 or c.MIN_LATE_CONFIRM_SPAN_DAYS < 0:
        raise ValueError("Invalid spatial neighbor count or confirmation span.")
    if not c.ITSA_NOISE_FLOOR <= 1 or not 0 <= c.ITSA_BASELINE_MINIMUM <= 1:
        raise ValueError("Invalid coherence-unit noise floor or baseline limit.")
    if c.SOURCE_VALIDITY_CONFIRMED and not c.SOURCE_VALIDITY_NOTES.strip():
        raise ValueError("SOURCE_VALIDITY_CONFIRMED requires SOURCE_VALIDITY_NOTES.")
    if not c.EXPORT_RAW_CLASSIFICATION:
        raise ValueError("Raw pixel classification is required for transparent comparison.")
    if c.CLASSIFICATION_MODE == "spatially_regularized_itsa" and not c.ENABLE_SPATIAL_REGULARIZATION:
        raise ValueError("Regularized mode requires ENABLE_SPATIAL_REGULARIZATION.")
    if c.ZERO_IS_INVALID and not (c.SOURCE_VALIDITY_CONFIRMED and c.SOURCE_VALIDITY_NOTES.strip()):
        raise ValueError("Excluding exact zeros requires verified source semantics and notes.")


def paint_result(ids, records, stack, mask, c):
    available = np.isfinite(stack).any(0) & mask
    labels = paint(ids, records, "classes", "uint8", 0)
    labels[(ids == 0) & available] = 5
    reason = paint(ids, records, "reasons", "uint8", 1)
    reason[(ids == 0) & available] = 2
    reason[~mask] = 11
    quality = paint(ids, records, "quality", "uint16", 0)
    quality[~mask] = 65535
    gates = {}
    if records:
        for key in next(iter(records.values()))[0].gates:
            gates[key] = paint(ids, records, "gate:"+key, "uint8", 0).astype(bool)
    return Assessment(labels, quality, reason, {key: paint(ids, records, key) for key in METRICS}, None, gates)


def object_window(stack, raw, win, ref, mask, barriers, obs, c, pixel_area, writers, diagnostic):
    cell, domain = cell_sizes(pixel_area, c)
    ids = np.zeros(mask.shape, dtype="uint32")
    merged = ids.copy()
    records, merged_records = {}, {}
    pre_indices = np.flatnonzero([o.phase == "pre" for o in obs])
    pre = stack[pre_indices]
    change = pre_change_only(pre, day_offsets(obs, c.EVENT_ORIGIN)[pre_indices])
    for r in range(0, stack.shape[1], domain):
        for s in range(0, stack.shape[2], domain):
            sl = np.s_[r:min(r+domain, stack.shape[1]), s:min(s+domain, stack.shape[2])]
            values = stack[:, sl[0], sl[1]]
            local = segment(pre[:, sl[0], sl[1]], None, change[sl], mask[sl], barriers[sl], cell,
                            int(win.row_off)+r, int(win.col_off)+s, ref.width, c)
            rec = region_assessments(values, local, obs, c, pixel_area)
            ids[sl] = local
            records.update(rec)
            write_records(writers["raw"], local, rec, pixel_area, obs, "object")
            reg, rrec = (regularize(values, local, rec, barriers[sl], obs, c, pixel_area)
                         if c.ENABLE_SPATIAL_REGULARIZATION else (local, rec))
            merged[sl] = reg
            merged_records.update(rrec)
            write_records(writers["regularized"], reg, rrec, pixel_area, obs, "regularized")
    obj = paint_result(ids, records, stack, mask, c)
    reg = paint_result(merged, merged_records, stack, mask, c)
    if c.ENABLE_THRESHOLD_SENSITIVITY_REPORT and records:
        series = np.column_stack([v[1] for v in records.values()])
        keys, counts = np.unique(ids[ids > 0], return_counts=True)
        count_map = dict(zip(keys, counts))
        areas = np.array([count_map[key] for key in records])
        codes = np.array([int(v[0].classes[0]) for v in records.values()])
        diagnostic.sweep(series, obs, c, areas, codes)
    return obj, reg, ids, merged, records


def run(c):
    validate_settings(c)
    obs = validate(c)
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    layout = create_layout(c, f"itsa_{obs[0].start:%Y%m%d}_{obs[-1].end:%Y%m%d}_{stamp}")
    out, figures, rasters = layout.tables, layout.graphs, layout.rasters
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s",
                        handlers=[logging.FileHandler(out / "processing.log", encoding="utf-8"), logging.StreamHandler()], force=True)
    logging.warning("Operational ITSA: thresholds provisional; source semantics confirmed=%s; no spatial filling.", c.SOURCE_VALIDITY_CONFIRMED)
    diagnostic, sample = Diagnostics(), Sample(c.SAMPLE_PIXELS, c.RANDOM_SEED)
    study_count = 0
    observation_counts = np.zeros(len(obs), dtype='int64')
    object_counts = {"object": np.zeros(6, dtype="int64"), "regularized": np.zeros(6, dtype="int64")}
    with ExitStack() as context:
        srcs = [context.enter_context(rasterio.open(o.path)) for o in obs]
        ref = srcs[0]
        pd.DataFrame(inventory(srcs, obs, c)).to_csv(out / "raster_inventory.csv", index=False)
        unit_factor = ref.crs.linear_units_factor[1]
        pixel_area = abs(ref.transform.a * ref.transform.e) * unit_factor**2
        if pixel_area > c.MAX_REGION_AREA_M2:
            raise ValueError("MAX_REGION_AREA_M2 is smaller than one reference pixel.")
        if ref.width * ref.height >= np.iinfo("uint32").max:
            raise ValueError("Reference grid exceeds available region-ID precision.")
        zones = read_zones(c.Study_Area, ref.crs, c.PROJECT_ROOT) if c.Study_Area else []
        geoms = [g for _, g in zones]
        barangays, boundary_report = load_barangays(c.BARANGAY_BOUNDARY, ref.crs, c)
        (out / 'barangay_boundary_audit.json').write_text(json.dumps(boundary_report, indent=2), encoding='utf-8')
        builtup_source = BuiltupLayer(context, ref, c)
        zonal = {name: ZonalSummary(name) for name, _ in zones}
        barrier_src = None
        if c.ANCILLARY_REGION_BARRIERS:
            barrier_path = local_path(c.ANCILLARY_REGION_BARRIERS, c.PROJECT_ROOT)
            barrier_src = context.enter_context(rasterio.open(barrier_path))
            if (barrier_src.crs != ref.crs or barrier_src.transform != ref.transform or barrier_src.shape != ref.shape):
                raise ValueError("Ancillary barriers require the exact reference grid; no automatic label interpolation.")
            if not np.issubdtype(np.dtype(barrier_src.dtypes[0]), np.integer):
                raise ValueError("Ancillary barriers must be positive integer categorical IDs; zero is excluded.")
        vrts = [context.enter_context(aligned(src, ref, c)) for src in srcs]
        outputs = create_outputs(context, rasters, ref, c, obs)
        writers = {name: [None, context.enter_context((out / f"{name}_region_statistics.csv").open("w", newline="", encoding="utf-8"))]
                   for name in ("raw", "regularized")}
        _, domain = cell_sizes(pixel_area, c)
        # I/O chunks are multiples of physical domains: no split objects at chunk seams.
        chunk = domain * max(1, c.WINDOW_SIZE // domain)
        for win in windows(ref.width, ref.height, chunk):
            stack = read_stack(vrts, srcs, win, c)
            shape = stack.shape[1:]
            mask = mask_for(geoms, shape, window_transform(win, ref.transform), c.ALL_TOUCHED)
            study_count += int(mask.sum())
            barriers = (barrier_src.read(1, window=win, masked=True).filled(0)
                        if barrier_src else np.ones(shape, dtype="uint8"))
            analysis_mask = mask & (barriers > 0)
            landcover = builtup_source.read(win, window_transform(win, ref.transform))
            strata = combine_strata(barriers, landcover, c.ENABLE_BUILTUP_PRIORITY)
            strata[~analysis_mask] = 0
            if strata.max() > np.iinfo('uint32').max:
                raise ValueError('Combined ancillary IDs exceed uint32; recode the external labels.')
            ownership = owners(barangays, shape, window_transform(win, ref.transform))
            ownership[~mask] = 0
            stack[:, ~analysis_mask] = np.nan
            observation_counts += np.isfinite(stack).sum(axis=(1,2))
            raw = assess(stack, obs, c)
            raw.reasons[~analysis_mask] = 11
            raw.quality[~analysis_mask] = 65535
            legacy = calculate_fraction(stack, obs, c).classes
            legacy = np.where(legacy == -9999, 0, legacy).astype("uint8")
            obj, reg, ids, merged, records = object_window(stack, raw, win, ref, analysis_mask,
                                                         strata, obs, c, pixel_area, writers, diagnostic)
            primary = {"pixel_itsa": raw, "object_itsa": obj, "spatially_regularized_itsa": reg}[c.CLASSIFICATION_MODE]
            products = dict(pixel=raw.classes, legacy=legacy, object=obj.classes, regularized=reg.classes)
            if any(not np.isin(values, range(6)).all() for values in products.values()):
                raise RuntimeError("Undefined recovery class before export.")
            for name, values in (("ITSA_Recovery_Classes", primary.classes),
                                 ("Raw_Pixel_ITSA_Classes", raw.classes), ("Legacy_Observed_Classes", legacy),
                                 ("Object_ITSA_Classes", obj.classes), ("Spatially_Regularized_ITSA_Classes", reg.classes),
                                 ("Spatial_Object_ID", ids), ("Regularized_Object_ID", merged),
                                 ("ITSA_Quality_Flags", primary.quality), ("ITSA_Exclusion_Reasons", primary.reasons),
                                 ("Pixel_Object_Classification_Change", np.where(mask, raw.classes != obj.classes, 255)),
                                 ("Regularization_Classification_Change", np.where(mask, obj.classes != reg.classes, 255))):
                if name in outputs:
                    write(outputs[name], values, win)
            write(outputs["ITSA_Metrics"], np.stack([primary.metrics[key] for key in METRICS]), win)
            write(outputs['Original_ITSA_Recovery_Classification'], primary.classes, win)
            write(outputs['Spatial_Recovery_Regions'], merged if c.CLASSIFICATION_MODE == 'spatially_regularized_itsa' else ids, win)
            write(outputs['Original_Evidence_Confidence'], np.where(primary.classes > 0, evidence_score(primary, c), np.nan), win)
            write(outputs['Raw_ITSA_Metrics'], np.stack([raw.metrics[k] for k in METRICS]), win)
            write(outputs['Raw_ITSA_Reasons'], raw.reasons, win)
            write(outputs['Raw_ITSA_Quality_Flags'], raw.quality, win)
            write(outputs['Spatial_Strata'], strata, win)
            write(outputs['Builtup_Status'], np.where(mask, landcover, -9999), win)
            write(outputs['Barangay_Owner'], ownership, win)
            write(outputs['Aligned_Coherence_Observations'], stack, win)
            write(outputs['Valid_Data_Coverage'], np.where(analysis_mask[None], np.isfinite(stack), 255), win)
            diagnostic.add(products, primary, mask)
            for name, geom in zones:
                zm = mask_for([geom], shape, window_transform(win, ref.transform), c.ALL_TOUCHED)
                zonal[name].add(primary.classes, primary.metrics["recovery_fraction"], zm & mask)
            for mode, path_ids, path_classes in (("object", ids, obj.classes), ("regularized", merged, reg.classes)):
                _, indices = np.unique(path_ids[path_ids > 0], return_index=True)
                object_counts[mode] += np.bincount(path_classes[path_ids > 0][indices], minlength=6)
            complete = np.isfinite(stack).all(0)
            if complete.any():
                rr, cc = np.where(complete)
                sample.add(stack[:, complete], rr+int(win.row_off), cc+int(win.col_off))
            logging.info("Completed input window row=%d col=%d; objects=%d", win.row_off, win.col_off, len(records))
        paths = {key: Path(ds.name) for key, ds in outputs.items()}
        geometry = (ref.crs, ref.transform, ref.shape)
    from .enhancement_reports import enhance_outputs
    enhanced_summary = enhance_outputs(paths, layout, c)
    from .barangay_outputs import export_barangays
    export_barangays(paths, barangays, layout, obs, c)
    # Closed-file header checks only; class/count QA was accumulated during writing.
    for key, path in paths.items():
        with rasterio.open(path) as check:
            if (check.crs, check.transform, check.shape) != geometry:
                raise RuntimeError(f"Output grid QA failed: {key}")
    rows = [row for name, counts in diagnostic.counts.items() for row in class_rows(name, counts, study_count, pixel_area)]
    if any(int(counts.sum()) != study_count for counts in diagnostic.counts.values()):
        raise RuntimeError("Class coverage totals do not equal the study-mask denominator.")
    pd.DataFrame(rows).to_csv(out / "class_area_statistics.csv", index=False)
    pd.DataFrame([dict(observation=o.name, phase=o.phase, master=o.master, slave=o.slave,
                      reference=o.reference, valid_pixels=int(n), study_pixels=study_count,
                      valid_fraction=n/study_count if study_count else np.nan)
                  for o,n in zip(obs,observation_counts)]).to_csv(out / 'Valid_Data_Coverage.csv', index=False)
    if zonal:
        pd.DataFrame([z.row() for z in zonal.values()]).to_csv(out / "zonal_recovery_statistics.csv", index=False)
    pd.DataFrame([dict(mode=mode, code=code, objects=int(n)) for mode, counts in object_counts.items()
                  for code, n in enumerate(counts)]).to_csv(out / "object_class_counts.csv", index=False)
    diagnostic.export(out, pixel_area)
    from .itsa_figures import render_itsa
    render_itsa(sample, obs, figures, c, out)
    from .enhancement_reports import render_comparison
    render_comparison(enhanced_summary, figures, c)
    source_hashes = {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in Path(__file__).parent.glob("*.py")}
    manifest = dict(settings={k: v for k, v in vars(c).items() if k.isupper() or k in ("Study_Area", "Output_Folder")},
                    dates=[dict(name=o.name, phase=o.phase, path=str(o.path), master=str(o.master), slave=str(o.slave),
                                reference=str(o.reference)) for o in obs], classes=CLASSES, quality=QUALITY,
                    reasons=REASONS, code_sha256=source_hashes, study_pixels=study_count, pixel_area_m2=pixel_area,
                    regularization="same-class pooled object refit; separately exported connected-neighbor refinement requiring own pixel ITSA support",
                    segmentation="PRE-only four-connected bins in fixed physical domains; no inferred river/coast barriers",
                    input_count=len(obs), scientific_status="exploratory until source semantics, thresholds and external evidence validated")
    manifest.update(output_roots=dict(rasters=str(rasters), graphs=str(figures), csv=str(out)),
                    barangay_boundary_audit=boundary_report,
                    enhancement_threshold_adjustment=False,
                    builtup_status='external verified mask' if c.BUILTUP_MASK else 'unknown; not inferred from coherence',
                    protection_score='normalized heuristic; not a calibrated probability')
    (out / "run_manifest.json").write_text(json.dumps(manifest, indent=2, default=str), encoding="utf-8")
    logging.info("Finished: %s", out)
    return 0

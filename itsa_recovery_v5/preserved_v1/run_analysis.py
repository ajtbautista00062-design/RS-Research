"""Discover actual coherence pairs, audit coverage, export descriptive measurements."""
from contextlib import ExitStack
from datetime import datetime, timezone
from pathlib import Path
import json
import logging
import importlib.util
import math
import time
import shutil
import numpy as np
import pandas as pd
import rasterio
from rasterio.enums import Resampling

from sarflow.temporal import discover, local_path
from sarflow.spatial import (inventory, verify_source_values, aligned, read_stack,
                             read_zones, mask_for, windows, export_profile)
from sarflow.recovery import baseline_recovery, settings as recovery_settings
from sarflow.recovery_io import open_outputs, verify as verify_recovery, finalize as finalize_recovery
from sarflow.sampling import Sample
from sarflow.figures import observed_series
from sarflow.itsa import design, band_names, raster_bands, MODEL, LIMITATION
from sarflow.zonal import ZoneStats
from sarflow.correction import corrected_window


def load_config():
    path = Path(__file__).resolve().with_name("config.py")
    spec = importlib.util.spec_from_file_location("sar_config", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def write_csv(path, rows):
    pd.DataFrame(rows).to_csv(path, index=False, encoding="utf-8-sig")


def observation_rows(obs):
    return [dict(name=o.name, phase=o.phase, path=str(o.path),
                 first_filename_date=o.master.isoformat(), second_filename_date=o.slave.isoformat(),
                 start=o.start.isoformat(), end=o.end.isoformat(), reference=o.reference.isoformat(),
                 midpoint=(o.start+(o.end-o.start)/2).isoformat(),
                 included_primary_itsa=o.phase != "event", pair_duration_days=(o.end-o.start).total_seconds()/86400) for o in obs]


def verify_export(path, ref, obs, c):
    names = band_names(obs)
    counts = np.zeros(len(names), dtype="int64")
    flags_path = path.with_name(path.name.replace("Infanta_Coherence_ITSA_",
                                                 "Infanta_Coherence_Correction_Flags_"))
    with rasterio.open(path) as src, rasterio.open(flags_path) as flags:
        if flags.crs != ref.crs or flags.transform != ref.transform or flags.shape != ref.shape or flags.count != len(obs):
            raise RuntimeError("Correction flags do not match the analysis grid.")
        if (src.crs != ref.crs or src.transform != ref.transform or src.shape != ref.shape
                or src.count != len(names) or src.nodata != c.NODATA
                or any(dtype != "float32" for dtype in src.dtypes)
                or src.descriptions != names):
            raise RuntimeError("GeoTIFF metadata verification failed.")
        for window in windows(src.width, src.height, c.WINDOW_SIZE):
            data = src.read(window=window)
            valid = data != c.NODATA
            if not np.isfinite(data).all():
                raise RuntimeError("Unencoded nonfinite output values.")
            n = len(obs)
            provenance = flags.read(window=window)
            if not np.isin(provenance,[0,1,2]).all() or not np.array_equal(provenance > 0,valid[:n]):
                raise RuntimeError("Correction provenance/support verification failed.")
            if np.any(valid[:n] & ((data[:n] < 0) | (data[:n] > 1))):
                raise RuntimeError("Observed coherence outside [0,1].")
            original = np.where(valid[:n], data[:n], np.nan)
            x = design(obs, c)
            expected = raster_bands(original, x)
            if not np.array_equal(valid[n:], np.isfinite(expected[n:])):
                raise RuntimeError("ITSA model/count/diagnostic masks differ.")
            if not np.allclose(data[n:][valid[n:]], expected[n:][valid[n:]], atol=1e-6):
                raise RuntimeError("ITSA arithmetic verification failed.")
            # Independent least-squares check of up to three real pixel series per block.
            complete = np.isfinite(original[x[:,0] == 1]).all(axis=0)
            rr, cc = np.where(complete)
            for r, col in zip(rr[:3], cc[:3]):
                use = x[:,0] == 1
                beta = np.linalg.lstsq(x[use], original[use,r,col].astype(float), rcond=None)[0]
                if not np.allclose(data[n:n+4,r,col], beta, atol=1e-6):
                    raise RuntimeError("Independent pixel OLS verification failed.")
            counts += valid.sum(axis=(1, 2))
    return {name: int(n) for name, n in zip(names, counts)}


def export_windows(srcs, ref, obs, c, out, dates, resources):
    """No baseline cutoff or all-date coverage filter is applied to measurements."""
    vrts = [resources.enter_context(aligned(src, ref, c)) for src in srcs]
    footprints = [resources.enter_context(aligned(src, ref, c, Resampling.nearest)) for src in srcs]
    profile = export_profile(ref, c)
    names = band_names(obs)
    x = design(obs, c)
    profile.update(count=len(names))
    path = out / f"Infanta_Coherence_ITSA_{dates}.partial.tif"
    dst = resources.enter_context(rasterio.open(path, "w", **profile))
    for j, name in enumerate(names, start=1):
        dst.set_band_description(j, name)
        dst.update_tags(j, units="coherence/day" if "per_day" in name else "count" if "count" in name else "dimensionless" if "squared" in name else "coherence",
                        validity="Observed bands: per-date support; model bands: complete PRE/POST only; count allows partial support")
    dst.update_tags(event_origin_utc=c.EVENT_ORIGIN, pair_reference=c.REFERENCE_DATE,
                    model=MODEL, limitation=LIMITATION,
                    residual_degrees_of_freedom=str(len(obs)-5),
                    zero_is_invalid=str(c.ZERO_IS_INVALID),
                    comparability_confirmed=str(c.ACQUISITION_COMPARABILITY_CONFIRMED))
    provenance_path = out/f"Infanta_Coherence_Correction_Flags_{dates}.partial.tif"
    provenance_dst = resources.enter_context(rasterio.open(provenance_path,"w",
        **(profile | dict(count=len(obs),dtype="uint8",nodata=0,predictor=2))))
    for j,o in enumerate(obs,1):
        provenance_dst.set_band_description(j,o.name)
    provenance_dst.update_tags(code_0="missing or outside boundary",code_1="valid aligned input",
        code_2="nearest-neighbor spatial estimate",max_hole_pixels=str(c.FILL_MAX_HOLE_PIXELS),
        max_distance_pixels=str(c.FILL_RADIUS_PIXELS))
    dst.update_tags(correction="small enclosed holes filled before ITSA; see correction flags",
                    boundary=c.Study_Area)
    class_dst = fraction_dst = None
    if c.ENABLE_BASELINE_CLASSIFICATION:
        class_dst, fraction_dst = open_outputs(resources, profile, out, dates, c)
    per_observation = [Sample(c.SAMPLE_PIXELS, c.RANDOM_SEED+i, 1) for i in range(len(obs))]
    common = Sample(c.SAMPLE_PIXELS, c.RANDOM_SEED, len(obs))
    zones = read_zones(c.Study_Area, ref.crs, c.PROJECT_ROOT) if c.Study_Area else []
    zonal_names = names + ("baseline_recovery_fraction",) if c.ENABLE_BASELINE_CLASSIFICATION else names
    zonal = {name: ZoneStats(name, zonal_names) for name, _ in zones}
    audit = dict(reference_grid_pixels=ref.width*ref.height, pixels_in_analysis_mask=0,
                 any_observation_valid=0, common_all_observations=0,
                 historical_040_cutoff_would_retain_on_common_latest_pre=0,
                 nearest_filled_pixels_by_observation=[0]*len(obs))
    latest_pre = [i for i, o in enumerate(obs) if o.phase == "pre"][-1]
    total_windows = math.ceil(ref.width/c.WINDOW_SIZE) * math.ceil(ref.height/c.WINDOW_SIZE)
    started = last_report = time.monotonic()
    logging.info("Raster processing: %d windows, %d observations, %d output bands. Writing %s",
                 total_windows, len(obs), len(names), path.name)
    for block, window in enumerate(windows(ref.width, ref.height, c.WINDOW_SIZE), 1):
        if block == 1:
            logging.info("Reading and aligning first window; inputs remain unchanged.")
        stack, correction_flags = corrected_window(vrts,srcs,ref,window,c,footprints,
                                                    [g for _,g in zones])
        provenance_dst.write(correction_flags,window=window)
        for j in range(len(obs)):
            audit["nearest_filled_pixels_by_observation"][j] += int((correction_flags[j] == 2).sum())
        shape = stack.shape[1:]
        transform = rasterio.windows.transform(window, ref.transform)
        mask = mask_for([g for _, g in zones], shape, transform, c.ALL_TOUCHED)
        stack[:, ~mask] = np.nan
        valid = np.isfinite(stack)
        complete = valid.all(axis=0)
        audit["pixels_in_analysis_mask"] += int(mask.sum())
        audit["any_observation_valid"] += int(valid.any(axis=0).sum())
        audit["common_all_observations"] += int(complete.sum())
        audit["historical_040_cutoff_would_retain_on_common_latest_pre"] += int((complete & (stack[latest_pre] >= .4)).sum())
        for j, sample in enumerate(per_observation):
            r, col = np.where(valid[j])
            if r.size:
                sample.add(stack[j:j+1, valid[j]], r+int(window.row_off), col+int(window.col_off))
        r, col = np.where(complete)
        if r.size:
            common.add(stack[:, complete], r+int(window.row_off), col+int(window.col_off))
        codes = fraction = None
        if c.ENABLE_BASELINE_CLASSIFICATION:
            codes, fraction = baseline_recovery(stack, obs, c)
        bands = raster_bands(stack, x)
        output_block = np.ascontiguousarray(np.where(
            np.isfinite(bands), bands, c.NODATA), dtype="float32")
        try:
            dst.write(output_block, window=window)
        except (Exception, KeyboardInterrupt):
            logging.exception("Raster write stopped at window %d/%d, row=%d col=%d; array shape=%s",
                              block, total_windows, window.row_off, window.col_off, output_block.shape)
            raise
        if class_dst is not None:
            class_dst.write(codes, 1, window=window)
            fraction_dst.write(np.where(np.isfinite(fraction), fraction, c.NODATA).astype("float32"),
                               1, window=window)
        zonal_bands = np.concatenate((bands, fraction[None]), axis=0) if fraction is not None else bands
        for name, geom in zones:
            zonal[name].add(zonal_bands, mask_for([geom], shape, transform, c.ALL_TOUCHED), codes)
        now = time.monotonic()
        if block == 1 or block == total_windows or now-last_report >= 10:
            logging.info("Raster progress: %d/%d windows (%.1f%%), elapsed %.1f min",
                         block, total_windows, 100*block/total_windows, (now-started)/60)
            last_report = now
    if audit["pixels_in_analysis_mask"] == 0:
        raise ValueError("Boundary does not overlap the reference grid; check CRS/extent.")
    logging.info("Flushing and closing compressed GeoTIFF; this may take additional time.")
    dst.close()
    provenance_dst.close()
    if class_dst is not None:
        class_dst.close()
        fraction_dst.close()
    return path, per_observation, common, zonal, audit


def summaries(per_observation, common, obs, out, c, dates):
    rows = []
    for sample, o in zip(per_observation, obs):
        row = sample.table([o.name]).to_dict("records")[0] if sample.count else dict(
            observation=o.name, complete_case_pixels=0, sample_pixels=0,
            mean_exact=np.nan, q25_sample=np.nan, median_sample=np.nan, q75_sample=np.nan)
        row["valid_pixels"] = row.pop("complete_case_pixels")
        rows.append(row)
    table = pd.DataFrame(rows)
    table.to_csv(out / "timeseries_available_support.csv", index=False)
    observed_series(table, obs, out, c, f"Coherence_AvailableSupport_{dates}",
                    "Per-observation available support; changing coverage can change spatial summaries")
    if common.count:
        common_table = common.table([o.name for o in obs]).rename(columns={"complete_case_pixels": "valid_pixels"})
        common_table.to_csv(out / "timeseries_common_support.csv", index=False)
        observed_series(common_table, obs, out, c, f"Coherence_CommonSupport_{dates}",
                        "Fixed common spatial support across all unique observations; no recovery filtering")
    else:
        logging.warning("No common pixels: common-support figure omitted; available-support data retained.")


def process_run(c, obs, duplicates, out, dates):
    graphs = out.parent / "Time Series Graphs" / out.name
    rasters = out.parent / "Recovery Map Raster File" / out.name
    graphs.mkdir(parents=True)
    rasters.mkdir(parents=True)
    logging.warning(LIMITATION)
    if c.ENABLE_BASELINE_CLASSIFICATION:
        logging.warning("Baseline recovery uses operational thresholds: %s", recovery_settings(c))
        if not c.ACQUISITION_COMPARABILITY_CONFIRMED:
            logging.warning("Acquisition comparability remains unverified: recovery classes are exploratory.")
    logging.info("Unique observations: %s", {p: sum(o.phase == p for o in obs) for p in ("pre", "event", "post")})
    write_csv(out/"observation_dates.csv", observation_rows(obs))
    write_csv(graphs/"pilot_metadata.csv", observation_rows(obs))
    pd.DataFrame(duplicates, columns=["file", "kept", "reason"]).to_csv(out/"duplicate_files.csv", index=False)
    logging.warning("Zero policy retained: %s; nodata=None does not establish that zero is invalid.", c.ZERO_IS_INVALID)
    with ExitStack() as resources:
        srcs = [resources.enter_context(rasterio.open(o.path)) for o in obs]
        reference_path = local_path(c.REFERENCE_RASTER, c.PROJECT_ROOT) if c.REFERENCE_RASTER else next(o.path for o in obs if o.phase == "event")
        ref = resources.enter_context(rasterio.open(reference_path))
        metadata = inventory(srcs, obs, c, ref)
        for index, (row, src) in enumerate(zip(metadata, srcs), 1):
            logging.info("Source validation %d/%d: %s", index, len(srcs), obs[index-1].name)
            row.update(verify_source_values(src, c))
            logging.info("Validated source: %s policy-valid pixels", format(row["valid_after_zero_policy"], ","))
        write_csv(out/"source_inventory_and_validity.csv", metadata)
        path, samples, common, zones, audit = export_windows(srcs, ref, obs, c, rasters, dates, resources)
        audit["exported_valid_pixels_by_band"] = verify_export(path, ref, obs, c)
        if c.ENABLE_BASELINE_CLASSIFICATION:
            audit["recovery_class_counts"] = verify_recovery(rasters, dates, path, ref, obs, c)
        grid = dict(crs=str(ref.crs), transform=tuple(ref.transform), shape=ref.shape, reference_file=str(reference_path))
    finished_path = path.with_name(path.name.replace(".partial.tif", ".tif"))
    path.rename(finished_path)
    flags = rasters/f"Infanta_Coherence_Correction_Flags_{dates}.partial.tif"
    flags.rename(flags.with_name(flags.name.replace(".partial.tif",".tif")))
    if c.ENABLE_BASELINE_CLASSIFICATION:
        finalize_recovery(rasters, dates)
    logging.info("Raster integrity checks passed: %s", finished_path.name)
    logging.info("Exporting time-series CSVs and PNG/PDF/SVG figures.")
    summaries(samples, common, obs, graphs, c, dates)
    from sarflow.pilot_exports import export_maps
    export_maps(finished_path, rasters, graphs, obs, c, dates)
    if zones:
        write_csv(out/"zonal_continuous_statistics.csv", [row for z in zones.values() for row in z.rows()])
        if c.ENABLE_BASELINE_CLASSIFICATION:
            write_csv(out/"zonal_class_counts.csv", [z.classes_row() for z in zones.values()])
    manifest = dict(observations=observation_rows(obs), duplicate_files=duplicates, grid=grid,
                    boundary_file=c.Study_Area,
                    correction=dict(enabled=c.FILL_SMALL_HOLES,max_area_pixels=c.FILL_MAX_HOLE_PIXELS,
                                    max_distance_pixels=c.FILL_RADIUS_PIXELS,
                                    stage="after grid alignment, before all temporal analysis"),
                    coverage=audit, band_definitions=dict(enumerate(band_names(obs), 1)),
                    itsa_model=MODEL, itsa_limitation=LIMITATION,
                    classification_exported=bool(c.ENABLE_BASELINE_CLASSIFICATION),
                    baseline_recovery_settings=recovery_settings(c),
                    zero_is_invalid=c.ZERO_IS_INVALID, resampling=c.RESAMPLING,
                    source_footprint="nearest-neighbor validity constrains interpolated data",
                    comparability_confirmed=c.ACQUISITION_COMPARABILITY_CONFIRMED,
                    comparability_notes=c.COMPARABILITY_NOTES,
                    unsupported=["validated statistical significance", "physical reconstruction inference"])
    (out/"run_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    (out/"pilot_results_summary.md").write_text(
        "# Real-data coherence pilot" + chr(10) + LIMITATION + chr(10) +
        "Primary regression excludes the event pair; plots retain it. Regional medians/IQR "
        "are estimated from the seeded spatial sample. Pixel fits use actual PRE/POST values. "
        "Recovery classes use the PRE median. Small-hole nearest-neighbor estimates enter all analyses; correction flags identify them. "
        "Thresholds and acquisition comparability require validation. No structural recovery "
        "is established. Source audits and class counts are in run_manifest.json." + chr(10) +
        f"Common spatial support: {common.count} pixels. Graphs: {graphs}. Rasters: {rasters}.",
        encoding="utf-8")
    shutil.copy2(out/"run_manifest.json", rasters/"run_manifest.json")
    logging.info("Finished. Graphs: %s; Rasters: %s", graphs, rasters)


def main():
    c = load_config()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s",
                        handlers=[logging.StreamHandler()])
    logging.info("Discovering input dates and checking duplicate files; hashing reads each source.")
    obs, duplicates = discover(c)
    design(obs, c)  # Validate temporal design before creating outputs.
    dates = f"{min(o.start for o in obs):%Y%m%d}_{max(o.end for o in obs):%Y%m%d}"
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    out = local_path(c.Output_Folder, c.PROJECT_ROOT) / f"series_{dates}_{stamp}"
    out.mkdir(parents=True, exist_ok=False)
    handler = logging.FileHandler(out/"processing.log", encoding="utf-8")
    handler.setFormatter(logging.Formatter("%(asctime)s %(levelname)s %(message)s"))
    logging.getLogger().addHandler(handler)
    state_path = out/"run_state.json"
    state_path.write_text(json.dumps(dict(status="running", output_directory=str(out)), indent=2), encoding="utf-8")
    try:
        process_run(c, obs, duplicates, out, dates)
    except (Exception, KeyboardInterrupt) as exc:
        status = "interrupted" if isinstance(exc, KeyboardInterrupt) else "failed"
        logging.exception("Run %s: %s: %s. Incomplete files are not validated outputs.",
                          status, type(exc).__name__, exc)
        state_path.write_text(json.dumps(dict(status=status, error_type=type(exc).__name__,
                                              error=str(exc)), indent=2), encoding="utf-8")
        raise
    else:
        state_path.write_text(json.dumps(dict(status="complete", output_directory=str(out)), indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()


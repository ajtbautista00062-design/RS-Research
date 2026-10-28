"""Run the configurable, windowed Infanta SAR coherence analysis.

Edit config.py first. This script analyzes existing coherence magnitude rasters;
it does not create interferograms or infer missing acquisition dates.
"""
from __future__ import annotations

import importlib.util
import json
import logging
from pathlib import Path
import sys

import numpy as np
import pandas as pd
import rasterio

from sarflow.temporal import validate
from sarflow.spatial import (aligned, export_profile, inventory, mask_for,
                             read_stack, read_zones, verify_source_values, windows)
from sarflow.metrics import calculate, calculate_fraction
from sarflow.sampling import Sample
from sarflow.zonal import ZoneStats
from sarflow.figures import portfolio

def load_config(path: Path):
    spec = importlib.util.spec_from_file_location("sar_config", path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot load configuration: {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def write_csv(path: Path, rows: list[dict]) -> None:
    if rows:
        pd.DataFrame(rows).to_csv(path, index=False, encoding="utf-8-sig")


def legacy_main() -> int:
    here = Path(__file__).resolve().parent
    c = load_config(here / "config.py")
    obs = validate(c)
    out = Path(c.Output_Folder).expanduser().resolve()
    out.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s %(levelname)s %(message)s",
                        handlers=[logging.FileHandler(out / "processing.log", encoding="utf-8"),
                                  logging.StreamHandler(sys.stdout)])
    logging.info("Validated %d rasters and pair dates.", len(obs))
    srcs = [rasterio.open(o.path) for o in obs]
    try:
        metadata = inventory(srcs, obs, c)
        for src in srcs:
            verify_source_values(src, c)
        write_csv(out / "raster_inventory.csv", metadata)
        pd.DataFrame([{"name": o.name, "phase": o.phase, "master_utc": o.master.isoformat(),
                       "slave_utc": o.slave.isoformat(), "pair_start_utc": o.start.isoformat(),
                       "pair_end_utc": o.end.isoformat(), "reference_utc": o.reference.isoformat(),
                       "temporal_baseline_days": (o.end-o.start).total_seconds()/86400}
                      for o in obs]).to_csv(out / "observation_dates.csv", index=False)
        ref = srcs[0]
        zones = read_zones(c.Study_Area, ref.crs, c.PROJECT_ROOT) if c.Study_Area else []
        zone_acc = {name: ZoneStats(name) for name, _ in zones}
        sample = Sample(c.SAMPLE_PIXELS, c.RANDOM_SEED)
        zone_samples = {name: Sample(c.SAMPLE_PIXELS, c.RANDOM_SEED + i + 1)
                        for i, (name, _) in enumerate(zones)}
        profile = export_profile(ref, c)
        profile.update(count=2, dtype="float32", nodata=c.NODATA, predictor=3)
        output_raster = out / "Infanta_Coherence_Recovery_Fraction_and_Class.tif"
        vrt_stack = [aligned(src, ref, c) for src in srcs]
        zone_geoms = [g for _, g in zones]
        with rasterio.open(output_raster, "w", **profile) as dst:
            dst.set_band_description(1, "end-state coherence recovery fraction (dimensionless; unclipped)")
            dst.set_band_description(2, "recovery class code")
            dst.update_tags(1, units="dimensionless", formula="(median(last 3 post coherence pairs) - event coherence) / (pre-event median - event coherence)",
                            valid_classes="1,2,3 only; band 1 is NoData for class 4 and invalid pixels",
                            end_state_post_pairs=c.ENDSTATE_POST_PAIRS,
                            nodata_meaning="invalid/unstable baseline, missing event, insufficient end-state pairs, or outside mask")
            dst.update_tags(2, units="class code", class_1="Recovered: final two consecutive post pairs meet pre median minus tolerance",
                            class_2="Partially recovered: detectable loss, not recovered, recovery fraction > 0",
                            class_3="Not recovered: detectable loss, not recovered, recovery fraction <= 0",
                            class_4="No detectable event-associated loss", nodata="-9999",
                            loss_threshold=f"max({c.LOSS_ABSOLUTE}, {c.LOSS_MAD_MULTIPLIER} * max(pre MAD sigma, {c.NOISE_FLOOR}))",
                            recovery_tolerance=c.RECOVERY_TOLERANCE,
                            recovery_confirmations=c.RECOVERY_CONFIRMATIONS,
                            reference_date=c.REFERENCE_DATE, event_origin_utc=c.EVENT_ORIGIN,
                            interpretation="Coherence response proxy; not proof of physical reconstruction")
            for win in windows(ref.width, ref.height, c.WINDOW_SIZE):
                stack = read_stack(vrt_stack, srcs, win, c)
                h, w = stack.shape[1:]
                transform = rasterio.windows.transform(win, ref.transform)
                mask = mask_for(zone_geoms, (h, w), transform, c.ALL_TOUCHED)
                stack[:, ~mask] = np.nan
                fm = calculate_fraction(stack, obs, c)
                recovery_fraction = fm.fraction.copy()
                recovery_class = fm.classes.astype("float32")
                recovery_fraction[fm.classes == 4] = c.NODATA
                recovery_fraction[fm.classes == -9999] = c.NODATA
                recovery_class[fm.classes == -9999] = c.NODATA
                dst.write(recovery_fraction.astype("float32"), 1, window=win)
                dst.write(recovery_class.astype("float32"), 2, window=win)
                complete = np.isfinite(stack).all(axis=0)
                if complete.any():
                    rows, cols = np.where(complete)
                    sample.add(stack[:, complete], rows + int(win.row_off), cols + int(win.col_off))
                    for (name, geom), zsample in zip(zones, zone_samples.values()):
                        zm = mask_for([geom], (h, w), transform, c.ALL_TOUCHED) & complete
                        if zm.any():
                            zr, zc = np.where(zm)
                            zsample.add(stack[:, zm], zr + int(win.row_off), zc + int(win.col_off))
                for (name, geom), acc in zip(zones, zone_acc.values()):
                    zm = mask_for([geom], (h, w), transform, c.ALL_TOUCHED)
                    acc.add(fm.fraction, fm.classes, zm)
        for vrt in vrt_stack:
            vrt.close()
        if sample.count:
            sample.table([o.name for o in obs]).to_csv(out / "timeseries_summary.csv", index=False)
        if zones:
            write_csv(out / "zonal_recovery_statistics.csv", [a.summary() for a in zone_acc.values()])
            # Zonal time series use fixed complete-case pixel support for comparability.
            for zone_index, ((name, _), s) in enumerate(zip(zones, zone_samples.values()), start=1):
                if s.count:
                    s.table([o.name for o in obs]).assign(zone=name).to_csv(
                        out / f"zone_timeseries_{zone_index:03d}.csv", index=False)
        # Sensitivity outputs count sustained returns under declared alternative thresholds.
        sensitivity = []
        for tolerance in c.SENSITIVITY_TOLERANCES:
            for k in c.SENSITIVITY_CONFIRMATIONS:
                if k <= 8 and sample.count:
                    sm = calculate_fraction(sample.values, obs, c, tolerance=tolerance, confirmations=k)
                    stat = sm.classes
                    assessed = np.isin(stat, [1, 2, 3]).sum()
                    sensitivity.append({"tolerance": tolerance, "confirmations": k,
                                        "sampled_loss_eligible": int(assessed),
                                        "sampled_recovered_fraction": float((stat == 1).sum()/assessed)
                                        if assessed else np.nan})
        write_csv(out / "recovery_sensitivity.csv", sensitivity)
        sample_fraction = calculate_fraction(sample.values, obs, c).fraction if sample.count else np.array([])
        portfolio(sample, {n: s for n, s in zip(zone_samples, zone_samples.values())},
                  obs, out, c, sample_fraction)
        # Dataset-level QA and reproducibility manifest.
        with rasterio.open(output_raster) as check:
            if check.crs != ref.crs or check.transform != ref.transform or check.shape != ref.shape:
                raise RuntimeError("Output georeferencing differs from the reference grid.")
            if check.count != 2 or check.dtypes != ("float32", "float32") or check.nodata != c.NODATA:
                raise RuntimeError("Output band count, dtype, or NoData failed verification.")
            counts = {str(v): 0 for v in (-9999, 1, 2, 3, 4)}
            for win in windows(check.width, check.height, c.WINDOW_SIZE):
                a = check.read(2, window=win)
                f = check.read(1, window=win)
                if not np.isin(a, [-9999, 1, 2, 3, 4]).all():
                    raise RuntimeError("Band 2 contains an undefined class code.")
                eligible_classes = np.isin(a, [1, 2, 3])
                if not np.isfinite(f).all() or not np.array_equal(f == c.NODATA, ~eligible_classes):
                    raise RuntimeError("Band 1 NoData locations do not match the class-band definition.")
                for code in counts:
                    counts[code] += int((a == int(code)).sum())
            manifest = {"input_count": 17, "output": output_raster.name,
                        "bands": {"1": "unclipped recovery fraction, classes 1-3 only", "2": "class code 1-4"},
                        "crs": str(check.crs), "transform": tuple(check.transform),
                        "shape": check.shape, "dtype": check.dtypes[0], "nodata": check.nodata,
                        "special_value_counts": counts, "settings": {
                            key: getattr(c, key) for key in ("REFERENCE_DATE", "EVENT_ORIGIN", "MIN_PRE",
                                "BASELINE_MIN", "BASELINE_MAD_MAX", "LOSS_ABSOLUTE", "RECOVERY_TOLERANCE",
                                "ENDSTATE_POST_PAIRS", "RECOVERY_CONFIRMATIONS", "MAX_CONFIRM_GAP_DAYS", "RESAMPLING")}}
            (out / "run_manifest.json").write_text(json.dumps(manifest, indent=2, default=str), encoding="utf-8")
        logging.info("Finished. Primary raster: %s", output_raster)
        logging.info("QA categories: %s", counts)
        return 0
    finally:
        for src in srcs:
            src.close()


def main() -> int:
    c = load_config(Path(__file__).resolve().parent / "config.py")
    if getattr(c, "RECOVERY_MODEL", "legacy_observed") == "legacy_observed":
        return legacy_main()
    if c.RECOVERY_MODEL != "itsa":
        raise ValueError("RECOVERY_MODEL must be itsa or legacy_observed.")
    from sarflow.research_pipeline import run
    return run(c)


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        logging.exception("Analysis failed: %s", exc)
        raise

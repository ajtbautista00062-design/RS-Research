"""Redraw observed series from a completed dated run; never process rasters."""
import argparse
import json
from pathlib import Path
import pandas as pd
from run_analysis import load_config
from sarflow.temporal import Observation, local_path, timestamp
from sarflow.figures import observed_series


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-directory", type=Path)
    args = parser.parse_args()
    c = load_config()
    root = local_path(c.Output_Folder, c.PROJECT_ROOT)
    available = sorted(p for p in root.glob("series_*") if p.is_dir())
    out = local_path(args.run_directory, c.PROJECT_ROOT) if args.run_directory else (available[-1] if available else None)
    if out is None:
        raise ValueError("No completed dated run is available.")
    if json.loads((out / "run_state.json").read_text())["status"] != "complete":
        raise ValueError("Selected run is not complete.")
    graphs = root / "Time Series Graphs" / out.name
    frame = pd.read_csv(out / "observation_dates.csv")
    obs = [Observation(row["name"], Path(row["path"]), row["phase"],
            timestamp(row["first_filename_date"]), timestamp(row["second_filename_date"]),
            timestamp(row["start"]), timestamp(row["end"]), timestamp(row["reference"]))
           for row in frame.to_dict("records")]
    dates = f"{min(o.start for o in obs):%Y%m%d}_{max(o.end for o in obs):%Y%m%d}"
    for support in ("available", "common"):
        csv = graphs / f"timeseries_{support}_support.csv"
        if csv.is_file():
            table = pd.read_csv(csv)
            observed_series(table, obs, graphs, c, f"Coherence_{support.title()}Support_{dates}",
                            f"{support.title()} spatial support; observed pairs only")


if __name__ == "__main__":
    main()


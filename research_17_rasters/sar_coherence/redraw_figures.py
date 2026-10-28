"""Regenerate journal figures from exact cached plotting inputs, without reprocessing rasters."""
import argparse
import hashlib
import json
import os
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
os.environ.setdefault('MPLCONFIGDIR', str(ROOT/'.mplconfig'))
import matplotlib
matplotlib.use('Agg')
import config as c
from sarflow import metrics
from sarflow.sampling import Sample
from sarflow.temporal import validate, local_path
from sarflow.journal.cache import read_cache
from sarflow.journal.portfolio import render_portfolio

STEMS = ('01_coherence_timeseries','02_coherence_anomaly','03_phase_distributions',
         '04_single_post_response','04_recovery_trajectory','05_zone_heatmap',
         '06_recovery_fraction_distribution','07_zone_trajectories',
         'fig_coherence_timeseries_20220925','fig_coherence_timeseries_20220925_single_column',
         'fig_coherence_synthesis_20220925')


def main(only_main=False):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--replace', action='store_true', help='Replace allowlisted figure exports only')
    parser.add_argument('--single-column', action='store_true')
    parser.add_argument('--check', action='store_true', help='Check configuration/cache without rendering')
    args = parser.parse_args()
    try:
        obs = validate(c)
        out = local_path(c.Output_Folder, ROOT)
        sample, zones, fraction = read_cache(out,c,obs,metrics,Sample)
    except (ValueError, FileNotFoundError) as exc:
        print(json.dumps({'status':'BLOCKED','reason':str(exc)},indent=2))
        return 2
    if args.check:
        print(json.dumps({'status':'READY','observations':len(obs),'cache':'publication_plot_cache.npz'}))
        return 0
    numerical = {p:hashlib.sha256(p.read_bytes()).hexdigest() for p in out.glob('*.csv')}
    if args.replace:
        stems = STEMS[:1]+STEMS[-3:-1] if only_main else STEMS
        for stem in stems:
            for suffix in ('.pdf','.svg','.png','_caption.md','_data.json'):
                path = (out/(stem+suffix)).resolve()
                if path.parent != out.resolve():
                    raise ValueError('Deletion target outside output directory')
                if path.is_file():
                    path.unlink()
    result = render_portfolio(sample,zones,obs,out,c,fraction,metrics,
                              only_main=only_main,single_column=args.single_column)
    assert all(hashlib.sha256(p.read_bytes()).hexdigest()==digest for p,digest in numerical.items())
    print(json.dumps({'status':'PASS','csv_files_unchanged':len(numerical),**result},indent=2))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())

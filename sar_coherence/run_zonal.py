"""Multi-scale zonal recovery classification for a completed study run.

Writes Zonal_* files into the run's existing A. RRI and B. ITSA folders.

Usage:  python run_zonal.py                          (latest completed run)
        python run_zonal.py --run vscode_run1
        python run_zonal.py --run vscode_run1 --replace   (regenerate Zonal_* files)
"""
import argparse
import json
import config
from sarflow.zonal_recovery import run_zonal


def latest_complete():
    runs = []
    for p in config.VERSION_ROOT.iterdir():
        state = p / 'run_state.json'
        if state.is_file() and json.loads(state.read_text(encoding='utf-8')).get('status') == 'complete':
            runs.append(p)
    if not runs:
        raise SystemExit('No completed run found; run run_study.py first.')
    return max(runs, key=lambda p: (p / 'run_state.json').stat().st_mtime)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--run', default=None, help='Run folder name under ' + str(config.VERSION_ROOT))
    parser.add_argument('--replace', action='store_true', help='Delete and regenerate existing Zonal_* files')
    args = parser.parse_args()
    root = config.VERSION_ROOT / args.run if args.run else latest_complete()
    out, overall = run_zonal(root, config, replace=args.replace)
    print('Zonal outputs written into A. RRI and B. ITSA of:', out)
    for k in ('population', 'I0_original_indeterminate', 'I1_zonal_indeterminate', 'indeterminate_reduction_percent',
              'C0_original_classes_1to3', 'C1_zonal_classes_1to3', 'recovery_class_increase_percent'):
        print(f'  {k}: {overall[k]}')


if __name__ == '__main__':
    main()

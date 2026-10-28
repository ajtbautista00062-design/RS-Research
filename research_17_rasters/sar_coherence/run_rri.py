"""Run separately AFTER a completed ITSA analysis; never refits the raster stack."""
import argparse
from pathlib import Path
import rri_config


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source',required=True,type=Path,help='Completed ITSA raster run directory')
    parser.add_argument('--manifest',type=Path,help='Matching run_manifest.json; inferred from CSV files/run ID by default')
    args=parser.parse_args()
    from sarflow.rri_pipeline import run
    run(args.source,args.manifest,rri_config)


if __name__=='__main__':
    main()

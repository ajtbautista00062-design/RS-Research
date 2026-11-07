"""Run the full ITSA/RRI recovery study using the rasters listed in config.py.

Steps: validate inputs -> align rasters -> pixel and PRE-object ITSA -> RRI ->
per-barangay exports and verification -> multi-scale zonal recovery product.

Usage:  python run_study.py --run-name run01
"""
import argparse
import config
from sarflow.study_pipeline import run


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--run-name', default=None,
                        help='New, short results subfolder name (existing runs are never replaced)')
    args = parser.parse_args()
    root = run(config, args.run_name)
    print('Run complete:', root)


if __name__ == '__main__':
    main()

# Active research implementation

`sar_coherence/` is the only implementation in this repository.

Entry points:
- `run_study.py`: full ITSA/RRI study plus multi-scale zonal recovery product
- `run_zonal.py`: zonal recovery product only, for an existing run
- `run_analysis.py`: alias of `run_study.py`

Users edit only `config.py` (Pre1..Pre8, Event_Straddler/During, Post1..Post8).
Read `sar_coherence/README.md` before changing code. Run the tests with
`python -m unittest discover -s tests -p "test_*.py"` from `sar_coherence/`.
Earlier versions are in the git history, not in the working tree.

"""Backwards-compatible main-figure command, using the shared publication system."""
from redraw_figures import main
if __name__ == '__main__':
    raise SystemExit(main(only_main=True))

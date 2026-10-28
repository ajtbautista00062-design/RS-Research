"""Journal graphics adapters; analytical modules remain unchanged."""
from pathlib import Path
import os
os.environ.setdefault("MPLCONFIGDIR", str(Path(__file__).resolve().parents[2]/".mplconfig"))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from . import metrics
from .journal.cache import write_cache
from .journal.portfolio import render_portfolio
from .journal.interactive import interactive


def save(fig, out, name, c):
    # Existing non-time-series categorical map export remains available.
    for ext in c.FIGURE_FORMATS:
        fig.savefig(out/f"{name}.{ext}", dpi=c.FIGURE_DPI, bbox_inches="tight")
    plt.close(fig)


def portfolio(sample, zones, obs, out, c, recovery_fraction):
    if not sample.count:
        import logging
        logging.warning("No common valid pixels; publication portfolio unavailable")
        return
    write_cache(out,sample,zones,obs,c,recovery_fraction,metrics)
    render_portfolio(sample,zones,obs,out,c,recovery_fraction,metrics)

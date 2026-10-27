"""Monochrome observed coherence and exploratory interrupted time series."""
from pathlib import Path
import os
os.environ.setdefault("MPLCONFIGDIR", str(Path(__file__).resolve().parents[2] / ".mplconfig"))
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.lines import Line2D
from .temporal import timestamp
from .itsa import design, fit_complete, MODEL, LIMITATION
from .pilot_statistics import export_statistics

MARKERS = {"pre": "o", "event": "D", "post": "s"}

def save(fig, out, name, c):
    for ext in c.FIGURE_FORMATS:
        fig.savefig(out / f"{name}.{ext}", dpi=c.FIGURE_DPI,
                    bbox_inches="tight", facecolor="white")
    plt.close(fig)

def observed_series(table, obs, out, c, stem, support_label):
    # Fit only a fixed population: varying coverage can produce spurious trends.
    fit = "CommonSupport" in stem or "common" in stem.lower()
    x = design(obs, c)
    y = table.median_sample.to_numpy(dtype=float)
    beta = fit_complete(y, x) if fit and np.isfinite(y[x[:,0] == 1]).all() else None
    dates = mdates.date2num([o.reference for o in obs])
    event = mdates.date2num(timestamp(c.EVENT_ORIGIN))
    with plt.rc_context({"font.family": "serif", "font.size": 10,
                         "axes.spines.top": False, "axes.spines.right": False,
                         "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none"}):
        fig, (ax, pair_ax) = plt.subplots(2, 1, figsize=(10, 6.5), sharex=True,
            gridspec_kw={"height_ratios": [4, 1.1]}, layout="constrained")
        ax.plot(dates, y, color="0.70", lw=.8, zorder=1)
        for j, o in enumerate(obs):
            row = table.iloc[j]
            ax.vlines(dates[j], row.q25_sample, row.q75_sample, color="0.65", lw=2)
            ax.plot(dates[j], y[j], marker=MARKERS[o.phase], ms=6,
                    color="#245A8D", mfc="white" if o.phase == "event" else "#245A8D")
            a, b = mdates.date2num([o.start, o.end])
            lane = {"pre": 2, "event": 1, "post": 0}[o.phase]
            pair_ax.plot([a, b], [lane, lane], color="0.25", lw=1.5)
            pair_ax.plot(dates[j], lane, marker=MARKERS[o.phase], color="black", ms=4)
        if beta is not None:
            pre = np.array([o.phase == "pre" for o in obs])
            post = np.array([o.phase == "post" for o in obs])
            for domain, after in (([dates[pre].min(), event], False),
                                  ([event, dates[post].max()], True)):
                grid = np.linspace(*domain, 100)
                t = grid-event
                predicted = beta[0]+beta[1]*t
                if after:
                    predicted += beta[2]+beta[3]*t
                ax.plot(grid, predicted, color="#245A8D", ls="--" if after else "-",
                        lw=1.8, zorder=3)
            grid = np.linspace(event, dates[post].max(), 100)
            ax.plot(grid, beta[0]+beta[1]*(grid-event), color="0.40", lw=1.3, ls=(0, (6, 3)), zorder=2)
            use, residual = export_statistics(y, x, beta,
                [o.reference.isoformat() for o in obs], out, stem)
            residual_fig, residual_ax = plt.subplots(figsize=(8,3), layout="constrained")
            residual_ax.axhline(0, color="#202020", lw=1)
            residual_ax.plot(dates[use], residual, "o-", color="#245A8D", lw=.8)
            residual_ax.set(xlabel="Observation date (UTC)", ylabel="Coherence residual",
                title=f"Exploratory residual diagnostics (n={int(use.sum())}, df={int(use.sum())-4})")
            residual_ax.xaxis.set_major_formatter(mdates.DateFormatter("%d %b %Y"))
            residual_fig.text(.5, -.03, "Short-series diagnostics have limited evidential value.",
                              ha="center", fontsize=8)
            save(residual_fig, out, stem+"_residuals", c)
        for axis in (ax, pair_ax):
            axis.axvspan(timestamp(c.EVENT_START), timestamp(c.EVENT_END), color="0.9", zorder=0)
            axis.axvline(event, color="black", ls="--", lw=1.2, zorder=2)
        # Do not hide unconstrained OLS predictions outside the physical data range.
        ax.set_ylim(min(0, ax.get_ylim()[0]), max(1, ax.get_ylim()[1]))
        ax.set_ylabel("Sentinel-1 coherence")
        ax.grid(axis="y", color="0.9", lw=.6)
        ax.set_title("Infanta, Quezon: exploratory interrupted time-series analysis" if beta is not None
                     else "Infanta, Quezon: observed coherence (available spatial support)",
                     fontsize=12, loc="left", pad=14)
        ax.annotate(
            "Karding / Noru interruption\n" +
            timestamp(c.EVENT_ORIGIN).strftime("%d %b %Y, %H:%M UTC"),
            xy=(event, .98), xycoords=("data", "axes fraction"),
            xytext=(6, -2), textcoords="offset points",
            ha="left", va="top", fontsize=8,
            bbox={"facecolor": "white", "edgecolor": "none", "alpha": .85})
        pair_ax.set_yticks([2, 1, 0], ["PRE", "EVENT-STRADDLING", "POST"])
        pair_ax.tick_params(axis="y", labelsize=8, length=0)
        pair_ax.set_ylim(-.5, 2.5)
        locator = mdates.AutoDateLocator(minticks=5, maxticks=9)
        pair_ax.xaxis.set_major_locator(locator)
        pair_ax.xaxis.set_major_formatter(mdates.ConciseDateFormatter(locator))
        pair_ax.set_xlabel(f"Acquisition intervals (UTC); markers at pair {c.REFERENCE_DATE}")
        handles = [Line2D([], [], color="#245A8D", marker=MARKERS[p], ls="none",
                          mfc="white" if p == "event" else "#245A8D",
                          label=p.upper() if p != "event" else "EVENT (excluded from fit)")
                   for p in MARKERS]
        handles.append(Line2D([], [], color="0.65", lw=2, label="Spatial IQR (not a confidence interval)"))
        handles += [Line2D([], [], color="black", ls="--", lw=1.2,
                           label="Interruption (configured event time)"),
                    Line2D([], [], color="0.70", lw=.8,
                           label="Observed medians (sample estimates)")]
        if beta is not None:
            handles += [
                Line2D([], [], color="#245A8D", lw=1.8, label="PRE regression trend"),
                Line2D([], [], color="#245A8D", lw=1.8, ls="--", label="POST regression trend"),
                Line2D([], [], color="0.40", ls=(0, (6, 3)),
                       label="Counterfactual (PRE-trend projection)")]
        fig.legend(handles=handles, loc="outside lower center", ncol=3, frameon=False, fontsize=8)
        save(fig, out, stem, c)
    caption = (support_label + ". Points are seeded-sample estimates of spatial medians; bars are sampled spatial IQR. "
               "Thin gray connectors are visual guides, not fitted trajectories. "
               "Lower bars show both acquisition dates. No observations are interpolated. "
               "Pair reference dates approximate interval measurements. Small enclosed missing areas were spatially estimated before analysis; see correction flags and run metadata. ")
    if beta is not None:
        caption += MODEL + ". " + LIMITATION + (
            " Blue solid/dashed lines show PRE/POST trends; gray dashed line is a modeled continuation "
            "of the PRE trend (counterfactual), not observed recovery. The unfilled event diamond is excluded from regression. "
            "The vertical dashed line marks the configured event origin. "
            "The POST/counterfactual gap is a model contrast, not an established causal effect. "
            "The level shift is extrapolated to the event origin; OLS predictions are not clipped.")
    else:
        caption += "No ITSA is fitted to changing spatial support."
    (out / f"{stem}_caption.md").write_text(caption+chr(10), encoding="utf-8")

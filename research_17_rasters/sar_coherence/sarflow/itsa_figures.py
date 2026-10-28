"""Formal observed-pair time plot and segmented fit on a fixed spatial sample."""
import numpy as np
import pandas as pd
from .recovery import assess
from .temporal import day_offsets


def render_itsa(sample, obs, folder, c):
    if not sample.count:
        import logging
        logging.warning("No complete common-pixel support: ITSA sample graph not generated; inspect coverage and per-date support.")
        return
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.dates as mdates
    series = np.median(sample.values, axis=1)
    q = np.quantile(sample.values, [.25, .75], axis=1)
    result = assess(series[:, None], obs, c)
    beta = result.fit.beta[:, 0]
    origin = pd.Timestamp(c.EVENT_ORIGIN)
    days = day_offsets(obs, c.EVENT_ORIGIN)
    phases = np.array([o.phase for o in obs])
    dates = [o.reference for o in obs]
    with plt.rc_context({"font.family": "DejaVu Sans", "font.size": 10,
                         "axes.spines.top": False, "axes.spines.right": False,
                         "axes.labelsize": 11, "legend.fontsize": 8}):
        fig, ax = plt.subplots(figsize=(10, 5.3), layout="constrained")
        for phase, marker in (("pre", "o"), ("event", "D"), ("post", "s")):
            indices = np.flatnonzero(phases == phase)
            for i in indices:
                ax.hlines(series[i], obs[i].start, obs[i].end, color=".65", lw=1)
            ax.errorbar([dates[i] for i in indices], series[indices],
                        yerr=np.vstack((series[indices]-q[0, indices], q[1, indices]-series[indices])),
                        marker=marker, color="black", mfc="white" if phase == "event" else "black",
                        linestyle="none", capsize=3,
                        label=f"{phase.upper()}: observed spatial median / IQR")
        if np.isfinite(beta).all():
            for phase, style in (("pre", "-"), ("post", "-")):
                selected = days[phases == phase]
                t = np.linspace(selected.min(), selected.max(), 100)
                model = beta[0] + beta[1]*t
                if phase == "post":
                    model += beta[2] + beta[3]*t
                ax.plot(origin + pd.to_timedelta(t, unit="D"), model, style, color="black",
                        lw=1.8, label=f"{phase.upper()} fitted ITSA segment")
            t = np.linspace(0, max(days), 100)
            ax.plot(origin + pd.to_timedelta(t, unit="D"), beta[0]+beta[1]*t,
                    "--", color=".4", label="PRE counterfactual (disturbance context)")
        baseline = float(result.metrics["baseline"][0])
        ax.axhline(baseline, color=".3", ls=":", label="Fixed observed PRE median")
        ax.axhline(baseline-c.RECOVERY_TOLERANCE, color=".5", ls="-.", label="Baseline return tolerance")
        ax.axvline(origin, color="black", ls="--", lw=1, label="Karding event origin")
        ax.set(title="Infanta coherence: interrupted time series and fixed-baseline recovery",
               xlabel="Acquisition-pair reference date (UTC; bars show pair intervals)",
               ylabel="Interferometric coherence (dimensionless)")
        ax.xaxis.set_major_formatter(mdates.ConciseDateFormatter(ax.xaxis.get_major_locator()))
        ax.grid(axis="y", color=".9", lw=.6)
        ax.legend(loc="upper center", bbox_to_anchor=(.5, -.20), ncol=2, frameon=False)
        fig.text(.01, -.17, f"Fixed common-pixel support: {sample.count:,}; sampled {sample.values.shape[1]:,}. "
                 "Spatial IQR is not model uncertainty. Event pair excluded from pure-phase OLS.\n"
                 "Aggregate fit is descriptive; pixel/object classifications use their own trajectories.", fontsize=8)
        for ext in c.FIGURE_FORMATS:
            fig.savefig(folder / f"ITSA_{obs[0].start:%Y%m%d}_{obs[-1].end:%Y%m%d}.{ext}",
                        dpi=c.FIGURE_DPI, bbox_inches="tight")
        plt.close(fig)
    pd.DataFrame({"date": dates, "phase": phases, "median": series, "spatial_q25": q[0],
                  "spatial_q75": q[1]}).to_csv(folder / "itsa_plot_observations.csv", index=False)

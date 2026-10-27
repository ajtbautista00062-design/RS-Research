"""Standalone coefficient rasters and a map layout from real output values."""
import numpy as np
import rasterio
from rasterio.enums import Resampling
from .itsa import COEFFICIENTS, DIAGNOSTICS
from .spatial import windows
from .recovery import LABELS, COLORS

def export_maps(source_path, rasters, graphs, obs, c, dates):
    from contextlib import ExitStack
    names = COEFFICIENTS+DIAGNOSTICS
    with ExitStack() as resources:
        src = resources.enter_context(rasterio.open(source_path))
        writers = []
        paths = []
        for j,name in enumerate(names):
            path = rasters/f"{name}_{dates}.partial.tif"
            dst = resources.enter_context(rasterio.open(path,"w",**(src.profile | dict(count=1))))
            dst.set_band_description(1,name)
            dst.update_tags(source_band=str(len(obs)+j+1),source_file=source_path.name,
                            units=src.tags(len(obs)+j+1).get("units",""))
            writers.append(dst)
            paths.append(path)
        for window in windows(src.width,src.height,c.WINDOW_SIZE):
            data=src.read(list(range(len(obs)+1,src.count+1)),window=window)
            for j,dst in enumerate(writers):
                dst.write(data[j],1,window=window)
    for path in paths:
        path.rename(path.with_name(path.name.replace(".partial.tif",".tif")))
    if not c.ENABLE_BASELINE_CLASSIFICATION:
        return
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap, BoundaryNorm
    from matplotlib.patches import Patch
    from .figures import save
    path=rasters/f"Infanta_Baseline_Recovery_Classes_{dates}.tif"
    with rasterio.open(path) as src:
        scale=max(src.width,src.height)/4000
        h,w=max(1,round(src.height/max(1,scale))),max(1,round(src.width/max(1,scale)))
        data=src.read(1,out_shape=(h,w),resampling=Resampling.nearest)
        bounds,crs=src.bounds,str(src.crs)
    cmap=ListedColormap([np.array(COLORS[i])/255 for i in range(6)])
    fig,ax=plt.subplots(figsize=(8,8),layout="constrained")
    ax.imshow(data,cmap=cmap,norm=BoundaryNorm(np.arange(-.5,6),6),
              extent=(bounds.left,bounds.right,bounds.bottom,bounds.top),interpolation="nearest")
    ax.set(xlabel=f"Easting ({crs})",ylabel=f"Northing ({crs})",
           title=f"Preliminary exploratory recovery classification â€” {len(obs)}-observation pilot")
    ax.ticklabel_format(style="plain",useOffset=False)
    fig.legend(handles=[Patch(facecolor=np.array(COLORS[i])/255,label=LABELS[i]) for i in LABELS],
               loc="outside lower center",ncol=2,fontsize=8,frameon=False)
    save(fig,graphs,f"pilot_recovery_map_{dates}",c)
    (graphs/f"pilot_recovery_map_{dates}_caption.md").write_text(
        "Operational baseline-coherence recovery, not physical reconstruction. "
        "Map preview uses nearest-neighbor downsampling; full-resolution TIFF is authoritative. "
        "NoData is transparent. Small-hole spatial estimates enter the analysis; see correction flags. "
        "Thresholds and acquisition comparability require validation.",
        encoding="utf-8")


    # Smooth continuous view is presentation only; numerical TIFF values stay unclipped.
    fraction_path=rasters/f"Infanta_Baseline_Recovery_Fraction_{dates}.tif"
    with rasterio.open(fraction_path) as src:
        fraction=src.read(1,out_shape=(h,w),masked=True,resampling=Resampling.bilinear)
    if fraction.count():
        from matplotlib.colors import LinearSegmentedColormap, Normalize
        ramp=LinearSegmentedColormap.from_list("recovery",["#b2182b","#ef8a62","#67a9cf","#2166ac"])
        fig,ax=plt.subplots(figsize=(8,8),layout="constrained")
        image=ax.imshow(fraction,cmap=ramp,norm=Normalize(-1,2),
            extent=(bounds.left,bounds.right,bounds.bottom,bounds.top),
            interpolation="bilinear")
        ax.set(xlabel=f"Easting ({crs})",ylabel=f"Northing ({crs})",
               title="Continuous baseline coherence recovery fraction")
        ax.ticklabel_format(style="plain",useOffset=False)
        bar=fig.colorbar(image,ax=ax,shrink=.7,extend="both")
        bar.set_label("Recovery fraction (dimensionless; display interpolation only)")
        save(fig,graphs,f"pilot_continuous_recovery_map_{dates}",c)

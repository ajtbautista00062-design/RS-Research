"""Deterministic object examples; actual pair support and fitted ITSA, no synthetic points."""
import numpy as np
import pandas as pd
from model import design


def example_trajectories(records,obs,out,c):
    from figures import save
    import matplotlib.pyplot as plt
    import matplotlib.dates as mdates
    _,_,time=design(obs,c.EVENT_ORIGIN)
    dates=mdates.date2num([o.reference for o in obs])
    event=mdates.date2num(pd.Timestamp(c.EVENT_ORIGIN))
    for key in sorted(records)[:6]:
        row=records[key]
        y=np.array([row[f'date_{j+1}_median'] for j in range(len(obs))])
        fig,ax=plt.subplots(figsize=(9,4.5),layout='constrained')
        for phase,marker in (('pre','o'),('event','D'),('post','s')):
            which=np.array([o.phase==phase for o in obs])
            ax.plot(dates[which],y[which],linestyle='none',marker=marker,color='black',
                    mfc='white' if phase=='event' else 'black',label=phase.upper())
        if np.isfinite(row['b0']):
            for phase in ('pre','post'):
                which=np.array([o.phase==phase for o in obs])
                t=np.linspace(time[which].min(),time[which].max(),60)
                fitted=row['b0']+row['b1']*t
                if phase=='post':
                    fitted+=row['b2']+row['b3']*t
                ax.plot(event+t,fitted,color='black',ls='-' if phase=='pre' else '--')
        ax.axhline(row['pre_baseline'],color='.3',lw=1,label='Fixed PRE median')
        ax.axhline(row['pre_baseline']-c.analysis_settings().RETURN_TOLERANCE,
                   color='.4',ls=':',label='Return threshold')
        ax.axvline(event,color='black',ls='-.',label='Noru / Karding')
        ax.set(ylim=(0,1),xlabel='Acquisition-pair reference date (UTC)',ylabel='Object median coherence',
               title=f"Exploratory object {key}: class {row['class_code']}; area {row['area_m2']:.0f} m²")
        locator=mdates.AutoDateLocator(minticks=4,maxticks=8)
        ax.xaxis.set_major_locator(locator); ax.xaxis.set_major_formatter(mdates.ConciseDateFormatter(locator))
        ax.legend(fontsize=8,ncol=3,frameon=False)
        ax.spines[['top','right']].set_visible(False)
        save(fig,out,f'Object_{key}_ITSA',c)

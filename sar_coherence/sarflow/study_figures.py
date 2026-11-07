"""Publication exports with full categorical legends and explicit analytical scope."""
import os
from pathlib import Path
os.environ.setdefault('MPLCONFIGDIR',str(Path(__file__).resolve().parents[1]/'.runtime_matplotlib'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.colors import ListedColormap, BoundaryNorm
from matplotlib.patches import Patch
import numpy as np
import pandas as pd
from .recovery import CLASSES
from .recovery_io import COLORS
from .rri import INTERPRETATION
from .temporal import day_offsets


def save(fig,directory,name,c):
    for ext in c.FIGURE_FORMATS:
        fig.savefig(directory/(name+'.'+ext),dpi=c.FIGURE_DPI,bbox_inches='tight',facecolor='white')
    plt.close(fig)


def itsa_graph(rows,result,obs,directory,title,c):
    days=day_offsets(obs,c.EVENT_ORIGIN)
    dates=[o.reference for o in obs]
    origin=pd.Timestamp(c.EVENT_ORIGIN)
    beta=result.fit.beta[:,0]
    med=np.array([r['median'] for r in rows])
    with plt.rc_context({'font.family':'DejaVu Sans','font.size':10,'pdf.fonttype':42,'svg.fonttype':'none'}):
        fig,ax=plt.subplots(figsize=(10,6),layout='constrained')
        for i,o in enumerate(obs):
            ax.plot([o.start,o.end],[med[i],med[i]],color='.7',lw=1)
            ax.vlines(o.reference,rows[i]['q25'],rows[i]['q75'],color='.6',lw=2)
            eligible=rows[i]['valid_fraction']>=c.MIN_PER_DATE_VALID_COVERAGE
            ax.plot(o.reference,med[i],marker='x' if not eligible else 'D' if o.phase=='event' else 'o' if o.phase=='pre' else 's',
                    color='#245a8d',mfc='white' if o.phase=='event' else '#245a8d')
        for phase in ('pre','post'):
            ix=np.array([o.phase==phase for o in obs])
            t=np.linspace(days[ix].min(),days[ix].max(),100)
            values=beta[0]+beta[1]*t+(beta[2]+beta[3]*t if phase=='post' else 0)
            ax.plot(origin+pd.to_timedelta(t,unit='D'),values,label=phase.upper()+' aggregate fit')
        if np.isfinite(beta).all():
            t=days[days>0]
            ax.plot(origin+pd.to_timedelta(t,unit='D'),beta[0]+beta[1]*t,':',color='.5',label='PRE counterfactual')
        baseline=float(result.metrics['baseline'][0])
        ax.axhline(baseline,color='black',lw=.8,label='Observed PRE median')
        ax.axhline(baseline-c.RECOVERY_TOLERANCE,color='.4',ls='-.',lw=.8,label='Recovery threshold')
        ax.axvline(origin,color='black',ls='--',lw=.8,label='Noru / Karding reference origin')
        ax.set(title=title+'\nSeparate descriptive fit to area medians',ylabel='Interferometric coherence',
               xlabel='Acquisition pair reference date (UTC)',ylim=(0,1))
        ax.xaxis.set_major_formatter(mdates.ConciseDateFormatter(ax.xaxis.get_major_locator()))
        ax.grid(axis='y',alpha=.15)
        ax.legend(fontsize=8,loc='upper center',bbox_to_anchor=(.5,-.16),ncol=3,frameon=False)
        fig.text(.02,-.14,'Event-straddling pair is displayed and excluded from OLS. Bars: spatial IQR, not confidence intervals.\n'
                 'Aggregate fit does not determine pixel or object recovery classes; source comparability remains unconfirmed.',fontsize=8)
        save(fig,directory,'ITSA_Time_Series',c)
        fig,ax=plt.subplots(figsize=(9,3.5),layout='constrained')
        use=np.array([o.phase!='event' for o in obs])
        x=np.column_stack((np.ones(len(obs)),days,np.array([o.phase=='post' for o in obs]),days*np.array([o.phase=='post' for o in obs])))
        ax.axhline(0,color='.4',lw=.8)
        use &= np.array([r['valid_fraction']>=c.MIN_PER_DATE_VALID_COVERAGE for r in rows])
        ax.plot(np.array(dates)[use],(med-x@beta)[use],'o-')
        ax.set(title=title+' aggregate residuals',ylabel='Coherence residual',xlabel='Pair reference date (UTC)')
        ax.xaxis.set_major_formatter(mdates.DateFormatter('%d %b'))
        save(fig,directory,'ITSA_Residuals',c)


def map_classes(array,member,directory,title,c,method='ITSA'):
    labels=CLASSES if method=='ITSA' else {0:'Undefined RRI',**{k:v for k,v in INTERPRETATION.items() if k}}
    colors=[np.array(COLORS[i])/255 for i in range(6)] if method=='ITSA' else ['#d8dde2','#c73a49','#8c78a3','#f4c95d','#167c80','#245a8d']
    fig,ax=plt.subplots(figsize=(8,7),layout='constrained')
    view=np.ma.array(array,mask=~member)
    ax.imshow(view,cmap=ListedColormap(colors),norm=BoundaryNorm(np.arange(-.5,6),6),interpolation='nearest')
    ax.set(title=title,xlabel='Reference-grid column offset',ylabel='Reference-grid row offset')
    fig.legend(handles=[Patch(facecolor=colors[i],label=labels[i]) for i in range(6)],
               loc='outside lower center',ncol=2,fontsize=8,frameon=False)
    save(fig,directory,method+'_Recovery_Map',c)


def rri_graph(values,classes,status,member,directory,title,c):
    valid=member&(status==1)&np.isfinite(values)
    fig,axes=plt.subplots(1,2,figsize=(12,5),layout='constrained')
    if valid.any():
        axes[0].hist(values[valid],bins=40,color='#245a8d',edgecolor='white')
    else:
        axes[0].text(.5,.5,'No eligible RRI values',ha='center',va='center',transform=axes[0].transAxes)
    for x in (-.05,.05,.95,1.05):
        axes[0].axvline(x,color='.5',ls=':',lw=.8)
    axes[0].set(xlabel='Final observation RRI (unclipped)',ylabel='Eligible raw pixels')
    counts=np.bincount(classes[valid],minlength=6)[1:]
    axes[1].bar(np.arange(5),counts,color=['#c73a49','#8c78a3','#f4c95d','#167c80','#245a8d'])
    axes[1].set_xticks(np.arange(5),['Further\ndegradation','No measurable\nrecovery','Partial\nrecovery','Baseline\nrecovery','Above baseline\nrecovery'],rotation=20,ha='right')
    axes[1].set(ylabel='Eligible raw pixels',xlabel='Independent RRI interpretation')
    fig.suptitle(title+f'\nRRI valid {valid.sum():,} of {member.sum():,} study pixels')
    save(fig,directory,'RRI_Distribution_and_Interpretation',c)

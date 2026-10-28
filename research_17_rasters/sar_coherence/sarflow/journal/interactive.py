"""Theme-free offline HTML using the same exact observed statistics."""
import numpy as np
from .style import DESIGN as D


def interactive(sample, obs, out, c):
    if not c.INTERACTIVE_HTML:
        return
    import logging
    try:
        from plotly.io import write_html
    except ImportError:
        logging.warning('Plotly unavailable: interactive HTML skipped')
        return
    q25, median, q75 = np.quantile(sample.values, [.25,.5,.75], axis=1)
    dates = [o.end.isoformat() for o in obs]
    symbols = [{'pre':'circle','event':'square','post':'triangle-up'}[o.phase] for o in obs]
    details = [[o.name, str(o.master), str(o.slave)] for o in obs]
    traces = [
        dict(type='scatter', x=dates, y=q25.tolist(), line={'width':0}, showlegend=False, hoverinfo='skip'),
        dict(type='scatter', x=dates, y=q75.tolist(), line={'width':0}, fill='tonexty',
             fillcolor='rgba(140,150,158,0.18)', name='Spatial IQR (P25–P75)', hoverinfo='skip'),
        dict(type='scatter', x=dates, y=(sample.sums/sample.count).tolist(), mode='lines+markers',
             name='Spatial mean', line={'color':D['mean'],'width':1.5},
             marker={'symbol':symbols,'size':7}, customdata=details,
             hovertemplate='%{customdata[0]}<br>Pair: %{customdata[1]} – %{customdata[2]}<br>γ=%{y:.6f}<extra></extra>'),
        dict(type='scatter', x=dates, y=median.tolist(), mode='lines+markers',
             name='Spatial median', line={'color':D['median'],'dash':'dash','width':1.2},
             marker={'symbol':[s+'-open' for s in symbols],'size':7})]
    span = float(q75.max()-q25.min())
    layout = {'paper_bgcolor':'white','plot_bgcolor':'white',
              'xaxis':{'type':'date','tickformat':'%Y-%m-%d','title':{'text':'Pair end date (UTC)'}},
              'yaxis':{'title':{'text':'Interferometric coherence (unitless)'},
                       'range':[max(0,float(q25.min())-.15*span),min(1,float(q75.max())+.15*span)],
                       'gridcolor':D['grid']},
              'legend':{'orientation':'h','y':-.25},'margin':{'b':105,'t':55},
              'shapes':[{'type':'line','xref':'x','yref':'paper','x0':D['event_date'],'x1':D['event_date'],
                         'y0':0,'y1':1,'line':{'color':D['reference'],'dash':'dash','width':1}}],
              'annotations':[{'x':D['event_date'],'y':1.1,'xref':'x','yref':'paper',
                              'text':'Typhoon Noru · 25 Sep 2022','showarrow':False,'font':{'size':11}}],
              'meta':{'limit':'no recovery trend inferred; expanded scale; theoretical range 0–1',
                      'common_pixels':sample.count,'sample_pixels':sample.values.shape[1]}}
    write_html({'data':traces,'layout':layout}, out/'interactive_timeseries.html',
               include_plotlyjs=True, validate=False)

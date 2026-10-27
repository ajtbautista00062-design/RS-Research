"""Synthetic in-memory plot rendering; no research graph or output files created."""
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import unittest
from unittest.mock import patch
import numpy as np
import pandas as pd
import config
from test_model_recovery import observations

class Sample:
    count=100
    def __init__(self,values):
        self.values=np.atleast_1d(values)
    def table(self,names):
        y=self.values
        return pd.DataFrame(dict(observation=names,complete_case_pixels=100,sample_pixels=100,
            mean_exact=y,median_sample=y,q25_sample=y-.02,q75_sample=y+.02))

class GraphTests(unittest.TestCase):
    def test_baseline_threshold_and_model_lines_render(self):
        import figures
        import matplotlib.pyplot as plt
        obs=observations()
        for j,o in enumerate(obs):
            o.name=f"SYNTHETIC_{j}"
        y=np.array([.7,.7,.7,.3,.5,.72,.94])
        rendered=[]
        def draw(fig,out,name,c):
            fig.canvas.draw()
            rendered.append((name,[line.get_label() for line in fig.axes[0].lines]))
            plt.close(fig)
        with patch.object(figures,"save",side_effect=draw),patch.object(pd.DataFrame,"to_csv"),patch.object(Path,"write_text"):
            figures.regional_graph(Sample(y),[Sample(v) for v in y],obs,config.EVENT_ORIGIN,
                config.analysis_settings(),Path(__file__).parent,config)
        self.assertEqual(len(rendered),2)
        self.assertIn("Fixed observed PRE median",rendered[0][1])
        self.assertIn("Baseline return threshold",rendered[0][1])

if __name__=="__main__":
    unittest.main()

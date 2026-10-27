"""Streaming comparison tables; no assumption that pixels are independent replicates."""
import json
from collections import defaultdict
import numpy as np
import pandas as pd
from recovery import LABELS,EXCLUSION,QUALITY,METRICS

class Summary:
    def __init__(self):
        self.inside=0; self.grid=0
        self.counts=defaultdict(lambda:np.zeros(6,dtype="int64"))
        self.reasons=np.zeros(max(EXCLUSION)+1,dtype="int64")
        self.quality={str(k):0 for k in QUALITY}
        self.transitions=defaultdict(lambda:np.zeros((6,6),dtype="int64"))
        self.filled=None; self.real_checks=0
        self.conditions=defaultdict(int); self.threshold_sensitive=0
    def add(self,name,classes,inside,primary=None):
        self.counts[name]+=np.bincount(classes[inside],minlength=6)
        if primary is not None:
            self.transitions[name]+=np.bincount((primary[inside]*6+classes[inside]),minlength=36).reshape(6,6)
    def export(self,out,metadata):
        out.mkdir(exist_ok=True)
        rows=[]
        for name,counts in self.counts.items():
            assessed=int(counts[1:5].sum()); eligible=int(counts[1:4].sum())
            rows.append(dict(scenario=name,inside_pixels=self.inside,assessed_pixels=assessed,
                recovered_fraction_eligible=counts[1]/eligible if eligible else np.nan,
                **{f"class_{i}":int(v) for i,v in enumerate(counts)},
                **{f"percent_{i}":100*v/self.inside for i,v in enumerate(counts)},
                changed_from_primary=int(self.transitions[name].sum()-np.trace(self.transitions[name])) if name in self.transitions else 0,
                independent_evidence_agreement="not available",false_positive_rate="not estimable"))
        pd.DataFrame(rows).to_csv(out/"classification_comparison.csv",index=False)
        pd.DataFrame([dict(comparison=k,primary_class=i,other_class=j,pixels=int(v))
            for k,m in self.transitions.items() for i,row in enumerate(m) for j,v in enumerate(row)]
            ).to_csv(out/"class_transition_counts.csv",index=False)
        pd.DataFrame([dict(code=i,reason=EXCLUSION[i],pixels=int(n)) for i,n in enumerate(self.reasons)]
            ).to_csv(out/"exclusion_counts_inside_boundary.csv",index=False)
        pd.DataFrame([dict(condition=k,accepted=n,rejected=self.inside-n)
                      for k,n in self.conditions.items()]).to_csv(out/"detection_conditions.csv",index=False)
        result=dict(metadata,grid_pixels=self.grid,inside_pixels=self.inside,outside_pixels=self.grid-self.inside,
                    independent_real_pixel_checks=self.real_checks,
                    threshold_sensitive_pixels=self.threshold_sensitive,
                    threshold_sensitive_percent=100*self.threshold_sensitive/self.inside,
                    quality_bit_counts=self.quality,spatial_filled_counts=self.filled,
                    classes=LABELS)
        (out/"validation_summary.json").write_text(json.dumps(result,indent=2,default=str),encoding="utf-8")
        text = [
            "# ITSA-derived baseline recovery run",
            "Class decisions use fitted trajectories and observed confirmation. "
            "Comparisons share the primary boundary/grid; masks and interpolation "
            "are explicitly labeled sensitivity assumptions.",
            f"Independent real-pixel checks: {self.real_checks}. "
            f"Grid pixels: {self.grid}; boundary pixels: {self.inside}.",
            "Zero semantics, comparability and baseline uncertainty remain governed by "
            "the manifest. Operational guards are not calibrated confidence intervals. "
            "Pixel counts are not independent-sample statistical tests. "
            "See METHOD.md and TECHNICAL_REPORT.md copied with this run."]
        (out/"runtime_report.md").write_text((chr(10)*2).join(text),encoding="utf-8")

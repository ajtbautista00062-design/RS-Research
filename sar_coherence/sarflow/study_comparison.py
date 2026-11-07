"""Fixed-population before/after comparisons, with all six class codes retained."""
from collections import Counter
import numpy as np
import pandas as pd
from .recovery import CLASSES, QUALITY, REASONS


class Comparison:
    def __init__(self):
        self.transitions = Counter()
        self.by_reason = Counter()
        self.gates = Counter()
        self.quality = Counter()
        self.reasons = Counter()
        self.sensitivity = Counter()
        self.rri_transitions = Counter()

    def add_rri(self,old,new,inside,owner):
        for zone in [0,*np.unique(owner[inside]).tolist()]:
            mask=inside if zone==0 else inside&(owner==zone)
            counts=np.bincount(old[mask].astype('int64')*7+new[mask],minlength=49).reshape(7,7)
            for i,j in np.ndindex(7,7):
                self.rri_transitions[(zone,i,j)]+=int(counts[i,j])

    def add(self, mode, original, corrected, inside, owner):
        for zone in [0, *np.unique(owner[inside]).tolist()]:
            member = inside if zone == 0 else inside & (owner == zone)
            pair = original.classes[member].astype('int64')*6+corrected.classes[member]
            counts = np.bincount(pair, minlength=36).reshape(6,6)
            for i,j in np.ndindex(6,6):
                self.transitions[(mode,zone,i,j)] += int(counts[i,j])
            triplet = np.stack([original.reasons[member], original.classes[member], corrected.classes[member]],axis=1)
            unique, numbers = np.unique(triplet,axis=0,return_counts=True)
            for (reason,old,new),n in zip(unique,numbers):
                self.by_reason[(mode,zone,int(reason),int(old),int(new))] += int(n)
            for code in REASONS:
                self.reasons[(mode,zone,code)] += int((corrected.reasons[member] == code).sum())
            for name,gate in corrected.gates.items():
                self.gates[(mode,zone,name,True)] += int(gate[member].sum())
                self.gates[(mode,zone,name,False)] += int((~gate[member]).sum())
            newly = member & (original.classes == 5) & np.isin(corrected.classes,(1,2,3,4))
            for bit in QUALITY:
                self.quality[(mode,zone,bit,'all')] += int(((corrected.quality[member]&bit)>0).sum())
                self.quality[(mode,zone,bit,'newly_resolved')] += int(((corrected.quality[newly]&bit)>0).sum())

    def add_sensitivity(self, name, classes, member, owner):
        for zone in [0,*np.unique(owner[member]).tolist()]:
            mask = member if zone == 0 else member & (owner == zone)
            counts = np.bincount(classes[mask],minlength=6)
            for code,n in enumerate(counts):
                self.sensitivity[(name,zone,code)] += int(n)

    def export(self, directory, names):
        directory.mkdir(parents=True,exist_ok=True)
        rows = [dict(mode=m,zone=names[z],original=i,corrected=j,pixels=n)
                for (m,z,i,j),n in self.transitions.items()]
        pd.DataFrame(rows).to_csv(directory/'Original_to_Corrected_Transitions.csv',index=False)
        summaries = []
        for mode,zone in sorted({(m,z) for m,z,_,_ in self.transitions}):
            matrix = np.array([[self.transitions[(mode,zone,i,j)] for j in range(6)] for i in range(6)])
            old,new = matrix.sum(1),matrix.sum(0)
            i0,i1,c0,c1 = int(old[5]),int(new[5]),int(old[1:4].sum()),int(new[1:4].sum())
            reduction = 100*(i0-i1)/i0 if i0 else None
            increase = 100*(c1-c0)/c0 if c0 else None
            summaries.append(dict(mode=mode,zone=names[zone],population=int(old.sum()),I0=i0,I1=i1,C0=c0,C1=c1,
                                  indeterminate_reduction_percent=reduction,resolved_recovery_increase_percent=increase,
                                  target_indeterminate_met=reduction is not None and reduction>=90,
                                  target_recovery_met=increase is not None and increase>=90,
                                  **{f'original_class_{k}':int(n) for k,n in enumerate(old)},
                                  **{f'corrected_class_{k}':int(n) for k,n in enumerate(new)}))
        pd.DataFrame(summaries).to_csv(directory/'Target_Assessment.csv',index=False)
        pd.DataFrame([dict(zone=names[z],original_status=i,corrected_status=j,pixels=n)
                      for (z,i,j),n in self.rri_transitions.items()]).to_csv(directory/'RRI_Validity_Transitions.csv',index=False)
        pd.DataFrame([dict(mode=m,zone=names[z],original_reason=r,original=i,corrected=j,pixels=n)
                      for (m,z,r,i,j),n in self.by_reason.items()]).to_csv(directory/'Changes_by_Original_Reason.csv',index=False)
        pd.DataFrame([dict(mode=m,zone=names[z],gate=g,accepted=a,pixels=n)
                      for (m,z,g,a),n in self.gates.items()]).to_csv(directory/'All_Gate_Outcomes.csv',index=False)
        pd.DataFrame([dict(mode=m,zone=names[z],quality_bit=b,label=QUALITY[b],population=p,pixels=n)
                      for (m,z,b,p),n in self.quality.items()]).to_csv(directory/'Quality_Status.csv',index=False)
        pd.DataFrame([dict(mode=m,zone=names[z],reason=r,label=REASONS[r],pixels=n)
                      for (m,z,r),n in self.reasons.items()]).to_csv(directory/'Remaining_Reasons.csv',index=False)
        pd.DataFrame([dict(scenario=s,zone=names[z],code=k,label=CLASSES[k],pixels=n)
                      for (s,z,k),n in self.sensitivity.items()]).to_csv(directory/'Sensitivity.csv',index=False)
        return summaries

"""Classical OLS diagnostics, explicitly conditional on unverified error assumptions."""
import numpy as np
import pandas as pd
from scipy.stats import t as student_t
from .itsa import COEFFICIENTS, LIMITATION

def export_statistics(y, x, beta, dates, out, stem):
    use = x[:,0] == 1
    xx, yy = x[use], y[use]
    fitted = xx@beta
    residual = yy-fitted
    n, k = xx.shape
    df = n-k
    sse = float(residual@residual)
    sst = float(np.sum((yy-yy.mean())**2))
    cov = (sse/df)*np.linalg.pinv(xx.T@xx)
    se = np.sqrt(np.maximum(0, np.diag(cov)))
    tval = np.divide(beta, se, out=np.full(k,np.nan), where=se > 0)
    ci = student_t.ppf(.975, df)*se
    table = pd.DataFrame(dict(term=COEFFICIENTS, estimate=beta, standard_error=se,
        ci95_lower=beta-ci, ci95_upper=beta+ci, p_value=2*student_t.sf(np.abs(tval),df),
        n=n, residual_df=df, r_squared=1-sse/sst if sst>1e-20 else np.nan,
        residual_rmse=np.sqrt(sse/n), uncertainty="Classical OLS; exploratory, error independence unverified"))
    table.to_csv(out/f"{stem}_coefficients.csv",index=False)
    fitted_all = np.full(len(y),np.nan)
    residual_all = fitted_all.copy()
    fitted_all[use], residual_all[use] = fitted,residual
    pd.DataFrame(dict(reference_date=dates, observed_median=y, included_primary_itsa=use,
                      fitted_median=fitted_all,residual=residual_all)).to_csv(out/f"{stem}_fit.csv",index=False)
    return use, residual

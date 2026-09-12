"""Two-stage robust anomaly detection.

Identical in structure to SSG-LUGIA: an Elliptic Envelope (Minimum Covariance Determinant
+ Mahalanobis distance) is fitted to all windows, then re-fitted on the windows the first
stage called native.  Cascading the two stages recovers alien windows whose signal was
masked by the contamination of the first fit.
"""
from __future__ import annotations

import numpy as np
from sklearn.covariance import EllipticEnvelope


def two_stage_envelope(X: np.ndarray, params: dict, random_state: int = 3):
    """Return ``(labels, scores)``; label -1 / score < 0 means *alien*.

    ``fit_stride`` (default 1, i.e. off) fits the robust covariance on every *k*-th window
    and then scores every window.  Consecutive windows overlap by 99% at the standard
    (10 kb, 100 bp) geometry, so decimating the *fitting* sample costs almost no
    information while making the Minimum-Covariance-Determinant search several times
    cheaper.  Scoring is unaffected: every window is always scored.
    """
    stride = int(params.get("fit_stride", 1))
    Xf = X[::stride] if stride > 1 else X

    e1 = EllipticEnvelope(contamination=params["contamination_model1"],
                          support_fraction=params["support_fraction_model1"],
                          random_state=random_state).fit(Xf)
    yp = e1.predict(X)
    ys = e1.decision_function(X)

    inlier = np.flatnonzero(yp == 1)
    if len(inlier) > X.shape[1] + 2:
        Xi = X[inlier]
        e2 = EllipticEnvelope(contamination=params["contamination_model2"],
                              support_fraction=params["support_fraction_model2"],
                              random_state=random_state).fit(Xi[::stride] if stride > 1 else Xi)
        ys[inlier] = e2.decision_function(Xi)
        yp[inlier] = e2.predict(Xi)
    return yp, ys

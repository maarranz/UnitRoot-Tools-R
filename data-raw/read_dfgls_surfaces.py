"""Minimal reader for inst/fixtures/dfgls_response_surfaces.json (Decision D11).

A reference for other implementations (e.g. the Python package UnitRoot-Tools):
evaluates DF-GLS critical values and p-values exactly as R's unitroottools
does, and checks the reference values stored in the file.
Run from the repository root:  python3 data-raw/read_dfgls_surfaces.py
"""
import json
from statistics import NormalDist

import numpy as np

SURFACES = json.load(open("inst/fixtures/dfgls_response_surfaces.json"))
_N = NormalDist()


def quantiles(trend, T, k):
    """Fitted quantiles at all tabulated probabilities; no extrapolation."""
    t = SURFACES["trends"][trend]
    kk = min(k, max(t["k"]))
    T = max(T, t["design"]["Tmin"][t["design"]["k"].index(kk)])
    slopes = np.array(t["slopes"][t["k"].index(kk)])
    return np.array(t["intercept"]) + slopes @ np.array([1 / T, 1 / T**2, 1 / T**3])


def critical_values(trend, T, k):
    """Critical values at the 1%, 5% and 10% levels."""
    q = quantiles(trend, T, k)
    p = np.array(SURFACES["trends"][trend]["probs"])
    return [float(q[np.argmin(abs(p - a))]) for a in (0.01, 0.05, 0.10)]


def p_value(stat, trend, T, k, npoints=11):
    """MacKinnon (1996): local cubic of qnorm(p) on the nearest fitted quantiles."""
    q = quantiles(trend, T, k)
    p = np.array(SURFACES["trends"][trend]["probs"])
    idx = np.argsort(abs(q - stat), kind="stable")[:npoints]
    qq = q[idx]
    X = np.column_stack([np.ones(npoints), qq, qq**2, qq**3])
    z = np.array([_N.inv_cdf(v) for v in p[idx]])
    b = np.linalg.lstsq(X, z, rcond=None)[0]
    return min(max(_N.cdf(float(b @ np.array([1, stat, stat**2, stat**3]))), 0.0), 1.0)


if __name__ == "__main__":
    worst = 0.0
    for r in SURFACES["reference_values"]["rows"]:
        cv = critical_values(r["trend"], r["T"], r["k"])
        worst = max(worst, abs(cv[0] - r["cv01"]), abs(cv[1] - r["cv05"]),
                    abs(cv[2] - r["cv10"]),
                    abs(p_value(r["stat"], r["trend"], r["T"], r["k"]) - r["p"]))
    print(f"max |difference| from R's reference values: {worst:.2e}")

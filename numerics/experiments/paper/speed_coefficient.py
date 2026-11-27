"""Table 1 of the article (label tab:speedcoef): first-order correction of the
flame speed.

Compares (c_eps - c_0)/eps, with c_eps the exact speed of the travelling wave
of (P_eps) (ODE shooting, combustion1d.exact.tw_speed_eps), against the
prediction of Proposition prop:speed (Proposition 2.2 in the draft in paper/;
the header printed by this script says "Prop 2.1", its number in an earlier
draft),

    c_eps = c_0 [1 + (eps/A)(1 - kappa (1 + w0))] + O(eps^2),
    kappa = (3 sqrt(3) - 2)/5          (for f(s) = (1 - s)^+),

for A = 1, w0 in {-0.5, 0, 1} and eps in {0.1, 0.05, 1e-2, 1e-3, 1e-4}.

Output: prints the table and writes results/paper/speed_coefficient.json =
{"A", "kappa", "rows": [{"w0", "c0", "prediction", "coef": {eps: value}}]}.

Run from numerics/:  python -m experiments.paper.speed_coefficient   (about 5 s)
"""

import numpy as np

from combustion1d.exact import tw_speed_eps, tw_speed_limit

from .common import save_json

A = 1.0
W0S = (-0.5, 0.0, 1.0)
EPS = (0.1, 0.05, 1e-2, 1e-3, 1e-4)
KAPPA = (3.0 * np.sqrt(3.0) - 2.0) / 5.0


def main():
    rows = []
    print(f"{'w0':>6} {'Prop 2.1':>11}" + "".join(f"{'eps=%g' % e:>12}" for e in EPS))
    for w0 in W0S:
        c0 = tw_speed_limit(A, w0)
        pred = c0 / A * (1.0 - KAPPA * (1.0 + w0))
        vals = [(tw_speed_eps(A, w0, e) - c0) / e for e in EPS]
        rows.append({"w0": w0, "c0": c0, "prediction": pred,
                     "coef": dict(zip(map(str, EPS), vals))})
        print(f"{w0:>6} {pred:>+11.6f}" + "".join(f"{v:>+12.6f}" for v in vals))
    save_json("speed_coefficient.json", {"A": A, "kappa": KAPPA, "rows": rows})


if __name__ == "__main__":
    main()

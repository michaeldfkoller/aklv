"""
PrivateCredit.py
================
Markov-chain valuation of private credit instruments using MarkovLib (MARKOVLV).

State Space
-----------
  S* ∪ {D}

where
  S* = {AAA, AA, A, BBB, BB, B, CCC, C}   (m = 8 S&P rating categories)
  D  = absorbing default state

Rating migrations follow the S&P long-run average annual transition matrix.
For intra-year frequencies the per-period matrix is obtained via eigendecomposition:
  P_period = P_annual^{1/freq}

Benefit structure (Convention A, consistent with MARKOVLV backward recursion)
------------------------------------------------------------------------------
  Cash-pay:
    t = 0 … T-1 : dSetPre(t,s,s) = nominal * r        [T coupons]
    t = T        : dSetPre(T,s,s) = nominal             [principal only]
    default      : dSetPost(t,s,D) = outstanding(t)*(1-LGD)

  PIK:
    t = 0 … T-1 : no cash benefit
    t = T        : dSetPre(T,s,s) = nominal*(1+r)^T    [compounded]
    default at t : dSetPost(t,s,D) = nominal*(1+r)^t*(1-LGD)

Extension (multi-tranche)
--------------------------
  ext_schedule = (f0, f1, ..., fK)  must sum to 1.0
    f0  fraction paid at T
    f1  fraction paid at T+1
    ...
    fK  fraction paid at T+K

  Example: (0.40, 0.40, 0.20) → 40% at T, 40% at T+1, 20% at T+2
"""

import os, sys
import numpy as np
from typing import List, Optional, Tuple, Union

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from markovlv import MARKOVLV, iSetLicense

_LICENSE_PATH = "./eth_summer_2027.lic"
_license_loaded = False


def _ensure_license():
    global _license_loaded
    if not _license_loaded:
        lic = open(_LICENSE_PATH).read().strip()
        if iSetLicense(lic) != 1:
            raise RuntimeError("markovlv: invalid or expired license")
        _license_loaded = True


# ---------------------------------------------------------------------------
# S&P long-run average annual transition matrix (1981-2023)
# Rows / cols: AAA, AA, A, BBB, BB, B, CCC, C, D
# ---------------------------------------------------------------------------

RATING_LABELS: List[str] = ["AAA", "AA", "A", "BBB", "BB", "B", "CCC", "C"]
M_STATES: int = len(RATING_LABELS)   # 8 live states
D_IDX:    int = M_STATES              # index 8 = Default

SP_ANNUAL: List[List[float]] = [
    # AAA    AA      A      BBB     BB      B      CCC     C       D
    [0.8744, 0.0834, 0.0060, 0.0006, 0.0010, 0.0000, 0.0000, 0.0000, 0.0000],
    [0.0058, 0.8825, 0.0783, 0.0069, 0.0008, 0.0006, 0.0002, 0.0000, 0.0001],
    [0.0005, 0.0210, 0.8822, 0.0541, 0.0072, 0.0026, 0.0001, 0.0000, 0.0006],
    [0.0002, 0.0021, 0.0389, 0.8624, 0.0547, 0.0103, 0.0017, 0.0000, 0.0026],
    [0.0003, 0.0007, 0.0019, 0.0549, 0.7561, 0.1065, 0.0098, 0.0000, 0.0108],
    [0.0000, 0.0007, 0.0023, 0.0032, 0.0558, 0.7352, 0.0544, 0.0000, 0.0522],
    [0.0000, 0.0000, 0.0027, 0.0039, 0.0129, 0.0934, 0.5158, 0.0000, 0.2787],
    [0.0000, 0.0000, 0.0010, 0.0020, 0.0080, 0.0600, 0.2000, 0.4500, 0.2790],
    [0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 1.0000],
]

# Historical year matrices used for bootstrap (AAA..CCC + D, 8 states)
SP_BOOTSTRAP: dict = {
    2000: [[0.8744,0.0834,0.0060,0.0006,0.0010,0.0000,0.0000,0.0000],
           [0.0058,0.8825,0.0783,0.0069,0.0008,0.0006,0.0002,0.0001],
           [0.0005,0.0210,0.8822,0.0541,0.0072,0.0026,0.0001,0.0006],
           [0.0002,0.0021,0.0389,0.8624,0.0547,0.0103,0.0017,0.0026],
           [0.0003,0.0007,0.0019,0.0549,0.7561,0.1065,0.0098,0.0108],
           [0.0000,0.0007,0.0023,0.0032,0.0558,0.7352,0.0544,0.0522],
           [0.0000,0.0000,0.0027,0.0039,0.0129,0.0934,0.5158,0.2787],
           [0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,1.0000]],
    2001: [[0.8700,0.0880,0.0070,0.0010,0.0005,0.0000,0.0000,0.0000],
           [0.0050,0.8770,0.0820,0.0090,0.0010,0.0005,0.0002,0.0002],
           [0.0004,0.0190,0.8780,0.0580,0.0080,0.0030,0.0001,0.0010],
           [0.0002,0.0018,0.0350,0.8560,0.0600,0.0120,0.0020,0.0040],
           [0.0002,0.0006,0.0015,0.0510,0.7430,0.1100,0.0110,0.0140],
           [0.0000,0.0005,0.0020,0.0030,0.0530,0.7200,0.0580,0.0680],
           [0.0000,0.0000,0.0020,0.0030,0.0110,0.0890,0.4900,0.3550],
           [0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,1.0000]],
    2004: [[0.8800,0.0800,0.0060,0.0005,0.0005,0.0000,0.0000,0.0000],
           [0.0060,0.8900,0.0760,0.0060,0.0007,0.0005,0.0001,0.0001],
           [0.0005,0.0220,0.8880,0.0510,0.0065,0.0022,0.0001,0.0003],
           [0.0002,0.0023,0.0400,0.8700,0.0510,0.0095,0.0015,0.0012],
           [0.0003,0.0007,0.0021,0.0570,0.7700,0.1020,0.0080,0.0070],
           [0.0000,0.0008,0.0025,0.0035,0.0580,0.7500,0.0520,0.0360],
           [0.0000,0.0000,0.0030,0.0042,0.0140,0.0970,0.5500,0.2350],
           [0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,1.0000]],
    2008: [[0.8620,0.0890,0.0120,0.0020,0.0010,0.0001,0.0000,0.0000],
           [0.0045,0.8620,0.0870,0.0130,0.0020,0.0010,0.0003,0.0010],
           [0.0004,0.0180,0.8650,0.0700,0.0110,0.0040,0.0002,0.0020],
           [0.0002,0.0016,0.0320,0.8420,0.0700,0.0150,0.0025,0.0060],
           [0.0002,0.0005,0.0013,0.0450,0.7200,0.1250,0.0160,0.0240],
           [0.0000,0.0004,0.0016,0.0025,0.0470,0.6900,0.0700,0.1050],
           [0.0000,0.0000,0.0018,0.0027,0.0100,0.0800,0.4500,0.4200],
           [0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,1.0000]],
    2009: [[0.8580,0.0910,0.0120,0.0020,0.0010,0.0001,0.0000,0.0000],
           [0.0042,0.8600,0.0880,0.0140,0.0022,0.0012,0.0003,0.0012],
           [0.0004,0.0170,0.8600,0.0740,0.0130,0.0048,0.0002,0.0028],
           [0.0002,0.0014,0.0300,0.8350,0.0740,0.0180,0.0030,0.0090],
           [0.0002,0.0004,0.0011,0.0400,0.7050,0.1350,0.0200,0.0350],
           [0.0000,0.0003,0.0014,0.0022,0.0440,0.6750,0.0750,0.1250],
           [0.0000,0.0000,0.0015,0.0024,0.0090,0.0750,0.4200,0.4700],
           [0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,1.0000]],
    2020: [[0.8650,0.0910,0.0110,0.0020,0.0010,0.0001,0.0000,0.0000],
           [0.0048,0.8680,0.0870,0.0130,0.0020,0.0010,0.0002,0.0008],
           [0.0004,0.0185,0.8700,0.0680,0.0105,0.0038,0.0001,0.0018],
           [0.0002,0.0016,0.0330,0.8480,0.0680,0.0155,0.0025,0.0075],
           [0.0002,0.0005,0.0013,0.0470,0.7280,0.1210,0.0155,0.0280],
           [0.0000,0.0004,0.0016,0.0026,0.0490,0.7020,0.0680,0.0950],
           [0.0000,0.0000,0.0017,0.0026,0.0095,0.0780,0.4600,0.4360],
           [0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,1.0000]],
    2023: [[0.8740,0.0836,0.0061,0.0007,0.0006,0.0001,0.0000,0.0000],
           [0.0057,0.8820,0.0786,0.0070,0.0009,0.0006,0.0002,0.0002],
           [0.0005,0.0212,0.8820,0.0542,0.0073,0.0026,0.0001,0.0006],
           [0.0002,0.0021,0.0382,0.8615,0.0550,0.0105,0.0017,0.0028],
           [0.0003,0.0007,0.0018,0.0544,0.7545,0.1068,0.0097,0.0115],
           [0.0000,0.0007,0.0022,0.0032,0.0554,0.7330,0.0545,0.0550],
           [0.0000,0.0000,0.0026,0.0038,0.0124,0.0922,0.5140,0.2900],
           [0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,1.0000]],
}
BOOTSTRAP_YEARS = sorted(SP_BOOTSTRAP.keys())


def _normalise(mat: np.ndarray) -> np.ndarray:
    mat = np.clip(mat, 0.0, None)
    rs  = mat.sum(axis=1, keepdims=True)
    return mat / np.where(rs == 0, 1.0, rs)


def _fractional_matrix(annual: np.ndarray, power: float) -> np.ndarray:
    vals, vecs = np.linalg.eig(annual)
    vals_p = np.power(vals.astype(complex), power).real
    result = (vecs @ np.diag(vals_p) @ np.linalg.inv(vecs)).real
    return _normalise(result)


# ---------------------------------------------------------------------------
# PrivateCredit
# ---------------------------------------------------------------------------

class PrivateCredit:
    """
    Private credit bullet loan valued via Markov-chain (MARKOVLV).

    Parameters
    ----------
    nominal       : face value
    coupon_rate   : ANNUAL coupon rate (decimal)
    n_years       : maturity in years
    freq          : payment frequency per year (1=annual, 4=quarterly)
    init_rating   : starting rating string ("BBB") or 0-based index
    lgd           : Loss Given Default (decimal)
    rf            : annual risk-free discount rate (decimal)
    pik           : True → PIK mode (no cash coupon; compounding at maturity)
    extension     : True → maturity extension active; see ext_schedule
    ext_schedule  : tuple of fractions (f0, f1, ..., fK) summing to 1.
                    f0 paid at T, f1 at T+1, ..., fK at T+K.
                    Default when extension=True: (0.40, 0.40, 0.20).
    trans_matrix  : optional (M+1)×(M+1) ANNUAL matrix to override S&P long-run.
    """

    _DEFAULT_EXT_SCHEDULE = (0.40, 0.40, 0.20)

    def __init__(
        self,
        nominal:      float = 100_000.0,
        coupon_rate:  float = 0.04,
        n_years:      int   = 5,
        freq:         int   = 1,
        init_rating:  Union[str, int] = "BBB",
        lgd:          float = 0.40,
        rf:           float = 0.045,
        pik:          bool  = False,
        extension:    bool  = False,
        ext_schedule: Optional[Tuple[float, ...]] = None,
        trans_matrix: Optional[List[List[float]]] = None,
    ):
        _ensure_license()

        self.nominal     = float(nominal)
        self.annual_rate = float(coupon_rate)
        self.n_years     = int(n_years)
        self.freq        = int(freq)
        self.lgd         = float(lgd)
        self.rf_annual   = float(rf)
        self.pik         = bool(pik)
        self.extension   = bool(extension)

        # Per-period rates
        self.r = coupon_rate / freq
        self.d = (1.0 + rf) ** (1.0 / freq) - 1

        # Total coupon periods
        self.T = n_years * freq

        # Rating index
        if isinstance(init_rating, str):
            self.init_state = RATING_LABELS.index(init_rating)
        else:
            self.init_state = int(init_rating)

        # Extension schedule
        if not self.extension:
            self._ext = (1.0,)
        elif ext_schedule is not None:
            self._ext = tuple(float(f) for f in ext_schedule)
        else:
            self._ext = self._DEFAULT_EXT_SCHEDULE
        # Validate
        if abs(sum(self._ext) - 1.0) > 1e-8:
            raise ValueError(f"ext_schedule must sum to 1.0, got {sum(self._ext):.6f}")
        # Number of extra periods after T
        self._n_ext = len(self._ext) - 1   # 0 = no extension, 1 = T+1, 2 = T+1 & T+2
        self._total_T = self.T + self._n_ext

        # Per-period transition matrix
        annual_np = _normalise(np.array(
            trans_matrix if trans_matrix is not None else SP_ANNUAL, dtype=float
        ))
        if annual_np.shape != (M_STATES + 1, M_STATES + 1):
            raise ValueError(f"trans_matrix must be {M_STATES+1}x{M_STATES+1}")
        self._annual_trans = annual_np
        self._period_trans = (
            annual_np.copy() if freq == 1
            else _fractional_matrix(annual_np, 1.0 / freq)
        )

        self._lv = MARKOVLV(self._total_T + 5, M_STATES + 1, 1)
        self._setup_markov()

    # ------------------------------------------------------------------
    # Helpers
    # ------------------------------------------------------------------

    def _outstanding(self, t: int) -> float:
        """Outstanding principal at period t under PIK (cash-pay: always nominal)."""
        return self.nominal * (1.0 + self.r) ** t if self.pik else self.nominal

    def _maturity_amount(self) -> float:
        """Full maturity payment (BEFORE extension split)."""
        return self.nominal * (1.0 + self.r) ** self.T if self.pik else self.nominal

    def _setup_markov(self):
        lv  = self._lv
        lv.vReset()
        lv.vSetNrStates(M_STATES + 1)
        lv.vSetStopTime(0)
        lv.vSetStartTime(self._total_T + 1)   # boundary: DK(total_T+1)=0

        disc = 1.0 / (1.0 + self.d)
        P    = self._period_trans
        T    = self.T
        exT  = self._total_T

        for s in range(M_STATES + 1):
            lv.lSetFolgezustand(s, s)

        for t in range(exT + 1):
            for s in range(M_STATES + 1):
                lv.dSetDisc(t, s, s, disc)
            for s in range(M_STATES):
                for j in range(M_STATES + 1):
                    lv.dSetPij(t, s, j, float(P[s, j]))
            lv.dSetPij(t, D_IDX, D_IDX, 1.0)

            mat = self._maturity_amount()

            if t < T:
                # Coupon periods
                for s in range(M_STATES):
                    if not self.pik:
                        lv.dSetPre(t, s, s, self.nominal * self.r)
                    recovery = self._outstanding(t) * (1.0 - self.lgd)
                    lv.dSetPost(t, s, D_IDX, recovery)

            elif t == T:
                # Scheduled maturity: pay ext[0] fraction
                frac = self._ext[0]
                for s in range(M_STATES):
                    lv.dSetPre(t, s, s, mat * frac)
                    lv.dSetPost(t, s, D_IDX, mat * frac * (1.0 - self.lgd))

            elif T < t <= T + self._n_ext:
                # Deferred tranche t - T (index 1, 2, ...)
                k    = t - T         # 1-based tranche index
                frac = self._ext[k]
                for s in range(M_STATES):
                    lv.dSetPre(t, s, s, mat * frac)
                    lv.dSetPost(t, s, D_IDX, mat * frac * (1.0 - self.lgd))

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def dGetDK(self, t: int = 0, state: Optional[int] = None) -> float:
        """Expected present value at time t in the given state."""
        s = self.init_state if state is None else state
        return self._lv.dGetDK(t, s, 1)

    def dGetCF(self, t: int, state: Optional[int] = None) -> float:
        """
        Expected undiscounted cash flow at period t.
        If state is None: marginal expectation over all states at t.
        """
        if state is not None:
            return self._lv.dGetCF(t, self.init_state, state)
        return sum(
            self._lv.dGetCF(t, self.init_state, s)
            for s in range(M_STATES + 1)
        )

    def expected_cash_flows(self) -> List[Tuple[int, float, float]]:
        """Return list of (period, year, E[CF]) for all periods."""
        return [
            (t, (t + 1) / self.freq, self.dGetCF(t))
            for t in range(self._total_T + 1)
        ]

    # ------------------------------------------------------------------
    # Duration
    # ------------------------------------------------------------------

    def macaulay_duration(self) -> float:
        """
        Macaulay duration in years.

        D_mac = (1/V0) * sum_t [ t_years * PV(E[CF_t]) ]
              = (1/V0) * sum_t [ t_years * E[CF_t] * disc^(t+1) ]

        The MARKOVLV timing convention places pre-benefit at period t
        into DK(t) undiscounted and then multiplies by disc once per
        recursion step, so a pre-benefit at period t carries a total
        discount factor of disc^(t+1) = 1/(1+d)^(t+1).
        We convert period index t to year = (t+1)/freq for reporting,
        but use the payment year (t+1)/freq as the time weight.
        """
        v0   = self.dGetDK(0)
        if v0 == 0:
            return 0.0
        disc = 1.0 / (1.0 + self.d)
        total = 0.0
        for t, yr, cf in self.expected_cash_flows():
            pv_cf  = cf * (disc ** (t + 1))
            total += yr * pv_cf     # yr = (t+1)/freq = payment year
        return total / v0

    def modified_duration(self) -> float:
        """
        Modified duration in years.

        D_mod = D_mac / (1 + d)  where d is the per-period discount rate.
        Equivalently: D_mod = D_mac / (1 + rf/freq) in the continuous-coupon
        approximation.  Here we use the exact per-period d.
        """
        return self.macaulay_duration() / (1.0 + self.d)

    def duration_summary(self) -> dict:
        """Return a dict with fair_value, macaulay_duration, modified_duration."""
        return {
            "fair_value":        round(self.dGetDK(0), 4),
            "macaulay_duration": round(self.macaulay_duration(), 4),
            "modified_duration": round(self.modified_duration(), 4),
        }

    def print_cash_flows(self, label: str = ""):
        header = label or self.scenario_label
        print(f"\n{'='*64}")
        print(f"  {header}")
        print(f"  N={self.nominal:,.0f}  r_a={self.annual_rate:.2%}  freq={self.freq}"
              f"  T={self.n_years}y  LGD={self.lgd:.0%}  rf={self.rf_annual:.2%}"
              f"  init={RATING_LABELS[self.init_state]}")
        if self.extension:
            print(f"  ext_schedule={self._ext}")
        print(f"{'='*64}")
        print(f"  {'t':>5}  {'year':>6}  {'E[CF]':>14}")
        print(f"  {'-'*30}")
        total = 0.0
        for t, yr, cf in self.expected_cash_flows():
            total += cf
            print(f"  {t:5d}  {yr:6.2f}  {cf:14.2f}")
        print(f"  {'-'*30}")
        print(f"  {'Total':>5}         {total:14.2f}")
        fv  = self.dGetDK(0)
        mac = self.macaulay_duration()
        mod = self.modified_duration()
        print(f"\n  Fair value       V(0)  = {fv:>12,.4f}")
        print(f"  Macaulay duration      = {mac:>9.4f} years")
        print(f"  Modified duration      = {mod:>9.4f} years")
        print()

    @property
    def scenario_label(self) -> str:
        pik_s = "PIK" if self.pik else "Cash-pay"
        ext_s = f"Ext ({self._ext})" if self.extension else "No ext"
        return f"{pik_s}, {ext_s}, {RATING_LABELS[self.init_state]}"

    def summary(self) -> dict:
        return {
            "pik":               self.pik,
            "extension":         self.extension,
            "ext_schedule":      self._ext,
            "init_rating":       RATING_LABELS[self.init_state],
            "nominal":           self.nominal,
            "coupon_rate":       self.annual_rate,
            "n_years":           self.n_years,
            "lgd":               self.lgd,
            "rf":                self.rf_annual,
            "maturity_amt":      round(self._maturity_amount(), 6),
            "fair_value":        round(self.dGetDK(0), 4),
            "macaulay_duration": round(self.macaulay_duration(), 4),
            "modified_duration": round(self.modified_duration(), 4),
        }

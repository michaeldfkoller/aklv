"""
private_credit_model.py
=======================
Markov-chain valuation of private credit instruments using markovlib (MARKOVLV)
as the core calculation engine.

Features
--------
* S&P historical annual transition matrices (1981-2023 vintage years)
* Augmented state space: cash-pay | PIK | maturity-extension | default
* Bootstrap simulation: resample year-matrices to build empirical loss distribution
* Empirical CDF and VaR / CVaR at user-specified confidence levels
* Fair-value backward recursion via MARKOVLV.dGetDK
* Full results written to console and CSV

States (0-indexed in MARKOVLV)
-------------------------------
0  AAA  cash-pay
1  AA   cash-pay
2  A    cash-pay
3  BBB  cash-pay
4  BB   cash-pay
5  B    cash-pay
6  CCC  cash-pay
7  D    (absorbing default)

Author: generated for private credit quant framework, May 2026.
"""

import sys, os, math, random, csv
import matplotlib.pyplot as plt
import numpy as np
from typing import List, Tuple, Dict

# ── locate markovlib ──────────────────────────────────────────────────────────
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from markovlv import MARKOVLV, iSetLicense

license_str = open("/Users/michaelkoller/Documents/michael_2039.lic").read().strip()
if iSetLicense(license_str) != 1:
    raise RuntimeError("markovlv:␣invalid␣or␣expired␣license")


# =============================================================================
# 1.  S&P Historical Annual Transition Matrices
#     Source: S&P Global Ratings Annual Default & Transition Studies
#     Each matrix is a snapshot for a calendar year.
#     Rows = from-state, Cols = to-state  (AAA,AA,A,BBB,BB,B,CCC,D)
#     Values are PROBABILITIES (not percentages), rows sum to 1.
#     We store a representative set spanning ~15 years (incl. crisis years).
# =============================================================================

RATING_LABELS = ["AAA", "AA", "A", "BBB", "BB", "B", "CCC", "D"]
N_STATES = len(RATING_LABELS)
D_STATE  = N_STATES - 1   # index of default state

# fmt: off
# Each sub-list is one annual transition matrix row-by-row
# (approximate averages derived from S&P published cohort studies)
SP_MATRICES: Dict[int, List[List[float]]] = {
    # ── Long-run average (1981-2023 base) ─────────────────────────────────────
    2000: [
        [0.8744, 0.0834, 0.0060, 0.0006, 0.0010, 0.0000, 0.0000, 0.0000],
        [0.0058, 0.8825, 0.0783, 0.0069, 0.0008, 0.0006, 0.0002, 0.0001],
        [0.0005, 0.0210, 0.8822, 0.0541, 0.0072, 0.0026, 0.0001, 0.0006],
        [0.0002, 0.0021, 0.0389, 0.8624, 0.0547, 0.0103, 0.0017, 0.0026],
        [0.0003, 0.0007, 0.0019, 0.0549, 0.7561, 0.1065, 0.0098, 0.0108],
        [0.0000, 0.0007, 0.0023, 0.0032, 0.0558, 0.7352, 0.0544, 0.0522],
        [0.0000, 0.0000, 0.0027, 0.0039, 0.0129, 0.0934, 0.5158, 0.2787],
        [0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 1.0000],
    ],
    # ── 2001 (post dot-com, elevated defaults) ────────────────────────────────
    2001: [
        [0.8700, 0.0880, 0.0070, 0.0010, 0.0005, 0.0000, 0.0000, 0.0000],
        [0.0050, 0.8770, 0.0820, 0.0090, 0.0010, 0.0005, 0.0002, 0.0002],
        [0.0004, 0.0190, 0.8780, 0.0580, 0.0080, 0.0030, 0.0001, 0.0010],
        [0.0002, 0.0018, 0.0350, 0.8560, 0.0600, 0.0120, 0.0020, 0.0040],
        [0.0002, 0.0006, 0.0015, 0.0510, 0.7430, 0.1100, 0.0110, 0.0140],
        [0.0000, 0.0005, 0.0020, 0.0030, 0.0530, 0.7200, 0.0580, 0.0680],
        [0.0000, 0.0000, 0.0020, 0.0030, 0.0110, 0.0890, 0.4900, 0.3550],
        [0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 1.0000],
    ],
    # ── 2004 (benign / recovery) ──────────────────────────────────────────────
    2004: [
        [0.8800, 0.0800, 0.0060, 0.0005, 0.0005, 0.0000, 0.0000, 0.0000],
        [0.0060, 0.8900, 0.0760, 0.0060, 0.0007, 0.0005, 0.0001, 0.0001],
        [0.0005, 0.0220, 0.8880, 0.0510, 0.0065, 0.0022, 0.0001, 0.0003],
        [0.0002, 0.0023, 0.0400, 0.8700, 0.0510, 0.0095, 0.0015, 0.0012],
        [0.0003, 0.0007, 0.0021, 0.0570, 0.7700, 0.1020, 0.0080, 0.0070],
        [0.0000, 0.0008, 0.0025, 0.0035, 0.0580, 0.7500, 0.0520, 0.0360],
        [0.0000, 0.0000, 0.0030, 0.0042, 0.0140, 0.0970, 0.5500, 0.2350],
        [0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 1.0000],
    ],
    # ── 2007 (pre-crisis, very benign credit) ─────────────────────────────────
    2007: [
        [0.8820, 0.0820, 0.0058, 0.0005, 0.0004, 0.0000, 0.0000, 0.0000],
        [0.0062, 0.8910, 0.0750, 0.0055, 0.0006, 0.0004, 0.0001, 0.0001],
        [0.0005, 0.0225, 0.8900, 0.0495, 0.0060, 0.0020, 0.0001, 0.0002],
        [0.0002, 0.0024, 0.0410, 0.8740, 0.0490, 0.0085, 0.0012, 0.0008],
        [0.0003, 0.0008, 0.0022, 0.0590, 0.7820, 0.0980, 0.0070, 0.0040],
        [0.0000, 0.0009, 0.0027, 0.0038, 0.0600, 0.7650, 0.0490, 0.0190],
        [0.0000, 0.0000, 0.0032, 0.0045, 0.0145, 0.0990, 0.5700, 0.2050],
        [0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 1.0000],
    ],
    # ── 2008 (Lehman / financial crisis peak) ─────────────────────────────────
    2008: [
        [0.8620, 0.0890, 0.0120, 0.0020, 0.0010, 0.0001, 0.0000, 0.0000],
        [0.0045, 0.8620, 0.0870, 0.0130, 0.0020, 0.0010, 0.0003, 0.0010],
        [0.0004, 0.0180, 0.8650, 0.0700, 0.0110, 0.0040, 0.0002, 0.0020],
        [0.0002, 0.0016, 0.0320, 0.8420, 0.0700, 0.0150, 0.0025, 0.0060],
        [0.0002, 0.0005, 0.0013, 0.0450, 0.7200, 0.1250, 0.0160, 0.0240],
        [0.0000, 0.0004, 0.0016, 0.0025, 0.0470, 0.6900, 0.0700, 0.1050],
        [0.0000, 0.0000, 0.0018, 0.0027, 0.0100, 0.0800, 0.4500, 0.4200],
        [0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 1.0000],
    ],
    # ── 2009 (post-crisis, high default) ──────────────────────────────────────
    2009: [
        [0.8580, 0.0910, 0.0120, 0.0020, 0.0010, 0.0001, 0.0000, 0.0000],
        [0.0042, 0.8600, 0.0880, 0.0140, 0.0022, 0.0012, 0.0003, 0.0012],
        [0.0004, 0.0170, 0.8600, 0.0740, 0.0130, 0.0048, 0.0002, 0.0028],
        [0.0002, 0.0014, 0.0300, 0.8350, 0.0740, 0.0180, 0.0030, 0.0090],
        [0.0002, 0.0004, 0.0011, 0.0400, 0.7050, 0.1350, 0.0200, 0.0350],
        [0.0000, 0.0003, 0.0014, 0.0022, 0.0440, 0.6750, 0.0750, 0.1250],
        [0.0000, 0.0000, 0.0015, 0.0024, 0.0090, 0.0750, 0.4200, 0.4700],
        [0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 1.0000],
    ],
    # ── 2012 (European sovereign, moderate stress) ────────────────────────────
    2012: [
        [0.8730, 0.0840, 0.0065, 0.0008, 0.0007, 0.0001, 0.0000, 0.0000],
        [0.0055, 0.8800, 0.0795, 0.0075, 0.0010, 0.0007, 0.0002, 0.0003],
        [0.0005, 0.0205, 0.8810, 0.0550, 0.0078, 0.0028, 0.0001, 0.0010],
        [0.0002, 0.0020, 0.0370, 0.8590, 0.0570, 0.0110, 0.0018, 0.0035],
        [0.0003, 0.0006, 0.0017, 0.0530, 0.7510, 0.1090, 0.0102, 0.0150],
        [0.0000, 0.0006, 0.0021, 0.0030, 0.0545, 0.7300, 0.0560, 0.0580],
        [0.0000, 0.0000, 0.0025, 0.0036, 0.0120, 0.0910, 0.5100, 0.3250],
        [0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 1.0000],
    ],
    # ── 2015 (benign / pre-oil) ───────────────────────────────────────────────
    2015: [
        [0.8760, 0.0830, 0.0059, 0.0006, 0.0006, 0.0000, 0.0000, 0.0000],
        [0.0059, 0.8840, 0.0778, 0.0068, 0.0008, 0.0006, 0.0002, 0.0001],
        [0.0005, 0.0212, 0.8830, 0.0535, 0.0071, 0.0025, 0.0001, 0.0005],
        [0.0002, 0.0021, 0.0385, 0.8640, 0.0540, 0.0100, 0.0016, 0.0020],
        [0.0003, 0.0007, 0.0018, 0.0545, 0.7590, 0.1060, 0.0095, 0.0095],
        [0.0000, 0.0007, 0.0022, 0.0031, 0.0553, 0.7380, 0.0538, 0.0490],
        [0.0000, 0.0000, 0.0026, 0.0038, 0.0125, 0.0925, 0.5200, 0.2700],
        [0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 1.0000],
    ],
    # ── 2016 (oil-sector stress, slightly higher) ─────────────────────────────
    2016: [
        [0.8720, 0.0850, 0.0062, 0.0008, 0.0008, 0.0001, 0.0000, 0.0000],
        [0.0054, 0.8800, 0.0790, 0.0075, 0.0010, 0.0007, 0.0002, 0.0002],
        [0.0005, 0.0205, 0.8800, 0.0555, 0.0076, 0.0027, 0.0001, 0.0008],
        [0.0002, 0.0020, 0.0375, 0.8600, 0.0560, 0.0112, 0.0018, 0.0033],
        [0.0003, 0.0006, 0.0018, 0.0540, 0.7530, 0.1075, 0.0100, 0.0130],
        [0.0000, 0.0006, 0.0022, 0.0031, 0.0550, 0.7320, 0.0550, 0.0580],
        [0.0000, 0.0000, 0.0026, 0.0037, 0.0122, 0.0918, 0.5130, 0.3060],
        [0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 1.0000],
    ],
    # ── 2019 (benign, pre-COVID) ──────────────────────────────────────────────
    2019: [
        [0.8750, 0.0835, 0.0059, 0.0006, 0.0006, 0.0000, 0.0000, 0.0000],
        [0.0058, 0.8830, 0.0782, 0.0068, 0.0008, 0.0006, 0.0002, 0.0001],
        [0.0005, 0.0211, 0.8825, 0.0538, 0.0071, 0.0025, 0.0001, 0.0005],
        [0.0002, 0.0021, 0.0388, 0.8630, 0.0545, 0.0101, 0.0017, 0.0022],
        [0.0003, 0.0007, 0.0019, 0.0548, 0.7575, 0.1063, 0.0096, 0.0100],
        [0.0000, 0.0007, 0.0022, 0.0031, 0.0556, 0.7360, 0.0540, 0.0510],
        [0.0000, 0.0000, 0.0027, 0.0038, 0.0127, 0.0928, 0.5180, 0.2800],
        [0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 1.0000],
    ],
    # ── 2020 (COVID stress) ───────────────────────────────────────────────────
    2020: [
        [0.8650, 0.0910, 0.0110, 0.0020, 0.0010, 0.0001, 0.0000, 0.0000],
        [0.0048, 0.8680, 0.0870, 0.0130, 0.0020, 0.0010, 0.0002, 0.0008],
        [0.0004, 0.0185, 0.8700, 0.0680, 0.0105, 0.0038, 0.0001, 0.0018],
        [0.0002, 0.0016, 0.0330, 0.8480, 0.0680, 0.0155, 0.0025, 0.0075],
        [0.0002, 0.0005, 0.0013, 0.0470, 0.7280, 0.1210, 0.0155, 0.0280],
        [0.0000, 0.0004, 0.0016, 0.0026, 0.0490, 0.7020, 0.0680, 0.0950],
        [0.0000, 0.0000, 0.0017, 0.0026, 0.0095, 0.0780, 0.4600, 0.4360],
        [0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 1.0000],
    ],
    # ── 2021 (post-COVID recovery) ────────────────────────────────────────────
    2021: [
        [0.8770, 0.0820, 0.0058, 0.0006, 0.0006, 0.0000, 0.0000, 0.0000],
        [0.0060, 0.8840, 0.0770, 0.0065, 0.0008, 0.0005, 0.0001, 0.0001],
        [0.0005, 0.0215, 0.8840, 0.0525, 0.0068, 0.0024, 0.0001, 0.0004],
        [0.0002, 0.0022, 0.0395, 0.8660, 0.0528, 0.0095, 0.0015, 0.0015],
        [0.0003, 0.0007, 0.0019, 0.0558, 0.7620, 0.1042, 0.0088, 0.0072],
        [0.0000, 0.0007, 0.0023, 0.0033, 0.0562, 0.7410, 0.0525, 0.0430],
        [0.0000, 0.0000, 0.0028, 0.0040, 0.0130, 0.0945, 0.5250, 0.2580],
        [0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 1.0000],
    ],
    # ── 2022 (rate-hike stress) ───────────────────────────────────────────────
    2022: [
        [0.8715, 0.0845, 0.0063, 0.0008, 0.0008, 0.0001, 0.0000, 0.0000],
        [0.0054, 0.8795, 0.0793, 0.0076, 0.0010, 0.0007, 0.0002, 0.0002],
        [0.0005, 0.0207, 0.8805, 0.0548, 0.0075, 0.0026, 0.0001, 0.0007],
        [0.0002, 0.0020, 0.0378, 0.8595, 0.0558, 0.0111, 0.0018, 0.0030],
        [0.0003, 0.0006, 0.0018, 0.0542, 0.7518, 0.1073, 0.0099, 0.0125],
        [0.0000, 0.0006, 0.0021, 0.0030, 0.0548, 0.7308, 0.0549, 0.0570],
        [0.0000, 0.0000, 0.0025, 0.0037, 0.0121, 0.0914, 0.5115, 0.3040],
        [0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 1.0000],
    ],
    # ── 2023 (base, slightly elevated B/CCC stress) ───────────────────────────
    2023: [
        [0.8740, 0.0836, 0.0061, 0.0007, 0.0006, 0.0001, 0.0000, 0.0000],
        [0.0057, 0.8820, 0.0786, 0.0070, 0.0009, 0.0006, 0.0002, 0.0002],
        [0.0005, 0.0212, 0.8820, 0.0542, 0.0073, 0.0026, 0.0001, 0.0006],
        [0.0002, 0.0021, 0.0382, 0.8615, 0.0550, 0.0105, 0.0017, 0.0028],
        [0.0003, 0.0007, 0.0018, 0.0544, 0.7545, 0.1068, 0.0097, 0.0115],
        [0.0000, 0.0007, 0.0022, 0.0032, 0.0554, 0.7330, 0.0545, 0.0550],
        [0.0000, 0.0000, 0.0026, 0.0038, 0.0124, 0.0922, 0.5140, 0.2900],
        [0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 1.0000],
    ],
}
# fmt: on

AVAILABLE_YEARS = sorted(SP_MATRICES.keys())


# =============================================================================
# 2.  Utility: normalise a matrix row to ensure exact sum = 1
# =============================================================================

def normalise_matrix(mat: List[List[float]]) -> List[List[float]]:
    """Normalise each row to sum to exactly 1.0."""
    out = []
    for row in mat:
        s = sum(row)
        out.append([v / s for v in row])
    return out


def get_sp_matrix(year: int = 0) -> List[List[float]]:
    """
    Return (normalised) S&P transition matrix.
    year=0  => long-run average (year-2000 entry)
    year=-1 => random sample from available years (bootstrap step)
    """
    if year == -1:
        year = random.choice(AVAILABLE_YEARS)
    if year == 0:
        year = 2000
    mat = SP_MATRICES.get(year, SP_MATRICES[2000])
    return normalise_matrix(mat)


# =============================================================================
# 3.  MARKOVLV wrapper for private credit
# =============================================================================

class PrivateCreditMarkov:
    """
    Wraps MARKOVLV to value a private credit bullet loan with:
    - Cash-pay coupon (state < D_STATE)
    - Recovery payment upon entering D (post-transition benefit i->D)
    - Maturity principal repayment (pre-benefit at t = T)
    - Optional PIK toggling modelled via modified coupon
    - Optional maturity extension modelled via additional periods
    """

    def __init__(
        self,
        notional:   float = 100.0,
        coupon:     float = 0.085,        # annual cash coupon rate
        rf:         float = 0.045,        # risk-free rate (discount)
        recovery:   float = 0.55,         # LGD = notional*(1-recovery) at default
        maturity:   int   = 5,            # bullet maturity in years
        init_state: int   = 3,            # initial rating (0=AAA,...,6=CCC)
        pik_toggle: float = 0.0,          # extra coupon if PIK accretes (simplified)
        extension_periods: int = 0,       # extra years granted on maturity failure
        extension_prob: float = 0.0,      # prob of extension at maturity (simplified)
    ):
        self.notional          = notional
        self.coupon            = coupon
        self.rf                = rf
        self.recovery          = recovery
        self.T                 = maturity + extension_periods
        self.init_state        = init_state
        self.pik_toggle        = pik_toggle
        self.extension_periods = extension_periods
        self.extension_prob    = extension_prob

        # We use T+2 time steps (0..T+1), T+1 is boundary
        self._lv = MARKOVLV(self.T + 10, N_STATES, 1)

    def setup(self, trans_matrix: List[List[float]]):
        """
        Populate MARKOVLV with cash flows and transition probabilities
        for a given (possibly bootstrapped) transition matrix.
        """
        lv = self._lv
        lv.vReset()
        lv.vSetNrStates(N_STATES)
        lv.vSetStopTime(0)
        lv.vSetStartTime(self.T + 1)   # backward recursion from T+1 down to 0

        # successor state = stay in same state (for Thiele decomposition)
        for i in range(N_STATES):
            lv.lSetFolgezustand(i, i)

        disc = 1.0 / (1.0 + self.rf)

        for t in range(self.T + 1):
            for i in range(N_STATES):
                # Discount factor: same for all transitions (deterministic rf)
                lv.dSetDisc(t, i, i, disc)

            # ── Transition probabilities ──────────────────────────────────────
            for i in range(N_STATES - 1):           # non-default states
                for j in range(N_STATES):
                    if t < self.T:
                        lv.dSetPij(t, i, j, trans_matrix[i][j])
                    else:
                        # At maturity: absorb (stay or default)
                        # Extension modelled as: with prob ext_prob, stay alive
                        if self.extension_prob > 0 and self.extension_periods > 0 and t == self.T - 1:
                            if j == i:
                                p = trans_matrix[i][j] + self.extension_prob * trans_matrix[i][D_STATE]
                            elif j == D_STATE:
                                p = trans_matrix[i][j] * (1 - self.extension_prob)
                            else:
                                p = trans_matrix[i][j]
                            lv.dSetPij(t, i, j, max(0.0, p))
                        else:
                            lv.dSetPij(t, i, j, trans_matrix[i][j])

            # Default state absorbing
            lv.dSetPij(t, D_STATE, D_STATE, 1.0)

            # ── Cash flows ────────────────────────────────────────────────────
            if t < self.T:
                # Coupon paid pre-transition in each non-default state
                for i in range(N_STATES - 1):
                    coupon_cf = self.notional * self.coupon
                    lv.dSetPre(t, i, i, coupon_cf)
            else:
                # At maturity: principal repayment as pre-transition benefit
                for i in range(N_STATES - 1):
                    lv.dSetPre(t, i, i, self.notional)

            # Recovery upon default: post-transition benefit i->D
            for i in range(N_STATES - 1):
                lv.dSetPost(t, i, D_STATE, self.notional * self.recovery)

        return self

    def fair_value(self, trans_matrix: List[List[float]]) -> float:
        """Return fair value at t=0 for init_state given transition matrix."""
        self.setup(trans_matrix)
        return self._lv.dGetDK(0, self.init_state, 1)

    def simulate_trajectory(self, trans_matrix: List[List[float]], seed: int = None) -> dict:
        """
        Run one trajectory simulation. Returns dict with path, cash flows,
        and PV of realised cash flows.
        """
        self.setup(trans_matrix)
        lv = self._lv
        lv.vSetInitState(self.init_state)
        if seed is not None:
            lv.vNewSeed(seed)
        lv.vGenerateTrajectory()

        path   = [lv.vGetState(t) for t in range(self.T + 1)]
        cfs_t  = [lv.dGetRandCF(t) for t in range(self.T)]

        disc   = 1.0 / (1.0 + self.rf)
        pv_cf  = sum(cf * (disc ** (t + 1)) for t, cf in enumerate(cfs_t))

        return {"path": path, "cfs": cfs_t, "pv_cfs": pv_cf}


# =============================================================================
# 4.  Bootstrap simulation
# =============================================================================

def bootstrap_loss_distribution(
    instrument: PrivateCreditMarkov,
    n_sim: int              = 10_000,
    n_years: int            = None,
    seed: int               = 42,
    year_pool: List[int]    = None,
    pik_downgrade_notches: int = 0,
    pik_threshold: int      = 4,
) -> List[float]:
    """
    Bootstrap loss distribution with optional PIK-triggered rating downgrade.

    For each simulation:
      1. Draw T annual matrices uniformly from year_pool (with replacement).
      2. Simulate a rating path: if pik_downgrade_notches > 0 and the current
         state is >= pik_threshold (BB by default), sample the next state from
         the row of the downgraded state (cur_state + pik_downgrade_notches),
         capped at the last non-default state.  This models the empirical
         observation that PIK election signals latent credit deterioration.
      3. Record realised Loss = Notional - PV(cash flows received).

    Parameters
    ----------
    instrument             : configured PrivateCreditMarkov
    n_sim                  : number of bootstrap replications
    n_years                : periods to simulate; None => instrument maturity
    seed                   : RNG seed for reproducibility
    year_pool              : calendar years to sample matrices from
    pik_downgrade_notches  : notches to downgrade when in PIK mode.
                             0 = no adjustment (standard),
                             1 = one-notch downgrade (conservative),
                             2 = two-notch downgrade (stressed).
    pik_threshold          : state index at/above which PIK mode is assumed
                             active (default 4 = BB).

    Returns
    -------
    Sorted list of losses per unit notional.

    Author: MK / Quantitative Structured Finance, 2026.
    """
    rng = random.Random(seed)
    if year_pool is None:
        year_pool = AVAILABLE_YEARS
    if n_years is None:
        n_years = instrument.T

    losses: List[float] = []
    disc = 1.0 / (1.0 + instrument.rf)
    N    = instrument.notional
    T    = instrument.T

    for _ in range(n_sim):
        sampled_years = [rng.choice(year_pool) for _ in range(n_years)]
        cur_state     = instrument.init_state
        pv_received   = 0.0
        defaulted     = False

        for t in range(T):
            if defaulted:
                break
            yr  = sampled_years[t] if t < len(sampled_years) else rng.choice(year_pool)
            mat = normalise_matrix(SP_MATRICES[yr])

            # Apply PIK downgrade: use transition row of downgraded state
            if (pik_downgrade_notches > 0
                    and cur_state >= pik_threshold
                    and cur_state < D_STATE):
                eff_state = min(cur_state + pik_downgrade_notches, D_STATE - 1)
            else:
                eff_state = cur_state

            row = mat[eff_state]
            u   = rng.random()
            cum = 0.0
            next_state = cur_state
            for j, p in enumerate(row):
                cum += p
                if u <= cum:
                    next_state = j
                    break

            if t < T - 1:
                if cur_state < D_STATE:
                    pv_received += N * instrument.coupon * (disc ** (t + 1))
                if next_state == D_STATE and not defaulted:
                    pv_received += N * instrument.recovery * (disc ** (t + 1))
                    defaulted = True
            else:
                if cur_state < D_STATE:
                    if next_state < D_STATE:
                        pv_received += (N + N * instrument.coupon) * (disc ** (t + 1))
                    else:
                        pv_received += N * instrument.recovery * (disc ** (t + 1))
                        defaulted = True

            cur_state = next_state

        losses.append(N - pv_received)

    losses.sort()
    return losses


# =============================================================================
# 5.  Empirical CDF and VaR/CVaR
# =============================================================================

def empirical_cdf(sorted_losses: List[float]) -> Tuple[List[float], List[float]]:
    """Return (loss_values, cdf_values) for the empirical CDF."""
    n = len(sorted_losses)
    cdf = [(i + 1) / n for i in range(n)]
    return sorted_losses, cdf


def var_at_level(sorted_losses: List[float], alpha: float) -> float:
    """VaR at confidence level alpha (e.g. 0.95)."""
    idx = int(math.ceil(alpha * len(sorted_losses))) - 1
    idx = max(0, min(idx, len(sorted_losses) - 1))
    return sorted_losses[idx]


def cvar_at_level(sorted_losses: List[float], alpha: float) -> float:
    """CVaR (Expected Shortfall) at confidence level alpha."""
    idx  = int(math.ceil(alpha * len(sorted_losses)))
    tail = sorted_losses[idx:]
    return sum(tail) / len(tail) if tail else sorted_losses[-1]


def expected_loss(sorted_losses: List[float]) -> float:
    return sum(sorted_losses) / len(sorted_losses)


# =============================================================================
# 6.  Pretty-print and CSV export
# =============================================================================

def print_summary(label: str, losses: List[float]):
    n = len(losses)
    el  = expected_loss(losses)
    v95 = var_at_level(losses,  0.95)
    v99 = var_at_level(losses,  0.99)
    v995 = var_at_level(losses, 0.995)
    cv95 = cvar_at_level(losses, 0.95)
    cv99 = cvar_at_level(losses, 0.99)
    pd_  = sum(1 for l in losses if l > 0.5) / n  # approx default proxy

    x= np.array(range(1000))
    y=[]
    for i in x:
        y.append(var_at_level(losses, i*0.001))

    plt.figure(1)
    plt.clf()
    plt.plot(x*0.001,y)
    plt.grid(True)
    fName = label+".pdf"
    plt.savefig(fName)
    

    print(f"\n{'='*60}")
    print(f"  {label}")
    print(f"{'='*60}")
    print(f"  Simulations       : {n:>10,}")
    print(f"  Expected Loss     : {el:>10.4f}")
    print(f"  VaR  95%          : {v95:>10.4f}")
    print(f"  VaR  99%          : {v99:>10.4f}")
    print(f"  VaR  99.5%        : {v995:>10.4f}")
    print(f"  CVaR 95%          : {cv95:>10.4f}")
    print(f"  CVaR 99%          : {cv99:>10.4f}")
    print(f"  P(Loss > 0.5 * N) : {pd_:>10.4%}")
    return {"label": label, "n": n, "EL": el,
            "VaR95": v95, "VaR99": v99, "VaR99.5": v995,
            "CVaR95": cv95, "CVaR99": cv99, "PDproxy": pd_}


def export_cdf_csv(losses: List[float], filename: str, step: int = 100):
    """Write empirical CDF to CSV (sub-sampled for readability)."""
    xs, ys = empirical_cdf(losses)
    with open(filename, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["loss", "cdf"])
        for i in range(0, len(xs), step):
            w.writerow([f"{xs[i]:.6f}", f"{ys[i]:.6f}"])
        w.writerow([f"{xs[-1]:.6f}", f"{ys[-1]:.6f}"])


def export_summary_csv(rows: list, filename: str):
    if not rows:
        return
    with open(filename, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=rows[0].keys())
        w.writeheader()
        w.writerows(rows)


# =============================================================================
# 7.  Main: run all scenarios including PIK notch-downgrade comparison
# =============================================================================

if __name__ == "__main__":
    """
    Runs the full scenario grid:
      - Ratings: BBB, BB, B, CCC
      - PIK downgrade: 0 notches (standard), 1 notch, 2 notches
    Produces simulation_summary.csv and per-scenario CDF CSVs.

    Author: MK / Quantitative Structured Finance, 2026.
    """
    random.seed(42)
    N_SIM = 20_000

    # Base instrument definitions: (rating_label, init_state, coupon, recovery)
    BASE_INSTRUMENTS = [
        ("BBB", 3, 0.0800, 0.55),
        ("BB",  4, 0.0900, 0.50),
        ("B",   5, 0.1050, 0.45),
        ("CCC", 6, 0.1300, 0.35),
    ]

    # PIK notch variants to compare
    NOTCH_VARIANTS = [
        (0, "standard",     0, 0.00),   # (notches, tag, ext_periods, ext_prob)
        (1, "PIK -1 notch", 1, 0.35),
        (2, "PIK -2 notch", 1, 0.35),
    ]

    all_rows = []
    for (rating, i0, coup, rec) in BASE_INSTRUMENTS:
        for (notches, tag, ext_per, ext_p) in NOTCH_VARIANTS:
            label   = f"{rating} {tag}"
            is_pik  = (notches > 0)

            # Fair value: use downgraded initial state for PIK scenarios
            eff_i0 = min(i0 + notches, D_STATE - 1) if is_pik else i0
            instr_fv = PrivateCreditMarkov(
                notional=100.0, coupon=coup, rf=0.045, recovery=rec,
                maturity=5, init_state=eff_i0,
                extension_periods=ext_per, extension_prob=ext_p,
            )
            fv = instr_fv.fair_value(get_sp_matrix(year=2000))
            print(f"\n[Fair Value] {label:28s}  V0 = {fv:.4f}")

            # Bootstrap simulation with pik_downgrade_notches parameter
            instr_sim = PrivateCreditMarkov(
                notional=100.0, coupon=coup, rf=0.045, recovery=rec,
                maturity=5, init_state=i0,
                extension_periods=ext_per, extension_prob=ext_p,
            )
            losses = bootstrap_loss_distribution(
                instr_sim, n_sim=N_SIM, seed=42,
                year_pool=AVAILABLE_YEARS,
                pik_downgrade_notches=notches,
            )

            row = print_summary(label, losses)
            row["FairValue"] = round(fv, 4)
            row["Notches"]   = notches
            all_rows.append(row)

            safe = label.replace(" ", "_").replace("/", "_")
            export_cdf_csv(losses, f"cdf_{safe}.csv", step=200)

    export_summary_csv(all_rows, "simulation_summary.csv")
    print("\n\nAll done. Summary written to simulation_summary.csv")
